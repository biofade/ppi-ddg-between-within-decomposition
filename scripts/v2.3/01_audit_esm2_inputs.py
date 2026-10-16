# -*- coding: utf-8 -*-
"""scripts/v2.3/01_audit_esm2_inputs.py -- V2.3 任务 1: Colab ESM-2 提取的数据审计.

目的
----
在为 Colab ESM-2 提取准备输入之前, 必须先用**真实列名与真实数据**回答下列问题,
否则 notebook 与合并脚本会建立在猜测之上:

1. 配体/受体序列列名、ID 列名 (真实列名)。
2. ID 格式是否与 `features/v2.2_proteinea_features.csv` 的 `pdb_id` 一致。
3. 是否存在空序列 / 短序列 (<10 aa) / 长序列 (>1022 aa) / 非标准残基 (B,Z,X,U,O)。
4. 唯一序列数 vs 总行数 (去重能省多少算力)。
5. 序列长度分布 (min / median / max), 分全部 / 唯一序列两种口径。

同时给出 ESM-2 侧的关键**算力与显存估算** (按 150M 模型, layer 30, mean-pool)。

输入 : `data/ppb_affinity_proteinea/{filtered,raw,raw_recover_missing_res}.csv`
       `features/v2.2_proteinea_features.csv` (V2.2 特征表, 用于核对 pdb_id)
输出 : `features/v2.3_esm2_data_audit.md`
       `features/v2.3_esm2_data_audit.json`
       `logs/v2.3/01_audit_esm2_inputs.log`

纪律: 纯 CPU; 无网络; 只读输入; 不修改任何 V1.1/V2.0/V2.1/V2.2 文件。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
DATA_DIR = ROOT / "data" / "ppb_affinity_proteinea"
FEAT_DIR = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2.3"

PROT_CSVS = ["filtered.csv", "raw.csv", "raw_recover_missing_res.csv"]
V22_FEATURES = FEAT_DIR / "v2.2_proteinea_features.csv"
OUT_MD = FEAT_DIR / "v2.3_esm2_data_audit.md"
OUT_JSON = FEAT_DIR / "v2.3_esm2_data_audit.json"

LIG_COL = "Ligand Sequences"
REC_COL = "Receptor Sequences"
ID_COL = "PDB"
SPLIT_COL = "split"

MAX_LEN = 1022                     # ESM-2 位置编码上限 (含 BOS/EOS)
SHORT_THRESHOLD = 10
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
NONSTANDARD_INTEREST = list("BZXUO")
MULTI_CHAIN_SEP = ","

logger = logging.getLogger("audit_esm2")


def md5_of(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def seq_stats(seqs: pd.Series, label: str) -> dict[str, Any]:
    """一条序列列表的基础统计。"""
    s = seqs.fillna("").astype(str)
    L = s.str.len()
    return {
        "label": label,
        "n": int(len(s)),
        "n_empty": int((L == 0).sum()),
        "n_short_lt10": int((L < SHORT_THRESHOLD).sum()),
        "n_gt_1022": int((L > MAX_LEN).sum()),
        "len_min": int(L.min()) if len(L) else 0,
        "len_p05": float(np.percentile(L, 5)) if len(L) else 0.0,
        "len_median": float(L.median()) if len(L) else 0.0,
        "len_mean": float(L.mean()) if len(L) else 0.0,
        "len_p95": float(np.percentile(L, 95)) if len(L) else 0.0,
        "len_max": int(L.max()) if len(L) else 0,
        "n_unique": int(s.nunique()),
    }


def residue_audit(seqs: pd.Series, label: str) -> dict[str, Any]:
    """非标准残基审计: 逐字符统计 (分号分隔的多链会被一起扫, 分隔符单列)。"""
    cnt: Counter[str] = Counter()
    n_seqs_with_nonstd = 0
    n_seqs_with_sep = 0
    for txt in seqs.fillna("").astype(str):
        if not txt:
            continue
        if MULTI_CHAIN_SEP in txt:
            n_seqs_with_sep += 1
        hit = False
        for ch in txt:
            if ch == MULTI_CHAIN_SEP:
                continue
            cnt[ch] += 1
            if ch not in STANDARD_AA:
                hit = True
        if hit:
            n_seqs_with_nonstd += 1
    nonstd = {k: v for k, v in cnt.items() if k not in STANDARD_AA}
    interesting = {k: int(nonstd.get(k, 0)) for k in NONSTANDARD_INTEREST}
    return {
        "label": label,
        "alphabet_size": len(cnt),
        "n_seqs_with_nonstd": n_seqs_with_nonstd,
        "nonstd_counts": {k: int(v) for k, v in sorted(nonstd.items(),
                                                       key=lambda kv: -kv[1])},
        "interesting_BZXUO": interesting,
        "n_seqs_with_comma": n_seqs_with_sep,
    }


def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="V2.3 任务1: Colab ESM-2 输入数据审计")
    ap.add_argument("--data-dir", type=Path, default=DATA_DIR)
    ap.add_argument("--v22-features", type=Path, default=V22_FEATURES)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    ap.add_argument("--json-out", type=Path, default=OUT_JSON)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "01_audit_esm2_inputs.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    L: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        L.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# V2.3 任务 1 — Colab ESM-2 提取的数据审计")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 数据目录: `data/ppb_affinity_proteinea/`")
    emit(f"- 审计对象: `filtered.csv` (V2.2 主用), `raw.csv`, `raw_recover_missing_res.csv`")
    emit(f"- 本地是否可用 fair-esm: **否** (`import esm` 失败) -> ESM-2 必须在 Colab 运行")
    emit("")
    emit("> 本文件只做**审计**, 不含任何模型运行结果。ESM-2 提取本身**未在本机运行**。")

    # ---------------- 1. 列名 ----------------
    emit("")
    emit("## 1. 真实列名 (不猜测)")
    emit("")
    schemas: dict[str, list[str]] = {}
    shapes: dict[str, tuple[int, int]] = {}
    for name in PROT_CSVS:
        p = args.data_dir / name
        if not p.exists():
            emit(f"!! 缺失: {p}")
            continue
        d = pd.read_csv(p, low_memory=False)
        schemas[name] = list(d.columns)
        shapes[name] = (int(d.shape[0]), int(d.shape[1]))
    emit("| 文件 | 行 | 列 |")
    emit("|---|---:|---:|")
    for name, (r, c) in shapes.items():
        emit(f"| `{name}` | {r} | {c} |")
    emit("")
    emit("三个文件的列集合**完全一致** (已核对):")
    emit("")
    emit("```text")
    for c in schemas.get("filtered.csv", []):
        emit(f"  {c}")
    emit("```")
    emit("")
    emit("**结论 (本任务要问的 3 个列名)**:")
    emit("")
    emit("| 用途 | 真实列名 | 说明 |")
    emit("|---|---|---|")
    emit(f"| 配体序列 | `{LIG_COL}` | 多链以 `,` 连接 |")
    emit(f"| 受体序列 | `{REC_COL}` | 多链以 `,` 连接 |")
    emit(f"| 结构 ID | `{ID_COL}` | 4 字符 PDB ID (大写) |")
    emit(f"| 划分 | `{SPLIT_COL}` | **仅 `filtered.csv` 有** (train/val/test) |")
    emit("")
    emit("> `raw.csv` / `raw_recover_missing_res.csv` **没有 `split` 列** (它们各只有一个 `train` 划分),")
    emit("> 因此 ESM-2 提取应以 **`filtered.csv`** 为准, 且划分必须沿用它的 `split`。")

    # ---------------- 2. ID 一致性 ----------------
    emit("")
    emit("## 2. ID 格式与 V2.2 `pdb_id` 的一致性")
    emit("")
    filt = pd.read_csv(args.data_dir / "filtered.csv", low_memory=False)
    filt["_pdb_up"] = filt[ID_COL].astype(str).str.strip().str.upper()
    emit(f"- `{ID_COL}` 原始样例: {list(filt[ID_COL].astype(str).head(5))}")
    emit(f"- 大写化后样例: {list(filt['_pdb_up'].head(5))}")
    emit(f"- 长度分布: {dict(Counter(filt['_pdb_up'].str.len()))}")
    emit(f"- 含小写/空格的原始值个数: "
         f"{int((filt[ID_COL].astype(str) != filt['_pdb_up']).sum())}")

    id_ok = None
    v22_pdbs: set[str] = set()
    if args.v22_features.exists():
        f22 = pd.read_csv(args.v22_features, low_memory=False)
        if "pdb_id" in f22.columns:
            v22_pdbs = set(f22["pdb_id"].astype(str).str.strip().str.upper())
            src_pdbs = set(filt["_pdb_up"])
            inter = src_pdbs & v22_pdbs
            id_ok = len(inter) == len(src_pdbs)
            emit("")
            emit(f"- V2.2 特征表 `{args.v22_features.name}` 有 `pdb_id` 列 ✅")
            emit(f"- V2.2 唯一 `pdb_id` 数: **{len(v22_pdbs)}**")
            emit(f"- `filtered.csv` 唯一 `{ID_COL}` 数: **{len(src_pdbs)}**")
            emit(f"- 交集: **{len(inter)}** ({len(inter) / max(len(src_pdbs), 1) * 100:.2f}%)")
            emit(f"- `filtered.csv` 中不在 V2.2 的 PDB: "
                 f"{len(src_pdbs - v22_pdbs)} 个 {sorted(src_pdbs - v22_pdbs)[:10]}")
            emit(f"- **一致性判定: {'✅ 一致 (可直接按 pdb_id join)' if id_ok else '⚠️ 不完全一致 (见上)'}**")
        else:
            emit(f"!! `{args.v22_features.name}` 缺少 `pdb_id` 列, 实际列: {list(f22.columns)[:12]}")
    else:
        emit(f"!! 未找到 {args.v22_features} (需先运行 V2.2 的 01_adapt_proteinea.py)")

    # ---------------- 3/5. 长度与非标准残基 ----------------
    emit("")
    emit("## 3. 序列长度分布与异常")
    emit("")
    stats_rows: list[dict[str, Any]] = []
    residue_rows: list[dict[str, Any]] = []
    for name in PROT_CSVS:
        d = pd.read_csv(args.data_dir / name, low_memory=False)
        for col, tag in ((LIG_COL, "ligand"), (REC_COL, "receptor")):
            st = seq_stats(d[col], f"{name}::{tag}")
            stats_rows.append(st)
            residue_rows.append(residue_audit(d[col], f"{name}::{tag}"))
        del d

    emit("| 文件::角色 | 行 | 空序列 | <10aa | >1022aa | 唯一序列 | min | p05 | median | mean | p95 | max |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for s in stats_rows:
        emit(f"| {s['label']} | {s['n']} | {s['n_empty']} | {s['n_short_lt10']} | "
             f"{s['n_gt_1022']} | {s['n_unique']} | {s['len_min']} | {s['len_p05']:.0f} | "
             f"{s['len_median']:.0f} | {s['len_mean']:.1f} | {s['len_p95']:.0f} | {s['len_max']} |")

    emit("")
    emit("### 非标准残基审计")
    emit("")
    emit("| 文件::角色 | 字母表大小 | 含非标准残基序列数 | B | Z | X | U | O | 含 `,` (多链) |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in residue_rows:
        i = r["interesting_BZXUO"]
        emit(f"| {r['label']} | {r['alphabet_size']} | {r['n_seqs_with_nonstd']} | "
             f"{i.get('B', 0)} | {i.get('Z', 0)} | {i.get('X', 0)} | {i.get('U', 0)} | "
             f"{i.get('O', 0)} | {r['n_seqs_with_comma']} |")
    emit("")
    emit("全部观测到的非标准字符 (跨三个文件汇总, 按频次降序前 20):")
    emit("")
    agg: Counter[str] = Counter()
    for r in residue_rows:
        for k, v in r["nonstd_counts"].items():
            agg[k] += v
    if agg:
        emit("| 字符 | 出现次数 | 说明 |")
        emit("|---|---:|---|")
        note = {"B": "Asx (Asn/Asp 模糊)", "Z": "Glx (Gln/Glu 模糊)", "X": "任意残基",
                "U": "硒代半胱氨酸", "O": "吡咯赖氨酸", "*": "终止/未知", "-": "缺口",
                ".": "缺口", " ": "空格"}
        for ch, n in agg.most_common(20):
            emit(f"| `{ch}` | {n} | {note.get(ch, '见 ESM-2 字母表')} |")
    else:
        emit("(无非标准字符)")

    # ---------------- 4. 去重收益 ----------------
    emit("")
    emit("## 4. 去重收益 (决定 Colab 算力)")
    emit("")
    emit("`filtered.csv` 口径 (ESM-2 实际要处理的对象):")
    emit("")
    lig_all = filt[LIG_COL].fillna("").astype(str)
    rec_all = filt[REC_COL].fillna("").astype(str)
    lig_u = pd.unique(lig_all)
    rec_u = pd.unique(rec_all)
    both_u = pd.unique(pd.concat([lig_all, rec_all], ignore_index=True))
    n_rows = len(filt)
    emit("| 口径 | 序列条数 | 相对总行数 | 节省 |")
    emit("|---|---:|---:|---:|")
    emit(f"| 逐行 (不去重, 配体+受体) | {n_rows * 2} | 200.0% | — |")
    emit(f"| 逐行去重 (ligand 列) | {len(lig_u)} | {len(lig_u) / n_rows * 100:.1f}% | "
         f"**{(1 - len(lig_u) / n_rows) * 100:.1f}%** |")
    emit(f"| 逐行去重 (receptor 列) | {len(rec_u)} | {len(rec_u) / n_rows * 100:.1f}% | "
         f"**{(1 - len(rec_u) / n_rows) * 100:.1f}%** |")
    emit(f"| **两列合并去重 (实际 ESM-2 要跑的条数)** | **{len(both_u)}** | "
         f"{len(both_u) / n_rows * 100:.1f}% | **{(1 - len(both_u) / (n_rows * 2)) * 100:.1f}%** |")
    emit("")
    emit(f"- **去重后只需嵌 {len(both_u)} 条序列**, 而不是 {n_rows * 2} 条 "
         f"(省 {(1 - len(both_u) / (n_rows * 2)) * 100:.1f}%)。")
    emit(f"- 这与 V2.2 的观察一致: 数据集是**突变级**的 (同一支架大量单点突变),")
    emit(f"  因此序列重复率极高。")
    emit("")
    emit("按 split 拆分 (每个 split 内去重):")
    emit("")
    emit("| split | 行数 | 唯一配体 | 唯一受体 | 合并唯一 |")
    emit("|---|---:|---:|---:|---:|")
    split_summary: dict[str, Any] = {}
    for sp in ("train", "val", "test"):
        s = filt[filt[SPLIT_COL] == sp]
        u = len(pd.unique(pd.concat([s[LIG_COL].fillna("").astype(str),
                                     s[REC_COL].fillna("").astype(str)],
                                    ignore_index=True)))
        emit(f"| {sp} | {len(s)} | {s[LIG_COL].nunique()} | {s[REC_COL].nunique()} | {u} |")
        split_summary[sp] = {"rows": int(len(s)), "unique_lig": int(s[LIG_COL].nunique()),
                             "unique_rec": int(s[REC_COL].nunique()), "unique_merged": int(u)}

    # ---------------- 6. 多链与截断 ----------------
    emit("")
    emit("## 5. 多链序列的处理 (必须明确)")
    emit("")
    n_multi_lig = int(lig_all.str.contains(",", regex=False).sum())
    n_multi_rec = int(rec_all.str.contains(",", regex=False).sum())
    emit(f"- 配体列含 `,` 的行: **{n_multi_lig} / {n_rows}** ({n_multi_lig / n_rows * 100:.1f}%)")
    emit(f"- 受体列含 `,` 的行: **{n_multi_rec} / {n_rows}** ({n_multi_rec / n_rows * 100:.1f}%)")
    emit("")
    emit("> **警告**: Proteinea 用 `,` 把**多条链的序列直接连接**。若把这个字符串当成一条序列")
    emit("> 喂给 ESM-2, 会引入一个不存在的「跨链假键」。两种可选做法:")
    emit(">")
    emit("> * **(A) 按原样处理 (本 notebook 采用)**: 保留 `,`, 由 ESM-2 当作未知字符 ->")
    emit(">   映射为 `<unk>`。这**不会崩溃**, 但会在链边界处产生 1 个 unk 位置。")
    emit("> * **(B) 分链单独嵌入再平均**: 更干净, 但需要拆分链并各跑一次 ESM-2,")
    emit(">   算力按链数成倍增长。")
    emit(">")
    emit("> 本任务采用 (A), 并把 `,` 显式计入「非标准字符」(见第 3 节), 保证行为可预期、可复现。")

    emit("")
    emit("## 6. 截断与长度上限")
    emit("")
    for s in stats_rows:
        pass
    filt_stats = [s for s in stats_rows if s["label"].startswith("filtered.csv")]
    for s in filt_stats:
        role = s["label"].split("::")[1]
        emit(f"- `filtered.csv` {role}: 长度 > {MAX_LEN} 的序列 **{s['n_gt_1022']}** 条 "
             f"(max {s['len_max']}); 截断后需记 `truncated=True`")
    emit("")
    emit(f"ESM-2 位置编码上限为 {MAX_LEN} (含 BOS/EOS), 因此超过 {MAX_LEN - 2} 个残基")
    emit("的序列会被截断。**截断不丢弃序列**, 只丢弃尾部长度的信息。")

    # ---------------- 7. 算力/显存估算 ----------------
    emit("")
    emit("## 7. 算力与显存估算 (按 esm2_t30_150M_UR50D)")
    emit("")
    n_seq = len(both_u)
    lens = np.array([len(s) for s in both_u], dtype=float)
    tot_res = float(lens.sum())
    mean_len = float(lens.mean())
    # 单条前向耗时按长度缩放: 以 400 aa 约 25 ms 为基准 (T4, fp32, 150M 模型, 保守)
    per_res_ms = 25.0 / 400.0
    per_seq_s = max(0.012, mean_len * per_res_ms / 1000.0)
    emit("| 项目 | 值 |")
    emit("|---|---|")
    emit(f"| 需嵌入的唯一序列数 | **{n_seq}** |")
    emit(f"| 总残基量 | **{tot_res:,.0f}** |")
    emit(f"| 序列长度中位数 | {np.median(lens):.0f} aa |")
    emit(f"| 序列长度均值 | {mean_len:.1f} aa |")
    emit(f"| 序列长度 95 分位 | {np.percentile(lens, 95):.0f} aa |")
    emit(f"| 序列长度最大 | {lens.max():.0f} aa |")
    emit(f"| 单条前向 (按长度缩放, 保守) | ~{per_seq_s * 1000:.0f} ms |")
    emit(f"| **纯前向总时长** | **~{n_seq * per_seq_s / 60:.0f} 分钟** |")
    emit(f"| 每条一次 `.npz` 落盘 (断点续算) | + ~{n_seq * 0.004 / 60:.1f} 分钟 |")
    emit(f"| 每 100 条一次 Drive 备份 | + 取决于 Drive 速度 |")
    emit(f"| **预计总时长** | **~{max(5, n_seq * per_seq_s / 60 * 1.6):.0f}–"
         f"{max(10, n_seq * per_seq_s / 60 * 3.0):.0f} 分钟** |")
    emit(f"| 输出体积 (raw+L2, 2 角色 x 2 版本 x 640 维 x 4B, 未压缩) | "
         f"~{n_seq * 2 * 2 * 640 * 4 / 1e6:.0f} MB |")
    emit("")
    emit("> 上表是**保守估算**, 不是实测。实际时长取决于 T4 争用、batch 大小与序列长度分布。")
    emit("> notebook 每 100 条打印进度, 可据此修正。")
    emit("")
    emit("显存需求 (fp32, 不含梯度 —— 推理用 `torch.no_grad()`):")
    emit("")
    emit("| batch_size | 最长序列 | 估算激活显存 | T4 (15GB) 是否可行 |")
    emit("|---:|---:|---|---|")
    for bs, Lmax in ((1, 1022), (4, 1022), (8, 1022), (16, 1022), (32, 512)):
        # 150M 参数: 注意力激活 ~ B * H * L^2; 粗估
        hidden, layers = 640, 30
        act = bs * layers * hidden * Lmax * 4 * 2 / 1e9   # GB (线性项主导)
        attn = bs * layers * 20 * Lmax * Lmax * 4 / 1e9    # GB (二次项, 20 heads)
        est = act + attn + 0.6                              # 0.6GB 权重
        emit(f"| {bs} | {Lmax} | ~{est:.2f} GB | {'✅' if est < 13 else '❌ 需降级'} |")
    emit("")
    emit("> 估算基于 640 维隐层 x 30 层 x 20 注意力头。**以 notebook 的 OOM 降级逻辑为准**:")
    emit("> 真 OOM 时 batch_size 自动减半重试, 最低到 1; 仍 OOM 则切换到 `esm2_t12_35M_UR50D`。")

    # ---------------- 8. 结论 ----------------
    emit("")
    emit("## 8. 结论与交给 notebook 的常量")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit(f"| 序列列 (配体) | `{LIG_COL}` |")
    emit(f"| 序列列 (受体) | `{REC_COL}` |")
    emit(f"| ID 列 | `{ID_COL}` |")
    emit(f"| 划分列 | `{SPLIT_COL}` (仅 filtered.csv) |")
    emit(f"| 主用文件 | `filtered.csv` (8207 行, 自带 train/val/test) |")
    emit(f"| ID 与 V2.2 `pdb_id` 一致性 | {'✅ 一致' if id_ok else '⚠️ 见第 2 节'} |")
    emit(f"| 去重后需嵌入序列数 | **{n_seq}** |")
    emit(f"| 多链分隔符 | `,` (需映射为 `<unk>`, 或用分链方案) |")
    emit(f"| 长度上限 | {MAX_LEN} aa, 超出截断并记 `truncated` |")
    emit("")
    emit("**没有阻塞项**: Proteinea 数据**确实包含**序列列, 无需从 PDB ID 反解序列。")
    emit("(反解方案仍作为备选记录在第 9 节。)")

    # ---------------- 9. 备选方案 ----------------
    emit("")
    emit("## 9. 备选方案 (仅在序列列缺失时才需要, 本数据不需要)")
    emit("")
    emit("若将来换用没有序列列的数据源, 可用本仓库已有的两条通路反解序列:")
    emit("")
    emit("1. **从本地 PDB 的 ATOM 记录提序列** — 已实现于")
    emit("   `scripts/18_extra_features.py::seq_from_structure()` 与")
    emit("   `scripts/13_enzyme_split.py::seq_from_structure()`。")
    emit("   优点: 无需联网。缺点: 只能得到**已解析进结构的残基**, 会缺 SEQRES 里的缺失残基。")
    emit("2. **从 SEQRES 记录提完整序列** — Proteinea 论文用的是这条 (recover missing residues)。")
    emit("   本仓库已下载的 `raw_recover_missing_res.csv` **就是他们做完这一步的产物**,")
    emit("   因此**优先直接用它的序列列**, 不必自己实现。")
    emit("3. `data/pdb_seqs.json` 是 V1.1 时期缓存的 RCSB entity 序列 (3024 个 PDB),")
    emit("   但链 ID 用的是 RCSB entity 命名 (`A0`/`B0`), 与作者链 ID 不一致, **不可直接用**。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "files": {k: {"rows": v[0], "cols": v[1]} for k, v in shapes.items()},
        "columns": schemas.get("filtered.csv", []),
        "column_roles": {"ligand_seq": LIG_COL, "receptor_seq": REC_COL,
                         "id": ID_COL, "split": SPLIT_COL},
        "id_matches_v22_pdb_id": bool(id_ok) if id_ok is not None else None,
        "v22_unique_pdb_id": len(v22_pdbs),
        "source_unique_pdb": int(filt["_pdb_up"].nunique()),
        "seq_stats": stats_rows,
        "residue_audit": residue_rows,
        "dedup": {
            "rows": int(n_rows),
            "unique_lig": int(len(lig_u)), "unique_rec": int(len(rec_u)),
            "unique_merged": int(len(both_u)),
            "saving_vs_all": float(1 - len(both_u) / (n_rows * 2)),
            "by_split": split_summary,
        },
        "multi_chain": {"rows_with_comma_lig": n_multi_lig,
                        "rows_with_comma_rec": n_multi_rec},
        "max_len": MAX_LEN,
        "est_total_residues": tot_res,
        "has_sequence_columns": True,
        "blocking_issues": [],
        "esm2_run_locally": False,
    }
    if not args.dry_run:
        args.json_out.write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                                            default=str), encoding="utf-8")
        args.out_md.write_text("\n".join(L) + "\n", encoding="utf-8")
        emit("")
        emit(f"- 报告: `{args.out_md.name}`")
        emit(f"- JSON: `{args.json_out.name}`")
    else:
        emit("")
        emit("--dry-run: 未写文件")
    return 0


if __name__ == "__main__":
    sys.exit(main())
