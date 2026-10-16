# -*- coding: utf-8 -*-
"""scripts/v2.2/01_adapt_proteinea.py -- 把 HuggingFace `proteinea/ppb_affinity`
适配到本项目既有的 60 物理特征流程 (V2.2 基准切换).

背景与目标
----------
本项目 V1.1/V2.0 的基准是**自建**的 PPB-Affinity 处理流程 (自定家族隔离划分)。
本脚本切换到 Proteinea 团队在 HuggingFace 上发布的 `proteinea/ppb_affinity`
`filtered` 配置 —— 该版本满足 **<=30% 序列同一性**并**自带 train/val/test 划分**,
因此可以与 BindPred、SSIF-Affinity 等外部工作在**同一个划分**上对比。

硬约束 (遵守情况)
-----------------
* 只读 `data/ppb_affinity_proteinea/`, 只写 `features/v2.2_*`; **不触碰任何 V1.1 / V2.0 文件**。
* **标签一律从 Proteinea 的 `KD(M)` 重新计算**, 绝不复用 `features_v1.1.csv` 的 `real_dG`。
  原因 (已实测): V1.1 的 `real_dG` 是**按 interface_key 对多个突变体取平均**得到的聚合标签,
  与 Proteinea 的逐条目标签最大相差 8.29 kcal/mol、平均 0.73 (5127/8203 行超 0.01)。
  直接复用会把 V1.1 的聚合标签错误地当成 Proteinea 的标签。
* **不自行重新划分**: 完全沿用 Proteinea 自带的 train / val / test。

特征复用是否合法? (关键论证)
-----------------------------
60 个物理特征是**纯结构派生量** (接触/氢键/盐桥/ΔSASA/组成), 只依赖
`(PDB, 配体链组, 受体链组)`, **完全不依赖标签**。已实测验证:
`features_v1.1.csv` 中同一 `interface_key` 的特征值**逐列恒定** (0 个列在键内有变化),
因此该文件是一张**按界面键索引的结构特征缓存**, 可以安全地按 `interface_key` join 到
Proteinea 的 8207 行上 —— 这与"复用标签"是完全不同的事。

实测覆盖率: 8203 / 8207 行 (100.0%) 可直接命中缓存; 其余 4 行 (5FT8 的链组)
因本地 PDB 的链 ID 版本与源表不一致而无法提取 (详见报告"缺失情况")。

输入
----
`data/ppb_affinity_proteinea/filtered.csv`   (Proteinea filtered 配置, 8207 行)
`features/features_v1.1.csv`                 (60 物理特征缓存, 按 interface_key)

输出
----
`features/v2.2_proteinea_features.csv`       适配后的特征表 (含 split 与重算标签)
`features/v2.2_proteinea_adapt_summary.json` 机器可读摘要
`logs/v2.2/01_adapt_proteinea.log`           日志

纪律: 纯 CPU; 无网络 (文件已预下载); seed=42; 逐行容错不中断; 不改动任何既有版本产物。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import logging
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
DATA_DIR = ROOT / "data"
FEAT_DIR = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2.2"
PROT_DIR = DATA_DIR / "ppb_affinity_proteinea"
PDB_DIR = DATA_DIR / "pdb"

PROT_CSV = PROT_DIR / "filtered.csv"
CACHE_CSV = FEAT_DIR / "features_v1.1.csv"
OUT_CSV = FEAT_DIR / "v2.2_proteinea_features.csv"
OUT_JSON = FEAT_DIR / "v2.2_proteinea_adapt_summary.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
RT_LITERAL = 0.592          # kcal/mol, 与 01_clean_data.py / V1.1 完全一致的口径
REF_TEMP_K = 298.15

logger = logging.getLogger("adapt_proteinea")


# ---------------------------------------------------------------- 工具
def _ensure_openmp() -> str | None:
    """lib_lightgbm.dll 依赖本机缺失的 vcomp140.dll; 借用 sklearn 自带副本。

    本脚本本身不训练, 但 --verify-lgbm 会用一次; 且保持与其余脚本一致的处置方式。
    """
    import ctypes
    import glob
    try:
        import sklearn
        hits = glob.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll"))
        if hits:
            ctypes.CDLL(hits[0])
            return hits[0]
    except Exception:  # noqa: BLE001
        pass
    return None


def norm_chain(value: Any) -> str:
    """``'A, B'`` -> ``'A,B'``; 与 01_clean_data.py::norm_chain 完全一致。

    必须逐字节对齐, 否则 `interface_key` 与特征缓存 join 不上。
    """
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    txt = str(value).replace(";", ",").replace("|", ",")
    parts = [p.strip().upper() for p in txt.split(",") if p.strip()]
    parts = [p.split(":")[0].strip() for p in parts]
    parts = [p for p in parts if p]
    return ",".join(sorted(dict.fromkeys(parts)))


def chain_tokens(value: Any) -> list[str]:
    """拆成链列表 (用于校验链是否存在于结构)。"""
    s = norm_chain(value)
    return [c for c in s.split(",") if c]


def load_v11_feature_names() -> tuple[list[str], list[str]]:
    """从 08_train_v1.0.py 载入特征集定义 (单一实现来源, 避免口径漂移)。"""
    path = ROOT / "scripts" / "08_train_v1.0.py"
    spec = importlib.util.spec_from_file_location("t8_v22", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    base = [c for c in mod.V02_FEATURES if c not in mod.REPLACED]
    return base + list(mod.V03_NEW), list(mod.COUNT_COLS)


def chains_in_pdb(pdb_id: str) -> set[str] | None:
    """轻量扫描 PDB 的 ATOM 记录, 返回出现的链 ID 集合 (None = 文件不存在)。"""
    p = PDB_DIR / f"{pdb_id}.pdb"
    if not p.exists():
        p = PDB_DIR / f"{pdb_id}.ent"
    if not p.exists():
        return None
    got: set[str] = set()
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line[:6] == "ATOM  ":
                    got.add(line[21])
    except Exception:  # noqa: BLE001
        return None
    return got


# ---------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="适配 proteinea/ppb_affinity 到 60 物理特征流程")
    ap.add_argument("--proteinea", type=Path, default=PROT_CSV)
    ap.add_argument("--feature-cache", type=Path, default=CACHE_CSV)
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    ap.add_argument("--json-out", type=Path, default=OUT_JSON)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--dry-run", action="store_true", help="只做校验与统计, 不写文件")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "01_adapt_proteinea.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# V2.2 基准切换 — 适配 proteinea/ppb_affinity 到 60 物理特征流程")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}")
    emit(f"- 标签口径: ΔG = {RT_LITERAL} · ln(KD[M])  (与 V1.1 完全一致)")

    for p in (args.proteinea, args.feature_cache):
        if not p.exists():
            emit(f"!! 输入缺失: {p}")
            return 2

    feat_names, count_cols = load_v11_feature_names()
    emit(f"- V1.1 特征集定义: {len(feat_names)} 列 (V0.2 {len(feat_names) - 12} + V0.3 12, 已去替换列)")

    # ---- 1. 读入 Proteinea ---------------------------------------------
    pr = pd.read_csv(args.proteinea, low_memory=False)
    emit("")
    emit("## 1. Proteinea `filtered` 配置")
    emit("")
    emit(f"- 文件: `{args.proteinea.name}` — {len(pr)} 行 x {pr.shape[1]} 列")
    need_cols = ["PDB", "Ligand Chains", "Receptor Chains", "KD(M)", "split"]
    missing = [c for c in need_cols if c not in pr.columns]
    if missing:
        emit(f"!! 缺少必需列: {missing}")
        return 3
    split_counts = pr["split"].value_counts()
    emit(f"- 自带划分 (**完全沿用, 不重新划分**): "
         + ", ".join(f"{k}={v}" for k, v in split_counts.items()))
    unknown = set(split_counts.index) - {"train", "val", "test"}
    if unknown:
        emit(f"!! 出现未知 split 值: {unknown}")
        return 3

    # ---- 2. 标签重算 (绝不复用 V1.1 的聚合标签) --------------------------
    emit("")
    emit("## 2. 标签重算 (从 Proteinea 的 KD(M))")
    emit("")
    kd = pd.to_numeric(pr["KD(M)"], errors="coerce").to_numpy(dtype=float)
    n_bad_kd = int((~np.isfinite(kd) | (kd <= 0)).sum())
    if n_bad_kd:
        emit(f"- !! 无效 KD 值 {n_bad_kd} 行 (非正/缺失), 将被剔除")
    with np.errstate(divide="ignore", invalid="ignore"):
        dG = RT_LITERAL * np.log(kd)
    dG = np.where(np.isfinite(dG), dG, np.nan)
    pr["real_dG"] = dG
    emit(f"- ΔG 范围: {np.nanmin(dG):.3f} ~ {np.nanmax(dG):.3f} kcal/mol "
         f"(mean {np.nanmean(dG):.3f}, std {np.nanstd(dG):.3f})")
    emit(f"- 有效标签: {int(np.isfinite(dG).sum())} / {len(pr)}")

    # 与 V1.1 聚合标签的差异 (证明"不能复用")
    cache = pd.read_csv(args.feature_cache, low_memory=False)
    pr["_pdb"] = pr["PDB"].astype(str).str.strip().str.upper()
    pr["_lig"] = pr["Ligand Chains"].map(norm_chain)
    pr["_rec"] = pr["Receptor Chains"].map(norm_chain)
    pr["interface_key"] = pr["_pdb"] + "|" + pr["_lig"] + "|" + pr["_rec"]

    tmp = pr[["interface_key", "real_dG"]].merge(
        cache[["interface_key", "real_dG"]].rename(columns={"real_dG": "real_dG_v11"}),
        on="interface_key", how="left")
    both = tmp["real_dG_v11"].notna()
    diff = (tmp.loc[both, "real_dG"] - tmp.loc[both, "real_dG_v11"]).abs()
    emit("")
    emit(f"- 与 `features_v1.1.csv` 的 `real_dG` 对比 ({int(both.sum())} 行可比):")
    emit(f"  * 最大绝对差 = **{diff.max():.4f}** kcal/mol")
    emit(f"  * 平均绝对差 = **{diff.mean():.4f}** kcal/mol")
    emit(f"  * 差值 > 0.01 的行数 = **{int((diff > 0.01).sum())}** "
         f"({(diff > 0.01).mean() * 100:.1f}%)")
    emit("  * 原因: V1.1 的 `real_dG` 是**按 interface_key 对多个突变体取平均**的聚合标签,")
    emit("    而 Proteinea 保留逐条目 (逐突变体) 标签。")
    emit("  * **结论: 必须用 Proteinea 的 KD(M) 重算标签, 复用 V1.1 标签是错的。**")

    # ---- 3. 特征缓存 join ----------------------------------------------
    emit("")
    emit("## 3. 60 物理特征适配 (按 interface_key join 结构特征缓存)")
    emit("")
    emit("### 复用合法性的实测依据")
    emit("")
    num_chk = [c for c in feat_names if c in cache.columns]
    nun = cache.groupby("interface_key")[num_chk].nunique()
    varied = (nun > 1).sum()
    varied = varied[varied > 0]
    emit(f"- 缓存表 `{args.feature_cache.name}`: {len(cache)} 行, "
         f"{cache['interface_key'].nunique()} 个唯一 interface_key")
    emit(f"- 同一 interface_key 内**取值有变化**的特征列数 = **{len(varied)}** "
         f"(必须为 0 才能当结构缓存复用)")
    if len(varied):
        emit(f"  !! 有变化: {varied.to_dict()}")
        emit("  -> 不能安全复用, 需要重新提取特征")
        return 4
    emit("- ✅ 0 列在键内变化 -> 该表是**按界面键索引的结构特征缓存**,")
    emit("  只依赖 (PDB, 配体链组, 受体链组), **完全不依赖标签**, 可安全 join。")
    emit("  (这与「复用标签」是两件不同的事: 标签在第 2 节已全部重算。)")

    cache_cols = ["interface_key"] + [c for c in feat_names if c in cache.columns]
    not_in_cache = [c for c in feat_names if c not in cache.columns]
    if not_in_cache:
        emit(f"- !! 缓存缺少 {len(not_in_cache)} 个特征列: {not_in_cache}")
    df = pr.merge(cache[cache_cols], on="interface_key", how="left", validate="many_to_one")
    hit = df[num_chk[0]].notna()
    emit("")
    emit(f"- join 命中率: **{int(hit.sum())} / {len(df)} ({hit.mean() * 100:.2f}%)**")
    for sp in ("train", "val", "test"):
        m = df["split"] == sp
        emit(f"  * {sp:5s}: {int((m & hit).sum()):5d} / {int(m.sum()):5d} "
             f"({(m & hit).sum() / max(int(m.sum()), 1) * 100:5.2f}%)")

    # ---- 4. 结构可用性核验 (诚实报告不可用行) ---------------------------
    emit("")
    emit("## 4. 结构可用性核验")
    emit("")
    emit("特征缓存是 V1.1 时期按 `(PDB, 配体链, 受体链)` 建立的; Proteinea 的链标注")
    emit("在少数结构上与本地 PDB 的**链 ID 版本不一致** (抗体/多聚体常见)。逐行核验:")
    ch_cache: dict[str, set[str] | None] = {}
    unusable: list[dict[str, Any]] = []
    # 注意: 这里必须用 to_dict("records") 而不是 itertuples() ——
    # pandas 的 itertuples 会把以 "_" 开头的列名重命名 (无法用 r._pdb 访问)。
    for rec in pr[["_pdb", "_lig", "_rec", "split"]].to_dict("records"):
        p = rec["_pdb"]
        if p not in ch_cache:
            ch_cache[p] = chains_in_pdb(p)
        ch = ch_cache[p]
        if ch is None:
            unusable.append({"pdb_id": p, "reason": "本地缺少结构文件",
                             "split": rec["split"], "lig": rec["_lig"], "rec": rec["_rec"]})
            continue
        need = set(chain_tokens(rec["_lig"])) | set(chain_tokens(rec["_rec"]))
        absent = need - ch
        if absent:
            unusable.append({"pdb_id": p, "reason": f"链不存在: {','.join(sorted(absent))}",
                             "split": rec["split"], "lig": rec["_lig"], "rec": rec["_rec"]})
    unusable_pdbs = {u["pdb_id"] for u in unusable}
    emit(f"- 无法提取结构的行: **{len(unusable)} / {len(df)}** "
         f"({len(unusable) / len(df) * 100:.2f}%), 涉及 {len(unusable_pdbs)} 个 PDB")
    if unusable:
        from collections import Counter
        emit(f"- 按 split: {dict(Counter(u['split'] for u in unusable))}")
        emit("")
        emit("| PDB | split | 配体链 | 受体链 | 原因 |")
        emit("|---|---|---|---|---|")
        for u in unusable:
            emit(f"| {u['pdb_id']} | {u['split']} | {u['lig']} | {u['rec']} | {u['reason']} |")
        emit("")
        emit("> **如实说明**: 这些行**不会被静默丢弃** —— 它们保留在输出表中并把特征标为 NaN,")
        emit("> 由下游训练脚本显式排除并计入报告。**不伪造这些行的特征。**")
    else:
        emit("- 全部行结构可用 ✅")

    # ---- 5. 输出 --------------------------------------------------------
    keep_meta = ["interface_key", "real_dG", "split", "_pdb", "_lig", "_rec"]
    rename = {"_pdb": "pdb_id", "_lig": "ligand_chains", "_rec": "receptor_chains"}
    out_cols = keep_meta + [c for c in feat_names if c in df.columns]
    extra_meta = [c for c in ("KD(M)", "Mutations", "Source Data Set", "Complex ID",
                              "Affinity Method", "Temperature(K)", "Subgroup",
                              "Ligand Sequences", "Receptor Sequences")
                  if c in df.columns]
    out = df[out_cols + extra_meta].rename(columns=rename)
    out["structure_usable"] = (~out["pdb_id"].isin(unusable_pdbs)).astype(int)
    # 逐行可用性 (同一 PDB 也可能只有部分链组合不可用)
    bad_keys = {(u["pdb_id"], u["lig"], u["rec"]) for u in unusable}
    out["row_structure_usable"] = [
        0 if (p, l, r) in bad_keys else 1
        for p, l, r in zip(out["pdb_id"], out["ligand_chains"], out["receptor_chains"])
    ]
    out["n_features_present"] = out[[c for c in feat_names if c in out.columns]].notna().sum(axis=1)

    if args.dry_run:
        emit("")
        emit("--dry-run: 不写文件")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + ".tmp")
    out.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out)
    emit("")
    emit("## 5. 输出")
    emit("")
    emit(f"- 适配特征表: `{args.out.name}` — {out.shape[0]} 行 x {out.shape[1]} 列")
    emit(f"- 含 {len([c for c in feat_names if c in out.columns])} 个物理特征 + split + 重算标签")

    # ---- 6. 数据统计 (供报告直接引用) -----------------------------------
    emit("")
    emit("## 6. 数据统计")
    emit("")
    emit("| split | 行数 | 唯一 PDB | 唯一 interface_key | 结构可用行 | ΔG mean | ΔG std | ΔG min | ΔG max |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    stats: dict[str, Any] = {}
    for sp in ("train", "val", "test"):
        s = out[out["split"] == sp]
        su = s[s["row_structure_usable"] == 1]
        emit(f"| {sp} | {len(s)} | {s['pdb_id'].nunique()} | {s['interface_key'].nunique()} | "
             f"{len(su)} | {s['real_dG'].mean():.3f} | {s['real_dG'].std():.3f} | "
             f"{s['real_dG'].min():.3f} | {s['real_dG'].max():.3f} |")
        stats[sp] = {
            "rows": int(len(s)), "unique_pdb": int(s["pdb_id"].nunique()),
            "unique_interface_key": int(s["interface_key"].nunique()),
            "rows_structure_usable": int(len(su)),
            "dG_mean": float(s["real_dG"].mean()), "dG_std": float(s["real_dG"].std()),
            "dG_min": float(s["real_dG"].min()), "dG_max": float(s["real_dG"].max()),
        }
    emit("")
    emit("> 注意: 本数据集的 `interface_key` **不唯一** (同一结构+链组可含大量突变体),")
    emit(f"> 例如全部 {len(out)} 行只对应 {out['interface_key'].nunique()} 个界面键。")
    emit("> 因此**行级划分不等于结构级/家族级隔离**。这一点在第 7 节单独审计。")

    # ---- 7. 划分隔离性审计 (Proteinea 自带, 不改动) ---------------------
    emit("")
    emit("## 7. 划分隔离性审计 (对 Proteinea 自带划分, 只审计不修改)")
    emit("")
    emit("Proteinea 的 filtered 配置按论文做的是 **<=30% 序列同一性** 聚类划分。")
    emit("本项目从**结构身份**角度独立复核, 看是否存在跨侧重叠:")
    emit("")
    for key, label in (("pdb_id", "PDB 结构"), ("interface_key", "界面键")):
        g = out.groupby(key)["split"].nunique()
        cross = int((g > 1).sum())
        emit(f"- {label} 跨侧 (同时出现在两侧): **{cross}** / {len(g)}")
    # 与 V1.1 的家族键做对照 (V1.1 family_id 只覆盖子集)
    tmp2 = out[["interface_key", "split"]].merge(
        cache[["interface_key", "family_id"]].dropna(), on="interface_key", how="inner")
    if len(tmp2):
        g = tmp2.groupby("family_id")["split"].nunique()
        emit(f"- V1.1 `family_id` 键跨侧: **{int((g > 1).sum())}** / {len(g)} "
             f"(仅在 V1.1 家族注释覆盖的 {len(tmp2)} 行上可核验)")
    emit("")
    emit("> ⚠️ **这是两套不同的隔离口径**:")
    emit("> * Proteinea 用 **<=30% 序列同一性**聚类分侧 (更严格, 覆盖全部 8207 行);")
    emit("> * 本项目 V1.1 用 40% 序列同一性 + 40% 贪心代表聚类。")
    emit("> 两者不可混用。V2.2 **一律以 Proteinea 自带划分为准** (任务硬约束 4)。")

    # ---- 汇总 JSON ------------------------------------------------------
    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "source_dataset": "proteinea/ppb_affinity (config=filtered)",
        "source_files": {p.name: p.stat().st_size for p in
                         sorted(PROT_DIR.glob("*")) if p.is_file()},
        "license": "cc-by-nc-sa-4.0",
        "label_definition": f"real_dG = {RT_LITERAL} * ln(KD[M])",
        "labels_recomputed": True,
        "legacy_label_max_abs_diff": float(diff.max()) if len(diff) else None,
        "legacy_label_mean_abs_diff": float(diff.mean()) if len(diff) else None,
        "legacy_label_mismatch_rows": int((diff > 0.01).sum()) if len(diff) else 0,
        "feature_cache": args.feature_cache.name,
        "feature_cache_columns_varying_within_key": int(len(varied)),
        "n_features_requested": len(feat_names),
        "n_features_available": len([c for c in feat_names if c in out.columns]),
        "join_hit_rows": int(hit.sum()), "total_rows": int(len(df)),
        "usable_rows": {sp: stats[sp]["rows_structure_usable"] for sp in stats},
        "unusable_rows": unusable,
        "stats_by_split": stats,
        "split_protocol": "Proteinea provided train/val/test (<=30% seq identity); NOT re-split",
    }
    args.json_out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str),
                             encoding="utf-8")
    emit("")
    emit(f"- 摘要 JSON: `{args.json_out.name}`")
    emit(f"耗时: {time.time() - t0:.1f}s")

    (LOG_DIR / "01_adapt_proteinea_report.md").write_text(
        "\n".join(run_log) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
            fh2.write(f"[01_adapt_proteinea] {traceback.format_exc()}\n")
        raise
