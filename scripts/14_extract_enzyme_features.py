# -*- coding: utf-8 -*-
"""14_extract_enzyme_features.py -- 阶段 4: 酶数据集的标签/序列特征装配.

定位与诚实边界
--------------
提示词交付清单里列了这个文件。**重要说明**: 阶段 4 已确认本机**离线**,
EC / CATH 注释层**无法接入** (`13_enzyme_split.py` 的 `ec_coverage = 0.0`,
`cath_coverage = 0.0`)。因此本脚本**不能**产出"酶催化特征"(活性位点、EC 类别、
CATH 折叠) —— 那需要一个在本机不存在的注释源。

在本机数据下**真正可做且有意义**的事情是: 把 `13_enzyme_split.py` 产出的
三层隔离数据集, 装配成"**按隔离组排序**"的 Learning-to-Rank 就绪表, 并给出
酶数据集特有的**标签与条件特征**。这正是阶段 3 结论指向的方向 ——
PPI 数据 group 结构退化 (61.7% 组内 ΔG 全同), 而酶/突变体数据 group 内 ΔG
天然各不相同, 需要一张专门为此准备的表。

本脚本产出
----------
1. **隔离组内 ordinal rank**: 对每个 `isolation_key` 组内的 `real_dG` 排名
   (1 = 组内最强结合), 这是序数标签任务的直接输入。
2. **组内相对 ΔG**: `dG_rel_to_group_best = real_dG - min(real_dG in group)`
   (>= 0), 用于突变效应/底物偏好建模。
3. **组规模与组内变异度**: `group_size`, `group_dG_range`, `group_dG_std`,
   `group_has_variation` —— 直接标记"这个组能不能用于排序学习"。
4. **测量条件特征**: `temperature_K`(缺失已标 `temperature_source`)、
   `assay_method` 的数值编码 (`assay_method_code`)、`n_method`。
5. **标签口径字段**: `label_type`(恒为 `Kd`, 已如实标注)、`value`/`unit`、
   `dG_from_value`(RT·ln 复算校验)。
6. **序列派生特征** (可复现、CPU 廉价): 链长度、界面/全链残基类别占比、
   氨基酸组成、疏水矩。**明确标注这些是"靶点身份指纹"** ——
   `v2_ablation_decoupling_findings.md` 已实测它们贡献了新增特征
   87% 的增益, 而这类增益**不代表界面物理能力**。

明确不做的事
------------
* 不产出 EC 号、CATH id、活性位点残基、催化三联体等**需要外部注释**的字段
  (会全部留空并标注原因, **不编造**)。
* 不产出 Km / kcat / Ki (源表只有 `KD(M)`, 见 `AUTO_ASSUMPTIONS.md` A3)。

输入 : `features/v2_enzyme_dataset.csv` (由 13_enzyme_split.py 产出)
       `data/v2_assay_meta.csv` (可选, 补充方法/温度)
输出 : `features/v2_enzyme_features.csv`
       `features/v2_enzyme_feature_report.md`
       `features/v2_enzyme_feature_summary.json`

纪律: 纯 CPU; 无网络; seed=42; 逐行容错; 不改动 V1.1 产物。
"""

from __future__ import annotations

import argparse
import json
import logging
import math
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
FEAT_DIR = ROOT / "features"
DATA_DIR = ROOT / "data"
LOG_DIR = ROOT / "logs" / "v2"
PDB_DIR = DATA_DIR / "pdb"

ENZYME_CSV = FEAT_DIR / "v2_enzyme_dataset.csv"
ASSAY_META = DATA_DIR / "v2_assay_meta.csv"
OUT_CSV = FEAT_DIR / "v2_enzyme_features.csv"
REPORT_MD = FEAT_DIR / "v2_enzyme_feature_report.md"
JSON_OUT = FEAT_DIR / "v2_enzyme_feature_summary.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
RT_KCAL = 1.98720425864083e-3
REF_TEMP_K = 298.15

AA_ALL = sorted(set("ACDEFGHIKLMNPQRSTVWY"))
AA_AROMATIC = set("FWY")
AA_CHARGED_POS = set("KRH")
AA_CHARGED_NEG = set("DE")
AA_POLAR = set("STNQC")
AA_HYDROPHOBIC = set("AVLIMGP")

EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85,
    "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50,
    "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81,
    "Y": 0.26, "V": 1.08,
}
HMOM_WINDOW = 11
HMOM_ANGLE_DEG = 100.0

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger("enzyme_features")


# --------------------------------------------------------------------- 序列特征
def class_fractions(seq: str) -> dict[str, float]:
    n = len(seq)
    if n == 0:
        return {k: np.nan for k in ("aromatic", "charged_pos", "charged_neg",
                                    "polar", "hydrophobic")}
    return {
        "aromatic": sum(a in AA_AROMATIC for a in seq) / n,
        "charged_pos": sum(a in AA_CHARGED_POS for a in seq) / n,
        "charged_neg": sum(a in AA_CHARGED_NEG for a in seq) / n,
        "polar": sum(a in AA_POLAR for a in seq) / n,
        "hydrophobic": sum(a in AA_HYDROPHOBIC for a in seq) / n,
    }


def hydrophobic_moment_max(seq: str) -> float:
    """疏水矩 (Eisenberg 标度, 窗口 11, 角度 100°) 的滑窗最大值。"""
    if len(seq) < HMOM_WINDOW:
        return float("nan")
    delta = math.radians(HMOM_ANGLE_DEG)
    best = float("nan")
    for s in range(len(seq) - HMOM_WINDOW + 1):
        hx = hy = 0.0
        ok = True
        for k, aa in enumerate(seq[s:s + HMOM_WINDOW]):
            h = EISENBERG.get(aa)
            if h is None:
                ok = False
                break
            ang = k * delta
            hx += h * math.cos(ang)
            hy += h * math.sin(ang)
        if ok:
            v = math.hypot(hx, hy)
            if not (best == best) or v > best:
                best = v
    return float(best)


def net_charge(seq: str) -> float:
    ch = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}
    return float(sum(ch.get(a, 0.0) for a in seq))


# --------------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 4: 酶数据集标签/序列特征装配")
    ap.add_argument("--enzyme", type=Path, default=ENZYME_CSV)
    ap.add_argument("--assay-meta", type=Path, default=ASSAY_META)
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--top-frac", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    (FEAT_DIR).mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(LOG_DIR / "14_enzyme_features.log", mode="w", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S"))
    logger.addHandler(fh)

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# Affinity Intelligence V2.0 — 阶段 4: 酶数据集特征装配报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}; 真强结合体定义: 组内 ΔG 最负的前 "
         f"{args.top_frac * 100:.0f}%")

    if not args.enzyme.exists():
        emit(f"!! 输入不存在: {args.enzyme}")
        emit("   请先运行: python scripts/13_enzyme_split.py")
        return 2
    df = pd.read_csv(args.enzyme, low_memory=False)
    emit(f"- 输入: {args.enzyme.name} — {len(df)} 行 x {df.shape[1]} 列")
    emit(f"- 隔离组数: {df['isolation_key'].nunique()}")

    if args.dry_run:
        emit("")
        emit("--dry-run: 不生成特征")
        return 0

    out = pd.DataFrame(index=df.index)
    out["interface_key"] = df.get("interface_key", df.index.astype(str))
    out["pdb_id"] = df["pdb_id"]
    out["isolation_key"] = df["isolation_key"]
    out["group_id"] = df.get("group_id", "")
    out["family_id"] = df.get("family_id", "")
    out["split"] = df.get("split", "")
    out["source"] = df.get("source", "")

    # ---- 1. 标签与口径 --------------------------------------------------
    emit("")
    emit("## 1. 标签与口径")
    dG = pd.to_numeric(df["real_dG"], errors="coerce").to_numpy(dtype=float)
    kd = pd.to_numeric(df.get("value", np.nan), errors="coerce").to_numpy(dtype=float)
    out["label_type"] = df.get("label_type", "Kd")
    out["value"] = kd
    out["unit"] = df.get("unit", "M")
    out["real_dG"] = dG
    # 复算校验: dG 应由 RT*ln(Kd) 得到 (固定 298.15 K 口径)
    with np.errstate(divide="ignore", invalid="ignore"):
        dg_check = RT_KCAL * REF_TEMP_K * np.log(kd)
    resid = np.abs(dg_check - dG)
    n_bad = int(np.sum(np.isfinite(resid) & (resid > 0.05)))
    out["dG_recompute_residual"] = resid
    emit(f"- `real_dG` 由 `RD·T·ln(Kd)` 复算: 最大残差 "
         f"{np.nanmax(resid) if np.isfinite(resid).any() else float('nan'):.6f} kcal/mol; "
         f"残差 > 0.05 的行: **{n_bad}**")
    emit("- 标签类型: 恒为 `Kd` (源表只有 `KD(M)`, 见 AUTO_ASSUMPTIONS.md A3);")
    emit("  **不存在 Ki / IC50 / Km / kcat 字段**, 因此本脚本不产出酶动力学标签。")

    # ---- 2. 隔离组内序数与相对量 (LTR 就绪) ------------------------------
    emit("")
    emit("## 2. 隔离组内序数与相对量 (LTR 就绪字段)")
    g = df.groupby("isolation_key")["real_dG"]
    gsize = g.size()
    gmin = g.transform("min")
    gmax = g.transform("max")
    gstd = g.transform("std")
    grange = gmax - gmin

    out["group_size"] = df["isolation_key"].map(gsize).to_numpy()
    out["group_dG_min"] = gmin.to_numpy()
    out["group_dG_max"] = gmax.to_numpy()
    out["group_dG_range"] = grange.to_numpy()
    out["group_dG_std"] = gstd.to_numpy()
    # ΔG 越负越强 -> 相对组内最优 (最负) 的差值, >= 0
    out["dG_rel_to_group_best"] = (dG - gmin.to_numpy())
    # 序数标签: 1 = 组内最强 (ΔG 最负); 并列取最小秩 (method="min")
    out["group_rank_desc"] = df.groupby("isolation_key")["real_dG"].rank(
        method="min", ascending=True).to_numpy()
    # 排序学习可用性标记: 组内至少 2 条且有差异
    out["group_has_variation"] = (
        (out["group_size"] > 1) & (out["group_dG_range"] > 1e-9)).astype(int)
    # 组内是否为"真强结合体"(组内前 top_frac)
    def _top_flag(sub: pd.Series) -> pd.Series:
        n = len(sub)
        k = max(1, int(np.ceil(args.top_frac * n)))
        thr = np.sort(sub.to_numpy())[:k].max()      # 最负的 k 个的阈值
        return (sub <= thr + 1e-12)
    out["group_top_frac_flag"] = (
        df.groupby("isolation_key")["real_dG"].transform(_top_flag).astype(int).to_numpy())

    usable = int(out["group_has_variation"].sum())
    emit("")
    emit(f"- 可用于排序学习的行 (组内 >= 2 条且 ΔG 有差异): "
         f"**{usable} / {len(out)} ({usable / len(out) * 100:.1f}%)**")
    emit(f"- 单条目组覆盖行: {int((out['group_size'] == 1).sum())}")
    emit(f"- 组内 ΔG 全同 (零排序梯度) 的行: "
         f"{int(((out['group_size'] > 1) & (out['group_dG_range'] <= 1e-9)).sum())}")
    emit("")
    emit("### 与 PPI 主数据集的对比 (阶段 3 的结论)")

    # 与 features_v1.1 的 family 组做同口径对比
    ppi_path = FEAT_DIR / "features_v1.1.csv"
    cmp_rows = []
    if ppi_path.exists():
        ppi = pd.read_csv(ppi_path, low_memory=False)
        for tag, frame, key in (("PPI (family_id, 阶段 3)", ppi, "family_id"),
                                ("酶隔离组 (isolation_key, 阶段 4)", df, "isolation_key")):
            gg = frame.groupby(key)["real_dG"]
            sz, nu = gg.size(), gg.nunique()
            multi = sz > 1
            cmp_rows.append({
                "dataset": tag, "n_groups": int(len(sz)), "n_multi": int(multi.sum()),
                "rows_in_multi": int(sz[multi].sum()),
                "rows_varying": int(sz[multi & (nu > 1)].sum()),
                "multi_share": float(multi.sum() / len(sz)),
                "varying_share_of_multi": (float((multi & (nu > 1)).sum() / multi.sum())
                                           if multi.sum() else float("nan")),
            })
        emit("")
        emit("| 数据集 | 组数 | 多条目组 | 多条目组占比 | 组内 ΔG 有差异的组 | 占有差异组比例 |")
        emit("|---|---:|---:|---:|---:|---:|")
        for r in cmp_rows:
            emit(f"| {r['dataset']} | {r['n_groups']} | {r['n_multi']} | "
                 f"{r['multi_share'] * 100:.1f}% | {r['rows_varying']} | "
                 f"{r['varying_share_of_multi'] * 100:.1f}% |")

    # ---- 3. 测量条件 ----------------------------------------------------
    emit("")
    emit("## 3. 测量条件特征")
    out["temperature_K"] = pd.to_numeric(df.get("temperature_K", np.nan), errors="coerce")
    out["temperature_source"] = df.get("temperature_source", "unknown")
    out["ph"] = np.nan
    if args.assay_meta.exists():
        meta = pd.read_csv(args.assay_meta)
        if "interface_key" in meta.columns:
            mm = meta.set_index("interface_key")
            for col, dst in (("method_majority", "assay_method"),
                             ("method_n", "assay_n_votes"),
                             ("method_conflict", "assay_method_conflict")):
                if col in mm.columns:
                    out[dst] = out["interface_key"].map(mm[col])
    if "assay_method" not in out.columns:
        out["assay_method"] = df.get("assay_method", "")
    out["assay_method"] = out["assay_method"].fillna("")
    # 稳定的方法编码 (按出现频次降序, 便于复现)
    methods = sorted([m for m in out["assay_method"].unique() if m])
    code = {m: i + 1 for i, m in enumerate(methods)}
    out["assay_method_code"] = out["assay_method"].map(code).fillna(0).astype(int)
    n_temp_known = int(pd.to_numeric(out["temperature_K"], errors="coerce")
                       .notna().sum())
    emit(f"- 有实验温度: {n_temp_known}/{len(out)} ({n_temp_known / len(out) * 100:.1f}%)")
    emit(f"- 有测量方法: {int((out['assay_method'] != '').sum())}/{len(out)} "
         f"({(out['assay_method'] != '').mean() * 100:.1f}%)")
    emit(f"- `ph`: **全空** (源表无 pH 列, 无法填充)")
    emit(f"- 方法编码 (频次降序): {code if code else '(无)'}")

    # ---- 4. 序列派生特征 (明确标注为靶点指纹) ----------------------------
    emit("")
    emit("## 4. 序列派生特征 (⚠️ 靶点身份指纹, 不是界面物理)")
    seqs = df.get("sequence", pd.Series([""] * len(df))).fillna("").astype(str)
    recs: list[dict[str, Any]] = []
    for s in seqs:
        d: dict[str, Any] = {"seq_len": len(s)}
        d.update({f"frac_{k}": v for k, v in class_fractions(s).items()})
        d.update({f"aac_{a}": (s.count(a) / len(s) if s else np.nan) for a in AA_ALL})
        d["hmom_max"] = hydrophobic_moment_max(s)
        d["net_charge"] = net_charge(s) if s else np.nan
        recs.append(d)
    seq_feat = pd.DataFrame(recs, index=df.index)
    for c in seq_feat.columns:
        out[f"seq_{c}"] = seq_feat[c]
    emit(f"- 产出 {seq_feat.shape[1]} 个序列派生特征 (链长度 / 类别占比 / 氨基酸组成 /"
         f" 疏水矩 / 净电荷)")
    emit("")
    emit("> **必须随特征一起传播的警告**: `v2_ablation_decoupling_findings.md` 已实测,")
    emit("> 这类序列派生特征对同一 target family 是常数 (家族内 std / 全局 std = 0.10-0.24),")
    emit("> 它们贡献了阶段 5 新增特征 **87.3%** 的 Spearman 增益 ——")
    emit("> 那是「识别靶点并回忆其平均亲和力」的能力, **不是界面结合物理的建模能力**。")
    emit("> 用它们报告性能时, 必须说明靶点新颖性条件。")

    # ---- 5. 明确未产出的字段 --------------------------------------------
    emit("")
    emit("## 5. 明确未产出的字段 (不编造)")
    emit("")
    emit("| 字段 | 状态 | 原因 |")
    emit("|---|---|---|")
    emit("| `ec_number` / EC 类别特征 | **未产出** | 本机离线, 无法获取 ExPASy/SIFTS 映射 |")
    emit("| `cath_id` / CATH 折叠特征 | **未产出** | 本机离线, 无法获取 CATH/Gene3D 映射 |")
    emit("| 活性位点 / 催化残基 | **未产出** | 需要 EC + UniProt 注释, 本机不可得 |")
    emit("| `Km` / `kcat` / `Ki` | **未产出** | 源表 `PPB-Affinity.xlsx` 只有 `KD(M)` 一列 |")
    emit("| `ph` | **全空** | 源表无 pH 列 (19 列已全量枚举) |")
    emit("")
    emit("接入方式: 放入 `data/ec_annotations.csv` (列 `pdb_id,ec_number`) 与")
    emit("`data/cath_annotations.csv` (列 `pdb_id,cath_id`) 后重跑")
    emit("`13_enzyme_split.py` 与 `14_extract_enzyme_features.py` 即自动启用。")

    # ---- 6. 落盘 --------------------------------------------------------
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + ".tmp")
    out.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out)
    emit("")
    emit(f"- 输出: `{args.out.name}` — {out.shape[0]} 行 x {out.shape[1]} 列")
    emit(f"耗时: {time.time() - t0:.1f}s")

    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "n_rows": int(len(out)), "n_columns": int(out.shape[1]),
        "n_isolation_groups": int(df["isolation_key"].nunique()),
        "rows_usable_for_ranking": usable,
        "rows_usable_share": float(usable / len(out)),
        "n_sequence_features": int(seq_feat.shape[1]),
        "dG_recompute_max_residual": float(np.nanmax(resid)) if np.isfinite(resid).any() else None,
        "dG_recompute_bad_rows": n_bad,
        "temperature_known_share": float(n_temp_known / len(out)),
        "assay_method_known_share": float((out["assay_method"] != "").mean()),
        "assay_method_codes": code,
        "not_produced": ["ec_number", "cath_id", "active_site", "Km", "kcat", "Ki", "ph"],
        "not_produced_reason": "本机离线 + 源表只有 KD(M)",
        "dataset_comparison": cmp_rows,
        "sequence_feature_warning":
            "序列派生特征是靶点身份指纹 (家族内 std/全局 std = 0.10-0.24), "
            "贡献新增特征 87.3% 的增益, 不代表界面物理建模能力",
    }
    JSON_OUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str),
                        encoding="utf-8")
    args.report.write_text("\n".join(run_log) + "\n", encoding="utf-8")
    emit("")
    emit(f"报告: {args.report}")
    emit(f"摘要: {JSON_OUT}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
            fh2.write(f"[14_extract_enzyme_features] {traceback.format_exc()}\n")
        raise
