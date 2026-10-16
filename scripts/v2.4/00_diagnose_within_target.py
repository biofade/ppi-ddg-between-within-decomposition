# -*- coding: utf-8 -*-
"""scripts/v2.4/00_diagnose_within_target.py -- V2.4 任务 0: within-target 可行性诊断.

为什么先做这个
--------------
V2.3.2 的错误分析给出 SKEMPI 子集 ρ=+0.5848, 但那是 **n=34** 的极小样本。
在写正式基准前必须回答三个决定性问题:

1. 数据量到底有多少? (SKEMPI 全表 vs 我们特征表里的 SKEMPI 子集)
2. 每个复合物有多少个突变体? (太少则 per-complex Spearman 无法计算)
3. **关键**: V2.2/V2.3.2 的嵌入是按 `pdb_id` join 的 (结构级平均) ——
   如果同一复合物内配体/受体序列**完全相同**, 那么同一复合物的所有突变体
   会拿到**完全相同的特征**, 此时 within-target Spearman 在数学上**必然为 NaN**
   (全并列)。必须先验证这一点, 否则整个 within-target 基准是空中楼阁。

输出: `features/v2.4_within_target_diagnosis.md`
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
DATA24 = ROOT / "data" / "v2.4"
SKEMPI_RAW = DATA24 / "skempi_v2_cleaned.csv"
V22 = FEAT / "v2.2_proteinea_features.csv"
CACHE = FEAT / "v2_esm2_cache.npz"
OUT = FEAT / "v2.4_within_target_diagnosis.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    emit("# V2.4 任务 0 — within-target 可行性诊断")
    emit("")

    # ---------- 1. SKEMPI 全表 ----------
    emit("## 1. SKEMPI 2.0 全表 (`skempi_v2_cleaned.csv`)")
    emit("")
    if SKEMPI_RAW.exists():
        sk = pd.read_csv(SKEMPI_RAW, sep=";", low_memory=False)
        # 首列名为 '#Pdb'
        pdbcol = sk.columns[0]
        emit(f"- 文件: `data/v2.4/skempi_v2_cleaned.csv` ({SKEMPI_RAW.stat().st_size} B)")
        emit(f"- 形状: **{sk.shape[0]} 行 x {sk.shape[1]} 列**")
        emit(f"- 复合物键列: `{pdbcol}` (格式 `PDB_配体链_受体链`)")
        sk["_complex"] = sk[pdbcol].astype(str)
        sk["_pdb"] = sk["_complex"].str.split("_").str[0].str.upper()
        vc = sk.groupby("_complex").size()
        emit("")
        emit("| 项 | 值 |")
        emit("|---|---:|")
        emit(f"| 复合物数 (唯一 `#Pdb`) | {sk['_complex'].nunique()} |")
        emit(f"| 唯一 PDB | {sk['_pdb'].nunique()} |")
        emit(f"| 每复合物突变体数: 中位 | {vc.median():.0f} |")
        emit(f"| 每复合物突变体数: 均值 | {vc.mean():.1f} |")
        emit(f"| 每复合物突变体数: 最大 | {vc.max()} |")
        emit(f"| 复合物数 (突变体 >= 3) | {int((vc >= 3).sum())} |")
        emit(f"| 复合物数 (突变体 >= 5) | {int((vc >= 5).sum())} |")
        emit(f"| 覆盖行数 (突变体 >= 3) | {int(vc[vc >= 3].sum())} |")
        emit(f"| 覆盖行数 (突变体 >= 5) | {int(vc[vc >= 5].sum())} |")
        emit("")
        emit(f"- `Affinity_mut_parsed` 非空: {int(sk['Affinity_mut_parsed'].notna().sum())}")
        emit(f"- `Affinity_wt_parsed` 非空: {int(sk['Affinity_wt_parsed'].notna().sum())}")
    else:
        emit("- ⚠️ SKEMPI 全表**不存在** (下载失败?), 只能用 V2.2 子集")
        sk = None

    # ---------- 2. V2.2 表里的 SKEMPI 子集 ----------
    emit("")
    emit("## 2. V2.2 特征表里的 SKEMPI 子集")
    emit("")
    d = pd.read_csv(V22, low_memory=False)
    sc = [c for c in d.columns if "Source" in c][0]
    s = d[d[sc].astype(str).str.contains("SKEMPI", case=False, na=False)].copy()
    emit(f"- 行数: **{len(s)}** | 唯一 `pdb_id`: {s['pdb_id'].nunique()} | "
         f"唯一 `interface_key`: {s['interface_key'].nunique()}")
    emit(f"- split 分布: {s['split'].value_counts().to_dict()}")
    emit("")
    vc2 = s.groupby("interface_key").size()
    emit(f"- 每 `interface_key` 突变体数: 中位 {vc2.median():.0f}, 均值 {vc2.mean():.1f}, "
         f"最大 {vc2.max()}")
    for th in (3, 5, 10, 20):
        emit(f"- `interface_key` 数 (突变体 >= {th}): **{int((vc2 >= th).sum())}** "
             f"(覆盖 {int(vc2[vc2 >= th].sum())} 行)")
    emit("")
    emit(f"- **test 划分里只有 {int((s['split'] == 'test').sum())} 行** —— "
         f"这是 V2.3.2 那个 ρ=+0.5848 (n=34) 的来源, "
         f"样本量不足以支撑 within-target 基准。")

    # ---------- 3. 关键: 复合物内特征是否真的变化 ----------
    emit("")
    emit("## 3. **关键诊断**: 同一复合物内, 特征是否会变化?")
    emit("")
    emit("> 这决定 within-target 基准是否**在数学上可能**。")
    emit("> V2.3.2 用 `pdb_id` join (结构级平均), 若同一结构的多个突变体序列完全相同,")
    emit("> 它们会拿到**完全相同的特征** -> 全并列 -> 组内 Spearman 必然 NaN。")
    emit("")

    # 3a. 直接用缓存里的序列 md5 (真实反映序列是否不同)
    if CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        pdb_ids = z["pdb_ids"].astype(str)
        lig_md5, rec_md5 = z["lig_md5"].astype(str), z["rec_md5"].astype(str)
        splits = z["splits"].astype(str)
        emit("### 3a. ESM-2 缓存层面 (序列 md5)")
        emit("")
        cdf = pd.DataFrame({"pdb_id": pdb_ids, "lig_md5": lig_md5,
                            "rec_md5": rec_md5, "split": splits})
        g = cdf.groupby("pdb_id")
        both_var = g.apply(lambda x: (x["lig_md5"].nunique() > 1)
                           and (x["rec_md5"].nunique() > 1), include_groups=False)
        lig_var = g["lig_md5"].nunique()
        rec_var = g["rec_md5"].nunique()
        multi = g.size() > 1
        emit(f"- 有多行的 pdb_id: {int(multi.sum())} / {len(g)}")
        emit(f"- 其中**配体序列有变化**的: {int((lig_var > 1).sum())}")
        emit(f"- 其中**受体序列有变化**的: {int((rec_var > 1).sum())}")
        emit(f"- 其中**两边都变**的: {int(both_var.sum())}")
        emit(f"- **两边都不变** (结构级平均后完全无法区分) 的: "
             f"{int((multi & ~both_var & (lig_var == 1) & (rec_var == 1)).sum())}")
        # SKEMPI 子集内的 test
        sk_pdbs = set(s[s["split"] == "test"]["pdb_id"].astype(str))
        sub = cdf[cdf["pdb_id"].isin(sk_pdbs)]
        gs = sub.groupby("pdb_id")
        emit("")
        emit(f"- **SKEMPI test 子集** ({len(sk_pdbs)} 个 pdb, {len(sub)} 行):")
        for pid, grp in list(gs)[:15]:
            lv, rv = grp["lig_md5"].nunique(), grp["rec_md5"].nunique()
            emit(f"    - {pid}: {len(grp)} 行 | 唯一配体序列 {lv} | 唯一受体序列 {rv} | "
                 f"{'✅ 可区分' if (lv > 1 or rv > 1) else '❌ 不可区分'}")
    else:
        emit("- (无 ESM-2 缓存)")

    # 3b. 特征表层面: 同一 interface_key 内特征是否变化
    emit("")
    emit("### 3b. 特征表层面 (V2.2 的 60 物理特征)")
    emit("")
    feat_cols = [c for c in s.columns
                 if c not in ("interface_key", "pdb_id", "split", "real_dG",
                              "KD(M)", "Mutations", "Source Data Set", "Subgroup",
                              "Affinity Method", "Temperature(K)")
                 and s[c].dtype.kind in "fi"]
    if s["interface_key"].nunique() < len(s):
        dup = s.groupby("interface_key").filter(lambda x: len(x) > 1)
        n_vary = 0
        for c in feat_cols:
            if dup.groupby("interface_key")[c].nunique().gt(1).any():
                n_vary += 1
        emit(f"- 有重复 `interface_key` 的行: {len(dup)}")
        emit(f"- 在 `interface_key` 内**有任何变化**的物理特征数: "
             f"**{n_vary} / {len(feat_cols)}**")
        if n_vary == 0:
            emit("- -> ❌ **同一界面内所有物理特征完全相同** (结构级特征), 无法区分突变体")
        else:
            emit(f"- -> 有 {n_vary} 个特征在界面内变化")
    else:
        emit("- 每个 `interface_key` 只有 1 行 (无重复), 无需检查")

    # ---------- 4. 结论 ----------
    emit("")
    emit("## 4. 诊断结论与可行方案")
    emit("")
    if CACHE.exists() and 'both_var' in dir():
        n_multi = int(multi.sum())
        n_disc = int(both_var.sum())
        emit(f"* ESM-2 层面: 多行 pdb 中 **{n_disc} / {n_multi}** 两边序列都变 -> ")
        emit(f"  这些结构的突变体**可以**被序列嵌入区分。")
        ling = int(((lig_var > 1) & (multi)).sum())
        emit(f"  另一部分 ({ling} 个) 只有配体变、受体固定 —— 受体侧嵌入恒定的前提下,")
        emit(f"  区分能力仅来自配体侧。")
    emit("")
    emit("* **V2.2 的 test 划分对 within-target 不可用**: SKEMPI 的 test 只有 32 行。")
    emit("* **可行方案**: 用 **复合物级 GroupKFold 交叉验证** 让**每个 SKEMPI 行都有"
         "折外 (out-of-fold) 预测**, 再在复合物内部排序。")
    emit("  * 这样可用上 SKEMPI 的 **5109 行 / 285 个 pdb** (而非 32 行);")
    emit("  * group = `pdb_id` (或 interface_key), 保证**同复合物的突变体不跨折** ->")
    emit("    与 cross-target 划分的精神一致, **不泄漏**;")
    emit("  * 每个复合物内部排序 -> per-complex Spearman -> 汇总中位数/IQR。")

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT} ({OUT.stat().st_size} B)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
