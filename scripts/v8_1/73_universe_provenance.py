# -*- coding: utf-8 -*-
"""scripts/v8_1/73_universe_provenance.py -- 判定 4,473 与 7,523 两个 universe 的关系。

背景（来自代码与记录，非推测）：
  * `build_pair_features` 产出 7,523 行真 ΔΔG forward 配对 / 230 复合物 / 112 family
    -> `matched_analysis_set.csv`，0 行被剔除（见 coverage_report.md）
  * `canonical_target_labels.csv` 由
    `scripts/v7_1_1/02_reconciliation_final.py` L99-106 生成：
    取 `block_L_local.csv` 的 pair_row 列表，到 canonical(7523) 里取
    `pdb_id` / `family_label`。它是 7,523 的一个**特征覆盖子集**。

本脚本用真实文件验证：
  1. 4,473 的 pair_row 是否确实是 7,523 的 row 的子集
  2. 两个 universe 的行数 / target 数 / family 数 / n>=k 分布
  3. 4,473 相对 7,523 缺了哪些 target / family
  4. 7,523 内部 (pdb, mutation) 重复观测的规模与 ddG 冲突规模
输出 machine-readable artifact。
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
MAS = R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv"
CAN = R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv"

mas = pd.read_csv(MAS, low_memory=False)
can = pd.read_csv(CAN, low_memory=False)
print(f"matched_analysis_set : {len(mas)} 行, {mas.pdb_id.nunique()} targets(pdb), "
      f"{mas.family_label.nunique()} families")
print(f"canonical_target_labels: {len(can)} 行, "
      f"{can.target_id_canonical.nunique()} targets, "
      f"{can.family_id_canonical.nunique()} families")
print()

# ---- 1) pair_row 子集关系 ----
rows_can = set(can.pair_row.tolist())
rows_mas = set(mas.row.tolist())
print("=== pair_row 关系 ===")
print(f"  canonical rows 全部存在于 7523 内 : {rows_can <= rows_mas}")
print(f"  canonical 行数 {len(rows_can)} / 7523 行数 {len(rows_mas)}")
print(f"  7523 \\ canonical = {len(rows_mas - rows_can)} 行")
print()

# ---- 2) n>=k 分布 ----
def dist(ids):
    s = ids.value_counts()
    return {f"n_ge{k}": int((s >= k).sum()) for k in (1, 2, 3, 5, 10)}


d_mas = dist(mas.pdb_id)
d_can = dist(can.target_id_canonical)
print("=== target n>=k 分布 ===")
print(f"  7523 集 : {d_mas}")
print(f"  4473 集 : {d_can}")
print()

# ---- 3) 缺失的 target / family ----
t_mas, t_can = set(mas.pdb_id.astype(str)), set(can.target_id_canonical.astype(str))
f_mas, f_can = set(mas.family_label.astype(str)), set(can.family_id_canonical.astype(str))
print("=== 4473 相对 7523 缺失 ===")
print(f"  target: {len(t_mas - t_can)} 个不在 4473 中")
print(f"  family: {len(f_mas - f_can)} 个不在 4473 中")
print()

# ---- 4) 7523 内重复观测 ----
mas["key"] = mas.pdb_id.astype(str) + "|" + mas.mutation.astype(str)
g = mas.groupby("key").agg(n=("ddG", "size"), nddg=("ddG", "nunique"))
dup_rows = int((g.n > 1).sum())
conflict = int((g.nddg > 1).sum())
print("=== 7523 集内部重复观测 ===")
print(f"  唯一 (pdb, mutation) 键     : {len(g)}")
print(f"  出现多行的键                : {dup_rows}")
print(f"  同一键 ddG 取值不一致的键   : {conflict}")
print(f"  唯一 state_id               : {mas.state_id.nunique()}")
print()

art = {
    "universes": {
        "U_MAS_7523": {
            "source": str(MAS).replace("\\", "/"),
            "producer": "scripts/v6_2/02b_canonical_set.py",
            "definition": ("build_pair_features 的真 ddG forward 配对全集；"
                           "保留四档 split 均可用的行；coverage_report 记录剔除 0 行"),
            "n_rows": int(len(mas)),
            "n_targets_pdb": int(mas.pdb_id.nunique()),
            "n_families": int(mas.family_label.nunique()),
            "n_unique_state_id": int(mas.state_id.nunique()),
            "n_unique_pdb_mutation_keys": int(len(g)),
            "n_keys_with_multiple_rows": dup_rows,
            "n_keys_with_conflicting_ddg": conflict,
            **d_mas,
        },
        "U_CANON_4473": {
            "source": str(CAN).replace("\\", "/"),
            "producer": "scripts/v7_1_1/02_reconciliation_final.py (L99-L106)",
            "definition": ("取 block_L_local.csv 的 pair_row 列表，映射到 "
                           "matched_analysis_set 的 pdb_id / family_label；"
                           "是 7523 的特征覆盖子集"),
            "n_rows": int(len(can)),
            "n_targets": int(can.target_id_canonical.nunique()),
            "n_families": int(can.family_id_canonical.nunique()),
            **d_can,
        },
    },
    "relation": {
        "canonical_rows_subset_of_mas_rows": bool(rows_can <= rows_mas),
        "rows_in_mas_not_in_canonical": int(len(rows_mas - rows_can)),
        "targets_in_mas_not_in_canonical": int(len(t_mas - t_can)),
        "families_in_mas_not_in_canonical": int(len(f_mas - f_can)),
    },
}
out = R / "v8_1" / "UNIVERSE_PROVENANCE.json"
out.write_text(json.dumps(art, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"wrote {out}")
