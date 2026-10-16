# -*- coding: utf-8 -*-
"""scripts/v8_1/143_recompute_feasibility.py -- 键级分组重训可行性评估。

问题：受控 split（7,523 rows）的 1,405 个 (pdb,mutation) key 带冲突 ΔΔG，
      且 2,011 个 key 重复。同一 key 的不同 row 可能落在不同 fold，
      造成 key-level leakage。要排除需「按 key 分组」重新划分并重算全部主结果。

本脚本不去重训，只**如实评估**：
  1. 需要重训的模型有哪些、能否在本机 CPU 上完成
  2. artifact 中是否保存了 fold 归属（决定能否只重划分而不重训）
  3. 按 key 分组后各 regime 的规模是否仍满足 n>=5 等门限
  4. 重划分对已有冻结数字的预期影响面
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
A = R / "v8_1"
out = {}

print("=" * 100)
print("1) 现有 artifact 是否保存 fold 归属 / 训练配置")
print("=" * 100)
cands = [
    "experiments/v7_1/phase2_rank/z_predictions.npz",
    "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
    "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
]
for c in cands:
    p = R / c
    if not p.exists():
        print(f"  [缺失] {c}"); continue
    if p.suffix == ".npz":
        z = np.load(p, allow_pickle=True)
        print(f"  {c}")
        print(f"    keys = {list(z.keys())}")
    else:
        d = pd.read_csv(p, low_memory=False, nrows=400)
        cols = list(d.columns)
        fold_like = [x for x in cols if any(k in x.lower()
                     for k in ("fold", "split", "regime", "seed", "partition"))]
        print(f"  {c}")
        print(f"    列数={len(cols)}  fold/split 相关列 = {fold_like or '无'}")

print()
print("=" * 100)
print("2) 主模型是什么、能否本地重训")
print("=" * 100)
reg = R / "scripts/v8/03_model_registry.py"
if reg.exists():
    txt = reg.read_text(encoding="utf-8", errors="replace")
    import re
    for m in re.finditer(r'"(model|name|family|backend|predictor)"\s*:\s*"([^"]{0,60})"', txt):
        print(f"    {m.group(1):12s} = {m.group(2)}")
print()
print("  注：主分析使用 experiments/v7_1/phase2_rank/z_predictions.npz 的 y 与各模型分数；")
print("      若该 npz 只存最终预测而不存 fold 归属，则无法「只重划分、不重训」。")

print()
print("=" * 100)
print("3) 按 key 分组重划分后的规模（用受控 split 实测）")
print("=" * 100)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
print(f"  受控 split 行数 = {len(mas):,}   列 = {list(mas.columns)[:12]}")
keycols = [c for c in mas.columns if "key" in c.lower()]
print(f"  key 相关列: {keycols or '无（需从 pdb+mutation 重建）'}")

# 重建 key
ms = pd.read_csv(R / "data/v5_0/mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
agg = ms.sort_values("pair_row").groupby("pair_row", as_index=False).agg(
    pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
    pos=("ann_pos1", "first"), mut=("mut_aa", "first"))
mas["pair_row"] = pd.to_numeric(mas.get("pair_row"), errors="coerce")
J = mas.merge(agg, on="pair_row", how="left")
J["key"] = (J.pdb4.astype(str) + "|"
            + pd.to_numeric(J.pos, errors="coerce").astype("float")
              .astype("Int64").astype(str) + "|"
            + J.mut.astype(str).str.upper())
n_keys = J.key.nunique()
dup = int((J.key.value_counts() > 1).sum())
print(f"  行数 = {len(J):,}   不同 key = {n_keys:,}   重复 key = {dup:,}")
out["rows"] = int(len(J)); out["distinct_keys"] = int(n_keys); out["dup_keys"] = dup

# 若按 key 分组，target 数会不会掉
tg = "target_id_canonical" if "target_id_canonical" in J.columns else None
if tg:
    print(f"  targets（行级） = {J[tg].nunique()}")
    gk = J.groupby("key")[tg].nunique()
    print(f"  同一 key 跨 target 的 key 数 = {int((gk > 1).sum()):,} "
          f"(这些 key 分组后无法同时属于两个 target)")
    out["keys_crossing_targets"] = int((gk > 1).sum())

json.dump(out, open(A / "RECOMPUTE_FEASIBILITY.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print()
print("wrote v8_1/RECOMPUTE_FEASIBILITY.json")
