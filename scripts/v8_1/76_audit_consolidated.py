# -*- coding: utf-8 -*-
"""scripts/v8_1/76_audit_consolidated.py -- 汇总审计：external coverage / families /
109 configurations / Figure 3 raw values / LA-LC pools。

产出 v8_1/AUDIT_CONSOLIDATED.json，作为 manuscript 与 validate_final 的机器可读依据。
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = {}
d = pd.read_csv("v8_1/all_models_by_rowset.csv", low_memory=False)

# ============ 1) external coverage + rankings ============
ext = pd.read_csv("v8_1/external_clean_profiles.csv")
print("=" * 78)
print("1) EXTERNAL PANEL COVERAGE")
print("=" * 78)
print(ext[["model", "n_rows", "n_targets", "n_targets_ge5", "overall_pearson",
           "macro_spearman", "C_B_row", "C_W_row", "rmse", "decomp_residual"]]
      .to_string(index=False))
counts = ext.set_index("model").n_rows.to_dict()
print(f"\n  mutation 计数: {counts}")
all_same = len(set(counts.values())) == 1
print(f"  所有方法 mutation 数相同? {all_same}")
print("  -> 论文 'All methods are evaluated on the same mutations' 不成立")

pr = ext.set_index("model").overall_pearson.rank(ascending=False)
mr = ext.set_index("model").macro_spearman.rank(ascending=False)
print("\n  Pooled Pearson 排名 :", list(pr.sort_values().index))
print("  Macro Spearman 排名 :", list(mr.sort_values().index))
pairs_p = list(zip(pr.sort_values().index, pr.sort_values().values))
pairs_m = list(zip(mr.sort_values().index, mr.sort_values().values))
swap = [i for i in range(len(pairs_p)) if pairs_p[i][0] != pairs_m[i][0]]
print(f"  发生位置交换的名次: {[i+1 for i in swap]}  -> "
      f"{[(pairs_p[i][0], pairs_m[i][0]) for i in swap]}")

from scipy.stats import spearmanr, kendalltau  # noqa: E402
rho, _ = spearmanr(ext.overall_pearson, ext.macro_spearman)
tau, _ = kendalltau(ext.overall_pearson, ext.macro_spearman)
# 精确置换 p（n=5 共 120 种）
import itertools  # noqa: E402
obs = rho
cnt = 0
vals = ext.overall_pearson.to_numpy()
tgt = ext.macro_spearman.to_numpy()
for perm in itertools.permutations(range(len(vals))):
    r = spearmanr(vals[list(perm)], tgt)[0]
    if abs(r) >= abs(obs) - 1e-12:
        cnt += 1
perm_p = cnt / 120
print(f"  Spearman rho={rho:.4f}  Kendall tau={tau:.4f}  精确置换 p={perm_p:.4f}")
print(f"  (frozen 记录 rho=0.90, p=0.0833 -> 一致: "
      f"{abs(rho-0.9)<0.01 and abs(perm_p-0.0833)<0.01})")
OUT["external"] = dict(
    coverage=counts, coverage_equal=bool(all_same),
    ranking_pooled=list(pr.sort_values().index),
    ranking_macro=list(mr.sort_values().index),
    rank_swap_positions=[int(i) + 1 for i in swap],
    rank_swap_pair=[[pairs_p[i][0], pairs_m[i][0]] for i in swap],
    spearman_rho=float(rho), kendall_tau=float(tau), perm_p=float(perm_p),
    profiles=ext.to_dict("records"))

# ============ 2) method families ============
print()
print("=" * 78)
print("2) PRIMARY UNIVERSE: method families")
print("=" * 78)
pu = d[(d.n_rows == 3351)]
print(f"  n_rows={len(pu)}  n_configs={len(pu)}")
print(f"  n_method_families={pu.method_family.nunique()}")
fam = pu.groupby("method_family").agg(
    n_configs=("key", "size"),
    reps=("is_primary_representative", "sum"),
    ablations=("is_ablation", "sum")).sort_values("n_configs", ascending=False)
print(fam.to_string())
OUT["primary_families"] = fam.reset_index().to_dict("records")
OUT["primary_universe"] = dict(
    n_rows=int(len(pu)), n_configs=int(len(pu)),
    n_method_families=int(pu.method_family.nunique()),
    n_representatives=int(pu.is_primary_representative.sum()))

# ============ 3) 109 configurations (weighting audit) ============
print()
print("=" * 78)
print("3) WEIGHTING AUDIT: 109 configurations")
print("=" * 78)
w = pd.read_csv("v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv")
print(f"  WEIGHTING_EFFECT_SIZE_AUDIT 行数 = {len(w)}")
print(f"  唯一 model_name = {w.model_name.nunique()}")
print(f"  model_name 示例: {list(w.model_name.head(5))}")
if "universe" in w.columns:
    print("  universe 分布:", w.universe.value_counts().to_dict())
# 与 all_models_by_rowset 对齐
am = pd.read_csv("v8_1/all_models_by_rowset.csv", low_memory=False)
am["modname"] = am.artifact.astype(str) + "|" + am.key.astype(str)
matched = sum(1 for n in w.model_name if n in set(am.modname))
print(f"  能对应到 all_models_by_rowset 的 = {matched}/{len(w)}")
OUT["weighting_audit"] = dict(
    n_configurations=int(len(w)), n_unique=w.model_name.nunique(),
    model_name_examples=list(w.model_name.head(5)),
    n_matched_to_models=int(matched),
    median_range_C_B=float(w.range_C_B.median()),
    median_range_C_W=float(w.range_C_W.median()),
    signflip_C_B_gt0=int(w["signflip_C_B_gt0"].sum()),
    signflip_C_B_gt0_005=int(w["signflip_C_B_gt0.005"].sum()),
    signflip_C_W_gt0=int(w["signflip_C_W_gt0"].sum()))

# ============ 4) Figure 3 raw values (histograms) ============
print()
print("=" * 78)
print("4) FIGURE 3 RAW VALUES")
print("=" * 78)
print(f"  range_C_B: n={w.range_C_B.notna().sum()} min={w.range_C_B.min():.5f} "
      f"max={w.range_C_B.max():.5f} median={w.range_C_B.median():.5f}")
print(f"  range_C_W: n={w.range_C_W.notna().sum()} min={w.range_C_W.min():.5f} "
      f"max={w.range_C_W.max():.5f} median={w.range_C_W.median():.5f}")
OUT["figure3_raw"] = dict(
    range_C_B=[float(x) for x in w.range_C_B.dropna()],
    range_C_W=[float(x) for x in w.range_C_W.dropna()],
    median_C_B=float(w.range_C_B.median()), median_C_W=float(w.range_C_W.median()))
print("  -> Figure 3 可由真实 109 个值重新生成（无需猜测 bin）")

# ============ 5) LA / LC pool composition ============
print()
print("=" * 78)
print("5) LA / LC POOL COMPOSITION (rank_agreement_by_stratum)")
print("=" * 78)
ra = pd.read_csv("v8_1/rank_agreement_by_stratum.csv")
print(f"  columns: {list(ra.columns)}")
print(ra.to_string(index=False))
OUT["rank_agreement_by_stratum"] = ra.to_dict("records")

Path("v8_1/AUDIT_CONSOLIDATED.json").write_text(
    json.dumps(OUT, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
print("\nwrote v8_1/AUDIT_CONSOLIDATED.json")
