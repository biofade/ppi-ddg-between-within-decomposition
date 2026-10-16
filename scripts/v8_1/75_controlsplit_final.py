# -*- coding: utf-8 -*-
"""scripts/v8_1/75_controlsplit_final.py -- 受控 split：median 口径复现 + paired attribution。

受控设计（来自 split4_report.md）：
  M0/M2 = pooled baseline；M1 = M0 + DGaux OOF
  模型定义/超参与 V5 Gate A 完全一致（300 树, PCA=128, 无早停）
  seeds = 42/43/44
  L1=random, L2=complex, L3=seq40, L4=family-OOD
  → 同一个 (model, seed) 在四档下都有结果，因此可做 **paired** 比较。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
d = pd.read_csv("v8_1/all_models_by_rowset.csv", low_memory=False)
# Table 2 / H4 的口径 = 仅 V6_2_split4（已核实可逐位复现 FROZEN H4 值）。
# V6_2_1_donly 是另一变体（D-only 控制），不参与受控 split 比较。
v = d[(d.n_rows == 7523) & (d.artifact == "V6_2_split4")].copy()
v["level"] = v.key.str.extract(r"^(L\d)")
v["model"] = v.key.str.extract(r"_(M\d)")
v["seed"] = v.key.str.extract(r"_(s\d+)$")
LEVEL = {"L1": "Random", "L2": "Complex", "L3": "seq40", "L4": "Family-OOD"}
v["regime"] = v.level.map(LEVEL)

print("=== 1) 每档 marginal medians（Table 2 口径）===")
med = {}
for lv in ("L1", "L2", "L3", "L4"):
    s = v[v.level == lv]
    mr, mb, mw = (s.overall_pearson.median(), s.C_B_row.median(), s.C_W_row.median())
    ms = s.macro_spearman.median(); mr_ = s.rmse.median()
    med[lv] = dict(n=len(s), median_r=mr, median_C_B=mb, median_C_W=mw,
                   macro=ms, rmse=mr_, sum_med=mb + mw, gap=mb + mw - mr)
    print(f"  {LEVEL[lv]:11s} n={len(s):2d}  r={mr:.4f}  C_B={mb:+.4f}  C_W={mw:+.4f}"
          f"  sum={mb+mw:.4f}  gap={mb+mw-mr:+.4f}  macro={ms:.4f}  rmse={mr_:.4f}")

print()
print("=== 2) per-configuration 恒等式 ===")
res = v.overall_pearson - v.C_B_row - v.C_W_row
print(f"  n={len(v)}  max|r-C_B-C_W|={np.abs(res).max():.3e}  "
      f"全部<=1e-12: {bool((np.abs(res)<=1e-12).all())}")

print()
print("=== 3) paired：Random(L1) -> Family-OOD(L4)，按 (model, seed, artifact) 配对 ===")
a = v[v.level == "L1"].set_index(["artifact", "model", "seed"])
b = v[v.level == "L4"].set_index(["artifact", "model", "seed"])
common = a.index.intersection(b.index)
print(f"  可配对 configuration 数 = {len(common)}")
print(f"  L1 独有={len(a.index.difference(b.index))}  L4 独有={len(b.index.difference(a.index))}")

pairs = []
for k in sorted(common):
    A, B = a.loc[k], b.loc[k]
    if isinstance(A, pd.DataFrame) or isinstance(B, pd.DataFrame):
        continue
    pairs.append(dict(artifact=k[0], model=k[1], seed=k[2],
                      r_L1=A.overall_pearson, r_L4=B.overall_pearson,
                      CB_L1=A.C_B_row, CB_L4=B.C_B_row,
                      CW_L1=A.C_W_row, CW_L4=B.C_W_row))
P = pd.DataFrame(pairs)
P["d_r"] = P.r_L4 - P.r_L1
P["d_CB"] = P.CB_L4 - P.CB_L1
P["d_CW"] = P.CW_L4 - P.CW_L1
P["paired_resid"] = P.d_r - P.d_CB - P.d_CW
print()
print(P[["artifact", "model", "seed", "r_L1", "r_L4", "d_r", "d_CB", "d_CW",
         "paired_resid"]].to_string(index=False))
print()
print(f"  paired 恒等式 max|d_r - d_CB - d_CW| = {P.paired_resid.abs().max():.3e}")
print()
print(f"  median d_r   = {P.d_r.median():+.4f}   mean = {P.d_r.mean():+.4f}")
print(f"  median d_CB  = {P.d_CB.median():+.4f}   mean = {P.d_CB.mean():+.4f}")
print(f"  median d_CW  = {P.d_CW.median():+.4f}   mean = {P.d_CW.mean():+.4f}")
r_paired = abs(P.d_CB.median()) / (abs(P.d_CB.median()) + abs(P.d_CW.median()))
r_marg = (abs(med['L4']['median_C_B'] - med['L1']['median_C_B']) /
          (abs(med['L4']['median_C_B'] - med['L1']['median_C_B']) +
           abs(med['L4']['median_C_W'] - med['L1']['median_C_W'])))
percfg = (P.d_CB.abs() / (P.d_CB.abs() + P.d_CW.abs())).median()
print()
print(f"  归因比例（marginal medians，论文现用）= {r_marg:.4f}  -> {100*r_marg:.1f}%")
print(f"  归因比例（paired medians）           = {r_paired:.4f}  -> {100*r_paired:.1f}%")
print(f"  归因比例（per-config median）        = {percfg:.4f}  -> {100*percfg:.1f}%")

# bootstrap
rng = np.random.default_rng(42)
B = 5000
idx = np.arange(len(P))
bs = {k: [] for k in ("d_r", "d_CB", "d_CW", "ratio")}
for _ in range(B):
    j = rng.choice(idx, len(idx), replace=True)
    q = P.iloc[j]
    bs["d_r"].append(q.d_r.median()); bs["d_CB"].append(q.d_CB.median())
    bs["d_CW"].append(q.d_CW.median())
    den = abs(q.d_CB.median()) + abs(q.d_CW.median())
    bs["ratio"].append(abs(q.d_CB.median()) / den if den else np.nan)
ci = {k: (float(np.nanpercentile(x, 2.5)), float(np.nanpercentile(x, 97.5)))
      for k, x in bs.items()}
print()
print("=== bootstrap (5,000, configuration 级, seed 42) ===")
print(f"  median d_r   95%CI = {ci['d_r'][0]:+.4f} .. {ci['d_r'][1]:+.4f}")
print(f"  median d_CB  95%CI = {ci['d_CB'][0]:+.4f} .. {ci['d_CB'][1]:+.4f}")
print(f"  median d_CW  95%CI = {ci['d_CW'][0]:+.4f} .. {ci['d_CW'][1]:+.4f}")
print(f"  归因比例     95%CI = {ci['ratio'][0]:.4f} .. {ci['ratio'][1]:.4f}")

out = dict(
    universe="U_MAS_7523 controlled split (V6.2)",
    controlled_design=("same model specs and hyperparameters (300 trees, PCA=128, "
                       "no early stopping); seeds 42/43/44; only the split changes"),
    per_configuration_identity=dict(
        n=int(len(v)), max_abs_residual=float(np.abs(res).max()),
        all_within_1e_12=bool((np.abs(res) <= 1e-12).all())),
    marginal_medians={LEVEL[k]: {kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else vv)
                                for kk, vv in mv.items()} for k, mv in med.items()},
    paired=dict(n_paired=int(len(P)),
                max_abs_paired_residual=float(P.paired_resid.abs().max()),
                delta_r_median=float(P.d_r.median()),
                delta_C_B_median=float(P.d_CB.median()),
                delta_C_W_median=float(P.d_CW.median()),
                attribution_marginal_median=float(r_marg),
                attribution_paired_median=float(r_paired),
                attribution_perconfig_median=float(percfg),
                bootstrap_ci=ci),
)
Path("v8_1/CONTROLSPLIT_IDENTITY.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print("\nwrote v8_1/CONTROLSPLIT_IDENTITY.json")
