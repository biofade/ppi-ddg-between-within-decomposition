# -*- coding: utf-8 -*-
"""scripts/v8_1/74_controlsplit_identity_and_paired.py -- 受控 split 的精确恒等式与 paired attribution。

回答三件事：
  1. 每个 configuration 上 r - C_B - C_W 的残差（应为 machine precision）
  2. Table 2 展示的 median 为何不满足 median(r)=median(C_B)+median(C_W)
  3. paired analysis：同一 configuration 从 Random -> Family-OOD 的
     delta_r / delta_C_B / delta_C_W，并检验 delta_r = delta_C_B + delta_C_W
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
d = pd.read_csv("v8_1/all_models_by_rowset.csv", low_memory=False)

# V6.2 受控 split：7,523 行 universe，四档 split_regime 各含同一批 configuration
v = d[d.n_rows == 7523].copy()
print(f"=== 7,523-row universe 的 configurations: {len(v)} ===")
print(v.groupby("split_regime").size().to_string())
print()

# ---- 1) per-configuration 恒等式 ----
res = v.overall_pearson - v.C_B_row - v.C_W_row
print("=== per-configuration 恒等式 r = C_B + C_W ===")
print(f"  n = {len(v)}")
print(f"  max |residual|            = {np.abs(res).max():.3e}")
print(f"  decomp_residual 列 max    = {v.decomp_residual.abs().max():.3e}")
print(f"  全部 <= 1e-12 ?           = {bool((np.abs(res) <= 1e-12).all())}")
print()

# ---- 2) marginal medians ----
print("=== 四档的 marginal medians（Table 2 展示口径）===")
rows = []
for regime in ("random", "complex", "seq40", "family_ood"):
    s = v[v.split_regime == regime]
    if not len(s):
        continue
    mr, mb, mw = s.overall_pearson.median(), s.C_B_row.median(), s.C_W_row.median()
    rows.append(dict(regime=regime, n_models=len(s),
                     median_r=mr, median_C_B=mb, median_C_W=mw,
                     sum_medians=mb + mw, gap=(mb + mw) - mr))
    print(f"  {regime:11s} median r={mr:.4f}  median C_B={mb:+.4f}  "
          f"median C_W={mw:+.4f}  sum={mb+mw:.4f}  gap={mb+mw-mr:+.4f}")
print()
print("  >>> gap != 0 是 median 的非线性所致，不是恒等式失效")
print()

# ---- 3) paired analysis（Random vs Family-OOD）----
rnd = v[v.split_regime == "random"].set_index("key")
ood = v[v.split_regime == "family_ood"].set_index("key")
common = rnd.index.intersection(ood.index)
print("=== paired analysis: Random -> Family-OOD ===")
print(f"  两档共同 configuration（按 key 配对）= {len(common)}")
print(f"  random 独有 = {len(rnd.index.difference(ood.index))}, "
      f"family_ood 独有 = {len(ood.index.difference(rnd.index))}")
print()

if len(common):
    dr = ood.loc[common, "overall_pearson"] - rnd.loc[common, "overall_pearson"]
    db = ood.loc[common, "C_B_row"] - rnd.loc[common, "C_B_row"]
    dw = ood.loc[common, "C_W_row"] - rnd.loc[common, "C_W_row"]
    pres = dr - db - dw
    print(f"  paired 恒等式 max |delta_r - delta_C_B - delta_C_W| = "
          f"{np.abs(pres).max():.3e}")
    print()
    print(f"  delta_r     : median={dr.median():+.4f}  mean={dr.mean():+.4f}  "
          f"[min={dr.min():+.4f}, max={dr.max():+.4f}]")
    print(f"  delta_C_B   : median={db.median():+.4f}  mean={db.mean():+.4f}  "
          f"[min={db.min():+.4f}, max={db.max():+.4f}]")
    print(f"  delta_C_W   : median={dw.median():+.4f}  mean={dw.mean():+.4f}  "
          f"[min={dw.min():+.4f}, max={dw.max():+.4f}]")
    print()
    # attribution ratio（per-configuration）
    denom = np.abs(db) + np.abs(dw)
    ratio = np.where(denom > 0, np.abs(db) / denom, np.nan)
    print("  per-configuration 归因比例 |dCB|/(|dCB|+|dCW|):")
    print(f"    median = {np.nanmedian(ratio):.4f}   mean = {np.nanmean(ratio):.4f}")
    print()
    # 与 marginal 做法的对比
    med_ratio = abs(v[v.split_regime=='family_ood'].C_B_row.median()
                    - v[v.split_regime=='random'].C_B_row.median())
    med_ratio_d = med_ratio + abs(v[v.split_regime=='family_ood'].C_W_row.median()
                                  - v[v.split_regime=='random'].C_W_row.median())
    print(f"  marginal-median 口径归因比例 = {med_ratio/med_ratio_d:.4f}  "
          f"(= 论文中的 81.1%)")
    print(f"  paired-median 口径归因比例   = "
          f"{abs(db.median())/(abs(db.median())+abs(dw.median())):.4f}")
    print(f"  per-config median 归因比例   = {np.nanmedian(ratio):.4f}")

# ---- bootstrap CI（configuration 级别重采样）----
if len(common) >= 5:
    rng = np.random.default_rng(42)
    B = 5000
    idx = np.arange(len(common))
    bs_ratio, bs_dr, bs_db, bs_dw = [], [], [], []
    for _ in range(B):
        j = rng.choice(idx, len(idx), replace=True)
        b_, w_ = db.to_numpy()[j], dw.to_numpy()[j]
        r_ = dr.to_numpy()[j]
        bs_dr.append(np.median(r_)); bs_db.append(np.median(b_)); bs_dw.append(np.median(w_))
        dd = np.abs(np.median(b_)) + np.abs(np.median(w_))
        bs_ratio.append(np.abs(np.median(b_)) / dd if dd else np.nan)
    q = lambda a: (float(np.nanpercentile(a, 2.5)), float(np.nanpercentile(a, 97.5)))
    print()
    print("=== bootstrap（5,000 次，configuration 级重采样，seed 42）===")
    print(f"  median delta_r    95% CI = {q(bs_dr)[0]:+.4f} .. {q(bs_dr)[1]:+.4f}")
    print(f"  median delta_C_B  95% CI = {q(bs_db)[0]:+.4f} .. {q(bs_db)[1]:+.4f}")
    print(f"  median delta_C_W  95% CI = {q(bs_dw)[0]:+.4f} .. {q(bs_dw)[1]:+.4f}")
    print(f"  归因比例          95% CI = {q(bs_ratio)[0]:.4f} .. {q(bs_ratio)[1]:.4f}")
    bs = dict(bs_dr=q(bs_dr), bs_db=q(bs_db), bs_dw=q(bs_dw), bs_ratio=q(bs_ratio))
else:
    bs = None

out = {
    "universe": "U_MAS_7523 controlled-split configurations",
    "n_configurations_total": int(len(v)),
    "identity": {
        "max_abs_residual": float(np.abs(res).max()),
        "all_within_1e-12": bool((np.abs(res) <= 1e-12).all()),
        "note": "r_pooled = C_B + C_W 在每个 configuration 上精确成立",
    },
    "marginal_medians": rows,
    "paired": {
        "n_paired": int(len(common)),
        "paired_identity_max_abs_residual": float(np.abs(pres).max()) if len(common) else None,
        "delta_r_median": float(dr.median()) if len(common) else None,
        "delta_C_B_median": float(db.median()) if len(common) else None,
        "delta_C_W_median": float(dw.median()) if len(common) else None,
        "attribution_marginal_median": float(med_ratio / med_ratio_d) if len(common) else None,
        "attribution_paired_median": (float(abs(db.median()) /
                                           (abs(db.median()) + abs(dw.median())))
                                      if len(common) else None),
        "attribution_perconfig_median": float(np.nanmedian(ratio)) if len(common) else None,
        "bootstrap": bs,
    },
}
Path("v8_1/CONTROLSPLIT_IDENTITY.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print("\nwrote v8_1/CONTROLSPLIT_IDENTITY.json")
