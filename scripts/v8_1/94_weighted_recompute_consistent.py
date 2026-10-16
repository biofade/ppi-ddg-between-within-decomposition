# -*- coding: utf-8 -*-
"""scripts/v8_1/94_weighted_recompute_consistent.py -- 用口径一致的加权分解重算 109 个配置。

问题（已定位）：scripts/v8/09_weighting_audit.weighted_metrics 中
    sigma 由 **行加权** 方差给出   vy = sum_ci w_ci (y_ci - m_y)^2 / W
    但 C_B 用 **target 权重** w_t 求和，却减 **行加权** 均值 m_y
两者不是同一加权测度，因此 r = C_B + C_W 只在 ROW scheme 精确成立。

正确做法（本脚本）：统一使用同一个加权测度
    a_c  = scheme_weights(n_c)          # ROW n_c / EQUAL 1 / SQRT sqrt(n_c)
                                        # LOG log1p(n_c) / CAP30 min(n_c,30)
    w_ci = a_c                          # target 内所有行同权（与代码一致）
    W    = sum_ci w_ci = sum_c n_c a_c
    m_x  = sum_ci w_ci x_ci / W
    sigma_x^2 = sum_ci w_ci (x_ci - m_x)^2 / W
    By_c = sum_i w_ci y_ci / sum_i w_ci     # 组内加权均值
    C_B  = sum_c  wbar_c (By_c - m_y)(Bp_c - m_p) / (W sigma_y sigma_p)
           where wbar_c = (sum_i w_ci) = n_c a_c
    C_W  = sum_ci w_ci (y_ci - By_c)(p_ci - Bp_c) / (W sigma_y sigma_p)
则 r_w = C_B,w + C_W,w 对每个 scheme、每个 configuration 精确成立。

输出：
  v8_1/WEIGHTING_CONSISTENT.csv        新口径的逐配置结果
  v8_1/WEIGHTING_CONSISTENT_AUDIT.csv  新口径的 range / signflip 汇总
  v8_1/WEIGHTING_CONSISTENT.json       校验信息
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]


def scheme_weights(n, scheme):
    c = np.asarray(n, float)
    return {"ROW": c, "EQUAL": np.ones_like(c), "SQRT": np.sqrt(c),
            "LOG": np.log1p(c), "CAP30": np.minimum(c, 30.0)}[scheme]


def consistent(y, p, tgt, scheme):
    """口径一致的加权分解。返回 (r, C_B, C_W)。"""
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv).astype(float)
    a = scheme_weights(n_t, scheme)          # a_c
    w = a[inv]                               # w_ci = a_c
    W = w.sum()
    my = np.sum(w * y) / W
    mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W
    vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 3
    sd = np.sqrt(vy * vp)
    # 组内加权均值 与 组内权重和 wbar_c = n_c * a_c
    By = np.empty(len(uq)); Bp = np.empty(len(uq)); wbar = np.empty(len(uq))
    for k in range(len(uq)):
        msk = inv == k
        wbar[k] = w[msk].sum()
        By[k] = np.sum(w[msk] * y[msk]) / wbar[k]
        Bp[k] = np.sum(w[msk] * p[msk]) / wbar[k]
    C_B = np.sum(wbar * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return float(C_B + C_W), float(C_B), float(C_W)


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
zcan = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y_CAN = np.asarray(zcan["y"], float)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
TGT_MAS = mas.pdb_id.astype(str).to_numpy()
S = pd.read_csv(R / "experiments/v8/BENCHMARK_SCOREBOARD.csv", low_memory=False)
old = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv",
                  low_memory=False).drop_duplicates("model_name").set_index("model_name")

cache, recs = {}, []
LAB = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
       "z_raw", "y_cal", "family", "seed", "masks"}
for r in S.itertuples():
    rel = r.prediction_file
    if rel not in cache:
        try:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        except Exception:  # noqa: BLE001
            cache[rel] = None
    zz = cache[rel]
    if zz is None or r.model_name not in zz.files or r.model_name in LAB \
            or r.model_name.lower().endswith("_y"):
        continue
    a = np.asarray(zz[r.model_name])
    if a.ndim != 1:
        continue
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        continue
    if r.join_mode.startswith("canonical"):
        y, tgt = Y_CAN, TGT
    else:
        y, tgt = Y_MAS, TGT_MAS
    if len(a) != len(y):
        continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        continue
    rec = {"model_name": r.model_name, "universe": r.universe}
    for s in SCHEMES:
        rr, cb, cw = consistent(yy, pp, tt, s)
        rec[f"pearson_{s}"] = rr; rec[f"C_B_{s}"] = cb; rec[f"C_W_{s}"] = cw
    # 记录旧口径用于对照
    if r.model_name in old.index:
        for s in SCHEMES:
            rec[f"old_C_B_{s}"] = old.loc[r.model_name, f"C_B_{s}"]
            rec[f"old_C_W_{s}"] = old.loc[r.model_name, f"C_W_{s}"]
    recs.append(rec)

D = pd.DataFrame(recs)
cbcols = [f"C_B_{s}" for s in SCHEMES]; cwcols = [f"C_W_{s}" for s in SCHEMES]
print(f"重算配置数 = {len(D)}")

# ---- 恒等式校验 ----
print()
print("=" * 78)
print("恒等式校验（新口径）")
print("=" * 78)
resid = {}
for s in SCHEMES:
    e = (D[f"pearson_{s}"] - D[f"C_B_{s}"] - D[f"C_W_{s}"]).abs().dropna()
    resid[s] = float(e.max())
    print(f"  {s:6s}  max|r - C_B - C_W| = {e.max():.3e}   "
          f"全部 <= 1e-12: {bool((e <= 1e-12).all())}")

# ---- 新 audit ----
A = pd.DataFrame({"model_name": D.model_name})
A["range_C_B"] = D[cbcols].max(axis=1) - D[cbcols].min(axis=1)
A["range_C_W"] = D[cwcols].max(axis=1) - D[cwcols].min(axis=1)
A["max_abs_C_B"] = D[cbcols].abs().max(axis=1)
A["max_abs_C_W"] = D[cwcols].abs().max(axis=1)
for lo in (0, 0.001, 0.005, 0.01):
    A[f"signflip_C_B_gt{lo}"] = (((D[cbcols] > lo).sum(axis=1) > 0) &
                                 ((D[cbcols] < -lo).sum(axis=1) > 0))
    A[f"signflip_C_W_gt{lo}"] = (((D[cwcols] > lo).sum(axis=1) > 0) &
                                 ((D[cwcols] < -lo).sum(axis=1) > 0))
A.to_csv(R / "v8_1/WEIGHTING_CONSISTENT_AUDIT.csv", index=False, encoding="utf-8")
D.to_csv(R / "v8_1/WEIGHTING_CONSISTENT.csv", index=False, encoding="utf-8")

print()
print("=" * 78)
print("C_B / C_W 跨 scheme 极差（新 vs 旧）")
print("=" * 78)
oldA = pd.read_csv(R / "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv")
print(f"  range_C_B  中位数: 新 = {A.range_C_B.median():.5f}   旧 = {oldA.range_C_B.median():.5f}")
print(f"  range_C_W  中位数: 新 = {A.range_C_W.median():.5f}   旧 = {oldA.range_C_W.median():.5f}")
for lo in (0, 0.001, 0.005, 0.01):
    print(f"  signflip C_B >{lo:<6}: 新 = {int(A[f'signflip_C_B_gt{lo}'].sum()):3d}   "
          f"旧 = {int(oldA[f'signflip_C_B_gt{lo}'].sum()):3d}")
for lo in (0, 0.001, 0.005, 0.01):
    print(f"  signflip C_W >{lo:<6}: 新 = {int(A[f'signflip_C_W_gt{lo}'].sum()):3d}   "
          f"旧 = {int(oldA[f'signflip_C_W_gt{lo}'].sum()):3d}")

# ---- C_W 是否变化 ----
print()
print("=" * 78)
print("C_W 是否随 scheme 改变")
print("=" * 78)
for s in SCHEMES:
    v = D[f"C_W_{s}"].dropna()
    print(f"  {s:6s}  C_W 中位数 = {v.median():+.5f}   范围 = [{v.min():+.5f}, {v.max():+.5f}]")
rng_cw = A.range_C_W
print()
print(f"  逐配置 C_W 极差中位数 = {rng_cw.median():.5f}   >0 比例 = {(rng_cw>1e-12).mean():.1%}")

json.dump(dict(
    n_configs=int(len(D)),
    identity_residual_by_scheme=resid,
    identity_holds_all=bool(all(v <= 1e-12 for v in resid.values())),
    range_C_B_median_new=float(A.range_C_B.median()),
    range_C_B_median_old=float(oldA.range_C_B.median()),
    range_C_W_median_new=float(A.range_C_W.median()),
    range_C_W_median_old=float(oldA.range_C_W.median()),
    signflip_C_B_new={str(lo): int(A[f"signflip_C_B_gt{lo}"].sum())
                      for lo in (0, 0.001, 0.005, 0.01)},
    signflip_C_B_old={str(lo): int(oldA[f"signflip_C_B_gt{lo}"].sum())
                      for lo in (0, 0.001, 0.005, 0.01)},
    signflip_C_W_new={str(lo): int(A[f"signflip_C_W_gt{lo}"].sum())
                      for lo in (0, 0.001, 0.005, 0.01)},
    signflip_C_W_old={str(lo): int(oldA[f"signflip_C_W_gt{lo}"].sum())
                      for lo in (0, 0.001, 0.005, 0.01)},
), open(R / "v8_1/WEIGHTING_CONSISTENT.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)
print("\nwrote v8_1/WEIGHTING_CONSISTENT{,_AUDIT}.csv + WEIGHTING_CONSISTENT.json")
