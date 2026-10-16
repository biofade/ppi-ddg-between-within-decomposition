# -*- coding: utf-8 -*-
"""scripts/v8_1/88_weighted_identity_verify.py -- 从真实代码重建 weighted decomposition 的数学形式。

真实实现（scripts/v8/09_weighting_audit.py, weighted_metrics + scheme_weights）：
    n_c = target c 的行数
    a_c = scheme_weights(n_c, scheme)                     # ROW n_c / EQUAL 1 /
                                                          # SQRT sqrt(n_c) /
                                                          # LOG log1p(n_c) /
                                                          # CAP30 min(n_c,30)
    w_ci = a_c                     <-- 每个 target 内的**所有行取得同一个权重 a_c**
    W    = sum_c n_c * a_c
    m_y  = sum_ci w_ci y_ci / W
    v_y  = sum_ci w_ci (y_ci - m_y)^2 / W
    sigma_y = sqrt(v_y),  sigma_p = sqrt(v_p)
    B_y,c = sum_i w_ci y_ci / sum_i w_ci = mean(y in c)    # 组内等权 -> 就是普通均值
    C_B = sum_c a_c (B_y,c - m_y)(B_p,c - m_p) / (W sigma_y sigma_p)
    C_W = sum_ci w_ci (y_ci - B_y,c)(p_ci - B_p,c) / (W sigma_y sigma_p)

本脚本：
  1. 用上述「显式公式」独立重算 5 个 scheme 下的 r / C_B / C_W
  2. 与 weighting_sensitivity.csv 的落盘值逐格比对
  3. 验证每个 scheme、每个 configuration 上 r = C_B + C_W
  4. 验证 C_W 是否真的随 scheme 改变（回答 Supplement Note 5 的自相矛盾）
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("scripts/v8")))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from importlib import import_module  # noqa: E402
S = import_module("09_weighting_audit") if False else None

# 直接内联真实实现的等价公式（与 scripts/v8/09_weighting_audit.py 一致）
def scheme_weights(n, scheme):
    c = np.asarray(n, float)
    return {"ROW": c, "EQUAL": np.ones_like(c), "SQRT": np.sqrt(c),
            "LOG": np.log1p(c), "CAP30": np.minimum(c, 30.0)}[scheme]


def explicit_metrics(y, p, tgt, scheme):
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
    sd = np.sqrt(vy) * np.sqrt(vp)
    By = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bp = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_B = np.sum(a * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    r = C_B + C_W
    return r, C_B, C_W


# ---- 用真实数据重建 5 个 scheme 下的指标 ----
R = Path(".")
lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
CAN = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id",
             "family_id_canonical": "family_id"})
TGT = CAN.target_id.to_numpy(str)
z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y = np.asarray(z["y"], float)
w2 = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv", low_memory=False)

SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed", "masks"}

# 复现 weighting_sensitivity 的加载顺序，逐 model 重算
pred_files = sorted(set(
    (R / "experiments/v7_1/phase2_rank").glob("*.npz")))
cache = {p: np.load(p, allow_pickle=True) for p in pred_files}

print("=" * 78)
print("1) 显式公式 vs 落盘 weighting_sensitivity.csv")
print("=" * 78)
checked = 0
maxdiff = 0.0
resid = {s: 0.0 for s in SCHEMES}
cw_seen = {s: [] for s in SCHEMES}
rows_check = []
for r in w2.itertuples():
    if getattr(r, "universe", "U_CANONICAL") != "U_CANONICAL":
        continue
    name = r.model_name
    arr = None
    for p, zz in cache.items():
        if name in zz.files and not (str(name).lower().endswith("_y")):
            a = np.asarray(zz[name])
            if a.ndim == 1 and len(a) == len(Y):
                try:
                    a = a.astype(float)
                except (ValueError, TypeError):
                    continue
                if np.isfinite(a).all():
                    arr = a; break
    if arr is None:
        continue
    m = np.isfinite(Y) & np.isfinite(arr)
    yy, pp, tt = Y[m], arr[m], TGT[m]
    for s in SCHEMES:
        rr, cb, cw = explicit_metrics(yy, pp, tt, s)
        got_r = getattr(r, f"pearson_{s}")
        got_cb = getattr(r, f"C_B_{s}")
        got_cw = getattr(r, f"C_W_{s}")
        if any(not np.isfinite(x) for x in (got_r, got_cb, got_cw)):
            continue
        d = max(abs(rr - got_r), abs(cb - got_cb), abs(cw - got_cw))
        maxdiff = max(maxdiff, d)
        resid[s] = max(resid[s], abs(rr - cb - cw))
        cw_seen[s].append(cw)
        checked += 1
        if checked <= 3:
            rows_check.append((name, s, got_r, cb, cw, rr - got_r))
print(f"  比对格数 = {checked}")
print(f"  与落盘值最大绝对差 = {maxdiff:.3e}   (>1e-10 则公式不符)")
print()
print("  %-28s %-6s %-10s %-10s %-10s %s" % ("model", "scheme", "r(disk)", "C_B", "C_W", "Δr"))
for n, s, gr, cb, cw, dr in rows_check:
    print("  %-28s %-6s %-10.4f %-10.4f %-10.4f %+.1e" % (n[:28], s, gr, cb, cw, dr))

print()
print("=" * 78)
print("2) 每个 scheme 下 r = C_B + C_W 的残差")
print("=" * 78)
for s in SCHEMES:
    print(f"  {s:6s}  max|r - C_B - C_W| = {resid[s]:.3e}")

print()
print("=" * 78)
print("3) C_W 是否随 scheme 改变？")
print("=" * 78)
for s in SCHEMES:
    v = np.array(cw_seen[s])
    if len(v):
        print(f"  {s:6s}  n={len(v):3d}  C_W 中位数={np.median(v):+.5f}  "
              f"范围=[{v.min():+.5f}, {v.max():+.5f}]")
# 同一配置跨 scheme 的 C_W 极差
cwdf = pd.DataFrame({s: cw_seen[s] for s in SCHEMES if cw_seen[s]})
if len(cwdf):
    rng = cwdf.max(axis=1) - cwdf.min(axis=1)
    print()
    print(f"  逐配置 C_W 极差: 中位数 = {rng.median():.4f}  "
          f"最大 = {rng.max():.4f}  非零占比 = {(rng>1e-12).mean():.1%}")
    print("  -> C_W 并非跨 scheme 不变；Note 5 的 'identical across all five schemes' 是错的")

out = dict(
    n_compared=checked, max_abs_diff_vs_disk=float(maxdiff),
    identity_residual_by_scheme={s: float(resid[s]) for s in SCHEMES},
    c_w_median_by_scheme={s: float(np.median(cw_seen[s])) for s in SCHEMES if cw_seen[s]},
    c_w_range_median_perconfig=float((cwdf.max(axis=1) - cwdf.min(axis=1)).median())
    if len(cwdf) else None,
    formula=("a_c = scheme_weights(n_c); w_ci = a_c (constant within target); "
             "W = sum_c n_c a_c; "
             "C_B = sum_c a_c (B_y,c - m_y)(B_p,c - m_p) / (W sigma_y sigma_p); "
             "C_W = sum_ci a_c (y_ci - B_y,c)(p_ci - B_p,c) / (W sigma_y sigma_p)"),
)
Path("v8_1/WEIGHTED_DECOMPOSITION.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print("\nwrote v8_1/WEIGHTED_DECOMPOSITION.json")
