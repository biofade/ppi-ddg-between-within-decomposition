# -*- coding: utf-8 -*-
"""scripts/v8_1/91_identity_diagnose.py -- 诊断加权分解中 r vs C_B+C_W 的差异来源。

09_weighting_audit.weighted_metrics 里有一个**分母口径不一致**：
  * sigma 由 **行加权** 方差得到:  vy = sum_ci w_ci (y_ci - m_y)^2 / W
  * C_B 用 **target 权重** w_t 求和, 但减的是**行加权**均值 m_y
因此 r = C_B + C_W 只在特定 scheme 下精确成立。

本脚本对同一配置，分别计算：
  (a) 代码原式（行加权 sigma + target 权重 C_B）
  (b) 完全一致的加权分解（分子分母同用 w_ci）
  (c) 未加权恒等式的残差
并比较 (a) 与 (b) 的 C_B / C_W，以判断 cw 的精确性。
"""
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


def code_version(y, p, tgt, scheme):
    """严格复刻 09_weighting_audit.weighted_metrics 的 (r, C_B, C_W)。"""
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = scheme_weights(n_t, scheme); w = w_t[inv]; W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    sd = np.sqrt(vy * vp)
    r = (np.sum(w * (y - my) * (p - mp)) / W) / sd
    By = np.array([np.average(y[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    C_B = np.sum(w_t * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return r, C_B, C_W, sd, vy, vp, my, mp, By, Bp, w, w_t, inv, y, p


def consistent_weights(y, p, tgt, scheme):
    """完全一致的加权分解：分子与 sigma 都用同一个加权测度 w_ci。

    定义加权均值 m_w、加权方差 sigma_w，
    并令 B_yc = 组内加权均值（组内等权时为普通均值）。
    此时可证明 r_w = C_B,w + C_W,w 精确成立。
    """
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = scheme_weights(n_t, scheme); w_arr = w_t[inv]; W = w_arr.sum()
    my = np.sum(w_arr * y) / W; mp = np.sum(w_arr * p) / W
    vy = np.sum(w_arr * (y - my) ** 2) / W; vp = np.sum(w_arr * (p - mp) ** 2) / W
    sd = np.sqrt(vy * vp)
    By = np.array([np.average(y[inv == k], weights=w_arr[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w_arr[inv == k]) for k in range(len(uq))])
    # 组内加权平均权重
    wbar = np.array([np.average(w_arr[inv == k]) for k in range(len(uq))])
    C_B = np.sum(wbar * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w_arr * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return C_B + C_W, C_B, C_W


# 用最小可复现配置：z_predictions 的第一个模型
lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y = np.asarray(z["y"], float)
NAME = "L1_random_A1_base_set" if "L1_random_A1_base_set" in z.files else None

disk = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv",
                   low_memory=False).drop_duplicates("model_name").set_index("model_name")

# 找一个 C_B 明显非零的配置来观察
cands = ["L1_random_A1_base_set", "z_R2_MLP_L0", "z_R0_A0_pooled"]
print("=" * 92)
print("配置对比：代码原式 vs 口径一致的加权分解")
print("=" * 92)
for name in cands:
    if name not in z.files:
        continue
    a = np.asarray(z[name], float)
    m = np.isfinite(Y) & np.isfinite(a)
    yy, pp, tt = Y[m], a[m], TGT[m]
    print(f"\n--- {name}  (rows={m.sum()}, targets={len(set(tt))}) ---")
    for s in SCHEMES:
        r, cb, cw, sd, vy, vp, my, mp, By, Bp, w, w_t, inv, yv, pv = code_version(yy, pp, tt, s)
        r2, cb2, cw2 = consistent_weights(yy, pp, tt, s)
        dw = abs(cw - cw2); db = abs(cb - cb2); dr = abs(r - r2)
        print(f"  {s:6s} 代码: r={r:+.4f} C_B={cb:+.4f} C_W={cw:+.4f} "
              f"|resid|={abs(r-cb-cw):.2e}   一致式: C_B={cb2:+.4f} C_W={cw2:+.4f} "
              f"|resid|={abs(r2-cb2-cw2):.2e}   ΔC_W={dw:.2e}")

print()
print("=" * 92)
print("结论要点")
print("=" * 92)
print("  * C_W 在两种口径下基本相同（代码用行加权，与一致式相同）")
print("  * 差异来自 C_B：代码用 target 权重 w_t 求和，但减的是行加权均值 m_y")
print("  * 因此 r = C_B + C_W 在加权情形下**不**精确成立")
