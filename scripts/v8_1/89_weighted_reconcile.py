# -*- coding: utf-8 -*-
"""scripts/v8_1/89_weighted_reconcile.py -- 把显式公式与两个落盘文件逐配置对齐。

目的：确定 Note 5 应该写什么公式，并确认 C_B / C_W 的真实量级。
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


def explicit(y, p, tgt, scheme):
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv).astype(float)
    a = scheme_weights(n_t, scheme)
    w = a[inv]
    W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    sd = np.sqrt(vy) * np.sqrt(vp)
    By = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bp = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_B = np.sum(a * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return C_B + C_W, C_B, C_W


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
CAN = lab[["pair_row", "target_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id"})
TGT = CAN.target_id.to_numpy(str)
z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y = np.asarray(z["y"], float)
print(f"canonical: y={len(Y)} rows, {len(set(TGT))} targets")
print(f"z_predictions 可用数组: {[k for k in z.files if not k.endswith('_y')][:12]}")
print()

disk = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv",
                   low_memory=False)
audit = pd.read_csv(R / "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv", low_memory=False)
a_map = audit.set_index("model_name")

cache = {}
print("%-26s %-6s %10s %10s %10s %10s" %
      ("model", "scheme", "r(mine)", "C_B(mine)", "C_W(mine)", "C_W(disk)"))
n_cmp = 0
maxd = 0.0
for r in disk[disk.universe == "U_CANONICAL"].itertuples():
    name = r.model_name
    if name not in z.files:
        continue
    a = np.asarray(z[name])
    if a.ndim != 1 or len(a) != len(Y):
        continue
    a = a.astype(float)
    m = np.isfinite(Y) & np.isfinite(a)
    yy, pp, tt = Y[m], a[m], TGT[m]
    for s in SCHEMES:
        rr, cb, cw = explicit(yy, pp, tt, s)
        dw = getattr(r, f"C_W_{s}"); dr = getattr(r, f"pearson_{s}")
        d = max(abs(rr - dr), abs(cw - dw))
        maxd = max(maxd, d); n_cmp += 1
        if n_cmp <= 6:
            print("%-26s %-6s %10.4f %10.4f %10.4f %10.4f" %
                  (name[:26], s, rr, cb, cw, dw))
print()
print(f"  与 disk 的 C_W / r 最大差 = {maxd:.3e}   (比对 {n_cmp} 格)")

# 落盘 C_B 量级
cbcols = [f"C_B_{s}" for s in SCHEMES]
cwcols = [f"C_W_{s}" for s in SCHEMES]
uc = disk[disk.universe == "U_CANONICAL"]
print()
print("=== 落盘 weighting_sensitivity.csv (U_CANONICAL, n=%d) ===" % len(uc))
print("  |C_B| 最大值           = %.3e" % uc[cbcols].abs().max().max())
print("  C_B 跨 scheme 极差中位 = %.3e" % (uc[cbcols].max(axis=1) - uc[cbcols].min(axis=1)).median())
print("  C_W 跨 scheme 极差中位 = %.5f" % (uc[cwcols].max(axis=1) - uc[cwcols].min(axis=1)).median())
um = disk[disk.universe == "U_MAS"]
if len(um):
    print("=== 落盘 (U_MAS, n=%d) ===" % len(um))
    print("  |C_B| 最大值           = %.4f" % um[cbcols].abs().max().max())
    print("  C_B 跨 scheme 极差中位 = %.5f" % (um[cbcols].max(axis=1) - um[cbcols].min(axis=1)).median())
    print("  C_W 跨 scheme 极差中位 = %.5f" % (um[cwcols].max(axis=1) - um[cwcols].min(axis=1)).median())
print()
print("=== audit (WEIGHTING_EFFECT_SIZE_AUDIT, n=%d) ===" % len(audit))
print("  range_C_B  median = %.6f" % audit.range_C_B.median())
print("  range_C_W  median = %.6f" % audit.range_C_W.median())
