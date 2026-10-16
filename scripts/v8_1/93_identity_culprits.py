# -*- coding: utf-8 -*-
"""scripts/v8_1/93_identity_culprits.py -- 找出加权恒等式残差大的具体配置。"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]


def sw(n, s):
    c = np.asarray(n, float)
    return {"ROW": c, "EQUAL": np.ones_like(c), "SQRT": np.sqrt(c),
            "LOG": np.log1p(c), "CAP30": np.minimum(c, 30.0)}[s]


def metrics(y, p, tgt, s):
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = sw(n_t, s); w = w_t[inv]; W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    sd = np.sqrt(vy * vp)
    r = (np.sum(w * (y - my) * (p - mp)) / W) / sd
    By = np.array([np.average(y[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    C_B = np.sum(w_t * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return float(r), float(C_B), float(C_W)


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
zcan = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y = np.asarray(zcan["y"], float)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
TGT_MAS = mas.pdb_id.astype(str).to_numpy()
S = pd.read_csv(R / "experiments/v8/BENCHMARK_SCOREBOARD.csv", low_memory=False)

cache = {}
res = []
for r in S.itertuples():
    rel = r.prediction_file
    if rel not in cache:
        try:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        except Exception:  # noqa: BLE001
            cache[rel] = None
    zz = cache[rel]
    if zz is None or r.model_name not in zz.files:
        continue
    if r.model_name.lower().endswith("_y") or r.model_name in (
            "y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
            "z_raw", "y_cal", "family", "seed", "masks"):
        continue
    a = np.asarray(zz[r.model_name])
    if a.ndim != 1:
        continue
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        continue
    y = Y if r.join_mode.startswith("canonical") else Y_MAS
    tgt = TGT if r.join_mode.startswith("canonical") else TGT_MAS
    if len(a) != len(y):
        continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        continue
    mx = 0.0; worst = None
    for s in SCHEMES:
        rr, cb, cw = metrics(yy, pp, tt, s)
        e = abs(rr - cb - cw)
        if e > mx:
            mx, worst = e, (s, rr, cb, cw)
    res.append(dict(model=r.model_name, universe=r.universe, join=r.join_mode,
                    n=len(yy), nt=len(set(tt)), max_resid=mx,
                    scheme=worst[0], r=worst[1], C_B=worst[2], C_W=worst[3]))

D = pd.DataFrame(res).sort_values("max_resid", ascending=False)
print(f"评估配置数 = {len(D)}")
print()
print("残差最大的 12 个配置:")
print(D.head(12).to_string(index=False))
print()
print(f"  max_resid <= 1e-12 的配置数 = {int((D.max_resid<=1e-12).sum())} / {len(D)}")
print(f"  残差中位数 = {D.max_resid.median():.3e}")
print()
print("按 universe 分组:")
print(D.groupby("universe").max_resid.agg(["count", "median", "max"]).to_string())
