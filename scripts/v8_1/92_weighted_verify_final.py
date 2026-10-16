# -*- coding: utf-8 -*-
"""scripts/v8_1/92_weighted_verify_final.py -- 对全部 109 个配置严格复现加权分解并验证恒等式。

按 BENCHMARK_SCOREBOARD.join_mode 正确对齐：
  canonical_rowid / canonical_align -> canonical y/target（4,473 行）
  mas_align                        -> matched_analysis_set 行对齐（7,523 行）
预测一律取自 SCOREBOARD 的 prediction_file 列（数据集自带的连接方式）。
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


def code_metrics(y, p, tgt, scheme):
    """与 scripts/v8/09_weighting_audit.weighted_metrics 完全一致。"""
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = scheme_weights(n_t, scheme); w = w_t[inv]; W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 3
    sd = np.sqrt(vy * vp)
    r = (np.sum(w * (y - my) * (p - mp)) / W) / sd
    By = np.array([np.average(y[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    C_B = np.sum(w_t * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return float(r), float(C_B), float(C_W)


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
ROWS = lab.pair_row.to_numpy(int); TGT = lab.target_id_canonical.to_numpy(str)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
TGT_MAS = mas.pdb_id.astype(str).to_numpy()
zcan = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y_CAN_FULL = np.asarray(zcan["y"], float)

S = pd.read_csv(R / "experiments/v8/BENCHMARK_SCOREBOARD.csv", low_memory=False)
disk = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv",
                   low_memory=False).drop_duplicates("model_name").set_index("model_name")

cache, recs, skipped = {}, [], []
LABELS = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
          "z_raw", "y_cal", "family", "seed", "masks"}
for r in S.itertuples():
    rel = r.prediction_file
    if rel not in cache:
        try:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        except Exception:  # noqa: BLE001
            cache[rel] = None
    zz = cache[rel]
    if zz is None or r.model_name not in zz.files or r.model_name in LABELS \
            or str(r.model_name).lower().endswith("_y"):
        skipped.append((r.model_name, "unavailable")); continue
    a = np.asarray(zz[r.model_name])
    if a.ndim != 1:
        skipped.append((r.model_name, "ndim")); continue
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        skipped.append((r.model_name, "dtype")); continue

    if r.join_mode.startswith("canonical"):
        y, tgt = Y_CAN_FULL, TGT
        if len(a) != len(y):
            skipped.append((r.model_name, "len vs canonical")); continue
    else:
        y, tgt = Y_MAS, TGT_MAS
        if len(a) != len(y):
            skipped.append((r.model_name, "len vs mas")); continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        skipped.append((r.model_name, "few targets")); continue
    row = {"model_name": r.model_name, "universe": r.universe, "join_mode": r.join_mode}
    for s in SCHEMES:
        rr, cb, cw = code_metrics(yy, pp, tt, s)
        row[f"pearson_{s}"] = rr; row[f"C_B_{s}"] = cb; row[f"C_W_{s}"] = cw
    recs.append(row)

D = pd.DataFrame(recs)
cbcols = [f"C_B_{s}" for s in SCHEMES]; cwcols = [f"C_W_{s}" for s in SCHEMES]
print(f"重算配置数 = {len(D)}   (跳过 {len(skipped)})")
if skipped:
    print("跳过原因:", pd.Series([x[1] for x in skipped]).value_counts().to_dict())
print()

print("=" * 80)
print("1) 公式 vs 落盘 weighting_sensitivity.csv")
print("=" * 80)
common = [n for n in D.model_name if n in disk.index]
maxd = 0.0; over = 0; ncells = 0
for n in common:
    mine = D[D.model_name == n].iloc[0]; ref = disk.loc[n]
    for s in SCHEMES:
        for key in (f"pearson_{s}", f"C_B_{s}", f"C_W_{s}"):
            va, vb = mine[key], ref[key]
            if not (np.isfinite(va) and np.isfinite(vb)):
                continue
            ncells += 1
            d = abs(va - vb)
            maxd = max(maxd, d)
            if d > 1e-9:
                over += 1
print(f"  共同配置 = {len(common)}   比对格数 = {ncells}")
print(f"  最大绝对差 = {maxd:.3e}   超 1e-9 的格数 = {over}")
print("  => 公式与代码逐格一致" if maxd < 1e-9 else "  => 仍有差异")

print()
print("=" * 80)
print("2) 逐配置恒等式 r = C_B + C_W")
print("=" * 80)
resid = {}
for s in SCHEMES:
    e = (D[f"pearson_{s}"] - D[f"C_B_{s}"] - D[f"C_W_{s}"]).abs().dropna()
    resid[s] = float(e.max())
    print(f"  {s:6s}  max = {e.max():.3e}   n={len(e)}   全部<=1e-12: {bool((e<=1e-12).all())}")

print()
print("=" * 80)
print("3) C_W 是否跨 scheme 改变")
print("=" * 80)
for s in SCHEMES:
    v = D[f"C_W_{s}"].dropna()
    print(f"  {s:6s}  C_W 中位数 = {v.median():+.5f}")
rng_cw = (D[cwcols].max(axis=1) - D[cwcols].min(axis=1))
rng_cb = (D[cbcols].max(axis=1) - D[cbcols].min(axis=1))
print()
print(f"  逐配置 C_W 极差中位数 = {rng_cw.median():.5f}   (audit: 0.06245)")
print(f"  逐配置 C_B 极差中位数 = {rng_cb.median():.5f}   (audit: 0.01828)")
print(f"  C_W 极差 > 0 的比例 = {(rng_cw>1e-12).mean():.1%}")

json.dump(dict(
    n_recomputed=int(len(D)), n_compared=len(common), n_cells=ncells,
    max_abs_diff_vs_disk=float(maxd), n_cells_over_tol=int(over),
    identity_residual_by_scheme=resid,
    c_w_range_median=float(rng_cw.median()), c_b_range_median=float(rng_cb.median()),
    c_w_varies_all=bool((rng_cw > 1e-12).all()),
), open("v8_1/WEIGHTED_DECOMPOSITION.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)
print("\nwrote v8_1/WEIGHTED_DECOMPOSITION.json")
