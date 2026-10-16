# -*- coding: utf-8 -*-
"""scripts/v8_1/90_weighted_full_verify.py -- 对全部 109 个配置复现加权分解。

严格照抄 scripts/v8/09_weighting_audit.py 的实现：
  * 预测按 BENCHMARK_SCOREBOARD.csv 的 prediction_file 载入
  * U_CANONICAL 用 canonical y/tgt，U_MAS 用 matched_analysis_set 的 ddG/pdb_id
  * 以 |corr(a, y)| > 0.98 作为「疑似标签泄漏」跳过（与代码一致）

输出：
  * 与 weighting_sensitivity.csv 的逐格最大差
  * 每个 scheme 的 r = C_B + C_W 残差
  * 每 scheme 的 C_W 分布（回答 Note 5 的自相矛盾）
  * C_B / C_W 的跨 scheme 极差中位数（应复现 audit 的 0.0183 / 0.0624）
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


def weighted_metrics(y, p, tgt, scheme):
    """与 09_weighting_audit.weighted_metrics 等价（只保留 (r, C_B, C_W)）。"""
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = scheme_weights(n_t, scheme)
    w = w_t[inv]
    W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 3
    r = (np.sum(w * (y - my) * (p - mp)) / W) / np.sqrt(vy * vp)
    By = np.array([np.average(y[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    sd = np.sqrt(vy) * np.sqrt(vp)
    C_B = np.sum(w_t * (By - my) * (Bp - mp)) / (W * sd)
    Wy = y - By[inv]; Wp = p - Bp[inv]
    C_W = np.sum(w * Wy * Wp) / (W * sd)
    return float(r), float(C_B), float(C_W)


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
CAN = lab[["pair_row", "target_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id"})
ROWS = CAN.row_id.to_numpy(int); TGT = CAN.target_id.to_numpy(str)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
MAS_ROWS = mas.row.to_numpy(int); MAS_TGT = mas.pdb_id.astype(str).to_numpy()
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)

S = pd.read_csv(R / "experiments/v8/BENCHMARK_SCOREBOARD.csv", low_memory=False)
disk = pd.read_csv(R / "experiments/v8/weighting/weighting_sensitivity.csv",
                   low_memory=False)
disk = disk.drop_duplicates("model_name").set_index("model_name")

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed", "masks"}
cache, recs = {}, []
loaded_ok, skipped = 0, []
for r in S.itertuples():
    rel = r.prediction_file
    if rel not in cache:
        try:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        except Exception:  # noqa: BLE001
            cache[rel] = None
    zz = cache[rel]
    if zz is None or r.model_name not in zz.files:
        skipped.append((r.model_name, "no array")); continue
    if str(r.model_name).lower().endswith("_y"):
        skipped.append((r.model_name, "label key")); continue
    a = np.asarray(zz[r.model_name])
    if a.ndim != 1:
        skipped.append((r.model_name, "ndim")); continue
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        skipped.append((r.model_name, "dtype")); continue

    if r.universe == "U_CANONICAL":
        z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz",
                    allow_pickle=True)
        y, tgt, base = np.asarray(z["y"], float), TGT, ROWS
    else:
        y, tgt, base = Y_MAS, MAS_TGT, MAS_ROWS
    if r.universe == "U_CANONICAL" and len(a) != len(base):
        skipped.append((r.model_name, "len mismatch")); continue
    if len(a) == len(y) and np.isfinite(a).all():
        if abs(np.corrcoef(a, y)[0, 1]) > 0.98:
            skipped.append((r.model_name, "leak-guard")); continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        skipped.append((r.model_name, "few targets")); continue
    row = {"model_name": r.model_name, "universe": r.universe}
    for s in SCHEMES:
        rr, cb, cw = weighted_metrics(yy, pp, tt, s)
        row[f"pearson_{s}"] = rr; row[f"C_B_{s}"] = cb; row[f"C_W_{s}"] = cw
    recs.append(row); loaded_ok += 1

D = pd.DataFrame(recs)
print(f"重算配置数 = {len(D)}   (跳过 {len(skipped)})")
print("跳过原因分布:", pd.Series([x[1] for x in skipped]).value_counts().to_dict())
print()

cbcols = [f"C_B_{s}" for s in SCHEMES]; cwcols = [f"C_W_{s}" for s in SCHEMES]

print("=" * 78)
print("1) 显式公式 vs 落盘 weighting_sensitivity.csv")
print("=" * 78)
common = [n for n in D.model_name if n in disk.index]
maxd = 0.0; ndiff = 0; ncells = 0
for n in common:
    mine = D[D.model_name == n].iloc[0]; ref = disk.loc[n]
    for s in SCHEMES:
        for key in (f"pearson_{s}", f"C_B_{s}", f"C_W_{s}"):
            a, b = mine[key], ref[key]
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            ncells += 1
            d = abs(a - b)
            if d > maxd: maxd = d
            if d > 1e-9: ndiff += 1
print(f"  共同配置 = {len(common)}   比对格数 = {ncells}")
print(f"  最大绝对差 = {maxd:.3e}   超过 1e-9 的格数 = {ndiff}")
print("  => 公式与代码一致" if maxd < 1e-9 else "  => 仍有差异，需继续核对")

print()
print("=" * 78)
print("2) 每个 scheme 下 r = C_B + C_W 的残差（逐配置）")
print("=" * 78)
resid = {}
for s in SCHEMES:
    e = (D[f"pearson_{s}"] - D[f"C_B_{s}"] - D[f"C_W_{s}"]).abs()
    resid[s] = float(np.nanmax(e))
    print(f"  {s:6s}  max|r - C_B - C_W| = {resid[s]:.3e}   "
          f"全部 <= 1e-12: {bool((e.dropna() <= 1e-12).all())}")

print()
print("=" * 78)
print("3) C_W 是否跨 scheme 改变？（Note 5 的关键）")
print("=" * 78)
for s in SCHEMES:
    v = D[f"C_W_{s}"].dropna()
    print(f"  {s:6s}  C_W 中位数 = {v.median():+.5f}   范围 = [{v.min():+.5f}, {v.max():+.5f}]")
rng_cw = (D[cwcols].max(axis=1) - D[cwcols].min(axis=1))
rng_cb = (D[cbcols].max(axis=1) - D[cbcols].min(axis=1))
print()
print(f"  逐配置 C_W 跨 scheme 极差: 中位数 = {rng_cw.median():.5f}  "
      f"(audit 记录 0.06245)")
print(f"  逐配置 C_B 跨 scheme 极差: 中位数 = {rng_cb.median():.5f}  "
      f"(audit 记录 0.01828)")
print(f"  C_W 极差 > 0 的配置占比 = {(rng_cw > 1e-12).mean():.1%}")

out = dict(
    n_recomputed=int(len(D)), n_compared=len(common), n_cells=ncells,
    max_abs_diff_vs_disk=float(maxd), n_cells_over_tol=int(ndiff),
    identity_residual_by_scheme=resid,
    c_w_range_median=float(rng_cw.median()),
    c_b_range_median=float(rng_cb.median()),
    c_w_varies_in_all=(bool((rng_cw > 1e-12).all())),
    formula=("a_c = scheme_weights(n_c); w_ci = a_c (constant within target); "
             "W = sum_ci w_ci; m = sum w x / W; sigma^2 = sum w (x-m)^2 / W; "
             "C_B = sum_c a_c (B_yc - m_y)(B_pc - m_p) / (W sigma_y sigma_p); "
             "C_W = sum_ci a_c (y_ci - B_yc)(p_ci - B_pc) / (W sigma_y sigma_p)"),
    skipped_reasons=pd.Series([x[1] for x in skipped]).value_counts().to_dict(),
)
Path("v8_1/WEIGHTED_DECOMPOSITION.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print("\nwrote v8_1/WEIGHTED_DECOMPOSITION.json")
