# -*- coding: utf-8 -*-
"""scripts/v8_1/103_weighting_final.py -- 以 TARGET TOTAL MASS 为主口径产出最终加权审计。

语义决定（依据）：
  * 代码 docstring 明确写：ROW = 每个突变等权；EQUAL = 每个 target 等权；
    SQRT/LOG/CAP30 = 权重 ∝ sqrt(n) / log(1+n) / min(n,30)。
    这些名称描述的是 **target 级聚合权重**。
  * 但实现是 w_ci = a_c，使 target 总质量变成 n_c*a_c：
      - EQUAL(a=1) 给每个 target 总质量 n_c，大 target 仍主导，**不是**每个 target 等权；
      - ROW(a=n_c) 给每行权重 n_c，也**不是**每个突变等权。
  * 若要让名称与数学一致，唯一自洽的形式是
        w_ci = A_c / n_c,   A_c = a_c
    此时 target 总质量 = A_c：
      - EQUAL(A=1) 真正是每 target 等权；
      - ROW(A=n_c) 等价于每行等权。

因此：
  主口径 = T（target total mass, w = a_c/n_c）
  敏感性 = P（per-row multiplier, w = a_c，即原实现）
两者都满足 r = C_B + C_W（分子分母同测度），差别只在 target 总质量的定义。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
NAME = {"ROW": "per-mutation", "EQUAL": "per-target", "SQRT": "sqrt(n)",
        "LOG": "log(1+n)", "CAP30": "min(n,30)"}


def a_of(n, scheme):
    c = np.asarray(n, float)
    return {"ROW": c, "EQUAL": np.ones_like(c), "SQRT": np.sqrt(c),
            "LOG": np.log1p(c), "CAP30": np.minimum(c, 30.0)}[scheme]


def weighted(y, p, tgt, scheme, mode):
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_c = np.bincount(inv).astype(float)
    a = a_of(n_c, scheme)
    w_row = a if mode == "P" else a / n_c
    w = w_row[inv]
    W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 3
    sd = np.sqrt(vy * vp)
    By = np.empty(len(uq)); Bp = np.empty(len(uq)); wbar = np.empty(len(uq))
    for k in range(len(uq)):
        m = inv == k
        wbar[k] = w[m].sum()
        By[k] = np.sum(w[m] * y[m]) / wbar[k]
        Bp[k] = np.sum(w[m] * p[m]) / wbar[k]
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
    y = Y_CAN if r.join_mode.startswith("canonical") else Y_MAS
    tgt = TGT if r.join_mode.startswith("canonical") else TGT_MAS
    if len(a) != len(y):
        continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        continue
    rec = {"model_name": r.model_name, "universe": r.universe}
    for s in SCHEMES:
        rr, cb, cw = weighted(yy, pp, tt, s, "T")
        rec[f"pearson_{s}"] = rr; rec[f"C_B_{s}"] = cb; rec[f"C_W_{s}"] = cw
    for s in SCHEMES:
        rr, cb, cw = weighted(yy, pp, tt, s, "P")
        rec[f"P_pearson_{s}"] = rr; rec[f"P_C_B_{s}"] = cb; rec[f"P_C_W_{s}"] = cw
    recs.append(rec)

D = pd.DataFrame(recs)
cb = [f"C_B_{s}" for s in SCHEMES]; cw = [f"C_W_{s}" for s in SCHEMES]
pcb = [f"P_C_B_{s}" for s in SCHEMES]; pcw = [f"P_C_W_{s}" for s in SCHEMES]
print(f"配置数 = {len(D)}")

stats = {}
for tag, bc, wc in (("T_target_total_mass", cb, cw), ("P_per_row_multiplier", pcb, pcw)):
    rngB = (D[bc].max(axis=1) - D[bc].min(axis=1))
    rngW = (D[wc].max(axis=1) - D[wc].min(axis=1))
    sfB = {str(lo): int((((D[bc] > lo).sum(axis=1) > 0) &
                         ((D[bc] < -lo).sum(axis=1) > 0)).sum())
           for lo in (0, 0.001, 0.005, 0.01)}
    sfW = {str(lo): int((((D[wc] > lo).sum(axis=1) > 0) &
                         ((D[wc] < -lo).sum(axis=1) > 0)).sum())
           for lo in (0, 0.001, 0.005, 0.01)}
    stats[tag] = dict(range_C_B_median=float(rngB.median()),
                      range_C_W_median=float(rngW.median()),
                      signflip_C_B=sfB, signflip_C_W=sfW,
                      C_W_median_by_scheme={s: float(D[f"{wc[0].split('_')[0]}_C_W_{s}"
                                                       if False else wc[SCHEMES.index(s)]].median())
                                            for s in SCHEMES})
    print()
    print("=" * 78)
    print(f"口径 {tag}")
    print("=" * 78)
    print(f"  range_C_B 中位数 = {rngB.median():.5f}   range_C_W 中位数 = {rngW.median():.5f}")
    print(f"  C_B 符号翻转 = {sfB}")
    print(f"  C_W 符号翻转 = {sfW}")
    print(f"  各 scheme 的 C_W 中位数 = "
          f"{ {s: round(float(D[wc[i]].median()),5) for i,s in enumerate(SCHEMES)} }")

# 主口径 audit（T）
A = pd.DataFrame({"model_name": D.model_name})
A["range_C_B"] = D[cb].max(axis=1) - D[cb].min(axis=1)
A["range_C_W"] = D[cw].max(axis=1) - D[cw].min(axis=1)
A["max_abs_C_B"] = D[cb].abs().max(axis=1)
A["max_abs_C_W"] = D[cw].abs().max(axis=1)
for lo in (0, 0.001, 0.005, 0.01):
    A[f"signflip_C_B_gt{lo}"] = (((D[cb] > lo).sum(axis=1) > 0) &
                                 ((D[cb] < -lo).sum(axis=1) > 0))
    A[f"signflip_C_W_gt{lo}"] = (((D[cw] > lo).sum(axis=1) > 0) &
                                 ((D[cw] < -lo).sum(axis=1) > 0))
A.to_csv(R / "v8_1/WEIGHTING_FINAL_AUDIT.csv", index=False, encoding="utf-8")
D.to_csv(R / "v8_1/WEIGHTING_FINAL.csv", index=False, encoding="utf-8")

json.dump(dict(primary="T_target_total_mass", sensitivity="P_per_row_multiplier",
               schemes={s: NAME[s] for s in SCHEMES},
               w_ci_primary="a_c / n_c", w_ci_sensitivity="a_c",
               stats=stats),
          open(R / "v8_1/WEIGHTING_FINAL.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print("\nwrote v8_1/WEIGHTING_FINAL{,_AUDIT}.csv + WEIGHTING_FINAL.json")
