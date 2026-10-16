# -*- coding: utf-8 -*-
"""scripts/v8/11_selection_v2.py -- V8 PHASE 11：AXIS C v2 区分度检验（修复 §9 的缺陷）。

步骤：
  1. 在 U_CANONICAL 的**全部 160 个**可用模型上重算 v2 效用指标
  2. 比较 v1（Efficiency@10 / NDCG@10 / Hit@10）与 v2 的 **IQR / 极差 / sd**
  3. 判断 v2 是否显著改善区分度
  4. 若改善 -> 用 v2 重检 H3（含分 seed 稳健性）
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M      # noqa: E402
import v8selection2 as S2  # noqa: E402

OUT = R / "experiments" / "v8"
SEL = OUT / "selection"
SEL.mkdir(parents=True, exist_ok=True)

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
           allow_pickle=True)
Y = np.asarray(z["y"], float)

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))]
print(f"U_CANONICAL 注册模型 = {len(uc)}")

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

cache, rows = {}, []
for _, r in uc.iterrows():
    rel = r["prediction_file"]
    mname = str(r["model_name"])
    if rel not in cache:
        cache[rel] = np.load(R / rel, allow_pickle=True)
    zz = cache[rel]
    if mname not in zz.files or is_label(mname):
        continue
    a = np.asarray(zz[mname])
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        continue
    if a.ndim != 1 or len(a) != len(Y):
        continue
    if abs(np.corrcoef(a, Y)[0, 1]) > 0.98:
        continue
    m = np.isfinite(Y) & np.isfinite(a)
    if m.sum() < 100:
        continue
    U = M.prep_universe(Y[m], a[m], TGT[m])
    v2 = S2.axis_selection_v2(U)
    rec = {"model_name": mname, "model_family": r["model_family"],
           "seed": (mname.split("_s")[-1] if "_s" in mname else "NA")}
    for c in ["efficiency@10", "ndcg@10", "hit@10", "macro_spearman",
              "overall_pearson", "C_W", "rmse"]:
        rec[c] = float(r[c]) if c in r.index and pd.notna(r[c]) else np.nan
    rec.update(v2)
    rows.append(rec)

D = pd.DataFrame(rows)
D.to_csv(SEL / "selection_v2_metrics.csv", index=False, encoding="utf-8")
print(f"评价模型数 = {len(D)}")

V1 = ["efficiency@10", "ndcg@10", "hit@10"]
V2 = [c for c in S2.V2_COLS if c in D.columns]
print("\n=== 区分度对比（U_CANONICAL）===")
print(f"{'metric':<26}{'mean':>9}{'sd':>9}{'IQR':>9}{'range':>9}{'sd_norm':>9}")
disc = []
for c in V1 + V2:
    v = D[c].dropna().to_numpy(float)
    if len(v) < 10:
        continue
    q25, q75 = np.percentile(v, [25, 75])
    print(f"{c:<26}{v.mean():>9.4f}{v.std():>9.4f}{q75-q25:>9.4f}"
          f"{v.max()-v.min():>9.4f}{v.std()/max(1e-9,abs(v.mean())):>9.4f}")
    disc.append({"metric": c, "axis": "v1" if c in V1 else "v2",
                 "mean": v.mean(), "sd": v.std(), "iqr": q75 - q25,
                 "range": v.max() - v.min(),
                 "sd_norm": v.std() / max(1e-9, abs(v.mean()))})
DD = pd.DataFrame(disc)
DD.to_csv(SEL / "discrimination_v1_v2.csv", index=False, encoding="utf-8")

print("\n=== 汇总 ===")
for ax, g in DD.groupby("axis"):
    print(f"  {ax}: IQR 中位={g.iqr.median():.4f} | 极差中位={g['range'].median():.4f} "
          f"| sd_norm 中位={g.sd_norm.median():.4f}")
i1 = DD[DD.axis == "v1"].iqr.median()
i2 = DD[DD.axis == "v2"].iqr.median()
r1 = DD[DD.axis == "v1"]["range"].median()
r2 = DD[DD.axis == "v2"]["range"].median()
print(f"\n  IQR 提升 = {i2/max(1e-12,i1):.2f}×   |   极差提升 = {r2/max(1e-12,r1):.2f}×")
best = DD[DD.axis == "v2"].nlargest(3, "iqr")
print("  v2 中 IQR 最大的三个:", list(best.metric))

# ---------- 用 v2 重检 H3 ----------
print("\n" + "=" * 66)
print("H3 重检（v2 效用）")
print("=" * 66)
PRED_CAL = ["overall_pearson", "rmse", "target_mean_pearson"]
PRED_RANK = ["macro_spearman"]
for c in ("macro_kendall", "pairwise_d05"):
    pass
PRED_RANK = [c for c in ["macro_spearman", "C_W"] if c in D.columns]
h3 = []
for ut in V2:
    for lab_ax, cols in (("calibration", PRED_CAL), ("ranking", PRED_RANK)):
        v = [spearmanr(D[c], D[ut]).correlation for c in cols if c in D.columns]
        v = [x for x in v if np.isfinite(x)]
        if not v:
            continue
        h3.append({"utility": ut, "predictor_axis": lab_ax,
                   "median_spearman": float(np.median(v))})
H3 = pd.DataFrame(h3)
piv = H3.pivot_table(index="utility", columns="predictor_axis",
                     values="median_spearman")
piv["ranking_minus_calibration"] = piv.get("ranking", np.nan) - piv.get(
    "calibration", np.nan)
print(piv.round(4).to_string())
piv.to_csv(SEL / "h3_v2_alignment.csv", encoding="utf-8")
n_better = int((piv["ranking_minus_calibration"] > 0).sum())
print(f"\n  ranking 优于 calibration 的效用指标数 = {n_better} / {len(piv)}")

# 分 seed 稳健性（用 IQR 最大的 v2 效用）
ut_best = best.metric.iloc[0]
print(f"\n=== 分 seed 稳健性（效用 = {ut_best}）===")
rob = []
for sd, g in D.groupby("seed"):
    if len(g) < 8:
        continue
    for lab_ax, cols in (("calibration", PRED_CAL), ("ranking", PRED_RANK)):
        v = [spearmanr(g[c], g[ut_best]).correlation for c in cols if c in g.columns]
        v = [x for x in v if np.isfinite(x)]
        if v:
            rob.append({"seed": sd, "n": len(g), "axis": lab_ax,
                        "median_spearman": float(np.median(v))})
    print(f"  seed={sd} (n={len(g)}): "
          + "  ".join(f"{r['axis']}={r['median_spearman']:+.3f}"
                      for r in rob if r["seed"] == sd))
RB = pd.DataFrame(rob)
RB.to_csv(SEL / "h3_v2_seed_robustness.csv", index=False, encoding="utf-8")
if len(RB):
    rp = RB.pivot_table(index="seed", columns="axis", values="median_spearman")
    if "ranking" in rp and "calibration" in rp:
        win = int((rp["ranking"] > rp["calibration"]).sum())
        print(f"\n  ranking 胜出的 seed 数 = {win} / {len(rp)}")
        print(f"  符号一致性（ranking 列）: "
              f"{'稳定' if (rp['ranking'] > 0).all() or (rp['ranking'] < 0).all() else '不稳定'}")

print(f"\n写出 selection_v2_metrics.csv / discrimination_v1_v2.csv / "
      f"h3_v2_alignment.csv / h3_v2_seed_robustness.csv")
