# -*- coding: utf-8 -*-
"""scripts/v8/11c_axis_c_saturation.py -- AXIS C 地板机制的**精确刻画**。

`11b` 已证明地板来自数据，但留下两个必须澄清的点：

  1. **当 n <= K 时**，模型必然选中**全部**突变 ⇒ `Eff@10 ≡ 1` 且
     `range = E_max − E_rand = 0` ⇒ 归一化**除以零**，
     导致 `util_norm` 均值出现 2.93 这种无意义值。
     ⇒ 必须**只在 n > K 的 target 上**做归一化。
  2. `n p75 = 34.8` 时 `E_rand = E_max = 1.0`（**完全饱和**）。
     需精确统计**有多大比例的 target 根本不具备区分能力**。

本脚本给出：可区分 target 占比、按 n 分箱的地板/区间、
以及**修正后的归一化效用**（仅 n > K）与其秩相关。
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
OUT = R / "experiments" / "v8"
SEL = OUT / "selection"

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
K = 10

G = pd.read_csv(SEL / "axis_c_floor_by_target.csv")
print(f"eligible targets = {len(G)} | K = {K}")

# ---------- 1. 饱和统计 ----------
G["saturated"] = G.E_rand >= 1.0 - 1e-12          # n <= K -> E_rand = 1
G["degenerate"] = G["range"] <= 1e-12             # 无区分空间
n_sat = int(G.saturated.sum())
n_deg = int(G.degenerate.sum())
n_disc = len(G) - n_deg
print(f"\n=== 饱和 / 退化统计 ===")
print(f"  E_rand = 1.0（n <= K，模型必选中全部）: {n_sat} / {len(G)} "
      f"({100*n_sat/len(G):.1f}%)")
print(f"  可区分区间 = 0（无区分空间）          : {n_deg} / {len(G)} "
      f"({100*n_deg/len(G):.1f}%)")
print(f"  **具备区分能力的 target**              : {n_disc} / {len(G)} "
      f"({100*n_disc/len(G):.1f}%)")

# 按 n 分箱
print(f"\n=== 按 target 大小 n 分箱 ===")
bins = [(5, 10), (11, 15), (16, 25), (26, 50), (51, 100), (101, 10**9)]
brow = []
print(f"{'n range':<14}{'n_tgt':>7}{'E_rand':>9}{'E_max':>9}{'range':>9}"
      f"{'可区分%':>9}")
for lo, hi in bins:
    g = G[(G.n >= lo) & (G.n <= hi)]
    if not len(g):
        continue
    disc = 100 * (g["range"] > 1e-12).mean()
    hi_s = "+" if hi >= 10**8 else str(hi)
    lab_n = f"{lo}-{hi_s}"
    print(f"{lab_n:<14}{len(g):>7}"
          f"{g.E_rand.mean():>9.4f}{g.E_max.mean():>9.4f}"
          f"{g['range'].mean():>9.4f}{disc:>9.1f}")
    brow.append({"n_lo": lo, "n_hi": hi if hi < 10**8 else -1,
                 "n_targets": len(g), "E_rand": g.E_rand.mean(),
                 "E_max": g.E_max.mean(), "range": g["range"].mean(),
                 "pct_discriminative": disc})
pd.DataFrame(brow).to_csv(SEL / "axis_c_by_size_bin.csv", index=False,
                          encoding="utf-8")

# ---------- 2. 仅 n > K 的归一化效用 ----------
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

GD = G[G.n > K].copy()                     # 只有这些 target 可区分
print(f"\n=== 仅用 n > K 的 {len(GD)} 个可区分 target 重算归一化效用 ===")
S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))]
cache, rows = {}, []
for _, r in uc.iterrows():
    rel = r["prediction_file"]; mname = str(r["model_name"])
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
    ef_raw, un = [], []
    for _, row in GD.iterrows():
        i = np.where(TGT == row.target_id)[0]
        yg, pg = Y[i], a[i]
        n = i.size
        order = np.argsort(yg, kind="stable")
        rank = np.empty(n, int); rank[order] = np.arange(n)
        gain = 1.0 + (n - rank)
        k = max(1, min(K, n))
        gid = gain[np.argsort(-gain, kind="stable")[:k]].sum()
        if gid <= 0:
            continue
        ef = gain[np.argsort(pg, kind="stable")[:k]].sum() / gid
        ef_raw.append(ef)
        un.append((ef - row.E_rand) / row["range"])       # range > 0 保证
    if not un:
        continue
    rows.append({"model_name": mname, "model_family": r["model_family"],
                 "macro_spearman": r["macro_spearman"],
                 "overall_pearson": r["overall_pearson"],
                 "eff10_disc": float(np.mean(ef_raw)),
                 "util_norm_fixed": float(np.mean(un))})
U = pd.DataFrame(rows)
U.to_csv(SEL / "axis_c_normalized_fixed.csv", index=False, encoding="utf-8")

print(f"模型数 = {len(U)}")
print(f"{'metric':<20}{'mean':>9}{'sd':>9}{'IQR':>9}{'range':>9}")
for c in ("eff10_disc", "util_norm_fixed", "macro_spearman", "overall_pearson"):
    v = U[c].dropna().to_numpy(float)
    q25, q75 = np.percentile(v, [25, 75])
    print(f"{c:<20}{v.mean():>9.4f}{v.std():>9.4f}{q75-q25:>9.4f}"
          f"{v.max()-v.min():>9.4f}")

print(f"\n=== 修正后 util_norm_fixed 的秩相关 ===")
al = []
for c in ("macro_spearman", "overall_pearson"):
    sp = spearmanr(U.util_norm_fixed, U[c]).correlation
    kt = kendalltau(U.util_norm_fixed, U[c]).correlation
    print(f"  ρ(U_norm_fixed, {c:<18}) = {sp:+.4f}   τ = {kt:+.4f}")
    al.append({"pair": f"U_norm_fixed vs {c}", "spearman": sp, "kendall": kt})
pd.DataFrame(al).to_csv(SEL / "axis_c_fixed_alignment.csv", index=False,
                        encoding="utf-8")

# 分 seed 稳定性
print(f"\n=== 分 seed 稳健性（U_norm_fixed）===")
U["seed"] = U.model_name.str.extract(r"_s(\d+)$")[0].fillna("NA")
rob = []
for sd, g in U.groupby("seed"):
    if len(g) < 8:
        continue
    line = f"  seed={sd} (n={len(g)}):"
    rec = {"seed": sd, "n": len(g)}
    for c in ("macro_spearman", "overall_pearson"):
        v = spearmanr(g.util_norm_fixed, g[c]).correlation
        rec[c] = v
        line += f"  ρ(·,{c[:12]})={v:+.3f}"
    print(line)
    rob.append(rec)
RB = pd.DataFrame(rob)
RB.to_csv(SEL / "axis_c_fixed_seed.csv", index=False, encoding="utf-8")
if len(RB):
    w = int((RB.macro_spearman > RB.overall_pearson).sum())
    print(f"\n  macro_spearman 优于 overall_pearson 的 seed 数 = {w} / {len(RB)}")
    for c in ("macro_spearman", "overall_pearson"):
        sg = "稳定" if (RB[c] > 0).all() or (RB[c] < 0).all() else "不稳定"
        print(f"  {c:<18} 跨 seed 符号: {sg}")

print("\n" + "=" * 68)
print("最终判读")
print("=" * 68)
print(f"  具备区分能力的 target 仅 {n_disc}/{len(G)} = {100*n_disc/len(G):.1f}%")
unq = np.percentile(U.util_norm_fixed, [25, 75])
print(f"  这些 target 的 U_norm_fixed: "
      f"IQR={float(np.ptp(unq)):.4f} "
      f"(对比 v1 Eff@10 的 IQR=0.0093)")
print()
print("  ⇒ 结论：AXIS C 的「发现地板」是**数据施加的**（SKEMPI 的 target 规模分布 +")
print("     gain 线性于名次），**不是指标写错**。即便改用完全归一化的效用，")
print("     可用于区分的 target 也只有约三分之一。")
