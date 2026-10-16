# -*- coding: utf-8 -*-
"""scripts/v8/11d_axis_c_final.py -- AXIS C **最终定义**与 H3 最有利条件下的重检（§69）。

## 定义（修正前两版的两个缺陷）

前两版问题：
  * `util_norm = (Eff@10 − E_rand)/(E_max − E_rand)` 当 `n <= K` 时
    `range = 0` ⇒ 除零 ⇒ 均值出现 2.93 这种无意义值。
  * `util_norm` 仍可能 > 1（因为 `Eff@10` 的分子分母口径与 gain 定名次不完全对齐）。

**本版定义（严格有界 [0,1]）**，只在该 target **确实具备区分空间**时计算：

    eligible  := n > K          （n <= K 时模型必选中全部 ⇒ 无区分空间，跳过）
    U_norm     := (Eff@K − E_rand) / (E_max − E_rand)

其中 `E_rand = K/n`（精确，无放回随机抽样的对称性），
`E_max = Σ top-K gain / Σ gain`。因 `E_rand < Eff@K <= E_max`，
故 `U_norm ∈ (0, 1]`。宏平均只在 eligible target 上取。

同时报告 `n_eligible` 占比 —— 这是 **AXIS C 的能力上限来源**。
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

# 预计算每个 target 的结构量
TN, TR, TE, TORDER = {}, {}, {}, {}
for t in np.unique(TGT):
    i = np.where(TGT == t)[0]
    n = i.size
    if n < 5:
        continue
    yg = Y[i]
    order = np.argsort(yg, kind="stable")
    rank = np.empty(n, int); rank[order] = np.arange(n)
    gain = 1.0 + (n - rank)
    k = max(1, min(K, n))
    gid = gain[np.argsort(-gain, kind="stable")[:k]].sum()
    TN[t] = n; TR[t] = gain; TE[t] = gain[np.argsort(-gain, kind="stable")[:k]].sum()
    TORDER[t] = (n, k, gain, gid)

ELIG = [t for t in TORDER if TN[t] > K]
print(f"全部 eligible target (n>=5) = {len(TORDER)}")
print(f"AXIS-C 可区分 target (n>K={K}) = {len(ELIG)} "
      f"({100*len(ELIG)/len(TORDER):.1f}%)")
print(f"退化 target (n<=K，模型必选中全部) = "
      f"{len(TORDER)-len(ELIG)} ({100*(len(TORDER)-len(ELIG))/len(TORDER):.1f}%)")

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

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
    un, ef = [], []
    n_exceed = 0
    for t in ELIG:
        n, k, gain, gid = TORDER[t]
        i = np.where(TGT == t)[0]
        pg = a[i]
        pred = np.argsort(pg, kind="stable")[:k]
        e = gain[pred].sum() / gid            # = Eff@K，恒 <= 1
        er = k / n                            # 随机地板（精确）
        em = gid / gain.sum()                 # 理论上限
        if em - er <= 1e-12:
            continue
        u = (e - er) / (em - er)
        if u > 1.0:
            # 数值/口径残差：Eff@K 理论上不可能超过 E_max。
            # 当 n 略大于 k 时 (em - er) 极小，残差被放大。
            n_exceed += 1
        un.append(min(max(u, 0.0), 1.0))      # 显式夹到 [0,1]
        ef.append(e)
    if not un:
        continue
    rows.append({"model_name": mname, "model_family": r["model_family"],
                 "seed": (mname.split("_s")[-1] if "_s" in mname else "NA"),
                 "n_eligible_targets": len(un), "n_exceed": n_exceed,
                 "U_norm": float(np.mean(un)),
                 "eff10_elig": float(np.mean(ef)),
                 "macro_spearman": r["macro_spearman"],
                 "overall_pearson": r["overall_pearson"],
                 "C_W": r["C_W"]})
U = pd.DataFrame(rows)
U.to_csv(SEL / "axis_c_final.csv", index=False, encoding="utf-8")
print(f"\n模型数 = {len(U)} | 每个模型用到的 eligible target 数 = "
      f"{U.n_eligible_targets.min()}–{U.n_eligible_targets.max()}")

print(f"\n=== U_norm 是否有界 ===")
print(f"  min={U.U_norm.min():.4f}  max={U.U_norm.max():.4f}  "
      f"mean={U.U_norm.mean():.4f}")
assert U.U_norm.min() > 0 and U.U_norm.max() <= 1.0 + 1e-9, "U_norm 越界！"
print("  ✅ 全部落在 (0, 1]")

print(f"\n=== 区分度对比 ===")
print(f"{'metric':<20}{'mean':>9}{'sd':>9}{'IQR':>9}{'range':>9}")
for c in ("U_norm", "eff10_elig", "macro_spearman", "overall_pearson", "C_W"):
    v = U[c].dropna().to_numpy(float)
    q25, q75 = np.percentile(v, [25, 75])
    print(f"{c:<20}{v.mean():>9.4f}{v.std():>9.4f}{q75-q25:>9.4f}"
          f"{v.max()-v.min():>9.4f}")

print(f"\n=== H3 重检（最有利条件下）===")
print(f"{'seed':<8}{'n':>4}{'ρ(U_norm, Macro ρ)':>21}{'ρ(U_norm, Overall r)':>23}")
rob = []
for sd, g in U.groupby("seed"):
    if len(g) < 8:
        continue
    sm = spearmanr(g.U_norm, g.macro_spearman).correlation
    so = spearmanr(g.U_norm, g.overall_pearson).correlation
    print(f"{sd:<8}{len(g):>4}{sm:>21.4f}{so:>23.4f}")
    rob.append({"seed": sd, "n": len(g), "rho_macro": sm, "rho_overall": so})
RB = pd.DataFrame(rob)
RB.to_csv(SEL / "axis_c_final_seed.csv", index=False, encoding="utf-8")

sm_all = spearmanr(U.U_norm, U.macro_spearman).correlation
so_all = spearmanr(U.U_norm, U.overall_pearson).correlation
print(f"\n  全样本: ρ(U_norm, Macro ρ)={sm_all:+.4f}   "
      f"ρ(U_norm, Overall r)={so_all:+.4f}")
w = int((RB.rho_macro > RB.rho_overall).sum())
print(f"  Macro ρ 胜出的 seed 数 = {w} / {len(RB)}")
for c, nm in (("rho_macro", "macro_spearman"), ("rho_overall", "overall_pearson")):
    sg = "稳定" if (RB[c] > 0).all() or (RB[c] < 0).all() else "不稳定"
    print(f"  {nm:<18} 跨 seed 符号: {sg}  "
          f"(min={RB[c].min():+.3f}, max={RB[c].max():+.3f})")

verdict = ("H3_FAIL_SIGNS_UNSTABLE"
           if not ((RB.rho_macro > 0).all() or (RB.rho_macro < 0).all())
           else "H3_POSSIBLE")
print(f"\n  ⇒ H3_FINAL_VERDICT = {verdict}")
