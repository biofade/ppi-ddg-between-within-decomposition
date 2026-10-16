# -*- coding: utf-8 -*-
"""scripts/v8/11b_axis_c_floor.py -- AXIS C「发现地板」的**数学来源**诊断（§34 方法论）。

## 要回答的问题

`Efficiency@10` 在 U_CANONICAL 上 IQR 仅 0.0143。这是**指标定义的缺陷**，
还是**数据本身**给 Top-K 效用设了上限？

## 做法：解析计算每个 target 的「随机地板」与「理论上限」

对 target 内 n 条突变，真值 gain（有利度）`g_i = 1 + (n − rank_i)`，`rank` 按
`ΔΔG` 升序（越负越有利）。令 `G_total = Σ g_i`。

  * **随机地板**  `E_rand`：从 n 条里**无放回**随机取 k 条，其 gain 之和的期望
    占 `G_total` 的比例。因抽样对称，`E_rand = k/n`（精确，无需模拟）。
  * **理论上限**  `E_max`：取 gain 最大的 k 条，`Σ top-k g_i / G_total`。
  * **可区分区间**  `E_max − E_rand`：模型**最多**能比随机好多少。
  * **归一化效用**  `U_norm = (Eff@10 − E_rand) / (E_max − E_rand)`：
    0 = 等同随机，1 = 理论最优。

**若 `E_max − E_rand` 普遍很小**，则「地板」来自 **gain 的分布形状**（即数据），
而不是指标写错 —— 且 Top-K 效用**原则上无法区分模型**。反之则是指标问题。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
SEL = OUT / "selection"
SEL.mkdir(parents=True, exist_ok=True)

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)

K = 10
MIN_N = 5

# ---------- 1. 每个 target 的地板 / 上限 ----------
rows = []
for t in np.unique(TGT):
    i = np.where(TGT == t)[0]
    n = i.size
    if n < MIN_N:
        continue
    yg = Y[i]
    order = np.argsort(yg, kind="stable")        # 越小越有利
    rank = np.empty(n, int); rank[order] = np.arange(n)
    gain = 1.0 + (n - rank)                      # 1..n
    k = max(1, min(K, n))
    Gtot = gain.sum()
    Emax = np.sort(gain)[::-1][:k].sum() / Gtot
    Erand = k / n
    rows.append({"target_id": t, "n": n, "k": k,
                 "gain_max": gain.max(), "gain_sum": Gtot,
                 "E_rand": Erand, "E_max": Emax,
                 "range": Emax - Erand,
                 "Gini_gain": float(
                     1 - 2 * np.sum(np.sort(gain) * np.arange(1, n + 1)) /
                     (n * Gtot) + (n + 1) / n) if n > 1 else 0.0})
G = pd.DataFrame(rows)
G.to_csv(SEL / "axis_c_floor_by_target.csv", index=False, encoding="utf-8")

print(f"eligible targets (n>={MIN_N}) = {len(G)} | K = {K}")
print("\n=== 每个 target 的随机地板 / 理论上限 / 可区分区间 ===")
print(f"{'quantile':<12}{'n':>7}{'E_rand':>10}{'E_max':>10}{'range':>10}")
for q in (0, 10, 25, 50, 75, 90, 100):
    v = np.percentile(G.n, q)
    print(f"  n p{q:<8}{v:>7.1f}"
          f"{np.percentile(G.E_rand,q):>10.4f}"
          f"{np.percentile(G.E_max,q):>10.4f}"
          f"{np.percentile(G['range'],q):>10.4f}")
print()
print(f"  E_rand  : mean={G.E_rand.mean():.4f} median={G.E_rand.median():.4f}")
print(f"  E_max   : mean={G.E_max.mean():.4f} median={G.E_max.median():.4f}")
print(f"  range   : mean={G['range'].mean():.4f} median={G['range'].median():.4f}")
print(f"  range/1 : 即「模型最多能比随机好多少」")
print()
rn = G['range'] / (1 - G.E_rand)
print(f"  归一化可区分区间 (E_max-E_rand)/(1-E_rand): "
      f"mean={rn.mean():.4f} median={rn.median():.4f}")

# ---------- 2. 实测 Eff@10 相对地板的归一化 ----------
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))]
cache, U = {}, []
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
    # 逐 target 计算 Eff@10 与归一化效用
    efs, unorms = [], []
    for t in G.target_id:
        i = np.where(TGT == t)[0]
        yg, pg = Y[i], a[i]
        n = i.size
        order = np.argsort(yg, kind="stable")
        rank = np.empty(n, int); rank[order] = np.arange(n)
        gain = 1.0 + (n - rank)
        k = max(1, min(K, n))
        ideal = np.argsort(-gain, kind="stable")[:k]
        pred = np.argsort(pg, kind="stable")[:k]
        gid = gain[ideal].sum()
        if gid <= 0:
            continue
        ef = gain[pred].sum() / gid
        efs.append(ef)
        row = G[G.target_id == t].iloc[0]
        rng = row["range"]
        unorms.append((ef - row.E_rand) / rng if rng > 1e-12 else np.nan)
    U.append({"model_name": mname, "model_family": r["model_family"],
              "macro_spearman": r["macro_spearman"],
              "overall_pearson": r["overall_pearson"],
              "eff10_v1": float(np.mean(efs)),
              "util_norm": float(np.nanmean(unorms))})
U = pd.DataFrame(U)
U.to_csv(SEL / "axis_c_normalized_utility.csv", index=False, encoding="utf-8")

print(f"\n=== 归一化效用 U_norm（0=等同随机, 1=理论最优）n={len(U)} 模型 ===")
for c in ("eff10_v1", "util_norm", "macro_spearman"):
    v = U[c].dropna()
    q25, q75 = np.percentile(v, [25, 75])
    print(f"  {c:<16} mean={v.mean():.4f} sd={v.std():.4f} "
          f"IQR={q75-q25:.4f} range={v.max()-v.min():.4f}")
print()
print("=== U_norm 与其它指标的秩相关 ===")
for c in ("macro_spearman", "overall_pearson", "C_W" if "C_W" in U else "eff10_v1"):
    if c in U.columns:
        print(f"  ρ(U_norm, {c}) = "
              f"{spearmanr(U.util_norm, U[c]).correlation:+.4f}")

# ---------- 3. 判读 ----------
print("\n" + "=" * 68)
print("判读")
print("=" * 68)
med_range = G['range'].median()
med_rand = G.E_rand.median()
if med_range < 0.10:
    verdict = "FLOOR_IS_DATA_IMPOSED"
    msg = ("可区分区间中位数 < 0.10 ⇒ **Top-K 效用被数据本身托住**。\n"
           "  任何模型最多只能比随机好 ~%.1f 个百分点。\n"
           "  ⇒ AXIS C 的「发现地板」**来自 gain 分布形状（数据），不是指标定义错误**。"
           % (100 * med_range))
elif med_range < 0.25:
    verdict = "FLOOR_PARTIALLY_DATA_IMPOSED"
    msg = (f"可区分区间中位数 = {med_range:.4f} ⇒ 区间偏窄但非极端；\n"
           "  指标可改进，但数据也贡献了显著压缩。")
else:
    verdict = "FLOOR_IS_METRIC_DEFECT"
    msg = (f"可区分区间中位数 = {med_range:.4f} ⇒ 区间充足，"
           "「地板」更可能来自指标定义而非数据。")
print(f"  median(E_rand) = {med_rand:.4f}   median(range) = {med_range:.4f}")
print(f"  ⇒ AXIS_C_FLOOR_DIAGNOSIS = {verdict}")
print(f"  {msg}")
