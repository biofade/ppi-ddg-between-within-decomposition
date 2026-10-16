# -*- coding: utf-8 -*-
"""scripts/v8/12c_selection_diagnostic.py -- 选择效用指标的**区分度诊断**。

R4 出现一个看似矛盾的现象：
  * `Overall Pearson` 选中 `L1_random_A1_base_set`，Eff@10 regret 仅 **0.0141**（近乎最优）
  * 但模型层面 `Overall Pearson` 与 `Eff@10` 的秩相关只有 **−0.007**
  * 而 `Macro Spearman` 选的模型 regret 高达 **0.696**

若 Eff@10 的**值域极窄**，则 normalized regret 会被小数放大，
「regret 大」并不代表绝对差异大。本脚本量化这一点。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)

COLS = ["efficiency@5", "efficiency@10", "ndcg@5", "ndcg@10",
        "hit@1", "hit@5", "hit@10", "precision@10", "ef@1", "ef@5",
        "macro_spearman", "overall_pearson", "C_W"]
print("=== 选择/排序指标的区分度（按 universe）===")
rows = []
for uni, g in S.groupby("universe"):
    print(f"\n[{uni}]  n_models = {len(g)}")
    print(f"{'metric':<18}{'min':>9}{'p10':>9}{'median':>9}{'p90':>9}{'max':>9}"
          f"{'range':>9}{'iqr':>9}{'sd':>9}")
    for c in COLS:
        if c not in g.columns:
            continue
        v = g[c].dropna().to_numpy(float)
        if len(v) < 5:
            continue
        q10, q50, q90 = np.percentile(v, [10, 50, 90])
        q25, q75 = np.percentile(v, [25, 75])
        print(f"{c:<18}{v.min():>9.4f}{q10:>9.4f}{q50:>9.4f}{q90:>9.4f}"
              f"{v.max():>9.4f}{v.max()-v.min():>9.4f}{q75-q25:>9.4f}{v.std():>9.4f}")
        rows.append({"universe": uni, "metric": c, "n": len(v),
                     "min": v.min(), "p10": q10, "median": q50, "p90": q90,
                     "max": v.max(), "range": v.max() - v.min(),
                     "iqr": q75 - q25, "sd": v.std()})
D = pd.DataFrame(rows)
D.to_csv(OUT / "selection" / "metric_discrimination.csv", index=False,
         encoding="utf-8")

print("\n=== 绝对差异 vs normalized regret（U_CANONICAL）===")
g = S[S.universe == "U_CANONICAL"].dropna(
    subset=["efficiency@10", "ndcg@10", "overall_pearson"]).reset_index(drop=True)
for ut in ("efficiency@10", "ndcg@10"):
    best, worst = g[ut].max(), g[ut].min()
    rng = best - worst
    print(f"\n  {ut}: best={best:.4f} worst={worst:.4f} range={rng:.4f} "
          f"sd={g[ut].std():.4f}")
    for rule in ("overall_pearson", "macro_spearman"):
        i = int(g[rule].idxmax())
        got = g.loc[i, ut]
        print(f"    {rule:16s} -> abs_gap={best-got:.4f}  "
              f"norm_regret={(best-got)/rng:.4f}  "
              f"in_sd_units={(best-got)/max(1e-12, g[ut].std()):.3f}")

print("\n=== 判读 ===")
for r in D[D.metric.isin(["efficiency@10", "ndcg@10"])].itertuples():
    print(f"  {r.universe:12s} {r.metric:14s} range={r.range:.4f} sd={r.sd:.4f} "
          f"| range/sd={r.range/max(1e-12,r.sd):.2f}")
print()
print("  若 range 很小（<0.1）则 Top-K 类效用指标在本 benchmark 上**区分度有限**，")
print("  normalized regret 会被过度放大；此类指标不应作为唯一决策依据。")
