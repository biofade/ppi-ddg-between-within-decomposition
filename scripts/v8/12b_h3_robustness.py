# -*- coding: utf-8 -*-
"""scripts/v8/12b_h3_robustness.py -- H3 的稳健性检验（诚实版）。

初步结果出现**自相矛盾**：
  * 聚合上「排序轴 vs 传统全局」有利于排序轴（U_CANONICAL Δρ=+0.14）
  * 但 **`C_W`（属排序轴）对 Eff@10 是负相关（−0.214）**，
    且 **`Hit@10` 几乎与所有指标负相关**
  * 更关键：**按 Overall Pearson 选的模型在 Eff@10 上 regret 仅 0.014（近乎最优）**

⇒ 不能简单判 H3 PASS。本脚本做三项稳健性检验：
  R1 分 seed 重复（排除 seed 混杂）
  R2 只保留「所有指标都可用」的高覆盖子集
  R3 逐效用指标分别判断（不合并），并检查 Hit@10 是否退化
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
DA = OUT / "decision_alignment"

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
PRED_CAL = ["overall_pearson", "rmse", "target_mean_pearson"]
PRED_RANK = ["macro_spearman", "macro_kendall", "pairwise_d05"]
UTILS = ["efficiency@10", "ndcg@10", "hit@10"]

rows = []
for uni, g in S.groupby("universe"):
    g = g.dropna(subset=PRED_CAL + PRED_RANK + UTILS + ["C_W"]).reset_index(drop=True)
    if len(g) < 12:
        continue
    # R1 逐 seed
    g["seed"] = g.model_name.str.extract(r"_s(\d+)$")[0]
    for seed, gs in g.groupby("seed"):
        if len(gs) < 8:
            continue
        for label, cols in (("calibration", PRED_CAL), ("ranking", PRED_RANK)):
            for ut in UTILS:
                v = [spearmanr(gs[c], gs[ut]).correlation for c in cols]
                v = [x for x in v if np.isfinite(x)]
                rows.append({"universe": uni, "subset": f"seed={seed}",
                             "n": len(gs), "axis": label, "utility": ut,
                             "median_spearman": float(np.median(v)) if v else np.nan})
    # R2 全部模型（合并 seed）
    for label, cols in (("calibration", PRED_CAL), ("ranking", PRED_RANK),
                        ("C_W_alone", ["C_W"])):
        for ut in UTILS:
            v = [spearmanr(g[c], g[ut]).correlation for c in cols]
            v = [x for x in v if np.isfinite(x)]
            rows.append({"universe": uni, "subset": "all", "n": len(g),
                         "axis": label, "utility": ut,
                         "median_spearman": float(np.median(v)) if v else np.nan})

D = pd.DataFrame(rows)
D.to_csv(DA / "h3_robustness.csv", index=False, encoding="utf-8")
print("=== H3 稳健性（median Spearman vs utility）===")
for uni, g in D.groupby("universe"):
    print(f"\n[{uni}]")
    p = g.pivot_table(index=["subset", "axis"], columns="utility",
                      values="median_spearman")
    print(p.round(3).to_string())

# R3：Hit@10 是否退化（几乎所有模型都命中 top-1）
print("\n=== R3：Hit@10 分布（是否退化到接近常数）===")
for uni, g in S.groupby("universe"):
    v = g["hit@10"].dropna()
    if len(v):
        print(f"  {uni}: n={len(v)} mean={v.mean():.4f} sd={v.std():.4f} "
              f"min={v.min():.4f} max={v.max():.4f} "
              f"| 取值数={v.nunique()}")

# R4：直接检查 §42 regret——按 Overall Pearson 选 vs 按排序指标选
print("\n=== R4：selection regret 对比（oracle 为描述性上界，§43）===")
for uni, g in S.groupby("universe"):
    g = g.dropna(subset=["overall_pearson", "macro_spearman", "C_W",
                         "efficiency@10", "ndcg@10"]).reset_index(drop=True)
    if len(g) < 12:
        continue
    print(f"\n[{uni}] n={len(g)}")
    for rule, col, asc in (("Overall Pearson", "overall_pearson", False),
                           ("Macro Spearman", "macro_spearman", False),
                           ("C_W", "C_W", False)):
        i = int(g[col].idxmin() if asc else g[col].idxmax())
        line = f"  {rule:16s} -> {g.loc[i,'model_name'][:34]:34s}"
        for ut in ("efficiency@10", "ndcg@10"):
            best = g[ut].max(); got = g.loc[i, ut]
            rng = best - g[ut].min()
            line += f"  {ut}: got={got:.4f} regret={((best-got)/rng):.4f}"
        print(line)
    for ut in ("efficiency@10", "ndcg@10"):
        j = int(g[ut].idxmax())
        print(f"  oracle-best {ut:15s} -> {g.loc[j,'model_name'][:44]}")

print("\n写出 h3_robustness.csv")
