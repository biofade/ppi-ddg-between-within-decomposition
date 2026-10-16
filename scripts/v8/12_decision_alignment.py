# -*- coding: utf-8 -*-
"""scripts/v8/12_decision_alignment.py -- V8 PHASE 12 / 核心实验 F：决策对齐（§38/§39/§42/§43）。

**问题（H3）**：哪些模型层面的指标，能预测该模型在**湿实验候选筛选**上的实际效用？

做法：在**同一 universe** 内，把每个模型看成一个观测，
用模型**排名**之间的 Spearman / Kendall 关联来衡量：

  候选指标（predictor）: Overall Pearson / RMSE / target-mean Pearson /
                        Macro Spearman / Macro Kendall / Pairwise δ=0.5 / C_W
  效用（target）       : Efficiency@10 / NDCG@10 / Hit@10

术语纪律（§39）：这叫 **decision alignment**，**不是** external validation
（同一批实验标签参与了双方）。

另做 §42 **model-selection regret**：若按 Overall Pearson 选模型，
它在 Eff@10 / NDCG@10 上落后"事后最优"多少（oracle 仅作**描述性上界**，§43）。
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
DA.mkdir(parents=True, exist_ok=True)

PREDICTORS = {
    "overall_pearson": "calibration",
    "rmse": "calibration(neg)",
    "target_mean_pearson": "calibration",
    "macro_spearman": "ranking",
    "macro_kendall": "ranking",
    "pairwise_d05": "ranking",
    "C_W": "ranking",
}
UTILITIES = {"efficiency@10": "selection", "ndcg@10": "selection",
             "hit@10": "selection"}

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
rows, regret_rows = [], []
for uni, g in S.groupby("universe"):
    need = list(PREDICTORS) + list(UTILITIES)
    g = g.dropna(subset=need)
    if len(g) < 8:
        print(f"[skip] {uni}: n={len(g)}")
        continue
    g = g.reset_index(drop=True)
    print(f"\n=== {uni}: n={len(g)} 模型 ===")
    print(f"{'predictor':<22}{'axis':<18}"
          + "".join(f"{u:>16}" for u in UTILITIES))
    for pr, ax in PREDICTORS.items():
        line = f"{pr:<22}{ax:<18}"
        for ut in UTILITIES:
            sp = spearmanr(g[pr], g[ut]).correlation
            line += f"{sp:>+16.3f}"
            rows.append({"universe": uni, "n_models": len(g),
                         "predictor": pr, "predictor_axis": ax,
                         "utility": ut,
                         "spearman": float(sp),
                         "kendall": float(kendalltau(g[pr], g[ut]).correlation)})
        print(line)

    # §42 selection regret
    for sel_metric in ("overall_pearson", "macro_spearman", "C_W",
                       "target_mean_pearson"):
        i_sel = int(g[sel_metric].idxmax()) if sel_metric != "rmse" \
            else int(g["rmse"].idxmin())
        m_sel = g.loc[i_sel, "model_name"]
        rec = {"universe": uni, "n_models": len(g),
               "selection_rule": sel_metric, "selected_model": m_sel}
        for ut in UTILITIES:
            got = float(g.loc[i_sel, ut])
            best = float(g[ut].max())
            worst = float(g[ut].min())
            rng = best - worst
            rec[f"got_{ut}"] = got
            rec[f"oracle_best_{ut}"] = best
            rec[f"regret_abs_{ut}"] = best - got
            rec[f"regret_norm_{ut}"] = (best - got) / rng if rng > 0 else np.nan
            rec[f"is_oracle_{ut}"] = bool(abs(got - best) < 1e-12)
        regret_rows.append(rec)

D = pd.DataFrame(rows)
D.to_csv(DA / "decision_alignment.csv", index=False, encoding="utf-8")
G = pd.DataFrame(regret_rows)
G.to_csv(DA / "selection_regret.csv", index=False, encoding="utf-8")

print("\n=== 按轴聚合（Spearman 中位数）===")
piv = D.pivot_table(index=["universe", "predictor_axis"], columns="utility",
                    values="spearman", aggfunc="median")
print(piv.round(3).to_string())

print("\n=== H3 直接对照（每个 universe）===")
h3 = []
for uni, g in D.groupby("universe"):
    def med(ax):
        v = g[g.predictor_axis.str.startswith(ax)].spearman
        return float(v.median()) if len(v) else np.nan
    cal = med("calibration")
    rank = med("ranking")
    # 仅用整体 Pearson / RMSE 作为"传统全局"代表
    trad = float(g[g.predictor.isin(["overall_pearson", "rmse"])].spearman.median())
    h3.append({"universe": uni, "n_models": int(g.n_models.iloc[0]),
               "median_spearman_calibration_axis": cal,
               "median_spearman_ranking_axis": rank,
               "median_spearman_traditional_global": trad,
               "ranking_better_than_traditional": bool(rank > trad),
               "delta_rank_minus_trad": rank - trad})
H = pd.DataFrame(h3)
H.to_csv(DA / "h3_summary.csv", index=False, encoding="utf-8")
print(H.round(4).to_string(index=False))

print("\n=== §42 Selection regret（按 Overall Pearson 选模型）===")
if len(G):
    show = ["universe", "selection_rule", "selected_model",
            "regret_norm_efficiency@10", "regret_norm_ndcg@10",
            "is_oracle_efficiency@10", "is_oracle_ndcg@10"]
    print(G[[c for c in show if c in G.columns]].round(4).to_string(index=False))

lines = ["# V8 Decision-Alignment Report（§38/§39/§42/§43 / H3）\n",
         "> **术语（§39）**：这是 **decision alignment**，**不是** external validation ——",
         "> 同一批实验标签同时参与了预测指标与效用指标的计算。\n",
         "## 1. 模型层面：各指标与选择效用的秩关联（Spearman）\n"]
for uni, g in D.groupby("universe"):
    lines.append(f"\n### universe `{uni}`（n={int(g.n_models.iloc[0])} 模型）\n")
    lines.append("| predictor | axis | " + " | ".join(UTILITIES) + " |")
    lines.append("|---|---|" + "---:|" * len(UTILITIES))
    p2 = g.pivot_table(index=["predictor", "predictor_axis"], columns="utility",
                       values="spearman")
    for (pr, ax), r in p2.iterrows():
        lines.append(f"| `{pr}` | {ax} | "
                     + " | ".join(f"{r[u]:+.3f}" for u in UTILITIES) + " |")
lines.append("\n## 2. H3 对照：排序轴 vs 传统全局指标\n")
lines.append("| universe | n | 校准轴中位 ρ | **排序轴中位 ρ** | 传统全局中位 ρ | 排序轴更优? |")
lines.append("|---|---:|---:|---:|---:|:--:|")
for r in H.itertuples():
    lines.append(f"| `{r.universe}` | {r.n_models} | "
                 f"{r.median_spearman_calibration_axis:+.3f} | "
                 f"**{r.median_spearman_ranking_axis:+.3f}** | "
                 f"{r.median_spearman_traditional_global:+.3f} | "
                 f"{'✅' if r.ranking_better_than_traditional else '❌'} |")
lines.append("\n## 3. Model-selection regret（§42）\n")
lines.append("> oracle-best 为**事后描述性上界**，**不可部署**（§43）。\n")
if len(G):
    lines.append("```")
    lines.append(G[[c for c in ["universe", "selection_rule", "selected_model",
                                "regret_norm_efficiency@10", "regret_norm_ndcg@10",
                                "is_oracle_efficiency@10", "is_oracle_ndcg@10"]
                    if c in G.columns]].round(4).to_string(index=False))
    lines.append("```")
(OUT / "DECISION_ALIGNMENT_REPORT.md").write_text("\n".join(lines) + "\n",
                                                  encoding="utf-8")
print("\n写出 decision_alignment.csv / h3_summary.csv / selection_regret.csv / "
      "DECISION_ALIGNMENT_REPORT.md")
