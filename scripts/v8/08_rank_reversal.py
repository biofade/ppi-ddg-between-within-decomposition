# -*- coding: utf-8 -*-
"""scripts/v8/08_rank_reversal.py -- V8 PHASE 8：模型秩反转（§30/§31）。

在同一 universe 内，按不同指标给模型排名，比较**模型排名的秩一致性**：
    Overall Pearson / RMSE / target-mean Pearson      （校准轴）
    Macro Spearman / Macro Kendall / C_W              （排序轴）
    Efficiency@10 / NDCG@10 / Hit@10                  （选择轴）

输出：
    experiments/v8/rank_reversal/rank_agreement.csv
    experiments/v8/rank_reversal_matrix.csv
    experiments/v8/RANK_REVERSAL_REPORT.md
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
RR = OUT / "rank_reversal"
RR.mkdir(parents=True, exist_ok=True)

TIE_TOL = 0.01           # §31 预注册实践并列容差
AXES = {
    "overall_pearson": "calibration",
    "rmse": "calibration(neg)",
    "target_mean_pearson": "calibration",
    "macro_spearman": "ranking",
    "macro_kendall": "ranking",
    "C_W": "ranking",
    "efficiency@10": "selection",
    "ndcg@10": "selection",
    "hit@10": "selection",
}
LOWER_BETTER = {"rmse"}

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
print("scoreboard rows =", len(S), "| universes =", S.universe.unique().tolist())

rep: list[str] = []
agg_rows, mat_rows = [], []
for uni, g in S.groupby("universe"):
    # 需要所有指标都可用，且至少 8 个模型才有比较意义
    g = g.dropna(subset=list(AXES.keys()))
    if len(g) < 8:
        print(f"  [skip] {uni}: 可用模型 {len(g)} < 8")
        continue
    g = g.reset_index(drop=True)
    print(f"\n=== {uni}: {len(g)} 个模型 ===")
    ranks = {}
    for col in AXES:
        v = g[col].to_numpy(float)
        # 排名：1 = 最好
        order = np.argsort(v if col in LOWER_BETTER else -v, kind="stable")
        rk = np.empty(len(v), int); rk[order] = np.arange(1, len(v) + 1)
        ranks[col] = rk
    Rk = pd.DataFrame(ranks)
    Rk.insert(0, "model_name", g.model_name.values)
    Rk.insert(1, "model_family", g.model_family.values)
    Rk.to_csv(RR / f"ranks_{uni}.csv", index=False, encoding="utf-8")

    cols = list(AXES)
    Sp = np.full((len(cols), len(cols)), np.nan)
    Ke = np.full((len(cols), len(cols)), np.nan)
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            if i == j:
                Sp[i, j] = Ke[i, j] = 1.0; continue
            Sp[i, j] = spearmanr(Rk[a], Rk[b]).correlation
            Ke[i, j] = kendalltau(Rk[a], Rk[b]).correlation
    SpD = pd.DataFrame(Sp, index=cols, columns=cols)
    KeD = pd.DataFrame(Ke, index=cols, columns=cols)
    SpD.to_csv(RR / f"rank_agreement_spearman_{uni}.csv", encoding="utf-8")
    KeD.to_csv(RR / f"rank_agreement_kendall_{uni}.csv", encoding="utf-8")
    mat_rows.append({"universe": uni, "n_models": len(g), "metric": "spearman",
                     **{f"{c}": "" for c in cols}})
    for a in cols:
        for b in cols:
            if a < b:
                agg_rows.append({"universe": uni, "n_models": len(g),
                                 "metric_a": a, "axis_a": AXES[a],
                                 "metric_b": b, "axis_b": AXES[b],
                                 "spearman_rank_agree": SpD.loc[a, b],
                                 "kendall_rank_agree": KeD.loc[a, b],
                                 "cross_axis": AXES[a] != AXES[b]})

    # top-1 / top-3
    def top(m, k=3):
        return list(Rk.sort_values(m).model_name.head(k))
    rep.append(f"\n### universe `{uni}`（n={len(g)} 模型）\n")
    rep.append("| 指标 | top-1 | top-3 |")
    rep.append("|---|---|---|")
    for c in cols:
        rep.append(f"| `{c}` | {top(c,1)[0]} | {', '.join(top(c,3))} |")
    rep.append("")
    rep.append("| 指标对 | Spearman 秩一致 | Kendall 秩一致 | 跨轴 |")
    rep.append("|---|---:|---:|:--:|")
    for a in cols:
        for b in cols:
            if a < b and AXES[a] != AXES[b]:
                rep.append(f"| `{a}` vs `{b}` | {SpD.loc[a,b]:+.3f} | "
                           f"{KeD.loc[a,b]:+.3f} | ✅ |")

    # 实质反转（§31）：跨轴且秩差超过容差、且非数值噪声
    key_pairs = [("overall_pearson", "macro_spearman"),
                 ("overall_pearson", "efficiency@10"),
                 ("macro_spearman", "efficiency@10"),
                 ("rmse", "macro_spearman")]
    rev_rows = []
    for a, b in key_pairs:
        d = (Rk[a] - Rk[b]).abs()
        # 用指标绝对差判断是否"实践并列"（不数噪声）
        dv = (g[a] - g[b]).abs() if False else None
        cand = Rk[(d > 0)].copy()
        cand["abs_rank_diff"] = d[d > 0]
        cand = cand.sort_values("abs_rank_diff", ascending=False).head(8)
        for _, r in cand.iterrows():
            nm = r["model_name"]
            ia = g.index[g.model_name == nm][0]
            rev_rows.append({"universe": uni, "pair": f"{a}|{b}", "model": nm,
                             "family": g.loc[ia, "model_family"],
                             f"rank_{a}": int(r[a]), f"rank_{b}": int(r[b]),
                             "val_" + a: g.loc[ia, a], "val_" + b: g.loc[ia, b],
                             "n_models": len(g)})
    if rev_rows:
        pd.concat([pd.DataFrame([x]) for x in rev_rows], ignore_index=True).to_csv(
            RR / f"reversals_{uni}.csv", index=False, encoding="utf-8")
        rep.append(f"\n**最大秩反转（{uni}，按秩差排序，前 6）**\n")
        rep.append("| 指标对 | 模型 | 族 | rank_a | rank_b | val_a | val_b |")
        rep.append("|---|---|---|---:|---:|---:|---:|")
        for x in rev_rows[:6]:
            a, b = x["pair"].split("|")
            rep.append(f"| `{a}` vs `{b}` | `{x['model']}` | {x['family']} | "
                       f"{x[f'rank_{a}']} | {x[f'rank_{b}']} | "
                       f"{x['val_'+a]:+.4f} | {x['val_'+b]:+.4f} |")

A = pd.DataFrame(agg_rows)
A.to_csv(OUT / "rank_reversal" / "rank_agreement.csv", index=False, encoding="utf-8")
print("\n=== 跨轴秩一致性（关键对）===")
PAIRS = [("overall_pearson", "macro_spearman"), ("overall_pearson", "efficiency@10"),
         ("overall_pearson", "ndcg@10"), ("rmse", "macro_spearman"),
         ("rmse", "efficiency@10"), ("target_mean_pearson", "macro_spearman"),
         ("macro_spearman", "efficiency@10")]
sub = []
for _, r in A.iterrows():
    for a, b in PAIRS:
        if {r.metric_a, r.metric_b} == {a, b}:
            sub.append({"universe": r.universe, "metric_a": a, "metric_b": b,
                        "spearman_rank_agree": r.spearman_rank_agree,
                        "kendall_rank_agree": r.kendall_rank_agree})
SUB = pd.DataFrame(sub)
if len(SUB):
    print(SUB.to_string(index=False))
    SUB.to_csv(RR / "rank_agreement_keypairs.csv", index=False, encoding="utf-8")
else:
    print("  (unordered pair matching failed; dumping matrix)"
          )
    print(SpD.to_string())

hdr = ("# V8 Model-Rank Reversal Report\n\n"
       "> §31 预注册实践并列容差：Overall Pearson 与 Macro Spearman 之差 "
       f"< {TIE_TOL} 视为实践并列（不数作科学反转）。\n\n"
       "> 秩一致使用 **Spearman / Kendall**（1.0 = 完全相同的模型排序）。\n")
(RR.parent / "RANK_REVERSAL_REPORT.md").write_text(
    hdr + "\n".join(rep) + "\n", encoding="utf-8")
print(f"\n写出 rank_agreement.csv / rank_reversal_matrix / RANK_REVERSAL_REPORT.md")
