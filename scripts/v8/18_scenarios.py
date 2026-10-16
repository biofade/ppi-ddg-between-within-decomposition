# -*- coding: utf-8 -*-
"""scripts/v8/18_scenarios.py -- V8 PHASE 18：§41 三场景模型推荐 + §42 regret 定稿。

§41 三个**任务场景**（不做通用复合分）：
    A 定量 ΔΔG       -> 首选指标 RMSE / 校准（C_B）
    B 突变文库排序   -> 首选指标 Macro per-target Spearman
    C 湿实验短名单   -> 首选指标 U_norm（AXIS C 最终定义）

对每个场景给出 top-3 推荐 + **该场景下的**并列情况，并回答 §41 的核心问题：
**「按 Overall Pearson 选，会不会推荐同一个模型？」**
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
SC = OUT / "scenarios"
SC.mkdir(parents=True, exist_ok=True)

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
CX = OUT / "selection" / "axis_c_final.csv"
U = pd.read_csv(CX) if CX.exists() else None

u = S[(S.universe == "U_CANONICAL")].copy()
if U is not None:
    u = u.merge(U[["model_name", "U_norm"]].drop_duplicates("model_name"),
                on="model_name", how="left")
    print(f"合并 U_norm：{u.U_norm.notna().sum()} / {len(u)} 个模型有值")
print(f"U_CANONICAL 模型 = {len(u)}")

SCEN = {
    "A_quantitative_ddG": {
        "desc": "定量 ΔΔG 预测（需要绝对数值可信）",
        "primary": "rmse", "asc": True,
        "secondary": ["C_B_row", "overall_pearson", "attenuation_slope"],
    },
    "B_library_ranking": {
        "desc": "突变文库排序（只要序对，不要绝对值）",
        "primary": "macro_spearman", "asc": False,
        "secondary": ["macro_kendall", "pairwise_d05"],
    },
    "C_wetlab_shortlist": {
        "desc": "湿实验短名单（Top-K 命中）",
        "primary": "U_norm", "asc": False,
        "secondary": ["efficiency@10", "ndcg@10"],
    },
}
CONVENTIONAL = "overall_pearson"

rows, rec_rows = [], []
for name, cfg in SCEN.items():
    g = u.dropna(subset=[cfg["primary"]]).copy()
    if not len(g):
        print(f"[skip] {name}: 无 {cfg['primary']} 数据"); continue
    g = g.sort_values(cfg["primary"], ascending=cfg["asc"]).reset_index(drop=True)
    top = g.head(3)
    print(f"\n=== 场景 {name}：{cfg['desc']} ===")
    print(f"  按 `{cfg['primary']}` 选（{'越小' if cfg['asc'] else '越大'}越好）")
    for i, r in top.iterrows():
        print(f"    #{i+1} `{r['model_name']}` ({r['model_family']}) "
              f"{cfg['primary']}={r[cfg['primary']]:.4f} "
              f"| Overall r={r['overall_pearson']:+.4f} "
              f"| Macro ρ={r['macro_spearman']:+.4f} "
              f"| RMSE={r['rmse']:.4f}")

    # §41 核心问题：按 Overall Pearson 选的是不是同一个？
    gc = u.dropna(subset=[CONVENTIONAL, cfg["primary"]]).copy()
    i_conv = int(gc[CONVENTIONAL].idxmax())
    m_conv = gc.loc[i_conv, "model_name"]
    in_top3 = m_conv in set(top.model_name)
    rank_conv = int(g.reset_index(drop=True).index[
        g.model_name == m_conv].tolist()[0]) + 1 if m_conv in set(g.model_name) else -1
    print(f"  §41: 按 Overall Pearson 选 -> `{m_conv}`；"
          f"在该场景排名第 {rank_conv}；{'✅ 在 top-3' if in_top3 else '❌ 不在 top-3'}")

    # 该场景下前 10 的族分布
    fams = g.head(10).model_family.value_counts().to_dict()
    rows.append({
        "scenario": name, "description": cfg["desc"],
        "primary_metric": cfg["primary"],
        "top1": top.model_name.iloc[0], "top1_family": top.model_family.iloc[0],
        "top1_primary_value": float(top[cfg["primary"]].iloc[0]),
        "top1_overall_r": float(top.overall_pearson.iloc[0]),
        "top1_macro": float(top.macro_spearman.iloc[0]),
        "top1_rmse": float(top.rmse.iloc[0]),
        "overall_pearson_pick": m_conv,
        "overall_pearson_pick_rank_in_scenario": rank_conv,
        "overall_pearson_pick_in_top3": in_top3,
        "top10_family_counts": str(fams),
    })
    rec_rows.append({"scenario": name, "criterion": cfg["primary"],
                     "rank": 1, "model": top.model_name.iloc[0],
                     "value": float(top[cfg["primary"]].iloc[0])})
    for i, r in top.iloc[1:].iterrows():
        rec_rows.append({"scenario": name, "criterion": cfg["primary"],
                         "rank": i + 1, "model": r["model_name"],
                         "value": float(r[cfg["primary"]])})
    # 传统对照
    rec_rows.append({"scenario": name, "criterion": CONVENTIONAL, "rank": 1,
                     "model": m_conv, "value": float(gc.loc[i_conv, CONVENTIONAL])})

D = pd.DataFrame(rows)
D.to_csv(SC / "scenario_recommendations.csv", index=False, encoding="utf-8")
pd.DataFrame(rec_rows).to_csv(SC / "scenario_top3.csv", index=False,
                             encoding="utf-8")

print("\n" + "=" * 70)
print("§41 核心结果：三个场景 vs 惯例（Overall Pearson）")
print("=" * 70)
print(f"{'场景':<24}{'按场景首选':<26}{'按 Overall r':<26}{'同?'}")
for r in D.itertuples():
    same = "✅" if r.top1 == r.overall_pearson_pick else "❌"
    print(f"{r.scenario:<24}{r.top1[:24]:<26}{r.overall_pearson_pick[:24]:<26}{same}")
n_diff = int((D.top1 != D.overall_pearson_pick).sum())
print(f"\n  ⇒ 三个场景中，按 Overall Pearson 会给出**不同首选模型**的个数 = "
      f"{n_diff} / {len(D)}")
print(f"  ⇒ 按 Overall Pearson 选的模型进入该场景 top-3 的个数 = "
      f"{int(D.overall_pearson_pick_in_top3.sum())} / {len(D)}")

lines = ["# V8 Scenario-Specific Recommendations（§41/§42）\n",
         "> **不使用通用复合分**（§41）。每个场景只按其**任务相关**指标排序。\n",
         "| 场景 | 说明 | 主排序指标 | 首选模型 | 该指标值 | 按 Overall r 会选 | 场景内排名 | 在 top-3? |",
         "|---|---|---|---|---:|---|---:|:--:|"]
for r in D.itertuples():
    lines.append(f"| `{r.scenario}` | {r.description} | `{r.primary_metric}` | "
                 f"`{r.top1}` | {r.top1_primary_value:.4f} | "
                 f"`{r.overall_pearson_pick}` | "
                 f"{r.overall_pearson_pick_rank_in_scenario} | "
                 f"{'✅' if r.overall_pearson_pick_in_top3 else '❌'} |")
lines += ["\n## 核心问题（§41）\n",
          f"- 三场景中，按 Overall Pearson 会给出**不同首选模型**的个数 = "
          f"**{n_diff} / {len(D)}**",
          f"- 按 Overall Pearson 选的模型进入该场景 top-3 的个数 = "
          f"**{int(D.overall_pearson_pick_in_top3.sum())} / {len(D)}**",
          "\n> ⚠️ 场景 C 的效用轴已证明在数据层面缺乏功效（见 `AXIS_C_DIAGNOSIS.md`），",
          "> 因此场景 C 的推荐**不可作为强结论**，只作记录。\n"]
(OUT / "SCENARIO_RECOMMENDATIONS.md").write_text("\n".join(lines) + "\n",
                                                 encoding="utf-8")
print(f"\n写出 scenario_recommendations.csv / scenario_top3.csv / "
      f"SCENARIO_RECOMMENDATIONS.md")
