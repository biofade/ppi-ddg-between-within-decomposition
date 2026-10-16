# -*- coding: utf-8 -*-
"""scripts/v8/19_figures_tables.py -- V8 PHASE 19 / §79 七张图 + §80 八张表。

§79 七张图（全部由已有 CSV 生成，不重算模型）：
  fig1 rank-agreement 热图（跨轴秩一致性）
  fig2 C_B / C_W 分解散点（校准信号 vs 排序信号）
  fig3 split inflation：三轴随划分严格度的变化
  fig4 remedy 链：每个 remedy 的 ΔMacro ρ vs ΔOverall r
  fig5 AXIS C：随机地板 / 理论上限 随 target 大小 n 的变化
  fig6 weighting sensitivity：C_B、C_W 在 5 种加权方案下的极差
  fig7 capability landscape：三轴二维投影

§80 八张表 -> experiments/v8/tables/T1..T8
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# 中文字体
for cand in ("Microsoft YaHei", "SimHei", "DengXian", "SimSun"):
    try:
        font_manager.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [cand]
        print(f"中文字体 = {cand}")
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.bbox"] = "tight"

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
u = S[S.universe == "U_CANONICAL"].copy()
CXV = OUT / "selection" / "axis_c_final.csv"
if CXV.exists():
    U = pd.read_csv(CXV)
    u = u.merge(U[["model_name", "U_norm"]].drop_deduplicated if False
                else U[["model_name", "U_norm"]].drop_duplicates("model_name"),
                on="model_name", how="left")
print(f"U_CANONICAL 模型 = {len(u)}")

# ============ fig1 秩一致性热图 ============
RA = pd.read_csv(OUT / "rank_reversal" / "rank_agreement.csv")
rc = pd.read_csv(OUT / "rank_reversal" / "rank_agreement_spearman_U_CANONICAL.csv",
                 index_col=0)
fig, ax = plt.subplots(figsize=(9.2, 7.4))
M = rc.to_numpy(float)
im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(rc.columns)))
ax.set_xticklabels(rc.columns, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(rc.index)))
ax.set_yticklabels(rc.index, fontsize=9)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isfinite(M[i, j]):
            ax.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center",
                    fontsize=7.5,
                    color="white" if abs(M[i, j]) > 0.55 else "black")
ax.set_title("图1  模型排名的跨轴秩一致性（Spearman, U_CANONICAL, 160 模型）",
             fontsize=11.5, pad=12)
fig.colorbar(im, ax=ax, shrink=0.82, label="模型排名的 Spearman 秩一致")
fig.text(0.5, -0.055,
         "对角线为 1；非对角越接近 0 表示两指标对模型的排序越不一致。\n"
         "Overall Pearson 与 Macro per-target Spearman = +0.18；与 Efficiency@10 = −0.01。",
         ha="center", fontsize=8.5, color="#333")
fig.savefig(FIG / "fig1_rank_agreement_heatmap.png")
plt.close(fig)
print("  fig1 ✅")

# ============ fig2 C_B / C_W 分解 ============
fig, ax = plt.subplots(figsize=(8.4, 6.4))
fams = u.model_family.astype("category")
cats = list(fams.cat.categories)
cmap = plt.get_cmap("tab10")
for k, f in enumerate(cats):
    g = u[u.model_family == f]
    ax.scatter(g.C_B_row, g.C_W, s=26, alpha=0.75, label=f,
               color=cmap(k % 10), edgecolors="none")
lim = max(abs(u.C_B_row).max(), abs(u.C_W).max()) * 1.08
ax.plot([-lim, lim], [-lim, lim], "--", color="#999", lw=0.9)
ax.axhline(0, color="#bbb", lw=0.8); ax.axvline(0, color="#bbb", lw=0.8)
ax.set_xlabel("C_B（跨靶定位的协方差信号贡献）", fontsize=10)
ax.set_ylabel("C_W（靶内共变的协方差信号贡献）", fontsize=10)
ax.set_title("图2  校准信号 C_B 与排序信号 C_W 的分解（U_CANONICAL）",
             fontsize=11.5, pad=10)
ax.legend(fontsize=7.4, ncol=2, loc="upper left", framealpha=0.9)
ax.grid(alpha=0.25)
fig.text(0.5, -0.035,
         "满足 r_full ≈ C_B + C_W（实测残差 ≤ 3.3e-16）。虚线为 C_B = C_W。\n"
         "二者按 §29 应读作「协方差信号贡献」，不是独立的标准化能力。",
         ha="center", fontsize=8.3, color="#333")
fig.savefig(FIG / "fig2_CB_CW_decomposition.png")
plt.close(fig)
print("  fig2 ✅")

# ============ fig3 split inflation ============
si = OUT / "split_inflation" / "split_inflation_raw.csv"
if si.exists():
    D = pd.read_csv(si)
    order = ["L1_random", "L2_complex", "L3_seq40", "L4_family"]
    lab = ["random\n(最宽松)", "complex", "seq40", "family-OOD\n(最严格)"]
    agg = D.groupby("level")[["overall_pearson", "macro_spearman", "C_B_row",
                              "C_W", "efficiency@10"]].median().reindex(order)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))
    ax = axes[0]
    x = np.arange(len(order))
    for c, mk in (("overall_pearson", "o-"), ("macro_spearman", "s-"),
                  ("C_B_row", "^-"), ("C_W", "d-")):
        ax.plot(x, agg[c], mk, label=c, lw=1.7, ms=6)
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.6)
    ax.axhline(0, color="#bbb", lw=0.8)
    ax.set_ylabel("中位指标值", fontsize=10)
    ax.set_title("(a) 各划分下的中位指标", fontsize=10.5)
    ax.legend(fontsize=8, framealpha=0.9); ax.grid(alpha=0.25)
    ax = axes[1]
    base = agg.loc["L1_random"]
    for c in ("overall_pearson", "macro_spearman", "efficiency@10"):
        ax.plot(x, agg[c] / base[c], "o-", label=c, lw=1.7, ms=6)
    ax.axhline(1, color="#999", ls="--", lw=0.9)
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.6)
    ax.set_yscale("log")
    ax.set_ylabel("相对 random 的倍数（log）", fontsize=10)
    ax.set_title("(b) 膨胀幅度：Overall r 掉 4.5×，Eff@10 只掉 2.5%", fontsize=10.5)
    ax.legend(fontsize=8, framealpha=0.9); ax.grid(alpha=0.25, which="both")
    fig.suptitle("图3  划分膨胀高度不均匀（同一批 7523 行、同一批 V6.2 模型）",
                 fontsize=11.5, y=1.03)
    fig.savefig(FIG / "fig3_split_inflation.png")
    plt.close(fig)
    print("  fig3 ✅")

# ============ fig4 remedy chain ============
rp = OUT / "remedy" / "remedy_sequence.csv"
if rp.exists():
    RS = pd.read_csv(rp)
    b = RS[RS.label.str.startswith("V7.3 F1")]
    if len(b):
        b = b.iloc[0]
        r5 = RS[RS.step == "R5"].copy()
        r5["d_overall"] = r5.overall_r - b.overall_r
        r5["d_macro"] = r5.macro_spearman - b.macro_spearman
        fig, ax = plt.subplots(figsize=(8.6, 6.0))
        ax.axhline(0, color="#bbb", lw=0.9); ax.axvline(0, color="#bbb", lw=0.9)
        cols = ["#2a9d8f" if d > 0.02 else "#e76f51" for d in r5.d_macro]
        ax.scatter(r5.d_overall, r5.d_macro, s=95, c=cols, zorder=3,
                   edgecolors="white", linewidths=1.1)
        for _, r in r5.iterrows():
            nm = r.label.replace("V7.3 ", "")
            ax.annotate(nm, (r.d_overall, r.d_macro), fontsize=8,
                        xytext=(5, 5), textcoords="offset points")
        ax.set_xlabel("ΔOverall Pearson（相对 F1 mutation-only）", fontsize=10)
        ax.set_ylabel("ΔMacro per-target Spearman", fontsize=10)
        ax.set_title("图4  Benchmark-guided remedy：增益落在哪一轴", fontsize=11.5)
        ax.grid(alpha=0.25)
        fig.text(0.5, -0.045,
                 "绿 = 真增益（两轴同涨）；红 = 校准幻觉或负向。\n"
                 "F2 static：Overall r +0.028 但 Macro ρ −0.0015 ⇒ 纯校准幻觉。\n"
                 "F4 mutation×partner：+0.070 且 +0.044 ⇒ 真增益。",
                 ha="center", fontsize=8.3, color="#333")
        fig.savefig(FIG / "fig4_remedy_chain.png")
        plt.close(fig)
        print("  fig4 ✅")

# ============ fig5 AXIS C floor ============
fl = OUT / "selection" / "axis_c_floor_by_target.csv"
if fl.exists():
    G = pd.read_csv(fl)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))
    ax = axes[0]
    ax.scatter(G.n, G.E_rand, s=16, alpha=0.6, label="随机地板 $E_{rand}=K/n$")
    ax.scatter(G.n, G.E_max, s=16, alpha=0.6, label="理论上限 $E_{max}$")
    ax.set_xscale("log")
    ax.set_xlabel("target 突变数 n（log）", fontsize=10)
    ax.set_ylabel("Efficiency@K 的取值", fontsize=10)
    ax.set_title("(a) 地板与上限", fontsize=10.5)
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    ax = axes[1]
    ax.scatter(G.n, G["range"], s=18, alpha=0.7, color="#e76f51")
    ax.axvline(10, color="#333", ls="--", lw=1.0)
    ax.text(10.6, G["range"].max() * 0.86, "K = 10", fontsize=8.6)
    ax.set_xscale("log")
    ax.set_xlabel("target 突变数 n（log）", fontsize=10)
    ax.set_ylabel("可区分区间 $E_{max}-E_{rand}$", fontsize=10)
    ax.set_title("(b) 模型最多能比随机好多少（n ≤ 10 时恒为 0）", fontsize=10.5)
    ax.grid(alpha=0.25)
    fig.suptitle("图5  AXIS C 的「发现地板」来自数据：40.8% 的 target 零区分空间",
                 fontsize=11.5, y=1.03)
    fig.savefig(FIG / "fig5_axis_c_floor.png")
    plt.close(fig)
    print("  fig5 ✅")

# ============ fig6 weighting sensitivity ============
wp = OUT / "weighting" / "weighting_sensitivity.csv"
if wp.exists():
    W = pd.read_csv(wp)
    SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
    cb = [f"C_B_{s}" for s in SCHEMES]; cw = [f"C_W_{s}" for s in SCHEMES]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))
    for ax, cols, nm, cc in ((axes[0], cb, "C_B（跨靶定位）", "#e76f51"),
                             (axes[1], cw, "C_W（靶内共变）", "#2a9d8f")):
        A = W[cols].to_numpy(float)
        x = np.arange(len(SCHEMES))
        for i in range(0, len(A), max(1, len(A) // 40)):
            ax.plot(x, A[i], color=cc, alpha=0.18, lw=0.9)
        ax.plot(x, np.nanmedian(A, axis=0), "k-o", lw=2.2, ms=6,
                label="中位")
        ax.axhline(0, color="#bbb", lw=0.9)
        ax.set_xticks(x); ax.set_xticklabels(SCHEMES, fontsize=9)
        rng = np.nanmax(A, axis=1) - np.nanmin(A, axis=1)
        sgn = int(((A > 0).sum(axis=1) > 0).__and__(
            (A > 0).sum(axis=1) < len(SCHEMES)).sum())
        ax.set_title(f"{nm}\n极差中位={np.nanmedian(rng):.4f}  "
                     f"符号变化模型={sgn}/{len(A)}", fontsize=10.3)
        ax.set_ylabel("值", fontsize=10); ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("图6  目标加权敏感性：C_B 敏感、C_W 稳健（预测完全不变，只改聚合权重）",
                 fontsize=11.5, y=1.04)
    fig.savefig(FIG / "fig6_weighting_sensitivity.png")
    plt.close(fig)
    print("  fig6 ✅")

# ============ fig7 capability landscape ============
fig, ax = plt.subplots(figsize=(8.8, 6.4))
uu = u.dropna(subset=["overall_pearson", "macro_spearman"])
for k, f in enumerate(cats):
    g = uu[uu.model_family == f]
    ax.scatter(g.overall_pearson, g.macro_spearman, s=32, alpha=0.78,
               label=f, color=cmap(k % 10), edgecolors="none")
ax.axhline(0, color="#bbb", lw=0.9); ax.axvline(0, color="#bbb", lw=0.9)
# 标注两个极端
hi = uu.nlargest(1, "overall_pearson").iloc[0]
hm = uu.nlargest(1, "macro_spearman").iloc[0]
for r, c in ((hi, "#d62828"), (hm, "#003049")):
    ax.annotate(r.model_name, (r.overall_pearson, r.macro_spearman),
                fontsize=8, xytext=(6, 6), textcoords="offset points", color=c)
    ax.scatter([r.overall_pearson], [r.macro_spearman], s=150, facecolors="none",
               edgecolors=c, linewidths=1.8, zorder=4)
ax.set_xlabel("Overall Pearson（校准轴 A）", fontsize=10)
ax.set_ylabel("Macro per-target Spearman（排序轴 B）", fontsize=10)
ax.set_title("图7  能力景观：两个轴的模型排序几乎不相关（ρ = −0.26）",
             fontsize=11.5, pad=10)
ax.legend(fontsize=7.3, ncol=2, loc="upper left", framealpha=0.9)
ax.grid(alpha=0.25)
fig.text(0.5, -0.045,
         "红圈 = Overall Pearson 最高者；蓝圈 = Macro ρ 最高者。二者是不同模型。\n"
         "跨三个独立模型子集，ρ(Overall r, Macro ρ) 一致落在 −0.20 ~ −0.26。",
         ha="center", fontsize=8.3, color="#333")
fig.savefig(FIG / "fig7_capability_landscape.png")
plt.close(fig)
print("  fig7 ✅")

# ============ §80 八张表 ============
print("\n=== §80 生成表格 ===")
T = {}
T["T1_model_inventory"] = pd.read_csv(OUT / "model_registry.csv")
T["T2_headline_scoreboard"] = u.nlargest(30, "macro_spearman")[
    ["model_name", "model_family", "n_rows", "n_targets_ge5", "overall_pearson",
     "rmse", "macro_spearman", "C_B_row", "C_W", "efficiency@10", "ndcg@10"]]
T["T3_rank_agreement"] = RA
T["T4_CB_CW_decomposition"] = u[["model_name", "model_family", "overall_pearson",
                                 "C_B_row", "C_W_row", "C_B_equal",
                                 "decomp_residual"]].copy() if "decomp_residual" in u else \
    u[["model_name", "model_family", "overall_pearson", "C_B_row", "C_W_row",
       "C_B_equal"]]
T["T5_split_inflation"] = pd.read_csv(
    OUT / "split_inflation" / "split_inflation_raw.csv") if (
    OUT / "split_inflation" / "split_inflation_raw.csv").exists() else pd.DataFrame()
T["T6_weighting_sensitivity"] = pd.read_csv(
    OUT / "weighting" / "weighting_sensitivity.csv") if (
    OUT / "weighting" / "weighting_sensitivity.csv").exists() else pd.DataFrame()
T["T7_axis_c_floor"] = pd.read_csv(
    OUT / "selection" / "axis_c_floor_by_target.csv") if (
    OUT / "selection" / "axis_c_floor_by_target.csv").exists() else pd.DataFrame()
T["T8_scenarios"] = pd.read_csv(OUT / "scenarios" / "scenario_recommendations.csv") \
    if (OUT / "scenarios" / "scenario_recommendations.csv").exists() else pd.DataFrame()

for k, v in T.items():
    if isinstance(v, pd.DataFrame) and len(v):
        v.to_csv(TAB / f"{k}.csv", index=False, encoding="utf-8")
        print(f"  {k}: {v.shape[0]} 行 × {v.shape[1]} 列")
    else:
        print(f"  {k}: 空（跳过）")

# 表清单
idx = ["# V8 表格清单（§80）\n",
       "| 表 | 文件 | 行数 | 说明 |", "|---|---|---:|---|"]
DESC = {
    "T1_model_inventory": "模型清单（211 预测键：连接方式/覆盖/证据层级/溯源）",
    "T2_headline_scoreboard": "U_CANONICAL 上按 Macro ρ 前 30 的三轴读数",
    "T3_rank_agreement": "跨轴模型排名秩一致性（Spearman / Kendall）",
    "T4_CB_CW_decomposition": "C_B / C_W 协方差信号分解（含恒等式残差 QC）",
    "T5_split_inflation": "四种划分下每个模型的全部指标（36 行）",
    "T6_weighting_sensitivity": "五种加权方案下每个模型的 C_B / C_W / Pearson",
    "T7_axis_c_floor": "每个 target 的随机地板 / 理论上限 / 可区分区间",
    "T8_scenarios": "§41 三场景的首选模型与惯例对照",
}
for k, v in T.items():
    n = v.shape[0] if isinstance(v, pd.DataFrame) else 0
    idx.append(f"| {k.split('_')[0]} | `tables/{k}.csv` | {n} | {DESC.get(k,'')} |")
(TAB / "TABLE_INDEX.md").write_text("\n".join(idx) + "\n", encoding="utf-8")
print(f"\n表格索引写出 -> tables/TABLE_INDEX.md")
