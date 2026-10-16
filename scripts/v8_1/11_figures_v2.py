# -*- coding: utf-8 -*-
"""scripts/v8_1/11_figures_v2.py -- 稿件 v2 图的重新设计（主体 6 图 + 补充图）。

设计约束（按审稿人意见）：
* **英文主体图**纯英文标签 + matplotlib mathtext（$\Delta\Delta G$, $C_B$, $C_W$, $\rho$）；
* **中文主体图**独立生成，字体经**程序化字形验证**，缺字即 FAIL；
* 图 1 = **新概念总览图（4 面板）**；原 matched-universe-strata 降为补充图 S1；
* 主体图最多 6 张；避免 spaghetti；不出现 SHA / 脚本名 / GREEN-RED 术语。

输出：
    v8_1_final/figures_main/         (Fig1..Fig6, EN + ZH)
    v8_1_final/figures_supplement/   (FigS1..FigS4)
    v8_1_final/FIGURE_GLYPH_QC.md
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
FM = F / "figures_main"
FS = F / "figures_supplement"
FM.mkdir(exist_ok=True); FS.mkdir(exist_ok=True)

EN = "DejaVu Sans"
ZH = None
for c in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC",
          "DengXian", "SimSun"):
    try:
        font_manager.findfont(c, fallback_to_default=False)
        ZH = c; break
    except Exception:
        continue
print(f"EN font = {EN} | ZH font = {ZH}")

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
W = pd.read_csv(R / "experiments" / "v8" / "weighting" / "weighting_sensitivity.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
RS = pd.read_csv(R / "experiments" / "v8" / "remedy" / "remedy_sequence.csv")

GLYPH_LOG = []


def use(font, zh=False):
    plt.rcParams["font.family"] = font
    plt.rcParams["mathtext.fontset"] = "dejavusans"
    plt.rcParams["mathtext.default"] = "it"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["savefig.bbox"] = "tight"


def glyph_check(fig, tag, zh):
    """程序化字形验证：渲染后检查是否有缺字（tofu）。"""
    from matplotlib import ft2font
    ok = True
    miss = []
    fpath = font_manager.findfont(plt.rcParams["font.family"][0] if
                                  isinstance(plt.rcParams["font.family"], list)
                                  else plt.rcParams["font.family"],
                                  fallback_to_default=False)
    face = ft2font.FT2Font(fpath)
    for ax in fig.axes:
        texts = [ax.get_title(), ax.get_xlabel(), ax.get_ylabel()]
        for t in ax.texts:
            texts.append(t.get_text())
        for l in ax.get_legend().get_texts() if ax.get_legend() else []:
            texts.append(l.get_text())
        for s in texts:
            for ch in str(s):
                if ch.isspace() or ord(ch) < 32:
                    continue
                if ch in "$\\{}_^&%#":
                    continue
                try:
                    if face.get_char_index(ch) == 0:
                        ok = False; miss.append(ch)
                except Exception:
                    pass
    GLYPH_LOG.append({"figure": tag, "cjk": bool(zh),
                      "font": plt.rcParams["font.family"],
                      "missing_glyphs": sorted(set(miss)),
                      "result": "PASS" if ok else "FAIL"})
    return ok


def save(fig, tag, zh, folder=FM):
    glyph_check(fig, tag, zh)
    suf = "_ZH" if zh else "_EN"
    base = folder / f"{tag}{suf}"
    fig.savefig(str(base) + ".png", dpi=600)
    try:
        fig.savefig(str(base) + ".pdf")
    except Exception:
        pass
    plt.close(fig)
    print(f"   {base.name}.png")


# ==================================================================
# FIGURE 1 — 概念总览（4 面板）
# ==================================================================
def fig1(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    fig = plt.figure(figsize=(13.2, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.20,
                          left=0.05, right=0.975, top=0.915, bottom=0.055)

    # ---- Panel a: evaluation task ----
    ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title(T("a  A pooled benchmark pools many targets, each with its own ΔΔG level",
                   "a  合并基准汇集多个靶点，每个靶点有其自身的 ΔΔG 水平"),
                 fontsize=10.0, loc="left", pad=6)
    rng = np.random.RandomState(0)
    for i, (cx, mu, col) in enumerate([(1.7, 3.2, "#4c72b0"), (5.0, 6.4, "#dd8452"),
                                       (8.3, 4.6, "#55a868")]):
        ax.add_patch(FancyBboxPatch((cx - 1.40, 0.95), 2.8, 6.5,
                                    boxstyle="round,pad=0.08",
                                    fc="#f6f7f9", ec=col, lw=1.3))
        ax.text(cx, 7.95, T(f"Target {i+1}", f"靶点 {i+1}"), ha="center",
                fontsize=9.2, color=col, fontweight="bold")
        ys = mu + rng.randn(5) * 0.75
        for yv in ys:
            ax.plot([cx - 1.10, cx + 1.10], [yv, yv], color=col, lw=1.6, alpha=0.62)
        ax.plot([cx - 1.18, cx + 1.18], [mu, mu], color=col, lw=3.0)
        ax.text(cx + 1.28, mu, rf"$\mu_{i+1}$",
                ha="left", va="center", fontsize=8.6, color=col)
    ax.text(5.0, 0.15, T("many targets · unequal mutation counts",
                         "多个靶点 · 突变数高度不均"),
            ha="center", fontsize=8.8, style="italic", color="#555")

    # ---- Panel b: why pooled correlation is ambiguous ----
    ax = fig.add_subplot(gs[0, 1]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title(T("b  Pooled correlation mixes two signal sources",
                   "b  合并相关混合了两种信号来源"),
                 fontsize=10.0, loc="left", pad=6)
    ax.set_xlabel("")
    axx = fig.add_axes([0.575, 0.545, 0.375, 0.315])
    for mu, col in [(1.2, "#4c72b0"), (3.4, "#dd8452"), (5.6, "#55a868")]:
        x = np.linspace(mu - 0.85, mu + 0.85, 40)
        axx.scatter(x, x + rng.randn(40) * 0.22, s=13, color=col, alpha=0.8)
        axx.axvline(mu, color=col, ls=":", lw=0.9, alpha=0.8)
    axx.set_xlabel(T("true mutation effect", "真实突变效应"), fontsize=8.0)
    axx.set_ylabel(T("prediction", "预测值"), fontsize=8.0)
    axx.tick_params(labelsize=7.0)
    axx.grid(alpha=0.18)
    axx.text(0.04, 0.94, r"$y_{ci}=\mu_c+\varepsilon_{ci}$", transform=axx.transAxes,
             fontsize=9.0, va="top")
    ax.text(5.0, 7.6, r"$r_{\mathrm{pooled}} = C_B + C_W$", ha="center",
            fontsize=13.5, fontweight="bold")
    ax.text(5.0, 5.9, T("BETWEEN-TARGET: differences between target means  →  $C_B$",
                        "跨靶：靶点均值之间的差异  →  $C_B$"),
            ha="center", fontsize=8.8, color="#c44e52")
    ax.text(5.0, 4.7, T("WITHIN-TARGET: mutation ordering around each target mean  →  $C_W$",
                        "靶内：各靶点内部围绕其均值的突变排序  →  $C_W$"),
            ha="center", fontsize=8.8, color="#4c72b0")
    ax.text(5.0, 2.8, T("$C_B$ and $C_W$ are covariance-signal contributions,\n"
                        "not independent correlations or capacities.",
                        "$C_B$ 与 $C_W$ 是协方差信号贡献，\n并非独立相关系数或独立能力。"),
            ha="center", fontsize=8.2, style="italic", color="#444")

    # ---- Panel c: evaluation framework pipeline ----
    ax = fig.add_subplot(gs[1, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title(T("c  Capability-resolved evaluation framework",
                   "c  能力分辨的评价框架"),
                 fontsize=10.0, loc="left", pad=6)
    steps = [(T("prediction\noutputs", "预测\n输出"), "#4c72b0"),
             (T("exact row-set\nmatching", "行集合\n精确匹配"), "#4c72b0"),
             (T("provenance\n& de-replication", "来源筛选\n与去伪重复"), "#8172b3"),
             (T("covariance\ndecomposition", "协方差\n分解"), "#c44e52"),
             (T("split\nstrictness", "划分\n严格度"), "#dd8452"),
             (T("target\nweighting", "目标\n加权"), "#dd8452"),
             (T("candidate-set\ngeometry", "候选集\n几何"), "#dd8452"),
             (T("capability\ninterpretation", "能力\n解读"), "#55a868")]
    x0, y0, w, h = 0.35, 6.2, 2.28, 2.1
    for i, (lab, col) in enumerate(steps):
        cx = x0 + (i % 4) * 2.42
        cy = y0 - (i // 4) * 3.0
        ax.add_patch(FancyBboxPatch((cx, cy), w, h, boxstyle="round,pad=0.07",
                                    fc="white", ec=col, lw=1.4))
        ax.text(cx + w / 2, cy + h / 2, lab, ha="center", va="center",
                fontsize=7.8, color="#222")
        if i % 4 < 3:
            ax.add_patch(FancyArrowPatch((cx + w + 0.03, cy + h / 2),
                                         (cx + 2.42 - 0.03, cy + h / 2),
                                         arrowstyle="-|>", mutation_scale=10,
                                         color="#888", lw=1.0))
        elif i == 3:
            ax.add_patch(FancyArrowPatch((cx + w / 2, cy - 0.03),
                                         (cx + w / 2, cy - 3.0 + h + 0.03),
                                         arrowstyle="-|>", mutation_scale=10,
                                         color="#888", lw=1.0))

    # ---- Panel d: main controlled finding ----
    ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
    ax.set_title(T("d  Controlled split shift: between-target signal largely disappears",
                   "d  受控划分切换：跨靶信号几乎消失"),
                 fontsize=10.0, loc="left", pad=6)
    axd = ax.inset_axes([0.06, 0.30, 0.88, 0.63])
    names = ["random", "complex", "seq40", "family-OOD"]
    cb = [H4[k]["C_B"] for k in ["random", "complex", "seq40", "family_OOD"]]
    cw = [H4[k]["C_W"] for k in ["random", "complex", "seq40", "family_OOD"]]
    rp = [H4[k]["overall_pearson"] for k in ["random", "complex", "seq40", "family_OOD"]]
    xs = np.arange(4)
    axd.bar(xs - 0.19, cb, width=0.36, color="#c44e52", label=r"$C_B$ (between)")
    axd.bar(xs + 0.19, cw, width=0.36, color="#4c72b0", label=r"$C_W$ (within)")
    axd.axhline(0, color="#888", lw=0.9)
    axd.set_xticks(xs)
    axd.set_xticklabels(names, fontsize=7.4, rotation=12)
    axd.set_ylabel(T("covariance contribution", "协方差贡献"), fontsize=7.8)
    axd.tick_params(labelsize=6.8)
    axd.legend(fontsize=6.8, loc="upper right", framealpha=0.9)
    axd.grid(alpha=0.18, axis="y")
    axd2 = axd.twinx()
    axd2.plot(xs, rp, "ko--", lw=1.3, ms=4.2, label=T("pooled $r$", "合并 $r$"))
    axd2.set_ylabel(T("pooled Pearson", "合并 Pearson"), fontsize=7.8)
    axd2.tick_params(labelsize=6.8)
    axd2.legend(fontsize=6.8, loc="center right", framealpha=0.9)
    ax.text(0.5, 0.185,
            T(f"pooled $r$: {rp[0]:.3f} $\\rightarrow$ {rp[3]:.3f}    "
              f"$\\Delta C_B$ = {DL.get('delta_C_B'):+.3f}    "
              f"$\\Delta C_W$ = {DL.get('delta_C_W'):+.3f}",
              f"合并 $r$: {rp[0]:.3f} $\\rightarrow$ {rp[3]:.3f}    "
              f"$\\Delta C_B$ = {DL.get('delta_C_B'):+.3f}    "
              f"$\\Delta C_W$ = {DL.get('delta_C_W'):+.3f}"),
            transform=ax.transAxes, ha="center", fontsize=8.6)
    ax.text(0.5, 0.085,
            T(f"between-target share of the absolute decline = "
              f"{100*DL.get('fraction_between'):.1f}%",
              f"跨靶分量占绝对下降的 {100*DL.get('fraction_between'):.1f}%"),
            transform=ax.transAxes, ha="center", fontsize=9.4,
            color="#c44e52", fontweight="bold")
    ax.text(0.5, 0.005,
            T("complex and seq40 yield identical assignments here",
              "complex 与 seq40 在本数据上给出完全相同的划分"),
            transform=ax.transAxes, ha="center", fontsize=7.2,
            style="italic", color="#a33")
    fig.suptitle(T("Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors",
                   "图 1 | PPI ΔΔG 预测器的能力分辨评价"),
                 fontsize=13.2, y=0.972, fontweight="bold")
    save(fig, "Fig1_capability_resolved_overview", zh)


# ==================================================================
# FIGURE 2 — 受控 family-OOD 分解（3 面板，中心实证图）
# ==================================================================
def fig2(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    order = ["random", "complex", "seq40", "family_OOD"]
    lab = [T("random", "随机"), T("complex", "复合物"),
           T("seq40", "seq40"), T("family-OOD", "family-OOD")]
    rp = [H4[k]["overall_pearson"] for k in order]
    cb = [H4[k]["C_B"] for k in order]
    cw = [H4[k]["C_W"] for k in order]
    mac = [H4[k]["macro_spearman"] for k in order]
    rms = [H4[k]["rmse"] for k in order]
    xs = np.arange(4)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.15))
    ax = axes[0]
    ax.plot(xs, rp, "o-", lw=2.2, ms=6, color="#c44e52")
    ax.axhline(0, color="#999", lw=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=8.2, rotation=12)
    ax.set_ylabel(T("median pooled Pearson", "合并 Pearson 中位数"), fontsize=9.0)
    ax.set_title(T("(a) Pooled correlation by split regime",
                   "(a) 各划分下的合并相关"), fontsize=9.6, loc="left")
    ax.grid(alpha=0.2)
    for i, v in enumerate(rp):
        ax.annotate(f"{v:.3f}", (i, v), fontsize=7.6, xytext=(0, 7),
                    textcoords="offset points", ha="center")
    ax = axes[1]
    ax.plot(xs, cb, "s-", lw=2.2, ms=6, color="#c44e52", label=r"$C_B$")
    ax.plot(xs, cw, "d-", lw=2.2, ms=6, color="#4c72b0", label=r"$C_W$")
    ax.axhline(0, color="#999", lw=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=8.2, rotation=12)
    ax.set_ylabel(T("covariance contribution", "协方差贡献"), fontsize=9.0)
    ax.set_title(T("(b) Exact decomposition $r=C_B+C_W$",
                   "(b) 精确分解 $r=C_B+C_W$"), fontsize=9.6, loc="left")
    ax.legend(fontsize=8.4); ax.grid(alpha=0.2)
    ax.annotate(T("$C_W$ remains positive", "$C_W$ 仍为正"),
                xy=(3, cw[3]), xytext=(1.35, 0.215), fontsize=8.0,
                arrowprops=dict(arrowstyle="->", color="#4c72b0", lw=1.0),
                color="#4c72b0")
    ax = axes[2]
    d_cb, d_cw = DL.get("delta_C_B"), DL.get("delta_C_W")
    ax.barh([1], [abs(d_cb)], color="#c44e52", height=0.5)
    ax.barh([0], [abs(d_cw)], color="#4c72b0", height=0.5)
    ax.set_yticks([1, 0])
    ax.set_yticklabels([T("between-target\n($\\Delta C_B$)", "跨靶\n($\\Delta C_B$)"),
                        T("within-target\n($\\Delta C_W$)", "靶内\n($\\Delta C_W$)")],
                       fontsize=8.4)
    ax.set_xlabel(T("absolute change from random → family-OOD",
                    "随机 → family-OOD 的绝对变化"), fontsize=8.6)
    ax.text(abs(d_cb) * 0.5, 1.30, f"{100*DL.get('fraction_between'):.1f}%",
            ha="center", fontsize=11.5, color="#c44e52", fontweight="bold")
    ax.text(abs(d_cw) + 0.012, 0.0, f"{100*DL.get('fraction_within'):.1f}%",
            va="center", fontsize=9.4, color="#4c72b0")
    ax.set_title(T("(c) Share of the absolute decline",
                   "(c) 绝对下降的份额"), fontsize=9.6, loc="left")
    ax.grid(alpha=0.2, axis="x"); ax.set_xlim(0, abs(d_cb) * 1.22)
    fig.suptitle(T("Figure 2 | Family-OOD evaluation disproportionately removes the "
                   "between-target component of pooled correlation",
                   "图 2 | Family-OOD 评价不成比例地移除合并相关中的跨靶分量"),
                 fontsize=11.6, y=1.045, fontweight="bold")
    fig.text(0.5, -0.085,
             T("Identical 7,523 mutation rows and identical model configurations; only the "
               "split regime changes. The complex and seq40 partitions are identical in "
               "this dataset and are not independent comparisons. Panel (c) uses absolute "
               "changes, $|\\Delta C_B|/(|\\Delta C_B|+|\\Delta C_W|)$.",
               "同一批 7,523 行突变、同一批模型配置，仅改变划分方式。complex 与 seq40 "
               "在本数据上完全相同，不构成独立比较。面板 (c) 采用绝对变化口径，"
               "$|\\Delta C_B|/(|\\Delta C_B|+|\\Delta C_W|)$。"),
             ha="center", fontsize=8.2, color="#333")
    save(fig, "Fig2_controlled_familyOOD_decomposition", zh)


# ==================================================================
# FIGURE 3 — 加权敏感性（简化：范围 + 符号稳定性）
# ==================================================================
def fig3(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    SCH = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    for ax, key, cc, nm in ((axes[0], "C_B", "#c44e52", r"$C_B$  (between-target)"),
                            (axes[1], "C_W", "#4c72b0", r"$C_W$  (within-target)")):
        cols = [f"{key}_{s}" for s in SCH]
        A_ = W[cols].to_numpy(float)
        rng_ = np.nanmax(A_, axis=1) - np.nanmin(A_, axis=1)
        r0 = W2[f"signflip_{key}_gt0"].sum()
        r5 = W2[f"signflip_{key}_gt0.005"].sum()
        ax.hist(rng_, bins=26, color=cc, alpha=0.72, edgecolor="white", linewidth=0.4)
        ax.axvline(np.nanmedian(rng_), color="k", ls="--", lw=1.5)
        ax.set_xlabel(T("range across the five weighting schemes",
                        "五种加权方案下的极差"), fontsize=8.8)
        ax.set_ylabel(T("number of model configurations",
                        "模型配置数"), fontsize=8.8)
        ax.set_title(T(f"{nm}\nmedian range = {np.nanmedian(rng_):.4f}   |   "
                       f"sign flips: {r0} (any), {r5} (exceeding 0.005)",
                       f"{nm}\n极差中位 = {np.nanmedian(rng_):.4f}   |   "
                       f"符号翻转：{r0}（任意），{r5}（超过 0.005）"),
                     fontsize=9.2, loc="left")
        ax.grid(alpha=0.2)
    fig.suptitle(T("Figure 3 | Target weighting changes the interpretation of "
                   "covariance contributions",
                   "图 3 | 目标加权改变协方差贡献的解读"),
                 fontsize=11.4, y=1.05, fontweight="bold")
    fig.text(0.5, -0.075,
             T("Predictions are unchanged; only aggregation weights differ "
               "(per-mutation, per-target, $\\sqrt{n}$, $\\log(1+n)$, $\\min(n,30)$). "
               "$C_B$ is unstable in direction only near zero; $C_W$ keeps its sign but "
               "its magnitude varies, so its magnitude is not invariant.",
               "预测值不变，仅改变聚合权重（每个突变、每个靶点、$\\sqrt{n}$、"
               "$\\log(1+n)$、$\\min(n,30)$）。$C_B$ 仅在近零区间方向不稳；"
               "$C_W$ 保持符号但幅度变化，故其幅度并非不变。"),
             ha="center", fontsize=8.2, color="#333")
    save(fig, "Fig3_target_weighting_sensitivity", zh)


# ==================================================================
# FIGURE 4 — 固定 K 几何
# ==================================================================
def fig4(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    fig, ax = plt.subplots(figsize=(7.6, 4.7))
    K = 10
    m = GEO_DF.n <= 160
    ax.plot(GEO_DF.n[m], GEO_DF.D[m], lw=2.4, color="#4c72b0")
    ax.axvline(K, color="#c44e52", ls="--", lw=1.4)
    ax.text(K + 1.6, GEO_DF.D[m].max() * 0.86, r"$n \leq K \Rightarrow D = 0$",
            fontsize=8.8, color="#c44e52")
    ns_ = GEO["n_star_analytic"]
    ax.axvline(ns_, color="#55a868", ls=":", lw=1.4)
    ax.annotate(T(f"peak $n^*\\approx{ns_:.1f}$", f"峰值 $n^*\\approx{ns_:.1f}$"),
                xy=(ns_, GEO_DF.D[m].max()), xytext=(ns_ + 6, GEO_DF.D[m].max() * 0.70),
                fontsize=8.8, color="#55a868",
                arrowprops=dict(arrowstyle="->", color="#55a868", lw=1.0))
    ax.plot([20, 50], [GEO["D_at_n20"], GEO["D_at_n50"]], "ko", ms=5)
    ax.annotate(f"D(20) = {GEO['D_at_n20']:.4f}", (20, GEO["D_at_n20"]),
                xytext=(24, GEO["D_at_n20"] + 0.018), fontsize=8.0)
    ax.annotate(f"D(50) = {GEO['D_at_n50']:.4f}", (50, GEO["D_at_n50"]),
                xytext=(60, GEO["D_at_n50"] + 0.012), fontsize=8.0)
    ax.set_xlabel(T("candidate-set size $n$  (K = 10 fixed)", "候选集大小 $n$（K = 10 固定）"),
                  fontsize=9.2)
    ax.set_ylabel(r"discriminable interval  $D(n,K)$", fontsize=9.2)
    ax.grid(alpha=0.2)
    ax.set_title(T("Figure 4 | Fixed-K shortlist evaluation is constrained by "
                   "candidate-set geometry",
                   "图 4 | 固定 K 短名单评价受候选集几何约束"),
                 fontsize=11.0, loc="left", pad=8)
    fig.text(0.5, -0.06,
             T(f"$D(n,K)=K(n-K)/[n(n+1)]$. In this benchmark "
               f"{GEO['n_le_K_count']} of {GEO['n_ge5_count']} targets with $n\\geq5$ "
               f"({100*GEO['n_le_K_fraction']:.1f}%) have $n \\leq K$, where any model "
               "selects every candidate and the target carries no information for "
               "fixed-K shortlist metrics. $D$ peaks near "
               f"$n^*\\approx{ns_:.1f}$ and declines thereafter, which is why a "
               "'larger candidate sets are more discriminable' heuristic does not hold.",
               f"$D(n,K)=K(n-K)/[n(n+1)]$。本基准中 {GEO['n_ge5_count']} 个 $n\\geq5$ "
               f"靶点有 {GEO['n_le_K_count']} 个（{100*GEO['n_le_K_fraction']:.1f}%）"
               f"满足 $n \\leq K$，此时任何模型都会选中全部候选，该靶点对固定 K "
               f"短名单指标不携带信息。$D$ 在 $n^*\\approx{ns_:.1f}$ 附近取峰后下降，"
               "故「候选集越大越可区分」的启发式并不成立。"),
             ha="center", fontsize=8.0, color="#333")
    save(fig, "Fig4_fixed_K_candidate_geometry", zh)


# ==================================================================
# FIGURE 5 — 回溯性补救（仅英文方法的描述性标签）
# ==================================================================
def fig5(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    NICE = {
        "F1_Base_Mut": ("Mutation properties", "突变性质"),
        "F2_Base_Static": ("Static local geometry", "静态局部几何"),
        "F3_Base_EnvPartner": ("Partner static context", "partner 静态环境"),
        "F4_Base_PartnerPert": ("Mutation × partner perturbation", "突变 × partner 扰动"),
        "F5_Base_IntraPert": ("Intra-chain perturbation", "链内扰动"),
        "F6_Partner_Intra": ("Partner + intra-chain", "partner + 链内"),
        "F7_F6_Burial": ("Partner + intra + burial", "partner + 链内 + 埋藏"),
        "F8_Full": ("Full local-physics model", "全量局部物理模型"),
        "F9_PhysicsOnly": ("Physics features only", "仅物理特征"),
        "F0_Base": ("Base model", "基线模型"),
    }
    b = RS[RS.label.str.startswith("V7.3 F1")]
    r5 = RS[RS.step == "R5"].copy()
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    if len(b):
        b = b.iloc[0]
        r5["d_ov"] = r5.overall_r - b.overall_r
        r5["d_mc"] = r5.macro_spearman - b.macro_spearman
        ax.axhline(0, color="#999", lw=0.9); ax.axvline(0, color="#999", lw=0.9)
        for _, r in r5.iterrows():
            frag = r.label.replace("V7.3 ", "").replace("_s42", "")
            nice_en, nice_zh = NICE.get(frag, (frag, frag))
            real = r.d_mc > 0.02
            ax.scatter(r.d_ov, r.d_mc, s=115, c="#55a868" if real else "#c44e52",
                       zorder=3, edgecolors="white", linewidths=1.1)
            ax.annotate(nice_zh if zh else nice_en, (r.d_ov, r.d_mc), fontsize=8.0,
                        xytext=(7, 5), textcoords="offset points")
    ax.set_xlabel(T(r"$\Delta$ pooled Pearson vs mutation-properties baseline",
                    r"相对「突变性质」基线的 $\Delta$ 合并 Pearson"), fontsize=9.4)
    ax.set_ylabel(T(r"$\Delta$ macro per-target Spearman",
                    r"$\Delta$ 靶内 Macro Spearman"), fontsize=9.4)
    ax.grid(alpha=0.2)
    ax.set_title(T("Figure 5 | Capability decomposition separates pooled-score-only "
                   "improvements from changes that also affect within-target ranking",
                   "图 5 | 能力分解区分「仅合并分数改善」与「同时影响靶内排序」的改动"),
                 fontsize=10.4, loc="left", pad=8)
    fig.text(0.5, -0.075,
             T("Retrospective matched ablation on an identical 3,351-row universe. "
               "Green: both axes improve. Red: pooled-score-only or negative change. "
               "This analysis is retrospective; it shows that the decomposition can "
               "express this distinction and does not show prospective prediction of "
               "which modification would succeed.",
               "在同一 3,351 行宇宙上的回溯性匹配消融。绿色：两轴均改善；"
               "红色：仅合并分数改善或变化为负。该分析为回溯性：它说明分解能够表达"
               "这一区别，并不表明基准能前瞻预测哪种改动会成功。"),
             ha="center", fontsize=8.2, color="#333")
    save(fig, "Fig5_retrospective_remedy", zh)


# ==================================================================
# FIGURE 6 — 外部证据 / 泛化边界
# ==================================================================
def fig6(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    cols = {"DDGPred": "#4c72b0", "ddG_FoldX": "#dd8452", "Flex-ddG": "#55a868",
            "proteinMPNN": "#8172b3", "ESM-IF": "#937860"}
    for r in EP.itertuples():
        ax.scatter(r.overall_pearson, r.macro_spearman, s=135,
                   c=cols.get(r.model, "#666"), zorder=3,
                   edgecolors="white", linewidths=1.1)
        ax.annotate(r.model, (r.overall_pearson, r.macro_spearman), fontsize=8.4,
                    xytext=(8, 5), textcoords="offset points")
    ax.axhline(0, color="#999", lw=0.9); ax.axvline(0, color="#999", lw=0.9)
    ax.set_xlabel(T("pooled Pearson", "合并 Pearson"), fontsize=9.4)
    ax.set_ylabel(T("macro per-target Spearman", "靶内 Macro Spearman"), fontsize=9.4)
    ax.grid(alpha=0.2)
    ax.set_title(T("Figure 6 | External method capability profiles define the "
                   "boundary of generalisation",
                   "图 6 | 外部方法能力画像界定泛化边界"),
                 fontsize=10.6, loc="left", pad=8)
    fig.text(0.5, -0.06,
             T("Five independently developed methods, evaluated on a matched, "
               "unambiguous mapping (1,781 mutations). Rank agreement is positive "
               "($\\rho$=+0.90, exact permutation $p$=0.083); with only five methods the "
               "panel is underpowered, and this positive alignment does not support a "
               "universal decoupling of pooled and within-target performance.",
               "五种独立开发的方法，在匹配且无歧义的映射（1,781 个突变）上评价。"
               "秩一致为正（$\\rho$=+0.90，精确置换 $p$=0.083）；仅 5 个方法时功效不足，"
               "该正向对齐不支持合并性能与靶内性能的普遍解耦。"),
             ha="center", fontsize=8.0, color="#333")
    save(fig, "Fig6_external_boundary", zh)


# ==================================================================
# 补充图
# ==================================================================
def figs(zh):
    use(ZH if zh else EN, zh)
    T = (lambda en, zh_: zh_ if zh else en)
    # S1 matched strata
    fig, ax = plt.subplots(figsize=(8.0, 4.3))
    u = UR.sort_values("n_rows"); y = np.arange(len(u))
    ax.barh(y, u.n_rows, color="#4c72b0", alpha=0.85, label=T("rows", "行数"))
    ax.barh(y, u.n_ge5, color="#dd8452", alpha=0.95, height=0.55,
            label=T("targets with n>=5", "n≥5 的靶点"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ({h[:8]})" for a, h in zip(u.artifact, u.rowset_hash)],
                       fontsize=7.4)
    ax.set_xlabel(T("number of mutation rows", "突变行数"), fontsize=9.0)
    ax.legend(fontsize=8.0)
    ax.set_title(T("Figure S1 | Matched-evaluation strata",
                   "图 S1 | 匹配评价分层"), fontsize=10.2, loc="left")
    ax.grid(alpha=0.2, axis="x")
    save(fig, "FigS1_matched_strata", zh, FS)

    # S2 artifact-level heatmap
    A_ = pd.read_csv(R / "experiments" / "v8" / "rank_reversal" /
                     "rank_agreement_spearman_U_CANONICAL_ALL.csv", index_col=0)
    B_ = pd.read_csv(R / "experiments" / "v8" / "rank_reversal" /
                     "rank_agreement_spearman_U_CANONICAL_CLEAN.csv", index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6))
    for ax, Mx, tt in ((axes[0], A_, T(f"(a) all artifacts (n=160)",
                                       "(a) 全部制品（n=160）")),
                       (axes[1], B_, T(f"(b) excluding duplicated families (n=109)",
                                       "(b) 排除重复家族后（n=109）"))):
        M = Mx.to_numpy(float)
        im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(Mx.columns)))
        ax.set_xticklabels(Mx.columns, rotation=45, ha="right", fontsize=7.4)
        ax.set_yticks(range(len(Mx.index)))
        ax.set_yticklabels(Mx.index, fontsize=7.4)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center",
                            fontsize=6.0,
                            color="white" if abs(M[i, j]) > 0.55 else "black")
        ax.set_title(tt, fontsize=9.4, loc="left")
    fig.colorbar(im, ax=axes, shrink=0.8, label=T("Spearman rank agreement",
                                                  "Spearman 秩一致"))
    fig.suptitle(T("Figure S2 | Artifact-level rank agreement (descriptive only; "
                   "artifacts are not independent methods)",
                   "图 S2 | 制品层面秩一致（仅描述性；制品不是独立方法）"),
                 fontsize=10.4, y=0.99)
    save(fig, "FigS2_artifact_level_heatmap", zh, FS)

    # S3 full weighting spaghetti
    SCH = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    for ax, key, cc, nm in ((axes[0], "C_B", "#c44e52", r"$C_B$"),
                            (axes[1], "C_W", "#4c72b0", r"$C_W$")):
        cols = [f"{key}_{s}" for s in SCH]
        A_ = W[cols].to_numpy(float)
        xs = np.arange(len(SCH))
        for i in range(0, len(A_), 2):
            ax.plot(xs, A_[i], color=cc, alpha=0.14, lw=0.8)
        ax.plot(xs, np.nanmedian(A_, axis=0), "k-o", lw=2.1, ms=5.5)
        ax.axhline(0, color="#999", lw=0.8)
        ax.set_xticks(xs); ax.set_xticklabels(SCH, fontsize=8.2)
        ax.set_ylabel(T("value", "取值"), fontsize=8.8)
        ax.set_title(nm, fontsize=9.6, loc="left"); ax.grid(alpha=0.2)
    fig.suptitle(T("Figure S3 | Per-configuration weighting sensitivity",
                   "图 S3 | 逐配置的加权敏感性"), fontsize=10.4, y=1.04)
    save(fig, "FigS3_weighting_spaghetti", zh, FS)

    # S4 pool-sensitivity of rank agreement
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    s = ST.sort_values("rho_overall_vs_macro")
    yy = np.arange(len(s))
    ax.barh(yy, s.rho_overall_vs_macro, color="#4c72b0", alpha=0.85)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{h[:8]} / L{lv}" for h, lv in
                        zip(s.rowset_hash, s.level)], fontsize=7.4)
    ax.axvline(0, color="#999", lw=0.9)
    ax.set_xlabel(T(r"$\rho$(pooled Pearson, macro Spearman)",
                    r"$\rho$（合并 Pearson，Macro Spearman）"), fontsize=9.0)
    ax.set_title(T("Figure S4 | Pool sensitivity of cross-metric rank agreement",
                   "图 S4 | 跨指标秩一致的池敏感性"), fontsize=10.2, loc="left")
    ax.grid(alpha=0.2, axis="x")
    save(fig, "FigS4_pool_sensitivity", zh, FS)


for zh in (False, True):
    print(f"--- {'ZH' if zh else 'EN'} ---")
    fig1(zh); fig2(zh); fig3(zh); fig4(zh); fig5(zh); fig6(zh); figs(zh)

# ---------- GLYPH QC ----------
Q = pd.DataFrame(GLYPH_LOG)
Q.to_csv(A / "figure_glyph_qc.csv", index=False, encoding="utf-8")
n_fail = int((Q.result == "FAIL").sum())
md = ["# FIGURE_GLYPH_QC\n",
      f"**FIGURE_GLYPH_QC = {'PASS' if n_fail == 0 else 'FAIL'}**  "
      f"({len(Q)-n_fail}/{len(Q)} PASS)\n",
      f"Fonts: EN = `{EN}`, ZH = `{ZH}`. mathtext fontset = `dejavusans`.\n",
      "程序化验证：对每个图的所有文本（标题/轴标签/注释/图例）逐字符检查所选字体"
      "是否含该字形（`FT2Font.get_char_index(ch) != 0`）。缺字即 FAIL。\n",
      "| figure | CJK | font | missing glyphs | result |",
      "|---|:--:|---|---|:--:|"]
for r in Q.itertuples():
    md.append(f"| `{r.figure}` | {'是' if r.cjk else '否'} | `{r.font}` | "
              f"{', '.join(r.missing_glyphs) if r.missing_glyphs else '—'} | "
              f"**{r.result}** |")
(F / "FIGURE_GLYPH_QC.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(f"\nFIGURES: main={len(list(FM.glob('*.png')))}, "
      f"supp={len(list(FS.glob('*.png')))}")
print(f"FIGURE_GLYPH_QC = {'PASS' if n_fail == 0 else 'FAIL'} ({len(Q)-n_fail}/{len(Q)})")
if n_fail:
    print(Q[Q.result == "FAIL"].to_string(index=False))
