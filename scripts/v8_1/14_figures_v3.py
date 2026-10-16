# -*- coding: utf-8 -*-
"""scripts/v8_1/14_figures_v3.py -- v3 图：统一视觉语言 + 新 Figure 1 + 能力热图。

调色板（§9）：
    between-target C_B   = muted purple  #7A6BA8
    within-target C_W    = green         #4C9A5E
    pooled Pearson       = charcoal      #3A3A3A
    random / permissive  = warm orange   #D98A3D
    family-OOD / strict  = deep blue     #2F5C8A
    negative / failed    = muted red     #B5544E
    external             = desaturated blue-grey #6E7C8C
字体：EN = DejaVu Sans；ZH = Microsoft YaHei（经程序化字形验证）
面板标签 A/B/C/D 加粗、左上角。数学用 mathtext。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
FM = F / "figures_main"; FS = F / "figures_supplement"
FM.mkdir(exist_ok=True); FS.mkdir(exist_ok=True)

PAL = {"CB": "#7A6BA8", "CW": "#4C9A5E", "POOL": "#3A3A3A",
       "RAND": "#D98A3D", "OOD": "#2F5C8A", "NEG": "#B5544E",
       "EXT": "#6E7C8C", "GRID": "#D9D9D9", "INK": "#222222"}
EN = "DejaVu Sans"
ZH = None
for c in ("Microsoft YaHei", "Noto Sans CJK SC", "Source Han Sans SC", "SimHei"):
    try:
        font_manager.findfont(c, fallback_to_default=False); ZH = c; break
    except Exception:
        continue
GLY = []


def use(zh):
    plt.rcParams.update({
        "font.family": ZH if zh else EN,
        "mathtext.fontset": "dejavusans",
        "mathtext.default": "it",
        "axes.unicode_minus": False,
        "figure.dpi": 200,
        "savefig.bbox": "tight",
        "axes.edgecolor": "#555555",
        "axes.linewidth": 0.9,
        "axes.labelcolor": PAL["INK"],
        "text.color": PAL["INK"],
        "xtick.color": PAL["INK"],
        "ytick.color": PAL["INK"],
        "font.size": 9.5,
    })


def panel(ax, letter):
    ax.text(-0.085, 1.10, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left", color=PAL["INK"])


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=PAL["GRID"], lw=0.7, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def gcheck(fig, tag, zh):
    from matplotlib import ft2font
    fam = plt.rcParams["font.family"]
    fam = fam[0] if isinstance(fam, (list, tuple)) else fam
    face = ft2font.FT2Font(font_manager.findfont(fam, fallback_to_default=False))
    miss = set()
    for ax in fig.axes:
        txts = [ax.get_title(), ax.get_xlabel(), ax.get_ylabel()]
        txts += [t.get_text() for t in ax.texts]
        lg = ax.get_legend()
        if lg:
            txts += [t.get_text() for t in lg.get_texts()]
        for s in txts:
            for ch in str(s):
                if ch.isspace() or ord(ch) < 32 or ch in "$\\{}_^&%#":
                    continue
                try:
                    if face.get_char_index(ch) == 0:
                        miss.add(ch)
                except Exception:
                    pass
    GLY.append({"figure": tag, "lang": "ZH" if zh else "EN", "font": fam,
                "missing": "".join(sorted(miss)), "result": "PASS" if not miss else "FAIL"})
    return not miss


def save(fig, stem, zh, folder=FM):
    gcheck(fig, stem, zh)
    suf = "_ZH" if zh else "_EN"
    for ext in ("svg", "pdf", "png"):
        dpi = 600 if ext == "png" else None
        try:
            fig.savefig(str(folder / f"{stem}{suf}.{ext}"), dpi=dpi)
        except Exception:
            pass
    plt.close(fig)


FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
W = pd.read_csv(R / "experiments" / "v8" / "weighting" / "weighting_sensitivity.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
D = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
RS = pd.read_csv(R / "experiments" / "v8" / "remedy" / "remedy_sequence.csv")
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]


# ============================ FIG 1 ============================
def fig1(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(15.6, 4.6))
    gs = fig.add_gridspec(1, 4, wspace=0.30, left=0.035, right=0.985,
                          top=0.845, bottom=0.135)

    # ---- A: benchmark ----
    ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    panel(ax, "A")
    ax.set_title(T("A PPI mutation benchmark", "PPI 突变基准"), fontsize=10,
                 loc="left", pad=20, color=PAL["INK"])
    rng = np.random.RandomState(3)
    for i, (cx, mu) in enumerate([(1.75, 3.0), (5.0, 6.3), (8.25, 4.5)]):
        ax.add_patch(FancyBboxPatch((cx - 1.35, 1.1), 2.7, 6.3,
                                    boxstyle="round,pad=0.07",
                                    fc="#FBFBFC", ec="#9AA4AE", lw=1.2))
        ax.text(cx, 7.75, T(f"Target {i+1}", f"靶点 {i+1}"), ha="center",
                fontsize=8.8, color=PAL["INK"], fontweight="bold")
        for yv in mu + rng.randn(5) * 0.72:
            ax.plot([cx - 1.05, cx + 1.05], [yv, yv], color="#9AA4AE", lw=1.6,
                    alpha=0.75)
        ax.plot([cx - 1.13, cx + 1.13], [mu, mu], color=PAL["INK"], lw=3.0)
        ax.text(cx + 1.22, mu, rf"$\mu_{i+1}$", ha="left", va="center",
                fontsize=8.4, color=PAL["INK"])
    ax.text(5.0, 0.35, T("targets differ in mean ΔΔG and in internal spread",
                         "各靶点的平均 ΔΔG 与靶内离散度均不同"),
            ha="center", fontsize=8.2, style="italic", color="#666")

    # ---- B: two sources ----
    ax = fig.add_subplot(gs[0, 1]); ax.axis("off")
    panel(ax, "B")
    ax.set_title(T("Two sources of signal", "两种信号来源"), fontsize=10,
                 loc="left", pad=20)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(5.0, 8.55, r"$y_{ci} = \mu_c + \varepsilon_{ci}$", ha="center",
            fontsize=15)
    ax.text(5.0, 7.35, r"$\hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}$",
            ha="center", fontsize=15)
    ax.add_patch(FancyBboxPatch((0.35, 3.35), 4.3, 2.9,
                                boxstyle="round,pad=0.10", fc="#F3F0F8",
                                ec=PAL["CB"], lw=1.4))
    ax.text(2.5, 5.62, r"$\mu_c$", ha="center", fontsize=12, color=PAL["CB"])
    ax.text(2.5, 4.35, T("target-level mean\ndifferences between targets",
                         "靶点均值\n靶点之间的差异"), ha="center", fontsize=8.2,
            color=PAL["CB"])
    ax.add_patch(FancyBboxPatch((5.35, 3.35), 4.3, 2.9,
                                boxstyle="round,pad=0.10", fc="#EFF6F0",
                                ec=PAL["CW"], lw=1.4))
    ax.text(7.5, 5.62, r"$\varepsilon_{ci}$", ha="center", fontsize=12,
            color=PAL["CW"])
    ax.text(7.5, 4.35, T("within-target deviation\nordering inside a target",
                         "靶内偏差\n靶点内部的排序"), ha="center", fontsize=8.2,
            color=PAL["CW"])
    ax.add_patch(FancyArrowPatch((2.5, 3.15), (2.5, 2.25), arrowstyle="-|>",
                                 mutation_scale=11, color=PAL["CB"], lw=1.2))
    ax.add_patch(FancyArrowPatch((7.5, 3.15), (7.5, 2.25), arrowstyle="-|>",
                                 mutation_scale=11, color=PAL["CW"], lw=1.2))
    ax.text(2.5, 1.75, T("BETWEEN-TARGET", "跨靶"), ha="center", fontsize=8.6,
            fontweight="bold", color=PAL["CB"])
    ax.text(7.5, 1.75, T("WITHIN-TARGET", "靶内"), ha="center", fontsize=8.6,
            fontweight="bold", color=PAL["CW"])
    ax.text(5.0, 0.85, T("observed effect = target level + mutation deviation",
                         "观测效应 = 靶点水平 + 突变偏差"),
            ha="center", fontsize=8.0, style="italic", color="#666")

    # ---- C: exact decomposition ----
    ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
    panel(ax, "C")
    ax.set_title(T("Exact correlation decomposition",
                   "精确的相关分解"), fontsize=10, loc="left", pad=20)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(5.0, 8.6, T("Pooled Pearson correlation", "合并 Pearson 相关"),
            ha="center", fontsize=9.2, color=PAL["POOL"])
    ax.text(5.0, 7.55, r"$r_{\mathrm{pooled}}$", ha="center", fontsize=15,
            color=PAL["POOL"])
    ax.add_patch(FancyArrowPatch((5.0, 7.05), (5.0, 6.35), arrowstyle="-|>",
                                 mutation_scale=12, color="#888", lw=1.2))
    ax.text(5.0, 5.95, T("exactly splits into", "精确拆分为"), ha="center",
            fontsize=8.2, style="italic", color="#666")
    ax.add_patch(FancyBboxPatch((0.30, 3.05), 4.35, 2.4,
                                boxstyle="round,pad=0.09", fc="#F3F0F8",
                                ec=PAL["CB"], lw=1.5))
    ax.text(2.48, 4.45, r"$C_B$", ha="center", fontsize=13, color=PAL["CB"])
    ax.text(2.48, 3.62, T("Between-target\ncovariance contribution",
                         "跨靶\n协方差贡献"), ha="center", fontsize=7.9,
            color=PAL["CB"])
    ax.add_patch(FancyBboxPatch((5.35, 3.05), 4.35, 2.4,
                                boxstyle="round,pad=0.09", fc="#EFF6F0",
                                ec=PAL["CW"], lw=1.5))
    ax.text(7.52, 4.45, r"$C_W$", ha="center", fontsize=13, color=PAL["CW"])
    ax.text(7.52, 3.62, T("Within-target\ncovariance contribution",
                         "靶内\n协方差贡献"), ha="center", fontsize=7.9,
            color=PAL["CW"])
    ax.text(5.0, 2.42, r"$r_{\mathrm{pooled}} = C_B + C_W$", ha="center",
            fontsize=13.5)
    ax.text(5.0, 1.15, T("covariance-signal contributions on a common scale —\n"
                         "not independent correlations, skill scores or mechanisms",
                         "同一尺度上的协方差信号贡献 ——\n并非独立相关、能力分数或因果机制"),
            ha="center", fontsize=7.6, style="italic", color="#B5544E")

    # ---- D: main result ----
    ax = fig.add_subplot(gs[0, 3]); ax.axis("off")
    panel(ax, "D")
    ax.set_title(T("Controlled split change", "受控划分切换"), fontsize=10,
                 loc="left", pad=20)
    axd = ax.inset_axes([0.05, 0.30, 0.92, 0.60]); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    labs = [T("Random", "随机"), T("Family-OOD", "家族外推")]
    cbs = [H4["random"]["C_B"], H4["family_OOD"]["C_B"]]
    cws = [H4["random"]["C_W"], H4["family_OOD"]["C_W"]]
    rs = [r0, r3]
    xs = np.arange(2)
    axd.bar(xs - 0.185, cbs, width=0.35, color=PAL["CB"], label=r"$C_B$ (between)")
    axd.bar(xs + 0.185, cws, width=0.35, color=PAL["CW"], label=r"$C_W$ (within)")
    axd.axhline(0, color="#999", lw=0.9)
    axd.set_xticks(xs); axd.set_xticklabels(labs, fontsize=8.6)
    axd.set_ylabel(T("covariance contribution", "协方差贡献"), fontsize=8.4)
    axd.tick_params(labelsize=8.0)
    axd.legend(fontsize=7.4, loc="upper right", frameon=False)
    grid(axd)
    axd2 = axd.twinx()
    axd2.plot(xs, rs, "o--", color=PAL["POOL"], lw=1.4, ms=5)
    axd2.set_ylabel(T("pooled Pearson", "合并 Pearson"), fontsize=8.4)
    axd2.tick_params(labelsize=8.0)
    axd2.spines["top"].set_visible(False)
    for i, v in enumerate(rs):
        axd2.annotate(f"{v:.3f}", (i, v), fontsize=8.4, xytext=(0, 8),
                      textcoords="offset points", ha="center", color=PAL["POOL"])
    ax.text(0.5, 0.19, T(f"pooled r  {r0:.3f} → {r3:.3f}", f"合并 r  {r0:.3f} → {r3:.3f}"),
            transform=ax.transAxes, ha="center", fontsize=8.6, color=PAL["POOL"])
    ax.text(0.5, 0.055,
            T(f"{100*DL['fraction_between']:.1f}% of the absolute decrease is carried "
              f"by the between-target component",
              f"绝对下降的 {100*DL['fraction_between']:.1f}% 由跨靶分量承载"),
            transform=ax.transAxes, ha="center", fontsize=8.8,
            color=PAL["CB"], fontweight="bold")
    save(fig, "Fig1_overview", zh)


# ============================ FIG 2 ============================
def fig2(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    order = ["random", "complex", "seq40", "family_OOD"]
    lab_en = ["Random", "Complex", "seq40", "Family-OOD"]
    lab_zh = ["随机", "复合物", "seq40", "家族外推"]
    lab = lab_zh if zh else lab_en
    rp = [H4[k]["overall_pearson"] for k in order]
    cb = [H4[k]["C_B"] for k in order]
    cw = [H4[k]["C_W"] for k in order]
    xs = np.arange(4)
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.3))
    fig.subplots_adjust(wspace=0.30, left=0.055, right=0.975, top=0.80, bottom=0.17)
    # A
    ax = axes[0]; panel(ax, "A")
    ax.plot(xs, rp, "o-", color=PAL["POOL"], lw=2.0, ms=6)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=8.4)
    ax.set_ylabel(T("pooled Pearson", "合并 Pearson"), fontsize=9.2)
    ax.set_title(T("Pooled correlation by split regime",
                   "各划分下的合并相关"), fontsize=9.8, loc="left", pad=8)
    for i, v in enumerate(rp):
        ax.annotate(f"{v:.3f}", (i, v), fontsize=8.0, xytext=(0, 8),
                    textcoords="offset points", ha="center")
    grid(ax)
    # B
    ax = axes[1]; panel(ax, "B")
    ax.plot(xs, cb, "s-", color=PAL["CB"], lw=2.0, ms=6, label=r"$C_B$ (between)")
    ax.plot(xs, cw, "d-", color=PAL["CW"], lw=2.0, ms=6, label=r"$C_W$ (within)")
    ax.axhline(0, color="#999", lw=0.9)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=8.4)
    ax.set_ylabel(T("covariance contribution", "协方差贡献"), fontsize=9.2)
    ax.set_title(T("Exact decomposition  $r_{\\rm pooled}=C_B+C_W$",
                   "精确分解  $r_{\\rm pooled}=C_B+C_W$"), fontsize=9.8, loc="left", pad=8)
    ax.legend(fontsize=8.0, frameon=False)
    grid(ax)
    ax.annotate(T("$C_W$ stays positive", "$C_W$ 保持为正"), xy=(3, cw[3]),
                xytext=(1.1, 0.245), fontsize=8.0, color=PAL["CW"],
                arrowprops=dict(arrowstyle="->", color=PAL["CW"], lw=1.0))
    # C
    ax = axes[2]; panel(ax, "C")
    d_cb, d_cw = abs(DL["delta_C_B"]), abs(DL["delta_C_W"])
    ax.barh([1], [d_cb], color=PAL["CB"], height=0.48)
    ax.barh([0], [d_cw], color=PAL["CW"], height=0.48)
    ax.set_yticks([1, 0])
    ax.set_yticklabels([T("Between-target\n$\\Delta C_B$", "跨靶\n$\\Delta C_B$"),
                        T("Within-target\n$\\Delta C_W$", "靶内\n$\\Delta C_W$")],
                       fontsize=8.6)
    ax.set_xlabel(T("absolute change, Random → Family-OOD",
                    "绝对变化，随机 → 家族外推"), fontsize=8.8)
    ax.text(d_cb * 0.5, 1.32, f"{100*DL['fraction_between']:.1f}%", ha="center",
            fontsize=13, color=PAL["CB"], fontweight="bold")
    ax.text(d_cw * 0.5, -0.32, f"{100*DL['fraction_within']:.1f}%", ha="center",
            fontsize=10, color=PAL["CW"], fontweight="bold")
    ax.set_xlim(0, d_cb * 1.18); ax.set_ylim(-0.65, 1.75)
    ax.set_title(T("Attribution of the decline", "下降的归因"),
                 fontsize=9.8, loc="left", pad=8)
    grid(ax, axis="x")
    fig.suptitle(T("Figure 2 | Family-OOD evaluation disproportionately removes the "
                   "between-target component of pooled correlation",
                   "图 2 | 家族外推评价不成比例地移除合并相关中的跨靶分量"),
                 fontsize=11.4, y=0.955, fontweight="bold", x=0.055, ha="left")
    fig.text(0.055, 0.055,
             T("Identical 7,523 mutation rows and identical model configurations; only the "
               "split regime changes. The complex and seq40 partitions are identical in "
               "this dataset and are not independent comparisons. Panel C uses absolute "
               "changes, $|\\Delta C_B|/(|\\Delta C_B|+|\\Delta C_W|)$.",
               "同一批 7,523 行突变、同一批模型配置，仅改变划分方式。complex 与 seq40 "
               "在本数据上完全相同，不构成独立比较。面板 C 采用绝对变化口径，"
               "$|\\Delta C_B|/(|\\Delta C_B|+|\\Delta C_W|)$。"),
             fontsize=8.0, color="#555")
    save(fig, "Fig2_split_decomposition", zh)


# ============================ FIG 3 ============================
def fig3(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    SCH = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.1))
    fig.subplots_adjust(wspace=0.30, left=0.055, right=0.975, top=0.79, bottom=0.20)
    for ax, key, col, nm, pnl in ((axes[0], "C_B", PAL["CB"], r"$C_B$ (between-target)", "A"),
                                  (axes[1], "C_W", PAL["CW"], r"$C_W$ (within-target)", "B")):
        panel(ax, pnl)
        cols = [f"{key}_{s}" for s in SCH]
        A_ = W[cols].to_numpy(float)
        rng_ = np.nanmax(A_, axis=1) - np.nanmin(A_, axis=1)
        ax.hist(rng_, bins=26, color=col, alpha=0.80, edgecolor="white", lw=0.35)
        ax.axvline(np.nanmedian(rng_), color=PAL["POOL"], ls="--", lw=1.5)
        ax.text(np.nanmedian(rng_) * 1.06, ax.get_ylim()[1] * 0.90,
                T(f"median {np.nanmedian(rng_):.4f}", f"中位 {np.nanmedian(rng_):.4f}"),
                fontsize=8.0, color=PAL["POOL"])
        ax.set_xlabel(T("range across the five weighting schemes",
                        "五种加权方案下的极差"), fontsize=8.8)
        ax.set_ylabel(T("model configurations", "模型配置数"), fontsize=8.8)
        ax.set_title(nm, fontsize=9.8, loc="left", pad=8)
        grid(ax)
    ax = axes[2]; panel(ax, "C")
    thr = [0, 0.001, 0.005]
    cb_fl = [int(W2[f"signflip_C_B_gt{t}"].sum()) for t in thr]
    cw_fl = [int(W2[f"signflip_C_W_gt{t}"].sum()) for t in thr]
    xx = np.arange(3); w = 0.36
    ax.bar(xx - w / 2, cb_fl, width=w, color=PAL["CB"], label=r"$C_B$")
    ax.bar(xx + w / 2, cw_fl, width=w, color=PAL["CW"], label=r"$C_W$")
    ax.set_xticks(xx)
    ax.set_xticklabels(["0", "0.001", "0.005"], fontsize=8.6)
    ax.set_xlabel(T("sign-flip threshold on $|C|$", "$|C|$ 的符号翻转阈值"),
                  fontsize=8.8)
    ax.set_ylabel(T("configurations reversing sign", "出现符号翻转的配置数"),
                  fontsize=8.8)
    ax.set_title(T("Directional stability", "方向稳定性"), fontsize=9.8,
                 loc="left", pad=8)
    ax.legend(fontsize=8.0, frameon=False)
    grid(ax)
    for i, v in enumerate(cb_fl):
        ax.annotate(str(v), (i - w / 2, v), fontsize=8.0, xytext=(0, 4),
                    textcoords="offset points", ha="center", color=PAL["CB"])
    fig.suptitle(T("Figure 3 | Target weighting changes the interpretation of "
                   "covariance contributions",
                   "图 3 | 目标加权改变协方差贡献的解读"),
                 fontsize=11.4, y=0.955, fontweight="bold", x=0.055, ha="left")
    fig.text(0.055, 0.055,
             T("Predictions are unchanged; only aggregation weights differ "
               "(per-mutation, per-target, $\\sqrt{n}$, $\\log(1+n)$, $\\min(n,30)$). "
               "$C_B$ reverses direction only near zero, whereas $C_W$ keeps its sign at "
               "every threshold but its magnitude is not invariant.",
               "预测值不变，仅改变聚合权重（每个突变、每个靶点、$\\sqrt{n}$、"
               "$\\log(1+n)$、$\\min(n,30)$）。$C_B$ 仅在近零处翻转方向，"
               "而 $C_W$ 在任何阈值下均保持符号，但其幅度并非不变。"),
             fontsize=8.0, color="#555")
    save(fig, "Fig3_weighting", zh)


# ============================ FIG 4 ============================
def fig4(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.86, bottom=0.20)
    K = 10
    m = GEO_DF.n <= 160
    ax.plot(GEO_DF.n[m], GEO_DF.D[m], lw=2.4, color=PAL["OOD"])
    ax.axvline(K, color=PAL["RAND"], ls="--", lw=1.5)
    ax.text(K + 2.0, GEO_DF.D[m].max() * 0.90, r"$n \leq K \Rightarrow D = 0$",
            fontsize=9.0, color=PAL["RAND"])
    ns_ = GEO["n_star_analytic"]
    ax.axvline(ns_, color=PAL["CW"], ls=":", lw=1.5)
    ax.annotate(T(f"peak  $n^*\\approx{ns_:.1f}$", f"峰值  $n^*\\approx{ns_:.1f}$"),
                xy=(ns_, GEO_DF.D[m].max()),
                xytext=(ns_ + 7, GEO_DF.D[m].max() * 0.68), fontsize=9.0,
                color=PAL["CW"],
                arrowprops=dict(arrowstyle="->", color=PAL["CW"], lw=1.0))
    ax.plot([20, 50], [GEO["D_at_n20"], GEO["D_at_n50"]], "o", ms=6,
            color=PAL["POOL"], zorder=4)
    ax.annotate(f"D(20) = {GEO['D_at_n20']:.3f}", (20, GEO["D_at_n20"]),
                xytext=(24, GEO["D_at_n20"] + 0.018), fontsize=8.6, color=PAL["POOL"])
    ax.annotate(f"D(50) = {GEO['D_at_n50']:.3f}", (50, GEO["D_at_n50"]),
                xytext=(60, GEO["D_at_n50"] + 0.010), fontsize=8.6, color=PAL["POOL"])
    ax.set_xlabel(T("candidate-set size $n$  ($K = 10$)", "候选集大小 $n$（$K = 10$）"),
                  fontsize=9.4)
    ax.set_ylabel(r"discriminable interval  $D(n,K)$", fontsize=9.4)
    ax.set_title(T("Figure 4 | Fixed-K shortlist evaluation is constrained by "
                   "candidate-set geometry",
                   "图 4 | 固定 K 短名单评价受候选集几何约束"),
                 fontsize=11.2, loc="left", pad=10, fontweight="bold")
    grid(ax)
    fig.text(0.10, 0.055,
             T(f"$D(n,K)=0$ for $n\\leq K$ and $D(n,K)=K(n-K)/[n(n+1)]$ otherwise. "
               f"In this benchmark {GEO['n_le_K_count']} of {GEO['n_ge5_count']} targets "
               f"with $n\\geq5$ ({100*GEO['n_le_K_fraction']:.1f}%) satisfy $n\\leq K$, "
               f"where every model selects all candidates and the target carries no "
               f"information for fixed-K shortlist metrics.",
               f"$n\\leq K$ 时 $D(n,K)=0$，否则 $D(n,K)=K(n-K)/[n(n+1)]$。"
               f"本基准中 {GEO['n_ge5_count']} 个 $n\\geq5$ 靶点有 "
               f"{GEO['n_le_K_count']} 个（{100*GEO['n_le_K_fraction']:.1f}%）"
               f"满足 $n\\leq K$，此时所有模型都会选中全部候选，"
               f"该靶点对固定 K 短名单指标不携带信息。"),
             fontsize=8.0, color="#555")
    save(fig, "Fig4_candidate_geometry", zh)


# ============================ FIG 5 ============================
NICE = {"F0_Base": ("Base model", "基线模型"),
        "F1_Base_Mut": ("Mutation properties", "突变性质"),
        "F2_Base_Static": ("Static local geometry", "静态局部几何"),
        "F3_Base_EnvPartner": ("Partner static context", "partner 静态环境"),
        "F4_Base_PartnerPert": ("Mutation × partner perturbation", "突变 × partner 扰动"),
        "F5_Base_IntraPert": ("Intra-chain perturbation", "链内扰动"),
        "F6_Partner_Intra": ("Partner + intra-chain", "partner + 链内"),
        "F7_F6_Burial": ("Partner + intra + burial", "partner + 链内 + 埋藏"),
        "F8_Full": ("Full local-physics model", "全量局部物理模型"),
        "F9_PhysicsOnly": ("Physics features only", "仅物理特征")}


def fig5(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig, ax = plt.subplots(figsize=(9.8, 6.0))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.86, bottom=0.16)
    b = RS[RS.label.str.startswith("V7.3 F1")]
    r5 = RS[RS.step == "R5"].copy()
    if len(b):
        b = b.iloc[0]
        r5["d_ov"] = r5.overall_r - b.overall_r
        r5["d_mc"] = r5.macro_spearman - b.macro_spearman
        ax.axhline(0, color="#999", lw=1.0)
        ax.axvline(0, color="#999", lw=1.0)
        ax.axhspan(0, 0.20, xmin=0.5, xmax=1.0, color="#EFF6F0", alpha=0.55, zorder=0)
        for _, r in r5.iterrows():
            frag = r.label.replace("V7.3 ", "").replace("_s42", "")
            nm = NICE.get(frag, (frag, frag))
            real = r.d_mc > 0.02
            col = PAL["CW"] if real else PAL["NEG"]
            ax.scatter(r.d_ov, r.d_mc, s=150, c=col, zorder=4,
                       edgecolors="white", linewidths=1.3)
            dx, dy = (8, 7) if real else (8, -14)
            ax.annotate(nm[1] if zh else nm[0], (r.d_ov, r.d_mc), fontsize=8.6,
                        xytext=(dx, dy), textcoords="offset points", color=col)
        ax.text(0.055, 0.055, T("pooled-score-only change", "仅合并分数改善"),
                transform=ax.transAxes, fontsize=8.2, color=PAL["NEG"],
                style="italic")
        ax.text(0.60, 0.93, T("both axes improve", "两轴同时改善"),
                transform=ax.transAxes, fontsize=8.2, color=PAL["CW"],
                style="italic")
    ax.set_xlabel(T(r"$\Delta$ pooled Pearson vs mutation-properties baseline",
                    r"相对「突变性质」基线的 $\Delta$ 合并 Pearson"), fontsize=9.4)
    ax.set_ylabel(T(r"$\Delta$ macro within-target Spearman",
                    r"$\Delta$ 靶内 Macro Spearman"), fontsize=9.4)
    ax.set_title(T("Figure 5 | Capability decomposition separates pooled-score gains "
                   "from within-target gains",
                   "图 5 | 能力分解区分合并分数增益与靶内增益"),
                 fontsize=11.2, loc="left", pad=10, fontweight="bold")
    grid(ax, axis="both")
    fig.text(0.10, 0.045,
             T("Retrospective matched ablation on an identical 3,351-row universe. "
               "Green: improvement on both axes. Red: pooled-score-only or negative "
               "change. The analysis is retrospective and does not imply prospective "
               "prediction of which modification will succeed.",
               "在同一 3,351 行宇宙上的回溯性匹配消融。绿色：两轴同时改善；"
               "红色：仅合并分数改善或变化为负。该分析为回溯性，"
               "不代表可前瞻预测哪种改动会成功。"),
             fontsize=8.0, color="#555")
    save(fig, "Fig5_remedy", zh)


# ============================ FIG 6 ============================
def fig6(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig, axes = plt.subplots(1, 2, figsize=(15.6, 4.8))
    fig.subplots_adjust(wspace=0.26, left=0.055, right=0.975, top=0.80, bottom=0.20)
    ax = axes[0]; panel(ax, "A")
    for r in EP.itertuples():
        ax.scatter(r.overall_pearson, r.macro_spearman, s=155, color=PAL["EXT"],
                   zorder=3, edgecolors="white", linewidths=1.3)
        ax.annotate(r.model, (r.overall_pearson, r.macro_spearman), fontsize=8.8,
                    xytext=(9, 6), textcoords="offset points", color=PAL["INK"])
    ax.axhline(0, color="#999", lw=1.0); ax.axvline(0, color="#999", lw=1.0)
    ax.set_xlabel(T("pooled Pearson", "合并 Pearson"), fontsize=9.2)
    ax.set_ylabel(T("macro within-target Spearman", "靶内 Macro Spearman"),
                  fontsize=9.2)
    ax.set_title(T("Within- vs between-target performance", "靶内与跨靶性能"),
                 fontsize=9.8, loc="left", pad=8)
    grid(ax, axis="both")
    ax = axes[1]; panel(ax, "B")
    names = list(EP.model)
    cols = ["overall_pearson", "macro_spearman", "C_B_row", "C_W"]
    heads = [T("Pooled\nPearson", "合并\nPearson"), T("Macro\nSpearman", "Macro\nSpearman"),
             r"$C_B$", r"$C_W$"]
    M = EP[cols].to_numpy(float)
    Z = (M - M.mean(axis=0)) / np.where(M.std(axis=0) > 0, M.std(axis=0), 1)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(heads, fontsize=8.6)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8.8)
    for i in range(len(names)):
        for j in range(4):
            ax.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center", fontsize=8.0,
                    color="white" if abs(Z[i, j]) > 1.25 else PAL["INK"])
    ax.set_title(T("Method-level capability profile (column-standardised)",
                   "方法级能力画像（按列标准化）"), fontsize=9.8, loc="left", pad=8)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.80, label=T("within-column z-score",
                                                 "列内 z 分数"))
    fig.suptitle(T("Figure 6 | External methods define the boundary of generalisation",
                   "图 6 | 外部方法界定泛化边界"),
                 fontsize=11.4, y=0.955, fontweight="bold", x=0.055, ha="left")
    fig.text(0.055, 0.055,
             T("Five independently developed methods on a matched, unambiguous mapping "
               "(1,781 mutations). Rank agreement is positive "
               "($\\rho = +0.90$, exact permutation $p = 0.083$); the panel is "
               "underpowered with five methods, so this positive alignment does not "
               "support a universal decoupling of pooled and within-target performance. "
               "Colours are standardised within each metric for visual comparison only; "
               "statistical analyses use raw metrics.",
               "五种独立开发的方法，在匹配且无歧义的映射（1,781 个突变）上评价。"
               "秩一致为正（$\\rho = +0.90$，精确置换 $p = 0.083$）；"
               "仅 5 个方法时功效不足，该正向对齐不支持合并性能与靶内性能的普遍解耦。"
               "颜色仅按列标准化以便视觉比较；统计分析使用原始指标。"),
             fontsize=8.0, color="#555")
    save(fig, "Fig6_external", zh)


# ==================== CAPABILITY HEATMAP ====================
def heatmap(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    # 主分层由 V7.3 消融变体构成，故其代表性键未被标记（§4.2：不对消融做折叠）。
    # 这里取每个消融变体的 seed-42 配置作为**同一方法族的受控变体**，
    # 并在图注中如实说明这一点。
    g = D[(D.rowset_hash == PU["rowset_hash"]) &
          D.key.str.endswith("_s42")].dropna(
        subset=["overall_pearson", "macro_spearman", "C_B_row", "C_W"])
    g = g.sort_values("macro_spearman", ascending=False)
    NICE_LBL = {"F0_Base": "Base model", "F1_Base_Mut": "Mutation properties",
                "F2_Base_Static": "Static local geometry",
                "F3_Base_EnvPartner": "Partner static context",
                "F4_Base_PartnerPert": "Mutation × partner perturbation",
                "F5_Base_IntraPert": "Intra-chain perturbation",
                "F6_Partner_Intra": "Partner + intra-chain",
                "F7_F6_Burial": "Partner + intra + burial",
                "F8_Full": "Full local-physics model",
                "F9_PhysicsOnly": "Physics features only"}
    ylabels = [NICE_LBL.get(m, m) for m in g.method_family]
    fig, ax = plt.subplots(figsize=(9.6, max(3.4, 0.34 * len(g) + 1.5)))
    fig.subplots_adjust(left=0.28, right=0.90, top=0.84, bottom=0.16)
    cols = ["overall_pearson", "macro_spearman", "C_B_row", "C_W"]
    heads = [T("Pooled\nPearson", "合并\nPearson"),
             T("Macro\nSpearman", "Macro\nSpearman"),
             T("$C_B$\n(between)", "$C_B$\n（跨靶）"),
             T("$C_W$\n(within)", "$C_W$\n（靶内）")]
    M = g[cols].to_numpy(float)
    Z = (M - M.mean(axis=0)) / np.where(M.std(axis=0) > 0, M.std(axis=0), 1)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(heads, fontsize=9.0)
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels(ylabels, fontsize=8.2)
    for i in range(len(g)):
        for j in range(4):
            ax.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center", fontsize=7.6,
                    color="white" if abs(Z[i, j]) > 1.25 else PAL["INK"])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(T("Capability profiles of matched, provenance-clean method families",
                   "匹配且来源清洁的方法族能力画像"),
                 fontsize=10.6, loc="left", pad=12, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.75, label=T("within-column z-score",
                                                 "列内 z 分数"))
    fig.text(0.28, 0.045,
             T(f"Primary matched evaluation set: {PU['n_rows']:,} mutation rows, "
               f"{PU['n_ge5']} targets with $n\\geq5$. Rows are the {len(g)} controlled "
               f"ablation variants of a single local-physics method family (seed 42), "
               f"not independent methods; values are standardised within each metric for "
               f"visual comparison only and statistical analyses use raw metrics. No "
               f"cross-metric correlation is implied by this display.",
               f"主匹配评价集合：{PU['n_rows']:,} 行突变、{PU['n_ge5']} 个 $n\\geq5$ 靶点。"
               f"行为同一局部物理方法族的 {len(g)} 个受控消融变体（seed 42），"
               f"并非独立方法；数值仅按列标准化以便视觉比较，统计分析使用原始指标。"
               f"该展示不隐含任何跨指标相关。"),
             fontsize=7.8, color="#555")
    save(fig, "Capability_heatmap_provenance_clean", zh)


# ==================== SUPPLEMENT ====================
def supp(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig, ax = plt.subplots(figsize=(10.6, 4.8))
    fig.subplots_adjust(left=0.24, right=0.97, top=0.86, bottom=0.16)
    u = UR.sort_values("n_rows"); y = np.arange(len(u))
    ax.barh(y, u.n_rows, color=PAL["EXT"], alpha=0.9,
            label=T("rows", "行数"))
    ax.barh(y, u.n_ge5, color=PAL["RAND"], alpha=0.95, height=0.55,
            label=T("targets with n >= 5", "n≥5 的靶点"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ({h[:8]})" for a, h in zip(u.artifact, u.rowset_hash)],
                       fontsize=7.6)
    ax.set_xlabel(T("mutation rows", "突变行数"), fontsize=9.0)
    ax.legend(fontsize=8.2, frameon=False)
    ax.set_title(T("Figure S1 | Matched-evaluation strata and row-set fingerprints",
                   "图 S1 | 匹配评价分层与行集指纹"),
                 fontsize=10.4, loc="left", pad=8, fontweight="bold")
    grid(ax, axis="x")
    save(fig, "FigS1_matched_universes", zh, FS)

    fig, ax = plt.subplots(figsize=(10.6, 4.6))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.84, bottom=0.22)
    s = ST.sort_values("rho_overall_vs_macro")
    yy = np.arange(len(s))
    colors = [PAL["NEG"] if v < 0 else PAL["OOD"] for v in s.rho_overall_vs_macro]
    ax.barh(yy, s.rho_overall_vs_macro, color=colors, alpha=0.9)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{h[:8]} / L{lv}" for h, lv in
                        zip(s.rowset_hash, s.level)], fontsize=7.4)
    ax.axvline(0, color="#999", lw=1.0)
    ax.set_xlabel(T(r"$\rho$(pooled Pearson, macro within-target Spearman)",
                    r"$\rho$（合并 Pearson，靶内 Macro Spearman）"), fontsize=9.0)
    ax.set_title(T("Figure S4 | Pool sensitivity of cross-metric rank agreement",
                   "图 S4 | 跨指标秩一致的池敏感性"),
                 fontsize=10.4, loc="left", pad=8, fontweight="bold")
    grid(ax, axis="x")
    save(fig, "FigS4_capability_pool_sensitivity", zh, FS)


for zh in (False, True):
    print(f"--- {'ZH' if zh else 'EN'} ---")
    fig1(zh); fig2(zh); fig3(zh); fig4(zh); fig5(zh); fig6(zh); heatmap(zh); supp(zh)

Q = pd.DataFrame(GLY)
Q.to_csv(A / "figure_glyph_qc_v3.csv", index=False, encoding="utf-8")
nf = int((Q.result == "FAIL").sum())
print(f"\nFIGURES v3: main={len(list(FM.glob('*_EN.png')))} EN / "
      f"{len(list(FM.glob('*_ZH.png')))} ZH ; supp={len(list(FS.glob('*.png')))}")
print(f"FIGURE_GLYPH_QC = {'PASS' if nf == 0 else 'FAIL'} ({len(Q)-nf}/{len(Q)})")
if nf:
    print(Q[Q.result == "FAIL"].to_string(index=False))
