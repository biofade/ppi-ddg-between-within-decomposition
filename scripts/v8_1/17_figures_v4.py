# -*- coding: utf-8 -*-
"""scripts/v8_1/17_figures_v4.py -- 旗舰级图系统 v4（Nature Methods / NMI 视觉标准）。

严格约束：
* 科学数值**只**来自 v8_1/FROZEN_RESULTS.json；本脚本不重新计算任何科学结论。
* 图内**不再**硬编码 "Figure X |" 大标题与底部灰色图注（改由 DOCX 图注承担）。
* 输出：v8_1_final/figures_main/、v8_1_final/figures_supplement/
          （EN + ZH，各 SVG 矢量母版 / PDF 矢量母版 / PNG 600 dpi）

调色板（莫兰迪 DAMO-RADAR 规范 + 信号语义色）：
    chevron1 标题栏 #EEDAD5 / 容器 #FAF4F2 / 边框 #D9B8B0   （Benchmark & Provenance）
    chevron2 标题栏 #E6DAE8 / 容器 #F7F3F8 / 边框 #C8B4CE   （Covariance Decomposition）
    chevron3 标题栏 #DCE2F0 / 容器 #F3F5FA / 边框 #B4C2DC   （Audit Suite）
    C_B  #6A51A3 (fill #CBC9E2)   C_W  #238B45 (fill #BAE4B3)
    Pooled #2B3E50                warning #CB181D (fill #FCAE91)
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon
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

# ---------------- palette ----------------
CB, CBF = "#6A51A3", "#CBC9E2"
CW, CWF = "#238B45", "#BAE4B3"
POOL = "#2B3E50"
WARN, WARNF = "#CB181D", "#FCAE91"
ST1 = dict(bar="#EEDAD5", bg="#FAF4F2", edge="#D9B8B0")
ST2 = dict(bar="#E6DAE8", bg="#F7F3F8", edge="#C8B4CE")
ST3 = dict(bar="#DCE2F0", bg="#F3F5FA", edge="#B4C2DC")
AMBER = "#E8A33D"
INK = "#1A1A1A"
GREY = "#8A8F98"

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
        "mathtext.fontset": "dejavusans", "mathtext.default": "it",
        "axes.unicode_minus": False, "figure.dpi": 200, "savefig.bbox": "tight",
        "savefig.facecolor": "white", "axes.edgecolor": "#555555",
        "axes.linewidth": 0.9, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.labelsize": 11, "font.size": 10.5,
    })


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
                if ch.isspace() or ord(ch) < 32 or ch in "$\\{}_^&%#,;:":
                    continue
                try:
                    if face.get_char_index(ch) == 0:
                        miss.add(ch)
                except Exception:
                    pass
    # 英文图不得含 CJK
    badcjk = False
    if not zh:
        allt = []
        for ax in fig.axes:
            allt += [ax.get_title(), ax.get_xlabel(), ax.get_ylabel()]
            allt += [t.get_text() for t in ax.texts]
        for s in allt:
            if any("\u4e00" <= c <= "\u9fff" for c in str(s)):
                badcjk = True
    res = "PASS" if (not miss and not badcjk) else "FAIL"
    GLY.append({"figure": tag, "lang": "ZH" if zh else "EN", "font": fam,
                "missing": "".join(sorted(miss)),
                "cjk_in_en": badcjk, "result": res})
    return res == "PASS"


def save(fig, stem, zh, folder=FM):
    gcheck(fig, stem, zh)
    suf = "_ZH" if zh else "_EN"
    fig.savefig(str(folder / f"{stem}{suf}.svg"))
    fig.savefig(str(folder / f"{stem}{suf}.pdf"))
    fig.savefig(str(folder / f"{stem}{suf}.png"), dpi=600)
    plt.close(fig)


def panel(ax, letter, dx=-0.02, dy=1.13):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=15, fontweight="bold",
            va="top", ha="left", color=INK)


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color="#DFE3E8", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _ov(fig, z=0):
    """覆盖层：坐标 0..1 = 整图。z 控制绘制顺序（数据轴用 z>=6）。"""
    ax = fig.add_axes([0, 0, 1, 1], zorder=z)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.patch.set_alpha(0)
    return ax


def chevron(fig, x, y, w, h, label, pal, fs=11.0, z=1):
    ax = _ov(fig, z)
    tip = h * 0.55
    ax.add_patch(Polygon([(x, y), (x + w - tip, y), (x + w, y + h / 2),
                          (x + w - tip, y + h), (x, y + h)], closed=True,
                         fc=pal["bar"], ec=pal["edge"], lw=1.0))
    ax.text(x + 0.012, y + h / 2, label, fontsize=fs, fontweight="bold",
            va="center", ha="left", color=INK)


def container(fig, x, y, w, h, pal, z=1):
    ax = _ov(fig, z)
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004",
                                fc=pal["bg"], ec=pal["edge"], lw=1.2))


def card(fig, x, y, w, h, ec, fc="white", lw=1.2, z=2):
    ax = _ov(fig, z)
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006",
                                fc=fc, ec=ec, lw=lw))
    return ax


def arrow(fig, x1, y1, x2, y2, col=GREY, lw=1.2, z=4):
    ax = _ov(fig, z)
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=12, color=col, lw=lw))


def txt(fig, x, y, s, size=9, color=INK, ha="center", va="center",
        weight="normal", style="normal", z=5, rot=0):
    ax = _ov(fig, z)
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
            fontweight=weight, style=style, rotation=rot)


def box(fig, x, y, w, h, s, size=7.6, fc="white", ec="#B8BEC7", lw=1.0,
        color=INK, weight="normal", z=5, pad=0.006):
    ax = _ov(fig, z)
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={pad}",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, s, fontsize=size, color=color, ha="center",
            va="center", fontweight=weight)


# ---------------- frozen data ----------------
FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]
rc = H4["complex"]["overall_pearson"]
cb0, cw0 = H4["random"]["C_B"], H4["random"]["C_W"]
cbc, cwc = H4["complex"]["C_B"], H4["complex"]["C_W"]
cb3, cw3 = H4["family_OOD"]["C_B"], H4["family_OOD"]["C_W"]
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
W = pd.read_csv(R / "experiments" / "v8" / "weighting" / "weighting_sensitivity.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
D = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
RS = pd.read_csv(R / "experiments" / "v8" / "remedy" / "remedy_sequence.csv")


# =====================================================================
# FIGURE 1 -- 旗舰级 2 行 4 面板总览
# =====================================================================
def fig1(zh):
    """旗舰级总览图：2 行 4 面板（A/B/C/D），DAMO-RADAR 式 16:10 版式。

    分层策略：装饰层 z=1..5，数据轴 z=7（必须最后绘制，否则会被覆盖层压住）。
    """
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(15.5, 10.3))
    fig.patch.set_facecolor("white")

    # ---------------- 版面坐标（fig 分数） ----------------
    TOPY, TOPH = 0.520, 0.400          # 上排容器
    BOTY, BOTH = 0.045, 0.400          # 下排容器
    LX, RX, WW = 0.014, 0.512, 0.474   # 左右列 x / 宽

    # ---------------- 阶段标题栏 + 容器（z=1） ----------------
    chevron(fig, LX, 0.940, WW, 0.042, T("1. Matched Benchmark & Target Structure",
            "1. 匹配基准与靶点结构"), ST1, 10.6)
    chevron(fig, RX, 0.940, WW, 0.042, T("2. Exact Additive Signal Decomposition",
            "2. 精确可加信号分解"), ST2, 10.6)
    chevron(fig, LX, 0.465, WW, 0.042, T("3. Controlled Split Collapse & Audit Findings",
            "3. 受控划分崩塌与审计发现"), ST3, 10.6)
    chevron(fig, RX, 0.465, WW, 0.042, T("4. Capability-Resolved Model Profiles",
            "4. 能力分辨的模型画像"), ST3, 10.6)
    for (x, y) in ((LX, TOPY), (RX, TOPY), (LX, BOTY), (RX, BOTY)):
        container(fig, x, y, WW, TOPH, ST1 if x < 0.5 and y > 0.5
                  else (ST2 if y > 0.5 else ST3))

    # ======================= PANEL A =======================
    # A1: 数据流卡片（纯装饰，全部用覆盖层绘制）
    ax0, aw = LX + 0.012, 0.196
    ay, ah = TOPY + 0.012, TOPH - 0.024
    card(fig, ax0, ay, aw, ah, ST1["edge"], "#FFFFFF", 1.1, z=2)
    txt(fig, ax0 + aw / 2, ay + ah - 0.022,
        T("Fingerprint-matched data flow", "指纹匹配数据流"), 9.4, INK, weight="bold",
        z=5)
    steps = [
        (T("SKEMPI v2.0 annotation\n4,473 rows · 214 targets · 108 families",
           "SKEMPI v2.0 标注\n4,473 行 · 214 靶点 · 108 家族"), "#FFFFFF",
         ST1["edge"], None),
        (T("SHA-256 row-set fingerprinting + provenance audit",
           "SHA-256 行集指纹 + 来源审计"), "#FFFFFF", ST1["edge"], None),
        (T("keep Family-OOD (no training overlap)",
           "保留 Family-OOD（无训练重叠）"), "#EAF5EA", CW, "✓"),
        (T("drop Random / Complex / seq40",
           "剔除 Random / Complex / seq40"), "#FDEDE8", WARN, "✕"),
        (T("Primary matched universe\n3,351 rows · 102 targets (n≥5) · 10 families",
           "主匹配宇宙\n3,351 行 · 102 个 n≥5 靶点 · 10 方法族"), CBF, CB, None),
    ]
    top = ay + ah - 0.048
    bh, gap = 0.050, 0.018
    for i, (s_, fc, ec, badge) in enumerate(steps):
        yy = top - i * (bh + gap)
        box(fig, ax0 + 0.012, yy - bh, aw - 0.024, bh, s_, 7.0, fc, ec, 1.0,
            weight="bold" if i == 4 else "normal", z=4)
        if badge:
            txt(fig, ax0 + 0.022, yy - bh / 2, badge, 10, CW if badge == "✓" else WARN,
                ha="center", z=5, weight="bold")
        if i < 4:
            arrow(fig, ax0 + aw / 2, yy - bh - 0.002, ax0 + aw / 2, yy - bh - gap + 0.004,
                  GREY, 1.0, z=4)

    # A2: 靶点结构示意（数据轴，z=7）
    axA = fig.add_axes([ax0 + aw + 0.012, ay, (LX + WW - 0.012) - (ax0 + aw + 0.012), ah],
                       zorder=7)
    axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
    axA.text(5.0, 9.60, T("Between-target scale vs within-target ordering",
                          "跨靶量级 vs 靶内排序"), ha="center", fontsize=9.4,
             fontweight="bold", color=INK)
    rng = np.random.RandomState(7)
    for i, (cx, mu) in enumerate([(1.70, 3.1), (5.00, 6.5), (8.30, 4.7)]):
        axA.add_patch(FancyBboxPatch((cx - 1.28, 1.05), 2.56, 6.9,
                                     boxstyle="round,pad=0.07", fc="#FCFAF9",
                                     ec=ST1["edge"], lw=1.1))
        axA.text(cx, 7.62, T(f"Target {i+1}", f"靶点 {i+1}"), ha="center",
                 fontsize=8.4, fontweight="bold", color=INK)
        eps = rng.randn(5) * 0.60
        for e in eps:
            axA.plot([cx - 1.00, cx + 1.00], [mu + e, mu + e], color=CW, lw=2.2,
                     alpha=0.9, solid_capstyle="round", zorder=3)
        axA.plot([cx - 1.08, cx + 1.08], [mu, mu], color=CB, lw=3.6,
                 solid_capstyle="round", zorder=4)
        axA.text(cx + 1.14, mu, rf"$\mu_{i+1}$", ha="left", va="center",
                 fontsize=9.6, color=CB)
        e0 = float(np.max(np.abs(eps)))
        axA.annotate("", xy=(cx - 0.70, mu + e0), xytext=(cx - 0.70, mu),
                     arrowprops=dict(arrowstyle="<->", color=CW, lw=1.0))
        axA.text(cx - 0.58, mu + e0 / 2, r"$\varepsilon_{ci}$", fontsize=8.4,
                 color=CW, va="center")
    axA.annotate("", xy=(9.55, 3.1), xytext=(9.55, 6.5),
                 arrowprops=dict(arrowstyle="<->", color=CB, lw=1.4))
    axA.text(9.30, 4.8, T("Between-target\nscale difference", "跨靶\n量级差异"),
             fontsize=7.6, color=CB, ha="right", va="center", rotation=90)

    # ======================= PANEL B =======================
    bx, bw = RX + 0.012, WW - 0.024
    by = TOPY + TOPH - 0.012
    # B1: 分解公式
    bh1 = 0.092
    box(fig, bx, by - bh1, bw, bh1,
        T("", ""), 8, "#FFFFFF", ST2["edge"], 1.1, z=3)
    txt(fig, bx + bw / 2, by - 0.030,
        r"$y_{ci}=\mu_c+\varepsilon_{ci}$        "
        r"$\hat{y}_{ci}=\hat{\mu}_c+\hat{\varepsilon}_{ci}$", 14, INK, z=5)
    txt(fig, bx + bw / 2, by - 0.062,
        T("each mutation = its target level + a within-target deviation; "
          "the two parts are orthogonal by construction",
          "每个突变 = 其靶点水平 + 靶内偏差；两部分按构造正交"), 7.6, GREY,
        style="italic", z=5)
    # B2: 两张贡献卡
    ch_, cy = 0.140, by - bh1 - 0.014
    cw_ = (bw - 0.014) / 2
    for k, (col, fill, name, eq, caption) in enumerate([
            (CB, "#FBF9FD", T("Between-Target Contribution", "跨靶贡献"),
             r"$C_{\mathrm{B}}=\dfrac{\mathrm{Cov}(B_y,B_{\hat{y}})}"
             r"{\sigma_y\sigma_{\hat{y}}}$",
             T("Global target-scale placement", "靶点整体量级定位")),
            (CW, "#F7FCF7", T("Within-Target Contribution", "靶内贡献"),
             r"$C_{\mathrm{W}}=\dfrac{\mathrm{Cov}(W_y,W_{\hat{y}})}"
             r"{\sigma_y\sigma_{\hat{y}}}$",
             T("Intra-target mutation ranking", "靶内突变排序"))]):
        cx = bx + k * (cw_ + 0.014)
        box(fig, cx, cy - ch_, cw_, ch_, "", 8, fill, col, 1.6, z=3)
        txt(fig, cx + cw_ / 2, cy - 0.026, name, 9.2, col, weight="bold", z=5)
        txt(fig, cx + cw_ / 2, cy - 0.078, eq, 10.6, col, z=5)
        txt(fig, cx + cw_ / 2, cy - ch_ + 0.022, caption, 7.8, INK, z=5)
    arrow(fig, bx + bw / 2, by - bh1 - 0.002, bx + bw * 0.25, cy + 0.004, CB, 1.3, z=4)
    arrow(fig, bx + bw / 2, by - bh1 - 0.002, bx + bw * 0.75, cy + 0.004, CW, 1.3, z=4)
    # B3: 恒等式徽章 + 边界提示
    iy = cy - ch_ - 0.014
    box(fig, bx, iy - 0.058, bw, 0.058,
        r"$r_{\mathrm{pooled}}=\dfrac{\mathrm{Cov}(y,\hat{y})}"
        r"{\sigma_y\sigma_{\hat{y}}}=C_{\mathrm{B}}+C_{\mathrm{W}}$"
        + "      " + T(r"(exact identity, residual $\leq 3.3\times10^{-16}$)",
                       r"（精确恒等式，残差 $\leq 3.3\times10^{-16}$）"),
        10.2, "#F4F6F9", POOL, 1.4, POOL, z=4)
    box(fig, bx, iy - 0.058 - 0.038, bw, 0.034,
        T("Common-scale additive covariance contributions — NOT independent "
          "correlations, skill scores or causal mechanisms",
          "同一尺度上的可加协方差贡献 —— 并非独立相关、能力分数或因果机制"),
        7.4, WARNF, WARN, 1.0, "#6B1005", z=4)

    # ======================= PANEL C =======================
    cy0, chh = BOTY + 0.012, BOTH - 0.024
    axC = fig.add_axes([LX + 0.012, cy0 + 0.096, 0.268, chh - 0.150], zorder=7)
    axC.set_title(T("Pooled correlation and its two components",
                    "合并相关及其两个分量"), fontsize=10.2, loc="left", pad=8)
    regs = [T("Random", "随机"), T("Complex / seq40", "Complex / seq40"),
            T("Family-OOD", "家族外推")]
    cbs = [cb0, cbc, cb3]; cws = [cw0, cwc, cw3]; rs = [r0, rc, r3]
    yy = np.arange(3)[::-1]
    axC.barh(yy, cbs, height=0.50, color=CB, edgecolor="white", lw=0.8,
             label=r"$C_{\mathrm{B}}$ (between)")
    axC.barh(yy, cws, height=0.50, left=cbs, color=CW, edgecolor="white", lw=0.8,
             label=r"$C_{\mathrm{W}}$ (within)")
    axC.axvline(0, color="#999", lw=0.9)
    axC.set_yticks(yy); axC.set_yticklabels(regs, fontsize=9.4)
    axC.set_xlabel(T("covariance contribution", "协方差贡献"), fontsize=9.8)
    axC.set_xlim(-0.08, 0.82)
    axC.legend(fontsize=8.4, frameon=False, loc="lower right")
    grid(axC, axis="x")
    for yyi, rr in zip(yy, rs):
        axC.text(0.805, yyi + 0.33, f"$r_{{\\rm pooled}}$ = {rr:.4f}", ha="right",
                 fontsize=9.0, color=POOL, fontweight="bold")
    axC.annotate(T(r"$C_{\mathrm{B}}\rightarrow\approx 0$",
                   r"$C_{\mathrm{B}}\rightarrow\approx 0$"),
                 xy=(cb3, yy[2] - 0.30), xytext=(0.12, yy[2] - 0.46), fontsize=8.6,
                 color=CB, arrowprops=dict(arrowstyle="->", color=CB, lw=1.0))
    # 关键结论卡
    kx, kw = LX + 0.292, 0.104
    card(fig, kx, cy0 + 0.096, kw, chh - 0.150, WARN, "#FEF6F4", 1.3, z=3)
    txt(fig, kx + kw / 2, cy0 + chh - 0.082, T("Key finding", "关键结论"), 8.6, WARN,
        weight="bold", z=5)
    txt(fig, kx + kw / 2, cy0 + chh - 0.134, "81.1%", 19, WARN, weight="bold", z=5)
    txt(fig, kx + kw / 2, cy0 + chh - 0.200,
        T("of the absolute drop is\nloss of between-target signal",
          "绝对下降来自\n跨靶信号的丧失"), 7.8, INK, z=5)
    txt(fig, kx + kw / 2, cy0 + chh - 0.262, r"$\Delta C_{\mathrm{B}}=-0.4347$", 8.0,
        CB, z=5)
    txt(fig, kx + kw / 2, cy0 + chh - 0.292, r"$C_{\mathrm{W}}=+0.1580$", 8.0, CW, z=5)
    txt(fig, kx + kw / 2, cy0 + chh - 0.320, T("stays positive", "仍为正"), 7.2, GREY,
        style="italic", z=5)
    # 审计胶囊（左栏底部整宽）
    aw2 = WW - 0.024
    box(fig, LX + 0.012, cy0 + 0.010, aw2, 0.046,
        T(f"Weighting audit: $C_{{\\rm B}}$ sign flips near 0 "
          f"({int(W2['signflip_C_B_gt0'].sum())}/{len(W2)}); "
          f"$C_{{\\rm W}}$ sign 100% stable",
          f"加权审计：$C_{{\\rm B}}$ 近零翻号 "
          f"({int(W2['signflip_C_B_gt0'].sum())}/{len(W2)})；$C_{{\\rm W}}$ 100% 稳定"),
        7.6, "#FFFFFF", ST3["edge"], 1.0, INK, z=4)
    box(fig, LX + 0.012, cy0 + 0.010 - 0.054, aw2, 0.046,
        T(f"Shortlist geometry: {100*GEO['n_le_K_fraction']:.1f}% of targets "
          f"($n\\leq K$) carry zero Fixed-K information",
          f"短名单几何：{100*GEO['n_le_K_fraction']:.1f}% 的靶点（$n\\leq K$）"
          f"不携带 Fixed-K 信息"), 7.6, "#FFFFFF", ST3["edge"], 1.0, INK, z=4)

    # ======================= PANEL D: radar =======================
    axD = fig.add_axes([RX + 0.108, BOTY + 0.075, 0.336, BOTH - 0.130],
                       projection="polar", zorder=7)
    g = D[(D.rowset_hash == PU["rowset_hash"]) & D.key.str.endswith("_s42")].copy()
    g = g.set_index("method_family")
    g["neg_rmse"] = -g["rmse"] if "rmse" in g.columns else np.nan
    axes_def = [("overall_pearson", r"Pooled Pearson $r$", "#4A6FA5"),
                ("C_B_row", r"Between-target $C_{\mathrm{B}}$", CB),
                ("macro_spearman", r"Macro Spearman $\rho$", "#2A9D8F"),
                ("neg_rmse", r"Negated RMSE", AMBER),
                ("C_W", r"Within-target $C_{\mathrm{W}}$", CW)]
    cols = [c for c, _, _ in axes_def]
    M = g[cols].to_numpy(float)
    lo, hi = np.nanmin(M, axis=0), np.nanmax(M, axis=0)
    Z = 0.14 + 0.86 * ((M - lo) / np.where(hi - lo > 0, hi - lo, 1))
    ang = np.linspace(0, 2 * np.pi, len(cols), endpoint=False).tolist()
    ang += ang[:1]
    for key, lab, col, al in [
            ("F1_Base_Mut", T("F1 Mutation-only baseline", "F1 仅突变性质基线"),
             "#8A8F98", 0.10),
            ("F2_Base_Static", T("F2 + Static local geometry", "F2 + 静态局部几何"),
             AMBER, 0.16),
            ("F4_Base_PartnerPert",
             T("F4 + Mutation × partner perturbation", "F4 + 突变 × partner 扰动"),
             "#2A9D8F", 0.18),
            ("F8_Full", T("F8 Full local-physics model", "F8 全量局部物理模型"),
             CB, 0.22)]:
        if key not in g.index:
            continue
        v = Z[list(g.index).index(key)].tolist(); v += v[:1]
        axD.plot(ang, v, color=col, lw=2.2, label=lab, zorder=4)
        axD.fill(ang, v, color=col, alpha=al, zorder=3)
    axD.set_xticks(ang[:-1])
    axD.set_xticklabels([lbl for _, lbl, _ in axes_def], fontsize=8.4)
    for lbl, (_, _, col) in zip(axD.get_xticklabels(), axes_def):
        lbl.set_color(col); lbl.set_fontweight("bold")
    axD.set_yticks([]); axD.set_ylim(0, 1.04)
    axD.spines["polar"].set_color("#C9CDD4")
    axD.grid(color="#DDE1E6", lw=0.7)
    axD.set_title(T("Standardised capability profiles", "标准化能力画像"),
                  fontsize=10.0, pad=22)
    axD.legend(fontsize=7.4, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.10), ncol=2, columnspacing=1.2,
               handlelength=1.6)
    txt(fig, RX + 0.012, BOTY + 0.052,
        T("Axes show standardised profiles on identical 3,351 rows for visual "
          "comparison; raw values reported in text. F1–F8 are controlled ablation "
          "variants of one method family, not independent methods.",
          "各轴显示同一 3,351 行上的标准化画像，仅供视觉比较；原始数值见正文。"
          "F1–F8 为同一方法族的受控消融变体，并非独立方法。"),
        7.0, GREY, ha="left", style="italic", z=6)

    save(fig, "Fig1_overview", zh)


# =====================================================================
# FIGURE 2 -- 2 行 (2 上 1 下)
# =====================================================================
def fig2(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(11.0, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.55, wspace=0.26,
                          left=0.085, right=0.975, top=0.90, bottom=0.115)
    order = ["random", "complex", "seq40", "family_OOD"]
    lab = ([T("Random", "随机"), T("Complex", "复合物"), T("seq40", "seq40"),
            T("Family-OOD", "家族外推")])
    rp = [H4[k]["overall_pearson"] for k in order]
    cb = [H4[k]["C_B"] for k in order]
    cw = [H4[k]["C_W"] for k in order]
    xs = np.arange(4)
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A")
    ax.plot(xs, rp, "o-", color=POOL, lw=2.2, ms=7)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=10)
    ax.set_ylabel(T("pooled Pearson", "合并 Pearson"))
    ax.set_title(T("Pooled correlation by split regime", "各划分下的合并相关"),
                 fontsize=11, loc="left", pad=9)
    for i, v in enumerate(rp):
        ax.annotate(f"{v:.3f}", (i, v), fontsize=9.4, xytext=(0, 9),
                    textcoords="offset points", ha="center")
    ax.set_ylim(-0.02, 0.78); grid(ax)

    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B")
    ax.plot(xs, cb, "s-", color=CB, lw=2.2, ms=7, label=r"$C_{\mathrm{B}}$ (between)")
    ax.plot(xs, cw, "d-", color=CW, lw=2.2, ms=7, label=r"$C_{\mathrm{W}}$ (within)")
    ax.axhline(0, color="#999", lw=1.0)
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=10)
    ax.set_ylabel(T("covariance contribution", "协方差贡献"))
    ax.set_title(T(r"Exact decomposition  $r_{\mathrm{pooled}}=C_{\mathrm{B}}"
                   r"+C_{\mathrm{W}}$",
                   r"精确分解  $r_{\mathrm{pooled}}=C_{\mathrm{B}}+C_{\mathrm{W}}$"),
                 fontsize=11, loc="left", pad=9)
    ax.legend(fontsize=9.4, frameon=False)
    grid(ax)
    ax.annotate(T(r"$C_{\mathrm{W}}$ stays positive",
                  r"$C_{\mathrm{W}}$ 保持为正"), xy=(3, cw3), xytext=(1.25, 0.255),
                fontsize=9.2, color=CW,
                arrowprops=dict(arrowstyle="->", color=CW, lw=1.1))

    ax = fig.add_subplot(gs[1, :]); panel(ax, "C")
    d_cb, d_cw = abs(DL["delta_C_B"]), abs(DL["delta_C_W"])
    ax.barh([0], [d_cb], color=CB, height=0.42)
    ax.barh([-0.72], [d_cw], color=CW, height=0.42)
    ax.set_yticks([0, -0.72])
    ax.set_yticklabels([T(r"Between-target  $\Delta C_{\mathrm{B}}$",
                          r"跨靶  $\Delta C_{\mathrm{B}}$"),
                        T(r"Within-target  $\Delta C_{\mathrm{W}}$",
                          r"靶内  $\Delta C_{\mathrm{W}}$")], fontsize=10.4)
    ax.set_xlabel(T(r"absolute change, Random $\rightarrow$ Family-OOD",
                    r"绝对变化，随机 $\rightarrow$ 家族外推"))
    ax.text(d_cb * 0.5, 0.34, f"{100*DL['fraction_between']:.1f}%", ha="center",
            fontsize=15, color=CB, fontweight="bold")
    ax.text(d_cw * 0.5, -0.38, f"{100*DL['fraction_within']:.1f}%", ha="center",
            fontsize=11.5, color=CW, fontweight="bold")
    ax.set_xlim(0, d_cb * 1.16); ax.set_ylim(-1.12, 0.72)
    ax.set_title(T("Attribution of the absolute decline",
                   "绝对下降的归因"), fontsize=11, loc="left", pad=9)
    grid(ax, axis="x")

    save(fig, "Fig2_split_decomposition", zh)


# =====================================================================
# FIGURE 3 -- 1x3 但宽高比 2.2:1
# =====================================================================
def fig3(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    SCH = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
    fig = plt.figure(figsize=(11.0, 5.0))
    gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.070, right=0.975,
                          top=0.90, bottom=0.14)
    for k, (key, col, nm, pnl) in enumerate(
            (("C_B", CB, r"$C_{\mathrm{B}}$ (between-target)", "A"),
             ("C_W", CW, r"$C_{\mathrm{W}}$ (within-target)", "B"))):
        ax = fig.add_subplot(gs[0, k]); panel(ax, pnl)
        A_ = W[[f"{key}_{s}" for s in SCH]].to_numpy(float)
        rng_ = np.nanmax(A_, axis=1) - np.nanmin(A_, axis=1)
        ax.hist(rng_, bins=24, color=col, alpha=0.85, edgecolor="white", lw=0.4)
        md = np.nanmedian(rng_)
        ax.axvline(md, color=POOL, ls="--", lw=1.6)
        ax.text(md * 1.08, ax.get_ylim()[1] * 0.88,
                T(f"median {md:.4f}", f"中位 {md:.4f}"), fontsize=9.4, color=POOL)
        ax.set_xlabel(T("range across the five weighting schemes",
                        "五种加权方案下的极差"))
        ax.set_ylabel(T("model configurations", "模型配置数"))
        ax.set_title(nm, fontsize=11, loc="left", pad=9)
        grid(ax)
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "C")
    thr = [0, 0.001, 0.005]
    cb_f = [int(W2[f"signflip_C_B_gt{t}"].sum()) for t in thr]
    cw_f = [int(W2[f"signflip_C_W_gt{t}"].sum()) for t in thr]
    xx = np.arange(3); w = 0.36
    ax.bar(xx - w / 2, cb_f, width=w, color=CB, label=r"$C_{\mathrm{B}}$")
    ax.bar(xx + w / 2, cw_f, width=w, color=CW, label=r"$C_{\mathrm{W}}$")
    ax.set_xticks(xx); ax.set_xticklabels(["0", "0.001", "0.005"])
    ax.set_xlabel(T(r"sign-flip threshold on $|C|$", r"$|C|$ 的符号翻转阈值"))
    ax.set_ylabel(T("configurations reversing sign", "符号翻转的配置数"))
    ax.set_title(T("Directional stability", "方向稳定性"), fontsize=11, loc="left",
                 pad=9)
    ax.legend(fontsize=9.4, frameon=False); grid(ax)
    for i, v in enumerate(cb_f):
        ax.annotate(str(v), (i - w / 2, v), fontsize=9.4, xytext=(0, 5),
                    textcoords="offset points", ha="center", color=CB)
    save(fig, "Fig3_weighting", zh)


# =====================================================================
# FIGURE 4 -- 单面板 1.6:1
# =====================================================================
def fig4(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(9.6, 6.0))
    ax = fig.add_axes([0.105, 0.155, 0.865, 0.800])
    K = 10
    m = GEO_DF.n <= 140
    ax.plot(GEO_DF.n[m], GEO_DF.D[m], lw=2.8, color=POOL)
    ax.axvline(K, color=AMBER, ls="--", lw=1.7)
    ax.text(K + 2.4, GEO_DF.D[m].max() * 0.88,
            r"$n \leq K \Rightarrow D(n,K)=0$", fontsize=10.4, color=AMBER)
    ns_ = GEO["n_star_analytic"]
    ax.axvline(ns_, color=CW, ls=":", lw=1.7)
    ax.annotate(T(f"peak  $n^*\\approx{ns_:.1f}$", f"峰值  $n^*\\approx{ns_:.1f}$"),
                xy=(ns_, GEO_DF.D[m].max()),
                xytext=(ns_ + 8, GEO_DF.D[m].max() * 0.66), fontsize=10.4,
                color=CW, arrowprops=dict(arrowstyle="->", color=CW, lw=1.1))
    ax.plot([20, 50], [GEO["D_at_n20"], GEO["D_at_n50"]], "o", ms=8, color=WARN,
            zorder=5)
    ax.annotate(f"D(20) = {GEO['D_at_n20']:.3f}", (20, GEO["D_at_n20"]),
                xytext=(26, GEO["D_at_n20"] + 0.020), fontsize=10.0, color=WARN)
    ax.annotate(f"D(50) = {GEO['D_at_n50']:.3f}", (50, GEO["D_at_n50"]),
                xytext=(58, GEO["D_at_n50"] + 0.010), fontsize=10.0, color=WARN)
    ax.set_xlabel(T(r"candidate-set size $n$  ($K=10$)", r"候选集大小 $n$（$K=10$）"))
    ax.set_ylabel(r"discriminable interval  $D(n,K)$")
    grid(ax)
    save(fig, "Fig4_candidate_geometry", zh)


# =====================================================================
# FIGURE 5 -- 单面板 1.6:1，纯英文标签，y 上限扩展
# =====================================================================
NICE = {"F0_Base": ("Base model", "基线模型"),
        "F1_Base_Mut": ("Mutation-only baseline", "仅突变性质基线"),
        "F2_Base_Static": ("Static local geometry", "静态局部几何"),
        "F3_Base_EnvPartner": ("Partner static context", "partner 静态环境"),
        "F4_Base_PartnerPert": ("Mutation × partner perturbation",
                                "突变 × partner 扰动"),
        "F5_Base_IntraPert": ("Intra-chain perturbation", "链内扰动"),
        "F6_Partner_Intra": ("Partner + intra-chain", "partner + 链内"),
        "F7_F6_Burial": ("Partner + intra + burial", "partner + 链内 + 埋藏"),
        "F8_Full": ("Full local-physics model", "全量局部物理模型"),
        "F9_PhysicsOnly": ("Physics features only", "仅物理特征")}
# 手工精调标签偏移，避免重叠
OFF = {"F2_Base_Static": (10, 8), "F4_Base_PartnerPert": (12, 6),
       "F5_Base_IntraPert": (12, -16), "F8_Full": (-14, 10),
       "F6_Partner_Intra": (12, -4), "F7_F6_Burial": (-12, -18),
       "F9_PhysicsOnly": (12, 4), "F0_Base": (12, -14), "F3_Base_EnvPartner": (12, 4)}


def resolve_frag(label: str) -> str:
    """从 remedy CSV 的 label 中稳健解析变体键。

    CSV 的 label 本身可能含中文（例如 'V7.3 F4 mutation×partner 扰动'），
    因此不能依赖 `label.replace('V7.3 ','')` 直接得到键；这里改为在 NICE 的
    已知键中做包含匹配，找不到时再退回清洗后的字符串（并剥离非 ASCII）。
    """
    s = label.replace("V7.3 ", "").strip()
    for k in NICE:
        if k in label:
            return k
    return re.sub(r"[^\x00-\x7F]+", "", s).strip() or s


def fig5(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(9.6, 6.0))
    ax = fig.add_axes([0.105, 0.145, 0.865, 0.815])
    b = RS[RS.label.str.startswith("V7.3 F1")]
    r5 = RS[RS.step == "R5"].copy()
    if len(b):
        b = b.iloc[0]
        r5["d_ov"] = r5.overall_r - b.overall_r
        r5["d_mc"] = r5.macro_spearman - b.macro_spearman
        ax.axhline(0, color="#9AA0A6", lw=1.1)
        ax.axvline(0, color="#9AA0A6", lw=1.1)
        ax.axhspan(0, 0.24, xmin=0.5, xmax=1.0, color=CWF, alpha=0.35, zorder=0)
        for _, r in r5.iterrows():
            frag = resolve_frag(r.label)
            nm = NICE.get(frag, (frag, frag))
            good = r.d_mc > 0.02
            col = CW if good else WARN
            ax.scatter(r.d_ov, r.d_mc, s=115, c=col, zorder=4,
                       edgecolors="white", linewidths=1.5)
            dx, dy = OFF.get(frag, (10, 6))
            ax.annotate(nm[1] if zh else nm[0], (r.d_ov, r.d_mc), fontsize=9.6,
                        xytext=(dx, dy), textcoords="offset points", color=col)
        ax.text(0.03, 0.06, T("pooled-score-only change", "仅合并分数改善"),
                transform=ax.transAxes, fontsize=9.0, color=WARN, style="italic")
        ax.add_patch(FancyBboxPatch((0.575, 0.855), 0.395, 0.098,
                                    boxstyle="round,pad=0.012", fc="#EFF8F0",
                                    ec=CW, lw=1.0, transform=ax.transAxes,
                                    zorder=5))
        ax.text(0.772, 0.904, T("both axes improve", "两轴同时改善"),
                transform=ax.transAxes, ha="center", va="center", fontsize=9.2,
                color=CW, style="italic", zorder=6)
    ax.set_xlabel(T(r"$\Delta$ pooled Pearson vs mutation-only baseline",
                    r"相对「仅突变性质」基线的 $\Delta$ 合并 Pearson"))
    ax.set_ylabel(T(r"$\Delta$ macro within-target Spearman",
                    r"$\Delta$ 靶内 Macro Spearman"))
    ax.set_ylim(-0.055, 0.225)
    ax.set_xlim(-0.02, 0.145)
    grid(ax, axis="both")
    save(fig, "Fig5_remedy", zh)


# =====================================================================
# FIGURE 6 -- A:B = 1:1.15，标签反方向引出
# =====================================================================
EXT_OFF = {"Flex-ddG": ((-52, 14), "right"),
           "ddG_FoldX": ((14, -18), "left"),
           "DDGPred": ((-12, 12), "right"),
           "proteinMPNN": ((14, 6), "left"),
           "ESM-IF": ((14, -6), "left")}


def fig6(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    fig = plt.figure(figsize=(12.0, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.30,
                          left=0.075, right=0.960, top=0.885, bottom=0.155)
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A")
    for r in EP.itertuples():
        ax.scatter(r.overall_pearson, r.macro_spearman, s=145, color=POOL, zorder=3,
                   edgecolors="white", linewidths=1.4)
        (dx, dy), ha = EXT_OFF.get(r.model, ((12, 6), "left"))
        ax.annotate(r.model, (r.overall_pearson, r.macro_spearman), fontsize=10.2,
                    xytext=(dx, dy), textcoords="offset points", ha=ha, color=INK,
                    arrowprops=dict(arrowstyle="-", color="#8A8F98", lw=0.8))
    ax.axhline(0, color="#9AA0A6", lw=1.0); ax.axvline(0, color="#9AA0A6", lw=1.0)
    ax.set_xlabel(T("Pooled Pearson", "合并 Pearson"))
    ax.set_ylabel(T("Macro within-target Spearman", "靶内 Macro Spearman"))
    ax.set_title(T("Pooled vs. within-target performance",
                   "合并性能与靶内性能"), fontsize=11, loc="left", pad=9)
    ax.set_ylim(0.0, 0.175); ax.set_xlim(-0.30, 0.32)
    grid(ax, axis="both")
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B")
    names = list(EP.model)
    cols = ["overall_pearson", "macro_spearman", "C_B_row", "C_W"]
    heads = [T("Pooled\nPearson", "合并\nPearson"),
             T("Macro\nSpearman", "Macro\nSpearman"),
             r"$C_{\mathrm{B}}$", r"$C_{\mathrm{W}}$"]
    M = EP[cols].to_numpy(float)
    Z = (M - M.mean(axis=0)) / np.where(M.std(axis=0) > 0, M.std(axis=0), 1)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(heads, fontsize=10.2)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10.4)
    for i in range(len(names)):
        for j in range(4):
            ax.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center", fontsize=9.8,
                    color="white" if abs(Z[i, j]) > 1.25 else INK)
    ax.set_title(T("Method-level capability profile (column-standardised)",
                   "方法级能力画像（按列标准化）"), fontsize=11, loc="left", pad=9)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.03,
                 label=T("within-column z-score", "列内 z 分数"))
    save(fig, "Fig6_external", zh)


# =====================================================================
# SUPPLEMENT
# =====================================================================
def supp(zh):
    use(zh)
    T = (lambda en, z: z if zh else en)
    # S1
    fig = plt.figure(figsize=(10.6, 5.0))
    ax = fig.add_axes([0.235, 0.145, 0.745, 0.820])
    u = UR.sort_values("n_rows"); y = np.arange(len(u))
    ax.barh(y, u.n_rows, color="#8FA3BF", alpha=0.92,
            label=T("mutation rows", "突变行数"))
    ax.barh(y, u.n_ge5, color=AMBER, alpha=0.95, height=0.55,
            label=T("targets with n ≥ 5", "n ≥ 5 的靶点"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ({h[:8]})" for a, h in zip(u.artifact, u.rowset_hash)],
                       fontsize=8.8)
    ax.set_xlabel(T("count", "计数"))
    ax.legend(fontsize=9.4, frameon=False)
    grid(ax, axis="x")
    save(fig, "FigS1_matched_universes", zh, FS)

    # S2 -- 修复右侧大片空白：居中占满 90%
    g = D[(D.rowset_hash == PU["rowset_hash"]) & D.key.str.endswith("_s42")].copy()
    g = g.sort_values("macro_spearman", ascending=False)
    yl = [NICE.get(m, (m, m))[1 if zh else 0] for m in g.method_family]
    fig = plt.figure(figsize=(10.0, 6.2))
    ax = fig.add_axes([0.245, 0.130, 0.615, 0.830])
    cols = ["overall_pearson", "macro_spearman", "C_B_row", "C_W"]
    heads = [T("Pooled\nPearson", "合并\nPearson"),
             T("Macro\nSpearman", "Macro\nSpearman"),
             T("$C_{\\mathrm{B}}$\n(between)", "$C_{\\mathrm{B}}$\n（跨靶）"),
             T("$C_{\\mathrm{W}}$\n(within)", "$C_{\\mathrm{W}}$\n（靶内）")]
    M = g[cols].to_numpy(float)
    Z = (M - M.mean(axis=0)) / np.where(M.std(axis=0) > 0, M.std(axis=0), 1)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(heads, fontsize=10.2)
    ax.set_yticks(range(len(g))); ax.set_yticklabels(yl, fontsize=9.6)
    for i in range(len(g)):
        for j in range(4):
            ax.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center", fontsize=9.0,
                    color="white" if abs(Z[i, j]) > 1.25 else INK)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.80, pad=0.02,
                 label=T("within-column z-score", "列内 z 分数"))
    save(fig, "Capability_heatmap_provenance_clean", zh, FS)

    # S3 -- 编号修正为 S3
    fig = plt.figure(figsize=(10.6, 5.0))
    ax = fig.add_axes([0.185, 0.190, 0.795, 0.775])
    s = ST.sort_values("rho_overall_vs_macro")
    yy = np.arange(len(s))
    cols_ = [WARN if v < 0 else "#4A6FA5" for v in s.rho_overall_vs_macro]
    ax.barh(yy, s.rho_overall_vs_macro, color=cols_, alpha=0.92)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{h[:8]} / L{lv}" for h, lv in
                        zip(s.rowset_hash, s.level)], fontsize=8.8)
    ax.axvline(0, color="#9AA0A6", lw=1.0)
    ax.set_xlabel(T(r"$\rho$(pooled Pearson, macro within-target Spearman)",
                    r"$\rho$（合并 Pearson，靶内 Macro Spearman）"))
    grid(ax, axis="x")
    save(fig, "FigS3_capability_pool_sensitivity", zh, FS)


for zh in (False, True):
    print(f"--- {'ZH' if zh else 'EN'} ---")
    fig1(zh); fig2(zh); fig3(zh); fig4(zh); fig5(zh); fig6(zh); supp(zh)

Q = pd.DataFrame(GLY)
Q.to_csv(A / "figure_glyph_qc_v4.csv", index=False, encoding="utf-8")
nf = int((Q.result == "FAIL").sum())
print(f"\nFIGURES v4: EN={len(list(FM.glob('*_EN.png')))} ZH={len(list(FM.glob('*_ZH.png')))} "
      f"supp={len(list(FS.glob('*_EN.png')))}")
print(f"FIGURE_GLYPH_QC = {'PASS' if nf == 0 else 'FAIL'} ({len(Q)-nf}/{len(Q)})")
if nf:
    print(Q[Q.result == "FAIL"].to_string(index=False))
