#!/usr/bin/env python3
"""
CR-PPI-ddG publication figures — vector-first plotting code
============================================================

Outputs SVG + PDF + high-resolution PNG for:
  Figure 1, 2, 3, 4, 5, 6, S1, S2, S3

Notes
-----
1. All numerical values encoded below are taken from the current manuscript /
   supplementary figures and tables supplied in the working draft.
2. Figure 3 is reconstructed from the published histogram bar counts visible in
   the supplied figure (each histogram sums to N=109); the exact raw per-model
   observations were not available in the PDF. If you have the original raw
   configuration-level ranges, replace FIG3_CB_CENTERS/COUNTS and
   FIG3_CW_CENTERS/COUNTS with the raw arrays and call np.histogram instead.
3. SVG text is kept as text (`svg.fonttype = none`), and PDFs embed TrueType
   fonts (`pdf.fonttype = 42`) for editing in Illustrator/Inkscape.
"""

from pathlib import Path
import argparse
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, Wedge, FancyArrowPatch

OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)

# ---------- publication style ----------
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.0,
    "axes.titlesize": 11,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.2,
    "ytick.labelsize": 8.2,
    "legend.fontsize": 8.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "lines.linewidth": 1.8,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "mathtext.fontset": "dejavusans",
})

NAVY = "#172033"
PURPLE = "#7956A8"
PURPLE2 = "#A98CCB"
GREEN = "#2E8B62"
GREEN2 = "#7BC49A"
BLUE = "#3F7DB7"
ORANGE = "#D9913C"
CORAL = "#C84B3A"
GRAY = "#7B8492"
GRID = "#D9DEE8"
BG = "#F7F8FB"
LIGHT_PURPLE = "#EEE7F6"
LIGHT_GREEN = "#E7F3EC"
LIGHT_BLUE = "#EAF2F9"
LIGHT_CORAL = "#FAE9E6"
LIGHT_ORANGE = "#FFF3E4"

def save_all(fig, stem, dpi=420):
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def panel_box(ax, title, letter, header_color=LIGHT_PURPLE):
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch(
        (0.002, 0.002), 0.996, 0.996,
        boxstyle="round,pad=0.008,rounding_size=0.025",
        fc="white", ec="#E1E5EC", lw=1.0, zorder=-20
    ))
    ax.add_patch(FancyBboxPatch(
        (0.002, 0.895), 0.996, 0.103,
        boxstyle="round,pad=0.008,rounding_size=0.025",
        fc=header_color, ec="none", zorder=-10
    ))
    ax.text(0.035, 0.947, letter, va="center", ha="center",
            fontsize=17, fontweight="bold", color=NAVY)
    ax.text(0.085, 0.947, title, va="center", ha="left",
            fontsize=11.0, fontweight="bold", color=NAVY)

def rounded(ax, xy, wh, fc="white", ec=GRID, lw=1.0, r=0.02, z=0):
    x, y = xy; w, h = wh
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0.008,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z)
    ax.add_patch(p)
    return p

def arrow(ax, p0, p1, color=GRAY, lw=1.6, ms=12):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                                 color=color, lw=lw))

def abstract_complex(ax, center=(0.5, 0.5), scale=1.0):
    """Small schematic PPI complex; intentionally illustrative, not structural data."""
    cx, cy = center
    rng = np.random.default_rng(2)
    for color, dx in [(BLUE, -0.035*scale), (CORAL, 0.035*scale)]:
        pts = rng.normal(size=(18, 2))
        pts[:, 0] = cx + dx + pts[:,0]*0.022*scale
        pts[:, 1] = cy + pts[:,1]*0.030*scale
        ax.scatter(pts[:,0], pts[:,1], s=35*scale, c=color, alpha=0.75,
                   edgecolors="white", linewidths=0.3, zorder=3)

# ============================================================================
# FIGURE 1 — manuscript-optimized, larger-text 2×2 layout
# ============================================================================

def figure1():
    fig = plt.figure(figsize=(14.4, 9.0))
    gs = fig.add_gridspec(2, 2, wspace=0.035, hspace=0.055)

    # --- A ---
    ax = fig.add_subplot(gs[0,0])
    panel_box(ax, "Provenance-filtered, fingerprint-matched benchmark construction", "A", LIGHT_CORAL)

    # source
    rounded(ax, (0.03,0.61),(0.20,0.23), fc="#F1ECFA", ec="#D8C8EC")
    ax.add_patch(Ellipse((0.13,0.80),0.13,0.035,fc=PURPLE2,ec=PURPLE,lw=1))
    ax.add_patch(Rectangle((0.065,0.675),0.13,0.125,fc="#DDD0EF",ec=PURPLE,lw=1))
    ax.add_patch(Ellipse((0.13,0.675),0.13,0.035,fc="#D2C1E8",ec=PURPLE,lw=1))
    ax.text(0.13,0.735,"SKEMPI\nv2.0",ha="center",va="center",fontsize=12,fontweight="bold",color=PURPLE)
    ax.text(0.13,0.635,"experimentally measured\nPPI ΔΔG mutations",ha="center",va="top",fontsize=7.7,color=NAVY)

    # source counts
    cards = [(0.27,"4,473","rows",PURPLE),(0.43,"214","targets",BLUE),(0.59,"108","families",ORANGE)]
    for x, num, lab, c in cards:
        rounded(ax,(x,0.705),(0.135,0.105),fc="white",ec="#E2E5EA")
        ax.text(x+0.0675,0.765,num,ha="center",va="center",fontsize=12,fontweight="bold",color=c)
        ax.text(x+0.0675,0.725,lab,ha="center",va="center",fontsize=7.8,color=NAVY)

    abstract_complex(ax, (0.84,0.745), 1.5)
    ax.text(0.84,0.635,"PPI target + mutation\nmeasurements",ha="center",va="center",fontsize=8,color=NAVY)

    # predictions
    rounded(ax,(0.03,0.42),(0.55,0.145),fc=LIGHT_BLUE,ec="#D5E4F2")
    ax.text(0.045,0.54,"Predictions from diverse method families",ha="left",va="center",
            fontsize=9,fontweight="bold",color=NAVY)
    family_labels=[("Physics",BLUE),("Structure ML",PURPLE),("Protein LM",GREEN),("Graph NN",ORANGE)]
    for j,(lab,c) in enumerate(family_labels):
        x=0.055+j*0.125
        rounded(ax,(x,0.445),(0.105,0.06),fc="white",ec=c,lw=1)
        ax.text(x+0.0525,0.475,lab,ha="center",va="center",fontsize=7.1,color=c,fontweight="bold")
    ax.text(0.555,0.474,"…",fontsize=18,ha="center",va="center",color=GRAY)

    # filters
    rounded(ax,(0.62,0.42),(0.34,0.145),fc=LIGHT_CORAL,ec="#EBCBC6")
    ax.text(0.79,0.54,"Provenance & leakage audit",ha="center",va="center",
            fontsize=9,fontweight="bold",color=NAVY)
    ax.text(0.64,0.493,"✓ SHA-256 row-set fingerprints",fontsize=7.7,color=GREEN,ha="left")
    ax.text(0.64,0.462,"✓ retain provenance-clean Family-OOD",fontsize=7.7,color=GREEN,ha="left")
    ax.text(0.64,0.431,"✕ remove training overlap / permissive splits",fontsize=7.5,color=CORAL,ha="left")

    arrow(ax,(0.50,0.37),(0.50,0.31),PURPLE)

    # matched universe
    rounded(ax,(0.18,0.07),(0.64,0.22),fc=LIGHT_PURPLE,ec="#CDBBE4",lw=1.2)
    ax.text(0.50,0.255,"Primary matched evaluation universe",ha="center",va="center",
            fontsize=10.5,fontweight="bold",color=PURPLE)
    ax.text(0.29,0.165,"3,351",ha="center",fontsize=17,fontweight="bold",color=PURPLE)
    ax.text(0.29,0.122,"mutation rows",ha="center",fontsize=8,color=NAVY)
    ax.text(0.50,0.165,"102",ha="center",fontsize=17,fontweight="bold",color=BLUE)
    ax.text(0.50,0.122,"targets (n ≥ 5)",ha="center",fontsize=8,color=NAVY)
    ax.text(0.71,0.165,"10",ha="center",fontsize=17,fontweight="bold",color=GREEN)
    ax.text(0.71,0.122,"method families",ha="center",fontsize=8,color=NAVY)

    # --- B ---
    ax = fig.add_subplot(gs[0,1])
    panel_box(ax, "Exact covariance-signal decomposition", "B", LIGHT_PURPLE)

    rounded(ax,(0.06,0.72),(0.40,0.13),fc=LIGHT_PURPLE,ec="#D7C8E9")
    ax.text(0.26,0.815,"True ΔΔG",ha="center",fontsize=10,fontweight="bold",color=NAVY)
    ax.text(0.26,0.765,r"$y_{ci}=\mu_c+\epsilon_{ci}$",ha="center",fontsize=15,color=PURPLE)
    ax.text(0.26,0.73,r"$\mu_c$: target mean   •   $\epsilon_{ci}$: within-target deviation",
            ha="center",fontsize=7.3,color=NAVY)

    rounded(ax,(0.54,0.72),(0.40,0.13),fc=LIGHT_BLUE,ec="#D5E4F2")
    ax.text(0.74,0.815,"Predicted ΔΔG",ha="center",fontsize=10,fontweight="bold",color=NAVY)
    ax.text(0.74,0.765,r"$\hat y_{ci}=\hat\mu_c+\hat\epsilon_{ci}$",ha="center",fontsize=15,color=BLUE)
    ax.text(0.74,0.73,r"$\hat\mu_c$: predicted target mean   •   $\hat\epsilon_{ci}$: predicted deviation",
            ha="center",fontsize=7.1,color=NAVY)

    rounded(ax,(0.18,0.54),(0.64,0.11),fc="white",ec="#C9B5E2",lw=1.3)
    ax.text(0.50,0.595,r"$r_{\mathrm{pooled}}=C_B+C_W$",ha="center",va="center",
            fontsize=21,color=NAVY,fontweight="bold")

    rounded(ax,(0.06,0.24),(0.40,0.23),fc=LIGHT_PURPLE,ec="#D7C8E9")
    ax.text(0.26,0.435,r"$C_B$",ha="center",fontsize=17,fontweight="bold",color=PURPLE)
    ax.text(0.26,0.395,"between-target signal",ha="center",fontsize=9.5,fontweight="bold",color=NAVY)
    ax.text(0.26,0.36,"target-level scale placement",ha="center",fontsize=7.7,color=GRAY)
    # schematic scatter
    xs=np.array([0.11,0.16,0.22,0.28,0.35,0.40])
    ys=np.array([0.285,0.30,0.32,0.345,0.36,0.385])
    ax.plot(xs,ys,color=PURPLE,lw=1.5); ax.scatter(xs,ys,s=18,color=PURPLE)
    ax.text(0.26,0.265,"correlation of target means",ha="center",fontsize=7.3,color=PURPLE)

    rounded(ax,(0.54,0.24),(0.40,0.23),fc=LIGHT_GREEN,ec="#CBE4D5")
    ax.text(0.74,0.435,r"$C_W$",ha="center",fontsize=17,fontweight="bold",color=GREEN)
    ax.text(0.74,0.395,"within-target signal",ha="center",fontsize=9.5,fontweight="bold",color=NAVY)
    ax.text(0.74,0.36,"mutation ordering within a target",ha="center",fontsize=7.7,color=GRAY)
    xs2=np.array([0.59,0.64,0.70,0.76,0.83,0.89])
    ys2=np.array([0.285,0.30,0.335,0.325,0.37,0.39])
    ax.plot(xs2,ys2,color=GREEN,lw=1.5); ax.scatter(xs2,ys2,s=18,color=GREEN)
    ax.text(0.74,0.265,"correlation of within-target deviations",ha="center",fontsize=7.2,color=GREEN)

    rounded(ax,(0.18,0.07),(0.64,0.11),fc=LIGHT_BLUE,ec="#CADCEC")
    ax.text(0.50,0.135,r"Exact identity: residual $\leq 3.3\times10^{-16}$",
            ha="center",fontsize=9,fontweight="bold",color=BLUE)
    ax.text(0.50,0.095,"common-scale additive covariance contributions — not independent correlations",
            ha="center",fontsize=7.6,color=NAVY)

    # --- C ---
    ax = fig.add_subplot(gs[1,0])
    panel_box(ax, "Controlled split collapse on identical 7,523 rows", "C", LIGHT_BLUE)

    ax.text(0.50,0.855,"Same mutation rows and model configurations; only split strictness changes",
            ha="center",fontsize=8.8,color=NAVY)
    split_cards=[
        (0.05,"Random","permutation",0.6844,BLUE,LIGHT_BLUE),
        (0.36,"Complex / seq40","interface / sequence",0.1931,CORAL,LIGHT_CORAL),
        (0.67,"Family-OOD","between-family",0.1508,CORAL,LIGHT_CORAL),
    ]
    for x,title,sub,val,c,fc in split_cards:
        rounded(ax,(x,0.61),(0.26,0.20),fc=fc,ec="#DCE1E9")
        ax.text(x+0.13,0.775,title,ha="center",fontsize=10,fontweight="bold",color=NAVY)
        ax.text(x+0.13,0.743,f"({sub})",ha="center",fontsize=7.2,color=GRAY)
        ax.text(x+0.13,0.66,f"{val:.4f}",ha="center",fontsize=18,fontweight="bold",color=c)
        ax.text(x+0.13,0.625,"pooled Pearson",ha="center",fontsize=7.3,color=NAVY)

    ax.text(0.50,0.56,"Decomposition of performance drop: Random → Family-OOD",
            ha="center",fontsize=9.5,fontweight="bold",color=NAVY)

    # horizontal contribution bars
    ax.text(0.08,0.475,r"Between-target  $C_B$",fontsize=9.5,fontweight="bold",color=PURPLE)
    ax.add_patch(Rectangle((0.30,0.455),0.55,0.045,fc="#EEF0F3",ec="none"))
    ax.add_patch(Rectangle((0.30,0.455),0.45,0.045,fc=PURPLE2,ec="none"))
    ax.text(0.88,0.477,r"$\Delta C_B=-0.4347$",ha="right",fontsize=10,color=PURPLE,fontweight="bold")

    ax.text(0.08,0.385,r"Within-target  $C_W$",fontsize=9.5,fontweight="bold",color=GREEN)
    ax.add_patch(Rectangle((0.30,0.365),0.55,0.045,fc="#EEF0F3",ec="none"))
    ax.add_patch(Rectangle((0.30,0.365),0.105,0.045,fc=GREEN2,ec="none"))
    ax.text(0.88,0.387,r"$\Delta C_W=-0.1011$",ha="right",fontsize=10,color=GREEN,fontweight="bold")

    # donut
    center=(0.21,0.19)
    ax.add_patch(Wedge(center,0.105,90,90+360*0.811,fc=PURPLE2,ec="white"))
    ax.add_patch(Wedge(center,0.105,90+360*0.811,450,fc=GREEN2,ec="white"))
    ax.add_patch(Circle(center,0.055,fc="white",ec="white"))
    ax.text(0.21,0.19,"81.1%",ha="center",va="center",fontsize=12,fontweight="bold",color=PURPLE)
    ax.text(0.37,0.215,"81.1% of the absolute component decline",fontsize=9.5,
            fontweight="bold",color=NAVY)
    ax.text(0.37,0.17,"is carried by loss of between-target signal.",fontsize=8.7,color=NAVY)
    ax.text(0.37,0.125,r"$C_W$ remains positive under Family-OOD (+0.1580).",
            fontsize=8.3,color=GREEN)

    # --- D ---
    ax = fig.add_subplot(gs[1,1])
    panel_box(ax, "Retrospective ablation trajectory on the 3,351-row universe", "D", LIGHT_GREEN)
    ax.text(0.50,0.855,
            "Matched local-physics ablations; raw values from Supplementary Fig. S2",
            ha="center",fontsize=8.8,color=NAVY)

    rows=[
        ("1","Mutation-only baseline (F1)",0.361,0.070),
        ("2","+ Static local geometry (F2)",0.389,0.068),
        ("3","+ Mutation × partner perturbation (F4)",0.431,0.113),
        ("4","Full local-physics model (F8)",0.466,0.145),
    ]
    # table region
    ax.text(0.61,0.79,"Pooled Pearson",ha="center",fontsize=8.2,fontweight="bold",color=PURPLE)
    ax.text(0.84,0.79,"Within-target\nMacro Spearman ρ",ha="center",fontsize=8.0,fontweight="bold",color=GREEN)
    for i,(idx,lab,pv,rv) in enumerate(rows):
        yy=0.73-i*0.085
        rounded(ax,(0.05,yy-0.032),(0.90,0.063),fc=("#F4FAF6" if i==3 else "#FBFBFC"),
                ec=("#B8DCC6" if i==3 else "#E1E4E9"),r=0.012)
        ax.add_patch(Circle((0.085,yy),0.022,fc=(GREEN if i==3 else "#A7ADB8"),ec="none"))
        ax.text(0.085,yy,idx,ha="center",va="center",fontsize=7.5,color="white",fontweight="bold")
        ax.text(0.125,yy,lab,ha="left",va="center",fontsize=7.8,color=NAVY,
                fontweight=("bold" if i==3 else "normal"))
        ax.text(0.61,yy,f"{pv:.3f}",ha="center",va="center",fontsize=9,fontweight="bold",color=PURPLE)
        ax.text(0.84,yy,f"{rv:.3f}",ha="center",va="center",fontsize=9,fontweight="bold",color=GREEN)

    # inset trajectory — larger and readable
    ins=ax.inset_axes([0.13,0.08,0.74,0.30])
    x=np.array([0.361,0.389,0.431,0.466])
    y=np.array([0.070,0.068,0.113,0.145])
    ins.set_facecolor("#FBFCFE")
    ins.grid(True,color=GRID,lw=0.7,alpha=0.8)
    ins.plot(x,y,color=BLUE,lw=2)
    point_cols=[GRAY,ORANGE,BLUE,GREEN]
    for j,(xx,yy,c) in enumerate(zip(x,y,point_cols),1):
        ins.scatter(xx,yy,s=45,c=c,edgecolors="white",linewidths=0.8,zorder=3)
        ins.text(xx,yy+0.006,str(j),ha="center",va="center",fontsize=7,color="white",
                 fontweight="bold",bbox=dict(boxstyle="circle,pad=0.15",fc=c,ec="white",lw=0.6))
    ins.set_xlim(0.35,0.475); ins.set_ylim(0.055,0.158)
    ins.set_xlabel(r"Pooled Pearson $r_{\mathrm{pooled}}$",fontsize=7.8,color=PURPLE,fontweight="bold")
    ins.set_ylabel("Macro Spearman ρ",fontsize=7.8,color=GREEN,fontweight="bold")
    ins.tick_params(labelsize=6.8)
    ins.annotate("F2: pooled ↑,\nwithin-target ≈ unchanged",
                 xy=(0.389,0.068),xytext=(0.352,0.095),fontsize=6.5,color="#9A5C1D",
                 arrowprops=dict(arrowstyle="->",color="#9A5C1D",lw=0.8),
                 bbox=dict(boxstyle="round,pad=0.18",fc=LIGHT_ORANGE,ec="#E6B46E",lw=0.7))
    ins.annotate("F4/F8: both ↑",xy=(0.466,0.145),xytext=(0.417,0.151),
                 fontsize=6.6,color=GREEN,
                 arrowprops=dict(arrowstyle="->",color=GREEN,lw=0.8),
                 bbox=dict(boxstyle="round,pad=0.18",fc=LIGHT_GREEN,ec="#91C8A5",lw=0.7))

    fig.suptitle("Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors",
                 y=0.995, fontsize=14, fontweight="bold", color=NAVY)
    fig.subplots_adjust(left=0.018,right=0.988,top=0.958,bottom=0.025)
    save_all(fig,"Figure1_large_manuscript")

# ============================================================================
# FIGURE 2
# ============================================================================
def figure2():
    split=["Random\n(Permissive)","Complex\n(Interface)","seq40\n(Sequence)","Family-OOD\n(Strict Extrap.)"]
    pooled=np.array([0.6844,0.1931,0.1931,0.1508])
    cb=np.array([0.4275,0.0233,0.0233,-0.0072])
    cw=np.array([0.2591,0.1685,0.1685,0.1580])

    fig=plt.figure(figsize=(11.5,6.2))
    gs=fig.add_gridspec(2,2,height_ratios=[1,0.82],wspace=0.26,hspace=0.42)
    x=np.arange(4)

    ax=fig.add_subplot(gs[0,0])
    ax.plot(x,pooled,marker="o",color=NAVY,lw=2)
    for xi,v in zip(x,pooled):
        ax.text(xi,v+0.025,f"{v:.4f}",ha="center",fontsize=8,fontweight="bold")
    ax.axvspan(0.75,2.25,color="#EEF0F3",alpha=0.9,zorder=-2)
    ax.text(1.5,0.40,"identical partition\nin this dataset",ha="center",fontsize=8,color=GRAY)
    ax.set_xticks(x,split)
    ax.set_ylim(0,0.78); ax.set_ylabel(r"Pooled Pearson ($r_{\mathrm{pooled}}$)")
    ax.set_title("A   Pooled correlation by split regime",loc="left",fontweight="bold")
    ax.grid(axis="y",color=GRID,ls="--",lw=0.7)

    ax=fig.add_subplot(gs[0,1])
    ax.plot(x,cb,marker="s",color=PURPLE,label=r"$C_B$ (between-target)")
    ax.plot(x,cw,marker="D",color=GREEN,label=r"$C_W$ (within-target)")
    for xi,v in zip(x,cb):
        ax.text(xi,v+(0.025 if v>=0 else -0.045),f"{v:+.4f}",ha="center",fontsize=7.2,color=PURPLE,fontweight="bold")
    for xi,v in zip(x,cw):
        ax.text(xi,v+0.025,f"{v:+.4f}",ha="center",fontsize=7.2,color=GREEN,fontweight="bold")
    ax.axhline(0,color=GRAY,lw=0.8,ls=":")
    ax.set_xticks(x,split); ax.set_ylim(-0.07,0.50)
    ax.set_ylabel("Covariance contribution")
    ax.set_title(r"B   Exact decomposition: $r_{\mathrm{pooled}}=C_B+C_W$",loc="left",fontweight="bold")
    ax.legend(frameon=False,loc="upper right")
    ax.grid(axis="y",color=GRID,ls="--",lw=0.7)

    ax=fig.add_subplot(gs[1,:])
    vals=[0.4347,0.1011]; labs=[r"Between-target $|\Delta C_B|$",r"Within-target $|\Delta C_W|$"]
    bars=ax.barh([1,0],vals,color=[PURPLE,GREEN],height=0.5)
    ax.set_yticks([1,0],labs)
    ax.set_xlim(0,0.60)
    ax.set_xlabel("Absolute change in covariance component (Random → Family-OOD)")
    ax.set_title(r"C   Attribution of the absolute decline: $|\Delta C_B|/(|\Delta C_B|+|\Delta C_W|)$",
                 loc="left",fontweight="bold")
    ax.text(vals[0]+0.012,1,"0.4347  (81.1%)",va="center",color=PURPLE,fontweight="bold")
    ax.text(vals[1]+0.012,0,"0.1011  (18.9%)",va="center",color=GREEN,fontweight="bold")
    ax.grid(axis="x",color=GRID,ls="--",lw=0.7)

    fig.suptitle("Family-out-of-distribution evaluation disproportionately removes the between-target component",
                 fontsize=13,fontweight="bold",y=0.995)
    fig.subplots_adjust(left=0.09,right=0.98,top=0.91,bottom=0.10)
    save_all(fig,"Figure2")

# ============================================================================
# FIGURE 3
# ============================================================================
FIG3_CB_COUNTS=np.array([7,3,3,10,3,1,5,3,5,14,8,15,12,5,1,6,4,4])
FIG3_CB_CENTERS=np.array([0.0007,0.0036,0.0050,0.0066,0.0081,0.0109,0.0124,0.0138,0.0152,0.0166,0.0180,0.0195,0.0210,0.0224,0.0248,0.0262,0.0277,0.0292])
FIG3_CW_COUNTS=np.array([7,2,5,1,4,25,14,8,14,2,1,1,5,3,6,4,3,4])
FIG3_CW_CENTERS=np.array([0.007,0.021,0.028,0.042,0.049,0.056,0.063,0.070,0.077,0.084,0.091,0.098,0.105,0.112,0.119,0.126,0.133,0.140])

def figure3():
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.7),gridspec_kw={"width_ratios":[1,1,1.02]})
    ax=axes[0]
    ax.bar(FIG3_CB_CENTERS,FIG3_CB_COUNTS,width=0.00125,color=PURPLE,alpha=0.88)
    ax.axvline(0.0183,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\n= 0.0183",xy=(0.0183,20),xytext=(0.023,22.5),
                arrowprops=dict(arrowstyle="->",color=NAVY),
                bbox=dict(boxstyle="round,pad=.25",fc=LIGHT_PURPLE,ec=PURPLE,lw=0.9),
                fontsize=8.5,fontweight="bold")
    ax.set_title(r"A   Between-target $C_B$ range",loc="left",fontweight="bold")
    ax.set_xlabel(r"Range of $C_B$ across 5 weighting schemes")
    ax.set_ylabel("Model configurations (N = 109)")
    ax.set_xlim(-0.001,0.031); ax.set_ylim(0,28); ax.grid(axis="y",color=GRID,ls="--",lw=.7)

    ax=axes[1]
    ax.bar(FIG3_CW_CENTERS,FIG3_CW_COUNTS,width=0.0055,color=GREEN,alpha=0.9)
    ax.axvline(0.0624,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\n= 0.0624",xy=(0.0624,20),xytext=(0.087,22.5),
                arrowprops=dict(arrowstyle="->",color=NAVY),
                bbox=dict(boxstyle="round,pad=.25",fc=LIGHT_GREEN,ec=GREEN,lw=0.9),
                fontsize=8.5,fontweight="bold")
    ax.set_title(r"B   Within-target $C_W$ range",loc="left",fontweight="bold")
    ax.set_xlabel(r"Range of $C_W$ across 5 weighting schemes")
    ax.set_ylabel("Model configurations (N = 109)")
    ax.set_xlim(0,0.143); ax.set_ylim(0,28); ax.grid(axis="y",color=GRID,ls="--",lw=.7)

    ax=axes[2]
    cats=[r"$|C|>0$"+"\n(No threshold)",r"$|C|>0.001$"+"\n(Near-zero filter)",r"$|C|>0.005$"+"\n(Moderate filter)"]
    bx=np.arange(3); w=.34
    cb=[35,13,0]; cw=[0,0,0]
    ax.bar(bx-w/2,cb,width=w,color=PURPLE,label=r"$C_B$ (between-target)")
    ax.bar(bx+w/2,cw,width=w,color=GREEN,label=r"$C_W$ (within-target)")
    for xx,v in zip(bx-w/2,cb):
        ax.text(xx,v+1,f"{v}/109",ha="center",color=PURPLE,fontweight="bold")
    for xx,v in zip(bx+w/2,cw):
        ax.text(xx,1 if v==0 else v+1,f"{v}/109",ha="center",color=GREEN,fontweight="bold")
    ax.set_xticks(bx,cats); ax.set_ylim(0,40)
    ax.set_ylabel("Configurations reversing sign")
    ax.set_title("C   Directional stability across weighting schemes",loc="left",fontweight="bold")
    ax.legend(frameon=False,loc="upper right")
    ax.text(1.05,24.5,
            r"$C_W$ never reverses sign"+"\n"+r"$C_B$ stable once $|C_B|>0.005$",
            color="#136C49",fontsize=8.8,fontweight="bold",
            bbox=dict(boxstyle="round,pad=.3",fc=LIGHT_GREEN,ec=GREEN,lw=1))
    ax.grid(axis="y",color=GRID,ls="--",lw=.7)

    fig.suptitle("Target weighting changes the interpretation of covariance contributions",
                 fontsize=13,fontweight="bold",y=0.995)
    fig.subplots_adjust(left=.065,right=.985,top=.88,bottom=.17,wspace=.28)
    save_all(fig,"Figure3")

# ============================================================================
# FIGURE 4
# ============================================================================
def figure4():
    K=10
    n=np.linspace(1,140,800)
    D=np.where(n<=K,0,K*(n-K)/(n*(n+1)))
    nstar=K+np.sqrt(K*K+K)
    d20=K*(20-K)/(20*21)
    d50=K*(50-K)/(50*51)

    fig,ax=plt.subplots(figsize=(11.8,6.5))
    ax.set_facecolor("#F8FAFC")
    ax.axvspan(1,K,color=LIGHT_CORAL,alpha=0.95,label=r"Degenerate regime ($n\leq K$: $D=0$)")
    ax.plot(n,D,color=NAVY,lw=3,label=r"$D(n,K)=\frac{K(n-K)}{n(n+1)}$ for $n>K$")
    ax.axvline(K,color=CORAL,ls="--",lw=1.8)
    ax.axvline(nstar,color=GREEN,ls=":",lw=1.6)

    ax.scatter([20,50],[d20,d50],s=90,c=[CORAL,CORAL],edgecolors="white",zorder=4)
    ax.annotate(r"Global peak: $n^*=K+\sqrt{K^2+K}\approx20.5$"+"\n"+r"$D(20)=0.2381$  (maximum discriminability)",
                xy=(20,d20),xytext=(33,0.251),
                arrowprops=dict(arrowstyle="->",color=GREEN,lw=1.3),
                bbox=dict(boxstyle="round,pad=.35",fc=LIGHT_GREEN,ec=GREEN,lw=1.1),
                color="#146843",fontsize=9.5,fontweight="bold")
    ax.annotate(r"$D(50)=0.1569$"+"\nLarger candidate libraries\nshrink the discriminable interval",
                xy=(50,d50),xytext=(62,0.175),
                arrowprops=dict(arrowstyle="->",color=NAVY,lw=1.1),
                bbox=dict(boxstyle="round,pad=.3",fc="white",ec=GRAY,lw=1),
                color=NAVY,fontsize=9,fontweight="bold")
    ax.annotate(r"Degenerate regime ($n\leq K$)"+"\n"+r"$D(n,K)=0$ (every model selects all candidates)"+"\n53 / 130 targets (40.8% of $n\geq5$) carry zero info",
                xy=(7,0.004),xytext=(23,0.032),
                arrowprops=dict(arrowstyle="->",color=CORAL,lw=1.2),
                bbox=dict(boxstyle="round,pad=.35",fc="#FFF6F3",ec=CORAL,lw=1.1),
                color="#8F3426",fontsize=8.8,fontweight="bold")

    # inset moved to a clear, non-overlapping zone
    ins=ax.inset_axes([0.66,0.47,0.30,0.30])
    vals=[53,77]
    ins.bar([0,1],vals,color=[CORAL,GREEN],width=.48)
    ins.set_ylim(0,105); ins.set_yticks([0,25,50,75,100])
    ins.set_xticks([0,1],[r"Degenerate"+"\n"+r"($n\leq10$)",r"Informative"+"\n"+r"($n>10$)"])
    ins.set_ylabel("Targets",fontsize=8)
    ins.set_title(r"SKEMPI v2.0 target geometry ($N=130$)",fontsize=8.8,fontweight="bold")
    ins.tick_params(labelsize=7)
    ins.grid(axis="y",color=GRID,ls="--",lw=.6)
    ins.text(0,57,"53 targets\n(40.8%)",ha="center",fontsize=7.5,fontweight="bold")
    ins.text(1,82,"77 targets\n(59.2%)",ha="center",fontsize=7.5,fontweight="bold")

    ax.set_xlim(0,142); ax.set_ylim(-0.012,0.275)
    ax.set_xlabel(r"Candidate-set size $n$  (fixed shortlist quota $K=10$)",fontsize=11,fontweight="bold")
    ax.set_ylabel(r"Discriminable interval $D(n,K)=E_{\max}-E_{\mathrm{rand}}$",fontsize=11,fontweight="bold")
    ax.set_title(r"Fixed-$K$ Shortlist Evaluation is Constrained by Candidate-Set Geometry",
                 fontsize=15,fontweight="bold",pad=14)
    ax.legend(loc="upper right",frameon=True)
    ax.grid(color=GRID,ls="--",lw=.7,alpha=.85)
    fig.subplots_adjust(left=.10,right=.98,top=.89,bottom=.13)
    save_all(fig,"Figure4")

# ============================================================================
# FIGURE 5
# ============================================================================
def figure5():
    pts={
        "F1 Mutation-only baseline":(0.000,0.000),
        "F2 + Static local geometry":(0.028,-0.002),
        "F3 + Partner static context":(0.043,0.020),
        "F4 + Mutation × partner perturbation":(0.070,0.044),
        "F5 + Intra-chain perturbation":(-0.005,-0.026),
        "F8 Full local-physics model":(0.105,0.076),
    }
    fig,ax=plt.subplots(figsize=(8.2,6.0))
    ax.axvspan(-0.03,0.12,color="#F5FAF6",zorder=-3)
    ax.axhspan(-0.035,0,color=LIGHT_CORAL,alpha=.65,zorder=-2)
    ax.axvline(0,color=GRAY,ls="--",lw=.9)
    ax.axhline(0,color=GRAY,ls="--",lw=.9)
    ax.grid(color=GRID,ls="--",lw=.65)

    # trajectory: F1 -> F2 -> F3 -> F4 -> F8
    order=["F1 Mutation-only baseline","F2 + Static local geometry","F3 + Partner static context",
           "F4 + Mutation × partner perturbation","F8 Full local-physics model"]
    xx=[pts[k][0] for k in order]; yy=[pts[k][1] for k in order]
    ax.plot(xx,yy,color=GREEN,lw=1.8,alpha=.9)

    for name,(x,y) in pts.items():
        good=(x>0 and y>0)
        c=GREEN if good else CORAL
        ax.scatter(x,y,s=58,c=c,edgecolors="white",linewidths=.8,zorder=4)
        short=name.split(" ",1)[0]
        dx=0.004 if x>=0 else -0.002
        dy=0.005 if y>=0 else -0.012
        label=name.replace("F1 ","F1: ").replace("F2 ","F2: ").replace("F3 ","F3: ").replace("F4 ","F4: ").replace("F5 ","F5: ").replace("F8 ","F8: ")
        ax.annotate(label,xy=(x,y),xytext=(x+dx,y+dy),fontsize=7.5,color=c,fontweight="bold",
                    bbox=dict(boxstyle="round,pad=.2",fc="white",ec=c,lw=.8),
                    arrowprops=dict(arrowstyle="-",color=c,lw=.7))

    ax.text(0.033,0.104,"Synergistic improvement regime\nboth pooled score and within-target ranking improve",
            color="#176C48",fontsize=8.5,fontweight="bold",
            bbox=dict(boxstyle="round,pad=.3",fc=LIGHT_GREEN,ec=GREEN,lw=1))
    ax.text(-0.027,-0.031,"Pooled-score-only / negative regime",
            color="#93372A",fontsize=8,fontweight="bold",
            bbox=dict(boxstyle="round,pad=.25",fc="#FFF6F3",ec=CORAL,lw=1))

    ax.set_xlim(-0.03,0.125); ax.set_ylim(-0.036,0.125)
    ax.set_xlabel(r"Δ Pooled Pearson ($r_{\mathrm{pooled}}$) vs. mutation-only baseline",fontweight="bold")
    ax.set_ylabel(r"Δ Macro within-target Spearman ($\rho_{\mathrm{macro}}$)",fontweight="bold")
    ax.set_title("Capability decomposition separates pooled-score gains from within-target gains",
                 fontsize=12.5,fontweight="bold")
    fig.subplots_adjust(left=.13,right=.97,top=.90,bottom=.13)
    save_all(fig,"Figure5")

# ============================================================================
# FIGURE 6
# ============================================================================
def figure6():
    methods=["Flex-ddG","DDGPred","FoldX\n(ddG_FoldX)","ESM-IF","ProteinMPNN"]
    pooled=np.array([0.2028,0.2406,0.2089,-0.2250,0.1509])
    macro=np.array([0.1092,0.1287,0.1088,0.0213,0.0320])
    cb=np.array([0.0478,0.0715,0.0657,-0.1149,0.0699])
    cw=np.array([0.1550,0.1690,0.1433,-0.1101,0.0810])
    data=np.column_stack([pooled,macro,cb,cw])

    fig=plt.figure(figsize=(10.8,4.8))
    gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.12],wspace=.30)

    ax=fig.add_subplot(gs[0,0])
    cols=[BLUE,GREEN,ORANGE,CORAL,PURPLE]
    for m,x,y,c in zip(methods,pooled,macro,cols):
        ax.scatter(x,y,s=70,c=c,edgecolors="white",linewidths=.8,zorder=3)
        ax.annotate(f"{m}\n({x:+.3f}, {y:.3f})",xy=(x,y),xytext=(5,5),
                    textcoords="offset points",fontsize=7,color=c,fontweight="bold")
    ax.grid(color=GRID,ls="--",lw=.7)
    ax.set_xlim(-.30,.31); ax.set_ylim(0,.16)
    ax.set_xlabel(r"Pooled Pearson ($r_{\mathrm{pooled}}$)",fontweight="bold")
    ax.set_ylabel(r"Macro within-target Spearman ($\rho_{\mathrm{macro}}$)",fontweight="bold")
    ax.set_title("A   Pooled vs. within-target performance",loc="left",fontweight="bold")
    ax.text(-.285,.151,"Matched unambiguous mapping\n1,781 mutations; 67 targets (n≥5)\nrank agreement ρ = +0.90",
            fontsize=7.5,va="top",
            bbox=dict(boxstyle="round,pad=.25",fc="#F4F6F9",ec=GRID,lw=.8))

    ax=fig.add_subplot(gs[0,1])
    z=(data-data.mean(axis=0))/data.std(axis=0,ddof=0)
    im=ax.imshow(z,cmap="RdBu_r",aspect="auto",vmin=-2,vmax=2)
    ax.set_yticks(np.arange(len(methods)),methods)
    ax.set_xticks(np.arange(4),["Pooled\nPearson","Macro\nSpearman",r"$C_B$"+"\n(between)",r"$C_W$"+"\n(within)"])
    ax.set_title("B   Method-level capability profile (raw values printed)",loc="left",fontweight="bold")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val=data[i,j]
            ax.text(j,i,f"{val:+.3f}" if j>=2 or val<0 else f"{val:.3f}",
                    ha="center",va="center",fontsize=7.5,fontweight="bold",
                    color=("white" if abs(z[i,j])>1.15 else NAVY))
    cbar=fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    cbar.set_label("Within-column standardized z-score",fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.suptitle("External methods define the boundary of generalisation",
                 fontsize=13,fontweight="bold",y=.99)
    fig.subplots_adjust(left=.08,right=.96,top=.87,bottom=.16)
    save_all(fig,"Figure6")

# ============================================================================
# SUPPLEMENTARY FIGURE S1
# ============================================================================
def figureS1():
    labels=[
        "V6.2.1_donly  (SHA: 8d003e00)",
        "V6.2_split4 [Controlled Split]  (SHA: 8d003e00)",
        "V7.1_lsr  (SHA: 34b45937)",
        "V7.0_losses  (SHA: 34b45937)",
        "V7.0_ablation  (SHA: 34b45937)",
        "V7.1_z  (SHA: 34b45937)",
        "V7.2R1_full  (SHA: 2f5a5c2e)",
        "V7.2R2_trackB  (SHA: e81801ff)",
        "V7.2R3_noise  (SHA: 17487417)",
        "V7.2R3_ppi  (SHA: 17487417)",
        "V7.3_physics [Primary]  (SHA: f8c33c7e)",
        "V7.2_nested  (SHA: c08fe062)",
    ]
    rows=np.array([7523,7523,4473,4473,4473,4473,4233,4219,3452,3452,3351,1807])
    targ=np.array([172,172,130,130,130,130,116,115,87,87,102,12])
    y=np.arange(len(labels))
    fig,ax=plt.subplots(figsize=(10.4,6.0))
    ax.barh(y,rows,color=BLUE,alpha=.92,height=.55,label="Mutation rows (N)")
    # targets shown as orange scaled bar for visual comparison, same original convention
    ax.barh(y,targ*5,color=ORANGE,alpha=.95,height=.18,label="Targets with n ≥ 5 (×5 display scale)")
    ax.set_yticks(y,labels); ax.invert_yaxis()
    ax.set_xlabel("Count (mutation rows / qualifying targets)")
    ax.set_title("Figure S1 | Matched-evaluation strata grouped by SHA-256 row-set fingerprints",
                 fontweight="bold")
    ax.grid(axis="x",color=GRID,ls="--",lw=.65)
    for yi,v in zip(y,rows):
        ax.text(v+80,yi,f"{v:,} rows",va="center",fontsize=7,color=NAVY)
    for yi,v in zip(y,targ):
        ax.text(v*5+25,yi+0.16,f"{v} targets (n≥5)",va="center",fontsize=6.5,color="#9A5C1D")
    ax.legend(frameon=False,loc="lower right")
    fig.subplots_adjust(left=.31,right=.97,top=.90,bottom=.11)
    save_all(fig,"FigureS1")

# ============================================================================
# SUPPLEMENTARY FIGURE S2
# ============================================================================
def figureS2():
    rows=[
        "Full local-physics model (F8)",
        "Partner + intra-chain",
        "Mutation × partner perturbation (F4)",
        "Partner + intra + burial",
        "Partner static context (F3)",
        "Physics features only",
        "Mutation-only baseline (F1)",
        "Static local geometry (F2)",
        "Intra-chain perturbation (F5)",
        "Base model",
    ]
    raw=np.array([
        [0.466,0.145,0.142,0.325],
        [0.416,0.129,0.099,0.317],
        [0.431,0.113,0.119,0.312],
        [0.458,0.111,0.140,0.318],
        [0.404,0.090,0.111,0.293],
        [0.361,0.074,0.166,0.195],
        [0.361,0.070,0.063,0.298],
        [0.389,0.068,0.090,0.299],
        [0.356,0.044,0.062,0.294],
        [0.360,0.042,0.094,0.266],
    ])
    z=(raw-raw.mean(axis=0))/raw.std(axis=0,ddof=0)
    fig,ax=plt.subplots(figsize=(8.4,5.7))
    im=ax.imshow(z,cmap="RdBu_r",aspect="auto",vmin=-2,vmax=2)
    ax.set_yticks(np.arange(len(rows)),rows)
    ax.set_xticks(np.arange(4),["Pooled\nPearson", "Macro\nSpearman", r"$C_B$"+"\n(between-target)", r"$C_W$"+"\n(within-target)"])
    ax.set_title("Figure S2 | Capability profiles of matched, provenance-clean ablation variants (seed 42)",
                 fontweight="bold",fontsize=11)
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            ax.text(j,i,f"{raw[i,j]:.3f}",ha="center",va="center",fontsize=7.5,
                    color=("white" if abs(z[i,j])>1.2 else NAVY),fontweight="bold")
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025)
    cb.set_label("Within-column standardized z-score",fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.subplots_adjust(left=.30,right=.92,top=.88,bottom=.14)
    save_all(fig,"FigureS2")

# ============================================================================
# SUPPLEMENTARY FIGURE S3
# ============================================================================
def figureS3():
    labels=[
        "V7.3_physics  (f8c33c7e / LA)*",
        "V6.2_split4  (8d003e00 / LA)*",
        "V6.2_split4  (8d003e00 / LC)*",
        "V7.2_nested  (c08fe062 / LA)*",
        "V7.2R1_full  (2f5a5c2e / LA)*",
        "V7.2R2_trackB  (e81801ff / LA)*",
        "V7.2R3_ppi  (17487417 / LC)*",
        "V7.2R3_ppi  (17487417 / LA)*",
        "V7.1_z  (34b45937 / LC)*",
        "V7.1_z  (34b45937 / LA)*",
    ]
    vals=np.array([0.82,0.76,0.69,0.53,0.43,0.34,0.08,-0.02,-0.07,-0.26])
    y=np.arange(len(labels))
    fig,ax=plt.subplots(figsize=(9.0,5.7))
    colors=np.where(vals>=0,BLUE,CORAL)
    ax.barh(y,vals,color=colors,height=.58)
    ax.axvline(0,color=NAVY,lw=1)
    ax.set_yticks(y,labels); ax.invert_yaxis()
    ax.set_xlim(-.32,.88)
    ax.set_xlabel(r"Rank correlation ρ(Pooled Pearson, Macro Within-Target Spearman)")
    ax.set_title("Figure S3 | Pool sensitivity of cross-metric rank agreement across matched strata",
                 fontweight="bold")
    ax.grid(axis="x",color=GRID,ls="--",lw=.65)
    for yi,v in zip(y,vals):
        ax.text(v+(0.018 if v>=0 else -0.018),yi,f"ρ={v:+.2f}",
                va="center",ha=("left" if v>=0 else "right"),fontsize=7.5,
                color=(BLUE if v>=0 else CORAL),fontweight="bold")
    fig.subplots_adjust(left=.34,right=.97,top=.90,bottom=.13)
    save_all(fig,"FigureS3")

FIG_FUNCS={
    "1":figure1, "2":figure2, "3":figure3, "4":figure4, "5":figure5, "6":figure6,
    "S1":figureS1, "S2":figureS2, "S3":figureS3,
}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--figure",default="all",
                   help="1,2,3,4,5,6,S1,S2,S3 or all")
    args=p.parse_args()
    if args.figure.lower()=="all":
        for f in FIG_FUNCS.values(): f()
    else:
        key=args.figure.upper()
        if key not in FIG_FUNCS:
            raise SystemExit(f"Unknown figure: {args.figure}")
        FIG_FUNCS[key]()

if __name__=="__main__":
    main()
