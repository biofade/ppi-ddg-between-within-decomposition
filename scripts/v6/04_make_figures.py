# -*- coding: utf-8 -*-
"""scripts/v6/04_make_figures.py -- V6.0 PHASE 4: 诊断图表.

Figure 1: The Illusion Scatter —— Overall Pearson vs Macro per-PDB Spearman
Figure 2: Efficiency Gap —— within-complex Top-10 真实平均 ΔΔG vs Random/Oracle
Figure 3: Hit@10 vs Base Rate

**如实呈现**：若外部模型表现良好（局部 Spearman 高），图中如实显示，
不为了支持假设而调整任何阈值或筛选。
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT = C.ROOT / "experiments" / "v6_sota_audit"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 8, "axes.grid": True,
                     "grid.alpha": 0.3})
made: list[str] = []


def save(fig, name: str) -> None:
    p = FIG / name
    fig.tight_layout()
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    made.append(name)
    print(f"写出 {p}")


def main() -> int:
    m = pd.read_csv(OUT / "sota_audit_metrics.csv", low_memory=False)
    # 每个模型一行（V5 baseline 取 3 seed 均值）
    v5 = m[m["dataset"].str.startswith("V5/")]
    v5row = {"dataset": "V5 (ours)", "model": "V5 baseline",
             "overall_pearson": v5["overall_pearson"].mean(),
             "macro_pdb_spearman": v5["macro_pdb_spearman"].mean(),
             "hit10_lt0": v5["hit10_lt0"].mean(),
             "base_rate_lt0": v5["base_rate_lt0"].mean(),
             "efficiency_at_10": v5["efficiency_at_10"].mean(),
             "n_group": int(v5["n_group"].max()),
             "utility_model": v5["utility_model"].mean(),
             "utility_random": v5["utility_random"].mean(),
             "utility_oracle": v5["utility_oracle"].mean()}
    ext = m[~m["dataset"].str.startswith("V5/")].copy()
    plot = pd.concat([ext, pd.DataFrame([v5row])], ignore_index=True)
    plot["label"] = plot["model"] + "\n(" + plot["dataset"].str.split("/").str[0] + ")"

    # ---------------- Figure 1 ----------------
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    is_v5 = plot["model"] == "V5 baseline"
    is_cath = plot["dataset"].str.contains("CATH")
    ax.scatter(plot.loc[is_cath, "overall_pearson"],
               plot.loc[is_cath, "macro_pdb_spearman"],
               s=90, marker="o", color="C0", alpha=0.8,
               label="CATH-ddG official held-out CATH (n=21 groups)")
    ax.scatter(plot.loc[~is_cath & ~is_v5, "overall_pearson"],
               plot.loc[~is_cath & ~is_v5, "macro_pdb_spearman"],
               s=90, marker="s", color="C2", alpha=0.8,
               label="DSSA-PPI RBM case study (n=1 group)")
    ax.scatter(plot.loc[is_v5, "overall_pearson"],
               plot.loc[is_v5, "macro_pdb_spearman"],
               s=230, marker="*", color="red", zorder=6,
               label="V5 baseline (ours, n=61 groups)")
    for _, r in plot.iterrows():
        ax.annotate(str(r["model"])[:16], (r["overall_pearson"],
                                          r["macro_pdb_spearman"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, color="k", lw=0.9, ls="--")
    ax.axhline(0.15, color="gray", lw=0.8, ls=":", label="hypothesis threshold (0.15)")
    ax.axhline(0.30, color="green", lw=0.8, ls=":", label="refutation threshold (0.30)")
    ax.axvline(0.40, color="gray", lw=0.8, ls=":")
    ax.set_xlabel("Overall Pearson  (global) ")
    ax.set_ylabel("Macro per-PDB Spearman  (within-target ranking)")
    ax.set_title("Fig 1 — The Illusion Scatter\n"
                 "hypothesis predicted lower-right cluster; observed: external models sit UPPER area")
    ax.legend(fontsize=6.5, loc="upper left")
    save(fig, "fig1_illusion_scatter.png")

    # ---------------- Figure 2 ----------------
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    d = plot.sort_values("utility_model")
    xs = np.arange(len(d))
    ax.bar(xs, d["utility_model"], color=["red" if x == "V5 baseline" else "C0"
                                          for x in d["model"]], alpha=0.85,
           label="model selection (mean −ΔΔG of top-10; higher = better)")
    ax.plot(xs, d["utility_random"], "s--", color="gray", ms=4,
            label="Random choice baseline")
    ax.plot(xs, d["utility_oracle"], "^--", color="green", ms=4,
            label="Oracle (perfect)")
    ax.set_xticks(xs)
    ax.set_xticklabels(d["label"], fontsize=5.5, rotation=45, ha="right")
    ax.set_ylabel("within-target top-10 utility  (= −mean true $\\Delta\\Delta G$)")
    ax.set_title("Fig 2 — Efficiency Gap: how far each model is from Oracle vs Random\n"
                 "(utility reported in the −ΔΔG orientation: higher = better)")
    ax.legend(fontsize=6.5, loc="lower right")
    save(fig, "fig2_efficiency_gap.png")

    # ---------------- Figure 3 ----------------
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    d = plot.sort_values("hit10_lt0", ascending=False)
    xs = np.arange(len(d))
    colors = ["red" if r.hit10_lt0 < r.base_rate_lt0 else "C0" for r in d.itertuples()]
    colors = ["red" if c == "red" else ("red" if m == "V5 baseline" else "C0")
              for c, m in zip(colors, d["model"])]
    ax.bar(xs, d["hit10_lt0"], color=colors, alpha=0.85,
           label="Hit@10 (ΔΔG<0): fraction of predicted top-10 that are beneficial")
    ax.plot(xs, d["base_rate_lt0"], "k^--", ms=5,
            label="Base rate (random selection expectation)")
    ax.set_xticks(xs)
    ax.set_xticklabels(d["label"], fontsize=5.5, rotation=45, ha="right")
    ax.set_ylabel("Hit@10  /  Base rate")
    ax.set_ylim(0, 1.05)
    ax.set_title("Fig 3 — Hit@10 vs Base Rate\n"
                 "(red = Hit@10 below its own base rate, i.e. worse than random)")
    for i, r in enumerate(d.itertuples()):
        ax.annotate(f"base={r.base_rate_lt0:.2f}", (i, r.base_rate_lt0),
                    fontsize=5, color="k", ha="center", va="bottom")
    ax.legend(fontsize=6.5, loc="lower left")
    save(fig, "fig3_hit10_vs_baserate.png")

    # ---------------- Figure 4 (bonus): Efficiency@10 ----------------
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    d = plot.sort_values("efficiency_at_10", ascending=False)
    xs = np.arange(len(d))
    ax.bar(xs, d["efficiency_at_10"],
           color=["red" if x == "V5 baseline" else "C0" for x in d["model"]],
           alpha=0.85)
    ax.axhline(0, color="k", lw=0.9)
    ax.axhline(1.0, color="green", lw=0.9, ls="--", label="Oracle = 1.0")
    ax.set_xticks(xs)
    ax.set_xticklabels(d["label"], fontsize=5.5, rotation=45, ha="right")
    ax.set_ylabel("Efficiency@10 = (Random − Model) / (Random − Oracle)")
    ax.set_title("Fig 4 — Efficiency@10 across all audited models\n"
                 "(0 = no better than random; 1 = oracle-level)")
    ax.legend(fontsize=7)
    save(fig, "fig4_efficiency_at_10.png")

    print(f"\n生成 {len(made)} 张图: {made}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
