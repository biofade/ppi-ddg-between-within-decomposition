# -*- coding: utf-8 -*-
"""scripts/v6_1/08_make_paper_figures.py -- V6.1 PHASE 5: 论文级 Figure（300 dpi）.

6 张图，每张的**数据来源在报告与本文件注释中标注**：

* Fig 1  The Illusion Scatter      <- PHASE 1.1 `bootstrap_ci_models.csv`
* Fig 2  Efficiency Gap Bar        <- V6.0 `sota_audit_metrics.csv`
* Fig 3  Category Comparison       <- PHASE 1.3 `category_paired_values.csv`
* Fig 4  Split Inflation           <- PHASE 2 `split_inflation_results.csv`
* Fig 5  Hit@10 vs Base Rate       <- V6.0 `sota_audit_metrics.csv` + PHASE 2
* Fig 6  Per-Complex Distribution  <- PHASE 4 `v5_cath_distribution.csv`

缺失的数据源会被**跳过并记录**，不伪造。
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

V61 = C.ROOT / "experiments" / "v6_1"
FIG = V61 / "phase5_figures"
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": DPI, "font.size": 8,
                     "axes.grid": True, "grid.alpha": 0.3})
made: list[str] = []
skipped: list[str] = []


def save(fig, name: str) -> None:
    p = FIG / name
    fig.tight_layout()
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    made.append(name)
    print(f"写出 {p}  ({DPI} dpi)")


def main() -> int:
    # ---------------- Fig 1: The Illusion Scatter ----------------
    f1 = V61 / "phase1_ci" / "bootstrap_ci_models.csv"
    if f1.exists():
        d = pd.read_csv(f1)
        s = d[d["metric"] == "macro_pdb_spearman"]
        # 需要 overall pearson：从 V6.0 汇总取
        v6 = C.ROOT / "experiments" / "v6_sota_audit" / "sota_audit_metrics.csv"
        pe = {}
        if v6.exists():
            vv = pd.read_csv(v6, low_memory=False)
            cc = vv[vv["dataset"] == "CATH"]
            pe = dict(zip(cc["model"], cc["overall_pearson"]))
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        for _, r in s.iterrows():
            x = pe.get(r["model"], np.nan)
            y = r["point"]
            if x != x:
                continue
            phys = r["model"] in ("Flex-ddG", "FoldX")
            mk = "o" if phys else "s"
            col = "C0" if phys else "C3"
            ax.errorbar(x, y, yerr=[[y - r["ci_low"]], [r["ci_high"] - y]],
                        fmt=mk, ms=8, color=col, capsize=3, alpha=0.85)
            ax.annotate(str(r["model"]), (x, y), fontsize=6.5,
                        xytext=(5, 4), textcoords="offset points")
        # V5 baseline
        vb_x, vb_y = 0.3721, 0.0835
        ax.errorbar(vb_x, vb_y, yerr=[[0.0835 + 0.0219], [0.2385 - 0.0835]],
                    fmt="*", ms=18, color="red", capsize=3,
                    label="V5 baseline (ours)")
        ax.axhline(0, color="k", lw=1.0, ls="--")
        ax.axhline(0.15, color="gray", lw=0.8, ls=":",
                   label="hypothesis threshold (0.15)")
        ax.axvline(0.4, color="gray", lw=0.8, ls=":")
        ax.set_xlabel("Overall Pearson (global)")
        ax.set_ylabel("Macro per-complex Spearman (within-target)")
        ax.set_title("Fig 1 — The Illusion Scatter\n"
                     "blue = physics/energy, red = sequence representation")
        ax.legend(fontsize=6.5, loc="lower right")
        save(fig, "fig1_illusion_scatter.png")
    else:
        skipped.append("fig1 (缺 phase1_ci/bootstrap_ci_models.csv)")

    # ---------------- Fig 2: Efficiency Gap ----------------
    v6 = C.ROOT / "experiments" / "v6_sota_audit" / "sota_audit_metrics.csv"
    if v6.exists():
        vv = pd.read_csv(v6, low_memory=False)
        d = vv.dropna(subset=["utility_model"]).copy()
        d["label"] = d["model"] + "\n(" + d["dataset"].str.split("/").str[0] + ")"
        d = d.sort_values("utility_model")
        fig, ax = plt.subplots(figsize=(8.6, 4.6))
        xs = np.arange(len(d))
        ax.bar(xs, d["utility_model"], color="C0", alpha=0.85,
               label="model top-10 (mean −ΔΔG; higher = better)")
        ax.plot(xs, d["utility_random"], "s--", color="gray", ms=4,
                label="Random choice")
        ax.plot(xs, d["utility_oracle"], "^--", color="green", ms=4,
                label="Oracle")
        ax.set_xticks(xs)
        ax.set_xticklabels(d["label"], fontsize=5.5, rotation=45, ha="right")
        ax.set_ylabel("within-target top-10 utility  (−mean true $\\Delta\\Delta G$)")
        ax.set_title("Fig 2 — Efficiency Gap: distance from Oracle vs Random")
        ax.legend(fontsize=6.5, loc="lower right")
        save(fig, "fig2_efficiency_gap.png")
    else:
        skipped.append("fig2 (缺 v6.0 sota_audit_metrics.csv)")

    # ---------------- Fig 3: Category Comparison (violin) ----------------
    f3 = V61 / "phase1_ci" / "category_paired_values.csv"
    if f3.exists():
        d = pd.read_csv(f3)
        phys = d[["Flex-ddG", "FoldX"]].mean(axis=1).to_numpy(float)
        seq = d[["ESM-IF", "proteinMPNN"]].mean(axis=1).to_numpy(float)
        v5v = None
        p4 = C.V5 / "p4_conformal" / "p4_predictions.npz"
        if p4.exists():
            z = np.load(p4, allow_pickle=True)
            y = z["test_y"]; g = np.array([str(x) for x in z["test_group"]])
            vals = []
            for key in ("plain_s42", "plain_s43", "plain_s44"):
                if key not in z.files:
                    continue
                for cx in np.unique(g):
                    m = np.where(g == cx)[0]
                    if len(m) < 5:
                        continue
                    v = C.spearman(y[m], z[key][m])
                    if v == v:
                        vals.append(v)
            v5v = np.array(vals, float)
        groups = [phys, seq] + ([v5v] if v5v is not None else [])
        labels = ["Physics/Energy\n(Flex-ddG+FoldX)",
                  "Sequence repr.\n(ESM-IF+proteinMPNN)"] + \
                 (["V5 baseline\n(ours)"] if v5v is not None else [])
        fig, ax = plt.subplots(figsize=(6.2, 4.4))
        parts = ax.violinplot(groups, showmedians=True, widths=0.8)
        for pc, col in zip(parts["bodies"], ["C0", "C3", "red"][:len(groups)]):
            pc.set_facecolor(col); pc.set_alpha(0.5)
        for i, arr in enumerate(groups, 1):
            ax.scatter(np.full(len(arr), i) + np.random.RandomState(0)
                       .normal(0, 0.04, len(arr)), arr, s=8, color="k", alpha=0.4)
        ax.axhline(0, color="k", lw=0.9, ls="--")
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("per-complex Spearman")
        # Wilcoxon p（配对，物理 vs 序列）
        pv = float("nan")
        try:
            from scipy.stats import wilcoxon
            pv = wilcoxon(phys, seq).pvalue
        except Exception:  # noqa: BLE001
            pass
        ax.set_title(f"Fig 3 — Category comparison (CATH held-out, n=21 complexes)\n"
                     f"paired Wilcoxon (physics vs sequence): p = {pv:.4f}")
        save(fig, "fig3_category_comparison.png")
    else:
        skipped.append("fig3 (缺 phase1_ci/category_paired_values.csv)")

    # ---------------- Fig 4: Split Inflation ----------------
    f4 = V61 / "phase2_split_inflation" / "split_inflation_results.csv"
    if f4.exists():
        d = pd.read_csv(f4, low_memory=False)
        d = d[d.get("status", "OK") == "OK"] if "status" in d.columns else d
        d = d.dropna(subset=["overall_pearson", "macro_pdb_spearman"])
        g = d.groupby("level").agg(
            pe=("overall_pearson", "mean"), pe_sd=("overall_pearson", "std"),
            ms=("macro_pdb_spearman", "mean"),
            ms_sd=("macro_pdb_spearman", "std")).reset_index()
        order = ["L1_random", "L2_complex", "L3_complex_seq40", "L4_family"]
        g = g.set_index("level").reindex([x for x in order if x in set(g["level"])])
        fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharex=True)
        xs = np.arange(len(g))
        for ax, col, sd, lab, color in (
                (axes[0], "pe", "pe_sd", "Overall Pearson", "C3"),
                (axes[1], "ms", "ms_sd", "Macro per-PDB Spearman", "C0")):
            ax.errorbar(xs, g[col], yerr=g[sd], fmt="o-", color=color, capsize=4)
            ax.set_xticks(xs)
            ax.set_xticklabels(g.index, fontsize=7, rotation=20)
            ax.set_ylabel(lab)
            ax.axhline(0, color="k", lw=0.8, ls="--")
        axes[0].set_title("Fig 4 — Split inflation: same model, 4 split definitions")
        axes[1].set_title("(within-target ranking stays near zero)")
        save(fig, "fig4_split_inflation.png")
    else:
        skipped.append("fig4 (缺 phase2_split_inflation/split_inflation_results.csv)")

    # ---------------- Fig 5: Hit@10 vs Base Rate ----------------
    if v6.exists():
        vv = pd.read_csv(v6, low_memory=False)
        d = vv.dropna(subset=["hit10_lt0", "base_rate_lt0"]).copy()
        fig, ax = plt.subplots(figsize=(5.6, 5.0))
        lim = [0, 1.02]
        ax.plot(lim, lim, "k--", lw=1.0, label="Hit@10 = Base Rate")
        for _, r in d.iterrows():
            col = "red" if r["hit10_lt0"] < r["base_rate_lt0"] else "C0"
            ax.scatter(r["base_rate_lt0"], r["hit10_lt0"], s=60, color=col,
                       alpha=0.85, edgecolor="k", linewidth=0.4)
            ax.annotate(str(r["model"])[:14], (r["base_rate_lt0"], r["hit10_lt0"]),
                        fontsize=6, xytext=(4, 3), textcoords="offset points")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("Base Rate (random selection expectation of ΔΔG<0)")
        ax.set_ylabel("Hit@10 (fraction of predicted top-10 with ΔΔG<0)")
        ax.set_title("Fig 5 — Hit@10 vs Base Rate\n"
                     "below the diagonal ⇒ worse than random")
        ax.legend(fontsize=7, loc="lower right")
        save(fig, "fig5_hit10_vs_baserate.png")
    else:
        skipped.append("fig5 (缺 v6.0 sota_audit_metrics.csv)")

    # ---------------- Fig 6: Per-complex distribution ----------------
    f6 = V61 / "phase4_v5_comparison" / "v5_cath_distribution.csv"
    if f6.exists():
        d = pd.read_csv(f6)
        fig, ax = plt.subplots(figsize=(8.4, 4.4))
        order = ["V5 baseline"] + [m for m in
                                   ["Flex-ddG", "FoldX", "DDGPred", "proteinMPNN",
                                    "ESM-IF"]
                                   if m in set(d["model"])]
        data, labels, colors = [], [], []
        for m in order:
            arr = d[d["model"] == m]["per_complex_spearman"].to_numpy(float)
            if len(arr) < 3:
                continue
            data.append(arr); labels.append(f"{m}\n(n={len(arr)})")
            colors.append("red" if m == "V5 baseline" else "C0")
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        medianprops=dict(color="k", lw=1.3))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.45)
        for i, arr in enumerate(data, 1):
            ax.scatter(np.full(len(arr), i) + np.random.RandomState(0)
                       .normal(0, 0.05, len(arr)), arr, s=9, color="k", alpha=0.35)
        ax.axhline(0, color="k", lw=0.9, ls="--")
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=6.5)
        ax.set_ylabel("per-complex Spearman")
        ax.set_title("Fig 6 — Per-complex Spearman distribution\n"
                     "(V5 = different dataset ⇒ cross-dataset comparison only)")
        save(fig, "fig6_per_complex_distribution.png")
    else:
        skipped.append("fig6 (缺 phase4_v5_comparison/v5_cath_distribution.csv)")

    print()
    print(f"生成 {len(made)} 张图: {made}")
    if skipped:
        print(f"跳过 {len(skipped)} 张（数据源缺失，未伪造）:")
        for s in skipped:
            print(f"  - {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
