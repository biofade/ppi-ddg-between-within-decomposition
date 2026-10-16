# -*- coding: utf-8 -*-
"""scripts/v5/15_make_figures.py -- V5.0: 生成 V5 §98 要求的论文级图.

Figure 1  Overall Pearson vs per-structure Spearman（scaffold 区分 vs 组内排序）
Figure 2  Prediction error vs nearest-train similarity（OOD 失效点）
Figure 3  MIPA ablation —— **`NOT_AVAILABLE_RESOURCE_BLOCKED`**（映射覆盖不足）
Figure 4  Top-K enrichment curves（global vs within-complex）
Figure 5  Conformal calibration（nominal vs empirical coverage）
Figure 6  Split strictness vs performance
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

FIG = C.V5 / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 8,
                     "axes.grid": True, "grid.alpha": 0.3})
made: list[str] = []


def save(fig, name: str) -> None:
    p = FIG / name
    fig.tight_layout()
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    made.append(name)
    print(f"写出 {p}")


def main() -> int:
    # ---------------- Figure 1 ----------------
    res = pd.read_csv(C.V5 / "all_experiments.csv", low_memory=False)
    res = res[res["status"].astype(str).str.upper() == "OK"]
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    fams = res["family"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(fams)))
    for c, fam in zip(colors, fams):
        g = res[res["family"] == fam]
        ax.scatter(pd.to_numeric(g["overall_pearson"], errors="coerce"),
                   pd.to_numeric(g["macro_pdb_spearman"], errors="coerce"),
                   s=18, alpha=0.75, color=c, label=str(fam))
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(0, color="k", lw=0.8, ls="--")
    b = res[res["variant"] == "C+DGaux_fwd+rev"]
    if len(b):
        ax.scatter([pd.to_numeric(b["overall_pearson"], errors="coerce").mean()],
                   [pd.to_numeric(b["macro_pdb_spearman"], errors="coerce").mean()],
                   marker="*", s=170, color="red", zorder=5,
                   label="V4 baseline (3-seed mean)")
    ax.set_xlabel("Overall Pearson (cross-complex, scaffold discrimination)")
    ax.set_ylabel("Macro per-structure Spearman (within-complex ranking)")
    ax.set_title("Fig 1 — Overall vs within-complex ability\n"
                 "high overall correlation does NOT imply within-complex ranking")
    ax.legend(fontsize=6, loc="lower right")
    save(fig, "fig1_overall_vs_per_structure.png")

    # ---------------- Figure 2 ----------------
    z = np.load(C.V5 / "p0_benchmark" / "baseline_repro_predictions.npz",
                allow_pickle=True)
    d = C.build_pair_features(verbose=False)
    y, groups = d["y"], d["groups"]
    te_i = np.where(d["mask_test"])[0]
    # nearest-train similarity：用 pooled 受体+配体嵌入的余弦（train 为参照）
    Za = d["Z_all"]
    tr_rows = d["split_idx"]["train"]
    tr_idx = np.unique(np.concatenate([d["fwd_v22_wt"][tr_rows],
                                       d["fwd_v22_mut"][tr_rows]]))
    tr_idx = tr_idx[tr_idx < len(Za)]
    T = Za[tr_idx]
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-9)
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for key, lab in (("C+DGaux_fwd+rev_s42", "V4 baseline (seed42)"),):
        if key not in z.files:
            continue
        p = z[key]
        Q = Za[np.concatenate([d["fwd_v22_wt"][te_i],
                               d["fwd_v22_mut"][te_i]])]
        Q = Q[:len(te_i)]
        # 逐 test pair 取 WT 行嵌入算相似度
        q = Za[d["fwd_v22_wt"][te_i]]
        qn = q / (np.linalg.norm(q, axis=1, keepdims=True) + 1e-9)
        sim = (qn @ Tn.T).max(axis=1)
        err = np.abs(y[te_i] - p)
        order = np.argsort(sim)
        nb = max(8, len(sim) // 8)
        xs, ys = [], []
        for i in range(0, len(order) - nb + 1, nb):
            sl = order[i:i + nb]
            xs.append(sim[sl].mean())
            ys.append(err[sl].mean())
        ax.plot(xs, ys, "o-", label=lab, ms=4)
        qs = pd.qcut(pd.Series(sim), 4, labels=["Q1 (most OOD)", "Q2", "Q3",
                                                "Q4 (most similar)"])
        for qq in qs.cat.categories:
            m = (qs == qq).to_numpy()
            print(f"    Fig2 {lab} {qq}: n={int(m.sum())}, "
                  f"mean|err|={err[m].mean():.3f}, sim={sim[m].mean():.3f}")
    ax.axhline(np.mean(np.abs(y[te_i] - z["C+DGaux_fwd+rev_s42"])), color="gray",
               ls=":", lw=1)
    ax.set_xlabel("Nearest-train WT embedding cosine similarity")
    ax.set_ylabel("Mean |prediction error| (kcal/mol)")
    ax.set_title("Fig 2 — Error vs OOD similarity (quartile trend)")
    ax.legend(fontsize=7)
    save(fig, "fig2_error_vs_similarity.png")

    # ---------------- Figure 3: MIPA ablation（不可用）----------------
    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    ax.axis("off")
    ax.text(0.5, 0.60, "Fig 3 — MIPA ablation", ha="center", fontsize=10,
            weight="bold")
    ax.text(0.5, 0.34,
            "NOT_AVAILABLE_RESOURCE_BLOCKED\n"
            "GATE B: mutation→residue mapping coverage\n"
            "strict EXACT = 13.12% (987/7523 pairs; test 189/1129)\n"
            "loose v2 = 61.35% — both below the 70% gate",
            ha="center", fontsize=8, color="darkred")
    save(fig, "fig3_mipa_ablation_NOT_AVAILABLE.png")

    # ---------------- Figure 4: Top-K enrichment curves ----------------
    z4 = np.load(C.V5 / "p4_conformal" / "p4_predictions.npz", allow_pickle=True)
    yt = z4["test_y"]
    gt = np.array([str(x) for x in z4["test_group"]])
    rs = np.random.RandomState(C.SEED)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    # 左：global —— 平均选中的真实 ΔΔG（越低越好）
    ax = axes[0]
    ks = np.arange(5, 105, 5)
    for key, lab, col in (("plain_s42", "model", "C0"),
                          ("plain_s43", "model (s43)", "C0"),
                          ("plain_s44", "model (s44)", "C0")):
        if key not in z4.files:
            continue
        p = z4[key]
        curve = [np.sort(yt - p)[:k].mean() for k in ks]
        ax.plot(ks, curve, "o-", ms=3, color=col, alpha=0.8,
                label="model (global pool)" if key == "plain_s42" else None)
    rnd = [np.mean([np.sort(yt[rs.permutation(len(yt))[:k]]).mean() for _ in range(60)])
           for k in ks]
    ax.plot(ks, rnd, "s--", ms=3, color="C1", label="random")
    ax.plot(ks, [np.sort(yt)[:k].mean() for k in ks], "^--", ms=3, color="C2",
            label="oracle")
    ax.axhline(yt.mean(), color="gray", ls=":", lw=1)
    ax.set_xlabel("K (budget)")
    ax.set_ylabel("Mean selected true $\\Delta\\Delta G$ (lower = better)")
    ax.set_title("Fig 4a — Global candidate pool\n"
                 "(dominated by between-complex differences)")
    ax.legend(fontsize=6.5)
    # 右：within-complex regret（越高越差），对比 oracle
    ax = axes[1]
    ks2 = np.arange(5, 55, 5)
    for key, lab, col in (("plain_s42", "model", "C0"),
                          ("plain_s43", "model (s43)", "C0"),
                          ("plain_s44", "model (s44)", "C0")):
        if key not in z4.files:
            continue
        p = z4[key]
        pts = {k: [] for k in ks2}
        for k in ks2:
            for _ in range(30):
                pass
        for cx in np.unique(gt):
            m = np.where(gt == cx)[0]
            if len(m) < 20:
                continue
            yy, pp = yt[m], p[m]
            sy = np.sort(yy)
            for k in ks2:
                kk = min(int(k), len(m))
                sel = np.argsort(-(-pp), kind="mergesort")[:kk]
                pts[k].append(yy[sel].mean() - sy[:kk].mean())
        ax.plot(ks2, [np.mean(pts[k]) for k in ks2], "o-", ms=3, color=col,
                label="model (=0 would mean oracle-level)" if key == "plain_s42" else None)
    ax.axhline(0, color="C2", ls="--", lw=1.2, label="oracle (regret = 0)")
    ax.set_xlabel("K (budget)")
    ax.set_ylabel("Regret vs within-complex oracle (kcal/mol; higher = worse)")
    ax.set_title("Fig 4b — Within-complex selection (real use case)\n"
                 "regret stays large => no usable within-target ranking")
    ax.legend(fontsize=6.5)
    save(fig, "fig4_topk_enrichment.png")

    # ---------------- Figure 5: Conformal calibration ----------------
    cal = pd.read_csv(C.V5 / "p4_conformal" / "p4_calibration.csv", low_memory=False)
    gp = z["oof_main"] if "oof_main" in z.files else None
    ytr = z["train_resid"] if "train_resid" in z.files else None
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    if ytr is not None:
        nom, emp = [], []
        for a in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50):
            q = float(np.quantile(ytr, 1 - a))
            mtr = d["mask_train"]
            emp.append(float(np.mean(d["y"][mtr] <= gp[mtr] + q)))
            nom.append(1 - a)
        ax.plot([0.5, 1.0], [0.5, 1.0], "k--", lw=0.9, label="ideal")
        ax.plot(nom, emp, "o-", color="C0", label="empirical (train OOF)")
        for n_, e_ in zip(nom, emp):
            ax.annotate("", xy=(n_, e_), xytext=(n_, n_),
                        arrowprops=dict(arrowstyle="->", lw=0.5, color="gray"))
        ax.set_xlabel("Nominal coverage")
        ax.set_ylabel("Empirical coverage")
        ax.set_title("Fig 5 — One-sided conformal calibration\n"
                     "(arrows show deviation from ideal)")
        ax.legend(fontsize=7)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no OOF available", ha="center")
    save(fig, "fig5_conformal_calibration.png")

    # ---------------- Figure 6: Split strictness vs performance ----------------
    pts = []
    v22 = json.loads((C.V4 if hasattr(C, "V4") else C.V5).as_posix()) if False else None
    # 用本项目已实测的三个口径（全部 PROJECT_EVIDENCE）
    pts = [("random split\n(V2.2 recorded)", 0.7409),
           ("PDB-level\n(V2.0 recorded)", 0.5409),
           ("family-level ≤40%\n(V2.0/V1.1 recorded)", 0.5198),
           ("Proteinea ≤30% official\n(V2.2 recorded)", 0.2375),
           ("complex+≤40%+family\n(V4/V5 ΔΔG, our strongest)", 0.3861)]
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    xs = np.arange(len(pts))
    ax.bar(xs, [p[1] for p in pts], color=["#d62728", "#ff7f0e", "#bcbd22",
                                           "#2ca02c", "#1f77b4"])
    for i, (lab, v) in enumerate(pts):
        ax.text(i, v + 0.015, f"{v:.4f}", ha="center", fontsize=7)
    ax.set_xticks(xs)
    ax.set_xticklabels([p[0] for p in pts], fontsize=6.5)
    ax.set_ylabel("Spearman / Pearson (metric differs per bar — see labels)")
    ax.set_title("Fig 6 — Split strictness dominates visible performance\n"
                 "(bars are NOT mutually comparable: different datasets & metrics)")
    save(fig, "fig6_split_strictness.png")

    print(f"\n生成 {len(made)} 张图: {made}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
