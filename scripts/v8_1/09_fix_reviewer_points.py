# -*- coding: utf-8 -*-
"""scripts/v8_1/09_fix_reviewer_points.py -- 回应审稿人三处文字修改 + 生成英文版图。

## 三处修改

1. `EXTERNAL_AUDIT_FINAL.md`：§8 门槛**不是**统计学通用规则 ⇒ 改为
   「with only five methods, method-level rank-correlation inference is low-powered」
   + 「under our prespecified practical threshold of 8 methods」。
2. 「external direction 与 internal 相反」⇒ 改为
   「The external panel showed strong positive rank agreement, which does not support a
   universal cross-metric decoupling and falls within the heterogeneous range observed
   internally.」
3. `BENCHMARK_NOVELTY_MATRIX.md` 版本不一致 ⇒ 同步到 V8.1 最终状态。

## 额外：英文版图（`figures_en/`）

英文版稿件需要英文图与英文图注；本脚本从 `FROZEN_RESULTS.json` 重新绘制**纯英文**图。
"""

from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
V8 = R / "experiments" / "v8"
FE = F / "figures_en"
FE.mkdir(exist_ok=True)

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.bbox"] = "tight"

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
D = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
W = pd.read_csv(V8 / "weighting" / "weighting_sensitivity.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
RS = pd.read_csv(V8 / "remedy" / "remedy_sequence.csv")
print("英文图数据就绪")

# ---------- 更新 1 & 2：EXTERNAL_AUDIT_FINAL.md ----------
p = F / "EXTERNAL_AUDIT_FINAL.md"
t = p.read_text(encoding="utf-8")
t = t.replace(
    "> **低功效警告（显式）**：外部可评价方法仅 **n_methods = 5**\n"
    "> （`Flex-ddG`、`DDGPred`、`ddG_FoldX`、`ESM-IF`、`proteinMPNN`），\n"
    "> 远低于秩相关推断所需的 **≥8** 门槛。\n",
    "> **Low-power warning (explicit).** Only **five** methods are externally evaluable\n"
    "> (`Flex-ddG`, `DDGPred`, `ddG_FoldX`, `ESM-IF`, `proteinMPNN`).\n"
    "> With only five methods, method-level rank-correlation inference is **low-powered**.\n"
    "> Under our **prespecified practical threshold of 8 methods**, the external analysis\n"
    "> is classified as `UNDERPOWERED`. This is a study-specific operating choice, **not**\n"
    "> a general statistical law.\n")
t = t.replace(
    "n_methods = 5   （门槛 >= 8）",
    "n_methods = 5   (below our prespecified practical threshold of 8)")
t = t.replace(
    "**⇒ 外部方向与内部相反**：EXACT_UNIQUE 上 `rho(Overall, Macro)` 为**正**（+0.90），\n"
    "**不是解耦** —— 但 n=5、p=0.083，\n"
    "**不足以支持或反驳**任何普遍主张。",
    "**⇒ 外部结果不构成解耦证据**：EXACT_UNIQUE 上 `rho(Overall, Macro)` 为**正**（+0.90）。\n"
    "该值**落在内部观测到的异质范围内**（主池 +0.82；全库 −0.27 至 +0.82），\n"
    "因此既不支持也不反驳「普遍解耦」；且 n=5、p=0.083，功效不足。")
t = t.replace(
    "> External evidence was **compatible with heterogeneous capability profiles**\n"
    "> but was **underpowered** regarding a universal cross-metric decoupling.\n"
    "> （低功效：n_methods = 5 < 8。）",
    "> The external panel showed **strong positive rank agreement**, which **does not\n"
    "> support a universal cross-metric decoupling** and **falls within the heterogeneous\n"
    "> range observed internally**. With only five methods this inference is low-powered.")
t = t.replace("**⚠️ 不因外部模型被加入而提升任何判定（§20）。**",
              "**⚠️ 不因外部模型被加入而提升任何判定（§20）。**\n\n"
              "> **英文定稿句**（供稿件直接引用）：*The external panel showed strong positive\n"
              "> rank agreement, which does not support a universal cross-metric decoupling and\n"
              "> falls within the heterogeneous range observed internally.*")
p.write_text(t, encoding="utf-8")
print("① ② EXTERNAL_AUDIT_FINAL.md 已改（低功效措辞 + 外部方向表述）")

# 同步到 04 脚本生成的 json，避免再次覆盖
j = json.load(open(A / "hermes_semantics_audit.json", encoding="utf-8"))
json.dump(j, open(A / "hermes_semantics_audit.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

# ---------- 更新 3：BENCHMARK_NOVELTY_MATRIX.md ----------
EXC = EP.set_index("model")
nov = f"""# BENCHMARK NOVELTY MATRIX — V8.1 FINAL FREEZE

> Synced to the V8.1 frozen state (`FROZEN_RESULTS.json`). Supersedes the V8 draft matrix.

## 0. Honest coverage of prior work (audited in V8 PHASE 1)

| Claim | Prior work | Our position |
|---|---|---|
| Leakage-aware / family-OOD PPI ΔΔG evaluation | Bushuiev et al. 2024 ([arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)) quantifies structural leakage (PDB-code 65%, sequence-30% 30%, SKEMPI v2.0 PPI-code 56%, expert grouping 0%) | **NOT ours. Do not claim.** |
| CATH-superfamily (family-OOD) splits + hard-mutation definition | CATH-ddG, Bioinformatics 2025 ([btaf228](https://doi.org/10.1093/bioinformatics/btaf228)) | **NOT ours. Do not claim.** |
| Per-target / per-structure Spearman reporting | already present in prior work | **NOT ours. Do not claim.** |
| Top-K / NDCG-style shortlist reporting | already present in prior work | **NOT ours. Do not claim.** |
| Overall leakage-aware benchmarking as a concept | prior work | **NOT ours. Do not claim.** |

## 1. External methods actually audited (V8.1 state)

| Method | Row-level predictions | Mapping used | Provenance | Evidence status |
|---|---|---|---|---|
| `DDGPred` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | literature-only provenance | TIER 2 |
| `ddG_FoldX` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | proprietary executable; predictions analysed | TIER 2 |
| `Flex-ddG` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | literature-only provenance | TIER 2 |
| `proteinMPNN` | yes (HERMES repo + CATH-ddG table) | exact key | inverse-folding, not ΔΔG-supervised | TIER 2 |
| `ESM-IF` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | inverse-folding, not ΔΔG-supervised | TIER 2 |
| `HERMES` | yes (official CSV) | exact 1:1 | cDNA117k stability fine-tune; **not** SKEMPI-binding-supervised; RMSE/MAE = NOT_APPLICABLE | TIER 2 |
| `PPIformer` | no official row-level predictions | — | upstream CPU blocker | NOT_RUN |
| `RDE-Network` | no official row-level predictions | — | needs self-run inference | NOT_RUN |
| `Prompt-DDG` | no official row-level predictions | — | no scorer CLI | NOT_RUN |
| `CATH-ddG` (itself) | **not released** | — | 88% of rows are its own training set | NOT_RUN |
| `PreFold-dG` | **not released** | No weights, no scorer | — | NOT_RUN |

## 2. Capability-resolved components — coverage

| Component | Covered by prior work? | Our status |
|---|---|---|
| Exact **between-target / within-target covariance-signal decomposition** `r = C_B + C_W` with identity QC | not identified | **contribution** (residual ≤ 3.3e-16) |
| **Controlled split-strictness** comparison of ΔC_B vs ΔC_W | not identified | **contribution** (between-fraction = {DL.get('fraction_between', float('nan')):.3f}) |
| **Target-weighting sensitivity** of the decomposition | not identified | **contribution** (C_B near-zero instability; C_W sign-stable) |
| **Matched-universe rank-reversal** analysis | not identified | **contribution** (range −0.27 … +0.82 across 10 strata) |
| **Fixed-K candidate-set geometry** audit, `D(n,K)=K(n−K)/[n(n+1)]` | not identified | **contribution** (40.8% of n≥5 targets have n ≤ K) |
| **Retrospective benchmark-guided remedy** attribution | not identified | **contribution** (retrospective only) |
| Integration into one capability-resolved evaluation framework | not identified | **contribution** |

## 3. Forbidden claims (V8.1, §19)

* first leakage-aware / first family-OOD / first per-target / first Top-K benchmark
* all existing benchmarks are invalid
* Overall Pearson is useless
* universal calibration–ranking orthogonality
* all PLMs fail; physics generally beats learning
* target size causes prediction difficulty
* HERMES log-odds is calibrated kcal/mol ΔΔG (**RMSE/MAE = NOT_APPLICABLE**)
* candidate-set floor proves SKEMPI ΔΔG values are intrinsically too similar
* benchmark prospectively predicted the successful remedy
* prospective wet-lab validity; universal model recommender
* `STRONGLY_SUPPORTED` merely because external models were added
* "recommend n ≥ 50 because fixed-K discrimination increases with n"
* any external capability claim computed on chain-ambiguous mappings

## 4. Positioning statement

> To our knowledge we did not identify prior work that **jointly** reports an exact
> between-/within-target covariance decomposition under matched, provenance-clean,
> family-OOD evaluation together with split-strictness, target-weighting, rank-reversal,
> fixed-K geometry and retrospective remedy analyses for PPI ΔΔG prediction.
> We make **no "first" claim** for leakage-aware benchmarking, family-OOD splitting,
> per-target reporting or Top-K evaluation.
"""
(F / "BENCHMARK_NOVELTY_MATRIX.md").write_text(nov, encoding="utf-8")
shutil.copy2(F / "BENCHMARK_NOVELTY_MATRIX.md",
             V8 / "BENCHMARK_NOVELTY_MATRIX_V8_1.md")
print("③ BENCHMARK_NOVELTY_MATRIX.md 已同步到 V8.1 最终状态")

# ================= 英文图 =================
K = 10
fig, ax = plt.subplots(figsize=(8.4, 4.6))
u = UR.sort_values("n_rows")
y = np.arange(len(u))
ax.barh(y, u.n_rows, color="#4c72b0", alpha=0.85, label="rows")
ax.barh(y, u.n_ge5, color="#dd8452", alpha=0.95, height=0.55,
        label="targets with n>=5")
ax.set_yticks(y)
ax.set_yticklabels([f"{a} ({h[:8]})" for a, h in zip(u.artifact, u.rowset_hash)],
                   fontsize=7.6)
for i, (n, g) in enumerate(zip(u.n_rows, u.n_ge5)):
    ax.text(n + 60, i, f"{n} / {g}", va="center", fontsize=7.0)
ax.set_xlabel("Number of mutation rows", fontsize=9.6)
ax.legend(fontsize=8.2, loc="lower right")
ax.set_title("Figure 1. Matched-universe strata\n"
             f"{len(u)} artifacts fall into {UR.rowset_hash.nunique()} distinct "
             "rowset hashes; only artifacts sharing a hash are directly comparable.",
             fontsize=9.8)
ax.grid(alpha=0.25, axis="x")
fig.savefig(FE / "Fig1_matched_universe_strata_EN.png"); plt.close(fig)

prim = D[D.rowset_hash == PU["rowset_hash"]].dropna(subset=["C_B_row", "C_W_row"])
fig, ax = plt.subplots(figsize=(6.2, 5.4))
ax.scatter(prim.C_B_row, prim.C_W_row, s=38, c="#4c72b0", alpha=0.82)
ax.axhline(0, color="#aaa", lw=0.9); ax.axvline(0, color="#aaa", lw=0.9)
lim = max(prim.C_B_row.abs().max(), prim.C_W_row.abs().max()) * 1.1
ax.plot([-lim, lim], [-lim, lim], "--", color="#999", lw=0.9)
res = float(np.nanmax(np.abs(prim.overall_pearson - (prim.C_B_row + prim.C_W_row))))
ax.set_xlabel(r"$C_B$  (between-target covariance contribution)", fontsize=9.6)
ax.set_ylabel(r"$C_W$  (within-target covariance contribution)", fontsize=9.6)
ax.set_title("Figure 2. Exact covariance-signal decomposition\n"
             f"universe={PU['rowset_hash'][:8]}, rows={PU['n_rows']}, "
             f"targets n>=5={PU['n_ge5']}, artifacts={PU['n_models']}, "
             f"method families={PU['n_method_families']}, 100% GREEN;\n"
             rf"$r_{{pooled}} = C_B + C_W$, residual $\leq$ {res:.1e}",
             fontsize=9.0)
ax.grid(alpha=0.25)
fig.savefig(FE / "Fig2_CB_CW_decomposition_EN.png"); plt.close(fig)

names = [k for k in H4 if k != "delta"]
x = np.arange(len(names))
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
axes[0].plot(x, [H4[k]["overall_pearson"] for k in names], "o-", lw=2, color="#c44e52")
axes[0].set_xticks(x); axes[0].set_xticklabels(names, fontsize=8.4, rotation=15)
axes[0].set_ylabel("Median pooled Pearson", fontsize=9.4)
axes[0].axhline(0, color="#aaa", lw=0.8); axes[0].grid(alpha=0.25)
axes[0].set_title("(a) Pooled correlation", fontsize=9.6)
axes[1].plot(x, [H4[k]["C_B"] for k in names], "s-", lw=2, color="#4c72b0",
             label=r"$C_B$ (between)")
axes[1].plot(x, [H4[k]["C_W"] for k in names], "d-", lw=2, color="#55a868",
             label=r"$C_W$ (within)")
axes[1].set_xticks(x); axes[1].set_xticklabels(names, fontsize=8.4, rotation=15)
axes[1].axhline(0, color="#aaa", lw=0.8); axes[1].grid(alpha=0.25)
axes[1].set_ylabel("Covariance contribution", fontsize=9.4)
axes[1].legend(fontsize=8.2); axes[1].set_title("(b) Decomposition", fontsize=9.6)
fig.suptitle("Figure 3. Controlled split inflation (identical rows and configurations; "
             "only the split regime changes)\n"
             f"universe = V6_2_split4 (7523 rows);  "
             rf"$\Delta r$={DL.get('delta_r',float('nan')):+.4f}, "
             rf"$\Delta C_B$={DL.get('delta_C_B',float('nan')):+.4f}, "
             rf"$\Delta C_W$={DL.get('delta_C_W',float('nan')):+.4f},  "
             f"between-fraction={DL.get('fraction_between',float('nan')):.3f}",
             fontsize=9.2, y=1.10)
fig.text(0.5, -0.06, "Note: L2_complex and L3_seq40 are identical in every metric here "
                     "and must not be treated as independent evidence.",
         ha="center", fontsize=8.2, color="#a33")
fig.savefig(FE / "Fig3_controlled_split_inflation_EN.png"); plt.close(fig)

b = RS[RS.label.str.startswith("V7.3 F1")]
r5 = RS[RS.step == "R5"].copy()
fig, ax = plt.subplots(figsize=(7.2, 5.2))
if len(b):
    b = b.iloc[0]
    r5["d_ov"] = r5.overall_r - b.overall_r
    r5["d_mc"] = r5.macro_spearman - b.macro_spearman
    ax.axhline(0, color="#aaa", lw=0.9); ax.axvline(0, color="#aaa", lw=0.9)
    cols = ["#55a868" if d > 0.02 else "#c44e52" for d in r5.d_mc]
    ax.scatter(r5.d_ov, r5.d_mc, s=100, c=cols, zorder=3,
               edgecolors="white", linewidths=1.1)
    for _, r in r5.iterrows():
        ax.annotate(r.label.replace("V7.3 ", ""), (r.d_ov, r.d_mc), fontsize=7.6,
                    xytext=(5, 4), textcoords="offset points")
ax.set_xlabel(r"$\Delta$ pooled Pearson vs F1 (mutation-only)", fontsize=9.4)
ax.set_ylabel(r"$\Delta$ macro per-target Spearman", fontsize=9.4)
ax.set_title("Figure 4. Retrospective matched remedy analysis (V7.3, identical 3351 rows)\n"
             "Green = both axes improve; red = pooled-score-only or negative.",
             fontsize=9.2)
ax.grid(alpha=0.25)
fig.savefig(FE / "Fig4_retrospective_remedy_EN.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.6))
m = GEO_DF.n <= 200
ax.plot(GEO_DF.n[m], GEO_DF.D[m], lw=2.2, color="#4c72b0")
ax.axvline(K, color="#c44e52", ls="--", lw=1.3)
ax.text(K + 1.2, GEO_DF.D[m].max() * 0.85, "n = K = 10", fontsize=8.2, color="#c44e52")
ns_ = GEO["n_star_analytic"]
ax.axvline(ns_, color="#55a868", ls=":", lw=1.3)
ax.text(ns_ + 1.8, GEO_DF.D[m].max() * 0.60, f"n* ~ {ns_:.1f}\n(peak)", fontsize=8.2,
        color="#55a868")
ax.set_xlabel("Candidate-set size n  (K = 10 fixed)", fontsize=9.4)
ax.set_ylabel("Discriminable interval  D(n, K)", fontsize=9.4)
ax.set_title("Figure 5. Fixed-K candidate-selection geometry\n"
             r"$D(n,K) = K(n-K)/[n(n+1)]$;  $n \leq K \Rightarrow D = 0$ "
             "(structurally non-informative)\n"
             f"D(20)={GEO['D_at_n20']:.4f} > D(50)={GEO['D_at_n50']:.4f} "
             "so the 'n >= 50' heuristic is not supported.",
             fontsize=9.0)
ax.grid(alpha=0.25)
fig.savefig(FE / "Fig5_fixed_K_geometry_EN.png"); plt.close(fig)

SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
cb = [f"C_B_{s}" for s in SCHEMES]; cw = [f"C_W_{s}" for s in SCHEMES]
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
for ax, cols, nm, cc in ((axes[0], cb, r"$C_B$ (between-target)", "#c44e52"),
                         (axes[1], cw, r"$C_W$ (within-target)", "#55a868")):
    A_ = W[cols].to_numpy(float)
    xx = np.arange(len(SCHEMES))
    for i in range(0, len(A_), 3):
        ax.plot(xx, A_[i], color=cc, alpha=0.15, lw=0.8)
    ax.plot(xx, np.nanmedian(A_, axis=0), "k-o", lw=2.1, ms=5.5)
    ax.axhline(0, color="#aaa", lw=0.9)
    ax.set_xticks(xx); ax.set_xticklabels(SCHEMES, fontsize=8.4)
    ax.grid(alpha=0.25)
    key = "C_B" if cc == "#c44e52" else "C_W"
    sf0 = int(W2[f"signflip_{key}_gt0"].sum())
    sf5 = int(W2[f"signflip_{key}_gt0.005"].sum())
    ax.set_title(f"{nm}\nsign flips: {sf0}/{len(A_)} (no threshold), "
                 f"{sf5}/{len(A_)} (|.|>0.005)", fontsize=9.2)
fig.suptitle("Figure 6. Target-weighting sensitivity: "
             r"$C_B$ direction is unstable near zero; "
             r"$C_W$ keeps its sign but its magnitude varies",
             fontsize=9.6, y=1.06)
fig.savefig(FE / "Fig6_weighting_sensitivity_EN.png"); plt.close(fig)

g = prim.dropna(subset=["overall_pearson", "macro_spearman"])
rho = float(ST[(ST.rowset_hash == PU["rowset_hash"]) &
               (ST.level == PU["level"])].rho_overall_vs_macro.iloc[0])
fig, ax = plt.subplots(figsize=(6.6, 5.0))
ax.scatter(g.overall_pearson, g.macro_spearman, s=42, c="#4c72b0", alpha=0.82)
for _, r in g.nlargest(3, "macro_spearman").iterrows():
    ax.annotate(r.key, (r.overall_pearson, r.macro_spearman), fontsize=7.2,
                xytext=(5, 4), textcoords="offset points")
ax.axhline(0, color="#aaa", lw=0.9); ax.axvline(0, color="#aaa", lw=0.9)
ax.set_xlabel("Pooled Pearson (axis A)", fontsize=9.4)
ax.set_ylabel("Macro per-target Spearman (axis B)", fontsize=9.4)
ax.set_title("Figure 7. Capability landscape, primary matched provenance-clean pool\n"
             f"universe={PU['rowset_hash'][:8]}, rows={PU['n_rows']}, "
             f"targets n>=5={PU['n_ge5']}, artifacts={PU['n_models']}, "
             f"method families={PU['n_method_families']}\n"
             f"rho = {rho:+.3f} in this pool; across {len(ST)} strata it ranges "
             f"{ST.rho_overall_vs_macro.min():+.2f} to {ST.rho_overall_vs_macro.max():+.2f}",
             fontsize=8.6)
ax.grid(alpha=0.25)
fig.savefig(FE / "Fig7_capability_landscape_EN.png"); plt.close(fig)

print(f"\n英文图 -> {FE}（{len(list(FE.glob('*.png')))} 张）")
for p_ in sorted(FE.glob("*.png")):
    print(f"   {p_.name}  {p_.stat().st_size/1024:.0f} KB")
