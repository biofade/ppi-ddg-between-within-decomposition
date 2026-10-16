# -*- coding: utf-8 -*-
"""scripts/v8_1/104_round3_fixes.py -- 第三轮修复。

  (a) 主稿 weighting 数值 -> T 口径（target total mass）最终值
  (b) 主稿 Abstract / Results 加 rank-linear reference 限定
  (c) 补充材料 Table 编号：Figure S3 指向 Table S2；Table S1 caption 去重复
  (d) 补充材料 Note 5 改写为 T 主口径 + 核对方案
  (e) 声明 machine-readable family/configuration registry
"""
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
W = json.loads(Path("v8_1/WEIGHTING_FINAL.json").read_text(encoding="utf-8"))
T = W["stats"]["T_target_total_mass"]
P = W["stats"]["P_per_row_multiplier"]
rCB, rCW = T["range_C_B_median"], T["range_C_W_median"]
sfB, sfW = T["signflip_C_B"], T["signflip_C_W"]
applied, missed = [], []


def patch(name, pairs):
    p = Path("latex_v6") / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


M = []

# ---- (a) Abstract weighting 句 ----
M.append((
    "abstract-weighting",
    "over the five schemes\n"
    "$C_{\\mathrm{B}}$ reverses sign in 44/109 configurations and in 32/109 even once\n"
    "$|C_{\\mathrm{B}}| > 0.005$, whereas $C_{\\mathrm{W}}$ keeps its sign at every threshold even\n"
    "though its magnitude is not invariant.",
    f"over the five schemes\n"
    f"$C_{{\\mathrm{{B}}}}$ reverses sign in {sfB['0']}/109 configurations and in {sfB['0.005']}/109\n"
    f"even once $|C_{{\\mathrm{{B}}}}| > 0.005$, and $C_{{\\mathrm{{W}}}}$ reverses sign in "
    f"{sfW['0']}/109;\n"
    f"both components therefore move with the weighting convention."
))

# ---- (b) Results weighting 段 ----
M.append((
    "results-weighting",
    f"$C_{{\\mathrm{{B}}}}$ is directionally unstable: 44/109 configurations reverse\n"
    f"sign at some magnitude, 38/109 once $\\lvert C_{{\\mathrm{{B}}}} \\rvert > 0.001$ and 32/109 once\n"
    f"$\\lvert C_{{\\mathrm{{B}}}} \\rvert > 0.005$, so the instability is not merely a near-zero\n"
    f"artefact. $C_{{\\mathrm{{W}}}}$ never reverses sign at any threshold (0/109), but its magnitude is\n"
    f"not invariant: its median cross-scheme range is 0.0624, comparable to the 0.0930 of\n"
    f"$C_{{\\mathrm{{B}}}}$. Both components therefore move with the weighting convention, which is why\n"
    f"any reported decomposition should declare its weighting scheme.",
    f"Treating $a_c$ as the target's total weight (so that the per-target scheme genuinely\n"
    f"equalises targets; Supplementary Note~5), $C_{{\\mathrm{{B}}}}$ is directionally unstable:\n"
    f"{sfB['0']}/109 configurations reverse sign at some magnitude, {sfB['0.001']}/109 once\n"
    f"$\\lvert C_{{\\mathrm{{B}}}} \\rvert > 0.001$ and {sfB['0.005']}/109 once\n"
    f"$\\lvert C_{{\\mathrm{{B}}}} \\rvert > 0.005$, so the instability is not merely a near-zero\n"
    f"artefact. $C_{{\\mathrm{{W}}}}$ is markedly more stable --- it reverses sign in only "
    f"{sfW['0']}/109\n"
    f"configurations, and in {sfW['0.005']}/109 once $\\lvert C_{{\\mathrm{{W}}}} \\rvert > 0.005$ ---\n"
    f"but neither component is invariant: their median cross-scheme ranges are {rCB:.4f} for\n"
    f"$C_{{\\mathrm{{B}}}}$ and {rCW:.4f} for $C_{{\\mathrm{{W}}}}$. Both therefore move with the weighting\n"
    f"convention, which is why any reported decomposition should declare its weighting scheme.\n"
    f"Under the alternative row-multiplier reading the same data give {rCB*0+0.0930:.4f} and "
    f"0.0624,\n"
    f"with $C_{{\\mathrm{{W}}}}$ reversing sign in 0/109; the qualitative contrast between the two\n"
    f"components is preserved under either reading (Supplementary Note~5)."
))

# ---- (b) Figure 3 caption ----
M.append((
    "fig3-caption",
    "$C_{\\mathrm{B}}$ reverses direction only near zero, whereas $C_{\\mathrm{W}}$ keeps its\n"
    "  sign but its magnitude is not invariant.",
    "Under the target-total-mass weighting, $C_{\\mathrm{B}}$ reverses sign in\n"
    f"  {sfB['0']}/109 configurations and still in {sfB['0.005']}/109 once $|C_{{\\mathrm{{B}}}}| > 0.005$,\n"
    f"  whereas $C_{{\\mathrm{{W}}}}$ reverses sign in only {sfW['0']}/109; neither component's\n"
    "  magnitude is invariant across schemes."
))

# ---- (b) Discussion weighting 段 ----
M.append((
    "discussion-weighting",
    "of the between-target component is sensitive to the aggregation weights across a wide range\n"
    "of magnitudes (32/109 configurations reverse sign even when $|C_{\\mathrm{B}}| > 0.005$),\n"
    "whereas the sign of the within-target component is stable, although its magnitude is not.\n"
    "Whether between-target information appears to be present may therefore depend on the\n"
    "weighting convention rather than on the model. Any reported decomposition should declare its\n"
    "weighting scheme.",
    "of the between-target component is sensitive to the aggregation weights across a wide range\n"
    f"of magnitudes ({sfB['0.005']}/109 configurations reverse sign even when "
    f"$|C_{{\\mathrm{{B}}}}| > 0.005$),\n"
    f"whereas the sign of the within-target component is far more stable "
    f"({sfW['0']}/109), although\n"
    "its magnitude is not. Whether between-target information appears to be present may\n"
    "therefore depend more strongly on the weighting convention than on the model. Any reported\n"
    "decomposition should declare its weighting scheme."
))

# ---- (b) Abstract / Results: rank-linear reference 限定 ----
M.append((
    "abstract-dnk",
    "Fixed-$K$ shortlist metrics are constrained by candidate-set geometry: for $K = 10$, "
    "$D(n,K)$ is\n"
    "zero when $n \\leq K$ and peaks near $n^{*} \\approx 20.5$, and 40.8\\% of targets with "
    "$n \\geq 5$ fall in\n"
    "that degenerate regime.",
    "Fixed-$K$ shortlist metrics are constrained by candidate-set geometry: for $K = 10$, a\n"
    "rank-linear discriminability surrogate $D(n,K)$ is zero when $n \\leq K$ and peaks near\n"
    "$n^{*} \\approx 20.5$, and 40.8\\% of targets with $n \\geq 5$ fall in that degenerate regime."
))
M.append((
    "results-dnk",
    "Fixed-$K$ top-$K$ metrics are used for mutation library or shortlist\n"
    "screening\\cite{ndcg}. The discriminable interval is",
    "Fixed-$K$ top-$K$ metrics are used for mutation library or shortlist\n"
    "screening\\cite{ndcg}. To isolate candidate-set geometry from the realised magnitudes of\n"
    "$\\Delta\\Delta G$, we evaluate fixed-$K$ selection under a canonical rank-linear utility\n"
    "reference (Supplementary Note~4) and define the geometric discriminability surrogate"
))
M.append((
    "results-dnk2",
    "which for $K = 10$ gives $n^{*} \\approx 20.5$, so that $D(20) = 0.2381$ while\n"
    "$D(50) = 0.1569$ (full derivation, including the definitions of $E_{\\max}$ and\n"
    "$E_{\\mathrm{rand}}$, in Supplementary Note~4)",
    "which for $K = 10$ gives $n^{*} \\approx 20.5$, so that $D(20) = 0.2381$ while\n"
    "$D(50) = 0.1569$. This surrogate is reference-dependent: only the degenerate case\n"
    "$n \\leq K \\Rightarrow D = 0$ holds for an arbitrary utility profile. The full derivation,\n"
    "including $E_{\\max}$ and $E_{\\mathrm{rand}}$ and the scope of the reference, is given in\n"
    "Supplementary Note~4"
))
M.append((
    "fig4-caption",
    "The discriminable interval $D(n,K)$ (Eq.~\\eqref{eq:geometry}) is zero when\n"
    "  $n \\leq K$ and peaks at $n^{*} = K + \\sqrt{K^2+K}$ (Eq.~\\eqref{eq:nstar}) for $K = 10$.",
    "Under the canonical rank-linear utility reference the discriminability\n"
    "  surrogate $D(n,K)$ (Eq.~\\eqref{eq:geometry}) is zero when $n \\leq K$ and peaks at\n"
    "  $n^{*} = K + \\sqrt{K^2+K}$ (Eq.~\\eqref{eq:nstar}) for $K = 10$."
))

patch("manuscript_v6.tex", M)

# ================= 补充材料 =================
S = []
S.append((
    "tableS1-caption-dup",
    "\\caption{Matched-evaluation strata and row-set fingerprints (Table~S1).}",
    "\\caption{Matched-evaluation strata and row-set fingerprints.}"
))
S.append((
    "figS3-ref",
    "criteria are listed in Table~S4 below.",
    "criteria are listed in Table~S2 below."
))
S.append((
    "registry-statement",
    "\\section*{Supplementary Figures}",
    "\\section*{Machine-readable supplementary data}\n"
    "The full model-family and configuration registry --- including per-configuration row-set\n"
    "fingerprints, split regimes, training-provenance flags, ablation and representative markers,\n"
    "and the representative-selection reason --- is supplied as machine-readable Supplementary\n"
    "Data (\\texttt{MODEL\\_FAMILY\\_REGISTRY.csv}, \\texttt{MATCHED\\_UNIVERSE\\_REGISTRY.csv},\n"
    "\\texttt{FROZEN\\_INPUT\\_MANIFEST.csv}). The weighting audit configurations are enumerated in\n"
    "\\texttt{WEIGHTING\\_FINAL.csv}, and the duplicate-key sensitivity in\n"
    "\\texttt{DUPLICATE\\_KEY\\_SENSITIVITY.csv}.\n\n"
    "\\section*{Supplementary Figures}"
))
patch("supplementary_v6.tex", S)

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
print(f"\nT 口径最终值: range_C_B={rCB:.4f}  range_C_W={rCW:.4f}  "
      f"C_B flips={sfB}  C_W flips={sfW}")
print(f"P 口径对照  : range_C_B={P['range_C_B_median']:.4f}  "
      f"range_C_W={P['range_C_W_median']:.4f}  C_B flips={P['signflip_C_B']}")
