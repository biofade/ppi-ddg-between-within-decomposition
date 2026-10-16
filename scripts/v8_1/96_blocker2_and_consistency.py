# -*- coding: utf-8 -*-
"""scripts/v8_1/96_blocker2_and_consistency.py -- 修复第二轮审计的其余问题。

  BLOCKER 2 : D(n,K) 推导缺 rank-linear utility 假设
  MAJOR  4  : Figure 6 "matched 1,781" 措辞 + caption；Discussion "ordering agreed"
  MAJOR  3  : duplicate-key grouped sensitivity（新增 Supplementary Note 6）
  MINOR  5  : Figure S1/S2/S3 与 Table S3/S4 编号
  MINOR  6  : Abstract novelty / 81.1% 措辞
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")
W = json.loads(Path("v8_1/WEIGHTING_CONSISTENT.json").read_text(encoding="utf-8"))
applied, missed = [], []


def patch(name, pairs):
    p = OUT / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


# ===================== SUPPLEMENTARY =====================
S = []

# ---- BLOCKER 2: D(n,K) 明确 rank-linear reference ----
S.append((
    "dnk-rank-linear",
    "Let a target contain $n$ candidate mutations carrying utilities $u_i$ (here the true\n"
    "$|\\Delta\\Delta G|$ or any monotone transform of it), and let\n"
    "$E = \\sum_{i \\in S_K} u_i / \\sum_{i=1}^{n} u_i$ denote the share of the total utility captured\n"
    "by a selected set $S_K$ of $K$ candidates. Two references define the discriminable interval:",
    "To separate candidate-set geometry from the realised magnitudes of $\\Delta\\Delta G$, we\n"
    "evaluate fixed-$K$ selection against a \\emph{canonical rank-linear utility reference}.\n"
    "Rank the $n$ candidates of a target from most to least favourable and assign the $j$-th\n"
    "ranked candidate the reference utility\n"
    "\\begin{equation}\n"
    "  u_{(j)} = n - j + 1, \\qquad j = 1,\\dots,n,\n"
    "  \\label{eq:rankutil}\n"
    "\\end{equation}\n"
    "so that utilities are equally spaced and strictly decreasing in rank. This reference is a\n"
    "deliberate idealisation: it fixes the \\emph{shape} of the utility profile and asks how much\n"
    "of the total utility a size-$K$ shortlist can capture as a function of $n$ and $K$ alone.\n"
    "It is not a model of any particular $\\Delta\\Delta G$ distribution, and the resulting\n"
    "quantity is a geometric surrogate for shortlist discriminability rather than a property of\n"
    "real affinity data. Let $E = \\sum_{i \\in S_K} u_i / \\sum_{i=1}^{n} u_i$ denote the share of\n"
    "total utility captured by a selected set $S_K$ of $K$ candidates. Two references define the\n"
    "discriminable interval:"
))

S.append((
    "dnk-emax",
    "  \\item $E_{\\max}$ --- the share attained by the optimal (oracle) selection. Sorting $u$ in\n"
    "        decreasing order and taking the top $K$ gives\n"
    "        \\begin{equation}\n"
    "          E_{\\max} = \\frac{K\\,(2n-K+1)}{n\\,(n+1)} ;\n"
    "          \\label{eq:emax}\n"
    "        \\end{equation}",
    "  \\item $E_{\\max}$ --- the share attained by the optimal (oracle) selection. Under\n"
    "        Eq.~\\eqref{eq:rankutil} the top $K$ candidates have utilities\n"
    "        $n, n-1, \\dots, n-K+1$, and $\\sum_{j=1}^{n} u_{(j)} = n(n+1)/2$, so\n"
    "        \\begin{equation}\n"
    "          E_{\\max} = \\frac{\\sum_{j=1}^{K}(n-j+1)}{n(n+1)/2}\n"
    "                    = \\frac{K\\,(2n-K+1)}{n\\,(n+1)} .\n"
    "          \\label{eq:emax}\n"
    "        \\end{equation}\n"
    "        The closed form depends on the reference Eq.~\\eqref{eq:rankutil}; for an arbitrary\n"
    "        utility vector $E_{\\max} = \\sum_{j\\le K} u_{(j)} / \\sum_{j\\le n} u_{(j)}$ and is\n"
    "        not determined by $n$ and $K$ alone."
))

S.append((
    "dnk-emrand",
    "  \\item $E_{\\mathrm{rand}}$ --- the share expected under uniform random selection of $K$ of the\n"
    "        $n$ candidates. By symmetry of the uniform distribution over subsets,\n"
    "        \\begin{equation}\n"
    "          E_{\\mathrm{rand}} = \\frac{K}{n} .\n"
    "          \\label{eq:erand}\n"
    "        \\end{equation}",
    "  \\item $E_{\\mathrm{rand}}$ --- the share expected under uniform random selection of $K$ of the\n"
    "        $n$ candidates. By symmetry of the uniform distribution over subsets, each candidate\n"
    "        is selected with probability $K/n$ irrespective of its utility, so\n"
    "        \\begin{equation}\n"
    "          E_{\\mathrm{rand}} = \\frac{K}{n} .\n"
    "          \\label{eq:erand}\n"
    "        \\end{equation}\n"
    "        This reference is exact under any utility profile, including the rank-linear one."
))

S.append((
    "dnk-limits",
    "\\emph{Applicability and limits.} $D$ is a purely combinatorial property of candidate-set\n"
    "geometry under a uniform-utility reference. It is not a model-performance metric, it does not\n"
    "measure the $\\Delta\\Delta G$ utility range of any real dataset, and it does not imply that\n"
    "SKEMPI $\\Delta\\Delta G$ values are intrinsically similar. Establishing that would require an\n"
    "analysis using realised utility magnitudes, which we did not perform.",
    "\\emph{Applicability and limits.} Everything above is conditional on the canonical\n"
    "rank-linear reference Eq.~\\eqref{eq:rankutil}. $D(n,K)$ is therefore a geometric\n"
    "surrogate for shortlist discriminability, not a statement about arbitrary utility\n"
    "distributions: with any other utility profile $E_{\\max}$ changes and so does $D$.\n"
    "It is not a model-performance metric, it does not measure the $\\Delta\\Delta G$ utility range\n"
    "of any real dataset, and it does not imply that SKEMPI $\\Delta\\Delta G$ values are\n"
    "intrinsically similar. Establishing that would require an analysis using realised utility\n"
    "magnitudes, which we did not perform. What is reference-free is the degenerate case: for\n"
    "$n \\leq K$ a size-$K$ shortlist selects every candidate, so $E_{\\max} = E_{\\mathrm{rand}} = 1$\n"
    "and $D = 0$ under \\emph{any} utility profile. That statement, and the target counts it\n"
    "implies, does not depend on Eq.~\\eqref{eq:rankutil}."
))

# ---- MAJOR 4: Figure S3 编号与指向 ----
S.append((
    "figS-numbering-fix",
    "  \\caption{\\textbf{Figure S1.} Matched-evaluation strata grouped by SHA-256 row-set",
    "  \\caption{Matched-evaluation strata grouped by SHA-256 row-set"
))
S.append((
    "figS2-numbering-fix",
    "  \\caption{\\textbf{Figure S2.} Capability profiles of matched, provenance-clean ablation",
    "  \\caption{Capability profiles of matched, provenance-clean ablation"
))
S.append((
    "figS3-numbering-fix",
    "  \\caption{\\textbf{Figure S3.} Pool sensitivity of cross-metric rank agreement across",
    "  \\caption{Pool sensitivity of cross-metric rank agreement across"
))
S.append((
    "tableS4-label-fix",
    "(Table~S4).}\\label{tab:s4pools}",
    "Table~S4.}\\label{tab:s4pools}"
))
S.append((
    "tableS4-caption-fix",
    "  \\caption{Model pools used in the rank-agreement sensitivity analysis Table~S4.}",
    "  \\caption{Model pools used in the rank-agreement sensitivity analysis.}"
))

# ---- MINOR 6 / MAJOR 3: 主稿 ----
M = []
M.append((
    "abstract-novelty",
    "We introduce an exact covariance-signal decomposition that\n"
    "writes pooled Pearson correlation as $r = C_{\\mathrm{B}} + C_{\\mathrm{W}}$, and apply it to",
    "We formulate an exact between-/within-target covariance decomposition for\n"
    "PPI $\\Delta\\Delta G$ benchmarking, which writes pooled Pearson correlation as\n"
    "$r = C_{\\mathrm{B}} + C_{\\mathrm{W}}$, and apply it to"
))
M.append((
    "abstract-attribution",
    "absolute-change basis the between-target component accounts for 81.1\\% of the decline while\n"
    "the within-target component remains positive.",
    "absolute-change basis the between-target component accounts for 81.1\\% of the summed\n"
    "absolute component change while the within-target component remains positive."
))
M.append((
    "abstract-weighting",
    "Target weighting changes how these\n"
    "contributions should be read: $C_{\\mathrm{B}}$ reverses direction only near zero\n"
    "(35/109 configurations), whereas $C_{\\mathrm{W}}$ keeps its sign at every threshold\n"
    "although its magnitude is not invariant.",
    "Target weighting changes how these contributions should be read: over the five schemes\n"
    "$C_{\\mathrm{B}}$ reverses sign in 44/109 configurations and in 32/109 even once\n"
    "$|C_{\\mathrm{B}}| > 0.005$, whereas $C_{\\mathrm{W}}$ keeps its sign at every threshold even\n"
    "though its magnitude is not invariant."
))
M.append((
    "results-weighting",
    "capped at 30), $C_{\\mathrm{B}}$ shows directional instability concentrated near zero:\n"
    "35/109 configurations reverse sign with no threshold, whereas none do once\n"
    "$\\lvert C_{\\mathrm{B}} \\rvert$ exceeds 0.005 (0/109). $C_{\\mathrm{W}}$ never reverses sign at\n"
    "any threshold (0/109), but its magnitude is not invariant --- its median range is 0.0624,\n"
    "larger than the 0.0183 of $C_{\\mathrm{B}}$. Any reported decomposition should declare its\n"
    "weighting scheme.",
    "capped at 30), $C_{\\mathrm{B}}$ is directionally unstable: 44/109 configurations reverse\n"
    "sign at some magnitude, 38/109 once $\\lvert C_{\\mathrm{B}} \\rvert > 0.001$ and 32/109 once\n"
    "$\\lvert C_{\\mathrm{B}} \\rvert > 0.005$, so the instability is not merely a near-zero\n"
    "artefact. $C_{\\mathrm{W}}$ never reverses sign at any threshold (0/109), but its magnitude is\n"
    "not invariant: its median cross-scheme range is 0.0624, comparable to the 0.0930 of\n"
    "$C_{\\mathrm{B}}$. Both components therefore move with the weighting convention, which is why\n"
    "any reported decomposition should declare its weighting scheme."
))
M.append((
    "discussion-weighting",
    "The weighting-sensitivity audit yields a conclusion of direct practical relevance: the sign\n"
    "of the between-target component is sensitive to aggregation weights near zero, whereas the\n"
    "sign of the within-target component is stable. When a reported $C_{\\mathrm{B}}$ is close to\n"
    "zero, the conclusion about whether between-target information is present may depend on the\n"
    "weighting convention rather than on the model. Any reported decomposition should declare its\n"
    "weighting scheme.",
    "The weighting-sensitivity audit yields a conclusion of direct practical relevance: the sign\n"
    "of the between-target component is sensitive to the aggregation weights across a wide range\n"
    "of magnitudes (32/109 configurations reverse sign even when $|C_{\\mathrm{B}}| > 0.005$),\n"
    "whereas the sign of the within-target component is stable, although its magnitude is not.\n"
    "Whether between-target information appears to be present may therefore depend on the\n"
    "weighting convention rather than on the model. Any reported decomposition should declare its\n"
    "weighting scheme."
))
M.append((
    "fig6-caption",
    "(A) Pooled\n"
    "versus within-target performance for five independently developed methods on a matched,\n"
    "unambiguous mapping.",
    "(A) Pooled versus within-target performance for five independently developed methods, each\n"
    "evaluated on its own uniquely mapped subset of a common unambiguous annotation universe."
))
M.append((
    "discussion-order",
    "Across five\n"
    "independently developed methods the ordering of methods agreed between pooled and\n"
    "within-target performance, a strongly positive agreement that neither supports nor refutes a\n"
    "universal decoupling of the two axes.",
    "Across five independently developed methods the rankings were strongly positively aligned\n"
    "between pooled and within-target performance (Spearman $\\rho = +0.90$), with a single\n"
    "adjacent swap; this agreement neither supports nor refutes a universal decoupling of the two\n"
    "axes."
))

patch("supplementary_v6.tex", S)
patch("manuscript_v6.tex", M)

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
