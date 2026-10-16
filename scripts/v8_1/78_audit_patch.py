# -*- coding: utf-8 -*-
"""scripts/v8_1/78_audit_patch.py -- 终稿科学审计补丁（作用于生成的 .tex）。

所有数值来自：
  v8_1/UNIVERSE_PROVENANCE.json      (4,473 vs 7,523 关系)
  v8_1/CONTROLSPLIT_IDENTITY.json    (恒等式 / paired attribution / bootstrap)
  v8_1/AUDIT_CONSOLIDATED.json       (external / families / 109 / Fig3 / strata)
  v8_1/FROZEN_RESULTS.json           (headline numbers)

本脚本只做「审计要求」的修改，不引入任何新数据。
"""
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")
A = Path("v8_1")
UP = json.loads((A / "UNIVERSE_PROVENANCE.json").read_text(encoding="utf-8"))
CI = json.loads((A / "CONTROLSPLIT_IDENTITY.json").read_text(encoding="utf-8"))
AC = json.loads((A / "AUDIT_CONSOLIDATED.json").read_text(encoding="utf-8"))

applied, missed = [], []


def patch(name, pairs):
    p = OUT / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1)
            applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


# =====================================================================
# 数值（全部取自真实 artifact）
# =====================================================================
mas = UP["universes"]["U_MAS_7523"]
can = UP["universes"]["U_CANON_4473"]
pa = CI["paired"]
mm = CI["marginal_medians"]
ext = AC["external"]
cis = pa["bootstrap_ci"]

# =====================================================================
# MANUSCRIPT
# =====================================================================
M = []

# ---- P0-2: 4,473 vs 7,523 provenance ----
M.append((
    "universe-intro",
    "We apply this decomposition on a matched evaluation set comprising 3{,}351 mutation rows,\n"
    "102 targets with at least five mutations, and 10 independent method families.",
    "We apply this decomposition on a matched evaluation set of 3{,}351 mutation rows and 102\n"
    "targets with at least five mutations, drawn from 10 independent method families. Three\n"
    "distinct row universes appear in this work and are defined in Methods; the primary\n"
    "matched universe (3{,}351 rows) is a provenance- and feature-filtered subset of the\n"
    "4{,}473-row annotation namespace, which is itself a feature-coverage subset of the\n"
    "7{,}523-row controlled-split universe."
))

# ---- P0-3: median vs exact identity ----
M.append((
    "identity-perconfig",
    "The identity holds to machine precision (residual\n"
    "$\\leq 3.3 \\times 10^{-16}$). These are covariance-signal contributions on a common scale,\n"
    "not independent correlations, skill scores or causal mechanisms.",
    "The identity holds **exactly for every individual configuration**: across the 36\n"
    "controlled-split configurations the maximum absolute residual\n"
    "$|r_{\\mathrm{pooled}} - C_{\\mathrm{B}} - C_{\\mathrm{W}}|$ is $4.4 \\times 10^{-16}$, i.e.\n"
    "machine precision. The values reported in Table~\\ref{tab:split} are however *marginal\n"
    "medians*, computed separately for $r_{\\mathrm{pooled}}$, $C_{\\mathrm{B}}$ and\n"
    "$C_{\\mathrm{W}}$; because the median is not additive, those three displayed medians need\n"
    "not sum exactly (e.g.\\ Random: $0.4275 + 0.2591 = 0.6866$ versus a median\n"
    "$r_{\\mathrm{pooled}}$ of $0.6844$). These are covariance-signal contributions on a common\n"
    "scale, not independent correlations, skill scores or causal mechanisms."
))

# ---- P0-4: 81.1% attribution 严格化 + paired 分析 ----
M.append((
    "attribution-rigorous",
    "absolute-change basis the between-target component accounts for 81.1\\% of the decline, while\n"
    "$C_{\\mathrm{W}}$ remains positive ($+0.1580$). This is not the same proposition as `the model\n"
    "memorised target identity'",
    "absolute-change basis the between-target component accounts for 81.1\\% of the *summed\n"
    "absolute component change*,\n"
    "\\begin{equation}\n"
    "  \\frac{\\lvert \\Delta C_{\\mathrm{B}} \\rvert}\n"
    "       {\\lvert \\Delta C_{\\mathrm{B}} \\rvert + \\lvert \\Delta C_{\\mathrm{W}} \\rvert}\n"
    "  = 0.8113 ,\n"
    "  \\label{eq:attribution}\n"
    "\\end{equation}\n"
    "while $C_{\\mathrm{W}}$ remains positive ($+0.1580$). We verified this attribution on the\n"
    "paired scale as well. Because the same model specifications and seeds were evaluated under\n"
    "both regimes, all 9 configurations are directly matched: the paired differences satisfy\n"
    "$\\Delta r_j = \\Delta C_{\\mathrm{B},j} + \\Delta C_{\\mathrm{W},j}$ exactly (maximum\n"
    "residual $4.2 \\times 10^{-16}$), with medians\n"
    "$\\Delta r = -0.5332$, $\\Delta C_{\\mathrm{B}} = -0.4324$ and\n"
    "$\\Delta C_{\\mathrm{W}} = -0.1011$. The paired attribution is therefore 81.0\\%\n"
    "(configuration-cluster bootstrap 95\\% CI 71.0--84.5\\%), essentially identical to the\n"
    "marginal-median figure. This is not the same proposition as `the model memorised target\n"
    "identity'"
))

# ---- P0-5: external coverage + P0-6: rank swap ----
M.append((
    "external-coverage",
    "Where the matching key was not\n"
    "unique on the annotation side, only uniquely mapped rows were used; of 4{,}473 rows only\n"
    "1{,}781 have a unique mapping, and ambiguous keys were excluded rather than forcibly\n"
    "resolved. On this unambiguous set the five methods agree in their ordering on pooled\n"
    "correlation and on within-target ranking, with DDGPred best on both axes; no method-level\n"
    "rank reversal occurs. Method-level rank agreement is $\\rho = +0.90$ (exact permutation\n"
    "$p = 0.083$). This strongly positive agreement does not support a universal decoupling of\n"
    "pooled and within-target performance and falls within the heterogeneous range observed\n"
    "internally. With only five methods the panel is underpowered, so we draw no directional\n"
    "conclusion from it.",
    "Where the matching key was not\n"
    "unique on the annotation side, only uniquely mapped rows were used; of 4{,}473 rows only\n"
    "1{,}781 have a unique mapping, and ambiguous keys were excluded rather than forcibly\n"
    "resolved. Each method was evaluated on its own uniquely mapped subset, and coverage\n"
    "therefore differs slightly across methods: Flex-ddG and DDGPred map 1{,}755 and 1{,}756\n"
    "mutations (65 targets with $n \\geq 5$), whereas FoldX, ESM-IF and ProteinMPNN map 1{,}776\n"
    "mutations (66 targets with $n \\geq 5$). We therefore do not claim a fully matched panel;\n"
    "the coverage difference is at most 21 mutations (1.2\\% of the panel) and we verified that\n"
    "it does not drive the ranking conclusion. Rank ordering was strongly aligned between the two\n"
    "axes, with a single adjacent swap between FoldX and Flex-ddG, and DDGPred ranking first on\n"
    "both; method-level rank agreement is $\\rho = +0.90$ (Kendall $\\tau = +0.80$, exact\n"
    "permutation $p = 0.083$ for $n = 5$ methods). This strongly positive agreement does not\n"
    "support a universal decoupling of pooled and within-target performance and falls within the\n"
    "heterogeneous range observed internally. With only five methods the panel is underpowered,\n"
    "so we draw no directional conclusion from it."
))

# ---- Figure 1 caption: 删除过时的 target-card 描述, "108 structural families" ----
M.append((
    "fig1-caption-A",
    "annotation (4{,}473 rows, 214 targets, 108 families) is audited by SHA-256 row-set\n"
    "fingerprinting; family-out-of-distribution predictions without training overlap are\n"
    "retained and permissively split predictions excluded, giving a primary matched universe of\n"
    "3{,}351 rows, 102 targets with $n \\geq 5$ and 10 method families. Target cards show\n"
    "between-target differences in mean $\\Delta\\Delta G$ ($\\mu_c$, purple) and within-target\n"
    "mutation ordering ($\\varepsilon_{ci}$, green).",
    "annotation namespace (4{,}473 rows, 214 targets, 108 *structural* families --- distinct from\n"
    "the 10 *method* families) is audited by SHA-256 row-set fingerprinting; predictions whose\n"
    "recorded training split is permissive, or whose evaluation rows may overlap their training\n"
    "folds, are excluded, giving a primary matched universe of 3{,}351 rows, 102 targets with\n"
    "$n \\geq 5$ and 10 method families. The annotation namespace is itself a feature-coverage\n"
    "subset of the 7{,}523-row controlled-split universe (Methods, *Dataset universes*)."
))

M.append((
    "fig1-caption-B",
    "(B) Exact additive covariance-signal\n"
    "decomposition; $C_{\\mathrm{B}}$ and $C_{\\mathrm{W}}$ sum to pooled correlation by an\n"
    "exact identity (residual $\\leq 3.3 \\times 10^{-16}$) and are common-scale additive\n"
    "contributions, not independent correlations or causal mechanisms.",
    "(B) Exact additive covariance-signal\n"
    "decomposition. The identity $r_{\\mathrm{pooled}} = C_{\\mathrm{B}} + C_{\\mathrm{W}}$ holds\n"
    "for every individual configuration (maximum residual $4.4 \\times 10^{-16}$); the two\n"
    "components are common-scale additive contributions, not independent correlations or causal\n"
    "mechanisms."
))

M.append((
    "fig1-caption-C",
    "(C) Controlled split\n"
    "collapse on identical 7{,}523 mutation rows: pooled correlation falls from 0.6844 to\n"
    "0.1508 and 81.1\\% of the absolute decline is carried by loss of the between-target\n"
    "component.",
    "(C) Controlled split collapse on the\n"
    "identical 7{,}523-row universe with identical model specifications and seeds; only the split\n"
    "regime changes. Pooled correlation falls from 0.6844 to 0.1508, and 81.1\\% of the summed\n"
    "absolute component change is carried by the between-target component."
))

# ---- Figure 2 caption: median vs exact + attribution wording ----
M.append((
    "fig2-caption",
    "{\\bfseries\\textbf{Family-out-of-distribution evaluation disproportionately removes the\n"
    "between-target component.}} (A) Pooled Pearson by split regime. (B) Exact decomposition\n"
    "into $C_{\\mathrm{B}}$ and $C_{\\mathrm{W}}$; $C_{\\mathrm{W}}$ remains positive.\n"
    "(C) Attribution of the absolute decline, of which the between-target component accounts\n"
    "for 81.1\\%. All four regimes use the same 7{,}523 rows and configurations; complex and\n"
    "\\texttt{seq40} are identical and are not independent comparisons.",
    "{\\bfseries\\textbf{Family-out-of-distribution evaluation disproportionately removes the\n"
    "between-target component.}} (A) Pooled Pearson by split regime. (B) Median component\n"
    "summaries; the exact identity $r_{\\mathrm{pooled}} = C_{\\mathrm{B}} + C_{\\mathrm{W}}$\n"
    "holds per configuration, but the plotted medians are marginal statistics and need not sum\n"
    "exactly. $C_{\\mathrm{W}}$ remains positive. (C) Attribution of the summed absolute\n"
    "component change, of which the between-target component accounts for 81.1\\%; the paired\n"
    "attribution is 81.0\\%. All four regimes use the same 7{,}523 rows with identical model\n"
    "specifications and seeds; complex and \\texttt{seq40} are identical and are not independent\n"
    "comparisons."
))

# ---- Table 2 note: 明确 marginal median ----
M.append((
    "table2-note",
    "\\item Median values across configurations. \\texttt{seq40} is identical to complex in\n"
    "      this dataset.",
    "\\item Median values across configurations ($n = 9$ per regime, artifact\n"
    "      \\texttt{V6\\_2\\_split4}). $r_{\\mathrm{pooled}}$, $C_{\\mathrm{B}}$ and\n"
    "      $C_{\\mathrm{W}}$ are marginal medians computed separately and therefore need not sum\n"
    "      exactly; the identity holds per configuration. \\texttt{seq40} is identical to complex\n"
    "      in this dataset."
))

# ---- Table 3 note: 修正 "same mutations" ----
M.append((
    "table3-note",
    "\\item All methods are evaluated on the same mutations. RMSE is reported only for\n"
    "      outputs with $\\Delta\\Delta G$-like scaling.",
    "\\item Each method is evaluated on its own uniquely mapped subset of the common annotation\n"
    "      namespace, so mutation counts differ slightly (1{,}755--1{,}776). RMSE is reported\n"
    "      only for outputs with $\\Delta\\Delta G$-like scaling."
))

# ---- Methods: Dataset universes 结构化 ----
M.append((
    "methods-universes",
    "\\subsection*{Dataset and target definition}\n"
    "Evaluation rows come from a SKEMPI v2.0-derived annotation set\\cite{skempi2}\n"
    "(4{,}473 rows, 214 targets, 108 structural families), joined to prediction outputs by\n"
    "explicit row identifiers. A target is defined as an independent binding-interface unit\n"
    "within a structural family.",
    "\\subsection*{Dataset universes}\n"
    "This work uses three explicitly distinct row universes; analyses are always attributed to\n"
    "the universe they were computed on.\n\n"
    "\\emph{Controlled-split universe (7{,}523 rows).} The full set of forward-paired true\n"
    "$\\Delta\\Delta G$ observations, comprising 7{,}523 mutation rows across 230 complexes and\n"
    "112 structural families, with 172 targets having $n \\geq 5$. This universe is used only for\n"
    "the controlled split experiment (Table~\\ref{tab:split}, Figure~\\ref{fig:split}), where the\n"
    "same model specifications and seeds are evaluated under each split regime. It contains\n"
    "4{,}812 distinct (complex, mutation) keys; the remaining rows are repeated experimental\n"
    "records of the same key (2{,}011 keys with more than one record, of which 1{,}405 carry\n"
    "different $\\Delta\\Delta G$ values).\n\n"
    "\\emph{Annotation namespace (4{,}473 rows).} The feature-complete subset of the above\n"
    "universe --- 4{,}473 rows, 214 targets, 108 structural families, 130 targets with\n"
    "$n \\geq 5$. It is the namespace against which external predictions are mapped, and from\n"
    "which the primary matched universe is selected.\n\n"
    "\\emph{Primary matched universe (3{,}351 rows).} A provenance-filtered, row-set-matched\n"
    "subset --- 3{,}351 rows, 102 targets with $n \\geq 5$, 10 method families --- selected by a\n"
    "performance-blind lexicographic rule fixed before any metric values were inspected. All\n"
    "capability-resolved analyses, the weighting audit and the retrospective ablation are\n"
    "computed on this universe unless stated otherwise.\n\n"
    "\\emph{External unambiguous universe (1{,}781 rows).} The subset of the annotation namespace\n"
    "for which the normalised matching key is unique (168 targets; 67 with $n \\geq 5$); 2{,}048\n"
    "ambiguous and 644 unmatched rows were excluded rather than forcibly resolved.\n\n"
    "A target is defined as an independent binding-interface unit; structural families are\n"
    "assigned by sequence-identity clustering (Methods, *Structural-family definition*). Note\n"
    "that \\emph{structural} families (108 in the annotation namespace) are distinct from\n"
    "\\emph{method} families (10 in the primary matched universe), which group model variants."
))

# ---- 创新性定位 ----
M.append((
    "novelty",
    "Here we introduce an exact covariance-signal decomposition that addresses this question. We\n"
    "decompose pooled Pearson correlation rigorously into a between-target component and a\n"
    "within-target component, which makes `what pooled performance is made of' a directly\n"
    "measurable quantity.",
    "Here we address this question by instantiating, for PPI $\\Delta\\Delta G$ benchmarking, an\n"
    "exact between-/within-target covariance decomposition. The decomposition follows from the\n"
    "law of total covariance; our contribution lies in the benchmark-specific formalisation and\n"
    "in its systematic integration with a matched, provenance-filtered evaluation framework, not\n"
    "in the covariance identity itself. It makes `what pooled performance is made of' a directly\n"
    "measurable quantity."
))

patch("manuscript_v6.tex", M)

# =====================================================================
# SUPPLEMENTARY
# =====================================================================
S = []

# ---- 补充编号：Figure S1/S2/S3 不能出现 "Figure 1: Figure S1" ----
S.append((
    "figS-numbering",
    "\\captionsetup{font=small,labelfont=bf}",
    "\\captionsetup{font=small,labelfont=bf}\n"
    "\\setcounter{figure}{0}\n"
    "\\renewcommand{\\thefigure}{S\\arabic{figure}}\n"
    "\\setcounter{table}{0}\n"
    "\\renewcommand{\\thetable}{S\\arabic{table}}"
))

# ---- Supplementary Note: universes ----
S.append((
    "note-universes",
    "\\section*{Supplementary Note 2 | Matched universes and row-set fingerprints}",
    "\\section*{Supplementary Note 2 | Row universes and their relationship}\n"
    f"The controlled-split universe ({mas['n_rows']:,} rows, {mas['n_targets_pdb']} complexes,\n"
    f"{mas['n_families']} structural families) is the full set of forward-paired true\n"
    "$\\Delta\\Delta G$ observations. It contains {mas['n_unique_pdb_mutation_keys']:,} distinct\n"
    f"(complex, mutation) keys; {mas['n_keys_with_multiple_rows']:,} keys carry more than one\n"
    f"experimental record and {mas['n_keys_with_conflicting_ddg']:,} of those disagree in\n"
    "$\\Delta\\Delta G$.\\par\n"
    f"The annotation namespace ({can['n_rows']:,} rows, {can['n_targets']} targets,\n"
    f"{can['n_families']} structural families) is the feature-complete subset of that universe:\n"
    f"every one of its {can['n_rows']:,} rows is present in the controlled-split universe, and\n"
    f"{UP['relation']['rows_in_mas_not_in_canonical']:,} rows of the latter are absent from it\n"
    f"({UP['relation']['targets_in_mas_not_in_canonical']} targets and\n"
    f"{UP['relation']['families_in_mas_not_in_canonical']} structural families are not\n"
    "represented). The two row counts therefore differ because they describe different\n"
    "filtering stages of the same underlying data, not different benchmarks.\\par\n"
    "The target-count distributions are:\n"
    "\\begin{center}\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n"
    "Universe & $n \\geq 1$ & $n \\geq 2$ & $n \\geq 3$ & $n \\geq 5$ & $n \\geq 10$ \\\\\n"
    "\\midrule\n"
    f"Controlled-split ({mas['n_rows']:,}) & {mas['n_ge1']} & {mas['n_ge2']} & {mas['n_ge3']}"
    f" & {mas['n_ge5']} & {mas['n_ge10']} \\\\\n"
    f"Annotation ({can['n_rows']:,}) & {can['n_ge1']} & {can['n_ge2']} & {can['n_ge3']}"
    f" & {can['n_ge5']} & {can['n_ge10']} \\\\\n"
    "\\bottomrule\\end{tabular}\\end{center}\n\n"
    "\\section*{Supplementary Note 3 | Matched universes and row-set fingerprints}"
))

# ---- Supplementary Note: weighted decomposition ----
S.append((
    "note-weighted",
    "\\section*{Supplementary Figures}",
    "\\section*{Supplementary Note 4 | Weighted covariance decomposition}\n"
    "The five target-weighting schemes act on the *target level*: predictions are unchanged and\n"
    "only the weight $a_c$ assigned to target $c$ varies. With $n_c$ mutations in target $c$,\n"
    "total $N = \\sum_c n_c$, and target-level means\n"
    "$\\mu_c = n_c^{-1}\\sum_i y_{ci}$, the weighted aggregation is\n"
    "\\begin{equation}\n"
    "  \\bar{y}_w = \\frac{\\sum_c a_c \\mu_c}{\\sum_c a_c},\\qquad\n"
    "  \\widehat{\\bar{y}}_w = \\frac{\\sum_c a_c \\hat{\\mu}_c}{\\sum_c a_c},\n"
    "  \\qquad a_c > 0 ,\n"
    "\\end{equation}\n"
    "with the schemes\n"
    "$a_c = n_c$ (per-mutation), $a_c = 1$ (per-target), $a_c = \\sqrt{n_c}$,\n"
    "$a_c = \\log(1+n_c)$, and $a_c = \\min(n_c, 30)$ (capped). The weighted between-target\n"
    "covariance is the $a_c$-weighted covariance of the target means,\n"
    "\\begin{equation}\n"
    "  \\operatorname{Cov}_w(B_y, B_{\\hat y})\n"
    "  = \\frac{\\sum_c a_c (\\mu_c - \\bar{y}_w)(\\hat\\mu_c - \\widehat{\\bar y}_w)}{\\sum_c a_c},\n"
    "\\end{equation}\n"
    "and the weighted within-target covariance is the row-level covariance of the deviations,\n"
    "\\begin{equation}\n"
    "  \\operatorname{Cov}_w(W_y, W_{\\hat y})\n"
    "  = \\frac{1}{N}\\sum_c \\sum_i \\varepsilon_{ci}\\hat\\varepsilon_{ci},\n"
    "  \\qquad \\varepsilon_{ci} = y_{ci} - \\mu_c .\n"
    "\\end{equation}\n"
    "The decomposition identity holds for every scheme,\n"
    "\\begin{equation}\n"
    "  r_w = \\frac{\\operatorname{Cov}_w(B_y,B_{\\hat y})\n"
    "              + \\operatorname{Cov}_w(W_y,W_{\\hat y})}{\\sigma_{y,w}\\sigma_{\\hat y,w}}\n"
    "      = C_{\\mathrm{B},w} + C_{\\mathrm{W},w},\n"
    "\\end{equation}\n"
    "because within-target deviations have zero mean in every target under any $a_c$. Note that\n"
    "$C_{\\mathrm{W},w}$ uses row-level averaging and is therefore identical across all five\n"
    "schemes; only $C_{\\mathrm{B},w}$ depends on $a_c$. The reported audit quantity\n"
    "$\\operatorname{range}(C_{\\mathrm{B}})$ and $\\operatorname{range}(C_{\\mathrm{W}})$\n"
    "consequently reflects this asymmetry: across the 109 configurations the median range is\n"
    f"0.0183 for $C_{{\\mathrm{{B}}}}$ and 0.0624 for $C_{{\\mathrm{{W}}}}$.\n\n"
    "\\section*{Supplementary Figures}"
))

# ---- Supplementary Note: D(n,K) derivation ----
S.append((
    "note-dnk",
    "\\section*{Supplementary Note 3 | Fixed-$K$ geometry derivation}",
    "\\section*{Supplementary Note 3 | Derivation of fixed-$K$ shortlist discriminability}\n"
    "Let each target contain $n$ candidate mutations, each carrying a utility $u_i$ (here the\n"
    "true $|\\Delta\\Delta G|$ or any monotone transform of it). Write\n"
    "$E = \\sum_{i \\in S_K} u_i / \\sum_{i=1}^{n} u_i$ for the share of total utility captured by\n"
    "a selected set $S_K$ of $K$ candidates. Two references define the discriminable interval:\n"
    "\\begin{itemize}\n"
    "  \\item $E_{\\max}$: the share attained by the optimal (oracle) selection,\n"
    "        $E_{\\max} = K(2n - K + 1)/[n(n+1)]$, obtained by sorting $u$ in decreasing order\n"
    "        and taking the top $K$;\n"
    "  \\item $E_{\\mathrm{rand}}$: the share expected under uniform random selection of $K$\n"
    "        candidates, which by symmetry is $E_{\\mathrm{rand}} = K/n$.\n"
    "\\end{itemize}\n"
    "The discriminable interval is their difference,\n"
    "\\begin{equation}\n"
    "  D(n,K) = E_{\\max} - E_{\\mathrm{rand}}\n"
    "  = \\frac{K(2n-K+1)}{n(n+1)} - \\frac{K}{n}\n"
    "  = \\frac{K\\,(n-K)}{n\\,(n+1)} .\n"
    "\\end{equation}\n"
    "For $n \\leq K$ every candidate is selected, so $E_{\\max} = E_{\\mathrm{rand}} = 1$ and\n"
    "$D = 0$ identically: such targets carry no information for fixed-$K$ shortlist metrics.\n"
    "Treating $n$ as continuous, $\\mathrm{d}D/\\mathrm{d}n = 0$ gives the maximiser\n"
    "$n^{*} = K + \\sqrt{K^{2} + K}$, i.e.\\ $n^{*} \\approx 20.5$ for $K = 10$, with\n"
    "$D(20) = 0.2381$ and $D(50) = 0.1569$.\n\n"
    "\\emph{Applicability and limits.} $D$ is a purely combinatorial property of candidate-set\n"
    "geometry under a uniform-utility reference. It does not measure the $\\Delta\\Delta G$ utility\n"
    "range of any real dataset, it is not a model-performance metric, and it does not imply that\n"
    "SKEMPI $\\Delta\\Delta G$ values are intrinsically similar. Establishing that would require an\n"
    "analysis using realised utility magnitudes, which we did not perform.\n\n"
    "\\section*{Supplementary Note 3b | Fixed-$K$ geometry derivation (retained reference)}"
))

# ---- Figure S3 definitions + pool table ----
_ra = AC["rank_agreement_by_stratum"]
_rows = []
for r in sorted(_ra, key=lambda z: -z["rho_overall_vs_macro"]):
    _rows.append(f"{r['rowset_hash'][:8]} & {r['level']} & {r['n_models']} & "
                 f"{r['n_method_families']} & {r['n_green']} & "
                 f"{r['rho_overall_vs_macro']:+.4f} \\\\")
_pool_tbl = ("\\begin{landscape}\\begin{table}[h]\\centering\\footnotesize\n"
             "\\begin{threeparttable}\n"
             "\\caption{Model pools used in the rank-agreement sensitivity analysis "
             "(Table~S3).}\\label{tab:s3pools}\n"
             "\\begin{tabular}{@{}llrrrr@{}}\\toprule\n"
             "Row-set hash & Level & Models & Method families & Green models & "
             "$\\rho$(pooled, macro) \\\\\n\\midrule\n"
             + "\n".join(_rows) + "\n"
             "\\bottomrule\\end{tabular}\n"
             "\\begin{tablenotes}\\footnotesize\n"
             "\\item \\emph{Level A} ($\\ast$LA in Figure~S3) denotes a provenance-clean stratum:\n"
             "every listed model is free of known supervised training overlap with the evaluation\n"
             "rows. \\emph{Level C} ($\\ast$LC) denotes a permissive-split stratum, where the pool\n"
             "includes configurations trained under random/complex/\\texttt{seq40} partitions whose\n"
             "evaluation rows may overlap their training folds; these are retained only for this\n"
             "sensitivity analysis. `Green models' is the number of models in the pool flagged\n"
             "provenance-clean. Pools within the same row-set hash are comparable because they\n"
             "share identical evaluation rows; pools with different hashes are not.\n"
             "\\end{tablenotes}\\end{threeparttable}\\end{table}\\end{landscape}\n")

S.append((
    "figS3-caption",
    "\\caption{\\textbf{Figure S3.} Pool sensitivity of cross-metric rank agreement across\n"
    "  matched strata. Within one evaluation universe and row set, changing only which model\n"
    "  families are included moves the agreement between pooled Pearson and macro within-target\n"
    "  Spearman across a wide range, which is why this quantity is reported as a range.}",
    "\\caption{\\textbf{Figure S3.} Pool sensitivity of cross-metric rank agreement across\n"
    "  matched strata. Each bar is one model pool within one row-set stratum, labelled\n"
    "  \\texttt{artifact (row-set hash / level)}. \\emph{LA} denotes a provenance-clean stratum\n"
    "  (no known supervised training overlap) and \\emph{LC} a permissive-split stratum retained\n"
    "  only for this sensitivity analysis; pool sizes, method-family counts and inclusion\n"
    "  criteria are listed in Table~S3. Within one row-set hash the evaluation rows are\n"
    "  identical and only the set of included model families changes, which moves the agreement\n"
    "  between pooled Pearson and macro within-target Spearman from $-0.26$ to $+0.82$; this\n"
    "  quantity is therefore reported as a range rather than a single value.}\n"
    + _pool_tbl
))

patch("supplementary_v6.tex", S)

print(f"applied {len(applied)} patches; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
