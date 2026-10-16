# -*- coding: utf-8 -*-
"""scripts/v8_1/53_build_latex_v6.py -- 把 v6（Word）逐字转成 LaTeX。

与 49_build_latex.py 的区别：**内容以 v6 PDF 为准**，不重新组织行文。
正文段落、表注、图注、Methods 小节标题、参考文献顺序全部照 v6。

四项格式要求：
  equation 环境 / 矢量 PDF 图 / booktabs 三线表 / BibTeX
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
EP = pd.read_csv(A / "external_clean_profiles.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]
sf0, sf5 = int(W2["signflip_C_B_gt0"].sum()), int(W2["signflip_C_B_gt0.005"].sum())
sfw = int(W2["signflip_C_W_gt0"].sum())
cw3 = H4["family_OOD"]["C_W"]
OUT = R / "latex_v6"
OUT.mkdir(exist_ok=True)


def esc(s) -> str:
    s = str(s)
    s = s.replace("\\_", "\x00")
    s = s.replace("\\", r"\textbackslash{}")
    for a, b in (("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#"),
                 ("$", r"\$"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s.replace("\x00", r"\_")


def row(cells):
    return " & ".join(esc(c) for c in cells) + r" \\"


# ---- Table 1 (v6 order / wording) ----
T1_ROWS = [
    ("Full annotation namespace", "4,473", "214", "130", "---", "Namespace definition"),
    ("Primary matched set", "3,351", "166", "102", "10",
     "Main capability-resolved analyses"),
    ("Controlled split set", "7,523", "230", "172", "12",
     "Split-strictness experiment"),
    ("External unambiguous set", "1,781", "168", "67", "5",
     "Boundary of generalisation"),
]
t1 = "\n".join(row(r) for r in T1_ROWS)

# ---- Table 2 ----
t2 = "\n".join(row(r) for r in [
    ("Random", "0.6844", "+0.4275", "+0.2591", "0.0707", "1.7678"),
    ("Complex", "0.1931", "+0.0233", "+0.1685", "0.0497", "2.6412"),
    ("seq40", "0.1961", "+0.0239", "+0.1721", "0.0631", "2.6911"),
    ("Family-OOD", "0.1508", "-0.0072", "+0.1580", "0.0377", "2.6952"),
])

# ---- Table 3 (v6 wording: FoldX (ddG_FoldX)) ----
NICE3 = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred",
         "ddG_FoldX": "FoldX (ddG\\_FoldX)", "ESM-IF": "ESM-IF",
         "proteinMPNN": "ProteinMPNN"}
t3 = "\n".join(row((
    NICE3.get(r.model, r.model), f"{int(r.n_rows):,}", f"{int(r.n_targets_ge5)}",
    f"{r.overall_pearson:.4f}", f"{r.macro_spearman:.4f}",
    f"{r.C_B_row:+.4f}", f"{r.C_W:+.4f}", f"{r.rmse:.4f}"))
    for _, r in EP.iterrows())

# ---- Table S1 (v6 order, full hash in two-column-ish form) ----
S1_ORDER = ["V6_2_split4", "V6_2_1_donly", "V7.1_z", "V7.0_ablation", "V7.0_losses",
            "V7.1_lsr", "V7.2R1_full", "V7.2R2_trackB", "V7.2R3_ppi",
            "V7.2R3_noise", "V7.3_physics", "V7.2_nested"]
_m = UR.set_index("artifact")
tS1 = "\n".join(row((
    _m.loc[a, "universe_id"], a, _m.loc[a, "rowset_hash"],
    f"{int(_m.loc[a,'n_rows'])}", f"{int(_m.loc[a,'n_targets'])}",
    f"{int(_m.loc[a,'n_families'])}", f"{int(_m.loc[a,'n_ge3'])}",
    f"{int(_m.loc[a,'n_ge5'])}", f"{int(_m.loc[a,'n_ge10'])}",
    f"{int(_m.loc[a,'n_models'])}")) for a in S1_ORDER)

TEX = r"""%=====================================================================
%  v6 -> LaTeX conversion
%  "Disentangling Between-Target and Within-Target Signal in PPI ddG Benchmarking"
%  XeLaTeX + fontspec | equation env | vector PDF figures |
%  booktabs three-line tables | BibTeX (unsrtnat, superscript)
%=====================================================================
\documentclass[11pt,a4paper]{article}
\usepackage[a4paper,top=2.4cm,bottom=2.4cm,left=2.3cm,right=2.3cm]{geometry}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{graphicx}
% 等比缩放：宽不超过文本宽，高不超过正文高的 0.82
\newcommand{\figfit}[2][\textwidth]{%
  \includegraphics[width=#1,height=0.82\textheight,keepaspectratio]{#2}}
\usepackage{longtable}
\usepackage{array}
\usepackage{tabularx}
\usepackage{pdflscape}
\usepackage{threeparttable}
\usepackage[super,square,sort&compress]{natbib}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{textcomp}
\usepackage{microtype}
\usepackage{parskip}

% ---- float placement: allow large figures to sit near their citation ----
\renewcommand{\topfraction}{0.92}
\renewcommand{\bottomfraction}{0.75}
\renewcommand{\textfraction}{0.06}
\renewcommand{\floatpagefraction}{0.75}
\setcounter{topnumber}{3}
\setcounter{bottomnumber}{2}
\setcounter{totalnumber}{4}

\graphicspath{{../figures_v5/}{./}}

\title{\bfseries Disentangling Between-Target and Within-Target Signal in
PPI $\Delta\Delta G$ Benchmarking}
\author{Author Name(s)\thanks{Affiliation(s). Corresponding author: [NAME] $\cdot$
Email: [EMAIL] $\cdot$ ORCID: [ORCID\_PLACEHOLDER].}}
\date{}

\begin{document}
\maketitle

\noindent\textbf{Keywords:} protein--protein interaction; binding affinity;
$\Delta\Delta G$; mutation effect prediction; benchmarking; out-of-distribution
generalisation

%=====================================================================
\section*{Abstract}
%=====================================================================
Predicting how mutations alter protein--protein binding affinity ($\Delta\Delta G$)
underpins antibody engineering, interface optimisation and the interpretation of
disease-associated variants. The pooled correlation commonly reported for such benchmarks
is difficult to interpret directly, because it depends both on whether a model places
different targets on the correct overall binding scale and on whether it orders mutations
correctly within each target. We introduce an exact covariance-signal decomposition that
writes pooled Pearson correlation as $r = C_{\mathrm{B}} + C_{\mathrm{W}}$, and apply it to
matched, de-replicated and provenance-filtered evaluation sets. In a controlled comparison
that changes only the split regime (identical 7{,}523 mutation rows and identical model
configurations), median pooled Pearson falls from 0.6844 under random splitting to 0.1508
under family-out-of-distribution evaluation, with
$\Delta C_{\mathrm{B}} = -0.4347$ and $\Delta C_{\mathrm{W}} = -0.1011$; on an
absolute-change basis the between-target component accounts for 81.1\% of the decline while
the within-target component remains positive. Target weighting changes how these
contributions should be read: $C_{\mathrm{B}}$ reverses direction only near zero
(35/109 configurations), whereas $C_{\mathrm{W}}$ keeps its sign at every threshold
although its magnitude is not invariant. Fixed-$K$ shortlist metrics are constrained by
candidate-set geometry: for $K = 10$, $D(n,K)$ is zero when $n \leq K$ and peaks near
$n^{*} \approx 20.5$, and 40.8\% of targets with $n \geq 5$ fall in that degenerate regime.
In a matched retrospective ablation the decomposition separated a change that improved
pooled scores only from one that also improved within-target ranking. External evidence
comprised five methods and is underpowered and positively aligned, so it does not support a
universal decoupling of pooled and within-target performance. We recommend that PPI
$\Delta\Delta G$ benchmarks report both components under an explicitly declared weighting,
together with candidate-set geometry.

%=====================================================================
\section{Introduction}
%=====================================================================
Predicting how mutations alter protein--protein binding affinity ($\Delta\Delta G$) is
central to antibody affinity maturation, interface optimisation and the interpretation of
disease-associated variants. As the number of deep-learning methods grows rapidly, progress
is increasingly judged by benchmark evaluation, and the properties of the evaluation
protocol itself --- what it measures, and under which extrapolation conditions --- determine
what the field can conclude.

Current practice relies largely on pooled metrics, computing Pearson or Spearman correlation
and RMSE over the entire mutation set\cite{biometry}, often supplemented by per-complex
or per-target metrics. Bushuiev et al.\ quantified structural leakage in PPI data: about
65\% of test examples have a near-duplicate training example under PDB-code splitting, about
30\% even under 30\% sequence-identity splitting, 56\% for SKEMPI
v2.0\cite{skempi2} split by PPI code, and 0\% under domain-expert
grouping\cite{leakage}. CATH-ddG adopted CATH homologous-superfamily splits and defined
hard mutations by a maximum TM-score below 0.6\cite{cathddg}; related methods likewise
report per-structure or per-target correlations\cite{ddgpred,prefolddg,ppiminer}. We emphasise that
leakage-aware splitting, family-out-of-distribution evaluation, per-target reporting and
shortlist (Top-$K$) evaluation were none of them introduced by this work.

A key question nevertheless remains unanswered directly: where does pooled benchmark
performance come from? Reporting pooled and per-target metrics side by side does not
quantify the additive contribution of between-target and within-target signal to pooled
performance. In benchmarks containing many targets with highly unequal mutation counts,
pooled correlation is additionally dominated by large targets, which amplifies the
conflation.

Here we introduce an exact covariance-signal decomposition that addresses this question. We
decompose pooled Pearson correlation rigorously into a between-target component and a
within-target component, which makes `what pooled performance is made of' a directly
measurable quantity. We quantify four things: which component is changed by split
strictness; how target weighting changes the interpretation of the two components; which
constraints candidate-set geometry places on fixed-$K$ shortlist evaluation; and whether
the decomposition can distinguish model modifications that improve pooled scores only from
those that also improve within-target ranking. We further bound the evidence using five
independently developed external methods.

%=====================================================================
\section{Results}
%=====================================================================

%% ---- Table 1 ----
\begin{table}[htbp]
  \centering
  \small
  \begin{threeparttable}
  \caption{Dataset and evaluation overview.}
  \label{tab:overview}
  \begin{tabular}{@{}lrrrrl@{}}
    \toprule
    Evaluation set & Rows & Targets & Targets ($n \geq 5$) & Method families & Purpose \\
    \midrule
__T1__
    \bottomrule
  \end{tabular}
  \begin{tablenotes}\footnotesize
    \item All rows are joined by explicit row identifiers.
  \end{tablenotes}
  \end{threeparttable}
\end{table}

\subsection*{Pooled PPI $\Delta\Delta G$ performance contains separable between-target and
within-target signal}

For a mutation $i$ belonging to target $c$, the observed $\Delta\Delta G$ can be written as a
target-level mean plus a mutation-specific deviation from that mean:
\begin{equation}
  y_{ci} = \mu_c + \varepsilon_{ci}
  \qquad\qquad
  \hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}
  \label{eq:decomp}
\end{equation}
With $\sigma_y$ and $\sigma_{\hat{y}}$ the pooled standard deviations, the pooled
correlation splits exactly as
\begin{equation}
  r_{\mathrm{pooled}}
    = \frac{\operatorname{Cov}(y,\hat{y})}{\sigma_y \sigma_{\hat{y}}}
    = C_{\mathrm{B}} + C_{\mathrm{W}},
  \label{eq:identity}
\end{equation}
where $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ are the between-target and within-target
covariance contributions. The identity holds to machine precision (residual
$\leq 3.3 \times 10^{-16}$). These are covariance-signal contributions on a common scale,
not independent correlations, skill scores or causal mechanisms.

We apply this decomposition on a matched evaluation set comprising 3{,}351 mutation rows,
102 targets with at least five mutations, and 10 independent method families.

%% ---- Figure 1 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig1_overview_EN.pdf}
  \caption{\textbf{Capability-resolved evaluation of PPI $\Delta\Delta G$ predictors.}
  (A) Provenance-filtered, fingerprint-matched benchmark construction: SKEMPI v2.0
  annotation (4{,}473 rows, 214 targets, 108 families) is audited by SHA-256 row-set
  fingerprinting; family-out-of-distribution predictions without training overlap are
  retained and permissively split predictions excluded, giving a primary matched universe of
  3{,}351 rows, 102 targets with $n \geq 5$ and 10 method families. Target cards show
  between-target differences in mean $\Delta\Delta G$ ($\mu_c$, purple) and within-target
  mutation ordering ($\varepsilon_{ci}$, green). (B) Exact additive covariance-signal
  decomposition; $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ sum to pooled correlation by an
  exact identity (residual $\leq 3.3 \times 10^{-16}$) and are common-scale additive
  contributions, not independent correlations or causal mechanisms. (C) Controlled split
  collapse on identical 7{,}523 mutation rows: pooled correlation falls from 0.6844 to
  0.1508 and 81.1\% of the absolute decline is carried by loss of the between-target
  component. (D) Retrospective matched ablation on the identical 3{,}351-row evaluation
  universe. Pooled Pearson and within-target Macro Spearman are shown for the mutation-only
  baseline (F1), static local geometry (F2), mutation $\times$ partner perturbation (F4),
  and the full local-physics model (F8). Static local geometry improves pooled correlation
  while leaving within-target ranking essentially unchanged, whereas F4 and F8 improve both
  axes. The capability trajectory is retrospective and is not intended as a causal or
  prospective model-selection analysis.}
  \label{fig:overview}
\end{figure}

\subsection*{Family-OOD evaluation disproportionately removes between-target signal}

To isolate the effect of split strictness, we compared identical model configurations on the
same 7{,}523 mutation rows while changing only the split regime. Median pooled Pearson falls
from 0.6844 to 0.1508. The decline is not uniform:
$\Delta C_{\mathrm{B}} = -0.4347$ and $\Delta C_{\mathrm{W}} = -0.1011$. On an
absolute-change basis the between-target component accounts for 81.1\% of the decline, while
$C_{\mathrm{W}}$ remains positive ($+0.1580$). This is not the same proposition as `the model
memorised target identity': a more defensible reading is that the pattern is consistent with
reliance on target-specific information that does not transfer across families. The
complex-level and \texttt{seq40} partitions produced identical assignments and were
therefore not treated as independent comparisons.

%% ---- Figure 2 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig2_controlled_split_collapse_EN.pdf}
  \caption{\textbf{Family-out-of-distribution evaluation disproportionately removes the
  between-target component.} (A) Pooled Pearson by split regime. (B) Exact decomposition
  into $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$; $C_{\mathrm{W}}$ remains positive.
  (C) Attribution of the absolute decline, of which the between-target component accounts
  for 81.1\%. All four regimes use the same 7{,}523 rows and configurations; complex and
  \texttt{seq40} are identical and are not independent comparisons.}
  \label{fig:split}
\end{figure}

%% ---- Table 2 ----
\begin{table}[htbp]
  \centering
  \small
  \begin{threeparttable}
  \caption{Controlled split comparison.}
  \label{tab:split}
  \begin{tabular}{@{}lrrrrr@{}}
    \toprule
    Split regime & Pooled Pearson & $C_{\mathrm{B}}$ & $C_{\mathrm{W}}$
      & Macro Spearman & RMSE \\
    \midrule
__T2__
    \bottomrule
  \end{tabular}
  \begin{tablenotes}\footnotesize
    \item Median values across configurations. \texttt{seq40} is reported separately; complex in
      this dataset.
  \end{tablenotes}
  \end{threeparttable}
\end{table}

\subsection*{Target aggregation changes the interpretation of covariance contributions}

Holding predictions fixed and varying only the aggregation weights across five schemes
(per-mutation, per-target, proportional to $\sqrt{n}$, proportional to $\log(1+n)$, and
capped at 30), $C_{\mathrm{B}}$ shows directional instability concentrated near zero:
35/109 configurations reverse sign with no threshold, whereas none do once
$\lvert C_{\mathrm{B}} \rvert$ exceeds 0.005 (0/109). $C_{\mathrm{W}}$ never reverses sign at
any threshold (0/109), but its magnitude is not invariant --- its median range is 0.0624,
larger than the 0.0183 of $C_{\mathrm{B}}$. Any reported decomposition should declare its
weighting scheme.

%% ---- Figure 3 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig3_target_weighting_sensitivity_EN.pdf}
  \caption{\textbf{Target weighting changes the interpretation of covariance
  contributions.} Predictions are unchanged; only aggregation weights differ. (A, B)
  Distribution of the range of $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ across the five
  schemes, with the median marked. (C) Configurations reversing sign at three thresholds.
  $C_{\mathrm{B}}$ reverses direction only near zero, whereas $C_{\mathrm{W}}$ keeps its
  sign but its magnitude is not invariant.}
  \label{fig:weighting}
\end{figure}

\subsection*{Fixed-$K$ shortlist evaluation is constrained by candidate-set geometry}

Fixed-$K$ top-$K$ metrics are used for mutation library or shortlist
screening\cite{ndcg}. The discriminable interval is
\begin{equation}
  D(n,K) = 0 \quad (n \leq K),
  \qquad
  D(n,K) = \frac{K\,(n-K)}{n\,(n+1)} \quad (n > K).
  \label{eq:geometry}
\end{equation}
When a target has no more than $K$ mutations, any model selects every candidate and the
target carries no information for fixed-$K$ shortlist metrics. In this benchmark, 53 of 130
targets with $n \geq 5$ (40.8\%) fall in that degenerate regime. $D$ is not monotonic in
$n$: it peaks at
\begin{equation}
  n^{*} = K + \sqrt{K^{2} + K},
  \label{eq:nstar}
\end{equation}
which for $K = 10$ gives $n^{*} \approx 20.5$, so that $D(20) = 0.2381$ while
$D(50) = 0.1569$. This geometric result does not imply that SKEMPI $\Delta\Delta G$ values
are intrinsically too similar; establishing that would require an analysis using actual
$\Delta\Delta G$ utility magnitudes, which we did not perform.

%% ---- Figure 4 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig4_fixed_K_geometry_EN.pdf}
  \caption{\textbf{Fixed-$K$ shortlist evaluation is constrained by candidate-set
  geometry.} The discriminable interval $D(n,K)$ is zero when $n \leq K$ and follows
  $K(n-K)/[n(n+1)]$ otherwise, peaking at $n^{*} = K + \sqrt{K^{2} + K} \approx 20.5$ for
  $K = 10$. In this benchmark 40.8\% of targets with $n \geq 5$ lie in the degenerate regime
  where every model selects all candidates.}
  \label{fig:geometry}
\end{figure}

\subsection*{Capability decomposition clarifies retrospective remedy effects}

In a matched ablation set sharing an identical 3{,}351-row universe, the decomposition
distinguishes modifications that raise pooled scores only from those that also improve
within-target ranking. Relative to a mutation-only baseline, adding static local geometry
raised pooled Pearson by about 0.028 while leaving within-target Macro Spearman essentially
unchanged ($-0.002$); mutation $\times$ partner perturbation improved both axes ($+0.070$
and $+0.044$); and a full local-physics model gave the largest improvement on both axes
($+0.105$ and $+0.076$). This analysis is retrospective and does not establish prospective
prediction of which modification will succeed.

%% ---- Figure 5 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig5_retrospective_remedy_EN.pdf}
  \caption{\textbf{Capability decomposition separates pooled-score gains from within-target
  gains.} Retrospective matched ablation on an identical 3{,}351-row universe; the baseline
  is the mutation-only model. Green markers indicate improvement on both axes; red markers
  indicate a pooled-score-only or negative change. The analysis is retrospective and does
  not imply prospective prediction of which modification will succeed.}
  \label{fig:remedy}
\end{figure}

\subsection*{External methods define the boundary of generalisation}

We built external evaluation sets from published official per-mutation predictions: Flex
ddG\cite{flexddg}, FoldX\cite{foldx5}, ProteinMPNN\cite{proteinmpnn},
ESM-IF\cite{esmif} and DDGPred\cite{ddgpred}. Where the matching key was not
unique on the annotation side, only uniquely mapped rows were used; of 4{,}473 rows only
1{,}781 have a unique mapping, and ambiguous keys were excluded rather than forcibly
resolved. On this unambiguous set the five methods agree in their ordering on pooled
correlation and on within-target ranking, with DDGPred best on both axes; no method-level
rank reversal occurs. Method-level rank agreement is $\rho = +0.90$ (exact permutation
$p = 0.083$). This strongly positive agreement does not support a universal decoupling of
pooled and within-target performance and falls within the heterogeneous range observed
internally. With only five methods the panel is underpowered, so we draw no directional
conclusion from it.

%% ---- Figure 6 ----
\begin{figure}[htbp]
  \centering
  \figfit{Fig6_external_methods_EN.pdf}
  \caption{\textbf{External methods define the boundary of generalisation.} (A) Pooled
  versus within-target performance for five independently developed methods on a matched,
  unambiguous mapping. (B) Method-level capability profile, standardised within each metric
  for visual comparison only; raw values are printed in each cell. Rank agreement is
  positive ($\rho = +0.90$, exact permutation $p = 0.083$) but the panel is underpowered
  with five methods.}
  \label{fig:external}
\end{figure}

%% ---- Table 3 ----
\begin{table}[htbp]
  \centering
  \small
  \begin{threeparttable}
  \caption{External method capability profiles on the unambiguous mapping set.}
  \label{tab:external}
  \begin{tabular}{@{}lrrrrrrr@{}}
    \toprule
    Method & Mutations & Targets ($n \geq 5$) & Pooled Pearson & Macro Spearman
      & $C_{\mathrm{B}}$ & $C_{\mathrm{W}}$ & RMSE \\
    \midrule
__T3__
    \bottomrule
  \end{tabular}
  \begin{tablenotes}\footnotesize
    \item All methods are evaluated on the same mutations. RMSE is reported only for
      outputs with $\Delta\Delta G$-like scaling.
  \end{tablenotes}
  \end{threeparttable}
\end{table}

%=====================================================================
\section{Discussion}
%=====================================================================
This study shows that the pooled correlation conventionally reported as a single number in
PPI $\Delta\Delta G$ benchmarking is composed of two qualitatively different predictive
signals: a between-target component that places different targets on the correct overall
scale, and a within-target component that orders mutation effects correctly inside each
target. Because our decomposition is an exact identity, `what pooled performance is made
of' becomes directly computable.

The relation to prior work should be delimited clearly. Structural leakage in PPI data and
its dependence on splitting strategy were quantified systematically by Bushuiev et al., who
also proposed interface-similarity splitting\cite{leakage}; CATH-ddG adopted CATH
homologous-superfamily splits and used structural similarity to define harder extrapolation
mutations\cite{cathddg}; and per-complex or per-target reporting together with
shortlist-style evaluation already appear in the $\Delta\Delta G$
literature\cite{ddgpred,prefolddg,ppiminer}. Our contribution lies not in any single one of these
ingredients but in the systematic integration of an exact between/within-target covariance
decomposition with four audits --- split strictness, target weighting, matched rank
agreement and candidate-set geometry --- inside one matched, provenance-filtered framework,
together with a retrospective attribution case study.

The interpretation of $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ has clear boundaries. They are
covariance contributions within one identity and share a common scale, so they sum to pooled
correlation; they are not two independently estimated skill scores and should not be read as
causal mechanisms. A model may obtain a substantial pooled correlation largely by placing
targets on approximately the right overall scale while its within-target ordering ability
remains limited --- precisely why pooled metrics can overstate practical screening value.

The weighting-sensitivity audit yields a conclusion of direct practical relevance: the sign
of the between-target component is sensitive to aggregation weights near zero, whereas the
sign of the within-target component is stable. When a reported $C_{\mathrm{B}}$ is close to
zero, the conclusion about whether between-target information is present may depend on the
weighting convention rather than on the model. Any reported decomposition should declare its
weighting scheme.

The candidate-set geometry result has direct implications for Top-$K$ practice. When a
target has no more than $K$ mutations, fixed-$K$ shortlist metrics cannot distinguish models
structurally; roughly four in ten targets with $n \geq 5$ fall in this regime in our
benchmark. Because the discriminable interval first rises and then falls as candidate sets
grow, enlarging a candidate set is not universally effective at improving discrimination. We
recommend reporting $n$ alongside $K$, labelling $n \leq K$ targets separately, or adopting
fraction-based evaluation.

The retrospective ablation illustrates the framework's explanatory value: it distinguishes
a change that raised pooled scores only from a change that also improved within-target
ranking. The retrospective nature of the case study must, however, be stated explicitly: it
does not establish prospective ability to predict which modifications will work.

Finally, external evidence bounds the generalisation of these conclusions. Across five
independently developed methods the ordering of methods agreed between pooled and
within-target performance, a strongly positive agreement that neither supports nor refutes a
universal decoupling of the two axes. The panel is underpowered, so we draw no directional
conclusion from it. Future work should broaden the external panel and test the relationship
between the decomposed components and realised screening yield on independent prospective
datasets.

%=====================================================================
\section{Methods}
%=====================================================================
\subsection*{Dataset and target definition}
Evaluation rows come from a SKEMPI v2.0-derived annotation set\cite{skempi2}
(4{,}473 rows, 214 targets, 108 structural families), joined to prediction outputs by
explicit row identifiers. A target is defined as an independent binding-interface unit
within a structural family.

\subsection*{Prediction sources and provenance criteria}
Prediction outputs were classified by their recorded supervised training split regime.
Configurations trained under a permissive (random / complex / \texttt{seq40}) partition were
excluded from deployment-oriented comparisons because their evaluation rows may overlap
their training folds; they were retained only in the controlled split experiment. Analyses
were performed at artifact, configuration and method-family level, with one
performance-blind representative per independent family.

\subsection*{Matched evaluation sets}
Each prediction output received a SHA-256 fingerprint of its sorted row identifiers;
cross-model comparison was performed only between outputs sharing a fingerprint
(Supplementary Note 2, Supplementary Fig.~S1). The primary matched set was selected by a
rule fixed before any metric values were inspected.

\subsection*{Between-/within-target signal decomposition}
Row-level weights are used and the identity residual is reported as a quality-control
quantity. Full derivation in Supplementary Note 1.

\subsection*{Target-weighting schemes}
Predictions are unchanged; only aggregation weights vary: per-mutation, per-target,
proportional to $\sqrt{n}$, proportional to $\log(1+n)$, and capped at 30.

\subsection*{Ranking and shortlist metrics}
Within-target ordering used per-target Spearman $\rho$ and Kendall $\tau$ correlations,
pairwise accuracy and the within-target covariance contribution; shortlist metrics used
fixed-$K$ and fraction-based hit rate, precision, recall, NDCG and enrichment
factor\cite{ndcg}.

\subsection*{Statistical analysis}
Rank agreement uses Spearman $\rho$ with Kendall $\tau$ as secondary; interval estimates use
target-cluster bootstrap (5{,}000 resamples, seed 42)\cite{bootstrap}; small external
panels use exact permutation tests. All cross-metric comparisons use higher-is-better
orientation\cite{biometry}.

\subsection*{External prediction mapping and ambiguity control}
Where the normalised matching key was not unique on the annotation side, only uniquely
mapped rows were used; ambiguous keys were excluded rather than resolved.

\subsection*{Reproducibility and consistency validation}
All headline numbers are generated from a single machine-readable file and verified by an
automated consistency checker. No new model was trained and no GPU was used.

%=====================================================================
\section{Limitations}
%=====================================================================
Our extrapolation setting is based on target-family partitioning, which does not exclude
the possibility that foundation-model pretraining encountered sequences or structures
similar to the evaluation families; we can control task-specific supervised overlap only.
The model pool is not a random sample of all PPI $\Delta\Delta G$ methods, and only five
methods are externally evaluable, so method-level inference is underpowered.

Shortlist and candidate-selection analyses rest on re-analysis of existing predictions and
include no prospective wet-lab validation; when a target has no more than $K$ mutations the
relevant metrics carry no information. The retrospective ablation likewise does not
establish prospective prediction.

$C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ are covariance-signal contributions rather than
causal mechanisms. We performed no identity-memorisation or target-shuffling intervention
and therefore make no causal claim about the origin of the between-target component.

No independent prospective experimental dataset was available, so no external dataset
validation was performed. This work does not provide a universal model recommender and does
not claim general superiority of any method family.

%=====================================================================
\section*{Data and code availability}
%=====================================================================
All analyses use existing predictions and published external predictions. Frozen input
hashes, the machine-readable source of every headline number and an automated consistency
checker are provided with the Supplementary Information. No new model was trained and no GPU
was used.

%=====================================================================
\section*{Author contributions}
%=====================================================================
Study conception, analysis design, interpretation of results and writing of the manuscript
were performed by the author(s). The use of AI-assisted tools is declared below.

%=====================================================================
\section*{Use of AI-assisted tools}
%=====================================================================
AI-assisted coding and language tools were used to support code drafting and debugging,
generation of plotting scripts, document automation, literature and repository retrieval,
and language editing. All metric definitions, analysis choices, scientific interpretations,
hypothesis decisions and final manuscript content were reviewed and approved by the authors,
who take full responsibility for the work. AI tools were not used to generate experimental
data or to modify upstream model predictions.

%=====================================================================
\section*{Acknowledgements}
%=====================================================================
[FUNDING\_PLACEHOLDER] We thank the developers of the publicly released prediction sets and
datasets used in this study.

% v6 reference order (for .bbl reordering):
%   1 skempi2  2 leakage  3 cathddg  4 ddgpred  5 prefolddg  6 flexddg
%   7 foldx5   8 proteinmpnn  9 esmif  10 ppiminer  11 ndcg  12 bootstrap
%   13 biometry
%=====================================================================
\bibliographystyle{plainnat}
\bibliography{references}

\end{document}
"""

TEX = TEX.replace("__T1__", t1).replace("__T2__", t2).replace("__T3__", t3)
(OUT / "manuscript_v6.tex").write_text(TEX, encoding="utf-8")
print(f"wrote manuscript_v6.tex ({len(TEX)/1024:.1f} KB)")

# ---------------- supplement（v6 只有 Note 1/2 + Table S1 + Fig S1-S3） ----------------
SI = r"""%=====================================================================
%  v6 Supplementary Information -> LaTeX
%=====================================================================
\documentclass[11pt,a4paper]{article}
\usepackage[a4paper,top=2.2cm,bottom=2.2cm,left=2.0cm,right=2.0cm]{geometry}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,longtable,array,tabularx}
\usepackage{graphicx}
%%FIGFIT_DEF%%
%%FIGFIT_DEF2%%
\usepackage{pdflscape}
\usepackage{threeparttable}
\usepackage[super,square,sort&compress]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{microtype}
\usepackage{parskip}

\graphicspath{{../figures_v5/}{./}}

\title{\bfseries Supplementary Information\\[6pt]
\large Disentangling Between-Target and Within-Target Signal in
PPI $\Delta\Delta G$ Benchmarking}
\author{}
\date{}

\begin{document}
\maketitle

\section*{Supplementary Note 1 | Derivation of the covariance decomposition}
Write the true value and prediction for mutation $i$ in target $c$ as
$y_{ci} = \mu_c + \varepsilon_{ci}$ and
$\hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}$. Expanding
$\operatorname{Cov}(y,\hat{y})$ gives a between-target term and a within-target term,
\begin{equation}
  \operatorname{Cov}(y,\hat{y})
  = \frac{1}{N}\sum_{c} n_c \left(\mu_c - \bar{\mu}\right)
    \left(\hat{\mu}_c - \bar{\hat{\mu}}\right)
  + \frac{1}{N}\sum_{c}\sum_{i} \varepsilon_{ci}\,\hat{\varepsilon}_{ci};
\end{equation}
all cross terms vanish identically because within-target deviations have zero mean in every
target. Dividing by $\sigma_y \sigma_{\hat{y}}$ gives $C_{\mathrm{B}}$ and
$C_{\mathrm{W}}$ and establishes
\begin{equation}
  r_{\mathrm{pooled}} = C_{\mathrm{B}} + C_{\mathrm{W}}
\end{equation}
exactly. The identity residual is reported as a quality-control quantity (observed maximum
$3.331 \times 10^{-16}$). Because $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ share the same
denominator they are covariance-signal contributions on a common scale, not independent
correlations.

\section*{Supplementary Note 2 | Matched universes and row-set fingerprints}
Each prediction output received a SHA-256 fingerprint of its sorted row identifiers.
Cross-model comparison was performed only between outputs sharing a fingerprint; raw metrics
were never compared across fingerprints.

\begin{landscape}
\begin{table}[h]
  \centering
  \scriptsize
  \begin{threeparttable}
  \caption{Matched-evaluation strata and row-set fingerprints (Table~S1).}
  \label{tab:s1}
  \begin{tabular}{@{}lllrrrrrrr@{}}
    \toprule
    Universe & Artifact & Row-set fingerprint & Rows ($N$) & Targets & Families
      & Targets ($n \geq 3$) & Targets ($n \geq 5$) & Targets ($n \geq 10$) & Models \\
    \midrule
__TS1__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}
\end{landscape}

\section*{Supplementary Figures}
\begin{figure}[htbp]
  \centering
  \figfit{FigS1_matched_universe_strata_EN.pdf}
  \caption{\textbf{Figure S1.} Matched-evaluation strata grouped by SHA-256 row-set
  fingerprints. Each bar is one prediction output; outputs sharing a fingerprint are
  directly comparable with each other and with no others.}
\end{figure}

\begin{figure}[htbp]
  \centering
  \figfit{FigS2_capability_heatmap_provenance_clean_EN.pdf}
  \caption{\textbf{Figure S2.} Capability profiles of matched, provenance-clean ablation
  variants (seed 42). Rows are controlled ablation variants of one local-physics method
  family, not independent methods. Values are standardised within each metric for visual
  comparison only; raw values are printed in each cell.}
\end{figure}

\begin{figure}[htbp]
  \centering
  \figfit{FigS3_pool_sensitivity_rank_agreement_EN.pdf}
  \caption{\textbf{Figure S3.} Pool sensitivity of cross-metric rank agreement across
  matched strata. Within one evaluation universe and row set, changing only which model
  families are included moves the agreement between pooled Pearson and macro within-target
  Spearman across a wide range, which is why this quantity is reported as a range.}
\end{figure}

\end{document}
"""
SI = SI.replace("%%FIGFIT_DEF%%", "\\newcommand{\\figfit}[2][\\textwidth]{%")
SI = SI.replace("%%FIGFIT_DEF2%%", "  \\includegraphics[width=#1,height=0.82\\textheight,keepaspectratio]{#2}}")
SI = SI.replace("__TS1__", tS1)
(OUT / "supplementary_v6.tex").write_text(SI, encoding="utf-8")
print(f"wrote supplementary_v6.tex ({len(SI)/1024:.1f} KB)")

# references.bib 复用
import shutil
shutil.copy2(R / "latex_v7" / "references.bib", OUT / "references.bib")
print("copied references.bib")

# ---- 终稿审计补丁（provenance / 统计口径 / external / Methods / 编号）----
import subprocess  # noqa: E402
for _name in ("78_audit_patch.py", "79_patch2.py", "80_patch3.py", "81_patch_fix.py", "82_fix_stray_backslash.py", "85_methods_expand.py", "86_tableS2.py", "87_tableS2_refine.py", "109_note5_final.py", "96_blocker2_and_consistency.py", "99_note6_and_final_fixes.py", "100_caption_fix.py", "104_round3_fixes.py", "105_round3_regex.py", "106_round3_last2.py", "117_round4_text_fixes.py", "118_abstract_dereplicated.py", "122_external_units_fix.py", "123_table3_rmse_fix.py", "124_external_semantics_methods.py", "126_round6_cleanup.py", "131_oriented_text_fixes.py", "132_supp_orientation.py", "135_fig6_caption_oriented.py", "136_round8_orientation_fix.py", "137_round8_fix2.py", "138_round8_supp.py", "139_om_symbol.py", "141_supp_stale_cleanup.py", "142_supp_notransform.py", "146_fig1_fullpage.py", "147_seq40_keyfix_propagate.py"):  # 148 已直接改图源，无需入链
    _patch = Path(__file__).resolve().parent / _name
    if not _patch.exists():
        continue
    _r = subprocess.run([sys.executable, str(_patch)], capture_output=True,
                        text=True, encoding="utf-8", errors="replace")
    print(f"[{_name}] " + ((_r.stdout or "").strip().replace("\n", " | ") or "no output"))
    if _r.returncode != 0:
        print(f"AUDIT PATCH {_name} FAILED:\n", (_r.stderr or "")[-2000:])
        raise SystemExit(1)
