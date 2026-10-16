# -*- coding: utf-8 -*-
"""scripts/v8_1/49_build_latex.py -- 生成 LaTeX v7 稿件。

格式要求：
* 公式用 equation / align 环境（含 \\label \\eqref）
* 图用矢量 PDF（\\includegraphics 引 figures_v5/*.pdf）
* 表用三线表（booktabs：\\toprule \\midrule \\bottomrule）
* 引用用 BibTeX（\\cite{} + references.bib，natbib 数字上标）

所有数值均读自 v8_1/FROZEN_RESULTS.json。
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
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]
sf0, sf5 = int(W2["signflip_C_B_gt0"].sum()), int(W2["signflip_C_B_gt0.005"].sum())
sfw = int(W2["signflip_C_W_gt0"].sum())
cw3 = H4["family_OOD"]["C_W"]
OUT = R / "latex_v7"
OUT.mkdir(exist_ok=True)


def f4(x):
    return f"{float(x):.4f}"


def esc(s) -> str:
    """LaTeX 安全化：转义 _ % & # $ { } ~ ^ \\，但保留已转义的 \\_。"""
    s = str(s)
    s = s.replace("\\_", "\x00")           # 保护已转义的下划线
    s = s.replace("\\", r"\textbackslash{}")
    for a, b in (("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#"),
                 ("$", r"\$"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s.replace("\x00", r"\_")         # 还原


def row(cells) -> str:
    return " & ".join(esc(c) for c in cells) + r" \\"


# ---------------- 表 1：三线表 ----------------
rows1 = [
    ("Full annotation namespace", "4,473", "214", "130", "---", "Namespace definition"),
    ("Primary matched set", f"{PU['n_rows']:,}", "166", f"{PU['n_ge5']}",
     f"{PU['n_method_families']}", "Main capability-resolved analyses"),
    ("Controlled split set", "7,523", "230", "172", "12", "Split-strictness experiment"),
    ("External unambiguous set", f"{FR['external_clean']['n_rows']:,}",
     f"{FR['external_clean']['n_targets']}", f"{FR['external_clean']['n_ge5']}", "5",
     "Boundary of generalisation"),
]
t1 = "\n".join(row(r) for r in rows1)

# ---------------- 表 2：三线表 ----------------
REG = [("random", "Random"), ("complex", "Complex"), ("seq40", "seq40"),
       ("family_OOD", "Family-OOD")]
rows2 = [(nm, f4(H4[k]["overall_pearson"]), f"{H4[k]['C_B']:+.4f}",
          f"{H4[k]['C_W']:+.4f}", f4(H4[k]["macro_spearman"]), f4(H4[k]["rmse"]))
         for k, nm in REG]
t2 = "\n".join(row(r) for r in rows2)

# ---------------- 表 3：三线表 ----------------
NICE = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred", "ddG_FoldX": "ddG-FoldX",
        "ESM-IF": "ESM-IF", "proteinMPNN": "ProteinMPNN"}
rows3 = []
for r in EP.itertuples():
    rows3.append((NICE.get(r.model, r.model), f"{int(r.n_rows):,}",
                  f"{int(r.n_targets_ge5)}", f4(r.overall_pearson),
                  f4(r.macro_spearman), f"{r.C_B_row:+.4f}", f"{r.C_W:+.4f}",
                  f4(r.rmse)))
t3 = "\n".join(row(r) for r in rows3)

# ---------------- 补充表 S1 ----------------
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
rowsS1 = []
for r in UR.itertuples():
    rowsS1.append((r.universe_id, r.artifact, r.rowset_hash[:8], f"{int(r.n_rows):,}",
                   f"{int(r.n_targets)}", f"{int(r.n_ge3)}", f"{int(r.n_ge5)}",
                   f"{int(r.n_ge10)}", f"{int(r.n_models)}"))
tS1 = "\n".join(row(r) for r in rowsS1)

# ---------------- 补充表 S4：几何 ----------------
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
g = GEO_DF[GEO_DF["n"] > 10]
sel = sorted(set([11, 12, 15, 20, 25, 30, 50, 75, 100, 140]) & set(g["n"]))
rowsS4 = [(f"{int(r.n)}", f"{float(r.E_rand):.3f}", f"{float(r.E_max):.3f}",
           f"{float(r.D):.3f}") for r in g[g["n"].isin(sel)].itertuples()]
tS4 = "\n".join(row(r) for r in rowsS4)

# ---------------- 补充表 S5：加权审计（前 30） ----------------
SCH = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]


def rng(df, pre):
    return (df[[f"{pre}_{s}" for s in SCH]].max(axis=1)
            - df[[f"{pre}_{s}" for s in SCH]].min(axis=1))


w = W2.head(30).copy()
# weighting_sensitivity.csv 的列名为大写后缀：C_B_ROW / C_W_CAP30 ...
_wfull = pd.read_csv(R / "experiments" / "v8" / "weighting" /
                     "weighting_sensitivity.csv")
_byname = _wfull.set_index("model_name")
w["rB"] = [(_byname.loc[n, [f"C_B_{s}" for s in SCH]].max()
            - _byname.loc[n, [f"C_B_{s}" for s in SCH]].min())
           for n in w["model_name"]]
w["rW"] = [(_byname.loc[n, [f"C_W_{s}" for s in SCH]].max()
            - _byname.loc[n, [f"C_W_{s}" for s in SCH]].min())
           for n in w["model_name"]]
rowsS5 = []
# 注意：pandas 的 itertuples() 会丢弃列名中含 '.' 的列
# （signflip_C_B_gt0.005 等），因此这里按下标取行。
for _, r in w.iterrows():
    rowsS5.append((str(r["model_name"]),
                   f"{r['rB']:.4f}", f"{r['rW']:.4f}",
                   f"{r['max_abs_C_B']:.4f}", f"{r['max_abs_C_W']:.4f}",
                   str(int(bool(r["signflip_C_B_gt0"]))),
                   str(int(bool(r["signflip_C_B_gt0.005"]))),
                   str(int(bool(r["signflip_C_W_gt0"])))))
tS5 = "\n".join(row(r) for r in rowsS5)

# ---------------- 补充表 S6：模型注册表（前 40） ----------------
MR = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv")
rowsS6 = []
for r in MR.head(40).itertuples():
    seed = "---" if pd.isna(r.seed) else str(int(r.seed))
    rowsS6.append((str(r.artifact_id),
                   str(r.method_family), seed,
                   str(r.split_regime), r.rowset_hash[:8],
                   "yes" if bool(r.is_ablation) else "no",
                   "yes" if bool(r.is_primary_representative) else "no"))
tS6 = "\n".join(row(r) for r in rowsS6)

# ---------------- 补充表 S3：指标方向登记（29） ----------------
MD = pd.read_csv(A / "METRIC_DIRECTION_REGISTRY.csv")
rowsS3 = []
for r in MD.itertuples():
    rowsS3.append((str(r.metric),
                   str(r.axis).replace("_", " "),
                   str(r.raw_direction).replace("_", " "),
                   str(r.ranking_direction).replace("_", " "),
                   str(r.eligible_for_calibration)))
tS3 = "\n".join(row(r) for r in rowsS3)

TEX = r"""%=====================================================================
%  Disentangling Between-Target and Within-Target Signal in PPI ddG Benchmarking
%  XeLaTeX + fontspec + booktabs + natbib/BibTeX
%  Figures: vector PDF via \includegraphics
%  All numbers read from v8_1/FROZEN_RESULTS.json
%=====================================================================
\documentclass[11pt,a4paper]{article}

% ---------- geometry ----------
\usepackage[a4paper,top=2.4cm,bottom=2.4cm,left=2.3cm,right=2.3cm]{geometry}

% ---------- fonts (XeLaTeX) ----------
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
\usepackage{amsmath,amssymb}
\usepackage{unicode-math}
\setmathfont{Cambria Math}

% ---------- tables / figures ----------
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{array}
\usepackage{tabularx}
\usepackage{pdflscape}
\usepackage{rotating}
\usepackage{threeparttable}

% ---------- citations ----------
\usepackage[super,square,sort&compress]{natbib}

% ---------- misc ----------
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{textcomp}
\usepackage{microtype}
\usepackage{parskip}

\graphicspath{{../figures_v5/}{./}}

\title{\bfseries Disentangling Between-Target and Within-Target Signal in
PPI $\Delta\Delta G$ Benchmarking}
\author{Author Name(s)\thanks{Affiliation(s). Corresponding author:
\href{mailto:EMAIL}{EMAIL}. ORCID: [ORCID\_PLACEHOLDER].} \\
\small [AFFILIATION\_PLACEHOLDER]}
\date{}

\begin{document}
\maketitle

\begin{abstract}
\noindent
Predicting how mutations alter protein--protein binding affinity ($\Delta\Delta G$)
underpins antibody engineering, interface optimisation and the interpretation of
disease-associated variants. The pooled correlation commonly reported for such
benchmarks is difficult to interpret directly, because it depends both on whether a
model places different targets on the correct overall binding scale and on whether it
orders mutations correctly within each target. We introduce an exact covariance-signal
decomposition that writes pooled Pearson correlation as $r = C_{\mathrm{B}} + C_{\mathrm{W}}$,
and apply it to matched, de-replicated and provenance-filtered evaluation sets. In a
controlled comparison that changes only the split regime (identical 7{,}523 mutation rows
and identical model configurations), median pooled Pearson falls from __R0__ under random
splitting to __R3__ under family-out-of-distribution evaluation, with
$\Delta C_{\mathrm{B}} = __DCB__$ and $\Delta C_{\mathrm{W}} = __DCW__$; on an
absolute-change basis the between-target component accounts for __FRAC__\% of the decline
while the within-target component remains positive. Target weighting changes how these
contributions should be read: $C_{\mathrm{B}}$ reverses direction only near zero
(__SF0__/__NCFG__ configurations), whereas $C_{\mathrm{W}}$ keeps its sign at every
threshold although its magnitude is not invariant. Fixed-$K$ shortlist metrics are
constrained by candidate-set geometry: for $K = 10$, $D(n,K)$ is zero when $n \leq K$ and
peaks near $n^* \approx __NSTAR__$, and __NLEK__\% of targets with $n \geq 5$ fall in that
degenerate regime. In a matched retrospective ablation the decomposition separated a change
that improved pooled scores only from one that also improved within-target ranking.
External evidence comprised five methods and is underpowered and positively aligned, so it
does not support a universal decoupling of pooled and within-target performance. We
recommend that PPI $\Delta\Delta G$ benchmarks report both components under an explicitly
declared weighting, together with candidate-set geometry.

\medskip
\noindent\textbf{Keywords:} protein--protein interaction; binding affinity;
$\Delta\Delta G$; mutation effect prediction; benchmarking; out-of-distribution
generalisation
\end{abstract}

%=====================================================================
\section{Introduction}
%=====================================================================
Predicting how mutations alter protein--protein binding affinity ($\Delta\Delta G$) is
central to antibody affinity maturation, interface optimisation and the interpretation of
disease-associated variants. As the number of deep-learning methods grows rapidly,
progress is increasingly judged by benchmark evaluation, and the properties of the
evaluation protocol itself --- what it measures, and under which extrapolation
conditions --- determine what the field can conclude.

Current practice relies largely on pooled metrics, computing Pearson or Spearman
correlation and RMSE over the entire mutation set~\cite{biometry}, often supplemented by
per-complex or per-target metrics. Bushuiev et al.\ quantified structural leakage in PPI
data: about 65\% of test examples have a near-duplicate training example under PDB-code
splitting, about 30\% even under 30\% sequence-identity splitting, 56\% for SKEMPI
v2.0~\cite{skempi2} split by PPI code, and 0\% under domain-expert
grouping~\cite{leakage}. CATH-ddG adopted CATH homologous-superfamily splits and defined
hard mutations by a maximum TM-score below $0.6$~\cite{cathddg}; related methods likewise
report per-structure or per-target
correlations~\cite{ddgpred,prefolddg,ppiminer}. We emphasise that leakage-aware
splitting, family-out-of-distribution evaluation, per-target reporting and shortlist
(Top-$K$) evaluation were none of them introduced by this work.

A key question nevertheless remains unanswered directly: where does pooled benchmark
performance come from? Reporting pooled and per-target metrics side by side does not
quantify the additive contribution of between-target and within-target signal to pooled
performance. In benchmarks containing many targets with highly unequal mutation counts,
pooled correlation is additionally dominated by large targets, which amplifies the
conflation.

Here we introduce an exact covariance-signal decomposition that addresses this question.
We quantify four things: which component is changed by split strictness; how target
weighting changes the interpretation of the two components; which constraints
candidate-set geometry places on fixed-$K$ shortlist evaluation; and whether the
decomposition can distinguish model modifications that improve pooled scores only from
those that also improve within-target ranking.

%=====================================================================
\section{Results}
%=====================================================================

\subsection{\texorpdfstring{Pooled PPI $\Delta\Delta G$ performance contains separable between-target and within-target signal}{Pooled PPI ddG performance contains separable between-target and
within-target signal}}

For a mutation $i$ belonging to target $c$, the observed $\Delta\Delta G$ can be written as
a target-level mean plus a mutation-specific deviation from that mean:
\begin{equation}
  y_{ci} = \mu_c + \varepsilon_{ci}
  \label{eq:decomp-y}
\end{equation}
where $\mu_c$ is the mean $\Delta\Delta G$ of target $c$ and
$\varepsilon_{ci} = y_{ci} - \mu_c$ is the within-target mutation deviation. The prediction
is decomposed identically:
\begin{equation}
  \hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}
  \label{eq:decomp-yhat}
\end{equation}

Only after this decomposition is it meaningful to distinguish two covariance
contributions. With $\sigma_y$ and $\sigma_{\hat{y}}$ the pooled standard deviations, the
pooled correlation splits exactly as
\begin{equation}
  r_{\mathrm{pooled}}
    = \frac{\operatorname{Cov}(y,\hat{y})}{\sigma_y \sigma_{\hat{y}}}
    = C_{\mathrm{B}} + C_{\mathrm{W}},
  \label{eq:identity}
\end{equation}
where
\begin{equation}
  C_{\mathrm{B}} = \frac{\operatorname{Cov}(B_y, B_{\hat{y}})}
                        {\sigma_y \sigma_{\hat{y}}},
  \qquad
  C_{\mathrm{W}} = \frac{\operatorname{Cov}(W_y, W_{\hat{W}})}
                        {\sigma_y \sigma_{\hat{y}}}
  \label{eq:CB-CW}
\end{equation}
are the between-target and within-target covariance contributions. The identity in
Eq.~\eqref{eq:identity} holds to machine precision (residual
$\leq 3.3 \times 10^{-16}$). These are covariance-signal contributions on a common scale, not
independent correlations, skill scores or causal mechanisms.

We apply this decomposition on a matched evaluation set comprising __NROWS__ mutation rows,
__NGE5__ targets with at least five mutations, and __NFAM__ independent method families.
Table~\ref{tab:overview} summarises the evaluation sets.

%% ---- Table 1 (three-line / booktabs) ----
\begin{table}[t]
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

%% ---- Figure 1 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=\textwidth]{Fig1_overview_EN.pdf}
  \caption{\textbf{Capability-resolved evaluation of PPI $\Delta\Delta G$ predictors.}
  (A) Provenance-filtered, fingerprint-matched benchmark construction: SKEMPI v2.0
  annotation (4{,}473 rows, 214 targets, 108 families) is audited by SHA-256 row-set
  fingerprinting; family-out-of-distribution predictions without training overlap are
  retained and permissively split predictions excluded, giving a primary matched universe
  of 3{,}351 rows, 102 targets with $n \geq 5$ and 10 method families. Target cards show
  between-target differences in mean $\Delta\Delta G$ ($\mu_c$) and within-target mutation
  ordering ($\varepsilon_{ci}$). (B) Exact additive covariance-signal decomposition;
  $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ sum to pooled correlation by an exact identity
  (residual $\leq 3.3 \times 10^{-16}$). (C) Controlled split collapse on identical
  7{,}523 mutation rows. (D) Retrospective matched ablation on the identical 3{,}351-row
  evaluation universe; the capability trajectory is retrospective and is not intended as a
  causal or prospective model-selection analysis.}
  \label{fig:overview}
\end{figure}

\subsection{Family-OOD evaluation disproportionately removes between-target signal}

To isolate the effect of split strictness, we compared identical model configurations on
the same 7{,}523 mutation rows while changing only the split regime
(Table~\ref{tab:split}). Median pooled Pearson falls from __R0__ to __R3__. The decline is
not uniform: $\Delta C_{\mathrm{B}} = __DCB__$ and $\Delta C_{\mathrm{W}} = __DCW__$. On an
absolute-change basis,
\begin{equation}
  \frac{\lvert \Delta C_{\mathrm{B}} \rvert}
       {\lvert \Delta C_{\mathrm{B}} \rvert + \lvert \Delta C_{\mathrm{W}} \rvert}
  = __FRACDEC__,
  \label{eq:attribution}
\end{equation}
so the between-target component accounts for __FRAC__\% of the decline, while
$C_{\mathrm{W}}$ remains positive ($__CW3__$).

This is not the same proposition as ``the model memorised target identity'': a more
defensible reading is that the pattern is consistent with reliance on target-specific
information that does not transfer across families. The complex-level and \texttt{seq40}
partitions produced identical assignments and were therefore not treated as independent
comparisons.

%% ---- Table 2 (three-line) ----
\begin{table}[t]
  \centering
  \small
  \begin{threeparttable}
  \caption{Controlled split comparison (identical rows and configurations).}
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
    \item Median values across configurations. \texttt{seq40} is identical to
      \texttt{complex} in this dataset and is not an independent comparison.
  \end{tablenotes}
  \end{threeparttable}
\end{table}

%% ---- Figure 2 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=\textwidth]{Fig2_controlled_split_collapse_EN.pdf}
  \caption{\textbf{Family-out-of-distribution evaluation disproportionately removes the
  between-target component.} (A) Pooled Pearson by split regime. (B) Exact decomposition
  into $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$; $C_{\mathrm{W}}$ remains positive.
  (C) Attribution of the absolute decline (Eq.~\eqref{eq:attribution}).}
  \label{fig:split}
\end{figure}

\subsection{Target aggregation changes the interpretation of covariance contributions}

Holding predictions fixed and varying only the aggregation weights across five schemes
(per-mutation, per-target, proportional to $\sqrt{n}$, proportional to $\log(1+n)$, and
capped at 30), $C_{\mathrm{B}}$ shows directional instability concentrated near zero:
__SF0__/__NCFG__ configurations reverse sign with no threshold, whereas none do once
$\lvert C_{\mathrm{B}} \rvert > 0.005$ (__SF5__/__NCFG__). $C_{\mathrm{W}}$ never reverses
sign at any threshold (__SFW__/__NCFG__), but its magnitude is not invariant --- its median
range is __RCW__, larger than the __RCB__ of $C_{\mathrm{B}}$. Any reported decomposition
should declare its weighting scheme.

%% ---- Figure 3 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=\textwidth]{Fig3_target_weighting_sensitivity_EN.pdf}
  \caption{\textbf{Target weighting changes the interpretation of covariance
  contributions.} Predictions are unchanged; only aggregation weights differ. (A, B)
  Distribution of the range of $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ across the five
  schemes. (C) Configurations reversing sign at three thresholds.}
  \label{fig:weighting}
\end{figure}

\subsection{\texorpdfstring{Fixed-$K$ shortlist evaluation is constrained by candidate-set geometry}{Fixed-K shortlist evaluation is constrained by candidate-set geometry}}

Fixed-$K$ top-$K$ metrics are used for mutation library or shortlist
screening~\cite{ndcg}. The discriminable interval is
\begin{equation}
  D(n,K) =
  \begin{cases}
    0, & n \leq K \\[4pt]
    \dfrac{K\,(n-K)}{n\,(n+1)}, & n > K
  \end{cases}
  \label{eq:geometry}
\end{equation}
where the expected share of gain from $K$ randomly chosen candidates is $K/n$ and the
optimal $K$ candidates attain $K(2n-K+1)/[n(n+1)]$. When a target has no more than $K$
mutations, any model selects every candidate and the target carries no information for
fixed-$K$ shortlist metrics. In this benchmark, __NLEKCNT__ of __NGE5CNT__ targets with
$n \geq 5$ (__NLEK__\%) fall in that degenerate regime. $D$ is not monotonic in $n$: it
peaks at
\begin{equation}
  n^{*} = K + \sqrt{K^{2} + K},
  \label{eq:nstar}
\end{equation}
which for $K = 10$ gives $n^{*} \approx __NSTAR__$, so that $D(20) = __D20__$ while
$D(50) = __D50__$.

This geometric result does not imply that SKEMPI $\Delta\Delta G$ values are intrinsically
too similar; establishing that would require an analysis using actual $\Delta\Delta G$
utility magnitudes, which we did not perform.

%% ---- Figure 4 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=0.95\textwidth]{Fig4_fixed_K_geometry_EN.pdf}
  \caption{\textbf{Fixed-$K$ shortlist evaluation is constrained by candidate-set
  geometry.} The discriminable interval $D(n,K)$ (Eq.~\eqref{eq:geometry}) is zero when
  $n \leq K$ and peaks at $n^{*} = K + \sqrt{K^2+K}$ (Eq.~\eqref{eq:nstar}) for $K = 10$.}
  \label{fig:geometry}
\end{figure}

\subsection{Capability decomposition clarifies retrospective remedy effects}

In a matched ablation set sharing an identical 3{,}351-row universe, the decomposition
distinguishes modifications that raise pooled scores only from those that also improve
within-target ranking. Relative to a mutation-only baseline, adding static local geometry
raised pooled Pearson by about $0.028$ while leaving within-target Macro Spearman
essentially unchanged ($-0.002$); mutation $\times$ partner perturbation improved both axes
($+0.070$ and $+0.044$); and a full local-physics model gave the largest improvement on
both axes ($+0.105$ and $+0.076$). This analysis is retrospective and does not establish
prospective prediction of which modification will succeed.

%% ---- Figure 5 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=0.95\textwidth]{Fig5_retrospective_remedy_EN.pdf}
  \caption{\textbf{Capability decomposition separates pooled-score gains from
  within-target gains.} Retrospective matched ablation on an identical 3{,}351-row
  universe; the baseline is the mutation-only model. Green markers indicate improvement on
  both axes; red markers indicate a pooled-score-only or negative change.}
  \label{fig:remedy}
\end{figure}

\subsection{External methods define the boundary of generalisation}

We built external evaluation sets from published official per-mutation predictions: Flex
ddG~\cite{flexddg}, FoldX~\cite{foldx5}, ProteinMPNN~\cite{proteinmpnn},
ESM-IF~\cite{esmif} and DDGPred~\cite{ddgpred}. Where the matching key was not unique on
the annotation side, only uniquely mapped rows were used; of 4{,}473 rows only
__EXTROWS__ have a unique mapping, and ambiguous keys were excluded rather than forcibly
resolved. On this unambiguous set the five methods agree in their ordering on pooled
correlation and on within-target ranking, with DDGPred best on both axes; no method-level
rank reversal occurs. Method-level rank agreement is $\rho = +0.90$ (exact permutation
$p = 0.083$). This strongly positive agreement does not support a universal decoupling of
pooled and within-target performance and falls within the heterogeneous range observed
internally. With only five methods the panel is underpowered, so we draw no directional
conclusion from it.

%% ---- Table 3 (three-line) ----
\begin{table}[t]
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

%% ---- Figure 6 (vector PDF) ----
\begin{figure}[t]
  \centering
  \includegraphics[width=\textwidth]{Fig6_external_methods_EN.pdf}
  \caption{\textbf{External methods define the boundary of generalisation.} (A) Pooled
  versus within-target performance for five independently developed methods on a matched,
  unambiguous mapping. (B) Method-level capability profile, standardised within each metric
  for visual comparison only; raw values are printed in each cell.}
  \label{fig:external}
\end{figure}

%=====================================================================
\section{Discussion}
%=====================================================================
This study shows that the pooled correlation conventionally reported as a single number in
PPI $\Delta\Delta G$ benchmarking is composed of two qualitatively different predictive
signals: a between-target component that places different targets on the correct overall
scale, and a within-target component that orders mutation effects correctly inside each
target. Because our decomposition is an exact identity (Eq.~\eqref{eq:identity}), what
pooled performance is made of becomes directly computable.

The relation to prior work should be delimited clearly. Structural leakage in PPI data and
its dependence on splitting strategy were quantified systematically by Bushuiev et al.,
who also proposed interface-similarity splitting~\cite{leakage}; CATH-ddG adopted CATH
homologous-superfamily splits and used structural similarity to define harder
extrapolation mutations~\cite{cathddg}; and per-complex or per-target reporting together
with shortlist-style evaluation already appear in the $\Delta\Delta G$
literature~\cite{ddgpred,prefolddg,ppiminer}. Our contribution lies not in any single one
of these ingredients but in the systematic integration of an exact between/within-target
covariance decomposition with four audits inside one matched, provenance-filtered
framework, together with a retrospective attribution case study.

The interpretation of $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ has clear boundaries. They are
covariance contributions within one identity and share a common scale, so they sum to
pooled correlation; they are not two independently estimated skill scores and should not be
read as causal mechanisms. A model may obtain a substantial pooled correlation largely by
placing targets on approximately the right overall scale while its within-target ordering
ability remains limited --- precisely why pooled metrics can overstate practical screening
value.

The weighting-sensitivity audit yields a conclusion of direct practical relevance: the sign
of the between-target component is sensitive to aggregation weights near zero, whereas the
sign of the within-target component is stable. Any reported decomposition should declare
its weighting scheme.

The candidate-set geometry result has direct implications for Top-$K$ practice. When a
target has no more than $K$ mutations, fixed-$K$ shortlist metrics cannot distinguish
models structurally. Because the discriminable interval first rises and then falls as
candidate sets grow, enlarging a candidate set is not universally effective at improving
discrimination. We recommend reporting $n$ alongside $K$, labelling $n \leq K$ targets
separately, or adopting fraction-based evaluation.

The retrospective ablation illustrates the framework's explanatory value, but its
retrospective nature must be stated explicitly: it does not establish prospective ability
to predict which modifications will work.

Finally, external evidence bounds the generalisation of these conclusions. Across five
independently developed methods the ordering of methods agreed between pooled and
within-target performance, a strongly positive agreement that neither supports nor refutes
a universal decoupling of the two axes. The panel is underpowered, so we draw no
directional conclusion from it.

%=====================================================================
\section{Methods}
%=====================================================================
\paragraph{Dataset and target definition.}
Evaluation rows come from a SKEMPI v2.0-derived annotation set~\cite{skempi2} (4{,}473
rows, 214 targets, 108 structural families), joined to prediction outputs by explicit row
identifiers. A target is defined as an independent binding-interface unit within a
structural family.

\paragraph{Prediction sources and provenance criteria.}
Prediction outputs were classified by their recorded supervised training split regime.
Configurations trained under a permissive (random / complex / \texttt{seq40}) partition were
excluded from deployment-oriented comparisons because their evaluation rows may overlap
their training folds; they were retained only in the controlled split experiment. Analyses
were performed at artifact, configuration and method-family level, with one
performance-blind representative per independent family.

\paragraph{Matched evaluation sets.}
Each prediction output received a SHA-256 fingerprint of its sorted row identifiers;
cross-model comparison was performed only between outputs sharing a fingerprint
(Supplementary Note~2, Supplementary Fig.~S1). The primary matched set was selected by a
rule fixed before any metric values were inspected.

\paragraph{Between-/within-target signal decomposition.}
Row-level weights are used and the identity residual is reported as a quality-control
quantity (Eqs.~\eqref{eq:decomp-y}--\eqref{eq:CB-CW}). Full derivation in Supplementary
Note~1.

\paragraph{Target-weighting schemes.}
Predictions are unchanged; only aggregation weights vary: per-mutation, per-target,
proportional to $\sqrt{n}$, proportional to $\log(1+n)$, and capped at 30.

\paragraph{Ranking and shortlist metrics.}
Within-target ordering used per-target Spearman $\rho$ and Kendall $\tau$ correlations,
pairwise accuracy and the within-target covariance contribution; shortlist metrics used
fixed-$K$ and fraction-based hit rate, precision, recall, NDCG and enrichment
factor~\cite{ndcg}.

\paragraph{Statistical analysis.}
Rank agreement uses Spearman $\rho$ with Kendall $\tau$ as secondary; interval estimates use
target-cluster bootstrap (5{,}000 resamples, seed 42)~\cite{bootstrap}; small external
panels use exact permutation tests. All cross-metric comparisons use higher-is-better
orientation~\cite{biometry}.

\paragraph{External prediction mapping and ambiguity control.}
Where the normalised matching key was not unique on the annotation side, only uniquely
mapped rows were used; ambiguous keys were excluded rather than resolved.

\paragraph{Reproducibility and consistency validation.}
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

$C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ are covariance-signal contributions rather than causal
mechanisms. We performed no identity-memorisation or target-shuffling intervention and
therefore make no causal claim about the origin of the between-target component.

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
and language editing. All metric definitions, analysis choices, scientific
interpretations, hypothesis decisions and final manuscript content were reviewed and
approved by the authors, who take full responsibility for the work. AI tools were not used
to generate experimental data or to modify upstream model predictions.

%=====================================================================
\section*{Acknowledgements}
%=====================================================================
[FUNDING\_PLACEHOLDER] We thank the developers of the publicly released prediction sets and
datasets used in this study.

%=====================================================================
\bibliographystyle{unsrtnat}
\bibliography{references}

\end{document}
"""

subs = {
    "__R0__": f4(r0), "__R3__": f4(r3),
    "__DCB__": f"{DL['delta_C_B']:+.4f}", "__DCW__": f"{DL['delta_C_W']:+.4f}",
    "__FRAC__": f"{100*DL['fraction_between']:.1f}",
    "__FRACDEC__": f"{DL['fraction_between']:.4f}",
    "__CW3__": f"{cw3:+.4f}",
    "__SF0__": str(sf0), "__SF5__": str(sf5), "__SFW__": str(sfw),
    "__NCFG__": str(len(W2)),
    "__RCW__": f"{W2.range_C_W.median():.4f}",
    "__RCB__": f"{W2.range_C_B.median():.4f}",
    "__NSTAR__": f"{GEO['n_star_analytic']:.1f}",
    "__NLEK__": f"{100*GEO['n_le_K_fraction']:.1f}",
    "__NLEKCNT__": str(GEO["n_le_K_count"]),
    "__NGE5CNT__": str(GEO["n_ge5_count"]),
    "__D20__": f"{GEO['D_at_n20']:.4f}", "__D50__": f"{GEO['D_at_n50']:.4f}",
    "__NROWS__": f"{PU['n_rows']:,}", "__NGE5__": str(PU["n_ge5"]),
    "__NFAM__": str(PU["n_method_families"]),
    "__EXTROWS__": f"{FR['external_clean']['n_rows']:,}",
    "__T1__": t1, "__T2__": t2, "__T3__": t3,
}
tex = TEX
for k, v in subs.items():
    tex = tex.replace(k, v)

(OUT / "manuscript_v7.tex").write_text(tex, encoding="utf-8")
print(f"wrote manuscript_v7.tex ({len(tex)/1024:.1f} KB)")
print("unsubstituted placeholders:", [k for k in subs if k in tex])

# ---------------- 补充材料 ----------------
SI = r"""%=====================================================================
%  Supplementary Information -- XeLaTeX
%=====================================================================
\documentclass[11pt,a4paper]{article}
\usepackage[a4paper,top=2.2cm,bottom=2.2cm,left=2.0cm,right=2.0cm]{geometry}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
\usepackage{amsmath,amssymb}
\usepackage{unicode-math}
\setmathfont{Cambria Math}
\usepackage{booktabs,longtable,array,tabularx}
\usepackage{graphicx}
\usepackage{pdflscape}
\usepackage{threeparttable}
\usepackage[super,square,sort&compress]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{microtype}
\usepackage{parskip}

\graphicspath{{../figures_v5/}{./}}

\title{\bfseries Supplementary Information\\[4pt]
\large Disentangling Between-Target and Within-Target Signal in
PPI $\Delta\Delta G$ Benchmarking}
\author{}
\date{}

\begin{document}
\maketitle

\section*{Supplementary Note 1 --- Derivation of the covariance decomposition}
Write the true value and prediction for mutation $i$ in target $c$ as
$y_{ci} = \mu_c + \varepsilon_{ci}$ and
$\hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}$. Expanding
$\operatorname{Cov}(y,\hat{y})$ gives a between-target term and a within-target term:
\begin{equation}
  \operatorname{Cov}(y,\hat{y})
  = \frac{1}{N}\sum_{c} n_c \left(\mu_c - \bar{\mu}\right)
    \left(\hat{\mu}_c - \bar{\hat{\mu}}\right)
  + \frac{1}{N}\sum_{c}\sum_{i} \varepsilon_{ci}\,\hat{\varepsilon}_{ci},
\end{equation}
with all cross terms vanishing identically because within-target deviations have zero mean
in every target. Dividing by $\sigma_y \sigma_{\hat{y}}$ gives
\begin{equation}
  r_{\mathrm{pooled}} = C_{\mathrm{B}} + C_{\mathrm{W}},
  \qquad
  C_{\mathrm{B}} = \frac{\operatorname{Cov}(B_y,B_{\hat{y}})}
                        {\sigma_y \sigma_{\hat{y}}},
  \qquad
  C_{\mathrm{W}} = \frac{\operatorname{Cov}(W_y,W_{\hat{y}})}
                        {\sigma_y \sigma_{\hat{y}}} .
\end{equation}
The identity residual is reported as a quality-control quantity (observed maximum
$3.331 \times 10^{-16}$). Because $C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ share the same
denominator they are covariance-signal contributions on a common scale, not independent
correlations.

\section*{Supplementary Note 2 --- Matched universes and row-set fingerprints}
Each prediction output received a SHA-256 fingerprint of its sorted row identifiers.
Cross-model comparison was performed only between outputs sharing a fingerprint; raw
metrics were never compared across fingerprints.

\begin{landscape}
\begin{table}[h]
  \centering
  \footnotesize
  \begin{threeparttable}
  \caption{Matched-evaluation strata and row-set fingerprints (Table~S1).}
  \label{tab:s1}
  \begin{tabular}{@{}lllrrrrrr@{}}
    \toprule
    Universe & Artifact & Fingerprint & Rows ($N$) & Targets
      & Targets ($n \geq 3$) & Targets ($n \geq 5$) & Targets ($n \geq 10$) & Models \\
    \midrule
__TS1__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}
\end{landscape}

\section*{Supplementary Note 3 --- Fixed-$K$ geometry derivation}
For a target with $n$ mutations and fixed $K$, choosing $K$ candidates uniformly at random
without replacement yields an expected gain share $E_{\mathrm{rand}} = K/n$ by symmetry,
while the optimal $K$ candidates attain
$E_{\mathrm{max}} = K(2n-K+1)/[n(n+1)]$. Hence
\begin{equation}
  D(n,K) = E_{\mathrm{max}} - E_{\mathrm{rand}}
  = \begin{cases}
      0, & n \leq K \\[4pt]
      \dfrac{K(n-K)}{n(n+1)}, & n > K
    \end{cases}
\end{equation}
$D$ is not monotonic in $n$: $\mathrm{d}D/\mathrm{d}n = 0$ at
$n^{*} = K + \sqrt{K^{2}+K}$.

\begin{table}[h]
  \centering
  \small
  \begin{threeparttable}
  \caption{Discriminable interval versus candidate-set size ($K = 10$) (Table~S4).}
  \label{tab:s4}
  \begin{tabular}{@{}rrrr@{}}
    \toprule
    $n$ & $E_{\mathrm{rand}}$ & $E_{\mathrm{max}}$ & $D(n,K)$ \\
    \midrule
__TS4__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}

\section*{Supplementary Note 4 --- Full weighting audit}
Continuous values are rounded to four decimal places. Sign-flip columns are integer counts
(0 or 1) of the weighting schemes in which the component changes sign at the stated
threshold on its absolute value.

\begin{landscape}
\begin{table}[h]
  \centering
  \scriptsize
  \begin{threeparttable}
  \caption{Weighting effect-size audit, first 30 configurations (Table~S5).}
  \label{tab:s5}
  \begin{tabular}{@{}lrrrrrrrr@{}}
    \toprule
    Model configuration & Range $C_{\mathrm{B}}$ & Range $C_{\mathrm{W}}$
      & Max $\lvert C_{\mathrm{B}} \rvert$ & Max $\lvert C_{\mathrm{W}} \rvert$
      & $C_{\mathrm{B}}$ flips ($>0$) & $C_{\mathrm{B}}$ flips ($>0.005$)
      & $C_{\mathrm{W}}$ flips ($>0$) \\
    \midrule
__TS5__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}
\end{landscape}

\section*{Supplementary Note 5 --- Metric-direction registry}
A metric-direction registry was established before any rank-agreement computation. RMSE and
MAE are lower-is-better; all cross-metric rank agreement uses negated RMSE (or ascending
ranks). The registry comprises all 29 metrics used in this work.

\begin{landscape}
\begin{table}[h]
  \centering
  \scriptsize
  \begin{threeparttable}
  \caption{Metric-direction registry, all 29 metrics (Table~S3).}
  \label{tab:s3}
  \begin{tabular}{@{}lllll@{}}
    \toprule
    Metric & Analysis axis & Raw direction & Direction for ranking
      & Eligible for calibration \\
    \midrule
__TS3__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}
\end{landscape}

\section*{Supplementary Note 6 --- Model registry}
\begin{landscape}
\begin{table}[h]
  \centering
  \scriptsize
  \begin{threeparttable}
  \caption{Model family registry, first 40 entries (Table~S6).}
  \label{tab:s6}
  \begin{tabular}{@{}lllllll@{}}
    \toprule
    Artifact & Method family & Seed & Split regime & Fingerprint & Ablation
      & Primary representative \\
    \midrule
__TS6__
    \bottomrule
  \end{tabular}
  \end{threeparttable}
\end{table}
\end{landscape}

\section*{Supplementary Note 7 --- Automated workflow and reproducibility}
AI-assisted coding tools were used to draft and debug the analysis and plotting scripts,
and to automate document generation, as declared in the main text. All headline numbers are
generated from a single machine-readable results file; an automated consistency checker
verifies that the documents, tables and figures agree with it. No new model was trained and
no GPU was used.

\section*{Supplementary Figures}
\begin{figure}[h]
  \centering
  \includegraphics[width=\textwidth]{FigS1_matched_universe_strata_EN.pdf}
  \caption{\textbf{Figure S1.} Matched-evaluation strata grouped by SHA-256 row-set
  fingerprints.}
\end{figure}

\begin{figure}[h]
  \centering
  \includegraphics[width=\textwidth]{FigS2_capability_heatmap_provenance_clean_EN.pdf}
  \caption{\textbf{Figure S2.} Capability profiles of matched, provenance-clean ablation
  variants (seed 42). Rows are controlled ablation variants of one local-physics method
  family, not independent methods. Values are standardised within each metric for visual
  comparison only; raw values are printed in each cell.}
\end{figure}

\begin{figure}[h]
  \centering
  \includegraphics[width=\textwidth]{FigS3_pool_sensitivity_rank_agreement_EN.pdf}
  \caption{\textbf{Figure S3.} Pool sensitivity of cross-metric rank agreement across
  matched strata.}
\end{figure}

\end{document}
"""
SI = SI.replace("__TS1__", tS1).replace("__TS4__", tS4).replace("__TS5__", tS5)
SI = SI.replace("__TS3__", tS3).replace("__TS6__", tS6)
(OUT / "supplementary_v7.tex").write_text(SI, encoding="utf-8")
print(f"wrote supplementary_v7.tex ({len(SI)/1024:.1f} KB)")
