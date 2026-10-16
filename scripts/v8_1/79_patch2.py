# -*- coding: utf-8 -*-
"""scripts/v8_1/79_patch2.py -- 补丁第二轮：处理 78 中因空白/顺序未命中的 5 处。

策略：用正则做「空白容忍」匹配，逐条独立应用；已应用的条目不会重复命中。
同样只做审计要求的修改，数值来自 v8_1/ 的真实 artifact。
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

done, miss = [], []


def sub_regex(path, tag, pattern, repl, flags=0):
    p = OUT / path
    s = p.read_text(encoding="utf-8")
    # 用 lambda 返回字面量，避免 re 把替换文本里的 \c 之类当转义解析
    new, n = re.subn(pattern, lambda m: repl, s, count=1, flags=flags)
    if n:
        p.write_text(new, encoding="utf-8"); done.append(tag)
    else:
        miss.append(tag)


# ---------------- Figure 1 caption：整体重写 A/B/C 三段 ----------------
FIG1_NEW = r"""  \caption{\textbf{Capability-resolved evaluation of PPI $\Delta\Delta G$ predictors.}
  (A) Provenance-filtered, fingerprint-matched benchmark construction: the SKEMPI v2.0
  annotation namespace (4{,}473 rows, 214 targets, 108 \emph{structural} families --- distinct
  from the 10 \emph{method} families) is audited by SHA-256 row-set fingerprinting;
  predictions whose recorded training split is permissive, or whose evaluation rows may
  overlap their training folds, are excluded, giving a primary matched universe of 3{,}351 rows,
  102 targets with $n \geq 5$ and 10 method families. The annotation namespace is itself a
  feature-coverage subset of the 7{,}523-row controlled-split universe (Methods,
  \emph{Dataset universes}). (B) Exact additive covariance-signal decomposition. The identity
  $r_{\mathrm{pooled}} = C_{\mathrm{B}} + C_{\mathrm{W}}$ holds for every individual
  configuration (maximum residual $4.4 \times 10^{-16}$); the two components are common-scale
  additive contributions, not independent correlations or causal mechanisms. (C) Controlled
  split collapse on the identical 7{,}523-row universe with identical model specifications and
  seeds; only the split regime changes. Pooled correlation falls from 0.6844 to 0.1508, and
  81.1\% of the summed absolute component change is carried by the between-target component.
  (D) Retrospective matched ablation on the identical 3{,}351-row evaluation universe. Pooled
  Pearson and within-target Macro Spearman are shown for the mutation-only baseline (F1),
  static local geometry (F2), mutation $\times$ partner perturbation (F4), and the full
  local-physics model (F8). Static local geometry improves pooled correlation while leaving
  within-target ranking essentially unchanged, whereas F4 and F8 improve both axes. The
  capability trajectory is retrospective and is not intended as a causal or prospective
  model-selection analysis.}"""

sub_regex("manuscript_v6.tex", "fig1-caption",
          r"  \\caption\{\\textbf\{Capability-resolved evaluation.*?model-selection analysis\.\}",
          FIG1_NEW, flags=re.S)

# ---------------- Figure 2 caption ----------------
FIG2_NEW = r"""  \caption{\textbf{Family-out-of-distribution evaluation disproportionately removes the
  between-target component.} (A) Pooled Pearson by split regime. (B) Median component
  summaries; the exact identity $r_{\mathrm{pooled}} = C_{\mathrm{B}} + C_{\mathrm{W}}$ holds
  per configuration, but the plotted medians are marginal statistics computed separately and
  therefore need not sum exactly. $C_{\mathrm{W}}$ remains positive. (C) Attribution of the
  summed absolute component change, of which the between-target component accounts for 81.1\%;
  the paired attribution is 81.0\%. All four regimes use the same 7{,}523 rows with identical
  model specifications and seeds; complex and \texttt{seq40} are identical and are not
  independent comparisons.}"""

sub_regex("manuscript_v6.tex", "fig2-caption",
          r"  \\caption\{\\textbf\{Family-out-of-distribution evaluation.*?independent comparisons\.\}",
          FIG2_NEW, flags=re.S)

# ---------------- Supplementary Note 3 (D(n,K)) ----------------
DNK_NEW = r"""\section*{Supplementary Note 3 | Derivation of fixed-$K$ shortlist discriminability}
Let each target contain $n$ candidate mutations carrying utilities $u_i$ (here the true
$|\Delta\Delta G|$ or any monotone transform of it). Write
$E = \sum_{i \in S_K} u_i / \sum_{i=1}^{n} u_i$ for the share of total utility captured by a
selected set $S_K$ of $K$ candidates. Two references define the discriminable interval:
\begin{itemize}
  \item $E_{\max}$: the share attained by the optimal (oracle) selection,
        $E_{\max} = K(2n - K + 1)/[n(n+1)]$, obtained by sorting $u$ in decreasing order and
        taking the top $K$;
  \item $E_{\mathrm{rand}}$: the share expected under uniform random selection of $K$
        candidates, which by symmetry is $E_{\mathrm{rand}} = K/n$.
\end{itemize}
The discriminable interval is their difference,
\begin{equation}
  D(n,K) = E_{\max} - E_{\mathrm{rand}}
  = \frac{K(2n-K+1)}{n(n+1)} - \frac{K}{n}
  = \frac{K\,(n-K)}{n\,(n+1)} .
\end{equation}
For $n \leq K$ every candidate is selected, so $E_{\max} = E_{\mathrm{rand}} = 1$ and $D = 0$
identically: such targets carry no information for fixed-$K$ shortlist metrics. Treating $n$ as
continuous, $\mathrm{d}D/\mathrm{d}n = 0$ gives the maximiser
$n^{*} = K + \sqrt{K^{2} + K}$, i.e.\ $n^{*} \approx 20.5$ for $K = 10$, with
$D(20) = 0.2381$ and $D(50) = 0.1569$.

\emph{Applicability and limits.} $D$ is a purely combinatorial property of candidate-set
geometry under a uniform-utility reference. It does not measure the $\Delta\Delta G$ utility
range of any real dataset, it is not a model-performance metric, and it does not imply that
SKEMPI $\Delta\Delta G$ values are intrinsically similar. Establishing that would require an
analysis using realised utility magnitudes, which we did not perform."""

sub_regex("supplementary_v6.tex", "note-dnk",
          r"\\section\*\{Supplementary Note 3 \| Fixed-\$K\$ geometry derivation\}.*?"
          r"(?=\\begin\{table\}\[h\]|\\section\*\{Supplementary Note 4)",
          DNK_NEW, flags=re.S)

print(f"applied {len(done)}: {done}")
if miss:
    print(f"MISSED {len(miss)}: {miss}")
