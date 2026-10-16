# -*- coding: utf-8 -*-
"""scripts/v8_1/80_patch3.py -- 第三轮：补齐 D(n,K) 推导、重排补充 Note 编号、加交叉引用。

数值来源：v8_1/FROZEN_RESULTS.json (axis_c_geometry) 与 v8_1/AUDIT_CONSOLIDATED.json。
"""
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")
A = Path("v8_1")
FR = json.loads((A / "FROZEN_RESULTS.json").read_text(encoding="utf-8"))
G = FR["axis_c_geometry"]
AC = json.loads((A / "AUDIT_CONSOLIDATED.json").read_text(encoding="utf-8"))

done, miss = [], []


def sub(path, tag, pattern, repl, flags=0):
    p = OUT / path
    s = p.read_text(encoding="utf-8")
    new, n = re.subn(pattern, lambda m: repl, s, count=1, flags=flags)
    if n:
        p.write_text(new, encoding="utf-8"); done.append(tag)
    else:
        miss.append(tag)


DNK = r"""\section*{Supplementary Note 4 | Derivation of fixed-$K$ shortlist discriminability}
Let a target contain $n$ candidate mutations carrying utilities $u_i$ (here the true
$|\Delta\Delta G|$ or any monotone transform of it), and let
$E = \sum_{i \in S_K} u_i / \sum_{i=1}^{n} u_i$ denote the share of the total utility captured
by a selected set $S_K$ of $K$ candidates. Two references define the discriminable interval:
\begin{itemize}
  \item $E_{\max}$ --- the share attained by the optimal (oracle) selection. Sorting $u$ in
        decreasing order and taking the top $K$ gives
        \begin{equation}
          E_{\max} = \frac{K\,(2n-K+1)}{n\,(n+1)} ;
          \label{eq:emax}
        \end{equation}
  \item $E_{\mathrm{rand}}$ --- the share expected under uniform random selection of $K$ of the
        $n$ candidates. By symmetry of the uniform distribution over subsets,
        \begin{equation}
          E_{\mathrm{rand}} = \frac{K}{n} .
          \label{eq:erand}
        \end{equation}
\end{itemize}
The discriminable interval is their difference, which simplifies in closed form:
\begin{equation}
  D(n,K) \;=\; E_{\max} - E_{\mathrm{rand}}
  \;=\; \frac{K(2n-K+1)}{n(n+1)} - \frac{K}{n}
  \;=\; \frac{K\,(n-K)}{n\,(n+1)} .
  \label{eq:dnk}
\end{equation}
For $n \leq K$ every candidate is selected, hence $E_{\max} = E_{\mathrm{rand}} = 1$ and
$D = 0$ identically: such targets carry no information for fixed-$K$ shortlist metrics.
Treating $n$ as continuous and differentiating Eq.~\eqref{eq:dnk},
\begin{equation}
  \frac{\mathrm{d}D}{\mathrm{d}n}
  = K\,\frac{(n+1) - (n-K)(2n+1)}{n^{2}(n+1)^{2}} = 0
  \;\Longrightarrow\; n - K = \frac{n(n+1)}{2n+1},
\end{equation}
whose positive root is
\begin{equation}
  n^{*} = K + \sqrt{K^{2} + K} .
  \label{eq:nstar2}
\end{equation}
For $K = 10$ this gives $n^{*} \approx 20.5$, with $D(20) = 0.2381$ and
$D(50) = 0.1569$; in the present benchmark $53$ of $130$ targets with $n \geq 5$ ($40.8\%$)
fall in the degenerate regime $n \leq K$.

\emph{Applicability and limits.} $D$ is a purely combinatorial property of candidate-set
geometry under a uniform-utility reference. It is not a model-performance metric, it does not
measure the $\Delta\Delta G$ utility range of any real dataset, and it does not imply that
SKEMPI $\Delta\Delta G$ values are intrinsically similar. Establishing that would require an
analysis using realised utility magnitudes, which we did not perform.

"""

# 1) 在 Note 3（matched universes）之后、Note 4（weighted）之前插入 D(n,K)
sub("supplementary_v6.tex", "insert-dnk",
    r"(\\section\*\{Supplementary Note 4 \| Weighted covariance decomposition\})",
    DNK + r"\section*{Supplementary Note 5 | Weighted covariance decomposition}", flags=0)

# 2) weighted 已在插入时直接命名为 Note 5，无需再重编号

# 3) 正文：Table 3 前加 Note 4 的交叉引用
sub("manuscript_v6.tex", "ref-note-dnk",
    r"which for \$K = 10\$ gives \$n\^\{\*\} \\approx 20\.5\$, so that \$D\(20\) = 0\.2381\$ while\n"
    r"\$D\(50\) = 0\.1569\$",
    r"which for $K = 10$ gives $n^{*} \approx 20.5$, so that $D(20) = 0.2381$ while\n"
    r"$D(50) = 0.1569$ (full derivation, including the definitions of $E_{\max}$ and\n"
    r"$E_{\mathrm{rand}}$, in Supplementary Note~4)")

# 4) 正文：weighting schemes 指向 Note 5
sub("manuscript_v6.tex", "ref-note-weighted",
    r"\\subsection\*\{Target-weighting schemes\}\nPredictions are unchanged; only aggregation "
    r"weights vary: per-mutation, per-target,\nproportional to \$\\sqrt\{n\}\$, proportional to "
    r"\$\\log\(1\+n\)\$, and capped at 30\.",
    r"\\subsection*{Target-weighting schemes}\n"
    r"Predictions are unchanged; only the target-level aggregation weight varies: per-mutation "
    r"($a_c = n_c$), per-target ($a_c = 1$), proportional to $\\sqrt{n_c}$, proportional to "
    r"$\\log(1+n_c)$, and capped at 30 ($a_c = \\min(n_c, 30)$). The corresponding weighted "
    r"decomposition, and the proof that $r_w = C_{\\mathrm{B},w} + C_{\\mathrm{W},w}$ holds for "
    r"every scheme, are given in Supplementary Note~5.")

# 5) Table S3 冲突：池表改名为 Table S4
sub("supplementary_v6.tex", "pooltab-rename",
    r"\(Table~S3\)\.\}\\label\{tab:s3pools\}",
    r"(Table~S4).}\\label{tab:s4pools}")
sub("supplementary_v6.tex", "pooltab-ref",
    r"criteria are listed in Table~S3\.",
    r"criteria are listed in Table~S4.")

print(f"applied {len(done)}: {done}")
if miss:
    print(f"MISSED {len(miss)}: {miss}")
