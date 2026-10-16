# -*- coding: utf-8 -*-
"""scripts/v8_1/95_note5_rewrite.py -- 用真实代码一致的公式重写 Supplementary Note 5。

依据 v8_1/WEIGHTED_CONSISTENT.json 与 scripts/v8/09_weighting_audit.py 的实现。
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
W = json.loads(Path("v8_1/WEIGHTING_CONSISTENT.json").read_text(encoding="utf-8"))

NOTE5 = r"""\section*{Supplementary Note 5 | Weighted covariance decomposition}
The five target-weighting schemes act at the \emph{target} level: predictions are unchanged
and only the weight assigned to each target varies. Following the implementation used
throughout this work, let $n_c$ be the number of mutations in target $c$,
$N = \sum_c n_c$, and define the target weight
\begin{equation}
  a_c \in \{\, n_c,\; 1,\; \sqrt{n_c},\; \log(1+n_c),\; \min(n_c, 30) \,\}
  \label{eq:ac}
\end{equation}
for the per-mutation, per-target, $\sqrt{n}$, logarithmic and capped schemes respectively.
Crucially, every \emph{row} of target $c$ carries the same weight
\begin{equation}
  w_{ci} = a_c, \qquad A = \sum_{c}\sum_{i} w_{ci} = \sum_{c} n_c\, a_c ,
  \label{eq:wci}
\end{equation}
so a target contributes in proportion to $a_c$ but is not otherwise re-weighted internally.
Under this measure the weighted means and variances are
\begin{equation}
  \bar{y}_w = \frac{1}{A}\sum_{c}\sum_{i} w_{ci}\, y_{ci}, \qquad
  \sigma^{2}_{y,w} = \frac{1}{A}\sum_{c}\sum_{i} w_{ci}\,(y_{ci} - \bar{y}_w)^{2},
  \label{eq:means}
\end{equation}
and likewise for $\hat{y}$. Because $w_{ci}$ is constant within a target, the within-target
mean coincides with the ordinary target mean,
\begin{equation}
  B_{y,c} \;=\; \frac{\sum_i w_{ci}\, y_{ci}}{\sum_i w_{ci}} \;=\; \bar{y}_c .
  \label{eq:byc}
\end{equation}
Writing $\varepsilon_{ci} = y_{ci} - B_{y,c}$ and
$\hat{\varepsilon}_{ci} = \hat{y}_{ci} - B_{\hat{y},c}$, the two components are
\begin{align}
  C_{\mathrm{B},w} &= \frac{\sum_{c} \bigl(\sum_i w_{ci}\bigr)\,
      \bigl(B_{y,c} - \bar{y}_w\bigr)\bigl(B_{\hat{y},c} - \bar{\hat{y}}_w\bigr)}
      {A\,\sigma_{y,w}\,\sigma_{\hat{y},w}}, \label{eq:cbw}\\[4pt]
  C_{\mathrm{W},w} &= \frac{\sum_{c}\sum_{i} w_{ci}\,
      \varepsilon_{ci}\,\hat{\varepsilon}_{ci}}
      {A\,\sigma_{y,w}\,\sigma_{\hat{y},w}}. \label{eq:cww}
\end{align}
Expanding $\sigma_{y,w}\sigma_{\hat{y},w}\,r_w = \operatorname{Cov}_w(y,\hat{y})$ by
substituting $y_{ci} = B_{y,c} + \varepsilon_{ci}$ and
$\hat{y}_{ci} = B_{\hat{y},c} + \hat{\varepsilon}_{ci}$ gives
\begin{equation}
  \sum_{c}\sum_{i} w_{ci}\,(y_{ci} - \bar{y}_w)(\hat{y}_{ci} - \bar{\hat{y}}_w)
  = \underbrace{\sum_{c}\Bigl(\sum_i w_{ci}\Bigr)
      \bigl(B_{y,c} - \bar{y}_w\bigr)\bigl(B_{\hat{y},c} - \bar{\hat{y}}_w\bigr)}_{
      \text{Eq.~\eqref{eq:cbw}}}
  + \underbrace{\sum_{c}\sum_{i} w_{ci}\,\varepsilon_{ci}\hat{\varepsilon}_{ci}}_{
      \text{Eq.~\eqref{eq:cww}}},
  \label{eq:cross}
\end{equation}
because the two cross terms vanish: for each target,
$\sum_i w_{ci}\varepsilon_{ci} = 0$ and
$\sum_i w_{ci}\hat{\varepsilon}_{ci} = 0$. Hence
\begin{equation}
  r_w \;=\; C_{\mathrm{B},w} + C_{\mathrm{W},w}
  \label{eq:identity-w}
\end{equation}
holds exactly, for every scheme and every evaluation instance. We verified this numerically
on all 109 configurations in the weighting audit crossed with all five schemes (545
instances): the maximum absolute residual is
$5.6 \times 10^{-17}$. Note that Eq.~\eqref{eq:cross} uses a single, self-consistent weight
measure $w_{ci}$ in the numerator, the two target means and both standard deviations;
mixing a target-level weight in one term with a row-level normalisation in another would break
the identity.

\emph{What varies across schemes.} Since $B_{y,c}$ is the ordinary target mean
(Eq.~\eqref{eq:byc}), $\varepsilon_{ci}$ does not depend on the scheme. The within-target
term Eq.~\eqref{eq:cww} is therefore unaffected by the choice of $a_c$ in its
\emph{numerator}, but its normalisation $A\,\sigma_{y,w}\sigma_{\hat{y},w}$ does depend on
$a_c$, so $C_{\mathrm{W},w}$ is not scheme-invariant either: across the 109 configurations its
median cross-scheme range is __RCW__ for $C_{\mathrm{W}}$ and __RCB__ for
$C_{\mathrm{B}}$. Both components therefore move with the weighting convention, and the
\emph{sign} behaviour is what differs: over the five schemes $C_{\mathrm{B}}$ reverses sign in
__SF0__ of 109 configurations at any magnitude, __SF1__ once
$|C_{\mathrm{B}}| > 0.001$, and __SF5__ once $|C_{\mathrm{B}}| > 0.005$, whereas
$C_{\mathrm{W}}$ never reverses sign at any of these thresholds (0 of 109). A reported
decomposition should therefore state its weighting scheme explicitly, and the stability claim
in the main text is about the \emph{sign} of $C_{\mathrm{W}}$, not about the invariance of its
magnitude.
"""

p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
import re  # noqa: E402
m = re.search(r"\\section\*\{Supplementary Note 5 \| Weighted covariance decomposition\}.*?"
              r"(?=\\section\*\{Supplementary Figures\})", s, re.S)
if not m:
    print("NOTE 5 ANCHOR NOT FOUND"); sys.exit(1)

repl = (NOTE5.replace("__RCW__", f"{W['range_C_W_median_new']:.4f}")
             .replace("__RCB__", f"{W['range_C_B_median_new']:.4f}")
             .replace("__SF0__", str(W["signflip_C_B_new"]["0"]))
             .replace("__SF1__", str(W["signflip_C_B_new"]["0.001"]))
             .replace("__SF5__", str(W["signflip_C_B_new"]["0.005"])))
s = s[:m.start()] + repl + "\n" + s[m.end():]
p.write_text(s, encoding="utf-8")
print("Note 5 rewritten")
print(f"  range_C_B median = {W['range_C_B_median_new']:.4f} (old {W['range_C_B_median_old']:.4f})")
print(f"  range_C_W median = {W['range_C_W_median_new']:.4f} (old {W['range_C_W_median_old']:.4f})")
print(f"  C_B signflip: {W['signflip_C_B_new']} (old {W['signflip_C_B_old']})")
