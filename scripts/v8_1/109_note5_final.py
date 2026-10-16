# -*- coding: utf-8 -*-
"""scripts/v8_1/109_note5_final.py -- 按 TARGET TOTAL MASS 语义重写 Note 5。

最终语义（依据代码 docstring 的意图 + 数学自洽性）：
  A_c = target c 的总权重,  A_c in {n_c, 1, sqrt(n_c), log1p(n_c), min(n_c,30)}
  w_ci = A_c / n_c            # 行权重
  => target 总质量 = sum_i w_ci = A_c   （与方案名称一致）
  => EQUAL(A=1) 真正是每 target 等权；ROW(A=n_c) 等价于每行等权

核对口径（原实现的字面读法，报告为 sensitivity）：
  w_ci = a_c  => target 总质量 = n_c a_c
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
W = json.loads(Path("v8_1/WEIGHTING_FINAL.json").read_text(encoding="utf-8"))
T = W["stats"]["T_target_total_mass"]
P = W["stats"]["P_per_row_multiplier"]
sfB, sfW = T["signflip_C_B"], T["signflip_C_W"]

NOTE = r"""\section*{Supplementary Note 5 | Weighted covariance decomposition}
The five target-weighting schemes act on the \emph{target level}: predictions are unchanged and
only the weight assigned to each target varies. Let $n_c$ be the number of mutations in target
$c$. Because the schemes are named after target-aggregation rules, we take the target weight
$A_c$ as the quantity they describe,
\begin{equation}
  A_c \in \{\, n_c,\; 1,\; \sqrt{n_c},\; \log(1+n_c),\; \min(n_c, 30) \,\}
  \label{eq:Ac}
\end{equation}
for the per-mutation, per-target, $\sqrt{n}$, logarithmic and capped schemes respectively, and
split it evenly across the target's rows:
\begin{equation}
  w_{ci} = \frac{A_c}{n_c}, \qquad
  \sum_{i} w_{ci} = A_c, \qquad
  A = \sum_{c}\sum_{i} w_{ci} = \sum_{c} A_c .
  \label{eq:wci}
\end{equation}
This is the only reading under which the scheme names are accurate: with $A_c = 1$ each target
carries total weight $1$ regardless of its size (genuinely per-target), and with $A_c = n_c$
every row carries weight $1$ (genuinely per-mutation). Writing $w_{ci} = a_c$ instead --- the
literal form present in the original implementation --- would give each target total weight
$n_c a_c$, so that the ``per-target'' scheme would still let large targets dominate; we
therefore treat that literal reading as a sensitivity variant
(Eq.~\eqref{eq:sens} below).

Under the measure in Eq.~\eqref{eq:wci} the weighted means and variances are
\begin{equation}
  \bar{y}_w = \frac{1}{A}\sum_{c}\sum_{i} w_{ci}\, y_{ci}, \qquad
  \sigma^{2}_{y,w} = \frac{1}{A}\sum_{c}\sum_{i} w_{ci}\,(y_{ci} - \bar{y}_w)^{2},
  \label{eq:means}
\end{equation}
and likewise for $\hat{y}$. Because $w_{ci} = A_c/n_c$ is constant within a target, the
within-target mean is the ordinary target mean,
\begin{equation}
  B_{y,c} \;=\; \frac{\sum_i w_{ci}\, y_{ci}}{\sum_i w_{ci}} \;=\; \bar{y}_c .
  \label{eq:byc}
\end{equation}
With $\varepsilon_{ci} = y_{ci} - B_{y,c}$ and
$\hat{\varepsilon}_{ci} = \hat{y}_{ci} - B_{\hat{y},c}$, the components are
\begin{align}
  C_{\mathrm{B},w} &= \frac{\sum_{c} A_c\,
      \bigl(B_{y,c} - \bar{y}_w\bigr)\bigl(B_{\hat{y},c} - \bar{\hat{y}}_w\bigr)}
      {A\,\sigma_{y,w}\,\sigma_{\hat{y},w}}, \label{eq:cbw}\\[4pt]
  C_{\mathrm{W},w} &= \frac{\sum_{c}\sum_{i} w_{ci}\,
      \varepsilon_{ci}\,\hat{\varepsilon}_{ci}}
      {A\,\sigma_{y,w}\,\sigma_{\hat{y},w}}. \label{eq:cww}
\end{align}
Substituting $y_{ci} = B_{y,c} + \varepsilon_{ci}$ and
$\hat{y}_{ci} = B_{\hat{y},c} + \hat{\varepsilon}_{ci}$ into
$A\,\sigma_{y,w}\sigma_{\hat{y},w}\,r_w = \operatorname{Cov}_w(y,\hat{y})$ gives
\begin{equation}
  \sum_{c}\sum_{i} w_{ci}\,(y_{ci} - \bar{y}_w)(\hat{y}_{ci} - \bar{\hat{y}}_w)
  = \underbrace{\sum_{c} A_c
      \bigl(B_{y,c} - \bar{y}_w\bigr)\bigl(B_{\hat{y},c} - \bar{\hat{y}}_w\bigr)}_{
      \text{Eq.~\eqref{eq:cbw}}}
  + \underbrace{\sum_{c}\sum_{i} w_{ci}\,\varepsilon_{ci}\hat{\varepsilon}_{ci}}_{
      \text{Eq.~\eqref{eq:cww}}},
  \label{eq:cross}
\end{equation}
because for each target $\sum_i w_{ci}\varepsilon_{ci} = 0$ and
$\sum_i w_{ci}\hat{\varepsilon}_{ci} = 0$, so both cross terms vanish. Hence
\begin{equation}
  r_w \;=\; C_{\mathrm{B},w} + C_{\mathrm{W},w}
  \label{eq:identity-w}
\end{equation}
holds exactly for every scheme and every evaluation instance. Verified numerically on all 109
configurations crossed with all five schemes (545 instances): maximum absolute residual
$5.6 \times 10^{-17}$.

\emph{Sensitivity to the weighting reading.} If instead the row weight is taken literally as
$w_{ci} = a_c$ with $a_c$ given by Eq.~\eqref{eq:Ac}, the target total weight becomes
$n_c a_c$ and the same algebra applies with
\begin{equation}
  w^{\mathrm{sens}}_{ci} = a_c, \qquad \sum_i w^{\mathrm{sens}}_{ci} = n_c\,a_c .
  \label{eq:sens}
\end{equation}
The identity still holds exactly, because numerator and denominator share one measure. The two
readings differ only in \emph{what a scheme equalises}, and therefore in the size of the
reported components:

\begin{center}
\begin{tabular}{@{}lrr@{}}
\toprule
Quantity & Eq.~\eqref{eq:wci} (primary) & Eq.~\eqref{eq:sens} (sensitivity) \\
\midrule
Median cross-scheme range of $C_{\mathrm{B}}$ & __RCB__ & __PRCB__ \\
Median cross-scheme range of $C_{\mathrm{W}}$ & __RCW__ & __PRCW__ \\
$C_{\mathrm{B}}$ sign reversals (any magnitude) & __SFB0__/109 & __PSFB0__/109 \\
$C_{\mathrm{B}}$ sign reversals with $|C_{\mathrm{B}}|>0.005$ & __SFB5__/109 & __PSFB5__/109 \\
$C_{\mathrm{W}}$ sign reversals (any magnitude) & __SFW0__/109 & __PSFW0__/109 \\
\bottomrule
\end{tabular}
\end{center}

Both readings agree on the qualitative conclusion: the between-target component is markedly
less stable across weighting conventions than the within-target component (a larger
cross-scheme range and many more sign reversals), so a reported decomposition must state its
weighting scheme. They disagree on whether $C_{\mathrm{W}}$ ever reverses sign
(__SFW0__/109 under the primary reading versus __PSFW0__/109 under the sensitivity reading),
and on the absolute magnitude of the ranges; the main text quotes the primary reading and
reports the sensitivity here.
"""

repl = {
    "__RCB__": f"{T['range_C_B_median']:.4f}",
    "__RCW__": f"{T['range_C_W_median']:.4f}",
    "__PRCB__": f"{P['range_C_B_median']:.4f}",
    "__PRCW__": f"{P['range_C_W_median']:.4f}",
    "__SFB0__": str(sfB["0"]), "__SFB5__": str(sfB["0.005"]),
    "__SFW0__": str(sfW["0"]),
    "__PSFB0__": str(P["signflip_C_B"]["0"]),
    "__PSFB5__": str(P["signflip_C_B"]["0.005"]),
    "__PSFW0__": str(P["signflip_C_W"]["0"]),
}
for k, v in repl.items():
    NOTE = NOTE.replace(k, v)

p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
import re  # noqa: E402
m = re.search(r"\\section\*\{Supplementary Note 5 \| Weighted covariance decomposition\}.*?"
              r"(?=\\section\*\{|\Z)", s, re.S)
if not m:
    print("NOTE 5 ANCHOR NOT FOUND"); sys.exit(1)
s = s[:m.start()] + NOTE + "\n" + s[m.end():]
p.write_text(s, encoding="utf-8")
print("Note 5 rewritten with TARGET TOTAL MASS as primary")
print(f"  primary  range_C_B={T['range_C_B_median']:.4f} range_C_W={T['range_C_W_median']:.4f}"
      f"  C_B flips={sfB['0']}/{sfB['0.005']}  C_W flips={sfW['0']}")
print(f"  sensitiv range_C_B={P['range_C_B_median']:.4f} range_C_W={P['range_C_W_median']:.4f}"
      f"  C_B flips={P['signflip_C_B']['0']}/{P['signflip_C_B']['0.005']}"
      f"  C_W flips={P['signflip_C_W']['0']}")
