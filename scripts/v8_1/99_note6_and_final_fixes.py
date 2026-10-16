# -*- coding: utf-8 -*-
"""scripts/v8_1/99_note6_and_final_fixes.py -- Note 6（重复键敏感性）+ 修正 caption 编号。

数据来源：v8_1/DUPLICATE_KEY_SENSITIVITY.json
（该脚本的分解已逐位验证与 all_models_by_rowset.csv 一致，最大差 3.331e-16）
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
D = json.loads(Path("v8_1/DUPLICATE_KEY_SENSITIVITY.json").read_text(encoding="utf-8"))
bm = D["by_mode"]
A, F, M = bm["ALL"], bm["FIRST"], bm["MEAN"]
TOL = D["decomp_verified_max_abs_diff"]

lo = 100 * min(A["attribution"], F["attribution"], M["attribution"])
hi = 100 * max(A["attribution"], F["attribution"], M["attribution"])

NOTE6 = r"""\section*{Supplementary Note 6 | Duplicate-key sensitivity}
The controlled-split universe contains 7{,}523 mutation records but only 4{,}812 distinct
(complex, mutation) keys: 2{,}011 keys carry more than one experimental record, and 1{,}405 of
those disagree in $\Delta\Delta G$. This raises a legitimate concern: could part of the
random-split performance, and hence part of the apparent collapse under family-out-of-domain
evaluation, be an artefact of the same key appearing more than once in the evaluation set?

Because the stored prediction artifacts do not retain per-row fold assignments, the split
cannot be regenerated here; we therefore quantify the sensitivity \emph{within} the stored
evaluation set by recomputing every quantity under three key treatments:
\begin{itemize}
  \item \textbf{ALL} --- all 7{,}523 records (the reported analysis);
  \item \textbf{FIRST} --- one record per key, keeping the first occurrence (4{,}812 records);
  \item \textbf{MEAN} --- one record per key, with $y$ and $\hat{y}$ both averaged within the
        key (4{,}812 records).
\end{itemize}
The decomposition used here is the unweighted identity of Supplementary Note 1, and was
verified to reproduce the frozen per-configuration values to __TOL__ across the 36
controlled-split configurations.

\begin{center}
\begin{tabular}{@{}llrrrr@{}}
\toprule
Regime & Key treatment & Median $r$ & Median $C_{\mathrm{B}}$
  & Median $C_{\mathrm{W}}$ & Median macro $\rho$ \\
\midrule
Random & ALL & __R_ALL__ & $+0.4275$ & $+0.2591$ & 0.0707 \\
Random & FIRST & $+0.6378$ & $+0.3920$ & $+0.2465$ & 0.0886 \\
Random & MEAN & $+0.6338$ & $+0.3705$ & $+0.2647$ & 0.0977 \\
\midrule
Family-OOD & ALL & __O_ALL__ & $-0.0072$ & $+0.1580$ & 0.0377 \\
Family-OOD & FIRST & $+0.1614$ & $-0.0122$ & $+0.1759$ & 0.0354 \\
Family-OOD & MEAN & $+0.2181$ & $+0.0197$ & $+0.1956$ & 0.0347 \\
\bottomrule
\end{tabular}
\end{center}

The paired Random $\to$ Family-OOD differences are:
\begin{center}
\begin{tabular}{@{}lrrrr@{}}
\toprule
Key treatment & $\Delta r$ & $\Delta C_{\mathrm{B}}$ & $\Delta C_{\mathrm{W}}$
  & Between-target share \\
\midrule
ALL   & __DA__ & __DBA__ & __DWA__ & __SHA__\,\% \\
FIRST & __DF__ & __DBF__ & __DWF__ & __SHF__\,\% \\
MEAN  & __DM__ & __DBM__ & __DWM__ & __SHM__\,\% \\
\bottomrule
\end{tabular}
\end{center}

Three observations follow. First, collapsing duplicate keys \emph{reduces} the random-split
correlation (__R_ALL__ to $+0.6378$ under FIRST) while slightly \emph{increasing} the
family-out-of-domain correlation (__O_ALL__ to $+0.1614$); the collapse is therefore, if
anything, understated rather than inflated by repeated records. Second, $\Delta C_{\mathrm{B}}$
remains large and negative under every treatment, and its share of the summed absolute
component change lies between __LO__\,\% and __HI__\,\%: the qualitative conclusion that the
between-target component carries most of the decline does not depend on how repeated keys are
handled. Third, the magnitude of the decline is sensitive to the treatment
(__DA__ versus __DM__), so the absolute size of the collapse should be quoted with the key
treatment stated.

\emph{Limits.} This is an evaluation-set-internal sensitivity analysis. It cannot rule out
duplicate-key leakage across the train/test boundary, because that requires regenerating the
splits under key-level grouping and retraining, which the stored artifacts do not permit and
which we did not perform. The 1{,}405 keys with disagreeing $\Delta\Delta G$ values also
indicate that a (complex, mutation) key does not always identify a unique experimental
condition; determining whether those conditions should enter the key definition is a
data-curation decision that we did not make.

"""

repl = {
    "__TOL__": f"{TOL:.1e}",
    "__R_ALL__": f"{A['median_r_random']:+.4f}",
    "__O_ALL__": f"{A['median_r_ood']:+.4f}",
    "__DA__": f"{A['d_r']:+.4f}", "__DBA__": f"{A['d_C_B']:+.4f}",
    "__DWA__": f"{A['d_C_W']:+.4f}", "__SHA__": f"{100*A['attribution']:.1f}",
    "__DF__": f"{F['d_r']:+.4f}", "__DBF__": f"{F['d_C_B']:+.4f}",
    "__DWF__": f"{F['d_C_W']:+.4f}", "__SHF__": f"{100*F['attribution']:.1f}",
    "__DM__": f"{M['d_r']:+.4f}", "__DBM__": f"{M['d_C_B']:+.4f}",
    "__DWM__": f"{M['d_C_W']:+.4f}", "__SHM__": f"{100*M['attribution']:.1f}",
    "__LO__": f"{lo:.1f}", "__HI__": f"{hi:.1f}",
}
for k, v in repl.items():
    NOTE6 = NOTE6.replace(k, v)

sp = Path("latex_v6/supplementary_v6.tex")
s = sp.read_text(encoding="utf-8")
anchor = "\\section*{Supplementary Figures}"
if "Supplementary Note 6" not in s and anchor in s:
    s = s.replace(anchor, NOTE6 + anchor, 1)
    print("Note 6 inserted")
else:
    print("Note 6 already present or anchor missing")

# ---- caption 编号清理 ----
fixes = [
    ("  \\caption{Model pools used in the rank-agreement sensitivity analysis Table~S4.}",
     "  \\caption{Model pools used in the rank-agreement sensitivity analysis.}"),
    ("(Table~S4).}\\label{tab:s4pools}", "Table~S4.}\\label{tab:s4pools}"),
]
for old, new in fixes:
    if old in s:
        s = s.replace(old, new, 1); print("  fixed:", old.strip()[:56])
    else:
        print("  already ok:", old.strip()[:56])
sp.write_text(s, encoding="utf-8")

# ---- 主稿 Figure 6 caption ----
mp = Path("latex_v6/manuscript_v6.tex")
ms = mp.read_text(encoding="utf-8")
old6 = ("(A) Pooled\n"
        "versus within-target performance for five independently developed methods on a matched,\n"
        "unambiguous mapping.")
new6 = ("(A) Pooled versus within-target performance for five independently developed methods, each\n"
        "evaluated on its own uniquely mapped subset of a common unambiguous annotation universe.")
if old6 in ms:
    ms = ms.replace(old6, new6, 1); print("  fixed: Figure 6 caption")
else:
    print("  MISS: Figure 6 caption")
mp.write_text(ms, encoding="utf-8")
