# -*- coding: utf-8 -*-
"""scripts/v8_1/138_round8_supp.py -- 补充材料第八轮修正。

  1) Note 7 两套表同步为「raw + oriented」双口径，并给出 oriented 的排名与 rho/tau/p
  2) Note 7 common-intersection 的能量型 RMSE 恢复（此前被误清空）
  3) Table S3 的 mu_c -> o_m，并加 oriented 列
  4) Note 8 改为诚实口径（o_m 由 evaluation correlation 判定）
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
O = json.loads(Path("v8_1/EXTERNAL_ORIENTATION.json").read_text(encoding="utf-8"))
RA, RB = O["raw"], O["oriented"]
CI = json.loads(Path("v8_1/EXTERNAL_COMMON_INTERSECTION.json").read_text(encoding="utf-8"))
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []

# ---------- 1) + 2) Note 7：重建两套表 ----------
VA = CI["coverage_variable"]; VI = CI["coverage_intersection"]
MULT = O["multipliers"]
# common intersection 的能量型 RMSE（原始重算值）
RMSE_INT = {"Flex-ddG": "2.1620", "DDGPred": "2.1265", "ddG_FoldX": "2.3623"}
ORDER = ["DDGPred", "ddG_FoldX", "Flex-ddG", "proteinMPNN", "ESM-IF"]
LABEL = {"ddG_FoldX": "ddG\\_FoldX"}


def row(rec, m, mu, rmse_ok):
    nm = LABEL.get(m, m)
    return (f"{nm} & {int(rec[m]['n_rows']):,} & {int(rec[m]['n_ge5'])} & "
            f"{mu*rec[m]['pearson']:+.4f} & {mu*rec[m]['macro']:+.4f} & "
            f"{mu*rec[m]['C_B']:+.4f} & {mu*rec[m]['C_W']:+.4f} & {rmse_ok} " + chr(92) * 2)


tblA_raw = "\n".join(row(VA, m, 1, f"{VA[m]['rmse']:.4f}" if m in RMSE_INT else "---")
                    for m in ORDER)
tblA_or = "\n".join(row(VA, m, MULT.get(m, 1),
                        f"{VA[m]['rmse']:.4f}" if m in RMSE_INT else "---")
                    for m in ORDER)
tblB_raw = "\n".join(row(VI, m, 1, RMSE_INT.get(m, "---")) for m in ORDER)
tblB_or = "\n".join(row(VI, m, MULT.get(m, 1), RMSE_INT.get(m, "---")) for m in ORDER)

NEW_TABLES = r"""\begin{center}
\begin{tabular}{@{}lrrrrrrr@{}}
\toprule
Method & Rows & Targets ($n \geq 5$) & Pooled Pearson & Macro Spearman
  & $C_{\mathrm{B}}$ & $C_{\mathrm{W}}$ & RMSE \\
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(A) variable coverage --- each method's own uniquely mapped subset, raw columns}} \\
__A_RAW__
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(A$'$) same subsets, exploratory sign-oriented values}} \\
__A_OR__
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(B) common intersection --- raw columns}} \\
__B_RAW__
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(B$'$) common intersection, exploratory sign-oriented values}} \\
__B_OR__
\bottomrule
\end{tabular}
\end{center}

RMSE is shown only for Flex-ddG, FoldX and DDGPred, whose columns are on a $\Delta\Delta G$-like
scale; ESM-IF and proteinMPNN return dimensionless mutation-effect scores for which an error in
kcal/mol is not defined (Supplementary Note~8).

The common intersection retains __NI__ of the 1{,}781 uniquely mappable rows (__PCT__\,\%),
covering __NT__ targets of which __NG__ have $n \geq 5$. Comparing the treatments:

\begin{itemize}
  \item on the \emph{raw} columns the two rankings are strongly positively aligned
        (Spearman $\rho = +__RA_RHO__$, exact permutation $p = __RA_P__$), as in the main text;
  \item after the exploratory sign-orientation the same ranks show no significant agreement
        ($\rho = +__RB_RHO__$, $p = __RB_P__$), and ESM-IF moves from last to second on the
        pooled axis while remaining last on the within-target axis;
  \item in both treatments DDGPred ranks first on both axes.
\end{itemize}

The coverage imbalance is therefore not what drives the difference between the two conclusions:
restricting all five methods to their common intersection leaves the raw panel aligned and the
oriented panel unaligned, exactly as in the full uniquely mapped subsets. The panel remains small
(five methods) and underpowered under either treatment.
"""
NEW_TABLES = (NEW_TABLES.replace("__A_RAW__", tblA_raw).replace("__A_OR__", tblA_or)
              .replace("__B_RAW__", tblB_raw).replace("__B_OR__", tblB_or)
              .replace("__NI__", f"{CI['n_intersection_rows']:,}")
              .replace("__PCT__", f"{100*CI['n_intersection_rows']/1781:.1f}")
              .replace("__NT__", str(CI["n_intersection_targets"]))
              .replace("__NG__", str(CI["n_intersection_ge5"]))
              .replace("__RA_RHO__", f"{RA['rho']:.2f}")
              .replace("__RA_P__", f"{RA['perm_p']:.3f}")
              .replace("__RB_RHO__", f"{RB['rho']:.2f}")
              .replace("__RB_P__", f"{RB['perm_p']:.3f}"))

import re  # noqa: E402
m = re.search(r"\\begin\{center\}\n\\begin\{tabular\}\{@\{\}lrrrrrrr@\{\}\}.*?"
              r"(?=\\section\*\{Supplementary Note 8)", s, re.S)
if m:
    s = s[:m.start()] + NEW_TABLES + "\n" + s[m.end():]
    done.append("note7-tables-rebuilt")
else:
    miss.append("note7-tables-anchor")

# 删除 Note 7 里残留的旧结论 bullet（已被新 bullet 取代）
for old_bullet in (
    "  \\item the pooled ranking is identical (DDGPred $>$ FoldX $>$ Flex-ddG $>$ proteinMPNN $>$\n"
    "        ESM-IF);",
    "  \\item the within-target ranking is identical (DDGPred $>$ Flex-ddG $>$ FoldX $>$ "
    "proteinMPNN\n        $>$ ESM-IF), so the single adjacent swap between FoldX and Flex-ddG "
    "persists;",
):
    if old_bullet in s:
        s = s.replace(old_bullet, "", 1); done.append("note7-old-bullet-removed")

p.write_text(s, encoding="utf-8")
print("note7:", done)
if miss:
    print("MISSED:", miss)
