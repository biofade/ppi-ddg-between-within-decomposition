# -*- coding: utf-8 -*-
"""scripts/v8_1/122_external_units_fix.py -- 修正 external RMSE 单位问题 + 新增 provenance 表。

依据：
  * v8_1/EXTERNAL_PROVENANCE.{csv,json} —— 原始量纲判定
  * scripts/v8_1/121_verify_raw_reproduces_paper.py —— 证实无隐藏变换
  * scripts/v8_1/04_external_audits.py —— 项目自身对同类打分的既有规则：
        "score only; NOT calibrated kcal/mol ddG"; RMSE = NOT_APPLICABLE
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
EP = pd.read_csv("v8_1/EXTERNAL_PROVENANCE.csv")
applied, missed = [], []


def patch(name, pairs):
    p = Path("latex_v6") / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


# ---------- 主稿：Table 3 数值行 ----------
# 表体由 78_audit_patch / 53 生成器写出；这里直接改已生成的 .tex 单元格。
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
RMSE_MAP = {"ESM-IF": "2.3587", "ProteinMPNN": "2.3199"}
changed = []
for name, val in RMSE_MAP.items():
    # 匹配 "ESM-IF & ... & 2.3587 \\" 形式的行
    import re
    pat = re.compile(r"(" + re.escape(name) + r"\s*&[^\n]*?)\s*&\s*" + re.escape(val) + r"\s*\\\\")
    s2, n = pat.subn(lambda m: m.group(1) + " & --- \\\\", s, count=1)
    if n:
        s = s2; changed.append(name)
    else:
        miss_pat = name
        missed.append(f"manuscript:rmse-{name}")
p.write_text(s, encoding="utf-8")
applied += [f"manuscript:rmse-{c}" for c in changed]

# ---------- 主稿：Table 3 表注 ----------
patch("manuscript_v6.tex", [(
    "table3-note",
    "\\item Each method is evaluated on its own uniquely mapped subset of the common annotation\n"
    "      namespace, so mutation counts differ slightly (1{,}755--1{,}776). RMSE is reported\n"
    "      only for outputs with $\\Delta\\Delta G$-like scaling.",
    "\\item Each method is evaluated on its own uniquely mapped subset of the common annotation\n"
    "      namespace, so mutation counts differ slightly (1{,}755--1{,}776), and their common\n"
    "      intersection is analysed in Supplementary Note~7. Prediction columns are used exactly\n"
    "      as published, with no sign flip and no recalibration (Supplementary Note~8). RMSE is\n"
    "      reported only for methods whose output is on a $\\Delta\\Delta G$-like scale;\n"
    "      ESM-IF and ProteinMPNN return dimensionless mutation-effect scores (raw standard\n"
    "      deviations 0.093 and 0.012) and are therefore shown as ``---'', since an error in\n"
    "      kcal/mol is not defined for them."
)])

# ---------- 主稿：Methods 增补一段 ----------
patch("manuscript_v6.tex", [(
    "methods-external-units",
    "\\subsection*{External prediction mapping and ambiguity handling}",
    "\\subsection*{External prediction semantics and units}\n"
    "The five external methods are drawn from a single published annotation table, and each\n"
    "prediction column is used exactly as supplied: no sign flip, no rescaling and no\n"
    "recalibration against the evaluation labels (verified by reproducing every reported\n"
    "correlation directly from the raw columns; Supplementary Note~8). Three of the columns\n"
    "(Flex-ddG, FoldX and DDGPred) are quantitative $\\Delta\\Delta G$ predictions on a\n"
    "kcal/mol-like scale, whereas ESM-IF and ProteinMPNN are inverse-folding sequence models\n"
    "whose published outputs are dimensionless mutation-effect scores rather than energies.\n"
    "We therefore report RMSE only for the first group, and we do not treat the score scale of\n"
    "the latter two as a calibrated binding energy. Correlation-based quantities, which are\n"
    "invariant to affine rescaling, are reported for all five methods.\n\n"
    "\\subsection*{External prediction mapping and ambiguity handling}"
)])

# ---------- 补充材料：Note 8 + Table S3 ----------
rows = []
for r in EP.itertuples():
    rmse = "yes" if bool(r.rmse_valid) else "**no**"
    rows.append(f"{r.method.replace('_', chr(92)+'_')} & \\texttt{{{r.raw_column.replace('_', chr(92)+'_')}}} "
                f"& {r.raw_quantity} & {r.units} & {r.transformation} & {rmse} " + chr(92)*2)
TBL = ("\\begin{landscape}\\begin{table}[h]\\centering\\scriptsize\n"
       "\\begin{threeparttable}\n"
       "\\caption{Provenance and score semantics of the external prediction columns "
       "(Table~S3).}\\label{tab:s3external}\n"
       "\\begin{tabular}{@{}llllll@{}}\\toprule\n"
       "Method & Raw column & Raw quantity & Units & Transformation applied here "
       "& RMSE meaningful? \\\\\n\\midrule\n"
       + "\n".join(rows) + "\n"
       "\\bottomrule\\end{tabular}\n"
       "\\begin{tablenotes}\\footnotesize\n"
       "\\item All five columns come from one published annotation table of SKEMPI-derived\n"
       "protein--protein interaction mutations. Values are used exactly as published; we apply\n"
       "no sign flip, no affine calibration and no recalibration against the evaluation labels.\n"
       "We verified this by recomputing every reported correlation directly from the raw columns,\n"
       "which reproduces the reported values to four decimals for all five methods. Raw standard\n"
       "deviations differ by more than an order of magnitude across columns (1.31, 1.55 and 1.92\n"
       "for the three $\\Delta\\Delta G$-like predictors versus 0.093 and 0.012 for the two\n"
       "inverse-folding scores), which is why an RMSE in kcal/mol is defined only for the former.\n"
       "\\end{tablenotes}\\end{threeparttable}\\end{table}\\end{landscape}\n")

NOTE8 = r"""\section*{Supplementary Note 8 | Provenance, orientation and units of the external predictions}
The five external methods are not separate downloads: they are five annotation columns of one
published, SKEMPI-derived table of protein--protein interaction mutations, which we read
directly (Table~S3). This matters for interpretation in three ways.

\emph{No transformation is applied.} For every method we use the published column exactly as
supplied. We apply no sign flip, no affine rescaling and no recalibration against the
evaluation labels, and therefore no evaluation label enters the construction of any prediction.
We verified this by recomputing each reported correlation directly from the raw column on the
same uniquely mapped subset: the values reproduce to four decimals for all five methods
(e.g.\ ESM-IF $-0.2250$ and ProteinMPNN $+0.1509$), so the reported correlations are properties
of the published scores and not of any post-processing.

\emph{Orientation.} All five columns are used with higher values meaning stronger binding, the
same orientation as the $-\Delta\Delta G$ convention adopted here, so no method required a sign
flip. We note explicitly that ESM-IF carries a negative pooled correlation
($-0.2250$): the sign is a property of that model's published score on this benchmark and not
an artefact of how we oriented it.

\emph{Units, and why two RMSE values are not reported.} Flex-ddG, FoldX and DDGPred are
quantitative $|\Delta\Delta G|$-like predictors whose columns span standard deviations of
1.31, 1.92 and 1.55 respectively, so an error in kcal/mol is meaningful for them. ESM-IF and
ProteinMPNN are inverse-folding sequence models: their published outputs are dimensionless
mutation-effect scores with standard deviations of 0.093 and 0.012, roughly two orders of
magnitude narrower than the experimental $\Delta\Delta G$ spread. An RMSE expressed in kcal/mol
is therefore not defined for those two columns, and we report ``---''. This follows the rule
already applied elsewhere in this work, where a dimensionless mutation-effect score is recorded
with \texttt{RMSE = NOT\_APPLICABLE} rather than converted into an apparent energy error.
Correlation-based quantities --- pooled Pearson, macro within-target Spearman,
$C_{\mathrm{B}}$ and $C_{\mathrm{W}}$ --- are invariant to affine rescaling and are reported for
all five methods.

\emph{Scope of the comparison.} Because the two score-based methods are not on the same scale
as the three energy-based ones, the RMSE column in Table~3 should not be used to compare across
those groups; the pooled and within-target correlations are the quantities that are comparable
across all five methods, and they are the basis of the rank-agreement analysis.

"""

p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
anchor = "\\section*{Machine-readable supplementary data}"
if "Supplementary Note 8" not in s and anchor in s:
    s = s.replace(anchor, TBL + "\n" + NOTE8 + "\n" + anchor, 1)
    p.write_text(s, encoding="utf-8")
    print("Note 8 + Table S3 inserted")
else:
    p.write_text(s, encoding="utf-8")
    print("Note 8 already present or anchor missing")

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
