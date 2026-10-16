# -*- coding: utf-8 -*-
"""scripts/v8_1/141_supp_stale_cleanup.py -- 清理补充材料残留的 "no sign flip" 与 Table~S4。"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []

REPL = [
    # Note 8：No transformation 段
    ("\\emph{No transformation is applied.} For every method we use the published column exactly as\n"
     "supplied. We apply no sign flip, no affine rescaling and no recalibration against the\n"
     "evaluation labels, and therefore no evaluation label enters the construction of any\n"
     "prediction.",
     "\\emph{The published columns are preserved unmodified.} For every method we ingest the\n"
     "published column exactly as distributed. We apply no affine rescaling and no recalibration\n"
     "against the evaluation labels. The only operation we apply anywhere is the single sign\n"
     "multiplier $o_m$ of Eq.~(orientation) in the main text, which is used solely for the\n"
     "exploratory orientation sensitivity and never modifies the raw columns reported here.\n"
     "Because $o_m$ is chosen from the evaluation data, that sensitivity is explicitly\n"
     "label-informed; the raw analysis requires no such choice and is the primary result."),
    # Table S3 表注
    ("protein--protein interaction mutations. Values are used exactly as published; we apply\n"
     "no sign flip, no affine calibration and no recalibration against the evaluation labels.",
     "protein--protein interaction mutations. Raw values are used exactly as published; we apply\n"
     "no affine calibration and no recalibration against the evaluation labels. The $o_m$ column\n"
     "records the sign multiplier used only for the exploratory orientation sensitivity."),
    # Table~S4 -> S3
    ("Table~S4 records each raw column with its raw", "Table~S3 records each raw column with its raw"),
]
for old, new in REPL:
    if old in s:
        s = s.replace(old, new, 1); done.append(old[:48])
    else:
        miss.append(old[:48])
p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
