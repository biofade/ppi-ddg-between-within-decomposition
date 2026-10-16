# -*- coding: utf-8 -*-
"""scripts/v8_1/142_supp_notransform.py -- 重写 Note 8 的 "No transformation" 段（精确文本）。"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")

OLD = ("\\emph{No transformation is applied.} For every method we use the published column exactly as\n"
       "supplied. We apply no sign flip, no affine rescaling and no recalibration against the\n"
       "evaluation labels, and therefore no evaluation label enters the construction of any prediction.\n"
       "We verified this by recomputing each reported correlation directly from the raw column on the\n"
       "same uniquely mapped subset: the values reproduce to four decimals for all five methods\n"
       "(e.g.\\ ESM-IF $-0.2250$ and ProteinMPNN $+0.1509$), so the reported correlations are properties\n"
       "of the publi")
NEW = ("\\emph{The published columns are preserved unmodified.} For every method we ingest the\n"
       "published column exactly as distributed. We apply no affine rescaling and no recalibration,\n"
       "and the raw values reproduce to four decimals when recomputed directly from the published\n"
       "columns on the same uniquely mapped subset (e.g.\\ ESM-IF $-0.2250$ and ProteinMPNN\n"
       "$+0.1509$). The only operation applied anywhere is the single sign multiplier $o_m$\n"
       "introduced in the main text Methods (\\emph{Score orientation and its limits}). It is used\n"
       "solely for the exploratory\n"
       "orientation sensitivity and never modifies the raw columns reported here; because $o_m$ is\n"
       "selected from the evaluation data, that sensitivity is explicitly label-informed, whereas\n"
       "the raw analysis requires no such choice and is the primary result. So the reported raw\n"
       "correlations are properties\n"
       "of the publi")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    p.write_text(s, encoding="utf-8")
    print("Note 8 'No transformation' paragraph rewritten")
else:
    print("MISS")
