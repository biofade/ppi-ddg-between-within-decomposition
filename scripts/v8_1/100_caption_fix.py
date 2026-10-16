# -*- coding: utf-8 -*-
"""scripts/v8_1/100_caption_fix.py -- 空白容忍地修正 Figure 6 caption 与 Table S4 caption。"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
done, miss = [], []

mp = Path("latex_v6/manuscript_v6.tex")
ms = mp.read_text(encoding="utf-8")
pat6 = re.compile(
    r"\(A\)\s+Pooled\s+versus within-target performance for five independently developed "
    r"methods on a matched,\s+unambiguous mapping\.",
    re.S)
new6 = ("(A) Pooled versus within-target performance for five independently developed\n"
        "  methods, each evaluated on its own uniquely mapped subset of a common unambiguous\n"
        "  annotation universe.")
ms2, n = pat6.subn(lambda m: new6, ms, count=1)
if n:
    mp.write_text(ms2, encoding="utf-8"); done.append("Figure 6 caption")
else:
    miss.append("Figure 6 caption")

sp = Path("latex_v6/supplementary_v6.tex")
ss = sp.read_text(encoding="utf-8")
pat = re.compile(r"\\caption\{Model pools used in the rank-agreement sensitivity analysis"
                 r"\s*Table~S4\.\}")
ss2, n2 = pat.subn(lambda m: r"\caption{Model pools used in the rank-agreement sensitivity "
                               r"analysis.}", ss, count=1)
if n2:
    sp.write_text(ss2, encoding="utf-8"); done.append("Table S4 caption")
else:
    miss.append("Table S4 caption")

print("fixed:", done)
if miss:
    print("MISSED:", miss)
