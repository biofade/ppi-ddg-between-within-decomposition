# -*- coding: utf-8 -*-
"""scripts/v8_1/118_abstract_dereplicated.py -- 修正 Abstract 的 "de-replicated" 措辞。"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
OLD = ("matched, de-replicated and provenance-filtered evaluation sets.")
NEW = ("matched and provenance-filtered evaluation sets, with explicit\n"
       "sensitivity analyses for repeated mutation keys.")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    p.write_text(s, encoding="utf-8")
    print("fixed: abstract de-replicated")
else:
    print("MISS")
