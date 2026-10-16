# -*- coding: utf-8 -*-
"""scripts/v8_1/135_fig6_caption_oriented.py -- Figure 6 caption 与图内/oriented 结论对齐。"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
O = json.loads(Path("v8_1/EXTERNAL_ORIENTATION.json").read_text(encoding="utf-8"))
RB, RA = O["oriented"], O["raw"]
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")

OLD = ("for visual comparison only; raw values are printed in each cell. Rank agreement is\n"
       "  positive ($\\rho = +0.90$, exact permutation $p = 0.083$) but the panel is underpowered\n"
       "  with five methods.")
NEW = (f"for visual comparison\n"
       f"  only; sign-oriented values are printed in each cell. After orientation the two axes show\n"
       f"  no significant rank agreement ($\\rho = +{RB['rho']:.2f}$, exact permutation\n"
       f"  $p = {RB['perm_p']:.3f}$); used raw the same numbers would give $\\rho = +{RA['rho']:.2f}$\n"
       f"  ($p = {RA['perm_p']:.3f}$). The panel is underpowered with five methods.")
if OLD in s:
    s = s.replace(OLD, NEW, 1); p.write_text(s, encoding="utf-8")
    print("Figure 6 caption aligned to oriented result")
else:
    print("MISS: caption anchor")
    i = s.find("Rank agreement is")
    print(repr(s[max(0, i - 120):i + 200]))
