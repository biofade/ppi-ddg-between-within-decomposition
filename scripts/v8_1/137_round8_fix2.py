# -*- coding: utf-8 -*-
"""scripts/v8_1/137_round8_fix2.py -- 补上 136 未命中的 3 处（换行/已改）。"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
O = json.loads(Path("v8_1/EXTERNAL_ORIENTATION.json").read_text(encoding="utf-8"))
RA, RB = O["raw"], O["oriented"]
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []

# ---- (2) Methods permutation p ----
OLD = ("For the five-method external panel, the exact permutation\n"
       "distribution of Spearman's $\\rho$ is enumerated over all $5! = 120$ relabellings, giving\n"
       "$p = 0.083$; asymptotic approximations are not used because the panel is far too small "
       "for\nthem to be reliable.")
NEW = ("For the five-method external panel, the exact permutation\n"
       "distribution of Spearman's $\\rho$ is enumerated over all $5! = 120$ relabellings, giving\n"
       f"$p = {RA['perm_p']:.3f}$ for the raw columns and $p = {RB['perm_p']:.3f}$ after the\n"
       "exploratory sign-orientation; asymptotic approximations are not used because the panel is\n"
       "far too small for them to be reliable.")
if OLD in s:
    s = s.replace(OLD, NEW, 1); done.append("methods-perm")
else:
    miss.append("methods-perm")

# ---- (3) Table 3 表注：删除 "no sign flip" ----
OLD3 = ("Prediction columns are used exactly\n"
        "      as published, with no sign flip and no recalibration (Supplementary Note~8). RMSE is\n"
        "      reported only for methods whose output is on a $\\Delta\\Delta G$-like scale;")
NEW3 = ("Raw prediction columns are\n"
        "      ingested and preserved exactly as published, with no rescaling or recalibration\n"
        "      (Supplementary Note~8). A sign multiplier is applied only for the orientation\n"
        "      sensitivity described in Methods. RMSE is reported only for methods whose output is\n"
        "      on a $\\Delta\\Delta G$-like scale;")
if OLD3 in s:
    s = s.replace(OLD3, NEW3, 1); done.append("table3-note")
else:
    miss.append("table3-note")

# ---- (3b) Methods 的 "no sign flip"（external prediction semantics 小节）----
OLD3b = "exactly as supplied: no sign flip, no rescaling and no recalibration against the evaluation\nlabels."
NEW3b = ("exactly as published, with no rescaling and no recalibration against the evaluation\n"
         "labels. A sign multiplier is applied only for the exploratory orientation sensitivity\n"
         "described below, and never to the raw columns themselves.")
if OLD3b in s:
    s = s.replace(OLD3b, NEW3b, 1); done.append("methods-semantics-nosignflip")
else:
    miss.append("methods-semantics-nosignflip")

p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
