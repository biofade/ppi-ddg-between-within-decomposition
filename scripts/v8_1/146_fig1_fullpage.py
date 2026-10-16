# -*- coding: utf-8 -*-
"""scripts/v8_1/146_fig1_fullpage.py -- Figure 1 改为满页 figure page。

图本身已改为接近竖版比例（0.871），因此可以真正占满一整页；
LaTeX 侧把 [htbp] 改为 [p]（专用浮动页）并去掉 0.82\\textheight 的高度上限，
使宽度约束（\\textwidth）成为唯一约束，从而得到最大可读尺寸。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []

OLD = ("%% ---- Figure 1 ----\n"
       "\\begin{figure}[htbp]\n"
       "  \\centering\n"
       "  \\figfit{Fig1_overview_EN.pdf}\n")
NEW = ("%% ---- Figure 1 ----（满页 figure page）\n"
       "\\begin{figure}[p]\n"
       "  \\centering\n"
       "  \\includegraphics[width=\\textwidth,height=0.90\\textheight,keepaspectratio]"
       "{Fig1_overview_EN.pdf}\n")
if OLD in s:
    s = s.replace(OLD, NEW, 1); done.append("fig1-fullpage")
else:
    miss.append("fig1-fullpage")

p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
    i = s.find("%% ---- Figure 1")
    print(repr(s[i:i + 200]))
