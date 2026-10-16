# -*- coding: utf-8 -*-
"""scripts/v8_1/82_fix_stray_backslash.py -- 清除补丁误插入的多余换行命令。

问题：补丁用非 raw 字符串书写替换文本，``\\`` 被 Python 解析为转义后写成了单个
反斜杠，落成 LaTeX 的换行命令 ``\\``，于是
  l.539 \\subsection*{Target-weighting schemes}
  l.208 ...(Table~S4).}\\label{tab:s4pools}
触发 "There's no line here to end"。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")
done = []

p = OUT / "manuscript_v6.tex"
s = p.read_text(encoding="utf-8")
if "\\\\subsection*{Target-weighting schemes}" in s:
    s = s.replace("\\\\subsection*{Target-weighting schemes}",
                  "\n\\subsection*{Target-weighting schemes}")
    p.write_text(s, encoding="utf-8"); done.append("manuscript subsection")

p = OUT / "supplementary_v6.tex"
s = p.read_text(encoding="utf-8")
if "(Table~S4).}\\\\label{tab:s4pools}" in s:
    s = s.replace("(Table~S4).}\\\\label{tab:s4pools}",
                  "(Table~S4).}\\label{tab:s4pools}")
    p.write_text(s, encoding="utf-8"); done.append("supplementary label")

print("fixed:", done if done else "nothing to fix")
