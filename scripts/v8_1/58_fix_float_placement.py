# -*- coding: utf-8 -*-
"""scripts/v8_1/58_fix_float_placement.py -- 让图紧跟其引用位置，不再全部漂到文末。

现象：v6 转 LaTeX 后 6 张主图全部落在第 9-12 页（正文只到第 6 页），
      因为图环境用 [t] 且图很大，LaTeX 把它们排到文末。

修复：
  1. 图环境改为 [htbp]（先试此处 -> 页顶 -> 页底 -> 单独浮动页）
  2. 放宽浮动参数：topfraction / bottomfraction / textfraction /
     floatpagefraction，允许页面容纳大图
  3. 表环境同样改为 [htbp]
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 1) 浮动参数：放在 geometry 之后
anchor = "\\graphicspath{{../figures_v5/}{./}}"
FLOATS = (
    "% ---- float placement: allow large figures to sit near their citation ----\n"
    "\\renewcommand{\\topfraction}{0.92}\n"
    "\\renewcommand{\\bottomfraction}{0.75}\n"
    "\\renewcommand{\\textfraction}{0.06}\n"
    "\\renewcommand{\\floatpagefraction}{0.75}\n"
    "\\setcounter{topnumber}{3}\n"
    "\\setcounter{bottomnumber}{2}\n"
    "\\setcounter{totalnumber}{4}\n\n"
    + anchor)
if anchor in s:
    s = s.replace(anchor, FLOATS, 1)
    print("float parameters injected")
else:
    print("WARNING: graphicspath anchor not found")

# 2) 图/表环境 [t] -> [htbp]
n1 = s.count("\\begin{figure}[t]")
n2 = s.count("\\begin{table}[t]")
s = s.replace("\\begin{figure}[t]", "\\begin{figure}[htbp]")
s = s.replace("\\begin{table}[t]", "\\begin{table}[htbp]")
# 补充材料里的 [h] 也放宽
n3 = s.count("\\begin{figure}[h]")
s = s.replace("\\begin{figure}[h]", "\\begin{figure}[htbp]")
print(f"float specifiers: figure[t] x{n1} -> htbp ; table[t] x{n2} -> htbp ; figure[h] x{n3} -> htbp")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
