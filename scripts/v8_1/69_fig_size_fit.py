# -*- coding: utf-8 -*-
"""scripts/v8_1/69_fig_size_fit.py -- 让图在保持比例的前提下尽可能大。

新图比例差异较大（Figure 1 为 2.33:1 横幅，Figure 5 为 1.60:1），
固定 width=\\textwidth 会让横幅图过小。改用「宽不超文本宽、高不超页高」的
等比上限，两轴同时约束，保证既填满页面又不溢出、且绝不改变长宽比。
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 在导言区加入 adjustbox / 高度上限宏
anchor = "\\usepackage{graphicx}"
if anchor in s and "maxheight" not in s:
    s = s.replace(
        anchor,
        anchor + "\n"
        "% 等比缩放：宽不超过文本宽，高不超过正文高的 0.82\n"
        "\\newcommand{\\figfit}[2][\\textwidth]{%\n"
        "  \\includegraphics[width=#1,height=0.82\\textheight,keepaspectratio]{#2}}",
        1)
    print("added \\figfit macro")

n = s.count("\\includegraphics[width=\\textwidth]{")
s = s.replace("\\includegraphics[width=\\textwidth]{", "\\figfit{")
print(f"converted {n} includegraphics calls to \\figfit")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
