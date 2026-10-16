# -*- coding: utf-8 -*-
"""scripts/v8_1/70_fix_si_macro.py -- 修正 SI 模板中 figfit 宏的反斜杠数量。

源码 L662/663 用了 4 个反斜杠（在 raw 模板里会原样输出为 4 个），
LaTeX 需要 2 个。改为 2 个。
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
B = open(P, "rb").read()

bad1 = b"\\\\\\\\newcommand{\\\\\\\\figfit}[2][\\\\\\\\textwidth]{%"
good1 = b"\\\\newcommand{\\\\figfit}[2][\\\\textwidth]{%"
bad2 = b"  \\\\\\\\includegraphics[width=#1,height=0.82\\\\\\\\textheight,keepaspectratio]{#2}}"
good2 = b"  \\\\includegraphics[width=#1,height=0.82\\\\textheight,keepaspectratio]{#2}}"

n1 = B.count(bad1)
n2 = B.count(bad2)
B = B.replace(bad1, good1).replace(bad2, good2)
open(P, "wb").write(B)
print(f"replaced: newcommand x{n1}, includegraphics x{n2}")

import ast  # noqa: E402
ast.parse(open(P, encoding="utf-8").read())
print("syntax OK")

# 复核源码里两处宏定义应完全一致
for i, l in enumerate(open(P, encoding="utf-8").read().split("\n"), 1):
    if "figfit" in l and "newcommand" in l:
        print(f"  L{i} {l!r}")
