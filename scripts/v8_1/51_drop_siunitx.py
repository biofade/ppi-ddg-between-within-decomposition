# -*- coding: utf-8 -*-
"""scripts/v8_1/51_drop_siunitx.py -- 移除 siunitx 依赖。

`\\SI{3.3e-16}{}` 在 XeLaTeX + 该字体组合下触发了
"Improper alphabetic constant \\mitDelta" 一类低层错误。
本文只用到一处科学计数法，直接用数学模式写更稳，也更简洁。
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/49_build_latex.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 1) 去掉 \usepackage{siunitx}
s = s.replace("\\usepackage{siunitx}\n", "")
# 2) 那一处 SI 改为普通数学
s = s.replace("(residual\n$\\leq \\SI{3.3e-16}{}$)",
              "(residual\n$\\leq 3.3 \\times 10^{-16}$)")
s = s.replace("$\\leq \\SI{3.3e-16}{}$", "$\\leq 3.3 \\times 10^{-16}$")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
print("siunitx refs left:", s.count("siunitx"), "| SI{ refs left:", s.count("\\SI{"))
