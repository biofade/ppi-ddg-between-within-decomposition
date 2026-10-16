# -*- coding: utf-8 -*-
"""scripts/v8_1/24_fix_fig1_text.py -- 修正上一步替换造成的 LaTeX 反斜杠丢失。"""
import io
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path("generate_publication_figures_v5.py")
s = io.open(P, encoding="utf-8").read()

# 找到被破坏的两行，用正确的 raw 字符串整体替换
bad_w = re.compile(r'"Weighting audit: \$C_\{\\mathrm\{B\}\}\$ flips[^"]*"')
bad_s = re.compile(r'"Shortlist geometry \(\$K=10\$\): 40\.8%[^"]*"')

good_w = (r'"Weighting audit: $C_{\mathrm{B}}$ flips only near zero '
          r'(35/109; 0/109 once $|C_{\mathrm{B}}|>0.005$); '
          r'$C_{\mathrm{W}}$ 100% sign-stable"')
good_s = (r'"Shortlist geometry ($K=10$): 40.8% of targets ($n\geq5$) are degenerate '
          r'($n\leq K$, $D=0$); peak at $n^*\approx20.5$"')

n1 = len(bad_w.findall(s))
n2 = len(bad_s.findall(s))
s = bad_w.sub(lambda m: good_w, s)
s = bad_s.sub(lambda m: good_s, s)

# 同时修复可能被破坏的警示条
s = re.sub(r'"Common-scale additive covariance contributions[^"]*"',
           r'"!  Common-scale additive covariance contributions — NOT independent '
           r'correlations or causal mechanisms"', s)

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
try:
    ast.parse(s)
    print("syntax OK")
except SyntaxError as e:
    print("SyntaxError line", e.lineno, ":", e.text)
    raise SystemExit(1)
print(f"replaced weighting-audit x{n1}, shortlist x{n2}")
for frag in (r"$C_{\mathrm{B}}$ flips only near zero",
             r"($n\leq K$, $D=0$)"):
    print(f"  present {frag!r}: {frag in s}")
