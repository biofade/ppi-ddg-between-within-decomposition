# -*- coding: utf-8 -*-
"""scripts/v8_1/26_fix_bad_escape.py -- 把两行的单反斜杠补成双反斜杠（raw 字符串）。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "generate_publication_figures_v5.py"
lines = io.open(P, encoding="utf-8").read().split("\n")

BAD_W = '"-  Weighting audit: $C_{\\mathrm{B}}$ flips only near zero (35/109; 0/109 once $|C_{\\mathrm{B}}|>0.005$); $C_{\\mathrm{W}}$ 100% sign-stable",'
BAD_S = '"-  Shortlist geometry ($K=10$): 40.8% of targets ($n\\geq5$) are degenerate ($n\\leq K$, $D=0$); peak at $n^*\\approx20.5$",'

GOOD_W = (r'"-  Weighting audit: $C_{\mathrm{B}}$ flips only near zero '
          r'(35/109; 0/109 once $|C_{\mathrm{B}}|>0.005$); '
          r'$C_{\mathrm{W}}$ 100% sign-stable",')
GOOD_S = (r'"-  Shortlist geometry ($K=10$): 40.8% of targets ($n\geq5$) are degenerate '
          r'($n\leq K$, $D=0$); peak at $n^*\approx20.5$",')

nw = ns = 0
for i, ln in enumerate(lines):
    st = ln.strip()
    if st.startswith('"-  Weighting audit:'):
        lines[i] = ln[:len(ln) - len(ln.lstrip())] + GOOD_W
        nw += 1
    elif st.startswith('"-  Shortlist geometry'):
        lines[i] = ln[:len(ln) - len(ln.lstrip())] + GOOD_S
        ns += 1

out = "\n".join(lines)
io.open(P, "w", encoding="utf-8").write(out)
import ast  # noqa: E402
ast.parse(out)
print(f"fixed weighting-audit={nw}, shortlist={ns}; syntax OK")
