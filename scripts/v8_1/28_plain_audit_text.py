# -*- coding: utf-8 -*-
"""scripts/v8_1/28_plain_audit_text.py -- 把两条审计胶囊文字改为纯 ASCII，规避 mathtext 解析。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "generate_publication_figures_v5.py"
lines = io.open(P, encoding="utf-8").read().split("\n")

GOOD_W = ('"-  Weighting audit: C_B flips only near zero (35/109; 0/109 once '
          '|C_B| > 0.005);  C_W 100% sign-stable",')
GOOD_S = ('"-  Shortlist geometry (K = 10): 40.8% of targets (n >= 5) are degenerate '
          '(n <= K, D = 0);  peak at n* ~ 20.5",')

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
print(f"plain-text audit captions: weighting={nw}, shortlist={ns}; syntax OK")
for i, ln in enumerate(out.split("\n"), 1):
    if "Weighting audit" in ln or "Shortlist geometry" in ln:
        print(f"{i:5d}| {ln.strip()[:150]}")
