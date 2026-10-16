# -*- coding: utf-8 -*-
"""scripts/v8_1/27_fix_approx.py -- 修复 \\approx 被解释成 BEL 控制字符的问题。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "generate_publication_figures_v5.py"
s = io.open(P, encoding="utf-8").read()

# 任何裸 \a（近似、以及被吃掉的 approx）都补成双反斜杠
before = s
s = s.replace("\\approx", "\\\\approx")          # 已是 \\ 的不再变
s = s.replace("pprox20.5", "pprox20.5")           # 占位保序
# 处理真实残留：BEL(0x07) + "pprox"
s = s.replace("\x07pprox", "\\\\approx")
# 修复被写坏的双反斜杠（避免 over-escape）
s = s.replace("\\\\\\\\approx", "\\\\approx")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
try:
    ast.parse(s)
    print("changed:", s != before, "| syntax OK")
except SyntaxError as e:
    print("SyntaxError", e.lineno, e.text)
    raise SystemExit(1)

lines = s.split("\n")
for i, ln in enumerate(lines, 1):
    if "Shortlist geometry" in ln or "peak at" in ln:
        print(f"{i:5d}| {ln.strip()[:180]}")
        print("      has BEL:", "\x07" in ln, "| has approx:", "approx" in ln)
