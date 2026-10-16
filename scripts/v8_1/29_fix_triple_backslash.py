# -*- coding: utf-8 -*-
"""scripts/v8_1/29_fix_triple_backslash.py -- 修正过度转义产生的 \\\\\\approx。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "generate_publication_figures_v5.py"
s = io.open(P, encoding="utf-8").read()
before = s

# 源码级：把 3 个反斜杠+approx 收敛为 2 个（源文件里的 "\\approx"）
s = s.replace("\\\\\\\\\\\\approx", "\\\\approx")   # 6 -> 2 (字面)
s = s.replace("\\\\\\approx", "\\\\approx")          # 3 -> 2 (字面)

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != before, "| syntax OK")

for i, ln in enumerate(s.split("\n"), 1):
    if "Global Peak" in ln:
        print(f"{i:5d}| {ln.strip()[:150]}")
        print("      literal backslash count before 'approx':",
              ln.count("\\", 0, ln.find("approx")))
