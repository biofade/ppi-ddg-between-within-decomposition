# -*- coding: utf-8 -*-
"""scripts/v8_1/25_locate_bad_escape.py -- 找出单反斜杠导致 ParseException 的行。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
lines = io.open("generate_publication_figures_v5.py", encoding="utf-8").read().split("\n")
pat = re.compile(r'"[^"]*\$[^"]*\\[a-zA-Z][^"]*"')
found = 0
for i, line in enumerate(lines, 1):
    if "Weighting audit" in line or "Shortlist geometry" in line or "Common-scale" in line:
        found += 1
        # 单反斜杠（不是 \\ ）即为坏行
        bad = re.search(r'(?<!\\)\\[a-zA-Z]', line)
        print(f"{i:5d} {'BAD ' if bad else 'ok  '}| {line.strip()[:170]}")
print("candidate lines:", found)
