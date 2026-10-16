# -*- coding: utf-8 -*-
"""scripts/v8_1/34_fix_fig1_call.py -- 修正 Fig1 与补充材料遗留的 add_picture 调用。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

# 主稿 Fig1
s = s.replace(
    '    D.add_picture(str(FIGS / "Fig1_overview_EN.png"), '
    'width=Inches(FIGW["Fig1_overview_EN"]))',
    '    addfig(D, "Fig1_overview_EN.png")')

# 补充材料：循环里的 stem 已含 "_EN"，补 .png
s = s.replace('addfig(D, f"{stem}.png")', 'addfig(D, stem)')

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("syntax OK")
print("add_picture occurrences:", s.count("add_picture") - 1)  # 减去 helper 内定义
for i, ln in enumerate(s.split("\n"), 1):
    if "add_picture" in ln or "addfig(D," in ln:
        print(f"{i:5d}| {ln.strip()[:118]}")
