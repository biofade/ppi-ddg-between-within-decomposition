# -*- coding: utf-8 -*-
"""scripts/v8_1/33_wire_addfig.py -- 把 v6 稿件的 add_picture 统一替换为 addfig。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

# 主稿 6 处
s = re.sub(
    r'D\.add_picture\(str\(FIGS / "(Fig[1-6]_[^"]+\.png)"\),\s*\n\s*width=Inches\(FIGW\["[^"]+"\]\)\)',
    lambda m: 'addfig(D, "%s")' % m.group(1), s)

# 补充材料 1 处（stem 变量）
s = re.sub(
    r'D\.add_picture\(str\(FIGS / stem[\s\S]{0,220}?width=Inches\(w\)\)',
    'addfig(D, f"{stem}.png")', s)

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("syntax OK")
print("addfig(D, ...) calls:", s.count("addfig(D,"))
print("remaining add_picture:", s.count("add_picture"))
for i, ln in enumerate(s.split("\n"), 1):
    if "add_picture" in ln or "addfig(D," in ln:
        print(f"{i:5d}| {ln.strip()[:120]}")
