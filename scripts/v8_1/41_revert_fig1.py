# -*- coding: utf-8 -*-
"""scripts/v8_1/41_revert_fig1.py -- 把 Figure 1 改回原来的单张整图安放（正文宽、纵向页）。

撤销内容：
* 撤销横向页面（landscape section）
* 撤销「切成两半 + 2 列并排」的版式
* 恢复为 `addfig(D, "Fig1_overview_EN.png")` 单张插入

保留：
* 引用上标（19 处）
* Figure 4 的新图
"""
import io
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path("scripts/v8_1/31_manuscript_v6.py")
s = io.open(P, encoding="utf-8").read()

start = s.find("    # ---- Figure 1 放大安放")
if start == -1:
    start = s.find("    # Figure 1 原图宽高比")
end = s.find('    cap(D, "Figure 1 |')
if start == -1 or end == -1:
    print("markers not found:", start, end)
    raise SystemExit(1)

s = s[:start] + '    addfig(D, "Fig1_overview_EN.png")\n' + s[end:]
io.open(P, "w", encoding="utf-8").write(s)

import ast  # noqa: E402
ast.parse(s)
print("Figure 1 reverted to single full-width portrait insertion; syntax OK")
print("landscape section present:", "WD_ORIENT.LANDSCAPE" in s)
print("half-image reference present:", "Fig1_half" in s)

# 清理不再使用的派生半图
FIGS = Path("figures_v5")
for name in ("Fig1_half1_EN.png", "Fig1_half2_EN.png",
             "Fig1_row1_EN.png", "Fig1_row2_EN.png", "Fig1_row3_EN.png",
             "Fig1_split.txt", "Fig1_rows.txt"):
    p = FIGS / name
    if p.exists():
        p.unlink()
        print("  removed derived file:", name)
