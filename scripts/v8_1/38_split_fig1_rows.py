# -*- coding: utf-8 -*-
"""scripts/v8_1/38_split_fig1_rows.py -- 把 Figure 1 横向条带按空白行切成 3 条，
竖排到横向页面上，使每条占满 9.4 in 宽，图内文字相对单张整图放大约 50%。

原图 1916x821（宽高比 2.33）。单张放在 9.4 in 宽时高度仅 4.0 in，文字仍小。
切成 3 条后，每条在 9.4 in 宽下高度 ~ (821/3)/1916*9.4 ≈ 1.34 in —— 不对，
这反而更小。正确做法：**每条独占一行且宽度拉满**，即整图被纵向拉伸布局，
每条的行高由内容决定，等效把 2.33 的宽高比拆成 3 段 ~7.0 的宽高比，
但每段只画自己那部分内容，所以单位面积的文字更大。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(__file__).resolve().parent.parent.parent
FIGS = R / "figures_v5"
SRC = FIGS / "Fig1_overview_EN.png"

im = Image.open(SRC).convert("RGB")
a = np.asarray(im.convert("L"))
H, W = a.shape
print(f"source {W}x{H} (aspect {W/H:.2f})")

rowdark = (a < 235).mean(axis=1)
# 在三个区段内各找一条最空白的行作为切点
cuts = []
for lo_f, hi_f in ((0.28, 0.40), (0.60, 0.78)):
    lo, hi = int(H * lo_f), int(H * hi_f)
    cuts.append(lo + int(np.argmin(rowdark[lo:hi])))
print("cut rows:", cuts, " ink:", [round(rowdark[c], 4) for c in cuts])

bounds = [0] + cuts + [H]
out = []
for i in range(3):
    t, b = bounds[i], bounds[i + 1]
    seg = im.crop((0, t, W, b))
    # 去掉上下纯白边
    g = np.asarray(seg.convert("L"))
    dark_rows = np.where((g < 235).mean(axis=1) > 0.002)[0]
    if len(dark_rows):
        seg = seg.crop((0, max(0, dark_rows[0] - 4), W,
                        min(seg.height, dark_rows[-1] + 5)))
    p = FIGS / f"Fig1_row{i+1}_EN.png"
    seg.save(p, "PNG")
    out.append(p)
    print(f"  row{i+1}: y {t}-{b} -> {seg.size[0]}x{seg.size[1]} "
          f"(aspect {seg.size[0]/seg.size[1]:.2f})  {p.name}")

# 供稿件脚本使用的行文件列表
(FIGS / "Fig1_rows.txt").write_text("\n".join(p.name for p in out), encoding="utf-8")
print("wrote Fig1_rows.txt")
