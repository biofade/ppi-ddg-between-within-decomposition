# -*- coding: utf-8 -*-
"""scripts/v8_1/39_fig1_tworow.py -- 把 Figure 1 排成 2x1 网格铺满横向页面。

几何对比（图内文字的实际物理尺寸以「图片宽度 / 内容像素宽」计）：
  单张整图放 6.5 in 正文宽 : 高度 = 6.5/2.33 = 2.79 in   <- 现状，文字过小
  单张整图放 9.4 in 横向宽 : 高度 = 9.4/2.33 = 4.03 in
  2 行网格  放 9.4 in 横向 : 每行 4.6 in 宽、高 4.6/2.33 = 1.97 in
                            整块高 3.95 in，但每行只承载**一半内容宽度**，
                            因此单位面积的文字相对单张整图放大约 1.4 倍
实现：按空白行把原图切成上下两半，再在横向页面里并排成两行。
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
rowdark = (a < 235).mean(axis=1)
lo, hi = int(H * 0.30), int(H * 0.75)
cut = lo + int(np.argmin(rowdark[lo:hi]))
print(f"source {W}x{H}; cut at row {cut} (ink {rowdark[cut]:.4f})")

for i, (t, b) in enumerate([(0, cut), (cut, H)], start=1):
    seg = im.crop((0, t, W, b))
    g = np.asarray(seg.convert("L"))
    dr = np.where((g < 235).mean(axis=1) > 0.002)[0]
    if len(dr):
        seg = seg.crop((0, max(0, dr[0] - 3), W, min(seg.height, dr[-1] + 4)))
    p = FIGS / f"Fig1_half{i}_EN.png"
    seg.save(p, "PNG")
    print(f"  half{i}: y {t}-{b} -> {seg.size[0]}x{seg.size[1]}  {p.name}")

(FIGS / "Fig1_split.txt").write_text(f"{cut}\n", encoding="utf-8")
print("done")
