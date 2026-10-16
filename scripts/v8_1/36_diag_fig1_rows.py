# -*- coding: utf-8 -*-
"""scripts/v8_1/36_diag_fig1_rows.py -- 找出 Figure 1 中可切分的空白行，用于竖向重排放大。"""
import sys
import numpy as np
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
im = Image.open("figures_v5/Fig1_overview_EN.png").convert("L")
a = np.asarray(im)
H, W = a.shape
print(f"size = {W} x {H}  (aspect {W/H:.2f})")

# 每行中「非白」像素占比
rowdark = (a < 235).mean(axis=1)
# 找中部区域最白的一行
lo, hi = int(H * 0.25), int(H * 0.75)
mid = lo + int(np.argmin(rowdark[lo:hi]))
print(f"最空白的中部行 = {mid}  (占比 {rowdark[mid]:.4f})")

# 打印若干候选切分行的空白度
print("\ncandidate split rows (row, nonwhite fraction):")
for f in np.linspace(0.30, 0.70, 17):
    r = int(H * f)
    print(f"  y={r:5d}  frac={rowdark[r]:.4f}")

# 上下两半的高度比
for frac in (0.40, 0.44, 0.46, 0.48):
    r = int(H * frac)
    print(f"\nsplit at {frac:.2f} -> top {r}px, bottom {H-r}px; "
          f"top aspect {W/r:.2f}, bottom aspect {W/(H-r):.2f}")
