# -*- coding: utf-8 -*-
"""scripts/v8_1/45_install_fig1.py -- 安装新 Figure 1，并裁掉底部纯文字说明区。

新图路径：F:\\dsh-home\\attachments\\v1\\files\\00\\00f0bdc6...\\figure1
处理：
  1. 读入 -> 统一 PNG
  2. 检测底部是否存在「纯文字说明带」（原图中 caption 类文字），若存在则裁掉；
     判据：从下往上扫描，若若干连续行中「墨迹主要集中在很窄的高度带且无图形色块」，
     则视为说明文字。
  3. 为避免误裁图元，仅当底部存在明显空白分隔带时才裁；否则保留整图并报告。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(__file__).resolve().parent.parent.parent
FIGS = R / "figures_v5"
SRC = Path(r"F:\dsh-home\attachments\v1\files\00"
           r"\00f0bdc6dce0f9a87cc9ab826286d0525de1adcd013f6ababfb677c26e1bbeae"
           r"\figure1")
DST = FIGS / "Fig1_overview_EN.png"

if not SRC.exists():
    print("SOURCE NOT FOUND:", SRC)
    raise SystemExit(1)

im = Image.open(SRC)
if im.mode not in ("RGB", "RGBA"):
    im = im.convert("RGB")
W, H = im.size
print(f"source: {W}x{H}  (aspect {W/H:.2f})  mode={im.mode}")

g = np.asarray(im.convert("L"))
rowdark = (g < 240).mean(axis=1)

# 从底部往上找「最后一行有内容」，并在其上寻找一条接近全白的分隔带
last = H - 1
while last > 0 and rowdark[last] < 0.004:
    last -= 1
print(f"last content row = {last} (ink {rowdark[last]:.4f})")

# 在 last 之上 40% 范围内找最白的一行作为潜在分隔带
lo = int(last * 0.55)
seg = rowdark[lo:last]
if len(seg):
    sep = lo + int(np.argmin(seg))
    print(f"best separator row in bottom region = {sep} (ink {rowdark[sep]:.4f})")
    # 只有确实存在明显空白带，且其下方内容高度较小（< 22% 图高）才裁
    below_h = last - sep
    if rowdark[sep] < 0.02 and 8 < below_h < H * 0.22:
        crop = im.crop((0, 0, W, sep))
        crop.save(DST, "PNG")
        print(f"TRIMMED bottom text band: {below_h}px removed -> {crop.size[0]}x{crop.size[1]}"
              f"  (aspect {crop.size[0]/crop.size[1]:.2f})")
    else:
        im.save(DST, "PNG")
        print(f"NO trim (separator ink {rowdark[sep]:.4f}, band {below_h}px) -> saved as-is")
else:
    im.save(DST, "PNG")
    print("NO trim (no bottom region) -> saved as-is")

final = Image.open(DST)
print(f"installed {DST.name}: {final.size[0]}x{final.size[1]}")
