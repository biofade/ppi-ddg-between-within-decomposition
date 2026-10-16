# -*- coding: utf-8 -*-
"""scripts/v8_1/67_check_vector_vs_v6.py -- 判定 figures_v5 里的 PDF 是否就是你 v6 里的那些图。

ground truth = v6 PDF 内嵌的位图（就是你确认过的图）。
把 figures_v5/*.pdf 渲染到同样尺寸后比对：
  r 高 -> PDF 即该图的矢量版本
  r 低 -> PDF 是别的东西，需要用你给的图重新生成矢量 PDF
"""
import sys
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

V6M = Path(r"F:\dsh-home\attachments\v1\files\ca"
           r"\cae30b6694f6f910c1e5e852eac21ffadeaaf91f748ab7d6e176a901ae2ec443"
           r"\CR-PPI-ddG_manuscript_v6_EN.pdf")
V6S = Path(r"F:\dsh-home\attachments\v1\files\e9"
           r"\e93b741cd80f50cfeb2577a0fa2586c9e0a425457d9856f1a3a06ce9500ce0cc"
           r"\CR-PPI-ddG_Supplementary_v6_EN.pdf")

NAMES = ["Fig1_overview_EN", "Fig2_controlled_split_collapse_EN",
         "Fig3_target_weighting_sensitivity_EN", "Fig4_fixed_K_geometry_EN",
         "Fig5_retrospective_remedy_EN", "Fig6_external_methods_EN",
         "FigS1_matched_universe_strata_EN",
         "FigS2_capability_heatmap_provenance_clean_EN",
         "FigS3_pool_sensitivity_rank_agreement_EN"]


def v6_images(path):
    d = pymupdf.open(str(path))
    out = []
    for pg in d:
        for im in pg.get_images(full=True):
            p = pymupdf.Pixmap(d, im[0])
            if p.n > 4:
                p = pymupdf.Pixmap(pymupdf.csRGB, p)
            out.append(Image.frombytes("RGB", (p.width, p.height), p.samples))
    d.close()
    return out


v6 = v6_images(V6M) + v6_images(V6S)
print(f"v6 内嵌图总数: {len(v6)}")
for i, im in enumerate(v6):
    print(f"  #{i+1} {im.size[0]}x{im.size[1]}")

print()
print("  %-46s %-10s %s" % ("figures_v5 PDF", "ratioDiff", "与 v6 对应图的相关系数"))
W = 800
todo = []
for i, n in enumerate(NAMES):
    if i >= len(v6):
        break
    ref = v6[i]
    h = int(round(ref.height * W / ref.width))
    A = np.asarray(ref.resize((W, h), Image.LANCZOS).convert("L")).astype(float)

    d = pymupdf.open("figures_v5/%s.pdf" % n)
    pg = d[0]
    pw, ph = pg.rect.width, pg.rect.height      # 先取出，避免 close 后失效
    m = pymupdf.Matrix(W / pw, h / ph)
    pix = pg.get_pixmap(matrix=m, colorspace=pymupdf.csGRAY)
    B = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    d.close()
    B = np.asarray(Image.fromarray(B).resize((W, h), Image.LANCZOS)).astype(float)
    r = np.corrcoef(A.ravel(), B.ravel())[0, 1]
    rd = abs((ref.width / ref.height) - (pw / ph))
    flag = "OK" if r > 0.97 else "不一致"
    print("  %-46s %-10.4f r=%.4f  %s" % (n, rd, r, flag))
    if r <= 0.97:
        todo.append(n)

print()
print("需要重新生成矢量 PDF 的图：", todo if todo else "无")
