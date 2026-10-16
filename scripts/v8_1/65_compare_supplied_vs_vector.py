# -*- coding: utf-8 -*-
"""scripts/v8_1/65_compare_supplied_vs_vector.py -- 把你这次给的图与工程内矢量 PDF 逐一比对。

判定方法：把矢量 PDF 渲染成与你给的图相同像素尺寸，算相关系数。
  r > 0.99  -> 同一张图，矢量 PDF 已是其矢量版本，无需更换
  r 偏低     -> 内容不同，需要用新图重新生成矢量版本
"""
import sys
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SUP = {
    "Fig1_overview_EN": r"F:\dsh-home\attachments\v1\files\00\00f0bdc6dce0f9a87cc9ab826286d0525de1adcd013f6ababfb677c26e1bbeae\figure1",
    "Fig4_fixed_K_geometry_EN": r"F:\dsh-home\attachments\v1\files\ff\ff942b146670da55db52fddf8b94bf35b42f7e76624618aa43e90d52ed9f7196\figure4",
    "Fig2_controlled_split_collapse_EN": r"F:\dsh-home\attachments\v1\objects\45\450ff2fe302b6afd3947e0906a592b6e099e54c3a40e01b1a4fcd04b60bc17c3",
    "Fig3_target_weighting_sensitivity_EN": r"F:\dsh-home\attachments\v1\objects\66\66902ce0ca7a9b93e4cd59af439e942a246dc56a8ef4a8515deef712e2b51afe",
    "Fig5_retrospective_remedy_EN": r"F:\dsh-home\attachments\v1\objects\67\6742ad4b421992695445ebe01321de4af535d3c12f21aa4a54b9f9339be0cbfc",
    "Fig6_external_methods_EN": r"F:\dsh-home\attachments\v1\objects\28\28aeb295edc0d6516628525ba8ca94ba406985d51a48acf43eacd963efb10399",
    "FigS1_matched_universe_strata_EN": r"F:\dsh-home\attachments\v1\objects\b9\b917a170b3334890da5a4424ab7d80d67bb4395e04d79c8034317aa04b97b030",
    "FigS2_capability_heatmap_provenance_clean_EN": r"F:\dsh-home\attachments\v1\objects\c2\c24657f1dcd9b01b74eca5bcb1ef41285eefabfcc879199d0c1a5ab4005ce8e8",
    "FigS3_pool_sensitivity_rank_agreement_EN": r"F:\dsh-home\attachments\v1\objects\25\258da5815eaebc7c329b747aefb15b7920d55407ed519e4a009e018f9e5b58e1",
}

print("  %-48s %-13s %-13s %s" % ("figure", "supplied", "vector PDF", "correlation"))
worst = []
for name, sp in SUP.items():
    if not Path(sp).exists():
        print("  %-48s 源文件不存在" % name)
        continue
    img = Image.open(sp).convert("RGB")
    W, H = img.size
    d = pymupdf.open("figures_v5/%s.pdf" % name)
    pg = d[0]
    m = pymupdf.Matrix(W / pg.rect.width, H / pg.rect.height)
    pix = pg.get_pixmap(matrix=m, colorspace=pymupdf.csRGB)
    n = min(W * H * 3, len(pix.samples))
    X = np.frombuffer(img.tobytes()[:n], dtype=np.uint8).astype(float)
    Y = np.frombuffer(pix.samples[:n], dtype=np.uint8).astype(float)
    r = np.corrcoef(X, Y)[0, 1] if X.size == Y.size else float("nan")
    verdict = "同一张图" if r > 0.99 else ("接近" if r > 0.95 else "内容不同")
    print("  %-48s %-13s %-13s r=%.5f  %s" % (
        name, "%dx%d" % (W, H), "%.0fx%.0f" % (pg.rect.width, pg.rect.height), r, verdict))
    if r <= 0.99:
        worst.append((name, r, W, H))
    d.close()

print()
if worst:
    print("需要处理的图：")
    for n, r, W, H in worst:
        print("  %-48s r=%.5f  (%dx%d)" % (n, r, W, H))
else:
    print("全部 9 张均为同一张图 —— 工程内矢量 PDF 已是这些图的矢量版本")
