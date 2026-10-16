# -*- coding: utf-8 -*-
"""scripts/v8_1/66_compare2.py -- 正确的比对：两侧都渲染成同一尺寸的 RGB 数组。"""
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

W = 900          # 统一渲染宽度
print("  统一渲染到宽度 %d px 后比较（灰度）\n" % W)
print("  %-48s %-8s %-8s %s" % ("figure", "r", "meanAbsDiff", "判定"))
todo = []
for name, sp in SUP.items():
    img = Image.open(sp).convert("L")
    h = int(round(img.height * W / img.width))
    A = np.asarray(img.resize((W, h), Image.LANCZOS)).astype(float)

    d = pymupdf.open("figures_v5/%s.pdf" % name)
    pg = d[0]
    m = pymupdf.Matrix(W / pg.rect.width, h / pg.rect.height)
    pix = pg.get_pixmap(matrix=m, colorspace=pymupdf.csGRAY)
    B = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    B = B[:h, :W].astype(float)
    if B.shape != A.shape:
        B = np.asarray(Image.fromarray(B.astype(np.uint8)).resize((W, h),
                       Image.LANCZOS)).astype(float)
    d.close()

    r = np.corrcoef(A.ravel(), B.ravel())[0, 1]
    mad = np.abs(A - B).mean()
    verdict = "同一张图" if r > 0.985 else ("接近" if r > 0.95 else "内容不同")
    print("  %-48s %-8.4f %-8.2f %s" % (name, r, mad, verdict))
    if r <= 0.985:
        todo.append(name)

print()
if todo:
    print("需要替换的图（需重新生成矢量版本）：")
    for t in todo:
        print("   -", t)
else:
    print("全部 9 张与矢量 PDF 为同一张图，无需替换。")
