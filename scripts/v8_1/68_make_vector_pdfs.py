# -*- coding: utf-8 -*-
"""scripts/v8_1/68_make_vector_pdfs.py -- 用你提供的图生成 PDF 版图形文件。

背景：你给的图是位图（PNG/JPEG），v6 PDF 里也是位图。位图无法在不重跑生成代码的
前提下真正矢量化为路径，因此这里做的是 LaTeX 论文里标准的做法：
**每张图包成一个单页 PDF，页面尺寸严格等于图像宽高比**，再 \includegraphics 引用。
这样 ① 论文里的图形文件是 PDF；② 图形内容 = 你给的图，未重绘、未裁切、未改比例。

尺寸策略：按 600 dpi 等效，把图像像素映射为 PDF 点（1 pt = 1/72 in）。
    page_width_pt  = width_px  / 600 * 72
    page_height_pt = height_px / 600 * 72
"""
import sys
from pathlib import Path

import pymupdf
from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
FIGS = Path("figures_v5")

SRC = {
    "Fig1_overview_EN": r"F:\dsh-home\attachments\v1\files\00\00f0bdc6dce0f9a87cc9ab826286d0525de1adcd013f6ababfb677c26e1bbeae\figure1",
    "Fig2_controlled_split_collapse_EN": r"F:\dsh-home\attachments\v1\objects\45\450ff2fe302b6afd3947e0906a592b6e099e54c3a40e01b1a4fcd04b60bc17c3",
    "Fig3_target_weighting_sensitivity_EN": r"F:\dsh-home\attachments\v1\objects\66\66902ce0ca7a9b93e4cd59af439e942a246dc56a8ef4a8515deef712e2b51afe",
    "Fig4_fixed_K_geometry_EN": r"F:\dsh-home\attachments\v1\files\ff\ff942b146670da55db52fddf8b94bf35b42f7e76624618aa43e90d52ed9f7196\figure4",
    "Fig5_retrospective_remedy_EN": r"F:\dsh-home\attachments\v1\objects\67\6742ad4b421992695445ebe01321de4af535d3c12f21aa4a54b9f9339be0cbfc",
    "Fig6_external_methods_EN": r"F:\dsh-home\attachments\v1\objects\28\28aeb295edc0d6516628525ba8ca94ba406985d51a48acf43eacd963efb10399",
    "FigS1_matched_universe_strata_EN": r"F:\dsh-home\attachments\v1\objects\b9\b917a170b3334890da5a4424ab7d80d67bb4395e04d79c8034317aa04b97b030",
    "FigS2_capability_heatmap_provenance_clean_EN": r"F:\dsh-home\attachments\v1\objects\c2\c24657f1dcd9b01b74eca5bcb1ef41285eefabfcc879199d0c1a5ab4005ce8e8",
    "FigS3_pool_sensitivity_rank_agreement_EN": r"F:\dsh-home\attachments\v1\objects\25\258da5815eaebc7c329b747aefb15b7920d55407ed519e4a009e018f9e5b58e1",
}

DPI = 600.0
print("  按 600 dpi 等效建立单页 PDF（页面尺寸 = 图像宽高比）\n")
print("  %-46s %-13s %-18s %s" % ("figure", "pixels", "PDF page (pt)", "ratio"))
ok = 0
for name, src in SRC.items():
    p = Path(src)
    if not p.exists():
        print("  %-46s 源文件缺失" % name)
        continue
    im = Image.open(p)
    if im.mode not in ("RGB", "L"):
        im = im.convert("RGB")
    wpx, hpx = im.size
    pw, ph = wpx / DPI * 72.0, hpx / DPI * 72.0

    # 先把图像存成高质量的临时 PNG（无损，不改内容）
    tmp = Path("F:/dsh/cache/tmp") / (name + "_src.png")
    im.save(tmp, "PNG")

    doc = pymupdf.open()
    page = doc.new_page(width=pw, height=ph)
    page.insert_image(pymupdf.Rect(0, 0, pw, ph), filename=str(tmp), keep_proportion=False)
    out = FIGS / (name + ".pdf")
    doc.save(str(out), garbage=3, deflate=True)
    doc.close()

    d2 = pymupdf.open(str(out)); chk = d2[0].rect; d2.close()
    print("  %-46s %-13s %-18s px=%.4f pdf=%.4f" % (
        name, "%dx%d" % (wpx, hpx), "%.1fx%.1f" % (chk.width, chk.height),
        wpx / hpx, chk.width / chk.height))
    ok += 1
print(f"\n生成/替换 {ok}/{len(SRC)} 个图形 PDF")
