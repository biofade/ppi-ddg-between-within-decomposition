# -*- coding: utf-8 -*-
"""scripts/v8_1/72_install_true_vectors.py -- 把绘图脚本产出的真矢量 PDF 装入论文。

来源：用户提供的 make_all_figures.py 运行结果（outputs/*.pdf，真矢量，无位图）。
目标：figures_v5/<论文用名>.pdf
备份：旧图（位图包装）移入 F:\\dsh\\cache\\tmp\\figs_raster_backup\\
"""
import shutil
import sys
from pathlib import Path

import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

GEN = Path(r"F:\dsh\cache\tmp\figuregen\outputs")
DST = Path("figures_v5")
BK = Path(r"F:\dsh\cache\tmp\figs_raster_backup")
BK.mkdir(parents=True, exist_ok=True)

MAP = {
    "Figure1_large_manuscript": "Fig1_overview_EN",
    "Figure2": "Fig2_controlled_split_collapse_EN",
    "Figure3": "Fig3_target_weighting_sensitivity_EN",
    "Figure4": "Fig4_fixed_K_geometry_EN",
    "Figure5": "Fig5_retrospective_remedy_EN",
    "Figure6": "Fig6_external_methods_EN",
    "FigureS1": "FigS1_matched_universe_strata_EN",
    "FigureS2": "FigS2_capability_heatmap_provenance_clean_EN",
    "FigureS3": "FigS3_pool_sensitivity_rank_agreement_EN",
}

print("  %-24s %-46s %-22s %s" % ("generated", "installed as", "PDF page pt", "内容"))
ok = 0
for stem, name in MAP.items():
    src = GEN / f"{stem}.pdf"
    if not src.exists():
        print("  %-24s 源缺失" % stem)
        continue
    # 备份旧图
    old = DST / f"{name}.pdf"
    if old.exists():
        shutil.copy2(old, BK / f"{name}.pdf")
    # 安装真矢量 PDF
    shutil.copy2(src, old)
    # 同步 SVG（矢量备份）与 PNG
    for ext in ("svg", "png"):
        s2 = GEN / f"{stem}.{ext}"
        if s2.exists():
            shutil.copy2(s2, DST / f"{name}.{ext}")

    d = pymupdf.open(str(old))
    pg = d[0]
    drw, img = len(pg.get_drawings()), len(pg.get_images(full=True))
    pw, ph = pg.rect.width, pg.rect.height     # close 前取出
    d.close()
    kind = "真矢量" if (drw > 10 and img == 0) else f"含 {img} 位图元素"
    print("  %-24s %-46s %-22s %s" % (
        stem, f"{name}.pdf", "%.1fx%.1f" % (pw, ph), kind))
    ok += 1

print(f"\n安装 {ok}/{len(MAP)} 个真矢量图")
print(f"旧图备份于 {BK}")
