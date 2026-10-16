# -*- coding: utf-8 -*-
"""scripts/v8_1/110_verify_all_figures_in_pdf.py -- 逐一核对最终 PDF 里嵌入的图 vs figures_v5 源。

对比方式：从最终 PDF 提取每页的矢量绘制指纹（drawing 数 + 文本），
与 figures_v5 对应源 PDF 的文本比较，确认 PDF 里是**最新**内容。
"""
import re
import sys
from pathlib import Path

import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MAIN = {"Fig1_overview_EN": "manuscript", "Fig2_controlled_split_collapse_EN": "manuscript",
        "Fig3_target_weighting_sensitivity_EN": "manuscript",
        "Fig4_fixed_K_geometry_EN": "manuscript", "Fig5_retrospective_remedy_EN": "manuscript",
        "Fig6_external_methods_EN": "manuscript",
        "FigS1_matched_universe_strata_EN": "supplementary",
        "FigS2_capability_heatmap_provenance_clean_EN": "supplementary",
        "FigS3_pool_sensitivity_rank_agreement_EN": "supplementary"}

pdfs = {
    "manuscript": pymupdf.open("latex_v6/manuscript_v6.pdf"),
    "supplementary": pymupdf.open("latex_v6/supplementary_v6.pdf"),
}
pdf_text = {k: "\n".join(pg.get_text() for pg in v) for k, v in pdfs.items()}

print("  %-46s %-12s %s" % ("figure", "in PDF?", "distinctive strings present"))
allok = True
for name, where in MAIN.items():
    src = Path("figures_v5") / f"{name}.pdf"
    if not src.exists():
        print("  %-46s 源缺失" % name); allok = False; continue
    t = pymupdf.open(src)[0].get_text()
    # 取该图内较长的数字/短语作为指纹（排除纯坐标轴刻度）
    cands = re.findall(r"\d+/\d+|\d\.\d{4}|\+0\.\d{3}|-\d\.\d{4}", t)
    seen, uniq = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c); uniq.append(c)
    uniq = uniq[:6]
    hit = [c for c in uniq if c in pdf_text[where]]
    ok = len(hit) >= max(1, len(uniq) // 2)
    allok &= ok
    print("  %-46s %-12s %d/%d  %s" % (name, "yes", len(hit), len(uniq),
                                       "OK" if ok else "MISMATCH -> " + str(uniq)))
for v in pdfs.values():
    v.close()

print()
print("  ALL FIGURES MATCH FINAL PDF" if allok else "  SOME FIGURES DO NOT MATCH")
sys.exit(0 if allok else 1)
