# -*- coding: utf-8 -*-
"""scripts/v8_1/111_verify_figs_v2.py -- 更稳的图形一致性核验。

对每张图，从 figures_v5 源 PDF 抽取若干**较长文本片段**（>=8 字符，去重、去纯刻度），
再检查这些片段是否出现在最终 PDF 中。片段命中率高 => PDF 内嵌的就是该图。
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

docs = {"manuscript": pymupdf.open("latex_v6/manuscript_v6.pdf"),
        "supplementary": pymupdf.open("latex_v6/supplementary_v6.pdf")}
blob = {k: re.sub(r"\s+", " ", "\n".join(pg.get_text() for pg in v)).lower()
        for k, v in docs.items()}

print("  %-46s %-6s %s" % ("figure", "score", "missed fragments (if any)"))
allok = True
for name, where in MAIN.items():
    t = re.sub(r"\s+", " ", pymupdf.open(f"figures_v5/{name}.pdf")[0].get_text()).lower()
    # 取较长的词串片段
    frags, seen = [], set()
    for m in re.finditer(r"[a-z][a-z0-9 \.\-\+%]{7,60}", t):
        f = m.group(0).strip()
        if len(f) >= 10 and not re.fullmatch(r"[\d\.\s%\-]+", f) and f not in seen:
            seen.add(f); frags.append(f)
    frags = frags[:12]
    if not frags:
        print("  %-46s %-6s %s" % (name, "n/a", "no text fragments (pure vector)"))
        continue
    hit = [f for f in frags if f in blob[where]]
    score = len(hit) / len(frags)
    ok = score >= 0.6
    allok &= ok
    miss = [f for f in frags if f not in blob[where]][:2]
    print("  %-46s %-6s %s" % (name, f"{score:.0%}", "OK" if ok else "MISS " + str(miss)))
for v in docs.values():
    v.close()
print()
print("  ALL FIGURES VERIFIED IN FINAL PDF" if allok else "  SOME FIGURES FAILED")
