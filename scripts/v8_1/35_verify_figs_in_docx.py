# -*- coding: utf-8 -*-
"""scripts/v8_1/35_verify_figs_in_docx.py -- 校验 DOCX 内嵌图与 figures_v5 源图一致。"""
from __future__ import annotations

import hashlib
import sys
import zipfile
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(__file__).resolve().parent.parent.parent
V6 = R / "v8_1_final" / "manuscript_v6"
FIGS = R / "figures_v5"

EXPECT = {
    "CR-PPI-ddG_manuscript_v6_EN.docx": [
        "Fig1_overview_EN.png", "Fig2_controlled_split_collapse_EN.png",
        "Fig3_target_weighting_sensitivity_EN.png", "Fig4_fixed_K_geometry_EN.png",
        "Fig5_retrospective_remedy_EN.png", "Fig6_external_methods_EN.png"],
    "CR-PPI-ddG_Supplementary_v6_EN.docx": [
        "FigS1_matched_universe_strata_EN.png",
        "FigS2_capability_heatmap_provenance_clean_EN.png",
        "FigS3_pool_sensitivity_rank_agreement_EN.png"],
}


def sha(b):
    return hashlib.sha256(b).hexdigest()[:16]


allok = True
for docx, figs in EXPECT.items():
    p = V6 / docx
    z = zipfile.ZipFile(p)
    media = {n: z.read(n) for n in z.namelist() if n.startswith("word/media/")}
    print(f"\n{docx}  ({len(media)} embedded images)")
    hashes = {sha(v): k for k, v in media.items()}
    for f in figs:
        src = FIGS / f
        if not src.exists():
            print(f"  MISSING source {f}")
            allok = False
            continue
        h = sha(src.read_bytes())
        found = h in hashes
        print(f"  {'OK  ' if found else 'DIFF'} {f:50s} sha={h}"
              f"{'  -> ' + hashes[h] if found else ''}")
        if not found:
            allok = False

print("\nFIGURE_EMBED_MATCH =", "PASS" if allok else "FAIL")
