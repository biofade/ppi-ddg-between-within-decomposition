# -*- coding: utf-8 -*-
"""scripts/v8_1/32_install_supplied_figures.py -- 把用户提供的成图安装到 figures_v5/。

映射：
    Fig1_overview_EN.png                        <- ChatGPT 图像 ...14_13_46-1.png
    Fig2_controlled_split_collapse_EN.png       <- Figure2.png
    Fig3_target_weighting_sensitivity_EN.png    <- Figure3.png
    Fig4_fixed_K_geometry_EN.png                <- Figure4.png
    Fig5_retrospective_remedy_EN.png            <- Figure5.png
    Fig6_external_methods_EN.png                <- Figure6.png
    FigS1_matched_universe_strata_EN.png        <- FigureS1.png
    FigS2_capability_heatmap_provenance_clean_EN.png <- FigureS2.png
    FigS3_pool_sensitivity_rank_agreement_EN.png     <- FigureS3.png
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
FIGS = R / "figures_v5"
ATT = Path(r"F:\dsh-home\attachments\v1\objects")

MAP = [
    # (目标文件名, 附件 sha256 前缀目录, 附件文件名)
    ("Fig1_overview_EN.png",
     "79/79b00a3f63d9c81ab4b141b84a943b60d879bfc5d45e0b1cb9f550d73e20c707",
     "ChatGPT 图像 2026年10月16日 14_13_46-1.png"),
    ("Fig2_controlled_split_collapse_EN.png",
     "45/450ff2fe302b6afd3947e0906a592b6e099e54c3a40e01b1a4fcd04b60bc17c3",
     "Figure2.png"),
    ("Fig3_target_weighting_sensitivity_EN.png",
     "66/66902ce0ca7a9b93e4cd59af439e942a246dc56a8ef4a8515deef712e2b51afe",
     "Figure3.png"),
    ("Fig4_fixed_K_geometry_EN.png",
     "df/df50363aa9c87b8354873847248543bed897904e8ea1e80855e2cbc97e511c81",
     "Figure4.png"),
    ("Fig5_retrospective_remedy_EN.png",
     "67/6742ad4b421992695445ebe01321de4af535d3c12f21aa4a54b9f9339be0cbfc",
     "Figure5.png"),
    ("Fig6_external_methods_EN.png",
     "28/28aeb295edc0d6516628525ba8ca94ba406985d51a48acf43eacd963efb10399",
     "Figure6.png"),
    ("FigS1_matched_universe_strata_EN.png",
     "b9/b917a170b3334890da5a4424ab7d80d67bb4395e04d79c8034317aa04b97b030",
     "FigureS1.png"),
    ("FigS2_capability_heatmap_provenance_clean_EN.png",
     "c2/c24657f1dcd9b01b74eca5bcb1ef41285eefabfcc879199d0c1a5ab4005ce8e8",
     "FigureS2.png"),
    ("FigS3_pool_sensitivity_rank_agreement_EN.png",
     "25/258da5815eaebc7c329b747aefb15b7920d55407ed519e4a009e018f9e5b58e1",
     "FigureS3.png"),
]

ok = miss = 0
for target, sub, _name in MAP:
    src = ATT / sub
    if not src.exists():
        print(f"  MISSING: {sub}")
        miss += 1
        continue
    files = sorted(p for p in src.iterdir() if p.is_file()) if src.is_dir() else [src]
    if not files:
        print(f"  MISSING file in {sub}")
        miss += 1
        continue
    s = files[0]
    dst = FIGS / target
    # 统一转为 PNG（附件可能是 jpeg 但扩展名标为 png）
    im = Image.open(s)
    im = im.convert("RGB") if im.mode not in ("RGB", "RGBA") else im
    im.save(dst, "PNG")
    print(f"  OK  {target:52s} {im.size[0]}x{im.size[1]}  <- {s.name}")
    ok += 1

print(f"\ninstalled {ok}/{len(MAP)}, missing {miss}")
if miss:
    raise SystemExit(1)
