# -*- coding: utf-8 -*-
"""scripts/v8_1/46_caption_fig1.py -- 用新版 (D) 描述替换 Figure 1 图注中的旧版描述。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

OLD = ('           "(D) Standardised capability profiles of matched ablation variants; axes are "\n'
       '           "standardised for visual comparison only.")')
NEW = ('           "(D) Retrospective matched ablation on the identical 3,351-row evaluation "\n'
       '           "universe. Pooled Pearson and within-target Macro Spearman are shown for the "\n'
       '           "mutation-only baseline (F1), static local geometry (F2), mutation × partner "\n'
       '           "perturbation (F4), and the full local-physics model (F8). Static local geometry "\n'
       '           "improves pooled correlation while leaving within-target ranking essentially "\n'
       '           "unchanged, whereas F4 and F8 improve both axes. The capability trajectory is "\n'
       '           "retrospective and is not intended as a causal or prospective model-selection "\n'
       '           "analysis.")')

n = s.count(OLD)
s = s.replace(OLD, NEW)
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print(f"Figure 1 caption (D) replaced x{n}; syntax OK")
print("old wording present:", "Standardised capability profiles of matched ablation" in s)
