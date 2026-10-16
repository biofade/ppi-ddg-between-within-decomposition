# -*- coding: utf-8 -*-
"""scripts/v8_1/47_caption_fig1_v2.py -- 替换 Figure 1 图注 (D)（按实际换行匹配）。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

OLD = ('           "between-target component. (D) Standardised capability profiles of matched "\n'
       '           "ablation variants; axes are standardised for visual comparison only.")')
NEW = ('           "between-target component. (D) Retrospective matched ablation on the identical "\n'
       '           "3,351-row evaluation universe. Pooled Pearson and within-target Macro Spearman "\n'
       '           "are shown for the mutation-only baseline (F1), static local geometry (F2), "\n'
       '           "mutation × partner perturbation (F4), and the full local-physics model (F8). "\n'
       '           "Static local geometry improves pooled correlation while leaving within-target "\n'
       '           "ranking essentially unchanged, whereas F4 and F8 improve both axes. The "\n'
       '           "capability trajectory is retrospective and is not intended as a causal or "\n'
       '           "prospective model-selection analysis.")')

n = s.count(OLD)
s = s.replace(OLD, NEW)
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print(f"(D) replaced x{n}; syntax OK")
print("old (D) wording still present:",
      "Standardised capability profiles of matched" in s)
