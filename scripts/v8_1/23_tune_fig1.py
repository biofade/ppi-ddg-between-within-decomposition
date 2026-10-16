# -*- coding: utf-8 -*-
"""scripts/v8_1/23_tune_fig1.py -- 仅调整 Figure 1 中造成重叠的文字长度与位置。

不重设计绘图逻辑，只做：
* 缩短 Panel B / Panel C 中过长、会溢出容器的说明文字
* Panel A 底部图例说明移出容器下边界
"""
import io
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path("generate_publication_figures_v5.py")
s = io.open(P, encoding="utf-8").read()
orig = s

# 1) Panel B 说明缩短
s = s.replace("Target centring orthogonalises between-target means",
              "target centring separates between-target means")

# 2) Panel B 警示条缩短
s = re.sub(r"Common-scale additive covariance contributions[^\"]*",
           "Common-scale additive covariance contributions - NOT independent "
           "correlations or causal mechanisms", s)

# 3) Panel C 两条审计胶囊文字缩短
s = re.sub(
    r"Weighting Audit \(5 Schemes\):[^\"]*",
    r"Weighting audit: $C_{\\mathrm{B}}$ flips only near zero "
    r"(35/109; 0/109 once $|C_{\\mathrm{B}}|>0.005$); $C_{\\mathrm{W}}$ 100% sign-stable",
    s)
s = re.sub(
    r"Shortlist Geometry \(\$K = 10\$\):[^\"]*",
    r"Shortlist geometry ($K=10$): 40.8% of targets ($n\\geq5$) are degenerate "
    r"($n\\leq K$, $D=0$); peak at $n^*\\approx20.5$",
    s)

# 4) Panel A 底部说明移到容器外
s = re.sub(
    r"ax_a\.text\(0\.73, 0\.115, \"Purple bars",
    "ax_a.text(0.73, 0.038, \"Purple bars", s)
s = s.replace(
    "ha='center', va='center', fontsize=8.8, fontweight='bold', "
    "color=PALETTE['muted_ink'])\n\n    # ------",
    "ha='center', va='center', fontsize=8.2, fontweight='bold', "
    "color=PALETTE['muted_ink'])\n\n    # ------")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig)
for frag in ("target centring separates",
             "NOT independent correlations or causal mechanisms",
             "100% sign-stable", "are degenerate"):
    print(f"  present {frag!r}: {frag in s}")
print("syntax OK")
