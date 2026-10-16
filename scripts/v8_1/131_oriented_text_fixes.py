# -*- coding: utf-8 -*-
"""scripts/v8_1/131_oriented_text_fixes.py -- 方向对齐后的正文/表格/补充材料同步。

核心变化（来自 v8_1/EXTERNAL_ORIENTATION.json）：
  ESM-IF 方向对齐后：pooled +0.2250、macro -0.0213、C_B +0.1149、C_W +0.1101
  方法级 rank agreement：rho 0.90 -> 0.30、tau 0.80 -> 0.20、p 0.083 -> 0.683
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
O = json.loads(Path("v8_1/EXTERNAL_ORIENTATION.json").read_text(encoding="utf-8"))
RA, RB = O["raw"], O["oriented"]
applied, missed = [], []


def patch(name, pairs):
    p = Path("latex_v6") / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


BS = chr(92) * 2
M = []

# ---- Table 3：ESM-IF 行改为 aligned 值 ----
M.append((
    "table3-esmif",
    "ESM-IF & 1,776 & 66 & -0.2250 & 0.0213 & -0.1149 & -0.1101 & --- " + BS,
    "ESM-IF$^{\\dagger}$ & 1,776 & 66 & +0.2250 & -0.0213 & +0.1149 & +0.1101 & --- " + BS
))
M.append((
    "table3-note-align",
    "      kcal/mol is not defined for them.",
    "      kcal/mol is not defined for them. $^{\\dagger}$The published ESM-IF column is\n"
    "      oriented oppositely (larger values mean stronger binding); it is sign-oriented to the\n"
    "      common performance orientation before comparison (Methods; Supplementary Note~8)."
))

# ---- Results: rank agreement 数字 ----
M.append((
    "results-rank",
    "Rank ordering was strongly aligned between the two\n"
    "axes, with a single adjacent swap between FoldX and Flex-ddG, and DDGPred ranking first on\n"
    "both; method-level rank agreement is $\\rho = +0.90$ (Kendall $\\tau = +0.80$, exact\n"
    "permutation $p = 0.083$ for $n = 5$ methods). This strongly positive agreement does not\n"
    "support a universal decoupling of pooled and within-target performance and falls within the\n"
    "heterogeneous range observed internally. With only five methods the panel is underpowered,\n"
    "so we draw no directional conclusion from it.",
    "After sign-orienting every column to a common performance orientation (Methods), the two\n"
    f"rankings show no significant agreement: Spearman $\\rho = +{RB['rho']:.2f}$\n"
    f"(Kendall $\\tau = +{RB['tau']:.2f}$, exact permutation $p = {RB['perm_p']:.3f}$ for "
    f"$n = 5$\n"
    "methods). DDGPred ranks first on both axes, but ESM-IF moves from last on the pooled axis\n"
    "to second while remaining last on the within-target axis, so three of the five positions\n"
    "differ between the two rankings. Used raw, the same numbers would give\n"
    f"$\\rho = +{RA['rho']:.2f}$ ($p = {RA['perm_p']:.3f}$); the sign convention therefore "
    "determines\n"
    "whether the panel looks aligned or unaligned. With only five methods the panel is\n"
    "underpowered either way, so we draw no directional conclusion from it."
))

# ---- Abstract ----
M.append((
    "abstract-external",
    "External evidence\n"
    "comprised five methods and is underpowered and positively aligned, so it does not support a\n"
    "universal decoupling of pooled and within-target performance.",
    "External evidence comprised five methods and is underpowered: once all prediction columns\n"
    "are sign-oriented to a common performance orientation, the two axes show no significant\n"
    "rank agreement, so the panel neither supports nor refutes a universal decoupling of pooled\n"
    "and within-target performance."
))

# ---- Discussion ----
M.append((
    "discussion-external",
    "Across five independently developed methods the rankings were strongly positively aligned\n"
    "between pooled and within-target performance (Spearman $\\rho = +0.90$), with a single\n"
    "adjacent swap; this agreement neither supports nor refutes a universal decoupling of the two\n"
    "axes.",
    f"Across five independently developed methods, and after sign-orienting every column to a\n"
    f"common performance orientation, the rankings on the two axes show no significant agreement\n"
    f"(Spearman $\\rho = +{RB['rho']:.2f}$, exact permutation $p = {RB['perm_p']:.3f}$); ESM-IF "
    f"ranks\n"
    "second on the pooled axis but last on the within-target axis. The panel is too small to\n"
    "settle the question, so it neither supports nor refutes a universal decoupling of the two\n"
    "axes."
))

# ---- Methods: higher-is-better 一句改写 + 方向对齐小节 ----
M.append((
    "methods-higher-better",
    "All cross-metric comparisons use higher-is-better\norientation\\cite{biometry}.",
    "All cross-metric comparisons use higher-is-better orientation\\cite{biometry}, and every\n"
    "prediction column is sign-oriented to a common performance orientation before any\n"
    "cross-method comparison (see \\emph{External prediction semantics and units})."
))
M.append((
    "methods-orientation",
    "\\subsection*{Missing prediction and ambiguity handling}",
    "\\subsection*{Orientation alignment for cross-method comparison}\n"
    "Holding the raw-score provenance fixed is not sufficient for a fair performance comparison:\n"
    "a score whose sign convention is opposite to the response would be penalised or rewarded for\n"
    "its convention rather than for its predictions. We therefore determine, for each prediction\n"
    "column, whether larger values mean a larger $\\Delta\\Delta G$ (weaker binding), and align it\n"
    "to the response by sign only:\n"
    "\\begin{equation}\n"
    "  s^{\\mathrm{aligned}} = \\mu_c \\cdot s^{\\mathrm{raw}},\n"
    "  \\qquad \\mu_c \\in \\{+1, -1\\}.\n"
    "  \\label{eq:orientation}\n"
    "\\end{equation}\n"
    "This is a direction normalisation: it uses no evaluation labels, involves no recalibration,\n"
    "and changes no information content, since it only flips the sign of the correlation-like\n"
    "quantities. We set $\\mu_c = -1$ for ESM-IF, whose published score increases with binding\n"
    "strength (within-complex correlation with $\\Delta\\Delta G$ of $-0.24$), and $\\mu_c = +1$ for\n"
    "Flex-ddG, FoldX, DDGPred and ProteinMPNN (within-complex correlations $+0.49$, $+0.49$,\n"
    "$+0.51$ and $+0.30$ respectively). The raw columns are reported unmodified in\n"
    "Supplementary Table~S4 alongside each $\\mu_c$, so the provenance and the aligned values are\n"
    "both auditable.\n\n"
    "\\subsection*{External prediction mapping and ambiguity handling}"
))

patch("manuscript_v6.tex", M)

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
print(f"oriented: rho={RB['rho']:.4f} tau={RB['tau']:.4f} p={RB['perm_p']:.4f}")
