# -*- coding: utf-8 -*-
"""scripts/v8_1/136_round8_orientation_fix.py -- 第八轮：方向选择的诚实处理 + 7 项一致性修正。

核心决定（回应 BLOCKER）：
  上游 CATH-ddG 仓库（_audit/x/cath/CATH-ddG-main）**不包含**生成那五个打分列的代码，
  发布表中也没有任何 score-definition 文档。因此 o_m **无法独立于当前 evaluation labels
  确定**。原 Methods 声称 "uses no evaluation labels" 是事实错误。

  改为诚实口径：
    * RAW（未改动列）为 **primary** —— 无需任何主观决定，provenance 驱动；
    * sign-oriented 为 **exploratory sensitivity**，明确声明其方向由 evaluation
      correlation 判定，因此是 label-informed；
    * 两者结论相反（rho +0.90 vs +0.30），因此 external panel **不用于支持任何方向性结论**。

其它修正：Note 7 同步、Methods 的 permutation p、删除 "no sign flip" 矛盾、
Table S4 -> S3、Figure 6 图内标题、common-intersection 能量型 RMSE 恢复、mu_c -> o_m。
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


# ============================================================
# 主稿
# ============================================================
M = []

# ---- (1) Methods：方向对齐小节改写为诚实口径 + o_m 符号 ----
M.append((
    "methods-orientation",
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
    "both auditable.\n\n",
    "\\subsection*{Score orientation and its limits}\n"
    "The five external score columns are distributed as part of a published annotation table and\n"
    "we were unable to locate, in the accompanying repository or documentation, a statement of the\n"
    "score definition or sign convention used to produce them. We therefore treat the published\n"
    "columns as \\(\\emph{raw}\\) and make the raw analysis the primary external analysis: it requires\n"
    "no orientation decision and is fully provenance-driven.\n\n"
    "Because a score whose sign convention is opposite to the response would be penalised for its\n"
    "convention rather than for its predictions, we additionally report an \\emph{exploratory}\n"
    "sign-orientation sensitivity. Writing\n"
    "\\begin{equation}\n"
    "  s_m^{\\mathrm{aligned}} = o_m \\, s_m^{\\mathrm{raw}}, \\qquad o_m \\in \\{+1, -1\\},\n"
    "  \\label{eq:orientation}\n"
    "\\end{equation}\n"
    "we determine \\(o_m\\) from the sign of the within-complex correlation between column \\(m\\)\n"
    "and the evaluation response \\(y = \\Delta\\Delta G\\) (which removes between-target\n"
    "differences and so isolates direction rather than predictive skill). This gives\n"
    "\\(o_m = -1\\) for ESM-IF (within-complex correlation \\(-0.24\\)) and \\(o_m = +1\\) for\n"
    "Flex-ddG, FoldX, DDGPred and proteinMPNN (\\(+0.49\\), \\(+0.49\\), \\(+0.51\\), \\(+0.30\\)).\n\n"
    "We state the limitation of this procedure explicitly: \\textbf{the orientation multiplier is\n"
    "selected using the evaluation labels}, because no label-independent source for it could be\n"
    "found. It is a single sign bit per method, it involves no rescaling or recalibration, and it\n"
    "preserves all information content; but it is a label-informed decision and we therefore\n"
    "report it as a sensitivity analysis rather than as the primary result. The raw columns are\n"
    "preserved unmodified in Supplementary Table~S3, which also records \\(o_m\\). The two\n"
    "orientations lead to opposite conclusions about cross-axis agreement, and we report both.\n\n"
))

# ---- (3) Methods：删除 "no sign flip" 矛盾 ----
M.append((
    "methods-nosgnflip",
    "each prediction column is used\n"
    "exactly as supplied: no sign flip, no rescaling and no recalibration against the evaluation\n"
    "labels. We verified this by recomputing every reported correlation directly from the raw\n"
    "columns, which reproduces the reported values to four decimals for all five methods\n"
    "(Supplementary Note~8).",
    "each prediction column is ingested and preserved\n"
    "exactly as published, with no rescaling and no recalibration; the raw values are reproduced\n"
    "directly from the published columns to four decimals for all five methods (Supplementary\n"
    "Note~8). A sign multiplier is applied only for the exploratory orientation sensitivity\n"
    "described below, and never to the raw columns themselves."
))

# ---- (2) Methods：permutation p ----
M.append((
    "methods-perm",
    "For the five-method external panel, the exact permutation distribution of Spearman's $\\rho$\n"
    "is enumerated over all $5! = 120$ relabellings, giving\n"
    "$p = 0.083$;",
    "For the five-method external panel, the exact permutation distribution of Spearman's $\\rho$\n"
    "is enumerated over all $5! = 120$ relabellings, giving $p = 0.083$ for the raw columns and\n"
    f"$p = {RB['perm_p']:.3f}$ after sign-orientation;"
))

# ---- Results：改为 raw primary + sensitivity ----
M.append((
    "results-rank",
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
    "underpowered either way, so we draw no directional conclusion from it.",
    "On the raw columns the pooled and within-target rankings are strongly positively aligned:\n"
    f"Spearman $\\rho = +{RA['rho']:.2f}$ (Kendall $\\tau = +{RA['tau']:.2f}$, exact permutation\n"
    f"$p = {RA['perm_p']:.3f}$ for $n = 5$ methods), with DDGPred first on both axes and a "
    "single\n"
    "adjacent swap between FoldX and Flex-ddG. This agreement is however not robust to score\n"
    "orientation. Because the published columns carry no documented sign convention\n"
    "(Methods), we also evaluated an exploratory sign-oriented panel, in which each column is\n"
    "aligned to the response by a single sign bit determined from the evaluation data; there the\n"
    f"same rankings show no significant agreement\n"
    f"($\\rho = +{RB['rho']:.2f}$, $\\tau = +{RB['tau']:.2f}$, $p = {RB['perm_p']:.3f}$), because "
    "ESM-IF\n"
    "moves from last to second on the pooled axis while remaining last on the within-target axis.\n"
    "The external evidence therefore does not support a claim in either direction: with five\n"
    "methods, and with the orientation of the published scores undetermined by their own\n"
    "documentation, the panel is both underpowered and orientation-dependent, and we draw no\n"
    "directional conclusion from it."
))

# ---- Abstract：改为两向都不支持 ----
M.append((
    "abstract-external",
    "External evidence comprised five methods and is underpowered: once all prediction columns\n"
    "are sign-oriented to a common performance orientation, the two axes show no significant\n"
    "rank agreement, so the panel neither supports nor refutes a universal decoupling of pooled\n"
    "and within-target performance.",
    "External evidence comprised five methods and is underpowered; its apparent cross-axis\n"
    "agreement is not robust to score orientation, so the panel neither supports nor refutes a\n"
    "universal decoupling of pooled and within-target performance."
))

# ---- Discussion ----
M.append((
    "discussion-external",
    f"Across five independently developed methods, and after sign-orienting every column to a\n"
    f"common performance orientation, the rankings on the two axes show no significant agreement\n"
    f"(Spearman $\\rho = +{RB['rho']:.2f}$, exact permutation $p = {RB['perm_p']:.3f}$); ESM-IF "
    f"ranks\n"
    "second on the pooled axis but last on the within-target axis. The panel is too small to\n"
    "settle the question, so it neither supports nor refutes a universal decoupling of the two\n"
    "axes.",
    f"Across five independently developed methods the raw score columns suggest strong cross-axis\n"
    f"agreement (Spearman $\\rho = +{RA['rho']:.2f}$), but that agreement disappears under an\n"
    f"exploratory sign-orientation ($\\rho = +{RB['rho']:.2f}$). Because the published columns\n"
    "document no sign convention, and because the panel contains only five methods, the external\n"
    "evidence neither supports nor refutes a universal decoupling of the two axes; we report both\n"
    "orientations rather than selecting one."
))

# ---- (3) Table 3 表注：删除 no sign flip ----
M.append((
    "table3-note",
    "Prediction columns are used exactly as published,\n"
    "      with no sign flip and no recalibration (Supplementary Note~8). RMSE is reported only\n"
    "      for methods whose output is on a $\\Delta\\Delta G$-like scale;",
    "Raw prediction columns are ingested and preserved\n"
    "      exactly as published, with no rescaling or recalibration (Supplementary Note~8). A\n"
    "      sign multiplier is applied only for the orientation sensitivity described in Methods.\n"
    "      RMSE is reported only for methods whose output is on a $\\Delta\\Delta G$-like scale;"
))

# ---- (4) Table S4 -> S3 ----
M.append(("tables4-ref", "Supplementary Table~S4", "Supplementary Table~S3"))

patch("manuscript_v6.tex", M)

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
print(f"RAW primary: rho={RA['rho']:+.2f} p={RA['perm_p']:.3f} | "
      f"oriented sensitivity: rho={RB['rho']:+.2f} p={RB['perm_p']:.3f}")
