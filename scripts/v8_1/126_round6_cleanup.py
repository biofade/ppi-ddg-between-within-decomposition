# -*- coding: utf-8 -*-
"""scripts/v8_1/126_round6_cleanup.py -- 第六轮终稿清理。

  1) ΔΔG 符号约定：从真实源确认（scripts/v5/10_p0_write_benchmark.py 的 sign 核实，
     以及 scripts/v8/03_model_registry.py 的 "ddG (positive = weaker binding)"）
  2) affine → positive affine rescaling
  3) 清除 Markdown * / ** 残留
  4) |ΔΔG|-like → ΔΔG-like
  5) Table S3 caption 去重复
  6) external 来源明确为 CATH-ddG 的 skempi_v2.csv，并去掉 official
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
applied, missed = [], []


def patch(name, pairs, regex=False):
    import re
    p = Path("latex_v6") / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if regex:
            s2, n = re.subn(old, lambda m: new, s, count=1, flags=re.S)
            if n:
                s = s2; applied.append(f"{name}:{tag}")
            else:
                missed.append(f"{name}:{tag}")
        elif old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


# ================= 主稿 =================
M = []

# (3) Markdown 残留 -> LaTeX
M.append(("md-bold", "**exactly for every individual configuration**",
          "\\textbf{exactly for every individual configuration}"))
M.append(("md-ital1", "*marginal\nmedians*", "\\emph{marginal\nmedians}"))
M.append(("md-ital2", "*summed\nabsolute component change*",
          "\\emph{summed absolute component change}"))
M.append(("md-ital3", "(Methods, *Structural-family definition*)",
          "(Methods, \\emph{Structural-family definition})"))

# (2) positive affine
M.append(("affine", "are invariant to affine rescaling and are reported for all five\nmethods.",
          "are invariant to positive affine rescaling (a positive change of scale and origin)\n"
          "and are reported for all five methods; a negative scaling would reverse their\n"
          "orientation and sign."))

# (1) ΔΔG 符号约定 —— 在 Methods 的 Dataset universes 首次定义数据处写明
M.append((
    "sign-convention",
    "\\emph{Controlled-split universe (7{,}523 rows).} The full set of forward-paired true\n"
    "$\\Delta\\Delta G$ observations,",
    "Throughout this work the analysed response is the binding-affinity change\n"
    "\\begin{equation}\n"
    "  y \\;=\\; \\Delta\\Delta G \\;=\\; \\Delta G_{\\mathrm{mut}} - \\Delta G_{\\mathrm{WT}},\n"
    "  \\qquad \\Delta G = 0.592\\,\\ln K_{\\mathrm{D}}[\\mathrm{M}],\n"
    "  \\label{eq:sign}\n"
    "\\end{equation}\n"
    "so that $y > 0$ denotes \\emph{weaker} binding and $y < 0$ \\emph{stronger} binding. This\n"
    "matches the SKEMPI v2.0 $\\Delta\\Delta G$ convention and the convention of the external\n"
    "annotation table used here, so no sign inversion is applied anywhere in the pipeline.\n"
    "Predictions are handled in the same orientation as they are published (Methods,\n"
    "\\emph{External prediction semantics and units}); we note that the published ESM-IF column\n"
    "is oriented oppositely, carrying a negative pooled correlation, which is a property of that\n"
    "score rather than of our processing.\n\n"
    "\\emph{Controlled-split universe (7{,}523 rows).} The full set of forward-paired true\n"
    "$\\Delta\\Delta G$ observations,"
))

# (6) external 来源写法
M.append((
    "external-source",
    "We built external evaluation sets from published official per-mutation predictions: Flex\n"
    "ddG\\cite{flexddg}, FoldX\\cite{foldx5}, ProteinMPNN\\cite{proteinmpnn},\n"
    "ESM-IF\\cite{esmif} and DDGPred\\cite{ddgpred}.",
    "We evaluated five published per-mutation prediction columns compiled in the SKEMPI-derived\n"
    "annotation table distributed with CATH-ddG\\cite{cathddg} (\\texttt{skempi\\_v2.csv}): Flex\n"
    "ddG\\cite{flexddg}, FoldX\\cite{foldx5}, ProteinMPNN\\cite{proteinmpnn},\n"
    "ESM-IF\\cite{esmif} and DDGPred\\cite{ddgpred}."
))

patch("manuscript_v6.tex", M)

# ================= 补充材料 =================
S = []

# (1) Note 8 符号方向修正（这是本轮发现的错误）
S.append((
    "note8-orientation",
    "\\emph{Orientation.} All five columns are used with higher values meaning stronger binding, the\n"
    "same orientation as the $-\\Delta\\Delta G$ convention adopted here, so no method required a sign\n"
    "flip. We note explicitly that ESM-IF carries a negative pooled correlation\n"
    "($-0.2250$): the sign is a property of that model's published score on this benchmark and not\n"
    "an artefact of how we oriented it.",
    "\\emph{Orientation.} The analysed response is $y = \\Delta\\Delta G = \\Delta G_{\\mathrm{mut}} -\n"
    "\\Delta G_{\\mathrm{WT}}$ (Methods), so $y > 0$ denotes weaker binding. We use every column\n"
    "exactly as published and apply no sign inversion to any method. Two consequences are worth\n"
    "stating. First, the Flex-ddG, FoldX and DDGPred columns are already oriented so that larger\n"
    "values mean weaker binding, matching $y$, and they correlate positively with it. Second, the\n"
    "published ESM-IF column is oriented the other way --- larger values mean stronger binding ---\n"
    "which is why it carries a negative pooled correlation ($-0.2250$). That negative sign is a\n"
    "property of the published ESM-IF score on this benchmark, not an artefact of our processing;\n"
    "we did not flip it, and we flag it here so that the sign is not misread."
))

# (4) |ΔΔG|-like -> ΔΔG-like
S.append(("abs-ddg", "quantitative $|\\Delta\\Delta G|$-like predictors whose columns span standard deviations of",
          "quantitative $\\Delta\\Delta G$-like predictors (signed binding-affinity changes) whose\n"
          "columns span standard deviations of"))

# (2) positive affine（补充材料）
S.append(("affine-sup",
          "--- are invariant to affine rescaling and are reported for\nall five methods.",
          "--- are invariant to positive affine rescaling and are reported for all five methods."))

# (3) Markdown **no** -> \textbf{no}
S.append(("md-no", "& **no** \\\\", "& \\textbf{no} \\\\"))

# (5) Table S3 caption 去重复
S.append(("tables3-cap",
          "\\caption{Provenance and score semantics of the external prediction columns (Table~S3).}",
          "\\caption{Provenance and score semantics of the external prediction columns.}"))

# (6) Note 8 来源明确 —— 用正文引文号，避免补充材料缺少 bibliography
S.append((
    "note8-source",
    "The five external methods are not separate downloads: they are five annotation columns of one\n"
    "published, SKEMPI-derived table of protein--protein interaction mutations, which we read\n"
    "directly (Table~S3).",
    "The five external methods are not separate downloads: they are five annotation columns of the\n"
    "SKEMPI-derived protein--protein interaction table distributed with CATH-ddG (main text\n"
    "reference~3; \\texttt{skempi\\_v2.csv}), which we read directly without modification\n"
    "(Table~S3)."
))

# (3) 其余 Markdown 残留（若有）
S.append(("md-generic", "**no**", "\\textbf{no}"))

patch("supplementary_v6.tex", S)

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
