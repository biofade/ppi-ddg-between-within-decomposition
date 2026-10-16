# -*- coding: utf-8 -*-
"""scripts/v8_1/132_supp_orientation.py -- 补充材料：Note 7 修正 + Table S4 方向列 + Note 8 更新。

  1) Note 7 两套表的 ESM-IF / ProteinMPNN RMSE 全部改 `---`
  2) Note 7 "evaluated against the same targets" 措辞修正
  3) Table S4（external provenance）增加 raw orientation / multiplier / evaluation orientation
  4) Note 8 增加方向对齐说明（与 Methods 一致）
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
O = json.loads(Path("v8_1/EXTERNAL_ORIENTATION.json").read_text(encoding="utf-8"))
applied, missed = [], []
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")

# ---------- 1) Note 7 表中的 RMSE ----------
RMSE_OLD = ["2.3587", "2.3199", "2.3496", "2.3111", "2.1620", "2.1265", "2.3623"]
for v in RMSE_OLD:
    if v in s:
        s = s.replace(v, "---"); applied.append(f"note7-rmse-{v}")
    else:
        missed.append(f"note7-rmse-{v} (not present)")

# ---------- 2) "same targets" 措辞 ----------
OLD_T = ("Because all five are evaluated against the same targets, we can remove this imbalance\n"
         "exactly by restricting every method to the rows for which \\emph{all five} supply a "
         "prediction.")
NEW_T = ("Because all five columns originate from the same annotation universe and share common\n"
         "row identifiers, we can remove the coverage imbalance by restricting every method to\n"
         "their common row intersection.")
if OLD_T in s:
    s = s.replace(OLD_T, NEW_T, 1); applied.append("note7-same-targets")
else:
    missed.append("note7-same-targets")

# ---------- 3) Table S4 增加方向列 ----------
EP = __import__("pandas").read_csv("v8_1/EXTERNAL_PROVENANCE.csv")
MULT = O["multipliers"]
ORIENT = {"Flex-ddG": "+1", "DDGPred": "+1", "ddG_FoldX": "+1",
          "ESM-IF": "-1", "proteinMPNN": "+1"}
BS = chr(92) * 2
US = chr(92) + "_"
NEWROWS = []
for r in EP.itertuples():
    m = r.method
    mu = ORIENT.get(m, "+1")
    raw_or = ("higher = larger $\\Delta\\Delta G$ (weaker binding)" if mu == "+1"
              else "higher = stronger binding")
    rhs = "yes" if r.rmse_valid else chr(92) + "textbf{no}"
    NEWROWS.append(
        f"{m.replace('_', US)} & {chr(92)}texttt{{{r.raw_column.replace('_', US)}}} "
        f"& {r.units} & {raw_or} & ${mu}$ & {rhs} " + BS)

OLD_TBL_START = s.find("\\caption{Provenance and score semantics of the external prediction columns.}")
if OLD_TBL_START > 0:
    # 找到该 table 的 tabular 主体
    i0 = s.find("\\midrule", OLD_TBL_START)
    i1 = s.find("\\bottomrule", i0)
    if i0 > 0 and i1 > i0:
        s = s[:i0] + "\\midrule\n" + "\n".join(NEWROWS) + "\n    " + s[i1:]
        applied.append("tableS4-orientation-columns")
    else:
        missed.append("tableS4-body")
else:
    missed.append("tableS4-caption")

# 表头也要改（去掉 Transformation，加方向列）
OLD_HDR = ("Method & Raw column & Raw quantity & Units & Transformation applied here "
           "& RMSE meaningful? " + chr(92) * 2)
NEW_HDR = ("Method & Raw column & Units & Raw orientation & $\\mu_c$ & RMSE meaningful? "
           + chr(92) * 2)
if OLD_HDR in s:
    s = s.replace(OLD_HDR, NEW_HDR, 1); applied.append("tableS4-header")
else:
    missed.append("tableS4-header")

# ---------- 4) Note 8 增加方向对齐段 ----------
ANCH = "\\emph{Scope of the comparison.}"
NEW_SEC = (
    "\\emph{Sign orientation for cross-method comparison.} The published columns are not all in\n"
    "the same direction. Using within-complex correlations against $y = \\Delta\\Delta G$ --- which\n"
    "removes between-target differences and therefore isolates direction rather than predictive\n"
    "skill --- Flex-ddG, FoldX, DDGPred and ProteinMPNN are aligned with $y$ ($+0.49$, $+0.49$,\n"
    "$+0.51$ and $+0.30$), whereas ESM-IF is anti-aligned ($-0.24$): its published score is larger\n"
    "for stronger binding. For the cross-method performance comparison we therefore apply a\n"
    "sign-only orientation, $s^{\\mathrm{aligned}} = \\mu_c\\, s^{\\mathrm{raw}}$ with\n"
    "$\\mu_{\\mathrm{ESM\\text{-}IF}} = -1$ and $\\mu_c = +1$ for the other four (Methods). This uses no\n"
    "evaluation labels and no recalibration, and it changes no information content: it only flips\n"
    "the sign of the correlation-like quantities. Table~S4 records each raw column with its raw\n"
    "orientation and $\\mu_c$, so the unmodified provenance and the aligned values are both\n"
    "auditable. The effect on the panel is material: used raw the two axes appear strongly\n"
    "aligned ($\\rho = +0.90$, $p = 0.083$), whereas after orientation they show no significant\n"
    "agreement ($\\rho = +0.30$, $p = 0.683$). We report the oriented values throughout and state\n"
    "the raw alternative explicitly, since with five methods neither reaches adequate power.\n\n")
if "Sign orientation for cross-method comparison" not in s and ANCH in s:
    s = s.replace(ANCH, NEW_SEC + ANCH, 1); applied.append("note8-orientation")
else:
    missed.append("note8-orientation")

p.write_text(s, encoding="utf-8")
print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
