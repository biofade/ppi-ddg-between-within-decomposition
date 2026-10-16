# -*- coding: utf-8 -*-
"""scripts/v8_1/124_external_semantics_methods.py -- 插入 Methods 小节 + 修正表编号。

  (a) Methods 新增 "External prediction semantics and units"（锚点用真实存在的小节）
  (b) 把 external provenance 表移到 Model pools 之后，使编号为 Table S3
      （Figure S3 与 Note 7 已分别引用 S2 / S3，此处保持一致）
"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
done, miss = [], []

# ---------- (a) Methods 小节 ----------
mp = Path("latex_v6/manuscript_v6.tex")
ms = mp.read_text(encoding="utf-8")
ANCHOR = "\\subsection*{Bootstrap and permutation inference}"
NEW = (
    "\\subsection*{External prediction semantics and units}\n"
    "The five external methods originate from a single published annotation table of\n"
    "SKEMPI-derived protein--protein interaction mutations, and each prediction column is used\n"
    "exactly as supplied: no sign flip, no rescaling and no recalibration against the evaluation\n"
    "labels. We verified this by recomputing every reported correlation directly from the raw\n"
    "columns, which reproduces the reported values to four decimals for all five methods\n"
    "(Supplementary Note~8). Three columns (Flex-ddG, FoldX and DDGPred) are quantitative\n"
    "$\\Delta\\Delta G$ predictions on a kcal/mol-like scale, whereas ESM-IF and ProteinMPNN are\n"
    "inverse-folding sequence models whose published outputs are dimensionless mutation-effect\n"
    "scores rather than energies. We therefore report RMSE only for the first group; for the\n"
    "latter two an error in kcal/mol is not defined and is shown as ``---''. Correlation-based\n"
    "quantities --- pooled Pearson, macro within-target Spearman, $C_{\\mathrm{B}}$ and\n"
    "$C_{\\mathrm{W}}$ --- are invariant to affine rescaling and are reported for all five\n"
    "methods.\n\n")
if "External prediction semantics and units" not in ms and ANCHOR in ms:
    ms = ms.replace(ANCHOR, NEW + ANCHOR, 1)
    mp.write_text(ms, encoding="utf-8")
    done.append("methods-external-units")
else:
    miss.append("methods-external-units")

# ---------- (b) 表编号：把 external 表移到 Model pools 之后 ----------
sp = Path("latex_v6/supplementary_v6.tex")
s = sp.read_text(encoding="utf-8")
m = re.search(r"\\begin\{landscape\}\\begin\{table\}\[h\].*?"
              r"external prediction columns.*?\\end\{landscape\}\n", s, re.S)
if m:
    tbl = m.group(0)
    s = s.replace(tbl, "", 1)
    # 插到 Model pools 表之后
    m2 = re.search(r"\\begin\{landscape\}\\begin\{table\}\[h\].*?"
                   r"rank-agreement sensitivity analysis\..*?\\end\{landscape\}\n", s, re.S)
    if m2:
        s = s[:m2.end()] + tbl + s[m2.end():]
        sp.write_text(s, encoding="utf-8")
        done.append("external-table-moved-after-pools")
    else:
        sp.write_text(s, encoding="utf-8")
        miss.append("pools-table-anchor")
else:
    miss.append("external-table-anchor")

print("fixed:", done)
if miss:
    print("MISSED:", miss)

# 报告当前表顺序
for i, mm in enumerate(re.finditer(r"\\caption\{([^}]{0,66})", s), 1):
    print(f"   caption {i}: {mm.group(1)}")
