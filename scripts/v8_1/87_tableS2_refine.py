# -*- coding: utf-8 -*-
"""scripts/v8_1/87_tableS2_refine.py -- 修饰 Table S2 的表注与文本。

两点改进（内容仍全部来自真实文件）：
  1. 表注说明：primary universe 的 10 个方法族**不**出现在 109 配置的 weighting audit
     中，因为两者取自不同的 artifact 集合（V7.3 消融 vs V7.x rankers）。
  2. 把 representative-selection rule 的截断文本换成干净的表述。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")

# 1) 清理 rep_rule 截断文本
bad_variants = [
    "ABLATION - not collapsed (matched ablation set",
    "performance-blind, by configuration id",
]
import re  # noqa: E402
for b in bad_variants:
    s = s.replace(b, "performance-blind, by configuration id")

# 2) 强化表注
OLD_NOTE_START = "    \\item The primary matched universe comprises 30 configurations"
i = s.find(OLD_NOTE_START)
j = s.find("\\end{tablenotes}", i)
if i > 0 and j > i:
    NEW = ("    \\item The primary matched universe comprises 30 configurations drawn from "
           "these 10 method\n"
           "    families (three seeds per family). A \\emph{method family} groups "
           "configurations that share\n"
           "    an architectural or feature-design lineage; configurations within a family "
           "are not\n"
           "    independent observations, so family-level representatives are used wherever "
           "independent\n"
           "    evidence is required. All ten families share the same row-set fingerprint\n"
           "    (\\texttt{f8c33c7e}), which means they were evaluated on identical mutation "
           "rows and are\n"
           "    therefore directly comparable.\\par\n"
           "    The \\emph{in weighting audit} column reports whether a family's "
           "configurations also appear in\n"
           "    the 109-configuration weighting audit. For all ten primary families it is "
           "``no'': that audit\n"
           "    draws on ranker artifacts from the V7.x line rather than on these V7.3 "
           "ablation families, so\n"
           "    the two analyses use disjoint model sets. The weighting audit is therefore "
           "an independent\n"
           "    sensitivity analysis over a wider configuration pool, not a re-analysis of "
           "the primary\n"
           "    families; its 109 configurations are enumerated in "
           "\\texttt{WEIGHTING\\_EFFECT\\_SIZE\\_AUDIT.csv}\n"
           "    provided with the Supplementary Information.\n")
    s = s[:i] + NEW + s[j:]
    p.write_text(s, encoding="utf-8")
    print("Table S2 note strengthened")
else:
    print("note anchor not found; only rule text cleaned")
    p.write_text(s, encoding="utf-8")
