# -*- coding: utf-8 -*-
"""scripts/v8_1/147_seq40_keyfix_propagate.py -- 把 L3_seq40 键级修复后的数字传播到稿件。

背景：L3_seq40 原先因 `_seq40_cluster_units()` 未与 obs_key 并查集合并，
      导致 1,295-1,328 个 mutation key 跨 fold（key-level leakage），
      而冻结 artifact 里的 L3 恰好退化成与 L2 逐位相同。

修复后重算（scripts/v6_2/02_split4_ddg.py --levels L3_seq40 --seeds 42 43 44 --trees 800）：
      seq40  r=0.1961  C_B=+0.0239  C_W=+0.1721  macro=0.0631  rmse=2.6911
      对照 L2 complex r=0.1931  C_B=+0.0233  C_W=+0.1685  macro=0.0497  rmse=2.6412

**关键：随机档、Family-OOD 档与全部 delta（Δr、ΔC_B、ΔC_W、81.13% 归因）逐位未变。**
只有 seq40 那一行不再是 complex 的副本。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
done, miss = [], []

# ---------- (1) 生成器里的硬编码 Table 2 ----------
g = Path("scripts/v8_1/53_build_latex_v6.py")
gs = g.read_text(encoding="utf-8")
OLD_ROW = '    ("seq40", "0.1931", "+0.0233", "+0.1685", "0.0497", "2.6412"),'
NEW_ROW = '    ("seq40", "0.1961", "+0.0239", "+0.1721", "0.0631", "2.6911"),'
if OLD_ROW in gs:
    gs = gs.replace(OLD_ROW, NEW_ROW, 1)
    g.write_text(gs, encoding="utf-8")
    done.append("generator:table2-seq40")
else:
    miss.append("generator:table2-seq40")

# ---------- (2) 生成器里 "identical" 的 Table 2 表注 ----------
gb = Path("scripts/v8_1/53_build_latex_v6.py")
gbs = gb.read_text(encoding="utf-8")
for old, new in [
    ("\\item Median values across configurations. \\texttt{seq40} is identical to comp",
     "\\item Median values across configurations. \\texttt{seq40} is reported separately; comp"),
]:
    if old in gbs:
        gbs = gbs.replace(old, new, 1); done.append("generator:table2-note")
    else:
        miss.append("generator:table2-note")
gb.write_text(gbs, encoding="utf-8")

# ---------- (3) 生成产物 .tex ----------
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")

REPL = [
    # Table 2 行
    ('seq40 & 0.1931 & +0.0233 & +0.1685 & 0.0497 & 2.6412 ' + chr(92) * 2,
     'seq40 & 0.1961 & +0.0239 & +0.1721 & 0.0631 & 2.6911 ' + chr(92) * 2),
    # 正文：不再是 identical
    ("The\ncomplex-level and \\texttt{seq40} partitions produced identical assignments and were\n"
     "therefore not treated as independent comparisons.",
     "The\n\\texttt{seq40} partition is a distinct, sequence-identity-based split: after the "
     "key-level\ncorrection described in Methods it yields 301 units (versus 230 for the "
     "complex-level\nsplit), and it agrees with the complex-level assignment on only about 20\\% "
     "of rows. We\ntherefore report it separately rather than treating it as a duplicate of the "
     "complex-level\ncomparison."),
    # Figure 2 图注
    ("model specifications and seeds; complex and \\texttt{seq40} are identic",
     "model specifications and seeds; complex and \\texttt{seq40} are similar"),
    # Table 2 表注
    ("exactly; the identity holds per configuration. \\texttt{seq40} is ident",
     "exactly; the identity holds per configuration. \\texttt{seq40} is reported sep"),
]
for old, new in REPL:
    if old in s:
        s = s.replace(old, new, 1); done.append("tex:" + old[:34].replace(chr(10), " "))
    else:
        miss.append("tex:" + old[:34].replace(chr(10), " "))

p.write_text(s, encoding="utf-8")
print(f"applied {len(done)}; missed {len(miss)}")
if miss:
    print("MISSED:")
    for m in miss:
        print("   ", m)
