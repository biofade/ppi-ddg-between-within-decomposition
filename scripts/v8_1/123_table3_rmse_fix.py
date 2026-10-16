# -*- coding: utf-8 -*-
"""scripts/v8_1/123_table3_rmse_fix.py -- 把 Table 3 中 ESM-IF / ProteinMPNN 的 RMSE 改为 ---。

依据 v8_1/EXTERNAL_PROVENANCE.csv：这两列是无量纲打分（raw sd 0.093 / 0.012），
RMSE(kcal/mol) 无定义。其余三列为 kcal/mol-like，保留。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []
BS = chr(92) * 2   # LaTeX 换行 \\
for name, val in (("ESM-IF", "2.3587"), ("ProteinMPNN", "2.3199")):
    old = f"{name} & "
    idx = s.find(old)
    if idx < 0:
        miss.append(name); continue
    end = s.find(BS, idx)
    line = s[idx:end]
    if val in line:
        newline = line.replace(val, "---")
        s = s[:idx] + newline + s[end:]
        done.append(name)
    else:
        miss.append(name)
p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
