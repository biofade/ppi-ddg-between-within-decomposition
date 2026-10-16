# -*- coding: utf-8 -*-
"""scripts/v8_1/116_fig_attribution_wording.py -- 统一图内归因措辞为「summed absolute component change」。

涉及：
  Figure 1C 文字 "81.1% of the absolute component decline"
  Figure 2C 标题 "Attribution of the absolute decline"
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")
done, miss = [], []

REPL = [
    ('"81.1% of the absolute component decline"',
     '"81.1% of the summed absolute\\ncomponent change"'),
    (r'ax.set_title(r"C   Attribution of the absolute decline: $|\Delta C_B|/(|\Delta C_B|+|\Delta C_W|)$",',
     r'ax.set_title(r"C   Attribution of the summed absolute component change: '
     r'$|\Delta C_B|/(|\Delta C_B|+|\Delta C_W|)$",'),
]
for old, new in REPL:
    if old in s:
        s = s.replace(old, new, 1); done.append(old[:56])
    else:
        miss.append(old[:56])
P.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
