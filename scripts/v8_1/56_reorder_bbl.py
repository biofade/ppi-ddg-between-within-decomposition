# -*- coding: utf-8 -*-
"""scripts/v8_1/56_reorder_bbl.py -- 把 BibTeX 生成的 .bbl 重排为 v6 的文献编号。

BibTeX 只能按 style 排序（unsrt=引用序 / plain=字母序），
而 v6 的文献表是固定顺序 1..13。natbib 会把编号写进 \\bibitem[...]{key}
的可选参数里，因此只要：
  1. 按 v6 顺序重排 \\bibitem 块
  2. 把每个块的可选参数 [..] 里的首字段改成新编号
就能得到与 v6 完全一致的编号，同时保留 BibTeX 生成的条目正文。
"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
V6 = Path("latex_v6")
BBL = V6 / "manuscript_v6.bbl"
ORDER = ["skempi2", "leakage", "cathddg", "ddgpred", "prefolddg", "flexddg",
         "foldx5", "proteinmpnn", "esmif", "ppiminer", "ndcg", "bootstrap",
         "biometry"]

txt = BBL.read_text(encoding="utf-8", errors="replace")

# 切出每个 \bibitem 块（到下一个 \bibitem 或 \end{thebibliography} 之前）
head, _, rest = txt.partition(r"\bibitem")
body, sep, tail = rest.rpartition(r"\end{thebibliography}")
blocks = re.split(r"(?=\\bibitem)", r"\bibitem" + body)
blocks = [b for b in blocks if b.strip().startswith(r"\bibitem")]
print("parsed bibitem blocks:", len(blocks))

by_key = {}
for b in blocks:
    m = re.search(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", b)
    if m:
        by_key[m.group(1)] = b
missing = [k for k in ORDER if k not in by_key]
extra = [k for k in by_key if k not in ORDER]
if missing or extra:
    print("  MISSING:", missing, " EXTRA:", extra)

# 重排 + 重编号
out = []
for n, key in enumerate(ORDER, start=1):
    b = by_key.get(key)
    if b is None:
        continue
    # 把可选参数 [..] 的首字段替换为 n
    def fix(m):
        inner = m.group(1)
        parts = inner.split("]", 1)
        first = parts[0]
        # 首字段形如 "Sokal and Rohlf(1995)"，只改开头的编号
        newfirst = re.sub(r"^\s*\d+\s*", "", first)
        rest2 = ("]" + parts[1]) if len(parts) > 1 else "]"
        return "[" + str(n) + newfirst + rest2
    b2 = re.sub(r"\[([^\]]*)\]", fix, b, count=1)
    out.append(b2.strip())

new_txt = (r"\begin{thebibliography}{99}" + "\n\n" + "\n\n".join(out) +
           "\n\n" + r"\end{thebibliography}" + "\n")
# 保证 \doi / \url 有定义（plainnat 的条目会用到 \doi）
PREAMBLE = ("\\providecommand{\\doi}[1]{doi: #1}\n"
            "\\providecommand{\\url}[1]{\\texttt{#1}}\n"
            "\\providecommand{\\href}[2]{#2}\n")
pre = txt.split(r"\bibitem")[0].strip()
if pre and not pre.startswith(r"\begin{thebibliography}"):
    new_txt = pre + "\n" + new_txt
# 前置宏定义必须在 \begin{thebibliography} 之前
new_txt = PREAMBLE + new_txt
BBL.write_text(new_txt, encoding="utf-8")
print("reordered .bbl ->", ORDER)
for i, k in enumerate(ORDER, 1):
    print(f"  [{i}] {k}")
