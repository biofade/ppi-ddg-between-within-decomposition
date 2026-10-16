# -*- coding: utf-8 -*-
"""scripts/v8_1/81_patch_fix.py -- 修正补丁自身引入的问题（保留既有补丁结果）。

修复项：
  1. 替换文本里误写成字面 ``\\n`` 而非真换行 -> LaTeX 报
     "There's no line here to end" / Undefined control sequence
  2. 补充材料 Row-universes 段落遗留未求值的占位符
     ``{mas['...']:,}`` -> 需用真实数值
  3. Figure S3 的 caption 里嵌了 landscape table -> "Not in outer par mode"
     改为把 table 移到 figure 之后
"""
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")
A = Path("v8_1")
UP = json.loads((A / "UNIVERSE_PROVENANCE.json").read_text(encoding="utf-8"))
mas = UP["universes"]["U_MAS_7523"]
can = UP["universes"]["U_CANON_4473"]

fixed = []

# ---- 1) manuscript: 把字面 \n 变成真换行 ----
p = OUT / "manuscript_v6.tex"
s = p.read_text(encoding="utf-8")
before = s.count("\\n")
# 只在明显的「行内连接」处替换：反斜杠+n 紧跟非字母 或 行首
s = s.replace("while\\n$D(50)", "while\n$D(50)")
s = s.replace("of $E_{\\max}$ and\\n$E_{\\mathrm{rand}}$",
              "of $E_{\\max}$ and $E_{\\mathrm{rand}}$")
s = s.replace("\\subsection*{Target-weighting schemes}\\nPredictions are unchanged",
              "\\subsection*{Target-weighting schemes}\nPredictions are unchanged")
p.write_text(s, encoding="utf-8")
fixed.append(f"manuscript literal-\\n: {before} -> {s.count(chr(92)+'n')}")

# ---- 2) supplementary: 替换未求值占位符 ----
p = OUT / "supplementary_v6.tex"
s = p.read_text(encoding="utf-8")
repl = {
    "{mas['n_unique_pdb_mutation_keys']:,}": f"{mas['n_unique_pdb_mutation_keys']:,}",
    "{mas['n_keys_with_multiple_rows']:,}": f"{mas['n_keys_with_multiple_rows']:,}",
    "{mas['n_keys_with_conflicting_ddg']:,}": f"{mas['n_keys_with_conflicting_ddg']:,}",
}
for k, v in repl.items():
    s = s.replace(k, v)
fixed.append("supplementary placeholders replaced: " +
             str({k: (k not in s) for k in repl}))
p.write_text(s, encoding="utf-8")

# ---- 3) supplementary: 把 landscape table 从 figure 内移到 figure 后 ----
p = OUT / "supplementary_v6.tex"
s = p.read_text(encoding="utf-8")
m = re.search(r"(\\begin\{landscape\}\\begin\{table\}\[h\].*?\\end\{landscape\}\n)", s, re.S)
if m:
    tbl = m.group(1)
    s = s.replace(tbl, "", 1)
    # 插到 Supplementary Figures 段之后（第一个 figure 之前）
    anchor = "\\section*{Supplementary Figures}"
    s = s.replace(anchor, anchor + "\n" + tbl, 1)
    fixed.append("moved pool table out of figure environment")
else:
    fixed.append("WARN pool table not found")

# caption 内改为指向 Table S4
s = s.replace("pool sizes, method-family counts and inclusion\n  criteria are listed in Table~S4.",
              "pool sizes, method-family counts and inclusion\n  criteria are listed in Table~S4 below.")
p.write_text(s, encoding="utf-8")

print(" | ".join(fixed))
