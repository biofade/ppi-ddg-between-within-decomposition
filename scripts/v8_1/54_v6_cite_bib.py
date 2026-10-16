# -*- coding: utf-8 -*-
"""scripts/v8_1/54_v6_cite_bib.py -- 把 v6 稿件的硬编码上标改为真正的 \\cite{}。

问题：初版生成器用了 \\textsuperscript{13} 之类硬编码上标，
      既不是 BibTeX 引用，也没有生成参考文献表（PDF 中 0 条）。
修复：建立 编号 -> BibTeX key 映射，替换为 \\cite{}，
      并由 unsrtnat + references.bib 生成文献表。
"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# v6 的文献顺序 1..13 -> BibTeX key
NUM2KEY = {
    1: "skempi2", 2: "leakage", 3: "cathddg", 4: "ddgpred", 5: "prefolddg",
    6: "flexddg", 7: "foldx5", 8: "proteinmpnn", 9: "esmif", 10: "ppiminer",
    11: "ndcg", 12: "bootstrap", 13: "biometry",
}


def to_cite(m):
    nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
    keys = [NUM2KEY[n] for n in nums if n in NUM2KEY]
    return r"\cite{" + ",".join(keys) + "}" if keys else m.group(0)


before = len(re.findall(r"\\textsuperscript\{[^}]*\}", s))
s = re.sub(r"\\textsuperscript\{([^}]*)\}", to_cite, s)
after = len(re.findall(r"\\textsuperscript\{[^}]*\}", s))
print(f"converted textsuperscript -> cite : {before} -> {after}")

# 主稿已有 \bibliographystyle/\bibliography；确认存在
assert "\\bibliographystyle{unsrtnat}" in s, "missing bibliographystyle"
assert "\\bibliography{references}" in s, "missing bibliography"
print("bibliography block present: yes")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
print("cite count now:", s.count("\\cite{"))
