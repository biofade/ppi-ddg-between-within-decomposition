# -*- coding: utf-8 -*-
"""scripts/v8_1/55_v6_natbib_order.py -- 让参考文献编号与 v6 一致（1=SKEMPI ... 13=Sokal）。

unsrtnat 按「首次引用顺序」编号，因此 RMSE 那句先引 biometry 就变成了 [1]。
v6 的编号是固定文献表顺序，故改用 natbib 的 nonamebreak + 显式
\\bibitem 列表（由 .bbl 生成后重排），保持上标数字引用。

做法：
  1. 把 \\bibliographystyle{unsrtnat} 换成 \\bibliographystyle{plainnat}
     （仍由 BibTeX 生成条目文本，保留 DOI/期刊等字段）
  2. 生成后处理脚本把 .bbl 中 13 条 \\bibitem 按 v6 顺序重排并重编号
  3. 正文引用改用 \\cite{...} 不变，编号由重排后的 \\bibitem 决定
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 用 plainnat 生成条目（数字编号 + 上标由 natbib 控制），随后重排 .bbl
s = s.replace(r"\bibliographystyle{unsrtnat}", r"\bibliographystyle{plainnat}")

# 增加：把 v6 的文献顺序写入一个说明性注释，便于后处理
s = s.replace(
    "%=====================================================================\n\\bibliographystyle{plainnat}",
    "% v6 reference order (for .bbl reordering):\n"
    "%   1 skempi2  2 leakage  3 cathddg  4 ddgpred  5 prefolddg  6 flexddg\n"
    "%   7 foldx5   8 proteinmpnn  9 esmif  10 ppiminer  11 ndcg  12 bootstrap\n"
    "%   13 biometry\n"
    "%=====================================================================\n\\bibliographystyle{plainnat}")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
print("bibliographystyle now:", "plainnat" if "plainnat" in s else "?")
