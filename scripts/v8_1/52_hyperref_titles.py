# -*- coding: utf-8 -*-
"""scripts/v8_1/52_hyperref_titles.py -- 章节标题中的数学改为 \\texorpdfstring。

原因：hyperref 生成 PDF 书签时会剥离数学，遇到 \\Delta\\Delta G 触发
"Improper alphabetic constant / \\mitDelta" 致命错误。
`\\texorpdfstring{排版}{书签纯文本}` 让 PDF 书签拿到纯文本，排版不受影响。
"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/49_build_latex.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 主稿：标题与副标题里的 $\Delta\Delta G$ / fixed-$K$ / Top-$K$
REPL = [
    (r"\\subsection{Pooled PPI $\Delta\Delta G$ performance contains separable between-target and\nwithin-target signal}",
     r"\\subsection{\texorpdfstring{Pooled PPI $\Delta\Delta G$ performance contains separable between-target and within-target signal}{Pooled PPI ddG performance contains separable between-target and within-target signal}}"),
    (r"\\section*{Supplementary Note 1 --- Derivation of the covariance decomposition}",
     r"\\section*{Supplementary Note 1 --- Derivation of the covariance decomposition}"),
]

# 通用做法：对 \section{...} / \subsection{...} 中含 $...$ 的标题自动包裹
def wrap(m):
    cmd, title = m.group(1), m.group(2)
    if "$" not in title:
        return m.group(0)
    plain = re.sub(r"\$[^$]*\$", lambda x: _plain(x.group(0)), title)
    title_flat = title.replace("\n", " ")
    return (r"\%s{\texorpdfstring{%s}{%s}}" % (cmd, title_flat, plain))


def _plain(math):
    body = math.strip("$")
    body = body.replace(r"\Delta\Delta G", "ddG")
    body = body.replace(r"\Delta", "D").replace(r"\rho", "rho")
    body = body.replace(r"\tau", "tau").replace(r"\mu", "mu")
    body = body.replace(r"\,", " ").replace("\\", "")
    body = body.replace("{", "").replace("}", "").replace("$", "")
    return body


# 只处理 \section{ / \subsection{ （带星号的不动，因为它们本就无编号书签问题）
s = re.sub(r"\\(section|subsection)\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}", wrap, s)

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
print("texorpdfstring count:", s.count("texorpdfstring"))
