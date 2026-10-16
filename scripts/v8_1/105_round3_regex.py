# -*- coding: utf-8 -*-
"""scripts/v8_1/105_round3_regex.py -- 用空白容忍正则补上 104 未命中的 3 处。"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []


def sub(tag, pat, new):
    global s
    s2, n = re.subn(pat, lambda m: new, s, count=1, flags=re.S)
    if n:
        s = s2; done.append(tag)
    else:
        miss.append(tag)


# 1) Abstract: D(n,K) 加 rank-linear surrogate 限定
sub("abstract-dnk",
    r"Fixed-\$K\$ shortlist metrics are constrained by candidate-set geometry: for \$K = 10\$,\s+"
    r"\$D\(n,K\)\$ is\s+zero when \$n \\leq K\$ and peaks near \$n\^\{\*\} \\approx 20\.5\$,",
    "Fixed-$K$ shortlist metrics are constrained by candidate-set geometry: for $K = 10$, a\n"
    "rank-linear discriminability surrogate $D(n,K)$ is zero when $n \\leq K$ and peaks near\n"
    "$n^{*} \\approx 20.5$,")

# 2) Results: 明确 D(n,K) 是 reference-dependent surrogate
sub("results-dnk2",
    r"which for \$K = 10\$ gives \$n\^\{\*\} \\approx 20\.5\$, so that \$D\(20\) = 0\.2381\$ while\s+"
    r"\$D\(50\) = 0\.1569\$ \(full derivation",
    "which for $K = 10$ gives $n^{*} \\approx 20.5$, so that $D(20) = 0.2381$ while\n"
    "$D(50) = 0.1569$. This surrogate is reference-dependent: only the degenerate case\n"
    "$n \\leq K \\Rightarrow D = 0$ holds for an arbitrary utility profile (full derivation")

# 3) Figure 4 caption
sub("fig4-caption",
    r"The discriminable interval \$D\(n,K\)\$ \(Eq\.~\\eqref\{eq:geometry\}\) is zero when\s+"
    r"\$n \\leq K\$ and peaks at",
    "Under the canonical rank-linear utility reference the discriminability\n"
    "  surrogate $D(n,K)$ (Eq.~\\eqref{eq:geometry}) is zero when $n \\leq K$ and peaks at")

p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("STILL MISSED:", miss)
