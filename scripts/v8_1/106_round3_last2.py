# -*- coding: utf-8 -*-
"""scripts/v8_1/106_round3_last2.py -- 补最后 2 处（Abstract D(n,K) 与 Figure 4 caption）。"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = Path("latex_v6/manuscript_v6.tex")
s = p.read_text(encoding="utf-8")
done, miss = [], []

A_OLD = ("Fixed-$K$ shortlist metrics are constrained by\n"
         "candidate-set geometry: for $K = 10$, $D(n,K)$ is zero when $n \\leq K$ and peaks near\n"
         "$n^{*} \\approx 20.5$,")
A_NEW = ("Fixed-$K$ shortlist metrics are constrained by\n"
         "candidate-set geometry: for $K = 10$, a rank-linear discriminability surrogate $D(n,K)$\n"
         "is zero when $n \\leq K$ and peaks near $n^{*} \\approx 20.5$,")
if A_OLD in s:
    s = s.replace(A_OLD, A_NEW, 1); done.append("abstract-dnk")
else:
    miss.append("abstract-dnk")

F_OLD = ("geometry.} The discriminable interval $D(n,K)$ is zero when $n \\leq K$ and follows\n"
         "  $K(n-K)/[n(n+1)]$ otherwise, peaking at $n^{*} = K + \\sqrt{K^{2} + K} \\approx 20.5$ for\n"
         "  $K = 10$.")
F_NEW = ("geometry.} Under the canonical rank-linear utility reference the discriminability\n"
         "  surrogate $D(n,K)$ is zero when $n \\leq K$ and follows $K(n-K)/[n(n+1)]$ otherwise,\n"
         "  peaking at $n^{*} = K + \\sqrt{K^{2} + K} \\approx 20.5$ for $K = 10$.")
if F_OLD in s:
    s = s.replace(F_OLD, F_NEW, 1); done.append("fig4-caption")
else:
    miss.append("fig4-caption")

p.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
