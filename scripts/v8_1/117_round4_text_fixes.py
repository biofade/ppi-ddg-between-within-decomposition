# -*- coding: utf-8 -*-
"""scripts/v8_1/117_round4_text_fixes.py -- 第四轮文字修复。

  (a) Discussion: 不再把 rank-linear surrogate 推成一般规律
  (b) Abstract: "de-replicated" 措辞修正
  (c) Results: coverage 声明改为已完成的 common-intersection 验证
  (d) 新增 Supplementary Note 7：external coverage / common intersection
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
X = json.loads(Path("v8_1/EXTERNAL_COMMON_INTERSECTION.json").read_text(encoding="utf-8"))
VA = X["coverage_variable"]; VI = X["coverage_intersection"]
applied, missed = [], []


def patch(name, pairs):
    p = Path("latex_v6") / name
    s = p.read_text(encoding="utf-8")
    for tag, old, new in pairs:
        if old in s:
            s = s.replace(old, new, 1); applied.append(f"{name}:{tag}")
        else:
            missed.append(f"{name}:{tag}")
    p.write_text(s, encoding="utf-8")


M = []

# (a) Discussion 过度泛化
M.append((
    "discussion-dnk",
    "Because the discriminable interval first rises and then falls as candidate sets\n"
    "grow, enlarging a candidate set is not universally effective at improving discrimination. We\n"
    "recommend reporting $n$ alongside $K$, labelling $n \\leq K$ targets separately, or adopting\n"
    "fraction-based evaluation.",
    "Under the canonical rank-linear reference used here, the geometric discriminability\n"
    "surrogate first rises and then falls with candidate-set size, which illustrates that\n"
    "candidate-set size alone does not guarantee monotonically greater geometric separation\n"
    "under that reference. The reference-free conclusion is narrower but important: when\n"
    "$n \\leq K$, fixed-$K$ shortlist metrics cannot distinguish models at all, because every\n"
    "candidate is selected. We therefore recommend reporting $n$ alongside $K$, labelling\n"
    "$n \\leq K$ targets separately, or adopting fraction-based evaluation --- recommendations\n"
    "that follow from the degenerate case rather than from the surrogate."
))

# (b) Abstract de-replicated
M.append((
    "abstract-dereplicated",
    "and apply it to matched, de-replicated and provenance-filtered evaluation sets.",
    "and apply it to matched and provenance-filtered evaluation sets, with explicit sensitivity\n"
    "analyses for repeated mutation keys."
))

# (c) Results coverage claim -> 已完成验证
M.append((
    "results-coverage",
    "We therefore do not claim a fully matched panel;\n"
    "the coverage difference is at most 21 mutations (1.2\\% of the panel) and we verified that\n"
    "it does not drive the ranking conclusion.",
    f"We therefore do not claim a fully matched panel. The coverage difference is at most 21\n"
    f"mutations (1.2\\% of the panel), and we tested its influence directly: restricting all five\n"
    f"methods to their common intersection ({X['n_intersection_rows']:,} rows, "
    f"{X['n_intersection_targets']} targets,\n"
    f"{X['n_intersection_ge5']} with $n \\geq 5$) leaves the pooled ranking, the "
    f"within-target ranking, the\n"
    f"single adjacent swap and the rank agreement unchanged "
    f"($\\rho = +{X['rank_intersection']['rho']:.2f}$,\n"
    f"exact permutation $p = {X['rank_intersection']['perm_p']:.3f}$; Supplementary Note~7)."
))
patch("manuscript_v6.tex", M)

# (d) Supplementary Note 7
def row(df, c, src):
    return (f"{c.replace('_', chr(92) + '_')} & {int(df[c]['n_rows']):,} & "
            f"{int(df[c]['n_ge5'])} & {df[c]['pearson']:.4f} & {df[c]['macro']:.4f} & "
            f"{df[c]['C_B']:+.4f} & {df[c]['C_W']:+.4f} & {df[c]['rmse']:.4f}" + chr(92) * 2)


orderA = sorted(VA, key=lambda k: -VA[k]["pearson"])
orderB = sorted(VI, key=lambda k: -VI[k]["pearson"])
tblA = "\n".join(row(VA, c, "A") for c in orderA)
tblB = "\n".join(row(VI, c, "B") for c in orderB)

NOTE7 = r"""\section*{Supplementary Note 7 | External coverage and common intersection}
The five external methods are annotation columns of the same published SKEMPI table, mapped to
the annotation namespace through the normalised key of
\texttt{scripts/v8\_1/05\_external\_clean.py}. Each method has its own pattern of missing
entries, so the uniquely mapped subsets differ slightly in size: Flex-ddG maps 1{,}755 rows,
DDGPred 1{,}756, and FoldX, ESM-IF and proteinMPNN 1{,}776, out of 1{,}781 uniquely mappable
rows. Because all five are evaluated against the same targets, we can remove this imbalance
exactly by restricting every method to the rows for which \emph{all five} supply a prediction.

\begin{center}
\begin{tabular}{@{}lrrrrrrr@{}}
\toprule
Method & Rows & Targets ($n \geq 5$) & Pooled Pearson & Macro Spearman
  & $C_{\mathrm{B}}$ & $C_{\mathrm{W}}$ & RMSE \\
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(A) variable coverage --- each method's own uniquely mapped subset}} \\
""" + tblA + r"""
\midrule
\multicolumn{8}{@{}l@{}}{\emph{(B) common intersection --- rows predicted by all five methods}} \\
""" + tblB + r"""
\bottomrule
\end{tabular}
\end{center}

The common intersection retains __NI__ of the 1{,}781 uniquely mappable rows (__PCT__\,\%),
covering __NT__ targets of which __NG__ have $n \geq 5$. Comparing the two panels:

\begin{itemize}
  \item the pooled ranking is identical (DDGPred $>$ FoldX $>$ Flex-ddG $>$ proteinMPNN $>$
        ESM-IF);
  \item the within-target ranking is identical (DDGPred $>$ Flex-ddG $>$ FoldX $>$ proteinMPNN
        $>$ ESM-IF), so the single adjacent swap between FoldX and Flex-ddG persists;
  \item DDGPred ranks first on both axes under either treatment;
  \item rank agreement is unchanged at Spearman $\rho = +0.90$ (Kendall
        $\tau = +0.80$, exact permutation $p = 0.083$).
\end{itemize}

Individual metrics move only slightly --- for example the within-target rank of FoldX changes
from __FXA__ to __FXB__, and of ESM-IF from __ESA__ to __ESB__ --- which is why we report the
variable-coverage panel in the main text and this intersection as a sensitivity analysis. The
conclusion that the two axes are strongly positively aligned in this panel is therefore not an
artefact of unequal coverage. The panel remains small (five methods) and is underpowered
regardless of which treatment is used.
"""

NOTE7 = (NOTE7
         .replace("__NI__", f"{X['n_intersection_rows']:,}")
         .replace("__PCT__", f"{100*X['n_intersection_rows']/1781:.1f}")
         .replace("__NT__", str(X["n_intersection_targets"]))
         .replace("__NG__", str(X["n_intersection_ge5"]))
         .replace("__FXA__", f"{VA['ddG_FoldX']['macro']:.4f}")
         .replace("__FXB__", f"{VI['ddG_FoldX']['macro']:.4f}")
         .replace("__ESA__", f"{VA['ESM-IF']['macro']:.4f}")
         .replace("__ESB__", f"{VI['ESM-IF']['macro']:.4f}"))

p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")
anchor = "\\section*{Machine-readable supplementary data}"
if "Supplementary Note 7" not in s and anchor in s:
    s = s.replace(anchor, NOTE7 + "\n" + anchor, 1)
    p.write_text(s, encoding="utf-8")
    print("Note 7 inserted")
else:
    p.write_text(s, encoding="utf-8")
    print("Note 7 already present or anchor missing")

print(f"applied {len(applied)}; missed {len(missed)}")
if missed:
    print("MISSED:", missed)
print(f"intersection: {X['n_intersection_rows']} rows ({100*X['n_intersection_rows']/1781:.1f}%), "
      f"{X['n_intersection_targets']} targets, {X['n_intersection_ge5']} n>=5")
