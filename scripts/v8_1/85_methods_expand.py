# -*- coding: utf-8 -*-
"""scripts/v8_1/85_methods_expand.py -- 扩充 Methods 至 benchmark 方法学论文所需覆盖度。

所有内容均有代码/文件依据：
  * artifact / configuration / method-family 三层 -> v8_1/MODEL_FAMILY_REGISTRY.csv
    列 artifact_id / key / method_family / configuration_id / seed / split_regime /
    training_provenance / rowset_hash / is_ablation / is_random_split_model /
    is_primary_representative / representative_reason
  * 结构族 -> canonical_target_labels.family_id_canonical (108 个)
  * 指纹 -> MATCHED_UNIVERSE_REGISTRY.rowset_hash (SHA-256 of sorted row ids)
  * 缺失预测处置 -> scripts/v8_1/05_external_clean.py (unmatched vs ambiguous)
  * 软件版本 -> 实测运行环境
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = Path("latex_v6")

METHODS_BLOCK = r"""\subsection*{Structural-family definition}
Structural families group complexes whose binding interfaces are homologous. They are
assigned by sequence-identity clustering of interface chains followed by union--find
transitive closure: two complexes join the same family if any of their chain clusters share
at least 40\% sequence identity. The annotation namespace comprises 108 structural families,
and family membership is frozen in the canonical target-label table used by every downstream
analysis. Structural families are the unit of the family-out-of-distribution split and are
unrelated to the \emph{method} families described below.

\subsection*{Prediction artifact hierarchy}
Prediction outputs are organised in three explicitly distinct levels, and every reported
comparison states the level at which it was computed.
\begin{itemize}
  \item \textbf{Artifact} --- one model-training run together with the row set it was
        evaluated on, identified by a row-set fingerprint. The registry contains 12
        artifacts spanning 7 distinct fingerprints.
  \item \textbf{Configuration} --- one trained model variant within an artifact, identified
        by a configuration key and a random seed. The full registry holds 202 configurations.
  \item \textbf{Method family} --- a group of configurations that share an architectural or
        feature-design lineage (for example the ten ablation families F0--F9 of the
        local-physics model). Families, not individual configurations, are the unit of
        independent evidence, because configurations within a family are not independent
        observations.
\end{itemize}
Cross-model comparison is only ever performed between configurations whose row-set
fingerprints are identical; metrics are never compared across fingerprints. The primary
matched universe contains 30 configurations drawn from 10 method families (three seeds per
family, with F0/F1/F2/F3/F4/F5/F6/F7/F8/F9 lineage labels).

\subsection*{Training-split provenance and known-overlap filtering}
Each configuration records its supervised training split regime as one of
\emph{random}, \emph{complex}, \emph{seq40} or \emph{family-OOD}. Configurations trained
under a permissive regime (random / complex / \texttt{seq40}) are excluded from all
deployment-oriented comparisons, because their evaluation rows may overlap their training
folds; they are retained only in the controlled split experiment, where the split regime is
itself the variable under study. Records of exclusions are retained in the leakage-exclusion
log shipped with the analysis. Throughout, we claim control only of \emph{task-specific
supervised} overlap; we cannot exclude the possibility that foundation-model pretraining
encountered similar sequences or structures.

\subsection*{Row-set fingerprinting}
Each prediction output receives a SHA-256 fingerprint computed over its sorted evaluation
row identifiers. The fingerprint is the sole basis for deciding that two configurations are
directly comparable, which makes the matching rule auditable and independent of any metric
value. Seven distinct fingerprints occur across the artifacts used here; Table~S1 lists each
stratum with its fingerprint, row count, target and family counts, and the number of models
available on it.

\subsection*{Primary matched universe selection}
The primary matched universe was selected by a performance-blind lexicographic rule fixed
before any metric values were inspected, ordering candidate strata by (1) provenance
cleanliness, (2) number of independent method families, (3) number of targets with
$n \geq 5$, and (4) number of rows. This yields the \texttt{f8c33c7e} stratum: 3{,}351 rows,
102 targets with $n \geq 5$, and 10 method families. Because the rule never references a
performance metric, the selection cannot be tuned towards a favourable result.

\subsection*{Representative selection within method families}
Analyses that require one entry per independent lineage use a single
\emph{performance-blind representative} per method family, chosen by configuration
identifier rather than by score. This rule is recorded per configuration in the model-family
registry, together with its stated reason. Reporting at the family level avoids treating
seed replicates as independent evidence, which would otherwise inflate the apparent number
of independent methods.

\subsection*{Missing prediction and ambiguity handling}
Two distinct failure modes are handled separately and are never conflated. A prediction is
\emph{unmatched} when its normalised matching key has no counterpart on the annotation side;
it is \emph{ambiguous} when the key maps to more than one annotation row. Of the 4{,}473
annotation rows, 644 are unmatched and 2{,}048 are ambiguous, leaving 1{,}781 uniquely
mapped rows. Ambiguous keys are excluded rather than forcibly resolved by an arbitrary
tie-break, because any deterministic tie-break would inject an unquantified bias. Each
external method is therefore evaluated on its own uniquely mapped subset, and coverage
differs slightly across methods (1{,}755--1{,}776 mutations).

\subsection*{Bootstrap and permutation inference}
Interval estimates for target-level quantities use a cluster bootstrap that resamples
targets, not rows, so that within-target dependence is preserved; 5{,}000 resamples are used
throughout with a fixed seed. For the controlled split, intervals are computed by resampling
the matched configuration pairs. For the five-method external panel, the exact permutation
distribution of Spearman's $\rho$ is enumerated over all $5! = 120$ relabellings, giving
$p = 0.083$; asymptotic approximations are not used because the panel is far too small for
them to be reliable.

\subsection*{Software and reproducibility}
Analyses were run under Python 3.11 with \texttt{numpy} 2.4, \texttt{pandas} 3.0,
\texttt{scipy} 1.17, \texttt{scikit-learn} 1.9 and \texttt{matplotlib} 3.11. No new model was
trained for this study and no GPU was used: all capability analyses re-analyse previously
produced prediction artifacts. Random seeds are fixed and recorded per configuration
(42/43/44 for the controlled split experiment; seed 42 for bootstrap resampling). Every
headline number is generated from a single machine-readable results file, and an automated
consistency checker re-derives the reported values from the underlying artifacts; the
checker, the frozen input manifest with hashes, and the machine-readable sources are
provided with the Supplementary Information.
"""

p = OUT / "manuscript_v6.tex"
s = p.read_text(encoding="utf-8")

# 插到 Statistical analysis 之后、External prediction mapping 之前
anchor = "\\subsection*{External prediction mapping and ambiguity control}"
if anchor in s:
    s = s.replace(anchor, METHODS_BLOCK.rstrip() + "\n", 1)
    # 旧的简短 ambiguity 小节并入上面的 Missing prediction 小节
    import re
    s = re.sub(r"\\subsection\*\{External prediction mapping and ambiguity control\}\n"
               r"Where the normalised matching key was not unique on the annotation side, only "
               r"uniquely\nmapped rows were used; ambiguous keys were excluded rather than "
               r"resolved\.\n", "", s, count=1)
    p.write_text(s, encoding="utf-8")
    print("Methods expanded: 9 subsections inserted before External prediction mapping")
    n = len(re.findall(r"\\subsection\*\{", s))
    print(f"Methods subsections now: {n}")
else:
    print("ANCHOR NOT FOUND — no change made")
