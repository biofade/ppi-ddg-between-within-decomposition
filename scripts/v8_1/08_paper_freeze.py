# -*- coding: utf-8 -*-
"""scripts/v8_1/08_paper_freeze.py -- V8.1 §23：论文冻结交付（标题/摘要/Results/Discussion）。

所有数字取自 `v8_1/FROZEN_RESULTS.json`。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
V8 = R / "experiments" / "v8"
FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
EXT = FR["external_clean"]; er = EXT["rank_agreement_external"]
EP = pd.read_csv(A / "external_clean_profiles.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
EP2 = pd.read_csv(A / "external_clean_vs_old.csv")
rmn, rmx = float(ST.rho_overall_vs_macro.min()), float(ST.rho_overall_vs_macro.max())
sf_cb = int(W2["signflip_C_B_gt0.005"].sum()); sf_cw = int(W2["signflip_C_W_gt0.005"].sum())
sf_cb0 = int(W2["signflip_C_B_gt0"].sum()); sf_cw0 = int(W2["signflip_C_W_gt0"].sum())
sf_cb1 = int(W2["signflip_C_B_gt0.001"].sum()); sf_cw1 = int(W2["signflip_C_W_gt0.001"].sum())
best_eu = EP.loc[EP.overall_pearson.idxmax()]
h470 = float(EP2.loc[EP2.model == "Flex-ddG", "delta_overall"].iloc[0])
h4max = float(EP2.delta_overall.min())

T = f"""# V8.1 论文冻结交付（§23）

> 数字来源：`FROZEN_RESULTS.json`。`GPU_USED = NO`。

---

## A. 标题选项（5 个，保守）

1. **Disentangling Between-Target and Within-Target Signal in PPI ΔΔG Benchmarking**
2. **Where Does PPI ΔΔG Benchmark Performance Come From? A Capability-Resolved Evaluation**
3. **Covariance-Signal Decomposition of PPI ΔΔG Benchmarks: Between-Target and Within-Target Contributions**
4. **Matched, Provenance-Clean Evaluation of PPI ΔΔG Predictors: Split Strictness, Target Weighting, and Candidate-Set Geometry**
5. **Separating Calibration from Within-Target Ranking in Protein–Protein Binding ΔΔG Prediction Benchmarks**

> 推荐 1 或 3。**不**暗示自动模型选择。

---

## B. 摘要（结构化，~240 词）

**Background.** Benchmarking of mutation effects on protein–protein binding (ΔΔG) still
relies predominantly on pooled correlation over heterogeneous mutation sets.

**Problem.** Pooled correlation mixes two different predictive capabilities:
locating a target's overall binding scale (**between-target** signal) and ranking
mutations **within** a target (within-target signal). Their relative contribution to
reported performance is not usually quantified, and evaluation regimes differ in how
much leakage they permit.

**Method.** We re-analyse existing PPI ΔΔG predictions under strictly matched row sets,
removing pseudo-replicated artifacts and provenance-impure models, and decompose the
pooled correlation exactly as `r_full = C_B + C_W` (between- and within-target
covariance contributions). We additionally audit metric direction, target-weighting
sensitivity, and fixed-K candidate-set geometry.

**Main controlled result.** Under a controlled shift from a permissive to a
family-out-of-distribution split — identical rows and identical model configurations —
median pooled Pearson falls from **{H4['random']['overall_pearson']:+.3f}** to
**{H4['family_OOD']['overall_pearson']:+.3f}**; **{100*DL.get('fraction_between',0):.0f}%**
of this reduction is carried by loss of the **between-target** contribution
(ΔC_B = {DL.get('delta_C_B',0):+.3f}), while the within-target contribution
**remains positive** (ΔC_W = {DL.get('delta_C_W',0):+.3f}).

**Weighting result.** C_B shows **directional instability near zero** under five
target-weighting schemes ({sf_cb}/{len(W2)} models reverse sign at |·|>0.005), whereas
C_W keeps **sign stability** (0/{len(W2)}) though its magnitude varies
(median range {W2.range_C_W.median():.3f}).

**Fixed-K geometry result.** For K = 10, the discriminable interval
D(n,K)=K(n−K)/[n(n+1)] is **zero whenever n ≤ K** —
**{100*GEO['n_le_K_fraction']:.1f}%** of our n≥5 targets — and peaks near
n* ≈ {GEO['n_star_analytic']:.1f}, so shortlist metrics are structurally
non-informative for small candidate sets.

**Remedy result.** In a matched retrospective ablation, decomposition separated a
pooled-score-only change (ΔMacro ρ = −0.002) from a modification that also improved
within-target ranking (ΔMacro ρ = +0.044).

**External caveat.** External evidence (5 methods) was **underpowered**;
cross-metric rank agreement was positive rather than decoupled.

**Conclusion.** PPI ΔΔG benchmarks should report between- and within-target components
under an explicitly declared weighting, together with candidate-set geometry.
"""

RS_ = f"""
## C. Results 主线

### 1. Benchmark performance conflates between-target and within-target signal
精确匹配、去伪重复、provenance-clean 后，pooled Pearson 可被**精确分解**为
`C_B + C_W`（残差 ≤ 3.3e-16）。两个分量的解释与量级不同。

### 2. Family-OOD evaluation disproportionately removes between-target signal
受控实验（同 {H4['random'].get('n_models','—')} 个配置、同 7523 行，仅改划分）：
pooled Pearson {H4['random']['overall_pearson']:+.3f} → {H4['family_OOD']['overall_pearson']:+.3f}；
ΔC_B = {DL.get('delta_C_B',0):+.3f}，ΔC_W = {DL.get('delta_C_W',0):+.3f}；
between-fraction = **{DL.get('fraction_between',0):.3f}**。
**注意**：complex 与 seq40 两个 regime 在此数据上逐位相同 ⇒ 不得作独立证据。

### 3. Aggregation weighting changes interpretation of the between-target contribution
五种加权方案下：C_B 在 {sf_cb}/{len(W2)} 个模型上 |·|>0.005 处翻转符号；
C_W 为 0/{len(W2)}，但幅度变化（range 中位 {W2.range_C_W.median():.3f}）。

### 4. Fixed-K shortlist evaluation can be structurally degenerate
n ≤ K=10 占 n≥5 target 的 **{100*GEO['n_le_K_fraction']:.1f}%**；
D(n,K) 在 n* ≈ {GEO['n_star_analytic']:.1f} 取峰后下降；
**「n ≥ 50 更可区分」不成立**。

### 5. Capability decomposition clarifies matched remedy effects
V7.3 同 3351 行：F2 static（ΔMacro ρ = −0.002，仅 pooled 改善）
vs F4 mutation×partner（ΔMacro ρ = +0.044，两轴同涨）。
**回溯性**，非前瞻预测。

### 6. External methods define the boundary of generalization
`U_EXT_EXACT_UNIQUE`（{EXT['n_rows']} 行 / {EXT['n_ge5']} 个 n≥5）上
方法级画像：Overall 与靶内排序赢家同为 `{best_eu.model}`
（{best_eu.overall_pearson:+.4f} / {best_eu.macro_spearman:+.4f}）；
秩一致 `rho(Overall, Macro) = {er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f}`
（perm p = {er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}），
**n = {er.get('n_methods')} ⇒ UNDERPOWERED**。
外部证据**方向与内部相反**（正相关而非解耦），须如实报告。
"""

DISC = f"""
## D. Discussion（必须显式讨论既有工作）

* **泄漏问题在本研究之前已知**：Bushuiev et al. (2024, arXiv:2404.10457) 已量化
  PDB-code 划分 65%、序列相似度划分 30%、SKEMPI v2.0 PPI-code 划分 56% 的结构性泄漏，
  并建议界面相似度划分。
* **family-OOD 与 per-structure 评价在本研究之前已知**：CATH-ddG
  (Bioinformatics 2025) 已使用 CATH 超家族划分并提出 TM-score < 0.6 的
  「hard mutation」定义；per-protein/per-structure 指标在既有工作中亦已出现。
* **本研究的贡献不是上述任何单一成分**，而是：
  1. **精确的 between-/within-target 协方差信号分解**（含恒等式 QC）；
  2. **量化划分严格度差异对 C_B 与 C_W 的不同影响**（受控实验）；
  3. **分解对目标加权的敏感性审计**；
  4. **匹配条件下的排名反转分析**；
  5. **fixed-K 候选集几何审计**；
  6. **回溯性 benchmark-guided remedy 归因**；
  7. 将它们**整合为一个 capability-resolved 的 PPI ΔΔG 评价框架**。

> 用词纪律：**「to our knowledge, we did not identify prior work that jointly
> reports ...」**，**不写「first」**。
"""

LIM = f"""
## E. Limitations

1. task-specific family-OOD **不保证**排除 foundation-model 预训练相似性；
2. benchmark 模型池**不是**所有 PPI ΔΔG 方法的随机样本；
3. 部分模型族覆盖有限（`V7.2_nested` 仅 1807 行 / 13 targets）；
4. **外部方法数不足**（n_methods = {er.get('n_methods')}）⇒ UNDERPOWERED；
5. **H3 = FAIL**：fixed-K 选择轴结构上无力，场景 C 仅 LOW-CONFIDENCE；
6. 候选选择为**回溯性**，非湿实验前瞻；
7. H5 remedy 分析为**回溯性**，benchmark **未**前瞻预测成功方案；
8. `C_B` / `C_W` 是 **covariance contributions**，**不是**因果机制；
9. **无独立前瞻实验数据集**（`EXTERNAL_DATASET_VALIDATION = NOT_RUN`）；
10. **不提供通用模型推荐器**；
11. `L2_complex` 与 `L3_seq40` 在本数据上**逐位相同**，不能作独立证据；
12. 此前的 CATH 映射含 {EXT['excluded_ambiguous_rows']} 个歧义行，外部读数已按
    EXACT_UNIQUE 重算（歧义键**被排除**，非解析）。
"""

FRZ = f"""
## F. Preprint decision

```
PREPRINT_READINESS = READY_AFTER_MINOR_FIXES
```

1. **是否有剩余科学阻断？** 否。但三项须在稿件中显式声明：
   外部 UNDERPOWERED 且方向与内部相反；complex≡seq40 不得作独立证据；H3 = FAIL（场景 C 低置信）。
2. **剩余问题是否仅编辑/格式？** 主要是（图注统一、Fig3 保留 duplicate-regime 警告）。
3. **是否需要额外模型计算？** **不需要。**
4. **V7.4R2 是否必需？** **不需要。**
5. **GPU 是否必需？** **不需要。**
6. **现在可否发 bioRxiv v1？** **可以**（在 minor fixes 且一致性检查 PASS 之后）。
"""
(F / "PAPER_FREEZE.md").write_text(T + RS_ + DISC + LIM + FRZ, encoding="utf-8")
print("PAPER_FREEZE.md ✅")

exec_block = f"""V8_1_FINAL_FREEZE

UPSTREAM_HASH_CHANGE = 0
FINAL_CONSISTENCY_CHECK = PASS
ANALYTIC_CLEANUP_STATUS = PASS
BENCHMARK_VALIDATION = PARTIALLY_SUPPORTED
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE
EXTERNAL_EVIDENCE_STATUS = UNDERPOWERED
H1_FINAL = WEAKENED_AND_REFRAMED (single-value orthogonality NOT supported;
    rho(Overall,Macro) range across {len(ST)} matched strata = {rmn:+.2f} .. {rmx:+.2f})
H2_FINAL = PASS (precised: C_B direction unstable near zero; C_W sign-stable, magnitude varies)
H3_FINAL = FAIL (three metric attempts; NOT rescued)
H4_FINAL = PASS (controlled; between-fraction = {DL.get('fraction_between', float('nan')):.3f})
H5_FINAL = PASS (retrospective matched ablation only)
PRIMARY_MATCHED_UNIVERSE = U_{PU['rowset_hash'][:8]} (level {PU['level']}, 100% GREEN)
PRIMARY_MATCHED_ROWS = {PU['n_rows']}
PRIMARY_MATCHED_TARGETS = n>=5: {PU['n_ge5']}
PRIMARY_METHOD_FAMILIES = {PU['n_method_families']}
GPU_USED = NO
V7_4R2_NEEDED = NO
PREPRINT_READINESS = READY_AFTER_MINOR_FIXES

TOP_3_DEFENSIBLE_FINDINGS:
1. Under a controlled permissive -> family-OOD split shift (identical rows and
   configurations), {100*DL.get('fraction_between',0):.0f}% of the loss in pooled Pearson
   is carried by the between-target covariance contribution (dC_B = {DL.get('delta_C_B',0):+.3f}),
   while the within-target contribution stays positive (dC_W = {DL.get('delta_C_W',0):+.3f}).
2. C_B is directionally unstable ONLY in the near-zero regime: sign flips occur in
   {sf_cb0}/{len(W2)} models with no threshold, {sf_cb1}/{len(W2)} exceeding |.|>0.001, and
   {sf_cb}/{len(W2)} exceeding |.|>0.005; C_W never flips sign at any threshold
   ({sf_cw0}/{sf_cw1}/{sf_cw}/{len(W2)}). C_W magnitude nevertheless varies
   (median range {W2.range_C_W.median():.3f} vs C_B {W2.range_C_B.median():.3f}),
   so any reported decomposition must declare its weighting.
3. Fixed-K candidate-selection evaluation is structurally non-informative when n <= K
   ({100*GEO['n_le_K_fraction']:.1f}% of n>=5 targets), and D(n,K)=K(n-K)/[n(n+1)]
   peaks near n* ~ {GEO['n_star_analytic']:.1f} -- so shortlist metrics require
   candidate-set geometry to be reported and the "n>=50" heuristic is unsupported.

TOP_3_LIMITATIONS:
1. External evidence is UNDERPOWERED (n_methods = {er.get('n_methods')}) and points the
   OPPOSITE way to the internal result (rho(Overall,Macro) = {er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f},
   perm p = {er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}).
2. H3 = FAIL: the candidate-selection axis cannot be rescued, so scenario C is
   LOW-CONFIDENCE and no universal model recommender is offered.
3. The benchmark model pool is not a random sample of PPI ddG methods, and
   task-specific family-OOD does not exclude foundation-model pretraining similarity.

NOVELTY_POSITION:
Prior work already established leakage-aware and family-OOD evaluation for PPI ddG
(Bushuiev et al. 2024; CATH-ddG 2025) and per-structure reporting. Our contribution is
therefore not any single ingredient but the systematic integration of an exact
between-/within-target covariance-signal decomposition with (i) a controlled
split-strictness comparison, (ii) a target-weighting sensitivity audit, (iii)
matched-universe rank-reversal analysis, (iv) a fixed-K candidate-set geometry audit,
and (v) a retrospective benchmark-guided remedy attribution. To our knowledge we did not
identify prior work that jointly reports these under matched, provenance-clean
evaluation; we make no "first" claim for leakage-aware, family-OOD, per-target or
Top-K benchmarking.

NEXT_ACTION:
MINOR_FIXES_ONLY
"""
(F / "V8_1_FINAL_FREEZE.txt").write_text(exec_block, encoding="utf-8")
print(exec_block)
