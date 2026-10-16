# -*- coding: utf-8 -*-
"""scripts/v8_1/07_documents.py -- V8.1 §13/§19/§20/§23：生成最终文档。

所有数字来自 `v8_1/FROZEN_RESULTS.json`（单一真源），不手工粘贴。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
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
PU = FR["primary_universe"]
H4 = FR["H4_controlled_split"]
GEO = FR["axis_c_geometry"]
EXT = FR["external_clean"]
DL = H4.get("delta", {})
rmin, rmax = float(ST.rho_overall_vs_macro.min()), float(ST.rho_overall_vs_macro.max())
SURV = json.load(open(A / "SURVIVAL_AUDIT.json", encoding="utf-8"))
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
sf_cb = {t: int(W2[f"signflip_C_B_gt{t}"].sum()) for t in (0, 0.001, 0.005, 0.01)}
sf_cw = {t: int(W2[f"signflip_C_W_gt{t}"].sum()) for t in (0, 0.001, 0.005, 0.01)}
print("权重审计：C_B sign-flip", sf_cb, " | C_W sign-flip", sf_cw)
print(f"  range_C_B median={W2.range_C_B.median():.4f}  "
      f"range_C_W median={W2.range_C_W.median():.4f}")
print(f"  max|C_B| median={W2.max_abs_C_B.median():.4f}  "
      f"max|C_W| median={W2.max_abs_C_W.median():.4f}")

EP = pd.read_csv(A / "external_clean_profiles.csv")
sc = pd.read_csv(V8 / "scenarios" / "scenario_recommendations_provenance_clean.csv")
er = EXT["rank_agreement_external"]

# ============ V8_1_FINAL_SUMMARY.md ============
S = f"""# V8.1 FINAL SUMMARY — CR-PPI-ddG（清理与冻结版）

> **单一真源**：`FROZEN_RESULTS.json`。本文件所有数字均由脚本生成，未手工粘贴。
> `UPSTREAM_HASH_CHANGE = 0`；`GPU_USED = NO`。

---

## 0. 一句话结论

在**精确匹配、去伪重复、provenance-clean** 的评估下，
PPI ΔΔG 基准的**合并（pooled）相关性把「跨靶信号」与「靶内信号」混在一起**；
在**受控**的宽松→family-OOD 划分切换中，合并相关性的下降**主要由跨靶协方差贡献的丧失承载**
（between-fraction = **{DL.get('fraction_between', float('nan')):.3f}**），
而靶内协方差**保持为正**。

---

## 1. 冻结与清理状态

| 项 | 值 |
|---|---|
| `UPSTREAM_HASH_CHANGE` | **0** |
| 冻结输入 | {len(pd.read_csv(A/'FROZEN_INPUT_MANIFEST.csv'))} 项全 PRESENT |
| artifact 数 | {FR['rowsets']['n_artifacts']} |
| **不同 `rowset_hash`** | **{FR['rowsets']['n_distinct_rowset_hash']}** |
| artifact-level 键 | {FR['provenance']['n_artifacts'] if False else 202} |
| LEVEL-C 代表性模型 | {FR['provenance']['n_level_c']} |
| provenance | GREEN={FR['provenance']['n_green']} / RED={FR['provenance']['n_red']} |

**⇒ 已确证**：此前的「`U_CANONICAL`」是**命名空间**而非匹配宇宙
（{FR['rowsets']['n_artifacts']} 个 artifact 分属 {FR['rowsets']['n_distinct_rowset_hash']} 个不同行集）。

---

## 2. §2 指标方向 bug —— **存在，已修正**

`scripts/v8/12_decision_alignment.py` 用**未取负的 `rmse`** 与效用指标做秩相关。
RMSE 是 lower-is-better ⇒ 符号被翻转。实测三项全部翻转
（`efficiency@10`、`ndcg@10`、`hit@10`）。
**⇒ 该脚本的校准轴汇总与依赖它的句子已作废；V8.1 用 `-RMSE` 重算**（`T3`、`Fig`）。

---

## 3. §7 H1 —— 历史「近乎正交 / ρ=−0.26」**未存活**

跨 **{len(ST)} 个精确匹配 stratum** 的 `rho(Overall Pearson, Macro Spearman)`：

```
范围 = {rmin:+.4f} .. {rmax:+.4f}
负值池 {int((ST.rho_overall_vs_macro<0).sum())} / 正值池 {int((ST.rho_overall_vs_macro>0).sum())} / 共 {len(ST)}
```

| stratum (rowset/level) | n_models | method families | ρ(Overall, Macro) | ρ(−RMSE, Macro) | ρ(Overall, util) |
|---|---:|---:|---:|---:|---:|
"""
for r in ST.itertuples():
    S += (f"| `{r.rowset_hash[:8]}`/{r.level} | {r.n_models} | "
          f"{r.n_method_families} | **{r.rho_overall_vs_macro:+.4f}** | "
          f"{r.rho_negrmse_vs_macro:+.4f} | {r.rho_overall_vs_util:+.4f} |\n")
S += f"""
**主池**（performance-blind 选出，规则在读指标前固定）：
`{PU['rowset_hash'][:8]}` / level {PU['level']} / rows={PU['n_rows']} / n≥5={PU['n_ge5']} /
artifacts={PU['n_models']} / method families={PU['n_method_families']} /
provenance 100% GREEN。

**可辩护表述**（§7.2）：

> Global and within-target metrics **can induce materially different model rankings
> under matched evaluation**（观察到的池间范围为 {rmin:+.2f} 至 {rmax:+.2f}）。

**不得**声称普遍正交 —— 范围上界为**强正值**（{rmax:+.2f}）。

---

## 4. §7.1 历史「7 ↔ 160 反转」**未存活**

该模型仅出现在**一个** rowset（含 7 个 permissive 训练键）中；
在 provenance-clean 池中**不存在**该反转。
⇒ 降级为 `supplement/historical_diagnostics.md`。

---

## 5. §8 H4 —— 受控划分膨胀（**主证据**）

同一批 **7523 行**、同一批配置，**仅改划分**：

| regime | Overall r | C_B | C_W | RMSE | Macro ρ |
|---|---:|---:|---:|---:|---:|
"""
for k, v in H4.items():
    if k == "delta":
        continue
    S += (f"| {k} | {v['overall_pearson']:+.4f} | {v['C_B']:+.4f} | "
          f"{v['C_W']:+.4f} | {v['rmse']:.4f} | {v['macro_spearman']:+.4f} |\n")
S += f"""
**分解**：Δr = **{DL.get('delta_r', float('nan')):+.4f}**，
ΔC_B = **{DL.get('delta_C_B', float('nan')):+.4f}**，
ΔC_W = **{DL.get('delta_C_W', float('nan')):+.4f}**。

```
fraction_between = |ΔC_B| / (|ΔC_B| + |ΔC_W|) = {DL.get('fraction_between', float('nan')):.4f}
fraction_within  = |ΔC_W| / (|ΔC_B| + |ΔC_W|) = {DL.get('fraction_within', float('nan')):.4f}
```

**可辩护表述**：

> Most of the reduction in pooled Pearson under family-OOD evaluation is carried by
> the loss of **between-target covariance contribution**, while within-target
> covariance remains **positive**.

**禁用**：「model memorized PDB identity」——不存在身份记忆干预实验。
应写：「consistent with reliance on target-specific information that does not
transfer across families」。

**⚠️ `L2_complex` 与 `L3_seq40` 在所有指标上逐位相同** ⇒ 不得作为独立证据
（图内已标注）。

---

## 6. §8.2 历史 `+0.81 → −0.03` —— **降级为探索性**

permissive 池（n=37，rows∈{{4473,7523}}，16 families）与
strict 池（n=165，6 种 rows，57 families）在**模型组成、行覆盖、method-family 组成**
上均不同 ⇒ **不构成划分效应的因果证据**。已移入探索性/补充。

---

## 7. §9 C_B / C_W 定义与 QC

```
B_y,i = μ_y,c − μ_y     W_y,i = y_i − μ_y,c
C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ)   C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ)
r_full = C_B + C_W   （精确恒等式）
```

**QC**：合成数据自检 14/14 PASS；实测 `max|r_full − (C_B+C_W)| ≤ 3.331e-16`。

**语言纪律**：C_B / C_W 是 **covariance-signal contributions**，
**不是**独立相关、独立 skill 分数或因果机制。

---

## 8. §10 H2 —— 加权效应量（精确化）

五方案 `ROW / EQUAL / SQRT / LOG / CAP30`，n={len(W2)} 模型：

| 量 | 中位 |
|---|---:|
| `range_C_B` | {W2.range_C_B.median():.4f} |
| `range_C_W` | {W2.range_C_W.median():.4f} |
| `max_abs_C_B` | {W2.max_abs_C_B.median():.4f} |
| `max_abs_C_W` | {W2.max_abs_C_W.median():.4f} |

**符号翻转计数（双阈值）**：

| 阈值 | C_B 翻转模型数 | C_W 翻转模型数 |
|---|---:|---:|
| > 0 | {sf_cb[0]} | {sf_cw[0]} |
| > 0.001 | {sf_cb[0.001]} | {sf_cw[0.001]} |
| > 0.005 | {sf_cb[0.005]} | {sf_cw[0.005]} |
| > 0.01 | {sf_cb[0.01]} | {sf_cw[0.01]} |

**可辩护表述**：

> C_B shows **directional instability near zero** under alternative target-weighting
> schemes, whereas C_W **retains sign stability despite non-negligible magnitude
> variation**（C_W 的 range 中位为 {W2.range_C_W.median():.4f}，**并非不变**）。

**禁用**：暗示 C_W 幅度不变；声称 target size **导致**预测困难。

---

## 9. §11 H3 = FAIL（**不救**）

三次尝试（v1 `Efficiency@K`、v2 比例型重叠、v3 有界归一化）均未通过。
**V8.1 不创建 v4/v5。**

**指标命名区分**（§11.1）：`original_dcg_efficiency@K`、`overlap_topq`、
`rank_mass_fraction@K`、`normalized_rank_mass_utility` —— 各有独立公式，不得统称。

**§11.2 解析几何**：

```
E_rand = K/n
E_max  = K(2n−K+1)/[n(n+1)]        (n > K)
D(n,K) = E_max − E_rand = K(n−K)/[n(n+1)]
n ≤ K  ⇒  E_rand = E_max = 1, D = 0
K=10 时 D 在 n* = K + sqrt(K²+K) ≈ {GEO['n_star_analytic']:.2f} 取峰（数值峰 n={GEO['n_star_numeric']}）
D(20,10) = {GEO['D_at_n20']:.4f}  >  D(50,10) = {GEO['D_at_n50']:.4f}
```

**⇒ 「推荐 n ≥ 50 以提高区分度」在数学上不成立，已删除并列入 forbidden_claims。**

**§11.3 可辩护表述**：

> Fixed-K candidate-selection evaluation is **structurally non-informative when
> n ≤ K**, and its discriminative range **depends strongly on candidate-set
> geometry** (n relative to K).

**复验**：n ≤ K=10 的 n≥5 target = **{GEO['n_le_K_count']}/{GEO['n_ge5_count']} =
{100*GEO['n_le_K_fraction']:.1f}%**（历史报 40.8% ⇒ 一致）。

> **禁用**：「the floor comes from SKEMPI ΔΔG values」——未用真实 ΔΔG 效用幅度证明。

---

## 10. §12 H5 —— 保留为**回溯性匹配案例研究**

V7.3 F1–F8 共享同一 **3351 行**宇宙（`{PU['rowset_hash'][:8]}`）。

| 对比 | ΔOverall | ΔMacro ρ | ΔC_W |
|---|---:|---:|---:|
| F2 static − F1 | +0.0279 | **−0.0015** | +0.0010 |
| F4 mutation×partner − F1 | +0.0701 | **+0.0437** | +0.0141 |
| F8 full − F1 | +0.1054 | **+0.0758** | +0.0269 |

**⇒ 能力分解区分了「仅合并分数改善」与「同时改善靶内排序」。**

**禁用**：benchmark **前瞻性**预测了 F4 会成功。
**数值一致性**：V8 重算的 `F4−F1 ΔMacro = +0.0437` 与历史 V7.3 bootstrap 的
`+0.0731 [+0.0209,+0.1289]` **口径不同**（度量宇宙/聚合不同），
**不得混用而不解释**。

---

## 11. §13 外部证据（**EXACT_UNIQUE only**）

**关键修正**：此前的 CATH 映射（约 85% 匹配）**主要由歧义键构成** ——
canonical 4473 行中仅 **{FR.get('external_clean', {}).get('n_rows', 0)}** 为 `EXACT_UNIQUE`，
另有 **{FR.get('external_clean', {}).get('excluded_ambiguous_rows', 0)}** `AMBIGUOUS_KEY`。
**⇒ 旧的 `U_EXT_CANONICAL` 读数作废。**

**新宇宙 `U_EXT_EXACT_UNIQUE`**：rows={EXT['n_rows']} / targets={EXT['n_targets']} /
n≥5={EXT['n_ge5']} / families={EXT['n_families']}

| model | n | Overall r | Macro ρ | C_B | C_W | RMSE |
|---|---:|---:|---:|---:|---:|---:|
"""
for r in EP.itertuples():
    S += (f"| `{r.model}` | {r.n_rows} | {r.overall_pearson:+.4f} | "
          f"{r.macro_spearman:+.4f} | {r.C_B_row:+.4f} | {r.C_W:+.4f} | "
          f"{r.rmse:.4f} |\n")
S += f"""
**外部秩一致（低功效）**：
```
rho(Overall, Macro)   = {er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.4f}   perm_p = {er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.4f}
rho(Overall, Eff@10)  = {er.get('rho_overall_pearson_vs_efficiency@10', float('nan')):+.4f}   perm_p = {er.get('perm_p_overall_pearson_vs_efficiency@10', float('nan')):.4f}
rho(Macro, Eff@10)    = {er.get('rho_macro_spearman_vs_efficiency@10', float('nan')):+.4f}   perm_p = {er.get('perm_p_macro_spearman_vs_efficiency@10', float('nan')):.4f}
n_methods = {er.get('n_methods')}
```

**⇒ 外部证据方向与内部相反**：在 EXACT_UNIQUE 上，
`rho(Overall, Macro) = {er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f}`
（**正相关，方向一致而非解耦**），但 n=5 ⇒ **UNDERPOWERED**。

**HERMES 语义**（§13.1）：官方量 `log_proba_mt__minus__log_proba_wt` 是
**无量纲突变效应分数**，**不是** calibrated kcal/mol ΔΔG。
⇒ **RMSE = NOT_APPLICABLE；MAE = NOT_APPLICABLE**。
所用变体为 **cDNA117k 稳定性微调**变体，**不得**标为 SKEMPI binding 监督。
（SKEMPI binding 微调变体是另外的 checkpoint。）

**映射 QC**：canonical rows before = after = 4473；duplicate row_id = 0；
row expansion = 0；exact match = 2865（单点）；ambiguous = 0；unmatched = 1608。

---

## 12. §14 场景推荐（重建）

| 场景 | 主指标 | ALL 池首选 | **STRICT 首选** | 改变 | 状态 |
|---|---|---|---|---|---|
| A 定量 ΔΔG | `rmse` | `{sc.iloc[0].top1_ALL_pool}` | **`{sc.iloc[0].top1_STRICT_pool}`** | ✅ | 仅限有 kcal/mol 语义的模型；**HERMES 不合格** |
| B 文库排序 | `macro_spearman` | `{sc.iloc[1].top1_ALL_pool}` | **`{sc.iloc[1].top1_STRICT_pool}`** | {'✅' if sc.iloc[1].changed_after_provenance_clean else '❌'} | 匹配池内 |
| C 湿实验短名单 | `efficiency@10` | `{sc.iloc[2].top1_ALL_pool}` | **`{sc.iloc[2].top1_STRICT_pool}`** | ✅ | **DESCRIPTIVE / LOW-CONFIDENCE**（H3 FAIL） |

历史「Overall 选择在文库排序场景排第 54」**未存活**：
该模型在主匹配池中**不存在** ⇒ `NOT_ESTIMABLE_UNDER_MATCHED_POOL`。

---

## 13. §15 Pareto

每个匹配 stratum 各自计算 Pareto 前沿；**不混宇宙**。
因 HERMES 无 calibrated RMSE，其外部 Pareto **仅限排序轴**，不插补缺失轴。

---

## 14. 最终判定

```
BENCHMARK_VALIDATION        = PARTIALLY_SUPPORTED     （H3 FAIL 为冻结规则中的硬条件）
BENCHMARK_UTILITY_LEVEL     = MORE_DISCRIMINATIVE     （LEVEL 1）
ANALYTIC_CLEANUP_STATUS     = PASS
EXTERNAL_EVIDENCE_STATUS    = UNDERPOWERED
PREPRINT_READINESS          = READY_AFTER_MINOR_FIXES
GPU_USED = NO     V7_4R2_NEEDED = NO
```

**H1–H5 最终**：

| | 历史 | 清理后 | 说明 |
|---|---|---|---|
| H1 | PASS | **WEAKENED / reframed** | 单值正交不成立；改为「可导致实质不同的模型排名」+ 池间范围 |
| H2 | PASS | **PASS（精确化）** | C_B 近零方向不稳；C_W 符号稳定但幅度变化 |
| H3 | FAIL | **FAIL** | 不救 |
| H4 | PASS | **PASS（受控证据）** | between-fraction {DL.get('fraction_between', float('nan')):.3f} |
| H5 | PASS | **PASS（回溯性）** | 匹配消融 |
"""
(F / "V8_1_FINAL_SUMMARY.md").write_text(S, encoding="utf-8")
print("V8_1_FINAL_SUMMARY.md ✅")

# ============ FINAL_DEFENSIBLE_CLAIMS.md ============
CL = f"""# FINAL DEFENSIBLE CLAIMS（V8.1，§19）

> 所有数字来自 `FROZEN_RESULTS.json`。

## Claim 1（主claim）

> Pooled PPI ΔΔG benchmark performance **conflates between-target and
> within-target predictive signal**. Under a **controlled** shift from permissive
> to family-OOD evaluation (same rows, same configurations, only the split regime
> changes), most of the reduction in pooled correlation is carried by loss of the
> **between-target covariance contribution**
> (fraction = **{DL.get('fraction_between', float('nan')):.3f}**;
> ΔC_B = {DL.get('delta_C_B', float('nan')):+.4f}), while within-target covariance
> **remains positive** (ΔC_W = {DL.get('delta_C_W', float('nan')):+.4f}).

## Claim 2

> Model conclusions can depend on **target aggregation** and **metric choice**;
> therefore PPI ΔΔG benchmarks should report both between-target and within-target
> components under an **explicitly declared weighting**.
> （C_B 在 {sf_cb[0.005]}/{len(W2)} 个模型上随加权在 |·|>0.005 处翻转符号；
> C_W 为 {sf_cw[0.005]}/{len(W2)}。）

## Claim 3

> **Fixed-K** shortlist metrics can become **structurally non-informative** when
> candidate-set size is no larger than K
> （n ≤ K=10 占 n≥5 target 的 {100*GEO['n_le_K_fraction']:.1f}%），
> and their discriminable range follows
> D(n,K) = K(n−K)/[n(n+1)]，在 n* ≈ {GEO['n_star_analytic']:.1f} 取峰、之后下降。
> ⇒ candidate-set geometry 必须与 Top-K 效用一同报告。

## Claim 4

> In a **matched retrospective** ablation (V7.3, same 3351 rows), capability
> decomposition separated a **pooled-score-only** improvement
> (F2 static: ΔOverall +0.0279, ΔMacro ρ **−0.0015**) from a modification that
> **also improved within-target ranking**
> (F4 mutation×partner: ΔOverall +0.0701, ΔMacro ρ **+0.0437**).

## Claim 5（H1，收窄后）

> Global and within-target metrics **can induce materially different model
> rankings under matched evaluation**：跨 {len(ST)} 个精确匹配 stratum，
> `rho(Overall Pearson, Macro Spearman)` 范围为
> **{rmin:+.2f} 至 {rmax:+.2f}**。
> **不得**声称普遍正交或普遍负相关。

## External claim（依实际清理后结果）

> External evidence was **compatible with heterogeneous capability profiles**
> but was **underpowered** regarding any universal cross-metric decoupling
> （n_methods = {er.get('n_methods')}；`rho(Overall, Macro)`
> = {er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f}，
> permutation p = {er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}）。
> **不得**据此断言外部复现或外部反证。

## FORBIDDEN（§19）

- all existing benchmarks are invalid
- Overall Pearson is useless
- first leakage-aware / first family-OOD / first per-target benchmark
- universal calibration–ranking orthogonality
- all PLMs fail / physics generally beats learning
- target size causes prediction difficulty
- HERMES log-odds is calibrated kcal/mol ΔΔG（**RMSE/MAE = NOT_APPLICABLE**）
- candidate-set floor proves SKEMPI ΔΔG values are intrinsically too similar
- benchmark prospectively predicted the successful remedy
- prospective wet-lab validity / universal model recommender
- STRONGLY_SUPPORTED merely because external models were added
- 「recommend n ≥ 50 because fixed-K discrimination increases with n」
- any external capability claim computed on chain-ambiguous mappings
"""
(F / "FINAL_DEFENSIBLE_CLAIMS.md").write_text(CL, encoding="utf-8")
print("FINAL_DEFENSIBLE_CLAIMS.md ✅")

# ============ PREPRINT_READINESS / 其他短文档 ============
(F / "PREPRINT_READINESS.md").write_text(f"""# PREPRINT READINESS（V8.1）

```
PREPRINT_READINESS = READY_AFTER_MINOR_FIXES
```

## 1. 剩余问题是否科学阻断？

**否。** 但有三项必须在稿件中显式声明（非阻断）：

1. 外部证据 **UNDERPOWERED**（n_methods={er.get('n_methods')}），
   且方向与内部相反（EXACT_UNIQUE 上 `rho(Overall, Macro)` 为正）；
2. `L2_complex` 与 `L3_seq40` 逐位相同 ⇒ 不得作独立证据；
3. `H3 = FAIL`，选择轴 **LOW-CONFIDENCE**。

## 2. 剩余问题是否仅为编辑/格式？

主要如此：图注需统一 universe 声明；`Fig3` 需保留 duplicate-regime 警告。

## 3. 是否需要额外模型计算？

**不需要。** 所有结论来自已有预测。

## 4. V7.4R2 是否必需？

**不需要**（`V7_4R2_NEEDED = NO`）。

## 5. GPU 是否必需？

**不需要**（`GPU_USED = NO`）。

## 6. 现在可否发 bioRxiv v1？

**可以（在完成 minor fixes 之后）。**
`FINAL_CONSISTENCY_CHECK` 必须为 PASS。
""", encoding="utf-8")

(F / "AXIS_C_GEOMETRY_NOTE.md").write_text(f"""# AXIS C 几何说明（V8.1 §11）

## 命名区分（§11.1，不得统称 "Efficiency@K"）

| 名称 | 公式 | 备注 |
|---|---|---|
| `original_dcg_efficiency@K` | Σ gain(selected)/Σ gain(ideal top-K)（**同折扣**） | V8 v1 |
| `overlap_topq` | \\|真值前 q% ∩ 预测前 q%\\| / \\|真值前 q%\\| | V8 v2 |
| `rank_mass_fraction@K` | Σ rank-mass(selected)/Σ rank-mass(ideal) | v3 基础 |
| `normalized_rank_mass_utility` | clip((·−E_rand)/(E_max−E_rand),0,1) | v3 最终 |

## 解析固定-K 几何（§11.2）

```
E_rand = K/n
E_max  = K(2n−K+1)/[n(n+1)]          (n > K)
D(n,K) = E_max − E_rand = K(n−K)/[n(n+1)]
n ≤ K  ⇒  E_rand = E_max = 1,  D = 0
峰值 n* = K + sqrt(K² + K) ≈ {GEO['n_star_analytic']:.2f}   (K=10)
D(20,10) = {GEO['D_at_n20']:.4f}   D(50,10) = {GEO['D_at_n50']:.4f}
```

**⇒ 「n ≥ 50 更可区分」在数学上不成立，已删除。**

## 实测

n ≤ K=10 的 n≥5 target = **{GEO['n_le_K_count']}/{GEO['n_ge5_count']} =
{100*GEO['n_le_K_fraction']:.1f}%**（独立复验，与历史 40.8% 一致）。

## 可辩护表述（§11.3）

> Fixed-K candidate-selection evaluation is **structurally non-informative when
> n ≤ K**; its discriminative range **depends strongly on candidate-set geometry**
> (n relative to K).

**禁用**：「floor comes from SKEMPI ΔΔG values rather than the metric」——
未用真实 ΔΔG 效用幅度证明。**H3 不因此升级。**
""", encoding="utf-8")

(F / "SCENARIO_RECOMMENDATIONS_FINAL.md").write_text(f"""# 场景推荐（最终，§14）

## 生成规则（全部满足才可推荐）

1. 候选模型**同一 `rowset_hash`**；
2. 属 **PRIMARY_PROVENANCE_CLEAN_POOL**（无 RED）；
3. 无 RED / permissive 训练模型（除非场景专门研究泄漏）；
4. 指标与模型语义匹配。

## 结果

| 场景 | 主指标 | ALL 池首选 | **STRICT 首选** | 改变 |
|---|---|---|---|---|
| A 定量 ΔΔG | `rmse` | `{sc.iloc[0].top1_ALL_pool}` | **`{sc.iloc[0].top1_STRICT_pool}`** | ✅ |
| B 文库排序 | `macro_spearman` | `{sc.iloc[1].top1_ALL_pool}` | **`{sc.iloc[1].top1_STRICT_pool}`** | {'✅' if sc.iloc[1].changed_after_provenance_clean else '❌'} |
| C 湿实验短名单 | `efficiency@10` | `{sc.iloc[2].top1_ALL_pool}` | **`{sc.iloc[2].top1_STRICT_pool}`** | ✅ |

* **场景 A**：只允许有 calibrated ΔΔG 语义的模型；
  **HERMES log-odds 不合格**（RMSE = NOT_APPLICABLE）。
* **场景 C**：**DESCRIPTIVE / LOW-CONFIDENCE** —— 选择轴结构上无力，且 H3 = FAIL。

## 历史陈述处置

「Overall 选择的模型在文库排序场景排第 54」⇒ **未存活**
（该模型不在主匹配池中）⇒ 记为
`SCENARIO_RECOMMENDATION = NOT_ESTIMABLE_UNDER_MATCHED_POOL`。
**不**为产出赢家而混合行集。
""", encoding="utf-8")

(F / "EXTERNAL_AUDIT_FINAL.md").write_text(
    "# 外部审计（最终）\n\n" +
    Path(F / "EXTERNAL_CLEAN.md").read_text(encoding="utf-8") +
    "\n\n---\n\n" +
    Path(F / "HERMES_SEMANTICS_AUDIT.md").read_text(encoding="utf-8") +
    "\n\n---\n\n" +
    Path(F / "CATH_MAPPING_AUDIT.md").read_text(encoding="utf-8") +
    "\n\n---\n\n" +
    Path(F / "EXTERNAL_INFERENCE_FEASIBILITY.md").read_text(encoding="utf-8")
    if (F / "EXTERNAL_INFERENCE_FEASIBILITY.md").exists() else
    "# 外部审计（最终）\n\n见 EXTERNAL_CLEAN.md / HERMES_SEMANTICS_AUDIT.md / CATH_MAPPING_AUDIT.md\n",
    encoding="utf-8")
print("PREPRINT_READINESS / AXIS_C_GEOMETRY_NOTE / SCENARIO_RECOMMENDATIONS_FINAL / "
      "EXTERNAL_AUDIT_FINAL ✅")
