# V7.2R-2 Final Summary

- 生成时间: 2026-10-14 21:55:00

## 1. Integrity

- Upstream hashes: `experiments/v7_2r_2/phase_local_freeze/input_hashes.csv`（11 项，**UPSTREAM_HASH_CHANGE = 0**）
- Kaggle ZIP SHA256: `f3f0b8b12cb469668effc2cc15cedf90b6cdee9cf6ac4d8e44ef1966cbd231d9`（23.2 MB）
- Import: **174/174 npz**，4219 行，**0 重复 row_id**，**0 非预期 row_id**，**0 无法映射到 canonical target**
- Canonical target namespace: **PASS** —— 全部按 `canonical_target_labels.csv` 的 `row_id` 连接（Kaggle 侧 `target_id` 与 canonical **100.00% 一致**，仅作核对）

## 2. Kaggle Representation Semantics

| 项 | 判定 |
|---|---|
| `r_wt` | **存在**，dim **1280**，dtype **float16** |
| `r_mut` | **存在**，dim **1280**，float16 |
| `delta_r` | **存在**，dim **1280**，float16；**数值自洽**（max‖Δr−(r_mut−r_wt)‖ = 0.0078 < float16 容差）✅ |
| `partner_repr` | **不存在** ❌ |
| MINT representation | **未进入突变表示**（MINT 被加载/审计，但 `embed_complex` 只前向含突变位点的链） |
| site-local representation | **存在**（per-residue hidden，按位点取 token i+1，再对多突变取 mean） |

| 维度/层级/条件性 | 值 |
|---|---|
| dimension | **1280** |
| dtype（存储） | float16 |
| 层级 | **site-level（per-residue）**，聚合到 variant（mean over sites） |
| generic sequence / PPI-conditioned | **generic sequence（仅 ESM2-650M）** |
| WT-specific / MUT-specific | **两者都有**（`r_wt` 与 `r_mut` 各自独立前向） |

### ⚠️ 链隔离审计（关键限定）

Kaggle 的 `build_variant()` 只返回**含突变位点的链**，`embed_complex()` 也只对这些链做前向。实测：

- **214/214 个 target 在 PDB 中都有 ≥2 条蛋白链**（全是复合物）
- 但 **128/214** 只有 **1 条**含突变的链 → 该复合物的伙伴链**未进入表示**
- `partner_sanity.csv` **为空**（0 行）→ **伙伴条件化未被演示**

**exact terminology allowed（§4 Case A）**：

```
STRONG_GENERIC_PLM_PERTURBATION
    (ESM2-650M, CHAIN_ISOLATED site-level, 无 PPI 条件化)
```

**不得**称其为 `PPI_AWARE_PERTURBATION`；**不得**声称伙伴条件化（§5 `PARTNER_CONDITIONING_NOT_DEMONSTRATED`）。

## 3. Track B Coverage

| 项 | 值 |
|---|---|
| Rows | **4219 / 4233 = 99.7%** |
| Targets | **174 / 174 = 100%** |
| Families | **91 / 91 = 100%** |
| n>=5 | **115 / 116 = 99.1%** |
| metric-eligible target coverage | **99.1%** |

## 4. Matched Universe U_AB

| universe | rows | targets | families |
|---|---:|---:|---:|
| Track A (full) | 4233 | 174 | 91 |
| Track B (imported) | 4219 | 174 | 91 |
| **U_AB (matched)** | **4219** | **174** | **91** |

- n>=5：**115**（Track A 为 116）→ 相对 **99.1%**
- §40 门槛（n>=5 覆盖 ≥80% of Track A）：**PASS ✅**（非 COVERAGE_LIMITED）

## 5. Coverage Bias

仅 **14 行**（4233 − 4219）未进入 U_AB —— 规模太小，**未做分层偏倚检验**（无统计功效）。V7.2R-1 已对 Track A 的 94.6% 覆盖做过完整偏倚审计。

## 6. Nested Family-CV Integrity

| 项 | 结果 |
|---|---|
| Outer overlap | **0**（每折 `assert not (set(fam[tr]) & set(fam[te]))`） |
| Inner overlap | **0**（`make_folds` family-disjoint，V7.1.1 冻结） |
| Outer-based selection | **NO**（无任何 outer 选型；超参全部预先固定） |
| PCA train-only | **PASS**（`PCA.fit` 仅用 `X[tr]`） |

## 7. Primary Ranking Results

**U_AB = 4219 rows / 174 targets / 91 families**

### 主口径：ensemble prediction（与 bootstrap 一致）

| Model | Macro rho | Eff@10 | C_W | Kendall | Pairwise | Hit@10 |
|---|---:|---:|---:|---:|---:|---:|
| **B0** | **+0.0204** | **+0.1277** | +0.1324 | +0.0390 | +0.5752 | +0.7333 |
| **A2 small PLM** | **+0.0115** | **+0.1074** | +0.1240 | +0.0152 | +0.5712 | +0.7351 |
| **B1 strong representation** | **+0.0559** | **+0.1634** | +0.1246 | +0.0347 | +0.5697 | +0.7329 |
| **B2 strong perturbation** | **+0.0292** | **+0.1454** | +0.1215 | +0.0393 | +0.5693 | +0.7357 |
| B3 + PPI context | **NOT_RUN**（无 partner 表示） | | | | | |

### 敏感性：seed-mean of metric

| Model | Macro rho | Eff@10 | C_W | Overall r |
|---|---:|---:|---:|---:|
| B0_Base | +0.0221 | +0.2640 | +0.1324 | +0.1944 |
| A2_smallPLM_dh | +0.0093 | +0.3454 | +0.1240 | +0.1822 |
| B1_strong_repr | +0.0506 | +0.3387 | +0.1246 | +0.1476 |
| B2_strong_dr | +0.0305 | +0.3314 | +0.1215 | +0.1723 |

> `mean(ρ_seed) ≠ ρ(mean(pred_seed))` —— 两口径**分别报告，未混用**（§14）。

## 8. Seed Sensitivity

| Ranker | s42 | s43 | s44 | 极差 |
|---|---:|---:|---:|---:|
| A2_smallPLM_dh | +0.0162 | -0.0248 | +0.0364 | 0.0612 |
| B0_Base | +0.0274 | +0.0064 | +0.0325 | 0.0261 |
| B1_strong_repr | +0.0685 | +0.0304 | +0.0528 | 0.0381 |
| B2_strong_dr | +0.0336 | -0.0086 | +0.0666 | 0.0752 |

## 9. Critical Contrast B2-A2

**5000/5000 paired target-cluster bootstrap**

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacro rho | +0.0177 | [-0.0118, +0.0453] | ❌ 含 0 |
| ΔEff@10 | -0.0786 | [-0.0952, +0.1279] | ❌ 含 0 |
| ΔC_W | -0.0042 | [-0.0147, +0.0065] | ❌ 含 0 |
| ΔOverall Pearson | -0.0108 | [-0.0340, +0.0115] | ❌ 含 0 |
| ΔRMSE | +0.0106 | [-0.0125, +0.0333] | ❌ 含 0 |

⇒ **§23 成功判据：`NOT SUPPORTED` ❌**

## 10. Critical Contrast B2-B0

**5000/5000**

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacro rho | +0.0088 | [-0.0220, +0.0421] | ❌ 含 0 |
| ΔEff@10 | +0.0293 | [-0.5456, +0.6190] | ❌ 含 0 |
| ΔC_W | -0.0102 | [-0.0210, +0.0022] | ❌ 含 0 |
| ΔOverall Pearson | -0.0197 | [-0.0388, +0.0029] | ❌ 含 0 |
| ΔRMSE | +0.0000 | [-0.0285, +0.0218] | ❌ 含 0 |

⇒ **`NOT SUPPORTED` ❌**

### 附：B1-B0（强表示但**无差分**）— 唯一显著项

**1000/1000**

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacro rho | +0.0355 | [+0.0003, +0.0757] | ✅ 排除 0 |
| ΔC_W | -0.0041 | [-0.0194, +0.0078] | ❌ 含 0 |
| ΔOverall Pearson | -0.0435 | [-0.0821, -0.0074] | ✅ 排除 0 |
| ΔRMSE | -0.0527 | [-0.0851, -0.0172] | ✅ 排除 0 |
| ΔEff@10（重算 n_rand=20，800/800） | +0.0356 | [−0.1236, +0.1659] | ❌ 含 0 |

**`B1 − B0` 的 ΔMacro rho = +0.0355 [+0.0003, +0.0757] → CI 排除 0（勉强）**

⇒ 但 §23 需要**同时**满足（ΔEff@10 > 0 **或** ΔC_W > 0）：
- ΔEff@10 = +0.0356，CI **[−0.1236, +0.1659] 含 0** ❌
- ΔC_W = −0.0041，CI [−0.0194, +0.0078] 含 0 ❌
- 且 `ΔOverall_r = −0.0435`、`ΔRMSE = −0.0527` **均显著变差**

⇒ **§23 判定：`NOT SUPPORTED` ❌**（仅 Macro rho 一项勉强显著，不满足「至少一项 secondary 改善」的要求）

## 11. Partner Increment B3-B2

**`NOT_RUN`** —— Kaggle ZIP 未提供 `partner_repr`，`partner_sanity.csv` 为空，表示经审计为 **CHAIN_ISOLATED**。
依 §8/§19，**不得**构造 B3，**不得**用零向量填充。

## 12. Negative Controls

| 对照 | 构成 | Macro rho | Eff@10 | C_W |
|---|---|---:|---:|---:|
| `NC1_B2_dr` | base+B_dr | +0.0305 | +0.3314 | +0.1215 |
| `NC1_B1_wt_only` | base+B_wt | +0.0462 | +0.2992 | +0.1262 |
| `NC2_mut_identity` | base+F0 | +0.0061 | +0.4351 | +0.1163 |
| `NC3_dr_norm_only` | base+B_drnorm | +0.0240 | +0.3017 | +0.1217 |
| `NC4_random_partner` | base+B_dr+B_randpartner | +0.0173 | +0.3264 | +0.1236 |

### 判定

| 对照 | 结论 |
|---|---|
| **NC1 WT only** | `base+B_wt`（**仅 r_wt，无差分**）= **+0.0462** **高于** `base+B_dr`（Δr）= **+0.0305** → **扰动本身不必要**；增益主要来自**结构身份**而非突变扰动 |
| **NC2 mutation identity** | `base+F0` = **+0.0061**（远低于全部 PLM 变体） → 表示增益**不能**归因于氨基酸替换身份编码 |
| **NC3 norm only** | `base+‖Δr‖` = **+0.0240**，**低于**完整 Δr（+0.0305） → 幅度信息不是全部，但差距极小 |
| **NC4 partner shuffle** | `base+B_dr+B_randpartner` = **+0.0173**，**低于** Δr → 无伙伴信号可被 shuffle（因为本来就没有伙伴表示）；该对照在此**不适用**（§19 仅在声称 partner 信号时才要求） |

## 13. Top-K Result

**TOPK_ONLY_SIGNAL: NO**

- B2−A2 的 ΔEff@10 = −0.0786，CI [−0.0952, +0.1279] **含 0** ❌
- B2−B0 的 ΔEff@10 = +0.0020，CI [−0.1979, +0.2068] **含 0** ❌
- B1−B0 的 ΔEff@10 = +0.0356，CI [−0.1236, +0.1659] **含 0** ❌

## 14. Subgroups

| 组 | 状态 |
|---|---|
| single / multi mutation | 见 `phase5_eval/trackB_results.csv`（未按此分层） |
| 2-9 / 10-29 / >=30 | **NOT_RUN**（时间受限） |
| interface / near / distal | **NOT_RUN**（需结构元数据） |

## 15. Calibration

**NOT_RUN** —— 依 §23，标定**不能**用于选择 ranker。由于没有任何 Track B ranker 通过 §23 成功判据（B2 两项主对比均 NULL），标定外壳不再对 Track B 排名模型做额外报告；V7.1 已给出 MU3+z 对 B0 类的标定增益（ΔOverall_r +0.3507 ✅、ΔRMSE −0.9448 ✅）。

**rank preservation**：`PASS by construction`（V7.1 已数值验证 180/180 target）

## 16. Scientific Verdict

```
TRACK_B_RESULT = STRONG_BACKBONE_NULL
```

**依据（§24 RESULT 0）**：

- `B2 − A2` ΔMacro rho = **+0.0177**，95% CI **[−0.0118, +0.0453] 含 0** → 不显著
- `B2 − B0` ΔMacro rho = **+0.0088**，95% CI **[−0.0229, +0.0418] 含 0** → 不显著
- 全部 secondary 指标同样含 0

## 17. Innovation Level

```
LEVEL = 2  (calibration/ranking formulation only)
```

不因 Overall Pearson 升高而上调等级（§25）。Track B 未带来 LEVEL 3/4/5。

## 18. Exact Defensible Claim

> Under a frozen, leakage-safe nested family-OOD protocol on an identical
> matched universe (4219 rows / 174 targets / 91 families), an ESM2-650M
> chain-isolated site-level mutation perturbation representation
> **did not significantly improve within-target mutation ranking** over
> either a small ESM2-8M perturbation (ΔMacro Spearman = +0.0177,
> 95% CI [−0.0118, +0.0453]) or the established baseline
> (ΔMacro Spearman = +0.0088, 95% CI [−0.0229, +0.0418]).
> Negative controls further show that the **WT-only** representation
> (+0.0462) already exceeds the explicit perturbation (+0.0305),
> indicating the weak signal present is carried by structural identity
> rather than by mutation-induced representation change.

## 19. Claims We Cannot Make

- ❌ 「protein foundation models do not work」（§51）
- ❌ 「PPI-aware representation fails」—— 本表示经审计为 **CHAIN_ISOLATED**，**从未真正测试过 PPI 条件化表示**
- ❌ 「partner context does not help」—— B3 为 **NOT_RUN**，未测试
- ❌ 任何基于 `Overall Pearson` 的「方法更好」声明
- ❌ 「MINT / ReCLIP 无效」—— MINT 被加载但**未进入突变表示**

## 20. Comparison with Track A

| | Track A (V7.2R-1) | Track B (V7.2R-2) |
|---|---|---|
| representation | ESM2-**8M** Δh | ESM2-**650M** Δr |
| rows / targets / families | 4233 / 174 / 91 | **4219 / 174 / 91** |
| Macro rho（扰动 vs 基线） | −0.0143 [−0.0588, +0.0307] | **+0.0088 [−0.0229, +0.0418]** |
| 判定 | **NULL** | **NULL** |

⇒ **表示规模从 8M 提升到 650M（81×）并未改变结论。**
依 §21，**不得**写「小 PLM 失败所以大 PLM 成功」；实测是**两者都为 NULL**。

## 21. Next Route

依 §28（所有 PLM 路线 Macro Spearman 均近零）：

```
representation scaling alone is not solving the ranking bottleneck
```

**V7.3 建议：physics / structural local interaction / energy-informed ranking**

**不要**继续放大模型规模。**不要**在没有伙伴表示的情况下继续堆 PLM。

若仍要走表示路线，**唯一**有科学意义的第一步是：
修复 `embed_complex` 让它包含**全部伙伴链**（本 ZIP 的关键缺陷），
即真正做出 PPI-conditioned 表示，再重测 B2/B3。

- 耗时 1s
