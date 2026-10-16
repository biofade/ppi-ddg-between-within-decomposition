# V7.1.1 Final Summary

**时间**：2026-10-14（约 18:0x → 18:3x，约 30 分钟；上限 60 分钟，未超时）
**约束**：未修改 V1–V7.1 任何文件；新命名空间 `v7_1_1`；CPU only；**未下载**；
**未训练任何新模型**；未开始 residue-level embedding；未开始 V7.2

**完成度**：PHASE 0/1/2/3/4/5/6/7 **全部 COMPLETED**

---

## 1. Provenance

- **hashes**：17 个输入冻结于 `phase0_provenance/input_manifest.csv`
  （V7.1 的 u_local_rows / u_full_rows / ranker_metrics / z_predictions /
  location_predictions / lsr_predictions / oracle_models / contrasts /
  stratification / final_summary；V7.0 的 matched_rows / block_L；
  canonical matched set；ddg_split family）
- **upstream unchanged**：**YES**（本任务只读 V7.0/V7.1 产物，未写入其目录）

---

## 2. 208 vs 214 Reconciliation

| universe | n_rows | n_targets |
|---|---:|---:|
| Canonical targets（`matched_analysis_set.pdb_id`） | 7523 | **230** |
| U_FULL | 7523 | 230 |
| U_LOCAL（`block_L*.pdb_id`，源自 `mutation_sites`） | 4473 | **208** |
| **ranker/LSR universe（`d['groups']`）** | 4473 | **214** |
| metric universe（Macro ρ，n≥5） | 4288 | **130** |

### Cause: **`TARGET_ID_NAMESPACE_MISMATCH`**

**三个已确证事实**（全部来自 row ID / target ID 的精确对账，无猜测）：

1. `d['groups']`（`build_pair_features`）vs canonical `pdb_id`：
   一致率 **100.00%**（7523/7523）→ 二者 unique 均为 **230**
2. `block_L_local.csv` / `block_L_rescued.csv` 的 `pdb_id` vs canonical `pdb_id`：
   一致率仅 **94.21%**（4214/4473），**259 行不一致**；
   且含 canonical 中**不存在**的 `1Y1K`（该行 pair_row=3791 在 canonical 中是 `1YCS`）
   → 该 `pdb_id` 源自 `mutation_sites.csv`，**不是** canonical 命名
3. 因此**同一 4473 行集合**按两种 `pdb_id` 聚合 → **208 vs 214**

**排除项**：不是 feature intersection、不是 metric eligibility、不是纯报告文本错误
—— 而是**两个真实不同的 target 命名空间**。

### Exact targets（按 row 集合差）

- `T_214 − T_208`：**6 个**（含 `1Y1K` 等，见 `target_membership_matrix.csv`）
- `T_208 − T_214`：**0 个**（208 是 214 在旧命名下的"合并"结果：259 行被归到相邻 PDB 名）

### 影响量化（`target_id_impact.csv`）

| 模型 · 指标 | (A) `d['groups']` 标签 | (B) 旧 `block_L.pdb_id` | Δ |
|---|---:|---:|---:|
| P0_raw · Macro ρ | **+0.1117** | +0.0945 | +0.0172 |
| P0_raw · target-mean r | **+0.1908** | +0.1224 | +0.0684 |
| P0_raw · C_W | +0.0441 | +0.0474 | −0.0033 |
| P0_raw · n_target | 214 | 208 | +6 |

### 受影响 / 未受影响

| 环节 | 受影响? | 说明 |
|---|---|---|
| LSR 的 target 标签（PHASE 3/4/5/6 聚合） | **未受影响** | 自始用 `d['groups']`（100% 对齐） |
| S/G 特征取回（`reindex(pdb_rows)`） | **未受影响** | S/G 索引同源于 canonical |
| BLOCK L 特征与 target 的对应 | **受影响（259 行）** | 该 259 行的 local 特征归属错 target |
| LSR / Oracle 的数值结论 | **未受影响** | 量级见上表（Δmacro ρ = +0.0172） |

---

## 3. Metric Eligibility

> **核心区分**：`DATA COVERAGE`（数据可用性）vs `METRIC ELIGIBILITY`（指标的最小 n 门槛）。
> **禁止**把后者写成模型 coverage loss。

数据：LSR predictions，**4473 行 / 214 targets**

| 阈值 | 满足的 target 数 |
|---|---:|
| n ≥ 1 | **214** |
| n ≥ 2 | **180** |
| n ≥ 3 | **164** |
| n ≥ 5 | **130** |
| n ≥ 10 | **85** |

| metric | min n | eligible targets | rows |
|---|---:|---:|---:|
| Overall Pearson | ≥1 | **214** | 4473 |
| RMSE | ≥1 | **214** | 4473 |
| `C_B` / `C_W` | ≥1 | **214** | 4473 |
| target-mean Pearson | ≥1 | **214** | 4473 |
| attenuation slope | ≥3 | **164** | 4407 |
| **Macro per-target Spearman** | **≥5** | **130** | 4288 |
| Kendall tau-b | ≥2 | **180** | 4439 |
| Pairwise Accuracy | ≥2 | **180** | 4439 |
| **Efficiency@10** | **≥5** | **130** | 4288 |
| Hit@10 / Precision@10 | ≥2 | **180** | 4439 |

### V7.1 数字的 target 基数（已复算确认）

| V7.1 报告值 | 复算值 | 使用的 target 数 |
|---|---:|---:|
| Macro ρ = 0.1117 | **+0.1117** ✅ | **130**（n≥5） |
| Eff@10 = 0.2522 | **+0.2522** ✅ | **130**（n≥5） |
| C_W = 0.0441 | **+0.0441** ✅ | **214**（全部） |

⇒ 同一个预测，`Macro ρ` 基于 **130** 个 target、`C_W` 基于 **214** 个 target。
**这些差异完全是 metric eligibility，不是模型 coverage loss。**

---

## 4. Outer / Inner Protocol

**Outer**
- unit: **family**
- folds: **5**（`GroupKFold` on family），seeds **42/43/44**
- 断言：`set(fam[train]) ∩ set(fam[test]) == ∅`
- 名称冻结为 **`OUTER_FOLDS_V7_2`**
- 历史 L4 可复现性：**YES**（`GroupKFold` 确定性，fold 只依赖 seed）

**Inner**
- unit: **family**
- folds: `k = min(3, max(2, n_unique_family_in_outer_train))`
- 断言：同 outer（已实测 family 零重叠 ✅）
- 规则**统一**，不随 fold 变化

**Selection metric**：`MacroSpearman`（每 target n≥5）

**Tie-break**：`Eff@10` → `C_W` → **simpler model**（`eps = 0.005`）

**`OUTER_LEAKAGE_CHECK`：PASS**
（协议要求所有选择只在 inner family CV 内进行；V7.1 的 outer-based 选择已作废）

---

## 5. Ranker Protocol

**Standardized z**（**唯一函数** `standardize_within_target()`）：

```
z_ci = (s_ci − mean_c(s)) / (SD_c(s) + eps)
```

**Collapsed target handling**：`SD_c(s) < 1e-6 → z = 0`，`PRED_SCORE_COLLAPSED = 1`
（各现有 ranker 的 collapse rate 实测见 `ranker_selection_protocol.md`）

**Final nested procedure**：

```
for each outer fold f:
    for each candidate ranker r (architecture family 事先固定):
        inner_cv_score(r) = mean over inner folds of MacroSpearman
    r* = argmax inner_cv_score   (tie-break: Eff@10 → C_W → simpler)
    retrain r* on full outer-train
    predict outer-test → OOF
concatenate OOF across folds
```

- **不同 outer fold 可选择不同超参/变体**；**architecture family 必须事先固定**
- **No global winner peeking**：V7.1 的 `R2_MLP_L0` 选择**正是此类泄漏**，
  故在 V7.2 中**作废**，仅作**描述性基线**保留

**允许比较的 ranker 家族（V7.1.1 不实现）**：
`R_BASE` / `R_RESIDUE_LOCAL` / `R_RESIDUE_LOCAL_PARTNER` /
`R_RESIDUE_PERTURBATION` / `R_RESIDUE_PERTURBATION_SET`

---

## 6. Calibration Shell

**Primary**：

```
y_hat_ci = mu_hat_c + alpha * z_ci
```

**mu model**：**MU3（target + set）**，且**必须在 outer-train 内重新拟合**

**alpha policy**（仅比较两个简单方案，按 **inner-family RMSE** 选）：

| alpha 方案 | inner-family RMSE |
|---|---:|
| `alpha = 1` | **2.5151** |
| `alpha = train 残差尺度` | 2.5078 |

- 差异 = **0.0072**（< 0.01 门限）→ 按「**差异极小则优先简单**」规则
- **冻结：`alpha = 1`**

**Secondary**：full sigma diagnostic only（`mu + sigma_set * z`），
**不作为 V7.2 主架构**

**永久 assert（已实测 PASS ✅）**：

```python
assert within_target_order(raw_z) == within_target_order(mu + alpha*z)
```

---

## 7. V7.2 Evaluation Contract

见 `experiments/v7_1_1/V7_2_EVALUATION_CONTRACT.md`（**FROZEN**）

- **Primary**：`Macro Spearman`（n≥5）
- **Secondary**：`Eff@10`、`C_W`、`Kendall tau-b`、`Pairwise acc (δ=0.5)`
- **Calibration**：`Pearson`、`RMSE`、`target mean`、`slope`
- **Bootstrap**：**5000** target clusters（paired）
- **Weighting**：**必须**同时报 row-weighted 与 equal-target `C_B`
- **Subsets**：all eligible / 2–9 / 10–29 / ≥30
- **Success**：`ΔMacroSpearman > 0` 且 95% CI 排除 0，
  且（`ΔEff@10 > 0` 或 `ΔC_W > 0`），且无灾难性 RMSE 退化

---

## 8. Corrections to V7.1

1. **`4473 rows / 208 targets` → `4473 rows / 214 targets`**
   （208 是旧 `block_L.pdb_id` 命名空间下的计数；权威命名空间是 `canonical.pdb_id`）
2. **target 覆盖率 `90.4%` → `93.0%`**（214/230）
3. **V7.1 的 ranker 选择（`R2_MLP_L0`，依据 outer L4 macro ρ）为泄漏**
   → 在 V7.2 中**作废**，降级为描述性基线；V7.2 必须用嵌套 inner-family 选择
4. **`OUTER_FOLDS_V7_2` 与 `protocol.json` 为新增冻结产物**（V7.1 无此协议）
5. V7.1 PHASE 1 的「救援未增加任何行」**结论成立**，但其 `208` 计数应为 `214`
6. **V7.1 §7 的 Oracle 解读须软化**（见 §9 下方"重解释规则"）：
   不得写「Oracle analysis proves ranking is intrinsically more important
   than calibration」，因为 `mu/sigma` 仿射标定**在构造上不可能**改变 within 排序，
   而 oracle z **平凡地**让排序变完美。

---

## 9. V7.2 Readiness

**`READY_FOR_V7_2`: YES**

**Blocking issues**：**无**

**但以下为已登记的先决条件（非阻塞，需在 V7.2 首个 PHASE 解决）**：

1. **residue-level 表示尚不存在** —— `features/v2_esm2_cache.npz` 只有
   **pooled 链级 640 维**，无 residue-level 嵌入；且 §13 禁止重训大型 PLM。
   ⇒ V7.2 第一步必须是**取得 residue-level 表示**（外部服务或小模型），
   否则 `R_RESIDUE_*` 家族无法实现。
2. **target 标签必须统一使用** `canonical_target_labels.csv`
   （不得再用 `block_L*.csv` 的 `pdb_id`）
3. **`alpha = 1` 已冻结**；`sigma` 方向不再投入（oracle 增益 ≈ 0）
4. V7.2 的 ranker 对比必须走**嵌套 family CV**，不得复用 V7.1 的 outer 选择

### ⚠️ 重解释规则（V7.1.1 §7 强制）

**不得**写：

```
Oracle analysis proves ranking is intrinsically more important than calibration.
```

因为：`mu/sigma` 仿射标定**在构造上不可能**改变靶点内排序，
而 oracle z **平凡地**使排序完美。**正确表述**：

> Existing deployable calibration has already achieved statistically supported
> gains in absolute prediction, while current family-OOD within-target ranking
> remains low; therefore ranking representation is the highest-priority
> engineering target.

---

## 10. Exact Next Step

**If READY（本轮判定为 YES）**：

```
Residue-level ranking model development.
```

具体顺序：

1. **取得 residue-level 表示**（不得重训大型 PLM；可用外部服务或小型结构感知模型）
2. 实现 `R_RESIDUE_PERTURBATION`（`Δh_m = h_mut_m − h_wt_m`）
   与 `R_RESIDUE_PERTURBATION_SET`（`Δh_m · c_m` 及 candidate-set 条件化）
3. 在 `OUTER_FOLDS_V7_2` + inner family CV 下做**嵌套 ranker 选择**
4. 用 `V7_2_EVALUATION_CONTRACT.md` 的 primary/secondary 指标与 5000 次
   paired cluster bootstrap 评估
5. calibration shell **固定为 `mu(MU3, outer-train 内重拟合) + 1.0 * z`**

---

## 11. 文件清单

### 脚本（`scripts/v7_1_1/`）
`01_universe_audit.py`、`02_reconciliation_final.py`、`03_protocol.py`

### 结果（`experiments/v7_1_1/`）
- `phase0_provenance/`：`input_manifest.csv`、`provenance.md`
- `phase1_universe/`：`universe_summary.csv`、`target_membership_matrix.csv`、
  `universe_reconciliation.md`、`reconciliation_final.md`、
  `target_id_impact.csv`、**`canonical_target_labels.csv`（权威标签）**
- `phase2_metrics/`：`metric_eligibility.csv`、`metric_subset_report.md`
- `phase3_nested_cv/`：`nested_cv_protocol.md`、`fold_manifest.csv`、`protocol.json`
- `phase4_ranker_protocol/`：`ranker_selection_protocol.md`
- `phase5_calibration/`：`calibration_shell.md`、`simple_scale_comparison.csv`、
  `calibration_choice.json`
- `V7_2_EVALUATION_CONTRACT.md`（**FROZEN**）
- `V7_1_1_FINAL_SUMMARY.md`（本文件）

### 日志（`logs/v7_1_1/`）
无 `failures.log` → **无功能性失败**
