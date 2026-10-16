# V7.2 EVALUATION CONTRACT

> **状态：FROZEN**（V7.1.1 冻结；V7.2 期间不得修改）
> **冻结时间**：2026-10-14
> **协议版本**：`V7_2_NESTED_FAMILY_OOD_v1`
> **权威 target 标签**：`experiments/v7_1_1/phase1_universe/canonical_target_labels.csv`

---

## 6.1 Primary scientific task

```
family-OOD mutation ranking
```

即：在 **未见蛋白家族** 上，对一个 target 的候选突变体做**靶点内排序**。

---

## 6.2 Primary metric

```
Macro per-target Spearman        （每 target n ≥ 5）
```

**禁止**用 `Overall Pearson` 作为 primary 或用于模型选择。

---

## 6.3 Key secondary metrics

```
Efficiency@10
C_W
Kendall tau-b                    （每 target n ≥ 2）
Pairwise accuracy (delta = 0.5)
```

---

## 6.4 Calibration secondary metrics

```
Overall Pearson
RMSE
target-mean Pearson
attenuation slope                （每 target n ≥ 3）
```

---

## 6.5 Weighting

**必须同时报告**：

```
row-weighted C_B
equal-target C_B
```

可选（**仅诊断**）：

```
sqrt(n)
capped30
log(1+n)
```

> `C_B` **不是模型固有能力**（V6.3/V7.0/V7.1 三次证明其符号可随加权口径翻转）。
> **禁止**用某一种 weighting 单独下结论。

---

## 6.6 Required target subsets

```
all eligible
2–9
10–29
>=30
```

---

## 6.7 Required uncertainty

Primary model contrasts：

```
paired cluster bootstrap
5000 iterations
cluster = target / complex
```

---

## 6.8 Primary success criterion

Residue-level ranker is **supported** only if：

```
Delta MacroSpearman > 0
AND
95% CI excludes 0
```

plus **at least one**：

```
Delta Efficiency@10 > 0
or
Delta C_W > 0
```

且 **no catastrophic RMSE degradation**（定义：ΔRMSE 的 CI 上界 < +0.5）。

---

## 6.9 Strong signal（**internal only**）

```
Delta MacroSpearman >= +0.05
```
且 bootstrap 支持正向增益。

---

## 6.10 Very strong signal（**internal only**）

```
Delta MacroSpearman >= +0.10
```
且：
```
relative Eff@10 gain >= 20%
```

---

## 6.11 Breakthrough candidate（**internal only**）

```
MacroSpearman_L4 >= 0.20
```
且相对基线的增益显著、Top-K 有改善。

> ⚠️ **该阈值仅为内部路线判断用，不得写成领域公认阈值。**

---

## 6.12 Metric eligibility 与 data coverage 的区分（**强制**）

**禁止**把 metric eligibility 写成模型 coverage loss。本项目的固定资格表：

| 指标 | 最小 target n | 合格 target 数（U_LOCAL） |
|---|---:|---:|
| Overall Pearson / RMSE / `C_B` / `C_W` / target-mean Pearson | 1 | **214** |
| attenuation slope | 3 | **164** |
| Macro per-target Spearman | 5 | **130** |
| Efficiency@10 | 5 | **130** |
| Kendall tau-b / Pairwise Accuracy / Hit@10 / Precision@10 | 2 | **180** |

---

## 6.13 冻结的嵌套协议（引用）

| 项 | 值 |
|---|---|
| Outer unit | **family**，5 folds，seeds 42/43/44 |
| Outer 断言 | `set(fam[train]) ∩ set(fam[test]) == ∅` |
| Inner unit | **family**，`k = min(3, max(2, n_family_in_outer_train))` |
| Selection | `MacroSpearman`；tie-break `Eff@10 → C_W → simpler`，`eps=0.005` |
| z 转换 | **唯一函数** `standardize_within_target()` |
| Collapsed | `SD < 1e-6 → z = 0`，`PRED_SCORE_COLLAPSED = 1` |
| 超参 | **仅**在 inner CV 预注册 |
| Outer peeking | **禁止**（V7.1 的 `R2_MLP_L0` 选择作废，仅作描述性基线） |
| Calibration shell | `y_hat = mu_hat(MU3, outer-train 内重拟合) + alpha * z`，**`alpha = 1`** |
| Rank assert | `within_target_order(z) == within_target_order(mu + alpha*z)` 必须成立 |

---

## 6.14 禁止事项（V7.2）

- 训练新的大型 ESM / 重训 PLM
- 使用 target-ID embedding 或 family-ID lookup 作为**预测特征**
- 使用测试标签做任何模型/超参/特征选择
- 在不同 row universe 之间直接比较
- 用单一 weighting 下的 `C_B` 下结论
- 声称因果（本项目无 causal design）

---

## 6.15 冻结的声明（V7.1.1 §8）

**Supported**
```
Location calibration improves absolute family-OOD prediction.
Scale is not currently a priority bottleneck.
Current family-OOD ranker remains weak.
Static V7.0 local geometry did not improve ranking.
```

**Not supported**
```
Mutation-centered representation generally fails.
Residue-level MIPA fails.
```

> **Residue-level MIPA has not yet been tested.**（`esm2_cache` 只有 pooled 链级嵌入）
