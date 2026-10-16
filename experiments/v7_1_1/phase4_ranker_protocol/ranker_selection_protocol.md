# V7.1.1 — V7.2 协议冻结

- 时间: 2026-10-14 18:12:51
- 权威 target 标签：`canonical_target_labels.csv`（`canonical.pdb_id` 按 `pair_row` 映射，与 `d['groups']` 100% 自洽）

## 4.1 允许比较的 ranker 家族（V7.1.1 不实现）

| 代号 | 说明 |
|---|---|
| `R_BASE` | 现有 pooled / D-only 基线 |
| `R_RESIDUE_LOCAL` | residue-level 局部环境 |
| `R_RESIDUE_LOCAL_PARTNER` | 局部 + partner 上下文 |
| `R_RESIDUE_PERTURBATION` | Δh_m = h_mut − h_wt 扰动 |
| `R_RESIDUE_PERTURBATION_SET` | 扰动 × candidate-set 条件化 |

## 4.2 z 转换（**唯一函数**）

所有 ranker 必须输出**原始分数** `score_ci`，然后**统一**转换：

```python
z, collapsed = standardize_within_target(score, target_id)
# z_ci = (s_ci - mean_c(s)) / (SD_c(s) + eps)
```

**禁止**每个 ranker 自行实现该转换。

## 4.3 Collapsed prediction 实测

| z 变量 | collapsed target 数 | collapsed 行占比 |
|---|---:|---:|
| `z_R0_A0_pooled` | 38 | 1.01% |
| `z_R1_A1_base_set` | 38 | 1.01% |
| `z_R2_MLP_L0` | 38 | 1.01% |
| `z_R3_Donly_M2` | 38 | 1.01% |
| `z_R4_hist_M0` | 38 | 1.01% |

## 4.4 嵌套选择流程（冻结）

```
for each outer fold f:
    for each candidate ranker r (architecture family 事先固定):
        inner_cv_score(r) = mean over inner folds of MacroSpearman
    r* = argmax inner_cv_score  (tie-break: Eff@10 -> C_W -> simpler)
    retrain r* on full outer-train
    predict outer-test  -> OOF predictions
concatenate OOF across folds
```

**不同 outer fold 可以选择不同的超参/变体**，
但**architecture family 必须事先固定**。

## 4.5 No global winner peeking（冻结为硬规则）

**禁止**：先看完整 L4 各模型结果，再选择一个 winner。
V7.1 的 `R2_MLP_L0` 选择**正是此类泄漏**，故在 V7.2 中**作废**，
`R2_MLP_L0` 仅作为**描述性基线**保留。

