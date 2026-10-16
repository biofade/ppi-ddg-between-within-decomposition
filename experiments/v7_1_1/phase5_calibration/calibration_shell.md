# V7.1.1 — V7.2 协议冻结

- 时间: 2026-10-14 18:12:51
- 权威 target 标签：`canonical_target_labels.csv`（`canonical.pdb_id` 按 `pair_row` 映射，与 `d['groups']` 100% 自洽）

## 5.1 Primary shell

```
y_hat_ci = mu_hat_c + alpha * z_ci
```

- `mu_hat` = **MU3（target+set）**，且**必须在 outer-train 内重新拟合**
- `alpha` = 见 §5.2

## 5.2 Scale policy：`alpha = 1` vs `alpha = train 残差尺度`

因 oracle σ 增益 ≈ 0（V7.1 PHASE 6），V7.2 **不再预测复杂 σ**。
仅比较两个简单方案，**只用现有预测**，按 **inner-family RMSE** 选：

| alpha 方案 | 定义 | inner-family RMSE |
|---|---|---:|
| `alpha=1` | — | **2.5151** |
| `alpha=train_scale` | — | **2.5078** |

- 两者 inner-family RMSE 差异 = **0.0072**
- 选择规则：差异极小（<0.01）→ **优先简单**
- **冻结：`alpha=1`**

## 5.3 Secondary full-LSR（diagnostic only）

保留 `mu + sigma_set * z` 作为**次级诊断**，**不作为 V7.2 主架构**。

## 5.4 Calibration 不得改变排序（**永久 assert**）

```python
assert within_target_order(raw_z) == within_target_order(mu + alpha*z)
```

- 实测（alpha > 0）：**PASS ✅**

