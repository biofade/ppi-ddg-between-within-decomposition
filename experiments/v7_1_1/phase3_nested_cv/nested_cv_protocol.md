# V7.1.1 — V7.2 协议冻结

- 时间: 2026-10-14 18:12:51
- 权威 target 标签：`canonical_target_labels.csv`（`canonical.pdb_id` 按 `pair_row` 映射，与 `d['groups']` 100% 自洽）

- 权威标签 vs `d['groups']`：**完全一致 ✅**
- rows **4473** / targets **214** / families **108**

## 3.1 Outer CV

沿用 historical L4 family-level KFold（`GroupKFold` on family，5 folds，seed 由 `make_folds` 决定），并冻结为 `OUTER_FOLDS_V7_2`。

- seed **42**：5 folds，family 零重叠 = **✅**；每折 test rows = [1342, 1014, 698, 483, 936]
- seed **43**：5 folds，family 零重叠 = **✅**；每折 test rows = [1006, 913, 595, 484, 1475]
- seed **44**：5 folds，family 零重叠 = **✅**；每折 test rows = [607, 1967, 518, 371, 1010]

- 历史 L4 可复现：**YES**（`GroupKFold` 是确定性的，fold 只依赖 seed）

## 3.2 Inner CV

对每一个 outer-train 集合，**再按 family** 做 `GroupKFold`：

```
inner_folds = min(3, max(2, n_unique_family_in_outer_train))
```

规则**统一**，不随 fold 随意变化（family 少于 3 个时降为 2）。

- outer fold 0：86 families → inner **3-fold**，family 零重叠 = **✅**
- outer fold 1：86 families → inner **3-fold**，family 零重叠 = **✅**
- outer fold 2：86 families → inner **3-fold**，family 零重叠 = **✅**
- outer fold 3：87 families → inner **3-fold**，family 零重叠 = **✅**
- outer fold 4：87 families → inner **3-fold**，family 零重叠 = **✅**

## 3.3–3.4 Selection 规则（冻结）

| 项 | 规则 |
|---|---|
| **Primary** | `MacroSpearman`（per-target Spearman 的 macro 均值，**每 target n ≥ 5**） |
| Tie-break 1 | 若 `|Δρ| < 0.005` → 取较高的 `Efficiency@10` |
| Tie-break 2 | 若仍接近（`|ΔEff| < 0.005`）→ 取较高的 `C_W` |
| Tie-break 3 | 仍相同 → **更简单的模型胜**（参数更少 / 无额外模块） |
| **禁止** | 用 `Overall Pearson` 选择 ranker |
| **禁止** | 人为组合任意权重 |

## 3.5 Hyperparameters

所有候选超参必须在 **inner CV 阶段预注册**；
**禁止**在看到 outer fold 结果后新增参数。

- 写出 `protocol.json` 与 `fold_manifest.csv`

