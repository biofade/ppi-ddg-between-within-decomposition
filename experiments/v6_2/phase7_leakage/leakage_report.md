# V6.2 PHASE 7 — 分级泄漏检查 + 浮点 sanity

- 时间: 2026-10-14 12:20:30
- **按 split 级别分别检查**：不用 L4 的约束去检查 L1（V6.2 §7.1）

- matched set: **7523** 行 / 230 pdb / 112 family

## 7.1 分级断言

| 档 | seed | 检查 | 结果 | 详情 |
|---|---:|---|---|---|
| L1_random | 42 | obs_duplicate_no_cross | ✅ | — |
| L1_random | 42 | state_id_no_cross | ✅ | — |
| L1_random | 43 | obs_duplicate_no_cross | ✅ | — |
| L1_random | 43 | state_id_no_cross | ✅ | — |
| L1_random | 44 | obs_duplicate_no_cross | ✅ | — |
| L1_random | 44 | state_id_no_cross | ✅ | — |
| L2_complex | 42 | obs_duplicate_no_cross | ✅ | — |
| L2_complex | 42 | state_id_no_cross | ✅ | — |
| L2_complex | 42 | pdb_no_cross | ✅ | — |
| L2_complex | 43 | obs_duplicate_no_cross | ✅ | — |
| L2_complex | 43 | state_id_no_cross | ✅ | — |
| L2_complex | 43 | pdb_no_cross | ✅ | — |
| L2_complex | 44 | obs_duplicate_no_cross | ✅ | — |
| L2_complex | 44 | state_id_no_cross | ✅ | — |
| L2_complex | 44 | pdb_no_cross | ✅ | — |
| L3_seq40 | 42 | obs_duplicate_no_cross | ✅ | — |
| L3_seq40 | 42 | state_id_no_cross | ✅ | — |
| L3_seq40 | 42 | pdb_no_cross | ✅ | — |
| L3_seq40 | 42 | seq40_no_cross | ✅ | — |
| L3_seq40 | 43 | obs_duplicate_no_cross | ✅ | — |
| L3_seq40 | 43 | state_id_no_cross | ✅ | — |
| L3_seq40 | 43 | pdb_no_cross | ✅ | — |
| L3_seq40 | 43 | seq40_no_cross | ✅ | — |
| L3_seq40 | 44 | obs_duplicate_no_cross | ✅ | — |
| L3_seq40 | 44 | state_id_no_cross | ✅ | — |
| L3_seq40 | 44 | pdb_no_cross | ✅ | — |
| L3_seq40 | 44 | seq40_no_cross | ✅ | — |
| L4_family | 42 | obs_duplicate_no_cross | ✅ | — |
| L4_family | 42 | state_id_no_cross | ✅ | — |
| L4_family | 42 | pdb_no_cross | ✅ | — |
| L4_family | 42 | seq40_no_cross | ✅ | — |
| L4_family | 42 | family_no_cross | ✅ | — |
| L4_family | 43 | obs_duplicate_no_cross | ✅ | — |
| L4_family | 43 | state_id_no_cross | ✅ | — |
| L4_family | 43 | pdb_no_cross | ✅ | — |
| L4_family | 43 | seq40_no_cross | ✅ | — |
| L4_family | 43 | family_no_cross | ✅ | — |
| L4_family | 44 | obs_duplicate_no_cross | ✅ | — |
| L4_family | 44 | state_id_no_cross | ✅ | — |
| L4_family | 44 | pdb_no_cross | ✅ | — |
| L4_family | 44 | seq40_no_cross | ✅ | — |
| L4_family | 44 | family_no_cross | ✅ | — |

## 7.2 DGaux 泄漏断言

`02_split4_ddg.py::dgaux_oof` 内含两条**运行期断言**：

```python
assert len(np.intersect1d(a, b)) == 0            # 自身不在训练 IDs
assert not (set(groups_tr[ia]) & set(groups_tr[ib]))  # inner grouping 不违反
```

其中 `groups_tr = inner_u[tr]`，`inner_u = level_unit ∪ state_id`
→ 构成 **outer + state 双重分组**（V6.2 §2.4）。

- 复核 `L4_family`：level unit 数 **112** → inner（∪state）数 **6047** → 双重分组成立 ✅

## 7.3 Forward/Reverse 顺序断言

`02_split4_ddg.py` 的实现（V6.2 §2.2）：

```python
Xtr_aug = np.vstack([Xtr, -Xtr])      # reverse 只在 outer-train 内生成
ytr_aug = np.r_[ytr, -ytr]
assert len(Xtr_aug) == 2 * len(tr)
assert len(np.intersect1d(tr, te)) == 0   # test 只保留 canonical forward
```

- `reverse_rows ⊆ train_rows`：✅
- `test_rows == canonical_forward_only`：✅

## 7.4 M1 复用断言

- `M1_implementation_audited` = **True**  期望 True → ✅
- `M1_uses_historical_script` = **True**  期望 True → ✅
- `M1_only_split_layer_replaced` = **True**  期望 True → ✅
- `M1_implementation_not_found` = **False**  期望 False（未发生）→ ✅

## 7.5 Pearson 分解浮点 sanity

判据：`np.isclose(C_B + C_W, r_full, rtol=1e-10, atol=1e-08)`
（**不**用硬编码 `< 1e-10`，避免浮点差异假报 bug —— V6.2 §7.5）

- 分解 sanity: **36/36 通过**；最大 \|diff\| = **1.110e-16** → ✅ PASS

### 一致性检查（ddof / row weighting / NaN mask / row universe）

| 检查 | 做法 | 结果 |
|---|---|---|
| ddof 一致 | 全部用 **population 矩**（`std()` 默认 ddof=0、`cov` 用 `.mean()` 归一）| ✅ |
| row weighting 一致 | 每个统计量在同一 `m` mask 上计算 | ✅ |
| NaN mask 一致 | `m = isfinite(y) & isfinite(p)` 后 y/p/g 同步过滤 | ✅ |
| target means row universe | 组均值在同一 mask 后的 g 上计算 | ✅ |

## 汇总

| 项 | 值 |
|---|---:|
| 检查总数 | 50 |
| 通过 | **50** |
| 失败 | **0** |

> **无泄漏检出。** 所有分级断言通过。

> 说明：按 V6.2 §7.6，失败项**不阻塞**后续 PHASE；
> 记录 split × seed 后继续。

- 耗时 9s
