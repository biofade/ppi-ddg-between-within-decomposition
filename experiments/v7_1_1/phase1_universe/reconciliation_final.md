# V7.1.1 PHASE 1 — 208 vs 214 定案

## 1. 三个已确证事实

**(1)** `d['groups']` vs canonical `pdb_id`：一致率 **100.00%**（7523/7523） → canonical 上二者 unique 均为 **230 / 230**
**(2)** `block_L.pdb_id` vs canonical `pdb_id`：一致率 **94.21%**，不一致 **259** 行；且含 canonical 中不存在的 `1Y1K`
**(3)** 同一 4473 行集合，按上述两种 pdb_id 聚合 → **208** vs **214** 个 target

## 2. 权威 target 标签（V7.2 唯一应使用者）

- 来源：**`canonical.pdb_id` 按 `pair_row` 映射**（唯一与 `d['groups']` 100% 自洽的命名）
- rows = **4473**，targets = **214**，families = **108**

| 阈值 | target 数 |
|---|---:|
| n ≥ 1 | **214** |
| n ≥ 2 | **180** |
| n ≥ 3 | **164** |
| n ≥ 5 | **130** |
| n ≥ 10 | **85** |

- `n == 1` 的 target：**34**（这些 target 的 within 类指标无法定义）

## 3. 错位对关键指标的**实际影响**

用 `lsr_predictions`（其 target 标签来自 `d['groups']`，100% 对齐）作为预测，
分别用 **(A) `d['groups']` 标签** 与 **(B) 旧 `block_L.pdb_id` 标签** 计算指标。

- ⚠️ `lsr_predictions.npz` 未含 `pair_row`，改用「pair_row → 旧 pdb_id」的按序拼接近似（已标注为近似）。
| 指标 | (A) canonical/`groups` 标签 | (B) 旧 `block_L.pdb_id` 标签 | Δ |
|---|---:|---:|---:|
| P0_raw · C_W | +0.0441 | +0.0474 | -0.0033 |
| P0_raw · macro_spearman | +0.1117 | +0.0945 | +0.0172 |
| P0_raw · target_mean_r | +0.1908 | +0.1224 | +0.0684 |
| P0_raw · n_target | 214 | 208 | +6 |
| P3_fullLSR · C_W | +0.1554 | +0.1614 | -0.0059 |
| P3_fullLSR · macro_spearman | +0.1117 | +0.1027 | +0.0090 |
| P3_fullLSR · target_mean_r | +0.2332 | +0.1840 | +0.0492 |
| P3_fullLSR · n_target | 214 | 208 | +6 |

## 4. 定案

| 判定 | 值 |
|---|---|
| **Cause** | **`TARGET_ID_NAMESPACE_MISMATCH`** |
| 208 的来源 | `block_L*.csv` 的 `pdb_id`（源自 `mutation_sites`，94.21% 对齐） |
| 214 的来源 | `d['groups']`（源自 `build_pair_features`，**100% 对齐 canonical**） |
| 行数 | 两者均为 **4473**（同一 row 集合） |
| 是否为 feature intersection | **否** |
| 是否为 metric eligibility | **否** |
| 是否为纯文本错误 | **否**（是两个真实不同的命名空间） |

### 受影响 / 未受影响

| 环节 | 是否受影响 | 说明 |
|---|---|---|
| LSR 的 target 标签（PHASE 3/4/5/6 聚合） | **未受影响** | 用 `d['groups']` |
| S/G 特征取回（`reindex(pdb_rows)`） | **未受影响** | S/G 索引同源自 canonical |
| BLOCK L 特征与 target 的对应 | **受影响（259 行）** | 该 259 行的 local 特征归错 target |
| LSR/Oracle 的数值结论 | **未受影响** | target 标签自始用 `d['groups']` |

### ⚠️ V7.1 报告须更正的表述

1. 「4473 rows / **208 targets**」应改为「4473 rows / **214 targets**」
   （208 是旧命名空间下的计数）
2. 「**90.4%** target 覆盖」应改为 **214/230 = 93.0%**
3. V7.1 PHASE 1 报告中「救援未增加任何行」的**结论仍成立**，
   但其 `208` 这个计数应替换为 `214`

- 耗时 10s
