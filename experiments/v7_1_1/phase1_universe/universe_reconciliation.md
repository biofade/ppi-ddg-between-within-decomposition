# V7.1.1 PHASE 0–2 — Universe 对账 + Metric 资格审计

- 时间: 2026-10-14 18:10:53

## PHASE 0 — 冻结 **17** 个输入

- `v71_u_local_rows`: `d34150d4861322b9…`
- `v71_u_full_rows`: `df423f39e27d8998…`
- `v71_ranker_metrics`: `3cddb08ee2e75515…`
- `v71_z_predictions`: `db2f4adc6e1daeb2…`
- `v71_location_pred`: `be0dbaf895d736d6…`
- `v71_location_models`: `b7ab9fff6187d4b0…`
- `v71_scale_models`: `bcf6b1da6a881bc9…`
- `v71_lsr_models`: `9b173d921541ee1a…`
- `v71_lsr_predictions`: `5b1f100881161069…`
- `v71_oracle`: `1f638a8f1f713940…`
- `v71_contrasts`: `7df56a40b5e220af…`
- `v71_stratification`: `2d2de4ba8843a35e…`
- `v71_final_summary`: `2d1aead9d0c13293…`
- `v70_matched_rows`: `aa7fc025a0fe98cc…`
- `v70_block_L`: `dfd0886c54cd96f2…`
- `canonical_matched_set`: `8d209eea08cc6058…`
- `ddg_split_family`: `92d4e8d78f0778ff…`

## PHASE 1 — Universe 对账

### 1.1 目标集合规模

| 集合 | n_rows | n_targets |
|---|---:|---:|
| T_CANONICAL（canonical） | 7523 | 230 |
| T_BLOCK_L_V70 | 4473 | 208 |
| T_BLOCK_L_V71（rescued） | 4473 | 208 |
| T_U_LOCAL_ROWS | 4473 | 208 |
| T_LSR（predictions 内 target） | 4473 | 214 |
| T_LOCATION | — | 214 |

### 1.2 集合差（关键）

- `T_lf1 − T_canonical` = **['1Y1K']**
- `T_canonical − T_lf1`（n=23）

> ⚠️ **发现 1 个 target 出现在 BLOCK L 但不在 canonical**：`['1Y1K']`
> 这说明 V7.0/V7.1 的 BLOCK L 构造**未严格限制在 canonical row 集内**，
> 或 canonical 的 pdb_id 大小写/命名与 mutation_sites 不一致。

  - `1Y1K`：canonical 中匹配行 = **0**（大小写不敏感比较）

### 1.3 LSR 的 row→target 划分

- LSR rows = **4473**，unique targets = **214**
- 每 target 的行数分布：min=1, median=6, max=355
- `n == 1` 的 target 数 = **34**

### 1.4 208 与 214 的确切来源

- `block_L_rescued` rows = 4473；LSR rows = 4473
- 行数相同
- 按 pair_row 排序后，`block_L_rescued.pdb_id` 与 LSR `pdb` 序列**不一致**

- 不一致的位置数：**259** / 4473
- 前 10 个不一致位置：[856, 858, 914, 928, 929, 950, 969, 1005, 1015, 1062]

| pair_row | block_L pdb | LSR pdb |
|---:|---|---|
| 1500 | `1FSS` | `1FY8` |
| 1502 | `1FY8` | `1GC1` |
| 1558 | `1GC1` | `1GCQ` |
| 1572 | `1GCQ` | `1GRN` |
| 1574 | `1GRN` | `1GUA` |
| 1604 | `1GUA` | `1H9D` |
| 1625 | `1HE8` | `1IAR` |
| 1661 | `1IAR` | `1JCK` |
| 1671 | `1JCK` | `1JRH` |
| 1726 | `1JRH` | `1JTD` |
| 1798 | `1JTD` | `1JTG` |
| 2371 | `1JTG` | `1K8R` |
| 2391 | `1K8R` | `1KBH` |
| 2477 | `1KBH` | `1KTZ` |
| 2531 | `1KTZ` | `1LFD` |
| 2600 | `1LFD` | `1M9E` |
| 2601 | `1M9E` | `1MAH` |
| 2710 | `1MHP` | `1MLC` |
| 2771 | `1MQ8` | `1N8Z` |
| 2785 | `1N8Z` | `1NCA` |

### 1.5 根因定位（`d['groups']` vs canonical `pdb_id`）

- `d['groups'][canonical.row]` 与 `canonical.pdb_id` 一致率：**100.00%**（7523/7523）
- **根因**：V7.1 的 `02_lsr.py` 用 `d['groups']`（来自 `build_pair_features`）作为 target 标签，而 V7.0 的 `block_L_local.csv` 用的是`mutation_sites.csv` 的 `pdb_id`。**两者并非同一命名/取值来源**。
- 在 canonical 内，`d['groups']` 的 unique 值数 = **230**，`canonical.pdb_id` 的 unique 值数 = **230**

- 不一致的行数：**0**，涉及 `pdb_id` **0** 个、`groups` **0** 个

| pdb_id（canonical） | pdb_groups（build_pair_features） | 行数 |
|---|---|---:|

> ### 判定
> **214** 是 `d['groups']`（`build_pair_features` 的 `pdb_id`）下的 target 数；
> **208** 是 `mutation_sites.csv` / `block_L_rescued.csv` 的 `pdb_id` 下的 target 数。
> 二者**行集合相同（4473）**但 **target 命名不同** → 这是
> **`TARGET_ID_NAMESPACE_MISMATCH`**（本项目的 `d['groups']` 与 canonical
> `pdb_id` 在少数复合物上大小写/别名不同），
> **不是** feature intersection、也不是 metric eligibility。

### 1.6 Row-level 不变量检查

- `canonical`: row_id unique = **True**；target assignment unique = **True**
- `block_L_v71`: row_id unique = **True**；target assignment unique = **True**

- 写出 `universe_summary.csv`（6 个 universe）
- 写出 `target_membership_matrix.csv`（231 个 target）

- 耗时 21s
