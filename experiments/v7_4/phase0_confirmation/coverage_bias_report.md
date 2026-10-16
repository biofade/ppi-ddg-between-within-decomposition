# V7.4 PHASE 0b — U_PHYS_STRICT 覆盖偏倚审计（§3）

- 时间: 2026-10-15 08:55:09

## 0b.1 行级覆盖

| 项 | 值 |
|---|---:|
| canonical rows | 4473 |
| **included** | **3351** |
| excluded | **1122** |
| included targets | **166** |
| excluded-only targets | **176** |

### 排除原因（行级，§10）

| reason | rows | 占比 |
|---|---:|---:|
| `ambiguous_mapping` | 720 | 64.2% |
| `wt_mismatch` | 257 | 22.9% |
| `missing_chain` | 83 | 7.4% |
| `idx_out_of_bounds` | 62 | 5.5% |
| **合计** | **1122** | 100% |

## 0b.2 逐 target 属性（included vs excluded-only targets）

| 组 | targets | rows |
|---|---:|---:|
| full_included | 38 | 386 |
| fully_excluded | 48 | 167 |
| partial | 128 | 3920 |

| 变量 | included 中位 | excluded 中位 | p | 类型 |
|---|---:|---:|---:|---|
| `n_mutations` | 17.500 | 5.000 | 8.135e-09 | DEPLOY (rb=0.594) |
| `n_unique_mutation_sites` | 10.000 | 3.000 | 2.506e-06 | DEPLOY (rb=0.485) |
| `single_mut_fraction` | 0.894 | 1.000 | 0.9835 | DEPLOY (rb=0.002) |
| `alanine_fraction` | 0.474 | 0.250 | 0.08999 | DEPLOY (rb=0.174) |
| `site_entropy` | 3.120 | 1.585 | 2.316e-06 | DEPLOY (rb=0.487) |
| `true_mean_ddg` | 0.888 | 0.885 | 0.8378 | LABEL (rb=0.021) |
| `ddg_std` | 1.217 | 0.965 | 0.07284 | LABEL (rb=0.185) |

## 0b.3 结构属性（PDB 链数、序列长度）

| 变量 | included 中位 | excluded 中位 | p | 类型 |
|---|---:|---:|---:|---|
| `n_pdb_chains` | 3.0 | 3.0 | 0.2258 | DEPLOY (rb=0.110) |
| `max_chain_len` | 239.0 | 221.0 | 0.2915 | DEPLOY (rb=0.100) |

## 0b.4 Family 分布

- included 涉及 family：**86** / 全体 **108**
- **完全无覆盖的 family：22** ['L0', 'L11', 'L113', 'L115', 'L120', 'L123', 'L38', 'L39', 'L5', 'L57']

## 0b.5 界面距离可得性

| | 有 BLOCK L（可算环境） | 无 BLOCK L |
|---|---:|---:|
| included | 3351 | 0 |
| excluded | 1122 | 0 |

## 0b.6 判定

- DEPLOYABLE 变量中 p<0.05 的个数：**3 / 7**
- **`COVERAGE_BIAS_PRESENT`**

> 标签派生量（**仅诊断，不用于选择**）：
>   - `true_mean_ddg`: included 中位 0.888 vs excluded 0.885（p=0.8378）
>   - `ddg_std`: included 中位 1.217 vs excluded 0.965（p=0.07284）

- 耗时 17s
