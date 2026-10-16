# V7.1 PHASE 1 — 全量覆盖审计 + 救援

- 时间: 2026-10-14 16:59:41

## PHASE 0 — 输入冻结

- 冻结 **14** 个输入 → `phase0_provenance/input_manifest.csv`

- canonical: **7523** 行 / **230** targets
- V7.0 U_LOCAL: **4473** 行 / **208** targets

### 1.2 排除机制分解

| 类别 | targets | 行数 |
|---|---:|---:|
| A 完全排除（0 个已验证位点） | **16** | **107** |
| B 部分排除（部分 variant 未通过） | **136** | **2943** |
| 合计排除 | 152 | **3050** |

> **关键**：主要损失来自 **B 类（部分排除）**，而非 16 个坏 target（仅 107 行）。

- 写出 `excluded_target_audit.csv`（16 个完全排除的 target）

### 1.3 救援（允许部分位点映射，仅用 WT 校验通过的位点）

- v2.2 行映射覆盖: 7523 / 7523

| 统计 | 值 |
|---|---:|
| n_total | 6738 |
| no_pdb | 0 |
| no_chain | 0 |
| idx_oob | 1632 |
| wt_mismatch | 2305 |
| verified | 4473 |
| no_site_ok | 2265 |

### 救援后 U_LOCAL
- **rows: 4473 / 7523 = 59.5%**（V7.0 为 4473 = 59.5%）
- **targets: 208 / 230**（V7.0 为 208）
- 其中 `partial_mapping=1` 的行: **546**（12.2%）

- target 覆盖 **90.4%**（目标 ≥90% 行覆盖：**未达成**）

### 1.5 两个 universe
- **U_FULL** = 7523 行（不需要 Local 的模型用这个）
- **U_LOCAL** = 4473 行（涉及 Local 的比较用这个）

### 1.6 Coverage-bias 检验（included vs excluded）

| 变量 | included 中位 | excluded 中位 | Mann-Whitney p | 类型 |
|---|---:|---:|---:|---|
| `n_mutations` | 14.000 | 6.000 | 0.002671 | DEPLOY (rb=0.419) |
| `n_unique_mutation_sites` | 10.000 | 4.000 | 0.0002248 | DEPLOY (rb=0.514) |
| `site_entropy` | 3.013 | 1.848 | 0.0003183 | DEPLOY (rb=0.502) |
| `alanine_fraction` | 0.454 | 0.333 | 0.8122 | DEPLOY (rb=0.033) |
| `single_mut_fraction` | 0.894 | 1.000 | 0.3744 | DEPLOY (rb=0.118) |
| `true_mean_ddg` | 0.882 | 0.928 | 0.5644 | LABEL (rb=0.081) |
| `ddg_std` | 1.169 | 1.060 | 0.505 | LABEL (rb=0.093) |

### 判定：**`COVERAGE_BIAS_PRESENT`**（3/7 个变量 p<0.05）

- 耗时 103s
