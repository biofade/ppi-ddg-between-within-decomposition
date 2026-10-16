- 读入 **3** 行 <- `experiments\v5_0\p0_benchmark\baseline_repro_results.csv`
- 读入 **12** 行 <- `experiments\v5_0\p2_state_potential\p2_results.csv`
- 读入 **9** 行 <- `experiments\v5_0\p3_balancing\p3_results.csv`
- 读入 **21** 行 <- `experiments\v5_0\p4_conformal\p4_results.csv`

## 写出 `all_experiments.csv`：**45** 行

| family | 行数 | variants |
|---|---:|---|
| P2_state | 12 | P2-0_pair_plain, P2-1_pair_dgaux, P2-2_state_absolute, P2-3_state_centered |
| P3_balancing | 9 | W0_equal, W1_inv_sqrt, W2_inv_n |
| P4_conformal | 21 | conformal_global_a10, conformal_global_a20, conformal_global_a25, mondrian_a10, mondrian_a20, mondrian_a25, plain |
| gateA_repro | 3 | C+DGaux_fwd+rev |

## 各 variant 的 3-seed 均值（test = DEVELOPMENT_BENCHMARK）

| family | variant | Pearson | Spearman | macro-PDB ρ | macro-PDB r | AUROC | RMSE | P@10 | EF@1% | Hit10(ΔΔG<0) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| P2_state | `P2-0_pair_plain` | +0.3535 | +0.1288 | +0.0933 | +0.1301 | 0.5572 | 2.2112 | 0.9000 | 8.4030 | 0.9000 |
| P2_state | `P2-1_pair_dgaux` | +0.3819 | +0.1876 | +0.1199 | +0.1441 | 0.5949 | 2.1513 | 1.0000 | 9.9035 | 1.0000 |
| P2_state | `P2-2_state_absolute` | +0.0594 | +0.1033 | -0.0055 | +0.0384 | 0.5436 | 2.9764 | 0.0000 | 0.0000 | 0.2000 |
| P2_state | `P2-3_state_centered` | -0.0820 | -0.1188 | -0.0298 | -0.0021 | 0.4314 | 2.5860 | 0.1000 | 0.9003 | 0.2000 |
| P3_balancing | `W0_equal` | +0.3862 | +0.1870 | +0.0956 | +0.1334 | 0.5935 | 2.1484 | 0.9667 | 9.3033 | 0.9667 |
| P3_balancing | `W1_inv_sqrt` | +0.3423 | +0.0999 | +0.0452 | +0.0717 | 0.5262 | 2.1985 | 0.9667 | 9.6034 | 0.9667 |
| P3_balancing | `W2_inv_n` | +0.2173 | +0.0383 | -0.0011 | +0.0023 | 0.4928 | 2.3387 | 0.9000 | 9.0032 | 0.9000 |
| P4_conformal | `conformal_global_a10` | +0.3721 | +0.1831 | +0.0834 | +0.1326 | 0.5933 | 7.2331 | 0.9667 | 9.6034 | 0.9667 |
| P4_conformal | `conformal_global_a20` | +0.3721 | +0.1831 | +0.0834 | +0.1326 | 0.5933 | 4.7910 | 0.9667 | 9.6034 | 0.9667 |
| P4_conformal | `conformal_global_a25` | +0.3721 | +0.1831 | +0.0834 | +0.1326 | 0.5933 | 4.0785 | 0.9667 | 9.6034 | 0.9667 |
| P4_conformal | `mondrian_a10` | +0.1558 | +0.1135 | +0.1044 | +0.1585 | 0.5515 | 5.7254 | 0.9667 | 9.6034 | 0.9667 |
| P4_conformal | `mondrian_a20` | +0.1451 | +0.1180 | +0.1121 | +0.1674 | 0.5552 | 4.6038 | 0.9667 | 9.3033 | 0.9667 |
| P4_conformal | `mondrian_a25` | +0.1461 | +0.1207 | +0.1144 | +0.1671 | 0.5572 | 4.2351 | 0.9000 | 9.0032 | 0.9000 |
| P4_conformal | `plain` | +0.3721 | +0.1831 | +0.0834 | +0.1326 | 0.5933 | 2.1640 | 0.9667 | 9.6034 | 0.9667 |
| gateA_repro | `C+DGaux_fwd+rev` | +0.3862 | +0.1870 | +0.0956 | +0.1334 | 0.5935 | 2.1484 | 0.9667 | 9.3033 | 0.9667 |

## 决策标签（规则来自 `v5_common.decide`，非人工判断）

- 基线 = `C+DGaux_fwd+rev`（Gate A 复现的 V4 最强配置）

| variant | ΔPearson | Δmacro-PDB ρ | ΔP@10 | KEEP(vs 全局判据) |
|---|---:|---:|---:|---|
| `P2-0_pair_plain` | -0.0327 | -0.0023 | -0.0667 | DROP/INCONCLUSIVE |
| `P2-1_pair_dgaux` | -0.0043 | +0.0244 | +0.0333 | DROP/INCONCLUSIVE |
| `P2-2_state_absolute` | -0.3268 | -0.1011 | -0.9667 | DROP/INCONCLUSIVE |
| `P2-3_state_centered` | -0.4681 | -0.1254 | -0.8667 | DROP/INCONCLUSIVE |
| `W0_equal` | +0.0000 | +0.0000 | +0.0000 | DROP/INCONCLUSIVE |
| `W1_inv_sqrt` | -0.0438 | -0.0503 | +0.0000 | DROP/INCONCLUSIVE |
| `W2_inv_n` | -0.1688 | -0.0967 | -0.0667 | DROP/INCONCLUSIVE |
| `conformal_global_a10` | -0.0141 | -0.0121 | +0.0000 | DROP/INCONCLUSIVE |
| `conformal_global_a20` | -0.0141 | -0.0121 | +0.0000 | DROP/INCONCLUSIVE |
| `conformal_global_a25` | -0.0141 | -0.0121 | +0.0000 | DROP/INCONCLUSIVE |
| `mondrian_a10` | -0.2303 | +0.0088 | +0.0000 | DROP/INCONCLUSIVE |
| `mondrian_a20` | -0.2411 | +0.0165 | +0.0000 | DROP/INCONCLUSIVE |
| `mondrian_a25` | -0.2401 | +0.0188 | -0.0667 | DROP/INCONCLUSIVE |
| `plain` | -0.0141 | -0.0121 | +0.0000 | DROP/INCONCLUSIVE |

## `paired_comparisons.csv`：**24** 行

| baseline | candidate | metric | Δ | p | BH q | decision |
|---|---|---|---:|---:|---:|---|
| `plain` | `conformal_global_a10` | p_at_10 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a10` | ef_1pct | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a10` | hit10_ddg_m05 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a20` | p_at_10 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a20` | ef_1pct | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a20` | hit10_ddg_m05 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a25` | p_at_10 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a25` | ef_1pct | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `conformal_global_a25` | hit10_ddg_m05 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a10` | p_at_10 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a10` | ef_1pct | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a10` | hit10_ddg_m05 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a20` | p_at_10 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a20` | ef_1pct | -0.3001 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a20` | hit10_ddg_m05 | +0.0000 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a25` | p_at_10 | -0.0667 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a25` | ef_1pct | -0.6002 | 1.0000 | 1.0000 | **DROP** |
| `plain` | `mondrian_a25` | hit10_ddg_m05 | -0.0667 | 1.0000 | 1.0000 | **DROP** |
| `P2-1_pair_dgaux` | `P2-0_pair_plain` | overall_pearson | -0.0284 | 0.1456 | 0.8736 | **DROP** |
| `P2-1_pair_dgaux` | `P2-0_pair_plain` | macro_pdb_spearman | -0.0266 | 0.4328 | 1.0000 | **DROP** |
| `P2-1_pair_dgaux` | `P2-2_state_absolute` | overall_pearson | -0.3225 | 0.0040 | 0.0480 | **DROP** |
| `P2-1_pair_dgaux` | `P2-2_state_absolute` | macro_pdb_spearman | -0.1255 | 0.2412 | 1.0000 | **DROP** |
| `P2-1_pair_dgaux` | `P2-3_state_centered` | overall_pearson | -0.4639 | 0.0017 | 0.0400 | **DROP** |
| `P2-1_pair_dgaux` | `P2-3_state_centered` | macro_pdb_spearman | -0.1498 | 0.1456 | 0.8736 | **DROP** |

## GATE 汇总（V5 §124）

| Gate | 问题 | 结论 |
|---|---|---|
| A | V4 baseline 复现成功？ | **PASS**（3 seeds 逐位一致，Δ≤1e-6） |
| B | mutation mapping 覆盖 ≥90%？ | **FAIL**（严格 EXACT 仅 13.12%；v2 宽松口径 61.35%） → `SOLVE_MAPPING_FIRST` |
| C | site Δh 有正信号？ | **NOT_TESTED**（ESM-2 不可用，见下） |
| D | state potential 有正收益？ | **NO**（ΔPearson −0.32，Δmacro −0.13） |
| E | Conformal Top-K 有改善？ | 见 P4 报告 |

## 资源阻塞清单（`NOT_RUN_RESOURCE_BLOCKED`）

| 项 | 阻塞原因 |
|---|---|
| P0 外部 SOTA 评估（官方 CATH split） | CPU-only 下需为 CATH 的 WT/mutant 序列重新做 ESM-2 推理；本机 `esm`/`fair-esm` **未安装**、无本地权重 |
| P1 MIPA full | 同上（需逐残基 ESM-2 表示）+ GATE B 覆盖率不足 |
| P1B inverse folding (ProteinMPNN) | 需 ProteinMPNN 权重与（建议）GPU |


