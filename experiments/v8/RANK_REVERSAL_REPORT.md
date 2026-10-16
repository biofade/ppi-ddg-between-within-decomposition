# V8 Model-Rank Reversal Report

> §31 预注册实践并列容差：Overall Pearson 与 Macro Spearman 之差 < 0.01 视为实践并列（不数作科学反转）。

> 秩一致使用 **Spearman / Kendall**（1.0 = 完全相同的模型排序）。

### universe `U_CANONICAL`（n=160 模型）

| 指标 | top-1 | top-3 |
|---|---|---|
| `overall_pearson` | L1_random_A1_base_set | L1_random_A1_base_set, L1_random_A4_base_set_local, L1_random_A2_base_global |
| `rmse` | L1_random_A1_base_set | L1_random_A1_base_set, L1_random_A4_base_set_local, L1_random_A2_base_global |
| `target_mean_pearson` | R0_Base_s42 | R0_Base_s42, R0_Base_s44, R3_Base_dh_s44 |
| `macro_spearman` | R1_Base_F0_s42 | R1_Base_F0_s42, R0_Base_s42, R2_Base_Residue_s42 |
| `macro_kendall` | R1_Base_F0_s42 | R1_Base_F0_s42, R0_Base_s42, R1_Base_F0_s44 |
| `C_W` | R0_Base_s42 | R0_Base_s42, R1_Base_F0_s42, F8_Full_s42 |
| `efficiency@10` | L1_random_A3_base_local | L1_random_A3_base_local, L1_random_A0_base, L1_random_A1_base_set |
| `ndcg@10` | L1_random_A0_base | L1_random_A0_base, L1_random_A1_base_set, L1_random_A3_base_local |
| `hit@10` | L4_family_L3_target_balanced | L4_family_L3_target_balanced, L4_family_L3_target_balanced_s43, z_R3_Donly_M2 |

| 指标对 | Spearman 秩一致 | Kendall 秩一致 | 跨轴 |
|---|---:|---:|:--:|
| `overall_pearson` vs `rmse` | +0.874 | +0.749 | ✅ |
| `rmse` vs `target_mean_pearson` | +0.415 | +0.288 | ✅ |
| `macro_spearman` vs `overall_pearson` | +0.181 | +0.123 | ✅ |
| `macro_spearman` vs `rmse` | -0.044 | -0.004 | ✅ |
| `macro_spearman` vs `target_mean_pearson` | -0.050 | -0.043 | ✅ |
| `macro_spearman` vs `ndcg@10` | +0.379 | +0.321 | ✅ |
| `macro_kendall` vs `overall_pearson` | +0.156 | +0.110 | ✅ |
| `macro_kendall` vs `rmse` | -0.065 | -0.016 | ✅ |
| `macro_kendall` vs `target_mean_pearson` | -0.068 | -0.057 | ✅ |
| `macro_kendall` vs `ndcg@10` | +0.356 | +0.300 | ✅ |
| `C_W` vs `overall_pearson` | +0.807 | +0.644 | ✅ |
| `C_W` vs `rmse` | +0.582 | +0.486 | ✅ |
| `C_W` vs `target_mean_pearson` | +0.320 | +0.217 | ✅ |
| `C_W` vs `efficiency@10` | -0.215 | -0.148 | ✅ |
| `C_W` vs `ndcg@10` | -0.052 | -0.019 | ✅ |
| `C_W` vs `hit@10` | -0.387 | -0.267 | ✅ |
| `efficiency@10` vs `overall_pearson` | -0.009 | -0.000 | ✅ |
| `efficiency@10` vs `rmse` | +0.166 | +0.133 | ✅ |
| `efficiency@10` vs `target_mean_pearson` | +0.070 | +0.091 | ✅ |
| `efficiency@10` vs `macro_spearman` | +0.188 | +0.154 | ✅ |
| `efficiency@10` vs `macro_kendall` | +0.159 | +0.133 | ✅ |
| `ndcg@10` vs `overall_pearson` | +0.093 | +0.079 | ✅ |
| `ndcg@10` vs `rmse` | +0.243 | +0.204 | ✅ |
| `ndcg@10` vs `target_mean_pearson` | -0.045 | -0.014 | ✅ |
| `hit@10` vs `overall_pearson` | -0.259 | -0.179 | ✅ |
| `hit@10` vs `rmse` | -0.114 | -0.089 | ✅ |
| `hit@10` vs `target_mean_pearson` | -0.225 | -0.161 | ✅ |
| `hit@10` vs `macro_spearman` | -0.189 | -0.117 | ✅ |
| `hit@10` vs `macro_kendall` | -0.192 | -0.125 | ✅ |

**最大秩反转（U_CANONICAL，按秩差排序，前 6）**

| 指标对 | 模型 | 族 | rank_a | rank_b | val_a | val_b |
|---|---|---|---:|---:|---:|---:|
| `overall_pearson` vs `macro_spearman` | `L1_random_A6_global_set_local` | V6_2_split4 | 7 | 160 | +0.6663 | -0.0626 |
| `overall_pearson` vs `macro_spearman` | `L1_random_A5_set_local` | V6_2_split4 | 6 | 158 | +0.6682 | -0.0527 |
| `overall_pearson` vs `macro_spearman` | `L4_family_L2_pairwise` | V6_2_split4 | 153 | 24 | -0.0984 | +0.1052 |
| `overall_pearson` vs `macro_spearman` | `s_R2_MLP_L0` | V7.1_rank_scaled | 150 | 26 | -0.0942 | +0.1051 |
| `overall_pearson` vs `macro_spearman` | `p_P0_raw_ranker` | V7.1_LSR | 152 | 28 | -0.0942 | +0.1051 |
| `overall_pearson` vs `macro_spearman` | `L4_family_L0_absolute` | V6_2_split4 | 151 | 27 | -0.0942 | +0.1051 |

### universe `U_MAS`（n=42 模型）

| 指标 | top-1 | top-3 |
|---|---|---|
| `overall_pearson` | L1_random_M1_s42 | L1_random_M1_s42, L1_random_M1_s43, L1_random_M1_s44 |
| `rmse` | L1_random_M1_s42 | L1_random_M1_s42, L1_random_M1_s43, L1_random_M1_s44 |
| `target_mean_pearson` | L1_random_M1_s42 | L1_random_M1_s42, L1_random_M1_s43, L1_random_M1_s44 |
| `macro_spearman` | L1_random_M0_s43 | L1_random_M0_s43, L1_random_M2_s43, L1_random_M2_s43 |
| `macro_kendall` | L1_random_M2_s43 | L1_random_M2_s43, L1_random_M0_s43, L1_random_M2_s43 |
| `C_W` | L1_random_M1_s42 | L1_random_M1_s42, L1_random_M1_s43, L1_random_M2_s42 |
| `efficiency@10` | L1_random_M1_s43 | L1_random_M1_s43, L1_random_M2_s43, L1_random_M0_s43 |
| `ndcg@10` | L1_random_M0_s43 | L1_random_M0_s43, L1_random_M2_s43, L1_random_M1_s43 |
| `hit@10` | L1_random_M1_s43 | L1_random_M1_s43, L1_random_M1_s42, L1_random_M0_s44 |

| 指标对 | Spearman 秩一致 | Kendall 秩一致 | 跨轴 |
|---|---:|---:|:--:|
| `overall_pearson` vs `rmse` | +0.961 | +0.884 | ✅ |
| `rmse` vs `target_mean_pearson` | +0.899 | +0.761 | ✅ |
| `macro_spearman` vs `overall_pearson` | +0.760 | +0.612 | ✅ |
| `macro_spearman` vs `rmse` | +0.706 | +0.533 | ✅ |
| `macro_spearman` vs `target_mean_pearson` | +0.694 | +0.517 | ✅ |
| `macro_spearman` vs `ndcg@10` | +0.823 | +0.668 | ✅ |
| `macro_kendall` vs `overall_pearson` | +0.583 | +0.461 | ✅ |
| `macro_kendall` vs `rmse` | +0.531 | +0.382 | ✅ |
| `macro_kendall` vs `target_mean_pearson` | +0.555 | +0.412 | ✅ |
| `macro_kendall` vs `ndcg@10` | +0.729 | +0.591 | ✅ |
| `C_W` vs `overall_pearson` | +0.902 | +0.779 | ✅ |
| `C_W` vs `rmse` | +0.853 | +0.738 | ✅ |
| `C_W` vs `target_mean_pearson` | +0.755 | +0.582 | ✅ |
| `C_W` vs `efficiency@10` | +0.798 | +0.640 | ✅ |
| `C_W` vs `ndcg@10` | +0.708 | +0.556 | ✅ |
| `C_W` vs `hit@10` | +0.273 | +0.224 | ✅ |
| `efficiency@10` vs `overall_pearson` | +0.892 | +0.731 | ✅ |
| `efficiency@10` vs `rmse` | +0.908 | +0.763 | ✅ |
| `efficiency@10` vs `target_mean_pearson` | +0.870 | +0.705 | ✅ |
| `efficiency@10` vs `macro_spearman` | +0.730 | +0.570 | ✅ |
| `efficiency@10` vs `macro_kendall` | +0.602 | +0.466 | ✅ |
| `ndcg@10` vs `overall_pearson` | +0.867 | +0.726 | ✅ |
| `ndcg@10` vs `rmse` | +0.862 | +0.703 | ✅ |
| `ndcg@10` vs `target_mean_pearson` | +0.876 | +0.728 | ✅ |
| `hit@10` vs `overall_pearson` | +0.302 | +0.231 | ✅ |
| `hit@10` vs `rmse` | +0.244 | +0.143 | ✅ |
| `hit@10` vs `target_mean_pearson` | +0.362 | +0.312 | ✅ |
| `hit@10` vs `macro_spearman` | +0.504 | +0.382 | ✅ |
| `hit@10` vs `macro_kendall` | +0.646 | +0.496 | ✅ |

**最大秩反转（U_MAS，按秩差排序，前 6）**

| 指标对 | 模型 | 族 | rank_a | rank_b | val_a | val_b |
|---|---|---|---:|---:|---:|---:|
| `overall_pearson` vs `macro_spearman` | `L4_family_M1_s43` | V6_2_split4 | 42 | 17 | +0.0800 | +0.0534 |
| `overall_pearson` vs `macro_spearman` | `L4_family_M2_s43` | V6_2_split4 | 41 | 23 | +0.0882 | +0.0503 |
| `overall_pearson` vs `macro_spearman` | `L4_family_M0_s43` | V6_2_split4 | 40 | 22 | +0.0882 | +0.0503 |
| `overall_pearson` vs `macro_spearman` | `L2_complex_M1_s42` | V6_2_split4 | 13 | 30 | +0.2464 | +0.0384 |
| `overall_pearson` vs `macro_spearman` | `L3_seq40_M1_s42` | V6_2_split4 | 14 | 31 | +0.2464 | +0.0384 |
| `overall_pearson` vs `macro_spearman` | `L1_random_M1_s44` | V6_2_split4 | 3 | 16 | +0.6881 | +0.0541 |
