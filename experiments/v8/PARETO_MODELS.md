# V8 Pareto-Optimal Models（§65）

> 避免过早使用复合评分；此处只给出**非支配集**。

> 目标：`RMSE` 越低、`Macro Spearman` 越高、`Efficiency@10` 越高。


## universe `U_CANONICAL`（12 个非支配模型）

| model | family | RMSE | Macro ρ | Eff@10 | Overall r |
|---|---|---:|---:|---:|---:|
| `R1_Base_F0_s42` | V7.2R1_representation | 2.9141 | +0.3388 | +0.7109 | +0.3241 |
| `R0_Base_s42` | V7.2R1_representation | 2.9853 | +0.3306 | +0.7203 | +0.3145 |
| `R2_Base_Residue_s44` | V7.2R1_representation | 2.8785 | +0.2190 | +0.6825 | +0.3246 |
| `F8_Full_s42` | V7.3_physics | 2.0178 | +0.1454 | +0.8235 | +0.4662 |
| `F3_Base_EnvPartner_s43` | V7.3_physics | 2.0547 | +0.1158 | +0.8262 | +0.4311 |
| `F6_Partner_Intra_s44` | V7.3_physics | 2.0417 | +0.1141 | +0.8236 | +0.4433 |
| `F7_F6_Burial_s43` | V7.3_physics | 2.0382 | +0.1074 | +0.8304 | +0.4461 |
| `L4_family_L2_pairwise` | V6_2_split4 | 3.4225 | +0.1052 | +0.8315 | -0.0984 |
| `p_P4_set_location` | V7.1_LSR | 2.4249 | +0.1051 | +0.8347 | +0.2636 |
| `L1_random_A0_base` | V6_2_split4 | 1.7714 | +0.0999 | +0.8534 | +0.6893 |
| `L1_random_A3_base_local` | V6_2_split4 | 1.7639 | +0.0875 | +0.8547 | +0.6922 |
| `L1_random_A1_base_set` | V6_2_split4 | 1.5464 | +0.0747 | +0.8517 | +0.7729 |

## universe `U_MAS`（5 个非支配模型）

| model | family | RMSE | Macro ρ | Eff@10 | Overall r |
|---|---|---:|---:|---:|---:|
| `L1_random_M0_s43` | V6_2_split4 | 1.7678 | +0.1003 | +0.8144 | +0.6844 |
| `L1_random_M2_s43` | V6_2_split4 | 1.7678 | +0.1003 | +0.8144 | +0.6844 |
| `L1_random_M2_s43` | V6_2_split4 | 1.7825 | +0.0993 | +0.8152 | +0.6822 |
| `L1_random_M1_s42` | V6_2_split4 | 1.7360 | +0.0876 | +0.8137 | +0.6990 |
| `L1_random_M1_s43` | V6_2_split4 | 1.7361 | +0.0765 | +0.8167 | +0.6987 |
