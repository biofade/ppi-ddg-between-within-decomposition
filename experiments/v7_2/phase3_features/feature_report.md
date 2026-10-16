# V7.2 PHASE 3–5 — 特征块 + ranker 家族 + 嵌套 family-CV

- 时间: 2026-10-14 18:54:11
- learner: **LightGBM**（`LGB_PARAMS`, trees=300）— 统一，不逐族调参
- PCA: **32**，**只在 train 行 fit**
- 协议: V7.1.1 冻结 `OUTER_FOLDS_V7_2` + inner family CV

- 残基嵌入：**13** 个 target 文件，**1807** 个可用 row
- 占 U_CONTRACT 4473 行的 **40.4%**

- 成功装配 h_wt/h_mut：**1807** 行（320 维）

- **U_CONTRACT 子集（有残基嵌入）**: rows **1807** / targets **13** / families **9**
- 其中 n≥5 的 target（Macro ρ 资格）: **12**

- **F0** mutation identity: 53 维
- **F1** h_wt,h_mut: 640 维
- **F2** Δh + ‖Δh‖ + cos: 322 维
- **F3** partner context: **NOT_RUN**（需结构映射；本轮未做）
- **F4** interaction: **NOT_RUN**（依赖 F3）
- **F5** candidate-set: 25 维

## PHASE 5 — 嵌套 family-CV（冻结协议）

| R0_Base | s42 | r=+0.3145 | C_W=+0.3288 | macro_ρ=+0.3272 | Eff@10=+0.2900 | collapse=0.2% |
| R1_Base_F0 | s42 | r=+0.3241 | C_W=+0.3268 | macro_ρ=+0.3363 | Eff@10=+0.2994 | collapse=0.2% |
| R2_Base_Residue | s42 | r=+0.2512 | C_W=+0.2933 | macro_ρ=+0.2511 | Eff@10=+0.1791 | collapse=0.0% |
| R3_Base_dh | s42 | r=+0.3022 | C_W=+0.3117 | macro_ρ=+0.2554 | Eff@10=+0.2845 | collapse=0.0% |
| R0_Base | s43 | r=+0.1010 | C_W=+0.2752 | macro_ρ=+0.2127 | Eff@10=+0.1068 | collapse=0.2% |
| R1_Base_F0 | s43 | r=+0.1783 | C_W=+0.2923 | macro_ρ=+0.2052 | Eff@10=+0.2474 | collapse=0.2% |
| R2_Base_Residue | s43 | r=+0.1938 | C_W=+0.3053 | macro_ρ=+0.2411 | Eff@10=+0.1146 | collapse=0.0% |
| R3_Base_dh | s43 | r=+0.1522 | C_W=+0.2867 | macro_ρ=+0.1441 | Eff@10=+0.1495 | collapse=0.0% |
| R0_Base | s44 | r=+0.2736 | C_W=+0.2927 | macro_ρ=+0.2559 | Eff@10=+0.1046 | collapse=0.2% |
| R1_Base_F0 | s44 | r=+0.3289 | C_W=+0.3019 | macro_ρ=+0.2547 | Eff@10=+0.0726 | collapse=0.2% |
| R2_Base_Residue | s44 | r=+0.3246 | C_W=+0.3016 | macro_ρ=+0.2104 | Eff@10=+0.2456 | collapse=0.0% |
| R3_Base_dh | s44 | r=+0.2861 | C_W=+0.2900 | macro_ρ=+0.2105 | Eff@10=+0.1038 | collapse=0.0% |

## 3-seed 均值

| Ranker | Overall r | C_W | Macro ρ | Eff@10 |
|---|---:|---:|---:|---:|
| R0_Base | +0.2297 | +0.2989 | +0.2653 | +0.1671 |
| R1_Base_F0 | +0.2771 | +0.3070 | +0.2654 | +0.2065 |
| R2_Base_Residue | +0.2566 | +0.3001 | +0.2342 | +0.1798 |
| R3_Base_dh | +0.2468 | +0.2961 | +0.2033 | +0.1793 |

- 耗时 6.6 分钟
