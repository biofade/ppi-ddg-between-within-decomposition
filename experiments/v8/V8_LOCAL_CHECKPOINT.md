==================================
V8 LOCAL CHECKPOINT
==================================

UPSTREAM_HASH_CHANGE = 0   (24 项冻结，含 V7.3/V7.4/V7.4R2 关键制品)
Historical prediction artifacts found = 167 文件 (prediction_like 63)
Internal models found = 211 预测键 / 12 制品 / 9 模型族
Primary-usable internal models = 211

External models audited = 0   (PHASE 1/13 尚未执行 —— 见 MISSING ASSETS)
External models with official row-level predictions = 0 (未审计)
External models with official processed/split artifacts = 0 (未审计)
External models requiring inference = 0 (未审计)

TIER_1 models = 内部 202 条评价键（U_CANONICAL 160 + U_MAS 42）
TIER_2 models = 0
TIER_3 models = 0（文献矩阵 PHASE 2 尚未执行）

--------------------------------
MATCHED UNIVERSE
--------------------------------
Universe ID = U_CANONICAL  (family-OOD; canonical_target_labels by row_id)
  rows = 4473   targets = 214   families = 108   n>=5 = 130
  models evaluated = 160
Universe ID = U_MAS  (row-aligned to matched_analysis_set; 含 random/complex/family 三种 split)
  rows = 7523   targets = 230   families = 112   n>=5 = 172
  models evaluated = 42

--------------------------------
PRELIMINARY SCOREBOARD  (U_CANONICAL, 按 Macro Spearman 前 10)
--------------------------------
| Model | Overall r | RMSE | Macro rho | C_W | Eff@10 | NDCG@10 |
|---|---:|---:|---:|---:|---:|---:|
| `R1_Base_F0_s42` | +0.3241 | 2.9141 | +0.3388 | +0.3268 | +0.7109 | +0.6911 |
| `R0_Base_s42` | +0.3145 | 2.9853 | +0.3306 | +0.3288 | +0.7203 | +0.7144 |
| `R2_Base_Residue_s42` | +0.2512 | 3.0723 | +0.2682 | +0.2933 | +0.6608 | +0.6363 |
| `R2_Base_Residue_s43` | +0.1938 | 3.0658 | +0.2577 | +0.3053 | +0.6510 | +0.6377 |
| `R3_Base_dh_s42` | +0.3022 | 2.9892 | +0.2564 | +0.3117 | +0.7029 | +0.6693 |
| `R1_Base_F0_s44` | +0.3289 | 2.9619 | +0.2561 | +0.3019 | +0.6531 | +0.6218 |
| `R0_Base_s44` | +0.2736 | 3.1243 | +0.2542 | +0.2927 | +0.6496 | +0.6482 |
| `R2_Base_Residue_s44` | +0.3246 | 2.8785 | +0.2190 | +0.3016 | +0.6825 | +0.6539 |
| `R3_Base_dh_s44` | +0.2861 | 3.0602 | +0.2083 | +0.2900 | +0.6481 | +0.6229 |
| `R0_Base_s43` | +0.1010 | 3.4462 | +0.2077 | +0.2752 | +0.6519 | +0.6243 |

同表按 Overall Pearson 前 5（用于对照）:
| Model | Overall r | RMSE | Macro rho | C_W | Eff@10 | NDCG@10 |
|---|---:|---:|---:|---:|---:|---:|
| `L1_random_A1_base_set` | +0.7729 | 1.5464 | +0.0747 | +0.2563 | +0.8517 | +0.7783 |
| `L1_random_A4_base_set_local` | +0.7701 | 1.5564 | +0.0617 | +0.2506 | +0.8475 | +0.7724 |
| `L1_random_A2_base_global` | +0.7626 | 1.5748 | +0.0743 | +0.2564 | +0.8505 | +0.7759 |
| `L1_random_A3_base_local` | +0.6922 | 1.7639 | +0.0875 | +0.2804 | +0.8547 | +0.7779 |
| `L1_random_A0_base` | +0.6893 | 1.7714 | +0.0999 | +0.2873 | +0.8534 | +0.7805 |

--------------------------------
PRELIMINARY RANK AGREEMENT
--------------------------------
[U_CANONICAL]  n_models = 160
  overall_pearson      vs macro_spearman     Spearman=+0.181  Kendall=+0.123
  overall_pearson      vs efficiency@10      Spearman=-0.009  Kendall=-0.000
  overall_pearson      vs ndcg@10            Spearman=+0.093  Kendall=+0.079
  rmse                 vs macro_spearman     Spearman=-0.044  Kendall=-0.004
  macro_spearman       vs efficiency@10      Spearman=+0.188  Kendall=+0.154
  target_mean_pearson  vs macro_spearman     Spearman=-0.050  Kendall=-0.043

[U_MAS]  n_models = 42
  overall_pearson      vs macro_spearman     Spearman=+0.760  Kendall=+0.612
  overall_pearson      vs efficiency@10      Spearman=+0.892  Kendall=+0.731
  overall_pearson      vs ndcg@10            Spearman=+0.867  Kendall=+0.726
  rmse                 vs macro_spearman     Spearman=+0.706  Kendall=+0.533
  macro_spearman       vs efficiency@10      Spearman=+0.730  Kendall=+0.570
  target_mean_pearson  vs macro_spearman     Spearman=+0.694  Kendall=+0.517

Potential material rank reversals:
  U_CANONICAL: `L1_random_A6_global_set_local` (V6_2_split4) rank(Overall)=7 vs rank(Macro)=160  (r=+0.666, rho=-0.063)
  U_CANONICAL: `L1_random_A5_set_local` (V6_2_split4) rank(Overall)=6 vs rank(Macro)=158  (r=+0.668, rho=-0.053)
  U_CANONICAL: `L4_family_L2_pairwise` (V6_2_split4) rank(Overall)=153 vs rank(Macro)=24  (r=-0.098, rho=+0.105)
  U_CANONICAL: `s_R2_MLP_L0` (V7.1_rank_scaled) rank(Overall)=150 vs rank(Macro)=26  (r=-0.094, rho=+0.105)
  U_CANONICAL: `p_P0_raw_ranker` (V7.1_LSR) rank(Overall)=152 vs rank(Macro)=28  (r=-0.094, rho=+0.105)

--------------------------------
TARGET-WEIGHTING (preliminary, H2)
--------------------------------
  C_B 极差: median=0.0183  max=0.0299  符号变化模型数=35/109
  C_W 极差: median=0.0624  max=0.1368  符号变化模型数=0/109
  [U_CANONICAL] 模型排名 ROW vs EQUAL: Spearman=+0.953  Kendall=+0.844  最大秩位移=20
  [U_MAS] 模型排名 ROW vs EQUAL: Spearman=+0.921  Kendall=+0.828  最大秩位移=12

--------------------------------
NOVELTY MATRIX（初判，PHASE 2 完整文献审计尚未执行）
--------------------------------
Existing work already covers (需 PHASE 2 逐一核实并标注):
  per-target / per-structure Spearman, family-OOD splits, Top-K prioritization,
  baseline/mutation decomposition, mutation microenvironment, partner-aware repr.
Potentially distinctive V8 components (待文献核实后再定):
  严格 task-specific family-OOD + matched-universe + 三轴能力分解 +
  目标加权审计 + 模型秩反转 + 决策对齐 + benchmark-guided remedy 的**系统组合**
Potential novelty already occupied = 未审计（不得预先声称任何 'first'）

--------------------------------
GPU
--------------------------------
GPU_REQUIRED_SO_FAR = NO
 理由：本轮全部结论均来自**已有内部预测**，CPU 即可；
       外部模型资产审计（PHASE 13）尚未开始，依 §75 不在 PHASE 13 之前决定 GPU。

--------------------------------
MISSING ASSETS
--------------------------------
  - PHASE 1 外部资源审计未执行（PPIformer / CATH-ddG / PreFold-dG / HERMES / Prompt-DDG / RDE-Network 的官方 row-level 预测）
  - PHASE 2 文献基准矩阵 + 新颖性矩阵未执行（tabulate 未安装，但已改用 to_string 规避）
  - PHASE 4 符号约定审计未执行（当前假设 +ddG = 结合变弱）
  - PHASE 10 split inflation 未定稿（U_MAS 已具备 random/complex/family 三种 split 预测）
  - PHASE 11–12 选择效用 / 决策对齐未定稿
  - PHASE 15 5000 target-cluster bootstrap 未执行
  - 特征工程外的外部方法（FoldX/Flex-ddG/DDGPred/ProteinMPNN/ESM-IF 预测文件）未在项目内找到

--------------------------------
NEXT_LOCAL_ACTION
--------------------------------
执行 PHASE 10（split inflation）：用同一批 V6.2 模型在 U_MAS 上
对比 random / complex(seq40) / family 三种 split 的 Overall r 与 Macro rho，
量化「宽松划分抬高表观性能」的幅度——这是 H4 的直接检验。

