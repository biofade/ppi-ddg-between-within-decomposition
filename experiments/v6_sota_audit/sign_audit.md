# V6.0 PHASE 2.1 — ΔΔG 符号一致性审计（手工可复核）

- 时间: 2026-10-13 22:45:14
- 物理常数: **R·T = 0.592485** kcal/mol（R = 1.98720425864083e-3, T = 298.15 K）
- 约定: `ΔG = RT·ln(KD[M])`；`ΔΔG = ΔG_mut − ΔG_WT`；**ΔΔG > 0 = 削弱结合，ΔΔG < 0 = 增强结合（beneficial）**

## 1. 独立复算：本地 SKEMPI v2.0 原始表（分号分隔，第三方真值）

- 读入 `data/v2.4/skempi_v2_cleaned.csv`: **6993 行 × 29 列**（分隔符 `;`）
- 可复算行（两个 KD 均为正）: **6519**

### 逐步手工复算（取 6 条，含 1 条典型 destabilizing）

| #Pdb | Mutation(cleaned) | KD_WT (M) | KD_mut (M) | ΔG_WT | ΔG_mut | ΔΔG_calc | ΔΔG<0? |
|---|---|---:|---:|---:|---:|---:|---|
| 1CSE_E_I | LI38G | 1.120e-12 | 5.260e-11 | -16.3038 | -14.0231 | **+2.2807** | 否（削弱） |
| 1CSE_E_I | LI38S | 1.120e-12 | 8.330e-12 | -16.3038 | -15.1150 | **+1.1888** | 否（削弱） |
| 1CSE_E_I | LI38P | 1.120e-12 | 1.020e-07 | -16.3038 | -9.5380 | **+6.7658** | 否（削弱） |
| 4NKQ_C_AB | RC11D | 5.800e-11 | 4.100e-11 | -13.9652 | -14.1707 | **-0.2055** | ✅ 是 |
| 1GC1_G_C | KC75A | 3.600e-09 | 4.700e-09 | -11.5193 | -11.3613 | **+0.1580** | 否（削弱） |
| 1FFW_A_B | TA86I | 1.200e-06 | 7.000e-07 | -8.0775 | -8.3968 | **-0.3193** | ✅ 是 |

**读法**：`LI38G` 的 KD 从 1.12e−12 升到 5.26e−11（结合**变弱**）→ ΔΔG = **正**值 ✅ 与「正值 = 削弱结合」一致。

## 2. 与 CATH-ddG 表 `ddG` 列交叉核对

- CATH-ddG 表: **7085 行 × 22 列**；`ddG` 非空 6798 行

- 与本地 SKEMPI 按 `(PDB, Mutation_cleaned)` 连接命中: **6502** 行
- **Pearson(ΔΔG_calc, ddG_CATH) = +nan**；**符号一致比例 = 94.99%**

| #Pdb | Mutation | ΔΔG_calc (本地复算) | ddG (CATH 表) | 符号一致? |
|---|---|---:|---:|---|
| 1CSE | LI38G | +2.2807 | +2.2806 | ✅ |
| 1CSE | LI38S | +1.1888 | +1.1888 | ✅ |
| 1ACB | LI38D | +6.8671 | +6.8667 | ✅ |
| 1IAR | NA89A | +1.5589 | +1.5588 | ✅ |
| 1JTG | EA143A,KB74A | +4.0622 | +4.0619 | ✅ |
| 1JTG | EA143A,FB142A | +2.5812 | +2.5810 | ✅ |

→ ⚠️ 相关仅 +nan，无法确证一致性；保守策略：**按符号一致比例判断**，并在报告中标注不确定性。

## 3. 模型预测列的符号判定（`needs_flip`）

> 判定规则（写死在代码里）：ΔΔG 真值**越小越好**，因此方向正确的模型应与真值**正相关**。若在 CATH 全表（train+test，最大样本量）上`Pearson(pred, ddG) < 0`，则判 `needs_flip = True`。
> **同时给出未反号与反号两套指标**，不隐藏反号带来的人为提升。

### 3a. CATH-ddG 表中的预测列（held-out CATH test 可用）

| 列 | 模型 | n 非空 | Pearson(全表) | 与真值同方向? | needs_flip |
|---|---|---:|---:|---|---|
| `Flex-ddG` | flex ddG | 6569 | +0.4970 | ✅ 同向 | **False** |
| `DDGPred` | DDGPred | 6583 | +0.5768 | ✅ 同向 | **False** |
| `ddG_FoldX` | FoldX | 6706 | +0.4965 | ✅ 同向 | **False** |
| `ESM-IF` | ESM-IF | 6706 | -0.2625 | ❌ 反向 | **True** |
| `proteinMPNN` | ProteinMPNN | 6706 | +0.3137 | ✅ 同向 | **False** |

### 3b. DSSA-PPI RBM 案例集（`EXP_fitness` 为真值）

- 行数 **1702**；列: ['wt', 'pos', 'mut', 'EXP_fitness', 'DSSA-PPI', 'PPIFormer', 'ESM-2-FFN', 'MuTon', 'MpbPPI', 'ProBASS', 'SAAMBE-SEQ', 'FoldX']
- ⚠️ 注意：该文件的真值列名为 `EXP_fitness`（实验**适应度**），**不是 ΔΔG**；其符号含义未在仓库中说明 → 按「预测应与真值同向（正相关）」判定 `needs_flip`，但**该数据集不进入 Hit@10（ΔΔG<0）口径**，只用于 Pearson/Spearman 与 within-position 排序对比。

| 列 | n | Pearson(全表) | needs_flip |
|---|---:|---:|---|
| `DSSA-PPI` | 1702 | +0.3081 | **False** |
| `PPIFormer` | 1702 | -0.1519 | **True** |
| `ESM-2-FFN` | 1702 | -0.0134 | **True** |
| `MuTon` | 1702 | +0.0009 | **False** |
| `MpbPPI` | 1702 | +0.0799 | **False** |
| `ProBASS` | 1702 | +0.0047 | **False** |
| `SAAMBE-SEQ` | 1702 | -0.1934 | **True** |
| `FoldX` | 1702 | -0.5151 | **True** |

## 4. 未能审计的数据集（如实记录）

| 目标 | 状态 | 原因 |
|---|---|---|
| SAMAffinity | `PREDICTION_FILE_NOT_FOUND` | 仓库 `data/test_fold/*/*.csv` 只含 `ddG` 真值与突变标识，**无任何预测值列**（已逐个检查 31 列） |
| Prompt-DDG | `PREDICTION_FILE_NOT_FOUND` | CATH 表中 `prompt-fold0/1/2` 三列只是 **split 标签**（train/val），无数值预测 |
| MINT | `PREDICTION_FILE_NOT_FOUND` | 仓库 68 个文件中仅 `mint/helpers/predict.py`（428 B）与一个 metrics notebook，无预测结果文件 |
| Pythia-PPI | `PREDICTION_FILE_NOT_FOUND` | 仓库 32 个文件，grep predict/result/test 无匹配结果文件 |

