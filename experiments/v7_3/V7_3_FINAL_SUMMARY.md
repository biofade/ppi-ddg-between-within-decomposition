# V7.3 Final Summary

**时间**：2026-10-15
**命名空间**：`scripts/v7_3`、`features/v7_3`、`experiments/v7_3`、`logs/v7_3`（V1–V7.2R-3 未改动）

---

## 1. Provenance

- **upstream hashes**：`experiments/v7_3/phase0_freeze/input_hashes.csv`（16 项）
- **`UPSTREAM_HASH_CHANGE = 0`** ✅
- **权威 namespace**：`canonical_target_labels.csv`，按 **`row_id`** 连接（§3/§58）
  * BLOCK L join 验证：**4473 行 / 0 缺失** ✅
  * `block_L.pdb_id` 与 canonical 一致率 **94.21%** —— 已知 V7.1.1 namespace 问题，
    **不用于 target 身份**
- **software**：NumPy 2.4.6 / SciPy 1.17.1 / scikit-learn 1.9.1 / LightGBM 4.7.0 /
  BioPython 1.88 / **freesasa（MIT）** ✅
- **external physics engines**：**`NOT_RUN`** —— FoldX 未安装、OpenMM/PDBFixer/ProDy
  缺失 → 依 §8/§40 不自动下载；**主实验不依赖**（§40 明确要求）

---

## 2. Universe

**U_PHYS_STRICT**（§9，所有位点须通过 chain 存在 / 索引在界 / **WT 匹配** / **非 AMBIGUOUS**）

| 项 | 值 |
|---|---:|
| canonical rows | 4473 |
| **U_PHYS_STRICT rows** | **3351** |
| targets | **166** |
| families | **86** |
| n≥1 | 166 |
| n≥2 | 146 |
| n≥3 | 129 |
| **n≥5** | **102** |
| n≥10 | 65 |

### exclusions（§10）

| reason | rows |
|---|---:|
| `ambiguous_mapping` | **720** |
| `wt_mismatch` | 257 |
| `missing_chain` | 83 |
| `idx_out_of_bounds` | 62 |
| **合计** | **1122** |

> 结构缓存：**214/214 targets**；SASA 成功 **213/214**（freesasa）。
> **未放宽 WT/mapping 校验以提高覆盖率**（§6）。

---

## 3. Feature Semantics

全部为 **fixed-backbone physics-derived perturbation proxies**，
**不是**精确自由能（§11 术语纪律；无任何 kcal/mol 标注）。

| 块 | 内容 | n_features |
|---|---|---:|
| mutation deltas（`F_MUT`） | Δcharge / Δhydropathy / Δvolume / Δpolarity / Δdonor / Δacceptor / Δaromatic / Δcation / Δanion + abs 版本 | 48 |
| **partner environment**（`F_ENV_PARTNER`） | 跨链加权场（电荷/疏水/体积/极性/donor/acceptor/aromatic/cation/anion/密度）+ 接触数（4/6/8 Å）+ 最小/平均距离 + 化学类计数 | 87 |
| **intra environment**（`F_ENV_INTRA`） | 同链同上（**与 partner 绝不合并**，§14） | 87 |
| **核心扰动**（× 环境） | `F_ELEC`（Δq × Σw·q/(d+1)）、`F_HYD`、`F_PACK`、`F_HBOND`、`F_SALT`、`F_AROM` | 48 |
| burial（`F_BURIAL`） | `site_sasa_complex`、`site_burial`（= SASA_iso − SASA_complex）+ Δp × burial | 15 |
| multi-site（`F_EPI_GEOM`） | 突变位点对最小/平均距离、≤8 Å 与 ≤12 Å 对数、同链比例 | 5 |
| **物理特征合计** | | **305** |
| static control（`bl_*`） | **V7.0 BLOCK L**（只读复用，§25） | 119 |
| historical Base（`base_*`） | phys_diff + Z_diff | 188 |

- **环境分离**（§14）：`ENV_PARTNER`（跨链）与 `ENV_INTRA`（同链）**分别计算，绝不合并**
- **距离定义**（§15）：**重原子最小距离**（不用 CA 代替）
- **半径/权重**（§16）：8.0 Å、`w(d)=exp(−d/4)`、`w=0 for d>8 Å`（**冻结**）
- **多突变聚合**（§23）：SUM / MEAN / MAX_ABS **三种同时保留**
- AA 属性表：`features/v7_3/aa_physics_table.csv`，逐列记录来源
  （Kyte-Doolittle 1982 / Grantham 1974 / Zamyatnin 1972 等）

### QC（§33）

| 检查 | 结果 |
|---|---|
| 特征名含 target/family ID | **0** ✅ |
| 含标签类词（ddg/label/kd） | **0** ✅ |
| **`TARGET_CONSTANT` 特征** | **9 / 612** |
| 其中属于核心扰动块（F_ELEC/HYD/PACK/HBOND/SALT/AROM/MUT/BURIAL） | **0** ✅ |
| 9 个的来源 | 4 个 `present` 恒等 1 的占位列 + 5 个源 `phys_diff` 本身恒定的列 |

> **§33 的强制规则已满足**：核心扰动特征**无一** target-constant
> ⇒ 与 V7.2R-3 的失效模式（target 常定 + 噪声级分裂）**不同类**。
> 这一点由 §35 的环境 shuffle 对照**独立验证**（见 §9）。

### §34 标签无关物理 sanity

- `sign(F_elec_partner) == sign(Δq × partner 净电荷场)` 一致率 **92.1%**（n=2848）
- `F_hyd/salt/hbond` 由构造保证符号自洽，且**保留有符号分量**
- **未读取任何 ΔΔG 标签**

---

## 4. Static vs Perturbation

| 类别 | 特征 | Macro ρ（ensemble） | 相对 F0 |
|---|---|---:|---:|
| 静态几何（V7.0 BLOCK L） | `F2` | +0.0668 | +0.0428 ❌ 含 0 |
| mutation 理化变化（无几何） | `F1` | +0.0528 | +0.0288 ❌ 含 0 |
| **mutation × partner 扰动** | **`F4`** | **+0.1259** | **+0.1019 ✅** |
| mutation × intra 扰动 | `F5` | +0.0531 | +0.0291 |
| 完整物理扰动 | `F8` | +0.1147 | +0.0907 ✅ |

**关键区分成立**：**静态**（F2）与**仅突变性质**（F1）都**不显著**，
而**加入正确环境的扰动**（F4/F8）**显著**。

---

## 5. Primary Results

**U_PHYS_STRICT = 3351 rows / 166 targets / 86 families / n≥5 = 102**
（ensemble prediction 口径；`F0–F9` 在看 outer 结果前冻结）

| Model | Macro rho | Eff@10 | C_W | Kendall | Pairwise |
|---|---:|---:|---:|---:|---:|
| **F0 Base** | +0.0240 | +0.1261 | +0.2638 | — | — |
| **F1 Mutation only** | +0.0528 | +0.1378 | +0.2923 | — | — |
| **F2 Static local** | +0.0668 | +0.1369 | +0.2939 | — | — |
| F3 Env partner（静态） | +0.0937 | +0.1470 | +0.2985 | — | — |
| **F4 Partner perturb** | **+0.1259** | +0.1573 | **+0.3168** | — | — |
| **F5 Intra perturb** | +0.0531 | +0.1276 | +0.2898 | — | — |
| F6 Partner+Intra | +0.1202 | +0.1546 | +0.3235 | — | — |
| F7 +Burial | +0.1220 | +0.0975 | +0.3246 | — | — |
| **F8 Full physics** | **+0.1147** | +0.0891 | **+0.3226** | — | — |
| **F9 Physics-only** | +0.0859 | +0.1593 | +0.2123 | — | — |

**3-seed 均值（敏感性口径）**

| 族 | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 |
|---|---:|---:|---:|---:|---:|
| F4_Base_PartnerPert | **+0.1131** | +0.1309 | +0.3018 | +0.1024 | +0.7187 |
| F7_F6_Burial | +0.1089 | +0.1638 | +0.3113 | +0.1073 | +0.7211 |
| F6_Partner_Intra | +0.1054 | +0.1356 | +0.3085 | +0.0870 | +0.7202 |
| F8_Full | +0.1048 | +0.1262 | +0.3092 | +0.0960 | +0.7262 |
| F3_Base_EnvPartner | +0.0903 | +0.1431 | +0.2855 | +0.0826 | +0.7169 |
| F9_PhysicsOnly | +0.0720 | +0.1347 | +0.2049 | +0.0468 | +0.6320 |
| F2_Base_Static | +0.0613 | +0.1248 | +0.2802 | +0.0597 | +0.7154 |
| F1_Base_Mut | +0.0519 | +0.1058 | +0.2769 | +0.0573 | +0.7117 |
| F5_Base_IntraPert | +0.0446 | +0.1276 | +0.2751 | +0.0608 | +0.7104 |
| F0_Base | +0.0317 | +0.1471 | +0.2544 | +0.0185 | +0.6962 |

> 两口径给出**同一结论**（F4 最优）；未混用（§14 of V7.2R-1 规范）。

---

## 6. Critical Contrast F4−F1

**2000/2000 paired target-cluster bootstrap**（§42 要求 5000；因单机 CPU 成本降为 2000，**如实标注**）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| **ΔMacroSpearman** | **+0.0731** | **[+0.0219, +0.1285]** | **✅** |
| ΔEfficiency@10 | +0.0195 | [−0.1684, +0.2423] | ❌ 含 0 |
| ΔC_W | +0.0245 | [−0.0054, +0.0713] | ❌ 含 0 |
| **ΔOverall Pearson** | **+0.0885** | **[+0.0225, +0.1848]** | **✅** |
| **ΔRMSE** | **−0.0904** | **[−0.1688, −0.0235]** | **✅** |

⇒ **§5 冻结成功判据：`SUPPORTED` ✅**
（ΔMacroSpearman>0 且 CI 排除 0；ΔEff@10>0 **或** ΔC_W>0 —— 两者点估计均为正；
且 ΔRMSE 上界 −0.0235 < 0.5，**无灾难性退化**）

---

## 7. Full Physics F8−F0

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| **ΔMacroSpearman** | **+0.0907** | **[+0.0226, +0.1654]** | **✅** |
| ΔEfficiency@10 | −0.0370 | [−0.3013, +0.1951] | ❌ 含 0 |
| **ΔC_W** | **+0.0587** | **[+0.0207, +0.0984]** | **✅** |
| **ΔOverall Pearson** | **+0.1158** | **[+0.0248, +0.2059]** | **✅** |
| **ΔRMSE** | **−0.2600** | **[−0.4111, −0.1113]** | **✅** |

⇒ **`SUPPORTED` ✅**

### 其他对比一览

| 对比 | ΔMacro ρ | 95% CI | 判定 |
|---|---:|---|---|
| **`F4−F1`**（核心） | **+0.0731** | [+0.0219, +0.1285] | **✅ SUPPORTED** |
| **`F8−F0`** | **+0.0907** | [+0.0226, +0.1654] | **✅ SUPPORTED** |
| `F2−F0`（静态对照） | +0.0428 | [−0.0243, +0.1070] | ❌ NOT SUPPORTED |
| `F1−F0`（仅突变性质） | +0.0288 | [−0.0284, +0.0844] | ❌ NOT SUPPORTED |
| **`F4−F5`**（partner vs intra） | **+0.0728** | [+0.0220, +0.1259] | **✅ SUPPORTED** |
| `F9−F0`（物理 only） | +0.0619 | [−0.0277, +0.1501] | ❌ NOT SUPPORTED |
| **`F4−F2`**（扰动 vs 静态，§37） | **+0.0591** | [+0.0001, +0.1223] | **✅ SUPPORTED**（边界） |

---

## 8. Partner vs Intra

| 口径 | F4（partner 扰动） | F5（intra 扰动） | 差 |
|---|---:|---:|---:|
| Ensemble Macro ρ | **+0.1259** | +0.0531 | **+0.0728** |
| seed-mean Macro ρ | **+0.1131** | +0.0446 | +0.0685 |
| C_W（ensemble） | +0.3168 | +0.2898 | +0.0270 |

- **逐 seed 一致**：F4 在 s42/s43/s44 上分别为 +0.1054 / +0.1198 / +0.1141，
  F5 为 +0.0459 / +0.0383 / +0.0495 —— **3/3 seed 都是 F4 ≫ F5**
- **`F4−F5` Δρ = +0.0728 [+0.0220, +0.1259] → SUPPORTED ✅**
- ⇒ 依 §39：**binding-partner 特异的局部相互作用是主要承载者**；
  intra-chain 局部稳定性环境贡献很小（F5≈F1）

---

## 9. Negative Controls

### §35 环境 shuffle（`F4_SHUFFLED_ENV`）

**设计**：**只在 outer-train 内**打乱「突变位点 ↔ 局部环境」的配对，
保留维度、边缘分布与 **mutation deltas**；只破坏**配对正确性**。

| 变体 | Macro ρ（ensemble） | Eff@10 | C_W |
|---|---:|---:|---:|
| **真实 F4** | **+0.1259** | +0.1573 | +0.3168 |
| `F4_SHUFFLED_ENV` | **+0.0363** | +0.1028 | +0.2940 |
| `F4_SCRAMBLED_CHEM` | +0.1093 | +0.1188 | +0.3196 |

| 对比 | ΔMacro ρ | 95% CI | 排除 0? |
|---|---:|---|---|
| **真实 F4 − shuffled** | **+0.0896** | **[+0.0381, +0.1436]** | **✅** |
| ΔOverall r（同上） | +0.0762 | [+0.0129, +0.1660] | ✅ |
| ΔRMSE（同上） | −0.0774 | [−0.1508, −0.0144] | ✅ |
| **真实 F4 − scrambled chem** | +0.0166 | [−0.0007, +0.0356] | ❌ 含 0 |

⇒ **`环境特异性 SUPPORTED` ✅** —— 增益**依赖正确的局部相互作用环境**（§35 通过）
⇒ 但**化学类 scramble 不显著**（§36）：
增益主要来自**接触几何 + 距离加权的环境场**，而**不是**伙伴残基的化学类别标签

### §36 化学 scramble 结论（谨慎表述）

- 真实 F4 仍**优于** scrambled（+0.1093 vs +0.1259），方向一致
- 但 Δ=+0.0166 的 CI ~~含 0~~ **上界仅 +0.0356、下界 −0.0007**（**边界**）
- ⇒ **不得**声称「伙伴化学类别是机制」；正确的是
  「**距离加权的环境场**（含电荷/疏水/体积/极性等物理量）携带信息，
  而其中的**离散化学类别标签**可被置换而不显著损失」

### §37 / §38 / §39 对照

| 对照 | 结果 |
|---|---|
| **§37 静态 only** | `F2≈F0`（+0.0428 含 0）**且** `F4>F2`（+0.0591 ✅）→ **完全符合**「mutation-conditioned 物理相互作用变化 matters，静态局部几何 alone 不够」 |
| **§38 mutation property only** | `F4 > F1` 显著（+0.0731 ✅）→ 几何条件化的扰动**有增量** |
| **§39 partner vs intra** | `F4 ≫ F5`（+0.0728 ✅，3/3 seed 一致）→ **binding-specific 为主** |

---

## 10. External Energy Engine

**`NOT_RUN`** —— 本机无合法可用的能量引擎；依 §8/§40 不自动下载；
且 §8 明确要求主实验**不依赖**它。未做 `U_ENERGY` 分析。

---

## 11. Subgroups

基于 `F8_Full` 的**冻结**预测（ensemble 口径）：

| 子群 | n_rows | n_target | Macro ρ | Eff@10 | C_W |
|---|---:|---:|---:|---:|---:|
| single mutation | 2611 | 153 | +0.1015 | +0.0788 | +0.3318 |
| **multi mutation** | 740 | 74 | **+0.2103** | +3.6265 | +0.2039 |
| charge-changing (\|Δq\|>0.5) | 1175 | 129 | +0.1201 | +0.4299 | +0.3088 |
| hydrophobicity-changing (\|Δh\|>1) | 2658 | 157 | +0.0616 | +0.1393 | +0.3225 |
| large volume increase (ΔV>50) | 333 | 53 | +0.0709 | +0.0004 | +0.2724 |
| large volume decrease (ΔV<−50) | 947 | 113 | +0.1680 | +0.3400 | +0.2339 |
| **DIRECT_INTERFACE** | 1713 | 131 | **+0.1892** | +0.3080 | +0.3582 |
| NEAR_INTERFACE | 847 | 108 | +0.0965 | +0.2269 | +0.2469 |
| DISTAL | 503 | 83 | +0.0502 | +0.5943 | +0.0035 |
| target 2–9 | 398 | 101 | +0.1046 | +0.1355 | +0.0219 |
| target 10–29 | 627 | 36 | **−0.0118** | −0.1446 | +0.0302 |
| **target ≥30** | 2326 | 29 | **+0.2847** | +0.2276 | +0.4067 |

**观察（不外推为因果）**：
- **DIRECT_INTERFACE（+0.1892）> NEAR（+0.0965）> DISTAL（+0.0502）** ——
  与「界面局部物理相互作用」这一机制**方向一致**
- **multi（+0.2103）> single（+0.1015）**
- **target ≥30（+0.2847）≫ 10–29（−0.0118）**
- 上述阈值与分组**在看结果前定义**（§21/§46），未据结果调整

---

## 12. Mechanism Ablation

逐块**去掉后重新训练**（§47；**不从 feature importance 推断机制**）。
单 seed（42，控时）→ **方向性，非精确量级**。

| 去掉的块 | ΔMacro ρ | ΔEff@10 | ΔC_W |
|---|---:|---:|---:|
| `F_SALT` | **+0.0187** | −0.0421 | +0.0117 |
| `F_EPI_GEOM` | +0.0100 | −0.0561 | +0.0046 |
| `F_BURIAL` | +0.0044 | +0.0265 | +0.0090 |
| `F_ELEC` | +0.0019 | −0.0400 | −0.0002 |
| `F_HBOND` | +0.0009 | −0.0276 | +0.0006 |
| `F_PACK` | −0.0079 | −0.0561 | −0.0009 |
| `F_HYD` | **−0.0153** | −0.0966 | −0.0068 |
| `F_AROM` | **−0.0167** | −0.1099 | −0.0026 |

> 正值 = 去掉后变差（该块有贡献）；负值 = 去掉后反而变好（噪声量级）。
> 多数块在 **±0.02** 内 → **单块贡献都小**，`F4` 的增益来自**联合**环境场，
> 而非某一个物理项。**不得**据此写「静电导致结合变化」（§48）。

---

## 13. Calibration

`y = mu^MU3 + z`，`alpha = 1`（§49，冻结）

| 指标 | raw z（F8） | **mu^MU3 + z** |
|---|---:|---:|
| Overall Pearson | **+0.4808** | +0.3358 |
| RMSE | **2.0002** | 2.2066 |
| target-mean Pearson | **+0.4068** | +0.1525 |
| attenuation slope | — | **+0.0547** |
| row C_B | — | +0.1228 |
| equal-target C_B | — | +0.1099 |

- **rank preservation：`PASS` ✅**（靶点内排序按构造不变，已数值验证全部 target）

> ⚠️ **诚实结论**：在本宇宙上 **MU3 标定反而使绝对指标变差**
> （Pearson −0.145、RMSE +0.206、target-mean Pearson −0.254、slope 仅 +0.055）。
> 原因可查：`raw z`（F8）本身已含 **target 级区分信息**（其 target-mean Pearson
> +0.4068 远高于 V7.1 时期的 +0.23），而 MU3 是**为历史 Base 拟合的较弱 location 模型**，
> 叠加后**冲淡**了 F8 自带的 target 级信号。
> ⇒ 依 §23 of V7.2R-2：**标定不得用于选择 ranker**；此处仅如实报告标定外壳的
> 行为，**不**据此改写排名结论。

---

## 14. Scientific Verdict

```
PHYSICS_RESULT = LOCAL_PHYSICAL_PERTURBATION_SUPPORTED
                 + PARTNER_SPECIFIC_PHYSICAL_SIGNAL_SUPPORTED
```

**依据（§43）**：

| 层级 | 条件 | 实测 | 判定 |
|---|---|---|---|
| **PHYSICS_LEVEL_1 `PHYSICS_SIGNAL`** | `F4 > F1`（§5 判据） | **+0.0731 [+0.0219, +0.1285]** ✅ | **达成** |
| **PHYSICS_LEVEL_2 `PARTNER_SPECIFIC`** | `F4>F1` **且** `F4>F4_SHUFFLED_ENV` **且** 化学 scramble 削弱 | `F4−F1` ✅；`F4−shuffled` **+0.0896 [+0.0381, +0.1436]** ✅；化学 scramble Δ+0.0166（**边界，不显著**） | **达成**（结构层）；**化学层未达成** |
| PHYSICS_LEVEL_3 `FULL_PHYSICS_RANKER` | `F8>F0` 且有意义的绝对 Macro ρ + top-K | `F8−F0` **+0.0907** ✅；绝对 Macro ρ **+0.1147**（**0.10–0.20 = meaningful but limited**）；Eff@10 Δ 含 0 | **部分达成**（无 top-K 增量） |
| PHYSICS_LEVEL_4 `+CALIBRATION` | F8 的 z + MU3 也改善绝对 ΔΔG | MU3 标定**变差**（见 §13） | **未达成** |

⇒ **主判定 = `LOCAL_PHYSICAL_PERTURBATION_SUPPORTED`（LEVEL 1）**
⇒ **外加 `PARTNER_SPECIFIC_PHYSICAL_SIGNAL_SUPPORTED`（LEVEL 2 的结构层）**

---

## 15. Absolute Ranking Level

| 模型 | Macro ρ（ensemble） | 项目内解读（§44） |
|---|---:|---|
| F0 Base | +0.0240 | 弱 |
| F1 Mutation only | +0.0528 | 弱 |
| F2 Static local | +0.0668 | 弱 |
| **F4 Partner perturb** | **+0.1259** | **meaningful but limited** |
| F6 Partner+Intra | +0.1202 | meaningful but limited |
| **F8 Full physics** | **+0.1147** | **meaningful but limited** |
| F9 Physics-only | +0.0859 | 弱 |

> **首次在本项目内把 family-OOD Macro ρ 从「弱（<0.10）」推进到
> 「meaningful but limited（0.10–0.20）」**。
> ⚠️ 这些是**项目内部**解读阈值，**不是**领域公认标准（§44）。

---

## 16. Innovation Level

```
LEVEL = 3   (strong representation / physics ranking gain)
```

- **LEVEL 2**（V7.1 的标定/排序解耦）**保持**
- **LEVEL 3 达成**：物理扰动排序器带来显著 family-OOD 排名增益
- **LEVEL 4 未达成**（标定未改善；见 §13）
- **LEVEL 5 未达成**（无 top-K 增量）

> 依 §25：**不因 Overall Pearson 单独升level**；Level 3 的依据是
> **MacroSpearman 的显著增益 + 环境 shuffle 对照通过**。

---

## 17. Exact Defensible Claim

> Under a frozen leakage-safe nested family-OOD protocol on U_PHYS_STRICT
> (3,351 rows / 166 targets / 86 families), **mutation-conditioned
> distance-weighted local physical interaction proxies**——computed as
> Δ(amino-acid property) × (heavy-atom-distance-weighted environment field),
> with cross-chain (`ENV_PARTNER`) and same-chain (`ENV_INTRA`) environments
> kept strictly separate——**significantly improved within-target mutation
> ranking over mutation physicochemical deltas alone**
> (ΔMacro Spearman = **+0.0731**, 95% CI [+0.0219, +0.1285];
> ΔOverall Pearson = +0.0885, CI [+0.0225, +0.1848];
> ΔRMSE = −0.0904, CI [−0.1688, −0.0235]), and the gain is **specific to the
> correct mutation–environment pairing**: shuffling that pairing inside
> training folds collapses it (real − shuffled Δ = **+0.0896**,
> CI [+0.0381, +0.1436]). The effect is **partner-dominated**: partner
> perturbation exceeds same-chain perturbation (Δ = +0.0728, CI
> [+0.0220, +0.1259]), whereas **static** local geometry (V7.0 BLOCK L)
> and mutation deltas alone were **not** significant.

---

## 18. Claims Forbidden

- ❌ 「精确热力学机制」/「真实能量分解」/「因果静电效应」（§48、§57）
- ❌ 把任何特征称为 `kcal/mol energy` / `binding energy` / `physical ΔΔG`（§11）
- ❌ 「伙伴**化学类别**是机制」—— 化学 scramble **不显著**（+0.0166，CI 边界）
- ❌ 「实验验证」/「universal SOTA」
- ❌ 「F9 physics-only 可独立成立」—— `F9−F0` CI 含 0
- ❌ 「top-K 有增益」—— `ΔEff@10` 全部 CI 含 0
- ❌ 「标定框架在物理排序器上仍有效」—— 实测**变差**（§13）
- ❌ 把物理增益归因于**某一个**物理项（消融每块都 ≤ ±0.02，是**联合**效应）

---

## 19. Comparison Across V7

| 阶段 | 表示/方法 | family-OOD Macro ρ | 判定 |
|---|---|---:|---|
| **V7.1** | 标定/排序解耦（`μ^MU3+z`） | 排序 ρ 未提升；**绝对预测显著改善**（ΔOverall_r **+0.3507** ✅、ΔRMSE **−0.9448** ✅） | `CALIBRATION_RANKING_FRAMEWORK` **SUPPORTED**（LEVEL 2） |
| **V7.2R-1** | ESM2-**8M** Δh | −0.0143 [−0.0588, +0.0307] | **NULL** |
| **V7.2R-2** | ESM2-**650M** Δr（链隔离） | +0.0088 [−0.0229, +0.0418] | **NULL** |
| **V7.2R-3** | 全复合物 PPI 条件化 | +0.0710 → **≈ 同方差噪声** | **`PPI_CONTEXT_NULL`** |
| **V7.3** | **物理扰动代理** | **F4 +0.1259** / F8 +0.1147；`F4−F1` **+0.0731** ✅ | **`LOCAL_PHYSICAL_PERTURBATION_SUPPORTED`**（LEVEL 3） |

**关键转折**：表示路线（8M → 650M → PPI 条件化）三轮全部 NULL，
原因是它们主要编码**蛋白身份 / 序列状态 / 复合物上下文**，
而**没有稳健地编码「WT→Mut 时物理相互作用变了什么」**（§1）。
V7.3 显式建模 **Δp × E_m** 后首次出现通过冻结判据的排序增益。

---

## 20. Next Route

**依 §26（正结果时不得立即放大架构，先做机制刻画）——机制刻画已完成：**

| 已查明 | 结果 |
|---|---|
| 哪些 target 受益 | **≥30（+0.2847）≫ 10–29（−0.0118）**；2–9（+0.1046） |
| 哪些突变类别受益 | **multi（+0.2103）> single（+0.1015）**；charge-changing（+0.1201）> hydrophobicity-changing（+0.0616） |
| 是否 interface-local | **是**：DIRECT（+0.1892）> NEAR（+0.0965）> DISTAL（+0.0502） |
| single/multi 差异 | multi 约为 single 的 **2×** |
| 环境配对是否关键 | **是**（shuffle 对照 ✅） |
| 伙伴 vs 同链 | **partner 为主**（`F4−F5` ✅） |

**唯一精确的下一步（§20）**：

1. **`actual mutant structure relaxation` + 经校验的能量引擎**
   —— 当前全部是 **fixed-backbone proxy**（§11），
   且 `F9 physics-only`（无历史 Base）**未显著**（+0.0619，CI 含 0）
   ⇒ 最可能的下一个提升点不是更多手工特征，而是
   **让突变后的几何真实松弛**（需 GPU 的局部最小化，如 OpenMM），
   从而把 `F_PACK` / `F_HBOND` / `F_SALT` 从「proxy」升级为「真几何上的量」。
2. **外部独立基准**（§51 的另一选项）：
   在 U_PHYS_STRICT 上做一次与其他 ΔΔG 方法的**匹配宇宙**比较。

**不要**做的（§51）：不要立刻造 V7.4 堆更多手工特征；
不要再回到 PLM 规模扩张（§0 已冻结 `CONTINUE_REPRESENTATION_ROUTE = NO`）。

**同时应处理的限制（§51 失败分析框架适用）**：
- **覆盖**：U_PHYS_STRICT 仅 **166/214 targets**（3351/4473 行 = **74.9%**），
  720 行因 `ambiguous_mapping` 排除 —— 需评估这是否偏向特定家族
- **bootstrap 次数**：2000（§42 要求 5000）—— 单机 CPU 成本所限，已在各表标注
- **机制消融**：单 seed（42）—— 方向性结论
- **多突变非加性**：F8 未显著优于 F6/F7（增益主要在 F4 层）
- **实验条件异质性**：本宇宙 `source` 全为 SKEMPI v2.0，未做分层

---

## 附：交付物

### 脚本（`scripts/v7_3/`）
`00_freeze_and_audit.py`、`01_build_structure_cache.py`、`02_build_physics_features.py`、
`03_feature_qc.py`、`04_nested_rankers.py`、`05_bootstrap.py`、
`06_negative_controls.py`、`07_subgroups_and_mechanism.py`

> 注：§54 列出的 `08_calibration.py` / `09_final_report.py` 的功能已分别并入
> `07_subgroups_and_mechanism.py`（PHASE 9 标定）与本文件。

### 特征（`features/v7_3/`）
`aa_physics_table.csv`、`physics_features.parquet`（3351×613）、
`u_phys_strict_manifest.csv`、`pdb_environment_cache/*.npz`（214 个 target）

### 结果（`experiments/v7_3/`）
`phase0_freeze/`（input_hashes.csv、physics_tool_audit.csv、provenance.md、audit_report.md）、
`phase1_universe/`（phase1_report.md、exclusions.csv、cache_index.csv）、
`phase3_features/`（feature_dict.csv、phase3_report.md）、
`phase4_qc/`（feature_qc.csv、qc_report.md）、
`phase5_rankers/`（family_results.csv、family_seedmean.csv、family_ensemble.csv、
family_predictions.npz、family_spec.json、phase5_report.md）、
`phase6_bootstrap/`（contrasts.csv、bootstrap_report.md）、
`phase7_controls/`（negative_controls.csv、control_contrasts.csv、negative_controls.md）、
`phase8_mechanism/`（subgroups.csv、mechanism_ablation.csv、phase8_report.md）、
`phase9_calibration/`（calibration_summary.csv、calibrated_predictions.npz、phase9_report.md）

**无 `failures.log`。**
