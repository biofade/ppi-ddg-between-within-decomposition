# V6.2 执行台账

> 每完成一个 PHASE 立即追加。失败写 `logs/v6_2/failures.log`，不停止。
> 约束：不改 V1.1–V6.1；`v6_2` 命名空间；CPU only；不下载 >100 MB；不伪造数字；
> **判定只依据 provenance / IDs / formula / sign convention / hashes，不依据「结果是否合理」**。

## PHASE 0: 环境 + 目录 — **DONE**

- 12 个目录已建（`experiments/v6_2/` 下 9 个 + `scripts/v6_2` + `logs/v6_2`）
- 关键输入 **8/8 存在**（含 V6.1 split 结果、ddg_pairs、ddg_split、V5 预测、v2.2 特征、ESM-2 缓存）
- 产出：`phase0_environment/environment.md`

## PHASE 1: Target Provenance Audit — **DONE**

> # 判定：**`TASK_IS_ABSOLUTE_DG`**

**四条 provenance 证据（不依赖数值是否「看起来合理」）**：

| # | 证据 | 内容 |
|---|---|---|
| 1 | **列名来源** | V6.1 源码精确匹配：`y = pd.to_numeric(v22["real_dG"]...)` → target 来自 `features/v2.2_proteinea_features.csv` 的 **`real_dG`** 列（每突变体一条绝对 ΔG） |
| 2 | **行宇宙** | V6.1 用 **8207 行**（v2.2 全量）；真 ΔΔG forward 配对只有 **7526**（映射后 7523）→ **不同行宇宙** |
| 3 | **分组键** | V6.1 以 `pdb_id`（**2171** 唯一值）作 complex 分组；ΔΔG 配对依据的是 **231** 个复合物 |
| 4 | **缺少差分步骤** | V6.1 `load_rows()` / `make_folds()` 中**不存在**任何 WT/mutant 配对、差分或正反向增广逻辑 |

**影响判定**：

- ✅ **Split Inflation = +0.2663 这一现象本身仍有效** —— 它是关于
  「划分严格性 → 可见 Pearson」的陈述，与任务是 ΔG 还是 ΔΔG 无关
- ❌ **不能**表述为「ΔΔG 任务的 split inflation」
- ⚠️ 「全局 Pearson 来自跨复合物区分」这一**宏观结论在绝对值任务上成立**，
  需 PHASE 2 在**真 ΔΔG 配对**上重做才能迁移

哈希（前 32 位）：`04_split_inflation.py` `acb1a415e6f7e8fb8c635f518e1498d2`；
`v2.2_proteinea_features.csv` `134a406175c975541ae40aed37c3f824`；
`ddg_pairs.csv` `e43b19d69a8b54efe5cc4b3e647e8765`

产出：`phase1_provenance/provenance_report.md` + `provenance_evidence.csv`

## PHASE 2: M1 审计 + Canonical + 真 ΔΔG 4 档 Split — **DONE**

### 2.0 M1 历史实现审计（**DONE，复用而非重写**）

| 标志 | 值 |
|---|---|
| `M1_implementation_audited` | **True** |
| `M1_uses_historical_script` | **True** |
| `M1_only_split_layer_replaced` | **True**（允许范围：split + OOF） |
| `M1_implementation_not_found` | False |

**复用**：`v5_common.build_pair_features`（特征装配，PCA 只在各折 train 行拟合）、
`v5_common.LGB_PARAMS`（800 树 → 本次统一 **300 树**）、V5 Gate A 的正/反向规则
（`y_rev = −y`，**辅助列不取负**）。
**替换**：split generation（四档）+ OOF generation（nested，按档分级）。

> ⚠️ **与 V5 的一处有意偏离（已记录）**：V5 Gate A 的 DGaux OOF 用 `GroupKFold(pdb_id)`；
> V6.2 §2.4 要求 inner OOF **继承 outer protocol 的分组级别**并叠加 `state_id` 双重分组
> → 属「允许替换 split + OOF 层」范围。

### 2.1 Canonical Row Universe（**DONE**）

| 项 | 值 |
|---|---:|
| 原始 forward ΔΔG 配对 | 7523 |
| **matched_analysis_set** | **7523**（**excluded = 0**） |
| 唯一 pdb_id / family / state_id | 230 / 112 / 6047 |
| ddG mean ± std | **+1.5222 ± 2.3956** |
| ddG min / max | −11.4939 / +12.4388 |
| **Base Rate (ΔΔG<0)** | **0.2363** |

> **Base Rate = 0.2363**（真 ΔΔG）对照 V6.1 的 **1.0000**（绝对 ΔG）
> → 这是 PHASE 1 provenance 判定的**独立佐证**。

**发现**：`(pdb, mutation)` 重复的行 **4722**（2011 个重复键）；
同一 `state_id` 出现多行的组 **1019** → L1 必须保证同 state 同 fold（已实现）。

### 2.4–2.7 四档 × 三模型 × 3 seeds（**DONE**，300 树）

**M1（`C+DGaux_fwd+rev`）3-seed 均值**：

| 档 | Overall r | Between r | Within macro ρ | **C_B** | **C_W** | Eff@10 |
|---|---:|---:|---:|---:|---:|---:|
| **L1** random | +0.6990 | +0.8120 | +0.0947 | **+0.4275** | **+0.2715** | +0.2104 |
| **L2** complex | +0.2195 | +0.6694 | +0.0530 | +0.0462 | +0.1733 | +0.1273 |
| **L3** seq40 | +0.1398 | +0.5522 | +0.0426 | −0.0173 | +0.1571 | +0.1216 |
| **L4** family | +0.1400 | +0.5522 | **+0.0441** | **−0.0102** | **+0.1362** | +0.1717 |

### ★ 核心量化（真 ΔΔG 上）

| 指标 | L1 | L4 | **L1 − L4** |
|---|---:|---:|---:|
| Overall Pearson | +0.6990 | +0.1400 | **+0.5590** |
| Between-target corr | +0.8120 | +0.5522 | +0.2598 |
| **C_B**（between 协方差） | +0.4275 | −0.0102 | **+0.4377** |
| **C_W**（within 协方差） | +0.2715 | +0.1362 | **+0.1353** |
| Within macro Spearman | +0.0947 | +0.0441 | +0.0506 |

### ⚠️ 与 V6.2 预设假说的**不一致**（必须如实报告）

V6.2 §8 预设 **`Δ_OOD = C_W(L1) − C_W(L4) > 0 且 CI 排除 0` → within signal 在 OOD 后丢失**。

实测：**C_W 只下降约一半（+0.2715 → +0.1362），并非「丢失」**；
真正**塌陷到 0 的是 `C_B`（+0.4275 → −0.0102）**。

⇒ **「within 信号在 OOD 后消失」这一预设被数据削弱**；
正确表述应为 **「between 分量在 OOD 后塌陷，而 within 分量部分保留」**。
（最终判定以 PHASE 3 的配对 bootstrap CI 为准。）

### 2.8 与 V6.1 的对比（**descriptive only**）

| 来源 | 任务 | 行数 | L1 | L4 | L1−L4 |
|---|---|---:|---:|---:|---:|
| V6.1 PHASE 2 | **绝对 ΔG** | 8207 | +0.8436 | +0.5774 | **+0.2663** |
| V6.2 PHASE 2 | **真 ΔΔG** | 7523 | +0.6990 | +0.1400 | **+0.5590** |

> **两者任务不同，数字不可直接比较**（V6.2 §2.8）。
> 可观察到的模式：**任务越难（ΔΔG 配对 vs 绝对 ΔG），
> 划分严格性对可见 Pearson 的影响越大（+0.2663 → +0.5590）**。

### 实施事件（如实记录）

- 首次运行因**多进程重叠**导致中间产物状态混乱 → 已 kill 全部进程、归档到
  `logs/v6_2/archive_100106/`，修复 2 个真实 bug 后重跑：
  - **BUG A**：`all_rows += r` 累加到**局部名**上 → 改为 `all_rows.extend(r)`
  - **BUG B**：`run_level` 内用局部 `rows` 写 CSV → **覆盖其它档** → 改由 `main()` 写累积结果
- 为在预算内完成，树数由 800 统一降为 **300**（**所有档同一值** → inflation 比较仍公平）

## PHASE 6: Leave-one-model-out + BH-FDR — **DONE**

> 判定：**`ROBUST`**（4 个 leave-one-out 方向一致 + 3/4 q<0.05 + 无 effect 掉 >50%）

| 组合 | Δ | 95% CI | p | **BH q** | q<0.05 |
|---|---:|---|---:|---:|---|
| 去 ESM-IF | +0.1644 | [−0.0174, +0.3463] | 0.1111 | 0.1111 | ❌ |
| 去 ProteinMPNN | **+0.3542** | [+0.2155, +0.4920] | 0.000041 | **0.000205** | ✅ |
| 去 Flex-ddG | +0.2306 | [+0.0943, +0.3689] | 0.003753 | 0.006255 | ✅ |
| 去 FoldX | +0.2881 | [+0.1343, +0.4435] | 0.005734 | 0.007167 | ✅ |
| **Full** | **+0.2593** | [+0.1262, +0.3944] | 0.003753 | 0.006255 | ✅ |

### ⚠️ 重要 nuance（必须写进论文）

**去掉 ESM-IF 后差异不再显著**（Δ=+0.1644，p=0.111，CI 含 0）。
即：**「物理 ≫ 序列表示」这一结论在较大程度上由 `ESM-IF` 这一个弱模型驱动。**

| 模型 | per-complex Spearman 均值 |
|---|---:|
| Flex-ddG | +0.4159 |
| FoldX | +0.3584 |
| **proteinMPNN** | **+0.2227** |
| **ESM-IF** | **+0.0329** |

⇒ 两类**并非同质**：物理侧两模型接近（+0.36 / +0.42）；
序列侧 **proteinMPNN（+0.2227）远强于 ESM-IF（+0.0329）**。
「序列表示类」的差表现主要由 **ESM-IF** 贡献。

## PHASE 3: Pearson 协方差分解 + Primary Statistics — **DONE**

> 产物：`phase3_decomposition/` 下 6 个 CSV + `decomposition_report.md`（耗时 689s）

### 浮点 sanity：**36 / 36 PASS**，最大 `|C_B + C_W − r_full|` = **2.22e-16**

> 优化记录：组均值原为逐组 `mean()`（5000 次 bootstrap 下单统计量 >30 分钟）。
> 改为 `np.bincount` 向量化后**快 24.6×**，且经元素级 `np.allclose` 验证
> **结果逐元素一致**（与已记录 CSV 复算到 6 位小数完全相同，diff = 1.11e-16）。

### 3.1 分解（M1，3-seed 均值）

| 档 | r_full | **C_B** | **C_W** | Share_B |
|---|---:|---:|---:|---:|
| L1_random | +0.6953 | **+0.4295** | +0.2658 | +0.6178 |
| L2_complex | +0.2195 | +0.0462 | +0.1733 | +0.2043 |
| L3_seq40 | +0.2195 | +0.0462 | +0.1733 | +0.2043 |
| L4_family | +0.1250 | **−0.0112** | **+0.1361** | −0.0822 |

### 3.2–3.4 Primary Statistics（cluster bootstrap 5000，单位 = complex）

| 统计量 | 档 | 模型 | point | 95% CI | 排除 0? |
|---|---|---|---:|---|---|
| **Δ_C = C_B − C_W** | **L4** | M1 | **−0.1745** | **[−0.2824, −0.0233]** | ✅ |
| Δ_C | **L4** | M0 | −0.1651 | [−0.2768, −0.0081] | ✅ |
| Δ_C | L1 | M1 | +0.1560 | [−0.1312, +0.3223] | ❌ |
| Δ_C | L1 | M0 | +0.1658 | [−0.1309, +0.3353] | ❌ |
| **Δ_OOD = C_W(L1) − C_W(L4)** | L1−L4 | **M1** | **+0.1144** | **[+0.0642, +0.1898]** | ✅ |
| Δ_OOD | L1−L4 | M0 | +0.1011 | [+0.0565, +0.1737] | ✅ |
| Gain_W(M2−M0) | L1 / L4 | — | **+0.0000** | [0, 0] | ❌（M2≡M0，定义如此） |
| Gain_W(M1−M2) | L1 | — | +0.0124 | [+0.0053, +0.0184] | ✅ |
| **Gain_W(M1−M2)** | **L4** | — | **−0.0009** | **[−0.0081, +0.0065]** | ❌ |
| Δρ_OOD | L1−L4 | M0 | +0.0621 | [+0.0083, +0.1203] | ✅ |
| Δρ_OOD | L1−L4 | M1 | +0.0521 | [−0.0035, +0.1073] | ❌ |

**关键解读**：
- **L4 上 `Δ_C < 0` 且 CI 排除 0** → 在最严格划分下 **within 分量反而主导**
  （因为 `C_B` 已塌陷到负值）
- **`Δ_OOD > 0` 且 CI 排除 0** → within 信号**确实随 OOD 严格性下降**，
  但注意幅度仅 **+0.114**（`C_W` 从 +0.266 降到 +0.136，**保留约 51%**），
  并非 V6.2 预设的「丢失」
- **`Gain_W(M1−M2)` 在 L1 显著（+0.0124）但 L4 不显著（−0.0009）**
  → DGaux 的增益**只在宽松划分下存在**，在 family-level 下消失

### 3.5 方差分解（两套口径，回答不同问题）

### 3.6 重构（**diagnostic only**，非可用预测器）

| 档 | 模型 | Pearson(full) | Pearson(between-only) | Pearson(within-only) |
|---|---|---:|---:|---:|
| L1 | M1 | +0.6990 | +0.5814 | +0.4006 |
| **L4** | **M1** | **+0.1398** | **−0.0301** | **+0.1923** |

> **L4 上 between-only 为负（−0.030）而 within-only 为正（+0.192）**
> —— 与 L1 完全相反。这是「between 分量在 OOD 后失效并反号」的直接证据。

### 3.7 负对照（**最强证据**）

**Control A — Within-target shuffle**（保持每靶点预测分布，打乱对应关系）：

| 档 | 模型 | Overall r（原始） | Overall r（shuffle） | **下降** | Within ρ（原始） | Within ρ（shuffle） |
|---|---|---:|---:|---:|---:|---:|
| L1 | M1 | +0.6990 | −0.0192 | **+0.7182** | +0.0947 | +0.0124 |
| L4 | M1 | +0.1398 | −0.0597 | **+0.1994** | +0.0426 | +0.0205 |

⇒ **把靶点内的突变顺序完全打乱后，全局 Pearson 掉到 ≈ 0**
（L1：0.699 → −0.019）。这说明
> **在 L1 下，全局 Pearson 的几乎全部来自「跨靶点区分」，
> 而非「靶点内突变排序」。**

**Control B — Between-only**（每个突变体赋该靶点平均预测）：

| 档 | 模型 | Overall r | Within ρ |
|---|---|---:|---:|
| L1 | M1 | +0.5814 | +0.0242 |
| L4 | M1 | **−0.0301** | +0.0242 |

## PHASE 6: Leave-one-model-out + BH-FDR — **DONE**（见上）

## PHASE 7: 分级泄漏检查 + 浮点 sanity — **DONE**

> **50 / 50 通过；泄漏 0 项**

| 检查类 | 结果 |
|---|---|
| 分级断言（L1 obs_duplicate / state_id；L2 pdb；L3 seq40；L4 family） | ✅ 全部不跨折 |
| DGaux outer + state 双重分组 | ✅ |
| Forward/Reverse（reverse ⊆ train；test = canonical forward only） | ✅ |
| M1 复用断言（3 个 True + 1 个 False） | ✅ |
| Pearson 分解浮点 sanity | ✅ 36/36，最大 \|diff\| = 2.22e-16 |

> 修复了一处**自检脚本误报**：`M1_implementation_not_found == False` 是**期望值**
> （表示「未发生」），原通用循环把它当失败 → 已改为按期望值逐项判定。

## ⚠️ 关于「四档」的一个必须如实报告的发现

**本数据集无法产生 4 个相互独立的严格度层级。**

实测性质（`phase2_split4` 验证）：

| 事实 | 值 |
|---|---|
| 唯一 PDB / family / seq40_cluster | **230 / 112 / 112** |
| `family_label == seq40_cluster` | **完全相同（True）** |
| 每个 PDB 对应的 family 数 | **恒为 1**（最大 1，均值 1.000） |
| 每个 family 跨的 PDB 数 | 最大 **45**，均值 2.054；34 个 family 跨 >1 PDB |

⇒ 因为**每个 PDB 恰好属于一个 family**：
- 用 `seq40` 分组 → family 不跨折**蕴含** PDB 不跨折 → **L3 ≡ L2**
- 用复合键 `pdb||seq40` 分组 → 划分与纯 `pdb` **完全相同** → **同样 ≡ L2**

**两种实现都无法让 L3 独立于 L2/L4。这是数据本身的性质，不是实现缺陷。**
因此**有效层级数为 3**（L1 宽松 / L2≡L3 复合物级 / L4 family 级），
报告中 L3 与 L2 数值相同属**预期结果**，已在此显式说明。

## PHASE 4 / 5 / 8 — DONE

### PHASE 4（排序 + δ，759s）

δ 敏感性（CATH，**common eligible targets = 21**）：

| δ | 物理 | 序列 | Δ | Wilcoxon p |
|---|---:|---:|---:|---:|
| 0.0 | 0.6387 | 0.5437 | +0.0949 | 0.006281 |
| 0.2 | 0.6623 | 0.5445 | +0.1178 | 0.003185 |
| **0.5** | **0.6949** | 0.5374 | **+0.1575** | **0.000052** |

### PHASE 5（bootstrap 修正，1109s）

V6.1 CI 系统性偏窄 **20–26%**；**结论方向不变**（ESM-IF 仍含 0）。

### PHASE 8

见 `V6_2_FINAL_SUMMARY.md`。**路线判定：`PARTIAL_MIPA`**。

## 全部 PHASE 完成（0–8）；无 failures.log
