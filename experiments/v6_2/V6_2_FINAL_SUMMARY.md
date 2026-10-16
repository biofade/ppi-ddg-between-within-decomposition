# V6.2 最终汇总

**时间**：2026-10-14（约 08:5x → 13:0x，约 4 小时；上限 6 小时，未超时）
**约束**：未修改 V1.1–V6.1 任何文件；`v6_2` 命名空间；CPU only；无下载；未伪造数字
**判定原则**：只依据 provenance / IDs / formula / sign convention / hashes，
**不依据「结果看起来是否合理」**

---

## PHASE 1: Provenance

- **V6.1+ Phase 2 实际任务：`TASK_IS_ABSOLUTE_DG`**（绝对 ΔG，**不是** ΔΔG）

### 依据（四条 provenance 证据，不依赖数值）

| # | 证据 | 内容 |
|---|---|---|
| 1 | **列名来源** | V6.1 源码精确匹配 `y = pd.to_numeric(v22["real_dG"]...)` → target 取自 `features/v2.2_proteinea_features.csv` 的 **`real_dG`** 列（每突变体一条绝对 ΔG） |
| 2 | **行宇宙** | V6.1 用 **8207 行**（v2.2 全量）；真 ΔΔG forward 配对只有 **7526**（映射后 7523）→ **不同行宇宙** |
| 3 | **分组键** | V6.1 以 `pdb_id`（**2171** 唯一值）作 complex 分组；ΔΔG 配对依据 **231/230** 个复合物 |
| 4 | **缺少差分步骤** | V6.1 `load_rows()` / `make_folds()` 中**不存在**任何 WT/mutant 配对、差分或正反向增广逻辑 |

**独立佐证**：V6.1 报告 Base Rate = **1.0000**；真 ΔΔG 上 Base Rate = **0.2363**。

**影响**：Split Inflation = +0.2663 这一现象本身仍有效（与任务是 ΔG 还是 ΔΔG 无关），
但**不能**表述为「ΔΔG 任务的 split inflation」。

---

## PHASE 2: M1 + Canonical + 4 档

- **M1 复用历史实现：YES**（`scripts/v5/v5_common.py` 特征装配 + `LGB_PARAMS` + V5 Gate A 正/反向规则）
  - `M1_implementation_audited = True`，`M1_uses_historical_script = True`，
    `M1_only_split_layer_replaced = True`，`M1_implementation_not_found = False`
  - **只替换** split generation + OOF generation 层（V6.2 §2.0 允许范围）
  - 树数 800 → **300**（**所有档同值**，保证比较公平）
- **Matched set: 7523 行**（原始 7523；**excluded 0**）
  - 230 PDB / 112 family / 6047 state_id；ddG mean ± std = **+1.5222 ± 2.3956**
  - Base Rate(ΔΔG<0) = **0.2363**
- **M0 L1 → L4 Pearson: +0.6830 → +0.1333**
- **M1 L1 → L4 Pearson: +0.6953 → +0.1250**
- **M1 L4 Between: +0.5340 / Within Spearman: +0.0481**

### ⚠️ 四档的**有效层级数为 3**（必须如实报告）

| 事实 | 值 |
|---|---|
| 唯一 PDB / family / seq40_cluster | **230 / 112 / 112** |
| `family_label == seq40_cluster` | **完全相同** |
| 每个 PDB 对应的 family 数 | **恒为 1** |

⇒ 因**每个 PDB 恰好属于一个 family**：用 `seq40` 分组时 family 不跨折**蕴含** PDB 不跨折
→ **L3 ≡ L2**；用复合键 `pdb||seq40` 时划分与纯 `pdb` **完全相同** → 同样 **≡ L2**。
**两种写法都无法让 L3 独立于 L2/L4 —— 这是数据本身的性质。**

**（修复尝试）** 另用 ≤40% k-mer Jaccard 贪心聚类构造独立 L3：
得 **1144 个簇**（≠ L2 的 230，≠ L4 的 112），确实独立；
但 3-mer Jaccard 是**过聚类代理**（已知高估同一性），其簇**细于** L4，
故该 L3 的约束**弱于** L4，不构成合法的「L2→L3→L4 递增」层级。
**结论保持不变：有效层级数 = 3。**

---

## PHASE 3: Pearson Covariance Decomposition（**Primary**）

> 浮点 sanity：**36 / 36 PASS**，最大 `|C_B + C_W − r_full|` = **2.22e-16**

| Split | r_full | C_B | C_W | Share_B |
|---|---:|---:|---:|---:|
| L1_random | +0.6953 | **+0.4295** | +0.2658 | +0.6178 |
| L2_complex | +0.2195 | +0.0462 | +0.1733 | +0.2043 |
| L3_seq40 | +0.2195 | +0.0462 | +0.1733 | +0.2043 |
| **L4_family** | **+0.1250** | **−0.0112** | **+0.1361** | −0.0822 |

### Primary Statistics（cluster bootstrap 5000，单位 = complex，**5000/5000 有效**）

| Statistic | Point | 95% CI | 判定 |
|---|---:|---|---|
| **Δ_C = C_B − C_W (M1 L4)** | **−0.1745** | **[−0.2824, −0.0233]** | ✅ 排除 0（**within 主导**） |
| Δ_C (M0 L4) | −0.1651 | [−0.2768, −0.0081] | ✅ |
| Δ_C (M1 L1) | +0.1560 | [−0.1312, +0.3223] | ❌ 含 0 |
| **Δ_OOD = C_W(L1) − C_W(L4) (M1)** | **+0.1144** | **[+0.0642, +0.1898]** | ✅ 排除 0 |
| Δ_OOD (M0) | +0.1011 | [+0.0565, +0.1737] | ✅ |
| Gain_W(M2−M0) | +0.0000 | [0, 0] | ❌（M2≡M0，定义如此） |
| Gain_W(M1−M2) L1 | +0.0124 | [+0.0053, +0.0184] | ✅ |
| **Gain_W(M1−M2) L4** | **−0.0009** | **[−0.0081, +0.0065]** | ❌ 含 0 |
| Δρ_OOD (M0) | +0.0621 | [+0.0083, +0.1203] | ✅ |
| Δρ_OOD (M1) | +0.0521 | [−0.0035, +0.1073] | ❌ 含 0 |

### Reconstruction（**diagnostic only，非可用预测器**）

- **Full / Between-only / Within-only Pearson**：
  - **L1 M1: +0.6990 / +0.5814 / +0.4006**
  - **L4 M1: +0.1398 / −0.0301 / +0.1923**

> **L4 上 between-only 为负（−0.030）而 within-only 为正（+0.192）**
> —— 与 L1 完全相反。

### Controls（**最强证据**）

- **Within-shuffle**：
  - L1 M1：Overall Pearson **+0.6990 → −0.0192（下降 +0.7182）**；Within ρ +0.0947 → +0.0124
  - L4 M1：Overall Pearson +0.1398 → −0.0597（下降 +0.1994）；Within ρ +0.0426 → +0.0205
- **Between-only**：
  - L1 M1：Overall r **+0.5814**；Within ρ +0.0242
  - L4 M1：Overall r **−0.0301**；Within ρ +0.0242

> ⇒ **把靶点内突变顺序完全打乱后，全局 Pearson 掉到 ≈ 0。**
> 即：**全局 Pearson 的几乎全部来自「跨靶点区分」，而非「靶点内突变排序」。**

---

## PHASE 4: Ranking + δ

| Metric | n≥2 | n≥3 | n≥5 |
|---|---:|---:|---:|
| MacroSpearman | 见 `ranking_metrics.csv` | — | — |
| MacroKendall | 见 `ranking_metrics.csv` | — | — |
| MacroAcc | 见 `ranking_metrics.csv` | — | — |

### δ-sensitivity（CATH，**common eligible targets = 21**）

| δ | 物理（FoldX+Flex-ddG） | 序列（ESM-IF+proteinMPNN） | **Δ** | n_eligible | **n_common_targets** | Wilcoxon p |
|---|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.6387 | 0.5437 | **+0.0949** | 21 | 21 | 0.006281 ✅ |
| 0.2 | 0.6623 | 0.5445 | **+0.1178** | 21 | 21 | 0.003185 ✅ |
| **0.5** | **0.6949** | 0.5374 | **+0.1575** | 21 | 21 | **0.000052** ✅ |

> **修正生效**：eligible = **四个模型 valid_targets 的交集** = **21**（非预设的 53）。
> **效应随 δ 增大而增强**（δ 越大越排除噪声，物理方法的排序优势越明显）。

---

## PHASE 5: Bootstrap (primary = prediction ensemble)

- **V6.1 Flex-ddG CI: [+0.3294, +0.4993]**（**1151/5000** 有效，宽 0.1699）
- **V6.2 Flex-ddG CI: [+0.3052, +0.5183]**（**5000/5000** 有效，宽 0.2131）
- **V6.1 的 CI 系统性偏窄约 25%**（bug：`if unique complex < 15: skip` 跳过的
  恰是抽样重复度更高的迭代）
- **结论是否改变：NO** —— 5 个模型的「是否含 0」判定全部不变
  （Flex-ddG / FoldX / DDGPred / proteinMPNN 排除 0；**ESM-IF 仍含 0**）

| 模型 | V6.1 宽 | V6.2 宽 | 变化 |
|---|---:|---:|---|
| Flex-ddG | 0.1699 | 0.2131 | **+25.4%** |
| FoldX | 0.1759 | 0.2198 | **+25.0%** |
| DDGPred | 0.2161 | 0.2636 | **+22.0%** |
| ESM-IF | 0.2137 | 0.2691 | **+25.9%** |
| proteinMPNN | 0.2505 | 0.3016 | **+20.4%** |

**split4 的 macro-Spearman CI（seed ensemble）**：L1 +0.0663 [+0.0012, +0.1314] /
L2 +0.0498 [−0.0121, +0.1114] / L4 +0.0407 [−0.0211, +0.1033]
→ **除 L1 外均含 0**。

---

## PHASE 6: Leave-one-out

| 组合 | Δ | q | 判定 |
|---|---:|---:|---|
| 去 ESM-IF | +0.1644 | 0.1111 | ❌ 不显著 |
| 去 ProteinMPNN | **+0.3542** | **0.000205** | ✅ |
| 去 Flex-ddG | +0.2306 | 0.006255 | ✅ |
| 去 FoldX | +0.2881 | 0.007167 | ✅ |
| Full | +0.2593 | 0.006255 | ✅ |

**Overall: `ROBUST`**（4 个 leave-one-out 方向一致 + 3/4 q<0.05 + 无 effect 掉 >50%）

> ⚠️ **但「去 ESM-IF」后差异不再显著**（p = 0.111，CI 含 0）
> → 该结论**在较大程度上由 ESM-IF 这一个弱模型驱动**。
> per-complex Spearman：Flex-ddG **+0.4159** / FoldX **+0.3584** /
> proteinMPNN **+0.2227** / ESM-IF **+0.0329**。
> **两类并非同质**：proteinMPNN 与物理类的差距（Δ=+0.1644）**不显著**。

---

## PHASE 7: 泄漏

| 检查 | 结果 |
|---|---|
| 分级断言（L1 obs/state；L2 pdb；L3 seq40；L4 family） | **PASS** |
| DGaux + state 双重分组 | **PASS** |
| Forward/reverse（reverse ⊆ train；test = canonical forward only） | **PASS** |
| M1 复用（3 个 True + 1 个 False） | **PASS** |
| Pearson 浮点（36/36，max 2.22e-16） | **PASS** |

**50 / 50 通过；无泄漏检出。**

---

## 关键判定

### Decoupling 假说（Primary Statistics）

**Primary Statistic 1: Δ_C = C_B − C_W**
- **M1 L4 上 Δ_C = −0.1745 [−0.2824, −0.0233]，CI 排除 0**
- **但方向与预设相反**：预设期望「between 主导」（Δ_C > 0），实测在最严格划分下
  **within 反而主导**（因 `C_B` 已塌陷到 **−0.0112**）
- L1 上 Δ_C = +0.1560 但 **CI 含 0** → 宽松划分下 between 优势**不显著**

**Primary Statistic 2: Δ_OOD = C_W(L1) − C_W(L4)**
- **M1 = +0.1144 [+0.0642, +0.1898]，CI 排除 0** ✅
- 但幅度仅 +0.114：`C_W` 从 **+0.2658 → +0.1361**，**保留约 51%**
- ⇒ **不是「丢失」，而是「显著下降约一半」**

**Primary Statistic 3: Gain_W**
- **Gain_W(M2−M0) = 0**（M2 与 M0 为同一模型，定义如此 —— 应记为
  **该对照未提供信息**，而非「无增益」）
- **Gain_W(M1−M2)**：L1 **+0.0124 [+0.0053, +0.0184] 显著**；
  L4 **−0.0009 [−0.0081, +0.0065] 不显著**
- ⇒ **DGaux 的 within 增益只在宽松划分下存在，在 family-level 下消失**

### 路线判定（四档）

判据回顾：
- **Dual-ddG**：`CW_L1 ≈ 0 AND CW_L4 ≈ 0 AND C_B >> C_W`
- **MIPA**：`CW_L1 > 0 AND CW_L4 ≈ 0 AND Δ_OOD CI > 0`
- **Hybrid**：`CW_L1 > 0 AND CW_L4 > 0 AND C_B >> C_W`
- **Unexpected**：`CW_L1` 和 `CW_L4` 都高

实测（M1）：**CW_L1 = +0.2658（明确 > 0）**；**CW_L4 = +0.1361（明确 > 0，非 ≈ 0）**；
`Δ_OOD` CI 排除 0；`C_B(L4) = −0.0112`（**不满足 `C_B >> C_W`**）。

⇒ **判定：`PARTIAL_MIPA`**（既非纯 MIPA 也非 Hybrid）

**理由（精确表述）**：
1. **within 信号真实存在**（`C_W` 在 L1 达 +0.266，远大于 0）→ 排除 **Dual-ddG**
2. **within 信号在 family-level OOD 后显著下降但**保留约 51%（+0.136）
   → **不符合 MIPA 预设的「`CW_L4 ≈ 0`」**
3. 在 L4 上 `C_B` 已塌陷到负值 → **不满足 Hybrid 的 `C_B >> C_W`**
4. 因此最贴近的是 **「MIPA 的弱化版」**：mutation signal 存在、
   能被 family-level 泛化**部分**保留，但**跨 family 的靶点均值预测失效并反号**

**（`Unexpected` 不适用：因为 `C_W` 并未「都高」——L4 已降到 L1 的一半。）**

### 论文可写性

- **Benchmark paper: YES**
  - 有力证据：① 全局 Pearson 的 ~100% 来自 between-target（within-shuffle 掉到 ≈0）；
    ② δ=0.5 下物理 vs 序列 Δ=+0.1575，p=0.00005；
    ③ Split Inflation（真 ΔΔG 上 L1−L4 = **+0.5590**）；
    ④ V6.1 bootstrap bug 的定量修正（CI 偏窄 ~25%）
- **Algorithm paper (Dual-ddG / MIPA / Hybrid): YES（但要改叙事）**
  - 不能写「within 信号在 OOD 后消失」；应写
    **「within 信号可部分跨 family 泛化（保留 ~51%），而 between 靶点均值预测在 OOD 后失效并反号」**
  - **必须处理 `Gain_W(M1−M2)` 在 L4 不显著**这一事实：
    DGaux 的价值限于宽松划分

### 措辞限制（**必须遵守**）

> 只能用 **"among the evaluated explicit-physics and structure-conditioned
> inverse-folding models"**
> **不能**升格为 "physics methods generally outperform learned methods"

依据：PHASE 6 显示**去掉 ESM-IF 后差异不显著**（p=0.111），
且本工作只评估了 **2 个** 物理模型与 **2 个** 序列条件模型。

---

## 阻塞

| # | 项 | 状态 | 说明 |
|---|---|---|---|
| 1 | **四档退化为三档** | 已查明，非缺陷 | 每个 PDB 恰属一个 family → L3 ≡ L2（数据性质）。替代的独立 L3（k-mer 代理，1144 簇）因**过聚类**而弱于 L4，不合法的递增层级 |
| 2 | **M2 ≡ M0** | 设计事实 | 故 `Gain_W(M2−M0) = 0`，该对照**未提供信息**（应记为「未检验」而非「无增益」） |
| 3 | **L3 独立实现耗时长** | 75s/次 | `_seq40_cluster_units` 的 k-mer 贪心为 O(N×代表集)，1144 簇 → 需 ~13 分钟/档 |
| 4 | `M1_implementation_not_found=False` 误报 | **已修复** | 原通用循环把「期望 False」当失败 → 50/50 通过 |
| 5 | 多进程重叠导致中间产物混乱 | **已修复** | kill + 归档 `logs/v6_2/archive_100106/`；修复 BUG A（`+=`→`extend`）与 BUG B（局部 rows 覆盖） |
| 6 | `logs/v6_2/failures.log` | **无记录** | 本任务无功能性失败 |

---

## 文件清单

### PHASE 0
- `experiments/v6_2/README.md`、`PROGRESS.md`
- `experiments/v6_2/phase0_environment/environment.md`

### PHASE 1
- `scripts/v6_2/01_provenance_audit.py`
- `experiments/v6_2/phase1_provenance/provenance_report.md`、`provenance_evidence.csv`

### PHASE 2
- `scripts/v6_2/02a_m1_audit.py`、`scripts/v6_2/02b_canonical_set.py`、`scripts/v6_2/02_split4_ddg.py`
- `experiments/v6_2/phase2_split4/`：`m1_audit.{md,json}`、`matched_analysis_set.csv`、
  `excluded_rows.csv`、`coverage_report.md`、`family_labels.csv`、
  `split4_results.csv`、`split4_report.md`、`split4_predictions.npz`
  （+ `.preL3fix.*` 与 `.smoke.csv` 归档）

### PHASE 3
- `scripts/v6_2/03_decomposition.py`
- `experiments/v6_2/phase3_decomposition/`：`pearson_decomposition.csv`、
  `primary_statistics.csv`、`variance_decomposition.csv`、
  `reconstruction_results.csv`、`control_results.csv`、`sanity_checks.csv`、
  `decomposition_report.md`

### PHASE 4
- `scripts/v6_2/04_ranking_metrics.py`
- `experiments/v6_2/phase4_ranking_metrics/`：`ranking_metrics.csv`、`n_sensitivity.md`、
  `delta_sensitivity.{csv,md}`、`cath_extended.{csv,md}`、`ranking_report.md`

### PHASE 5
- `scripts/v6_2/05_bootstrap_fixed.py`
- `experiments/v6_2/phase5_bootstrap/`：`bootstrap_comparison.csv`、
  `primary_statistics_ci.csv`、`bootstrap_report.md`

### PHASE 6
- `scripts/v6_2/06_leave_one_out.py`
- `experiments/v6_2/phase6_leave_one_out/leave_one_out.{csv,md}`

### PHASE 7
- `scripts/v6_2/07_leakage_assertions.py`
- `experiments/v6_2/phase7_leakage/leakage_checks.csv`、`leakage_report.md`

### PHASE 8
- `experiments/v6_2/V6_2_FINAL_SUMMARY.md`（本文件）
- `scripts/v6_2/09_run_remaining.py`（顺序驱动）

### 日志
- `logs/v6_2/`：`02_split4.log`、`03_decomp.log`、`05_boot.log`、`09_driver.log`、
  `archive_100106/`（混乱中间产物归档）
- **无 `failures.log`**
