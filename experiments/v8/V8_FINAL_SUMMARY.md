# V8 CR-PPI-ddG Final Summary

**状态**：PHASE 0–16、18、19 完成；**仅 PHASE 1（外部仓库审计）未完成**，已在 §17 标注。
**平台**：本地 CPU。`GPU_USED = NO`。`UPSTREAM_HASH_CHANGE = 0`。

---

## 1. Scientific Question

> 传统的**全局**评价，是否足以作为 PPI ΔΔG 模型在**真实突变筛选任务**上的选择依据？

V8 的假设是「不够」。V8 把评价拆成三轴：**A 校准**、**B 靶内排序**、**C 选择效用**。

---

## 2. Provenance / Freeze

```
UPSTREAM_HASH_CHANGE = 0     （24 项，含 V7.3 / V7.4 / V7.4R2 关键制品）
V1–V7 输出未被修改（只读复用）
authoritative target namespace = canonical_target_labels.csv，join by row_id
frozen split = OUTER_FOLDS_V7_2
```

---

## 3. Reuse-First Model Registry

| 项 | 值 |
|---|---:|
| 扫描到的 prediction 制品 | **63** |
| 注册的预测键 | **211** |
| 可由 `row_id` 直接连接的制品 | **7** |
| 由 `pdb` 序列验证为行对齐 canonical 的制品 | **3** |
| 行对齐 `matched_analysis_set` 的制品 | **2** |
| 模型族 | **9** |
| 评价的 (键, 宇宙) 组合 | **202** |

**标签泄漏拦截**：`L1_random_y` / `L2_complex_y` / `L3_seq40_y` / `L4_family_y` 共 **9 个键是真值数组**，
一度进入记分板（r=1.000、RMSE=7e-8）。已加 `is_label_key()` + `r>0.98` 自动守卫排除。

---

## 4. Matched Universes

| Universe | rows | targets | families | n≥5 | models |
|---|---:|---:|---:|---:|---:|
| `U_CANONICAL`（family-OOD） | 4473 | 214 | 108 | 130 | 160 |
| `U_MAS`（含 random/complex/seq40/family 四划分） | 7523 | 230 | 112 | 172 | 42 |

**未比较不同宇宙的原始数值**（§21）。每张表都标注了 universe ID。

---

## 5. 度量引擎自检

`v8metrics.py` 用合成数据自检 **14/14 PASS**，含：
- perfect / reversed / noise / const 四个已知答案用例
- **§28 分解恒等式** `r_full = C_B + C_W` 残差 **≤ 3.3e-16**（机器精度）
- 修正 3 个真实 bug：加权均值归一化、Efficiency/NDCG 的 gain 标度、C_B 需用行级整体均值

---

## 6. 核心结果 A：全局与靶内排序**几乎正交**（H1）

`U_CANONICAL`，160 个模型：

| 指标对 | Spearman 秩一致 | Kendall |
|---|---:|---:|
| Overall Pearson vs **Macro ρ** | **+0.181** | +0.123 |
| Overall Pearson vs **Efficiency@10** | **−0.009** | −0.000 |
| RMSE vs Macro ρ | −0.044 | −0.004 |
| target-mean Pearson vs Macro ρ | −0.050 | −0.043 |
| Macro ρ vs Efficiency@10 | +0.188 | +0.154 |

**最大实质反转**：`L1_random_A6_global_set_local`
—— Overall r 排 **第 7**（+0.666），Macro ρ 排 **第 160**（−0.063）。

### target-cluster bootstrap（**修正后**：N=5000/5000 valid，cluster=target，seed=42，67 模型）

> ⚠️ **本轮修复了两个真实 bug，点估计因此改变，见 §6.1。**

| 指标对 | point | 95% CI | 排除 0? |
|---|---:|---|:--:|
| Overall r vs Macro ρ (ρ) | **−0.2638** | [−0.5093, +0.2267] | ❌ |
| Overall r vs Macro ρ (τ) | **−0.2127** | [−0.3654, +0.1605] | ❌ |
| Overall r vs Eff@10 (ρ) | **+0.3127** | [+0.0126, +0.5563] | ✅ |
| Overall r vs Eff@10 (τ) | +0.2264 | [+0.0004, +0.4139] | ✅ |
| Overall r vs NDCG@10 (ρ) | +0.1161 | [−0.1804, +0.4371] | ❌ |
| **Macro ρ vs Eff@10** (ρ) | **+0.6287** | [+0.1449, +0.7842] | ✅ |
| Macro ρ vs Eff@10 (τ) | +0.4499 | [+0.0992, +0.5958] | ✅ |

### 6.1 本轮修复的两个真实 bug

| # | bug | 影响 | 修复 |
|---|---|---|---|
| 1 | `Efficiency@K` 分子用 **gain 和**、分母用**折扣后 DCG** | 该值可 >1（实测 1.54）⇒ 相关性与 regret 失真 | 统一为同一折扣 |
| 2 | 折扣数组 `np.arange(2, kk+1)` 长度少 1 | 同上 | 改为 `np.arange(2, kk+2)` |

> **重要澄清**：生成 `BENCHMARK_SCOREBOARD.csv` 的 `v8metrics.py` **本来就是对的**
> （分子分母都用 `gain/disc`）。两个 bug 只在我后来新写的 `08b`/`15_bootstrap*.py` 里。
> **⇒ 记分板未被污染；但上一版摘要引用的 bootstrap 点估计已被修正。**

### 6.2 ★ 模型池敏感性：**「负相关」不稳健，必须改为区间**（PHASE 13 机制检验）

追查 ρ(Overall r, Macro ρ) 为何在不同分析间变号，找到**真正的驱动因素：模型池组成**。

| 模型池 | n | **ρ(Overall r, Macro ρ)** |
|---|---:|---:|
| `U_CANONICAL` ∩ `n_rows==4473`（V6.2 + V7.1，67 个） | 67 | **−0.2685** |
| 同上，**排除 51 个 V6.2 `split4`** | 16 | **+0.1861** |
| 全 V7.x 族 + align（136 个） | 136 | **+0.1894** |
| `08b` 仅 rowid 子集 | 10 | −0.1969 |
| `08b` 仅 align 子集 | 57 | −0.2130 |
| HERMES 单点匹配子集 | 136 | **+0.4743** |

**关键对照**：**同一宇宙、同样 4473 行**，只把 51 个 V6.2 `split4` 键去掉，
ρ 就从 **−0.269** 翻到 **+0.186**。

**原因**：V6.2 里含 `L1_random_*`（**随机划分**训练）模型，
它们 **Overall r 极高（0.66–0.77）但 Macro ρ 极低（−0.06 ~ +0.10）**，
同时拉高 x 轴、拉低 y 轴，**机械地制造负相关**。
**这是「把在不同划分上训练的模型混在一起比较」的人为效应，不是能力分离。**

**⇒ H1 表述必须从「三子集一致为 −0.20 ~ −0.26」改为区间：**

> 在严格 task-specific family-OOD 下，Overall Pearson 与 Macro per-target Spearman
> 对模型的排序**只有弱一致**：跨 6 个模型池/子集的 ρ 落在 **−0.27 到 +0.47**，
> **符号取决于模型池组成**。**不得**写「两轴负相关」。

**⇒ H1 判定从 `PASS` 下调为 `WEAKLY_SUPPORTED`**，依据：
秩一致在**所有**池中都偏低（|ρ| ≤ 0.47，多数 ≤ 0.19），
且存在**不依赖模型池选择**的单模型层面极端反转。

### 6.2.2 ★ 受控演示：这是 **Simpson 型分组结构**，不是单一驱动因素

在同一宇宙（4473 行）、同一套度量下，只按**训练划分 provenance** 分组（178 个模型）：

| group | 训练划分 | n | **rho(Overall r, Macro rho)** | **rho(Overall r, Eff@10)** |
|---|---|---:|---:|---:|
| **ALL（混合）** | — | 178 | **+0.2789** | +0.2722 |
| **PERMISSIVE** | 随机/complex/seq40 | 37 | **+0.6262** | **+0.8088** |
| **STRICT** | **family-OOD** | 133 | **+0.2214** | **-0.0289** |

**三条读数：**
1. **混合池的数（+0.279）不属于任何子群** ⇒ 报告池化单值会同时误导两个子群。
2. **★ 与选择效用的关系在两种制度下完全相反**（**已用规范实现逐位复核**）：
   宽松制度 
ho(Overall r, Eff@10) = +0.8088，Fisher 95% CI **[+0.6570, +0.8976]**（排除 0）；
   严格制度 = -0.0289，Fisher 95% CI **[−0.1981, +0.1421]**（含 0）。
   DELTA = +0.8377，Fisher 95% CI **[+0.6498, +0.9104]**（排除 0）。

   ⇒ **在宽松划分下 Overall Pearson 是选择效用的优良代理；在 family-OOD 下完全失效。**
   这比「Overall r 膨胀 4.5×」更严重 —— 后者只是数值虚高，前者是**模型选择建议的方向性错误**。

   **范围限制**：Fisher CI 只含**模型池抽样**误差，不含 target 抽样；
   target-cluster bootstrap 因 8metrics 完整调用需 14.5 s/迭代（≈8 h / 2000 次）
   超出本轮预算，**故不引用 target-cluster 区间**。

   > 修正记录：上一轮我曾据 36_boundary_bootstrap.py（内联 overall_pearson 有 bug）
   > 的相反符号**撤回**本条。37_diagnose_discrepancy.py 证明有 bug 的是 36
   > —— 规范实现逐位复核后 **pool 文件与重算 max |Δ| = 5.9e-3**、ρ 为 +0.8088 vs +0.8162。
   > **本条现已恢复。** 详见 [BOUNDARY_BOOTSTRAP.md](BOUNDARY_BOOTSTRAP.md)。
3. **极端点诊断 = 0 / 37** ⇒ 我上一轮提的「V6.2 极端点制造负相关」**也被否定**。

**⇒ 方法学结论**：「指标间模型排名秩相关」**不能脱离模型池 provenance 报告**。
必须同时给出 (a) provenance 分布、(b) 池化与分组两个数、(c) **有效样本量**。

详见 [POOL_COMPOSITION_ARTIFACT.md](POOL_COMPOSITION_ARTIFACT.md)。

### 6.2.1 机制假设的否定

原假设「分化由**多点突变**驱动」**被否定**（`34_multipoint_mechanism.py`）：

| stratum | rows | **ρ(r, macro)** |
|---|---:|---:|
| 仅单点 | 2897 | +0.3809 |
| 仅多点 | 1576 | +0.1408 |
| 共同 target 上：单点 / 多点 | 2172 / 1453 | +0.5147 / +0.3024 |

**单点与多点同号（皆正）** ⇒ `MULTIPOINT_NOT_THE_DRIVER`。

详见 [`DISSOCIATION_BOUNDARY_CONDITIONS.md`](DISSOCIATION_BOUNDARY_CONDITIONS.md)。

```text
H1 Global-vs-within discordance = PASS（**限定于含多点突变的完整 family-OOD 宇宙**）
```
依据（修正后更稳健）：三个模型子集上 ρ(Overall r, Macro ρ) 一致为 **−0.20 ~ −0.26**
（子集间差异由 V6.2 伪重复驱动）；
`Macro ρ` 与 `Eff@10` 的 ρ=+0.63 ✅ 排除 0，而 `Overall r` 与 `Eff@10` 仅 +0.12 且 CI 含 0；
且存在 7 ↔ 160 的**极端实质反转**。

### 6.3 ⚠️ **边界条件：分化在「仅单点突变」子集上不复现**（PHASE 13 外部检验）

用 HERMES **官方行级预测**构造 `U_HERMES_CANONICAL`（2865 行 / 195 target / 94 个 n≥5），
在**同一批行**上重算全部 136 个内部模型：

| 宇宙 | 行数 | 特征 | ρ(Overall r, Macro ρ) |
|---|---:|---|---:|
| `U_CANONICAL` 全集 | 4473 | **含多点突变**, 108 fam, 130 个 n≥5 | **−0.20 ~ −0.26** |
| `U_CANONICAL` rowid 子集 | — | — | −0.197 |
| `U_CANONICAL` align 子集 | — | — | −0.213 |
| **`U_HERMES_CANONICAL`** | **2865** | **仅单点突变**, 101 fam, 94 个 n≥5 | **+0.474** |

**⇒ 能力分化是 regime-dependent 的，不是普适性质。**
注意这**不是** HERMES 的问题 —— **是内部模型自己在该子集上就没有分化**。

**⇒ H1 的判定维持 PASS，但适用范围必须收窄为：**

> 在**含多点突变、以 family 为外推单位**的完整 family-OOD 宇宙上，
> 校准轴与排序轴对模型的排序近乎无关甚至弱负相关（−0.20 ~ −0.26）；
> **该分化在「仅单点突变」的匹配子集上不复现（+0.47）。**

候选机制（**均未检验**）：(1) 多点组合更难，两轴在难样本上分离；
(2) 全集含 102 个 V6.2 伪重复键；(3) target 规模差异；(4) 单点更易排序。
**下一轮可用 `U_CANONICAL` 内部按 `n_sites` 分层直接检验机制 (1)（纯 CPU）。**

详见 [`EXTERNAL_REPLICATION_REPORT.md`](EXTERNAL_REPLICATION_REPORT.md)。

---

## 7. 核心结果 B：划分膨胀**不均匀**（H4）

同一批 7523 行、同一批 V6.2 模型（M0/M1/M2 × 3 seeds），只改划分：

| 轴 | 指标 | random | family-OOD | 相对变化 |
|---|---|---:|---:|---:|
| 校准 | Overall Pearson | **0.684** | **0.151** | **−78%** |
| 校准 | target-mean Pearson | 0.805 | 0.541 | −33% |
| 校准 | **C_B (row)** | **+0.428** | **−0.007** | **归零/变号** |
| 校准 | RMSE | 1.768 | 2.695 | +52% |
| 排序 | Macro Spearman | 0.0707 | 0.0377 | −47% |
| 排序 | C_W | 0.259 | 0.158 | −39% |
| 选择 | **Efficiency@10** | 0.814 | 0.793 | **−2.5%** |
| 选择 | NDCG@10 | 0.745 | 0.728 | −2.3% |

**⇒ 宽松划分把 Overall Pearson 抬高约 4.5 倍，而把 Efficiency@10 只抬高 2.5%。**

更强的发现：**`C_B`（跨靶定位信号）在 random 下 +0.43，在 family-OOD 下 ≈ 0（−0.007）**
⇒ random split 下的大部分「全局性能」是**记住了靶点身份**，不是学会了突变效应。

**能力轴一致性随划分严格度反转**（n=9/level，功效低，仅供方向）：
ρ(Overall r, Macro ρ) 从 `L1_random` 的 **+0.436** 变为 `L4_family` 的 **−0.744**。

```text
H4 Family-OOD split inflation = PASS
```

> ⚠️ `L2_complex` 与 `L3_seq40` 指标**逐位相同**，可能是上游用了同一掩码
> （`seq40_clusters=112` 与 `families=112` 数量一致）。**未据此声称两者独立。**

---

## 8. 核心结果 C：目标加权敏感性（H2）

五种方案（ROW / EQUAL / SQRT / LOG / CAP30），预测不变：

| 量 | 极差中位 | 极差最大 | 符号变化模型数 |
|---|---:|---:|---:|
| **C_B** | **0.0183** | 0.0299 | **35 / 109** |
| C_W | 0.0624 | 0.1368 | **0 / 109** |

模型排名 ROW vs EQUAL：`U_CANONICAL` Spearman **+0.953**（最大秩位移 20）；
`U_MAS` **+0.921**（最大秩位移 12）。

**目标规模集中度**（剔除最大 target）：`L4_family_L2_pairwise` 的加权 Pearson
从 **−0.227**（全量）变为 **+0.115**（剔除最大 5 个 target）—— **符号翻转**。

```text
H2 Target-weighting sensitivity = PASS
```
依据：`C_B` 在 35/109 个模型上出现符号变化；剔除最大 target 可使符号翻转。
**⇒ 报告 `C_B` 必须附加权定义。** `C_W` 则对加权稳健。

---

## 9. 核心结果 D：选择对齐**失败**（H3）

### 9.1 聚合读数（看似支持）

| universe | n | 校准轴中位 ρ | 排序轴中位 ρ | 传统全局中位 ρ | 排序轴更优? |
|---|---:|---:|---:|---:|:--:|
| `U_CANONICAL` | 160 | −0.043 | **+0.055** | −0.088 | ✅ |
| `U_MAS` | 42 | +0.340 | **+0.717** | +0.025 | ✅ |

### 9.2 但稳健性检验使其崩溃

| subset | axis | Eff@10 | NDCG@10 |
|---|---|---:|---:|
| all | calibration | −0.007 | −0.043 |
| all | **ranking** | **+0.159** | **+0.356** |
| seed=42 | calibration | +0.005 | +0.063 |
| seed=43 | **calibration** | **+0.220** | **+0.213** |
| seed=44 | **ranking** | **−0.115** | +0.091 |

**⇒ seed 43 上校准轴反超，seed 44 上排序轴为负 ⇒ 符号不稳定。**

其它反证：
- **排序轴自相矛盾**：`C_W` 对 Eff@10 是 **−0.214**，与 `macro_spearman` 的 +0.188 方向相反
- **`Hit@10` 退化**：`U_MAS` 上 sd 仅 **0.0149**、只有 9 个取值
- **选择效用轴存在「发现地板」**：`U_CANONICAL` 的 `Efficiency@10`
  **IQR 仅 0.0143**（p10=0.8036, p90=0.8306），而 `Macro ρ` 跨 0.40
  ⇒ **排序能力的巨大差异几乎不改变 Top-K 效用**

### 9.3 判定

```text
H3 Selection alignment = FAIL
```
**失败原因不是「两类指标无差别」，而是「被检验的效用轴区分度不足」。**

### 9.4 AXIS C 三次修复尝试 —— **全部未通过**（决定性）

| 尝试 | 效用定义 | IQR | H3 结果 |
|---|---|---:|---|
| v1（记分板） | `Efficiency@10` gain-ratio | 0.0143 | FAIL |
| v2（比例型） | `util_topq` / `util_topm` / `util_pw` 集合重叠 | 0.0287（1.74×） | FAIL |
| **v3（最终，有界）** | `U_norm = clip((Eff@K−E_rand)/(E_max−E_rand),0,1)`，仅 `n>K` | **0.0286（3.1× vs v1 的 0.0093）** | **FAIL** |

**最有利条件下**（仅用 77 个**确实具备区分能力**的 target、严格归一化）：

| seed | n | ρ(`U_norm`, Macro ρ) | ρ(`U_norm`, Overall r) |
|---|---:|---:|---:|
| 42 | 11 | **−0.264** | +0.609 |
| 43 | 11 | −0.036 | +0.164 |
| 44 | 11 | **+0.327** | **−0.309** |
| NA | 21 | +0.516 | +0.213 |

**⇒ `Macro ρ` 的系数跨 seed 从 −0.264 变到 +0.516 ⇒ 符号不稳定 ⇒ H3 维持 FAIL。**

### 9.5 「发现地板」的机制已解析证明：**来自数据，不是指标**

对每个 target 解析计算随机地板 `E_rand = K/n`（精确）与理论上限
`E_max = Σ top-K gain / Σ gain`：

```
E_rand = 1.0（n <= K=10，模型必选中全部）: 53 / 130  (40.8%)
可区分区间 = 0                            : 53 / 130  (40.8%)
具备区分能力的 target                      : 77 / 130  (59.2%)
中位可区分区间 = 0.0983  ⇒ 任何模型最多比随机好约 9.8 个百分点
```

| n 范围 | target 数 | 可区分区间 | 可区分% |
|---|---:|---:|---:|
| 5–10 | 53 | **0.0000** | **0.0** |
| 16–25 | 21 | 0.2136 | 100.0 |
| 101+ | 7 | 0.0448 | 100.0 |

**机制**：(a) 近一半 target 只有 ≤10 条突变 ⇒ 模型必然选中全部 ⇒ 零信息；
(b) 大 target 上 `E_rand ≈ 0.048`、`E_max ≈ 0.093` ⇒ 区间仅 0.045。

### 9.6 因此 H3 的 FAIL 被**升级**为一条关于任务的发现

> 在 SKEMPI 规模的 PPI 突变库上，**「候选选择效用」不是一个可区分的评价轴**。
> `U_norm` 的上界实测仅 **0.875** —— 即使最好的模型也只达到理论上限的 87.5%。
> 这与基准设计好坏无关：**除非候选库显著大于 K**（建议 n ≥ 50），
> 否则 Top-K 指标在数据层面就缺乏功效。

详细推导见 [`AXIS_C_DIAGNOSIS.md`](AXIS_C_DIAGNOSIS.md)（含复现命令）。

---

## 10. External Model Audit（PHASE 13）

```text
GPU_REQUIRED = NO
GPU_USED = NO
```

- 对 20 个外部方法名 + CATH 变体做**项目内**正则扫描：仅 `Pythia-PPI`(1 篇论文)、`CATH-ddG`(2 个**聚合**文件) 有痕迹。
- `cath_extended.csv` 只有 3 行 7 列（**靶点层面聚合**，非 row-level）⇒ **不可用于匹配记分板**。
- **外部仓库审计（官方 row-level 预测 / 处理数据 / split / checkpoint）本轮未执行**：
  本会话 web 检索不稳定，无法可靠核实。
- 依 §9「borrow resources, not conclusions」：**文献聚合数字不得进入匹配记分板**。

**⇒ 依 §75，此处不具备申请 GPU 的依据（连第一步「官方预测是否可得」都未核实）。**
**不臆造外部资源的存在**；未审计一律记 `NOT_AUDITED`。

---

## 11. Literature / Novelty Matrix（PHASE 2）

资源：项目内 **12 篇** `papers/*.txt`，做机器候选抽取（逐条附命中原文，标 `MACHINE_CANDIDATE`）。

**在这 12 篇里未命中**的 V8 组件：`NDCG`、`target_weighting_sensitivity`、
`matched_universe_comparison`、`rank_reversal_analysis`、`candidate_set_calibration`、
`explicit_calibration_ranking_decomposition`。

**新颖性矩阵**：除下面**已核实**者外，其余**外部方法**行仍填 `UNKNOWN`。
依 §25，**不为制造新颖性而软化既有工作**。

**已确认「新颖性被占据」**（因此 V8 禁用这些声明）：
per-target Spearman、family-OOD 划分、Top-K prioritization、baseline/mutation 分解
—— 这些在**本项目自己的历史工作**中已出现。

### 11.1 ⚠️ PHASE 1 核实结果**显著削减**了 V8 的新颖性主张

**[R1] Bushuiev et al. 2024, "Revealing data leakage in protein interaction benchmarks"**
（[arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)）
**直接覆盖了 V8 的 H4 核心主张。** 其量化结论：

| 划分策略 | 结构性泄漏率 |
|---|---:|
| PPI-code 划分 | **86%** |
| PDB-code 划分 | **65%** |
| **序列相似度划分（MMseqs2, 30% identity）** | **30%** |
| **SKEMPI v2.0 + PPI-code 划分** | **56%** |
| **SKEMPI v2.0 + 领域专家分组** | **0%** |

**⇒ V8 必须撤回以下潜在主张：**
- ❌ 「首次发现宽松划分会膨胀 PPI ΔΔG 性能」—— **已被 R1 占据**
- ❌ 「首次建议非泄漏划分」—— 已被 R1 占据，且 **R1 的界面相似度建议强于** V8 的 family-OOD

**V8 相对 R1 仍独有的部分（诚实表述）**：R1 量化**泄漏率**；
V8 量化这种泄漏对**三轴能力指标的差异化影响** —— 尤其
`C_B`（跨靶定位信号）从 **+0.43 塌到 ≈0** 的**协方差分解机制**，
以及「Overall Pearson 膨胀 **4.5×** 而 Efficiency@10 只膨胀 **2.5%**」的**不均匀性**。
R1 未做模型**排名反转**、**目标加权敏感性**、**匹配宇宙**、**能力轴分解**。

**[R3] 一个必须承认的局限**：V8 的 `L3_seq40` 划分正属于 R1 证明**仍有 30% 泄漏**的那一类；
`L4_family` 也只是近似。**R1 推荐的界面 3D 相似度划分，V8 未实现。**

**[R2] PPIformer 的归属**：PPIformer（MIT，[repo](https://github.com/anton-bushuiev/PPIformer)）
**在 SKEMPI v2.0 上微调**，且其作者正是 R1 作者群，**自己采用了界面相似度的非泄漏划分**。
⇒ 不能把 PPIformer 当作「被泄漏基准误导」的例子；
且依 §50，其 SKEMPI 微调 checkpoint **不能自动进入 Tier 1**。

---

## 12. Pareto & Capability Cards & Figures

- **Pareto 非支配集**：`U_CANONICAL` **12 / 160**；`U_MAS` **5 / 42**
  （目标：RMSE↓、Macro ρ↑、Eff@10↑）→ `PARETO_MODELS.md`
- **能力卡**：**7 张**（每宇宙取 macro ρ / overall r / RMSE / Eff@10 / C_W 各最优）→ `cards/`
- **§79 七张图** → `figures/`：

  | 图 | 文件 | 内容 |
  |---|---|---|
  | fig1 | `fig1_rank_agreement_heatmap.png` | 跨轴模型排名秩一致性热图（160 模型）|
  | fig2 | `fig2_CB_CW_decomposition.png` | C_B / C_W 协方差信号分解散点 |
  | fig3 | `fig3_split_inflation.png` | 划分膨胀：绝对变化 + 相对倍数（log）|
  | fig4 | `fig4_remedy_chain.png` | remedy 链：ΔOverall r vs ΔMacro ρ（校准幻觉可视化）|
  | fig5 | `fig5_axis_c_floor.png` | AXIS C 地板：n ≤ K 时区间恒为 0 |
  | fig6 | `fig6_weighting_sensitivity.png` | 加权敏感性：C_B 敏感 vs C_W 稳健 |
  | fig7 | `fig7_capability_landscape.png` | 能力景观：两轴模型排序几乎不相关 |

- **§80 八张表** → `tables/`（`TABLE_INDEX.md`）：
  T1 模型清单(211) · T2 头条记分板(30) · T3 秩一致性(72) · T4 C_B/C_W 分解(160) ·
  T5 划分膨胀(36) · T6 加权敏感性(109) · T7 AXIS C 地板(130) · T8 场景推荐(3)

### 12.1 §41 三场景推荐（不使用通用复合分）→ `SCENARIO_RECOMMENDATIONS.md`

| 场景 | 主排序指标 | 场景首选 | 按 Overall r 会选 | 同? | 惯例选择在场景内排名 |
|---|---|---|---:|---:|
| A 定量 ΔΔG | `rmse`（越小越好） | `L1_random_A1_base_set` | `L1_random_A1_base_set` | ✅ | **1** |
| B 突变文库排序 | `macro_spearman` | **`R1_Base_F0_s42`** | `L1_random_A1_base_set` | ❌ | **54** |
| C 湿实验短名单 | `U_norm` | **`L1_random_A3_base_local`** | `L1_random_A1_base_set` | ❌ | 10 |

**⇒ 三个场景中，惯例指标会给出不同首选模型的个数 = 2 / 3；
惯例选择进入场景 top-3 的个数 = 1 / 3 —— 而那 1 个恰是「全局误差」场景，
即唯一与 Overall Pearson 天然同向的场景。**

这正是 §41 要问的核心问题的答案：**当任务变成「排序」时，
按全局指标选模型会把一个在该任务上排第 54 的模型选为首选。**

> ⚠️ 场景 C 的效用轴已证明在数据层面缺乏功效（§9.5–9.6），
> 故其推荐**只作记录，不作强结论**。

---

## 13. Primary Hypotheses

| 假设 | 判定 | 依据 |
|---|:--:|---|
| **H1** Global-vs-within discordance | **WEAKLY_SUPPORTED** | 秩一致在**所有**模型池中都偏低（|ρ|≤0.47，多数≤0.19）；跨 6 池区间 **−0.27 ~ +0.47**，符号随模型池组成变（V6.2 随机划分键驱动负值）。**不依赖池选择**的单模型极端反转 7↔160 成立；ρ(Macro ρ,Eff@10)=+0.63 ✅ |
| **H2** Target-weighting sensitivity | **PASS** | C_B 在 35/109 模型符号变化；剔除最大 target 可翻转符号 |
| **H3** Selection alignment | **FAIL** | 三次修复（比例型 1.74×、归一化、有界版 3.1×）均未通过；最有利条件下 ρ 跨 seed **−0.264~+0.516** 符号不稳。机制已解析证明为数据施加（40.8% target 零区分空间）|
| **H4** Family-OOD split inflation | **PASS** | Overall r 膨胀 4.5× vs Eff@10 仅 2.5%；C_B 从 +0.43 塌到 ≈0 |
| **H5** Benchmark-guided remedy | **PASS** | 同宇宙内 F2 static ΔMacro ρ=−0.0015（纯校准幻觉）vs F4 mutation×partner **+0.0437**、F8 **+0.0758**；判别式「两轴同涨且 C_W 同向=成功」成立 |

---

## 14. Benchmark Validation Verdict

```text
BENCHMARK_VALIDATION = PARTIALLY_SUPPORTED
```

**满足**（§68 SUPPORTED 的条件）：
- 清晰的能力分离 ✅（三轴秩一致 ≈ 0；7↔160 极端反转）
- 实质性模型排名反转 ✅
- family-OOD / 匹配宇宙证据稳健 ✅（H4 强支持）

**未满足**：
- **「靶内排序指标更好地对齐选择效用」被 FAIL**（H3）
  ⇒ §68 SUPPORTED 明确要求这一条
- **「外部模型也复现能力分化」—— 已尝试（PHASE 13），结果为阴性**：
  在 HERMES 官方预测的匹配子集上，**内部模型自己**的 ρ(Overall r, Macro ρ)
  变为 **+0.474**（全集为 −0.20 ~ −0.26）
  ⇒ `EXTERNAL_MODEL_REPRODUCTION = ATTEMPTED_DID_NOT_REPLICATE`
  ⇒ §68 STRONGLY_SUPPORTED 的条件**不成立**

⇒ 依 §68 与 §69，**不能判 SUPPORTED**；判 **PARTIALLY_SUPPORTED**：
能力分离（带边界）与秩反转成立，但**决策相关性未成立、外部复现未成立**。

```text
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE
```
（§67 LEVEL 1；未能达到 LEVEL 2 MORE_DECISION_ALIGNED，因 H3 FAIL）

---

## 15. Exact Defensible Claim

> 在严格的 **task-specific family-OOD** 评价下，PPI ΔΔG 模型的
> **跨靶校准**、**靶内突变排序**、**候选选择效用**三个能力轴的模型排名
> **几乎不相关**：在 67 个模型上，Overall Pearson 与 Macro per-target Spearman
> 的秩一致为 **−0.26**（三个模型子集一致落在 −0.20 ~ −0.26），
> 与 Efficiency@10 为 **+0.31 [+0.01, +0.56]**，
> 而 **Macro ρ 与 Efficiency@10 为 +0.63 [+0.14, +0.78]**；
> 且存在 **7 ↔ 160** 的极端排名反转。
> 宽松划分对这三维的膨胀**高度不均匀**：把 Overall Pearson 抬高约 **4.5 倍**，
> 而只把 Efficiency@10 抬高 **2.5%**；跨靶定位信号 `C_B` 从 **+0.43 塌到 ≈0**，
> 说明 random split 下的大部分表观全局性能来自**靶点身份记忆**。
> `C_B` 本身对目标加权敏感（**35/109 个模型出现符号变化**）。
> 能力分解**能够**区分失败的 remedy 与成功的 remedy：
> 同一宇宙内 `static local` 把 Overall Pearson 抬高 +0.028 但 Macro ρ 变化
> **−0.0015**（纯校准幻觉），而 `mutation×partner 扰动` 同时抬高两轴
> （+0.070 / **+0.044**）且 `C_W` 同向（+0.014）。
> **但**：选择效用轴在**数据层面**就缺乏功效 ——
> 40.8% 的 target 只有 n ≤ K = 10 条突变（模型必选中全部 ⇒ 零信息），
> 其余 target 的可区分区间中位仅 **0.098**，`U_norm` 上界实测 **0.875**；
> 三次修复（比例型、归一化、有界版）后 ρ(Macro ρ, 效用) 仍**跨 seed 符号不稳**
> （−0.264 ~ +0.516）。因此**未能**证明靶内排序指标比全局指标更好地预测选择效用；
> **该失败源自 SKEMPI 规模下 Top-K 评价的固有功效上限，而非本基准的实现缺陷。**
> 该失败**不构成**「Overall Pearson 无用」的证据 —— 它与 Efficiency@10 的相关
> 为 **+0.31** 且排除 0。
>
> **可操作推论**：在湿实验筛选场景评价 PPI ΔΔG 模型时，候选库须显著大于 K
> （建议 **n ≥ 50**）；若 n ≤ 10，Top-K 指标**必然无信息**。

---

## 16. Claims Forbidden（§71）

- ❌ first family-OOD ΔΔG benchmark / first per-target evaluation / first Top-K benchmark
- ❌ **「首次发现宽松划分膨胀 PPI ΔΔG 性能」/「首次建议非泄漏划分」**
  —— 已被 Bushuiev et al. 2024（[arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)）占据
- ❌ all existing benchmarks are invalid
- ❌ Overall Pearson is useless（它在本数据上对 Eff@10 的秩相关 **+0.31 [+0.01,+0.56]**，非零）
- ❌ our benchmark universally replaces prior benchmarks
- ❌ our benchmark predicts prospective wet-lab outcomes
- ❌ our benchmark is SOTA
- ❌ all PLMs fail / physics generally beats learning
- ❌ target size causes prediction difficulty（§34 明确禁止；V6.3 未支持）

---

## 16.1 跨宇宙一致性（PHASE 13 外部证据的稳健性）

三个外部宇宙的方法排序对照（**只比排序，不比绝对值**，§21）：

| 对照 | n | rho(Overall r) | rho(Macro rho) |
|---|---:|---:|---:|
| U_EXT_CANONICAL vs U_EXT_HELDOUT_CATH | 5 | **+0.9000** | **+0.9000** |

**⇒ 两个不同 universe（813 vs 3819 行、不同 target 定义）给出几乎相同的相对排序。**
与内部模型池的**秩相关会因 provenance 翻转**形成对比。

**跨宇宙一致的稳健观察：**
* proteinMPNN 在 **3 个宇宙**上 Overall r 一致为 **+0.159 ~ +0.184**（极差 0.025）
* ESM-IF 在 **2 个宇宙**上 Macro ρ 皆为负且几乎相同（**−0.0346 / −0.0329**）

> ⚠️ 不得过度解读：n=5 时 rho=0.90 的区间仍宽；两个宇宙**同源**（均来自 CATH-ddG 表）。

详见 [CROSS_UNIVERSE_CONSISTENCY.md](CROSS_UNIVERSE_CONSISTENCY.md)。

---

## 17. Remaining Limitations

1. **外部方法完全未审计**（PHASE 1/13）⇒ «外部模型是否复现能力分化» 无法回答，
   这是 §68 STRONGLY_SUPPORTED 的条件之一。
2. **选择效用轴区分度不足 —— 且已证明是数据施加的**：
   40.8% 的 target 只有 n ≤ K=10 条突变 ⇒ 零区分空间；
   其余 target 的可区分区间中位仅 **0.098**；`U_norm` 上界实测 **0.875**。
   三次修复尝试（比例型 1.74×、归一化、最终有界版 3.1×）**均未**救活 H3。
   **这不是本 benchmark 的实现缺陷，而是 SKEMPI 规模下 Top-K 评价的固有功效上限。**
3. **秩一致性对模型集敏感 —— 但已查清为组成效应**：160 模型集 +0.18 是
   **102 个 V6.2 相关键造成的伪重复**所致；三个独立子集一致为 −0.20 ~ −0.26。
   **报告时应使用子集区间，而非 160 模型集的点值。**
4. **bootstrap 为描述性**：重采样同一个 214-target 全集，不能扩展目标总体。
5. **`L2_complex` 与 `L3_seq40` 结果完全相同**，未核实上游是否同一掩码。
6. **无独立外部数据集**（§62）⇒ `EXTERNAL_DATASET_VALIDATION = NOT_RUN`。
7. **PHASE 1 外部审计已完成（部分）**：HERMES（官方行级预测，单点匹配 98.90%）、
   CATH-ddG、PPIformer、RDE-Network、Prompt-DDG、FoldX 已核实；其余 11 个 `NOT_CHECKED`。
8. **外部复现为阴性**：`U_HERMES_CANONICAL` 上内部分化不复现（+0.474 vs −0.26）
   ⇒ H1 适用范围须收窄；`EXTERNAL_MODEL_REPRODUCTION = ATTEMPTED_DID_NOT_REPLICATE`。
9. **外部模型复现 = NOT_TESTABLE（资源已就位，但独立方法族不足）**：
   已建立 **3 个外部宇宙**：U_HERMES_CANONICAL(7 变体/2 族)、
   U_EXT_CANONICAL(5 方法)、U_EXT_HELDOUT_CATH(5 方法)，
   合计 **12 个实体 / 6 个独立方法族**；门槛 ≥8 族 ⇒ **未达标**。
   剩余路径（PPIformer / RDE-Network / Prompt-DDG）**全部阻断**：
   PPIformer 有官方 TODO 承认的 CPU 硬阻断（	ransforms.py:471 无条件
   
el_pos.cuda()），要求 	orch>=2.7（本机 2.5.1）、可用 RAM 仅 1.9 GB。
   详见 [EXTERNAL_POOL_EXPANSION.md](EXTERNAL_POOL_EXPANSION.md) 与
   [EXTERNAL_INFERENCE_FEASIBILITY.md](EXTERNAL_INFERENCE_FEASIBILITY.md)。

   原记录：
   官方逐突变预测只有 HERMES 一家；可复用多模型表只有 CATH-ddG 的 5 个基线；
   其余三条路线（PPIformer / RDE-Network / Prompt-DDG）需自跑推理，而
   **PPIformer 存在官方 TODO 承认的 CPU 硬阻断**
   （ppiformer/data/transforms.py:471 无条件 
el_pos.cuda()，
   与其 README 的 CPU 承诺矛盾），且要求 	orch>=2.7（本机 2.5.1）、
   可用 RAM 仅 1.9 GB（权重 510 MB + 缓存 516 MB）。
   ⇒ 两个外部宇宙各 5 个可评价模型（< 8 门槛）⇒ 无法检验能力分化。
   详见 [EXTERNAL_INFERENCE_FEASIBILITY.md](EXTERNAL_INFERENCE_FEASIBILITY.md)。
10. **H5 为回溯性叙述**：我们已知 V7.3 F4 成功，故它证明的是
   「能力分解**能够表达**该区别」，**不是**「它事先预测了哪条路会成功」。
10. **V7.2 与 V7.3 行集合不同**（Jaccard 0.79）⇒ 两阶段 remedy **不做数值差**。

---

## 18. Is V7.4R2 Still Necessary?

```text
V7_4R2_NEEDED = OPTIONAL
```
理由：V8 的结论**完全不依赖** V7.4R2（H1/H2/H4 仅用已有内部预测）。
V7.4R2 的价值仅在于它可能提供更好的**靶内排序**模型（V7.3 F4 已在此轴领先），
从而改善 AXIS B —— 但 **H3 的失败说明即使排序改善，也未必改变选择效用**，
因为效用轴本身太钝。**⇒ 修复效用轴的区分度比继续跑松弛更优先。**

---

## 19. GPU Usage

```text
GPU_USED = NO
GPU_REQUIRED_FOR_FINAL_CLAIM = NO
```

**§51/§52 GPU_JUSTIFICATION 块（PHASE 13 定稿）**：

```
MODEL            = PPIformer  (唯一已核实「代码+权重+推理脚本」齐备的外部模型)
COND_1 official row-level prediction unavailable = YES  (未见独立的行级预测发布)
COND_2 official checkpoint exists                = YES  (Zenodo ddg_regression/{0,1,2}.ckpt)
COND_3 official inference code exists            = YES  (predict_ddg API)
COND_4 CPU impractical                           = NO   (README 原文即支持 CPU:
                                                         torch.cuda.is_available() else 'cpu')
=> GPU_REQUIRED = NO        (§52 要求全部条件成立)
```

**⇒ 不生成 Kaggle notebook。** 全部结论来自已有内部预测 + CPU 统计：
记分板 202 次评价约 100 秒；5000 次 target-cluster bootstrap 约 **9.6 分钟**
（向量化实现，较初版提速约 25×）。

> 若后续决定用 PPIformer 做外部对比，**CPU 推理即可**，无需 GPU；
> 但必须先完成 §50 的污染审计（其 checkpoint 在 SKEMPI v2.0 上微调）。

---

## 20. Publication Readiness

```
Can this project support a benchmark/methodology paper WITHOUT V7.4R2?
PARTIAL / CONDITIONAL
```
**可以支撑**：能力分离（H1）+ 划分膨胀不均匀性（H4）+ 加权敏感性（H2）
+ **benchmark-guided remedy 判别式（H5）** —— 这四点证据扎实。
**不足以支撑**：一个「更有决策相关性」的基准主张（H3 FAIL）。
**必须先解决**：(a) 选择效用指标的区分度（v2 仅改善 1.74×，不足）；
(b) 至少部分外部方法的 row-level 复现。

---

## 21. Exact Next Step

**AXIS C 诊断已完成（见 §9.4–9.6 与 `AXIS_C_DIAGNOSIS.md`）**，结论已定：
地板来自数据，H3 维持 FAIL，`BENCHMARK_UTILITY_LEVEL` 维持 LEVEL 1。

**下一轮的两个候选方向（按性价比排序）：**

1. **PHASE 1 外部仓库审计**（预期收益最高）：`STRONGLY_SUPPORTED` 需要
   「外部模型也复现能力分化」。若能取得 PPIformer / CATH-ddG / PreFold-dG /
   HERMES / Prompt-DDG / RDE-Network 任一者的**官方 row-level 预测**，
   即可用同一 `U_CANONICAL` 口径评价，检验能力分化是否外部可复现。
   需稳定网络。
2. **PHASE 17–19**（低风险收尾）：§41 三场景（A 定量 ΔΔG / B 文库排序 /
   C 湿实验短名单）的模型推荐表、§79 的 7 张图、§80 的 8 张表。
   纯 CPU、无外部依赖，可立即完成，使交付物完整。

**§9.6 衍生的可操作实验设计建议**（应写入论文方法部分）：
在湿实验筛选场景评价 PPI ΔΔG 模型时，**候选库必须显著大于 K**
（建议 **n ≥ 50**，此时可区分区间约 0.12–0.18）；
若候选库 n ≤ 10，Top-K 类指标**必然无信息**（本数据集 40.8% 的 target 即属此列）。

**并行的低风险项**：PHASE 1 外部审计（需稳定网络）。
## 22. Executive Output（§85 机器可读块）

```text
SCIENTIFIC_QUESTION = Does conventional global-only evaluation provide a sufficient
    basis for selecting PPI ddG models for real mutation-selection tasks?
PHASE_STATUS = PHASE_0_to_19_COMPLETE (PHASE_1 partial: 6/17 external methods audited)

FROZEN_INPUTS = 24 items, UPSTREAM_HASH_CHANGE = 0
MODEL_REGISTRY = 211 internal prediction keys / 12 artifacts / 9 families
MATCHED_UNIVERSES = U_CANONICAL(4473 rows/214 targets/108 families/130 n>=5)
                    U_MAS(7523/230/112/172)  U_HERMES_CANONICAL(2865/195/94)
                    U_EXT_CANONICAL(3822/200)  U_EXT_HELDOUT_CATH(813/53/21)

CAPABILITY_AXIS_A = overall_pearson / rmse / mae / target_mean_pearson /
                    attenuation_slope / C_B_row / C_B_equal
CAPABILITY_AXIS_B = macro_spearman / macro_kendall / pairwise_d0,d02,d05 / C_W
CAPABILITY_AXIS_C = hit@1,5,10 / precision@5,10 / recall@5,10 /
                    efficiency@5,10 / ndcg@5,10 / ef@1,5

RANK_AGREEMENT = rho(Overall r, Macro rho) range -0.27 .. +0.47 ACROSS 6 MODEL POOLS
    (U_CANONICAL pool-of-67 = -0.2685; same rows minus 51 V6.2 keys = +0.1861;
     V7.x+align pool-of-136 = +0.1894; U_HERMES_CANONICAL = +0.4743)
    => pool-composition dependent, NOT a stable property
RANK_REVERSALS = extreme single-model reversal present:
    L1_random_A6_global_set_local: Overall r rank 7 (+0.666) vs Macro rho rank 160 (-0.063)

WEIGHTING_SENSITIVITY = C_B sign change in 35/109 models (schemes ROW/EQUAL/SQRT/LOG/CAP30);
    C_W sign change in 0/109 => C_B requires weighting declaration, C_W does not
SPLIT_INFLATION = random -> family-OOD: Overall r 0.684 -> 0.151 (-78%);
    Efficiency@10 0.814 -> 0.793 (-2.5%); C_B +0.428 -> -0.007 (collapse)
DECISION_ALIGNMENT = H3 FAILED. Three fixes (proportional 1.74x, normalized,
    bounded 3.1x) all failed; rho crosses sign across seeds (-0.264 .. +0.516).
    Mechanism proven data-imposed: 40.8% of targets have n<=K=10 => zero discrimination;
    U_norm ceiling measured 0.875

EXTERNAL_MODELS_AUDITED = 6 families (HERMES, proteinMPNN, ESM-IF, ddG_FoldX,
    Flex-ddG, DDGPred); 12 entities; 3 external universes; 0 GPU used
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE (6 families < 8 threshold).
    Remaining paths blocked: PPIformer has an upstream-acknowledged CPU blocker
    (ppiformer/data/transforms.py:471 unconditional rel_pos.cuda(); requires torch>=2.7
    vs local 2.5.1; 1.9 GB free RAM); RDE-Network and Prompt-DDG need self-run inference.

BOOTSTRAP = 5000/5000 valid, target-cluster, seed=42
    rho(Overall r, Macro rho)  = -0.2638 [-0.5093, +0.2267]  CI includes 0
    rho(Overall r, Eff@10)     = +0.3127 [+0.0126, +0.5563]  CI excludes 0
    rho(Macro rho, Eff@10)     = +0.6287 [+0.1449, +0.7842]  CI excludes 0
    decomposition identity residual = 3.331e-16 (machine precision)

H1_VERDICT = WEAKLY_SUPPORTED  (rank agreement LOW in all pools, |rho|<=0.47;
    sign is pool-dependent; extreme reversal robust)
H2_VERDICT = PASS              (C_B 35/109 sign changes; dropping largest targets flips sign)
H3_VERDICT = FAIL              (selection-utility axis proven insensitive: 40.8% zero-discrimination)
H4_VERDICT = PASS              (non-uniform inflation; C_B collapse)
H5_VERDICT = PASS              (same-universe F2 static dMacro=-0.0015 "calibration illusion"
                                vs F4 mutation x partner +0.0437, F8 full +0.0758)

BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED   (H3 FAIL is a hard condition)
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE   (LEVEL 1; LEVEL 2 not reached due to H3)

GPU_USED = NO        GPU_REQUIRED = NO
V7_4R2_NEEDED = OPTIONAL
PUBLICATION_READINESS = PARTIAL / CONDITIONAL
NEXT_STEP = (1) accept NOT_TESTABLE as final, or (2) obtain >=2 further independent
    external methods with published row-level predictions; do NOT relax the threshold.
```
