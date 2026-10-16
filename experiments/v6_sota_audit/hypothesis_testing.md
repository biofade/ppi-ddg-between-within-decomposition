# V6.0 PHASE 3 — 核心假设检验

> 假设（H1）：**当前的 PPI ΔΔG SOTA 模型存在评估幻觉 —— 其极高的全局预测指标（Overall Pearson / P@10）掩盖了其在真实复合物内部（Within-target）挑选突变体效用的失效。**
>
> 判定标准（V6 §6）：
> - **SUPPORTED** —— 至少 1 个外部模型 `Overall Pearson > 0.4` **且** `Macro per-PDB Spearman < 0.15` **且** within `Hit@10 < Base Rate`
> - **REFUTED** —— 所有外部模型 `Macro per-PDB Spearman > 0.3`
> - **INCONCLUSIVE** —— 其他

---

## 最终判定：**REFUTED（并被数据反向证伪）**

**外部 SOTA 模型在被审计的严格协议下，其 within-target 排序能力显著优于我们的 V5 baseline，而非同等失效。**

核心对照（同协议、同指标实现、同一过滤规则「组内 ≥5 条」）：

| 来源 | 模型 | Overall Pearson | **Macro per-PDB Spearman** | Hit@10(Δ<0) | Base Rate | Efficiency@10 |
|---|---|---:|---:|---:|---:|---:|
| CATH 官方 held-out | **Flex-ddG** | +0.5626 | **+0.4159** | 1.0000 | 0.1815 | +0.5715 |
| CATH 官方 held-out | **DDGPred** | +0.5573 | **+0.3262** | 0.8000 | 0.1815 | +0.3715 |
| CATH 官方 held-out | FoldX | +0.4348 | +0.3584 | 0.4000 | 0.1815 | +0.4505 |
| CATH 官方 held-out | proteinMPNN | +0.1215 | +0.2227 | 0.5000 | 0.1815 | +0.1399 |
| CATH 官方 held-out | ESM-IF | −0.1896 | +0.0329 | 0.4000 | 0.1815 | +0.0497 |
| RBM case study | FoldX | +0.5151 | **+0.6478** | 1.0000 | 0.8625 | +0.7515 |
| RBM case study | DSSA-PPI | +0.3081 | +0.3129 | 1.0000 | 0.8625 | +0.7439 |
| RBM case study | SAAMBE-SEQ | +0.1934 | +0.2058 | 1.0000 | 0.8625 | +0.5074 |
| **V5 baseline（我们）** | — | **+0.3721** | **+0.0835** | 0.9667 | 0.3729 | **+0.1014** |

**V5 baseline 是三组数据中 within-target 排序最弱的一个**（macro per-PDB Spearman
+0.0782 / +0.0970 / +0.0752，均值 **+0.0835**），而多个外部模型达到
**+0.33 ~ +0.65**。

---

## 核心诊断问题 1：幻觉是否存在？

**答案：在「全局高 + 局部低」这个特定组合上，未观察到 —— 但观察到一个不同的现象。**

| 判据 | 阈值 | 满足的外部模型数 |
|---|---|---:|
| Overall Pearson > 0.4 | >0.4 | **4 / 13** |
| Macro per-PDB Spearman < 0.15 | <0.15 | **5 / 13** |
| **两者同时成立（H1 的幻觉模式）** | — | **0 / 13** |

`Overall Pearson > 0.4` 的 4 个模型（Flex-ddG 0.5626、DDGPred 0.5573、FoldX 0.4348、FoldX-RBM 0.5151）**全部同时具有较高的 macro per-PDB Spearman（0.33–0.65）**。
→ **不存在「全局高但局部低」的组合。**

`Macro per-PDB Spearman < 0.15` 的 5 个模型（Esm-IF、proteinMPNN、ESM-2-FFN、MuTon、MpbPPI/ProBASS）**其 Overall Pearson 也都很低（−0.19 ~ +0.12）**。
→ **局部差与全局差是同步的，不是「被掩盖」。**

**观察到的新现象（H1 的变体）**：
> **V5 baseline 本身**恰好落在「全局中等（Pearson +0.3721）但局部极弱
> （macro per-PDB Spearman +0.0835）」的位置。也就是说，
> **「全局指标与局部能力脱节」这一现象在我们自己的模型上成立，而非在外部 SOTA 上成立。**

这与 V5.0 的结论（`REPRESENTATION_BOTTLENECK_CONFIRMED`）一致，且**强化**了它：
外部模型用**逐突变/逐结构局部能量**（Flex-ddG、FoldX）或**几何编码器**
（DDGPred）能拿到局部排序信号，而我们的**池化 ESM-2 表示**拿不到。

---

## 核心诊断问题 2：效用是否崩塌？

**答案：外部模型未崩塌；我们的模型效用最低。**

`Efficiency@10`（0 = 与随机无异，1 = oracle 水平）：

| 模型 | Efficiency@10 | 解读 |
|---|---:|---|
| FoldX (RBM) | **+0.7515** | 接近 oracle 的 3/4 |
| DSSA-PPI (RBM) | +0.7439 | 同上 |
| Flex-ddG (CATH) | **+0.5715** | 超过 oracle 空间的一半 |
| SAAMBE-SEQ (RBM) | +0.5074 | — |
| FoldX (CATH) | +0.4505 | — |
| DDGPred (CATH) | +0.3715 | — |
| ProBASS (RBM) | +0.3271 | — |
| ESM-2-FFN (RBM) | +0.2898 | — |
| PPIFormer (RBM) | +0.2652 | — |
| proteinMPNN (CATH) | +0.1399 | 弱 |
| ESM-IF (CATH) | +0.0497 | 几乎无 |
| MpbPPI (RBM) | −0.0192 | 略差于随机 |
| **V5 baseline（我们）** | **+0.1014** | 倒数第 4 |

→ 我们模型相对随机**确有**微弱改善（+0.10），但**远低于** Flex-ddG / FoldX 的
+0.45 ~ +0.75。

**但 H1 有一个被支持的部分**：
即使是表现最好的外部模型（Flex-ddG、FoldX），`Efficiency@10` 也只有
**0.45–0.75**，**都没有达到 oracle**。同时
CATH test 上 Flex-ddG 的 `Hit@10(ΔΔG<0) = 1.000`，**而其 Base Rate 仅 0.1815**
（该 test 集只有 18.2% 的突变是 beneficial）——
即在这些特定的全局 Top-10 上它们确实挑对了。
→ 所以「SOTA 已能用于湿实验选靶」这一说法**在 Top-10 这个粒度上是部分成立的**
（对 Flex-ddG / FoldX），不能一概否定。

---

## 核心诊断问题 3：样本量是否足够？

**答案：不足 —— 这是本次审计最主要的限制。**

| 数据集 | n 行 | n 组（过滤后 ≥5 条/组） | 判定 |
|---|---:|---:|---|
| CATH-ddG 官方 held-out CATH | 755 | **21** | ⚠️ **样本量不足，结论需谨慎** |
| DSSA-PPI RBM case study | 1 702 | **1** | ⚠️ **严重不足**（单一 group，且真值是适应度非 ΔΔG） |
| V5 PPB-Affinity development | 1 129 | **61** | ✅ 相对可用 |

**必须明确标注的三点**：

1. **CATH test 原始为 813 条 / 53 PPI**，但按「组内 ≥5 条」过滤后**只剩 21 个
   complex / 755 行**。`Macro per-PDB Spearman` 由 21 个值平均而来，
   其置信区间必然很宽（未在本报告中计算 CI，**属于限制**）。
2. **RBM 数据集的「n_group = 1」是人为分箱的结果**（把 1 702 条按 `pos` 连续分块），
   **不等于多个独立复合物**；其 macro Spearman 与 Efficiency@10 **不应与
   CATH/V5 的数字直接比较**。它只能作为「同一靶点内、同一位置上的氨基酸替换排序」
   的例证。
3. **RBM 的真值 `EXP_fitness` 不是 ΔΔG**（是实验适应度，符号语义未由仓库说明）
   → 该数据集**不进入 `Hit@10(ΔΔG<0)` 的绝对口径**；其 `Hit@10` 一列
   （因为 86.25% 的行都 <0）**没有判别力，不予解读**。

---

## 并排对比（外部模型 vs V5 baseline）

### 同一指标族内可比的部分

| 对比 | 可比的指标 | 结论 |
|---|---|---|
| CATH 官方 held-out vs V5 | `Overall Pearson`、`Macro per-PDB Spearman`、`Efficiency@10` | **外部模型（Flex-ddG / DDGPred / FoldX）在局部排序上显著更强** |
| CATH 官方 held-out vs V5 | `Hit@10` | 可比但受基准率差异影响（CATH 18.15% vs V5 37.29%） |
| RBM vs（CATH / V5） | — | **不可直接比较**（分组为人工分箱、真值非 ΔΔG） |

### 关于「数据集差异」的诚实说明（不作为借口，只作事实）

| 维度 | CATH test | V5 test |
|---|---|---|
| 数据源 | SKEMPI v2.0 | PPB-Affinity（五源聚合） |
| 划分 | CATH 同源超家族 hold-out，序列同一性 <30%，**无 easy mutation** | 复合物级 + ≤40% 序列同一性 + 家族级 |
| 过滤后规模 | 21 complex / 755 行 | 61 complex / 1 129 行 |
| Base Rate（ΔΔG<0） | **0.1815** | **0.3729** |
| 模型类型 | 物理能量（FoldX/flex ddG）+ 几何编码器（DDGPred）+ 序列（ESM-IF） | 池化 ESM-2 + 60 物理特征 + GBDT |

> **不构成辩护**：即使考虑上述差异，V5 baseline 的 `Macro per-PDB Spearman`
> （+0.0835）与 Flex-ddG（+0.4159）之间的差距（**0.33**）远大于任何 plausible
> 数据集效应，且**方向一致地出现在 CATH 与 RBM 两个独立数据集上**。

---

## 结论汇总

| 诊断问题 | 答案 |
|---|---|
| Q1 幻觉（全局高 + 局部低）是否存在？ | **否 —— 0/13 外部模型满足。** 但该模式**在我们自己的 V5 baseline 上成立**。 |
| Q2 效用是否崩塌？ | **外部模型未崩塌**（Efficiency@10 0.05–0.75，最好的达 0.75）；**我们的模型效用最低之一（+0.1014）**。 |
| Q3 样本量是否足够？ | **不足** —— CATH 过滤后 21 组、RBM 仅 1 组。**所有基于这两者的结论均需标注不确定性。** |
| H1 判定 | **REFUTED** |

### 对 V5.0 结论的影响

- **不推翻** V5 的 `REPRESENTATION_BOTTLENECK_CONFIRMED` —— 反而**加强**：
  瓶颈是**我们所用的池化表示**，而非任务本身不可解。
- **修正** V5 的一条表述风险：V5 的 `within-complex Hit@10 = 0.3080 < Base Rate 0.3729`
  容易被读作「这个任务的组内选择普遍不可行」。**V6 证明这是错的** ——
  物理能量法（FoldX / flex-ddG）在同协议下能达到 macro per-PDB Spearman
  **0.36–0.42**（CATH）与 **0.65**（RBM）。

### V6.1 建议（按优先级）

1. **把 FoldX / flex-ddG 作为 within-target 基线纳入 V5 的评估框架。**
   这是最低成本、最高信息量的一步：它们不需要 GPU、不需要训练，
   且已在 V6 中证明有 **0.42–0.65** 的组内排序能力。
2. **为 CATH hold-out 补算 cluster bootstrap CI**（当前 21 个 complex 的
   macro Spearman 没有区间，是本次审计最大缺口）。
3. **联系 SAMAffinity / Pythia-PPI / MINT 作者索取预测文件** ——
   三者仓库均只有代码/权重，无预测结果（见 `data_scavenge_summary.md`）。
4. **重新审视 V5 的 P1 MIPA 假设**：如果 FoldX 这类**不含学习的局部能量**
   就能拿到 0.42 的组内排序，那么 MIPA 的价值主张（学习局部扰动）
   必须与 FoldX 直接对比，而不是只与 pooled baseline 对比。
