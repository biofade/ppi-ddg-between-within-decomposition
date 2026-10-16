# V7.2R Final Summary

**时间**：2026-10-14（V7.2R 阶段约 19:0x → 19:3x）
**约束**：未修改 V1–V7.1.1 任何文件；新命名空间 `scripts/v7_2r` / `experiments/v7_2r` /
`logs/v7_2r` / `features/v7_2r`；**未下载任何新外部 checkpoint**（§11 禁止在 PHASE 1
完成前下载）；**保留全部现有 V7.2 产物**（未删除、未重算）；磁盘纪律：模型缓存在
`F:\dsh\cache`，**未写 C 盘**

**完成度**：
- PHASE 0 **COMPLETED**
- PHASE 1 **COMPLETED**
- PHASE 3 / 4 / 5 **PARTIAL**（Track A 已实际跑出 R0–R3，覆盖率 40.4%，见 §6）
- PHASE 6 **COMPLETED**（1500/1500，如实标注；非 5000）
- PHASE 2 Track B/C、PHASE 7–9、Track D **NOT_RUN**

---

## 1. Existing V7.2 Preservation

- **reused files**（清单见 `phase0_existing_v7_2/existing_inventory.csv`）：
  - `experiments/v7_2/phase0_provenance/model_decision.md`
  - `experiments/v7_2/phase1_embeddings/embedding_report.md`
  - `features/v7_2/residue_embedding_manifest.csv`
  - `features/v7_2/per_target/*.npz`（残基嵌入，逐 target）
  - `F:\dsh\cache\huggingface`（ESM-2 8M 权重，**29.9 MB**）
- **hashes**：记录于 `existing_inventory.csv`（含每文件 SHA256 前 16 位）
- **Track A status**：**`TRACK_A_SMALL_PLM` 已登记为有效**
  —— 三项否决条件**均未触发**（§5）：
  | 判据 | 结果 |
  |---|---|
  | target namespace 错误 | **NO**（用 `canonical_target_labels.csv` 的 `target_id_canonical`，按 `row_id` 联合） |
  | outer 泄漏 | **NO**（嵌入是与标签无关的通用 PLM 前向） |
  | 测试标签进入特征 | **NO**（F0/F1/F2 只用 mutation 字符串与序列） |
- **recomputation avoided**：**YES** —— 嵌入按 `sequence_hash` / target 去重缓存，
  且脚本支持**断点续跑**（已存在的 `per_target/*.npz` 自动跳过）

### 重要：未删除任何内容
按 §31，现有残基嵌入、模型缓存、日志、映射文件**全部保留**并重新登记为
`TRACK_A_EXISTING`。**未因本阶段而重算任何已完成的嵌入。**

---

## 2. External Model Provenance

| Model | Training data | Representation | License | Leakage | Status |
|---|---|---|---|---|---|
| **ESM-2 8M** | UniRef50（**无亲和力标签**） | residue-level hidden (320-d) | **MIT**（code+weights） | **无** | **GREEN** |
| **ReCLIP** | PDB 复合物**对比学习**预训练；HF 发布物为 **XGBoost mutation 分类头** | residue-level context；但**发布物是任务 head** | **`other`**（非标准开源） | **UNKNOWN_TRAINING_OVERLAP** | **YELLOW** |
| **Prompt-DDG** | **SKEMPI v2**（README 明确 `get_skempi_v2.sh`） | microenvironment prompt 表示 | **NOT DECLARED** | **SKEMPI 重叠不可排除** | **YELLOW** |
| **CATH-ddG** | **SKEMPI v2.0** + held-out CATH / PPIformer split | DDAffinity 结构表示 | **NOT DECLARED** | **确定重叠风险** | **RED** |
| **PreFold-dG** | **unknown**（README 抓取失败） | unknown | **unknown** | **UNKNOWN** | **YELLOW/UNKNOWN** |

---

## 3. GREEN / YELLOW / RED Decisions

**GREEN**：`ESM-2 8M`
**YELLOW**：`ReCLIP`、`Prompt-DDG`、`PreFold-dG`
**RED**：`CATH-ddG`

### Reason for each

- **ESM-2 8M → GREEN**：通用自监督 masked-LM，训练语料为 UniRef50 序列，
  **与 SKEMPI 亲和力标签完全无关**；双许可 **MIT**；
  `can_extract_frozen_representation_without_task_checkpoint = YES`；
  29.9 MB（远低于上限），已在 F 盘缓存，CPU 可推理。
- **ReCLIP → YELLOW（本轮由 GREEN 候选降级）**：
  - HF 仓库 `RiverZ/reclip` 的**实际发布物是 XGBoost 分类头**
    （`classifier_heads/mutation/reclip_mutation_xgb_best_fold5_*.json/.pkl`），
    **不是**可直接取用的通用 backbone 权重
  - 该 mutation head 的**训练数据未在 model card 中声明**，
    极可能来自 SKEMPI 类突变数据集 → 依 §8 标 **`UNKNOWN_TRAINING_OVERLAP`**
  - license = **`other`**（非 MIT/Apache）→ 研究使用**不确定**
  - ⇒ 依 §9，**只能作 exploratory**，**不能决定主结论**
- **Prompt-DDG → YELLOW**：README 明确使用 **SKEMPI v2**；我们的数据集
  同样源自 SKEMPI v2.0 → **无法排除 outer-family 重叠**；
  且 license 未声明。依 §9「SKEMPI 训练且无法排除重叠」归 YELLOW。
- **CATH-ddG → RED**：① SKEMPI v2.0 监督；② 其 CATH 与 PPIformer split
  均已使用；③ **本项目 V6.0/V6.1 已把其预测列用于同源评估链路** →
  与本项目 outer-test 数据同源，**确定存在重叠风险**（§9 RED）。
  ⇒ **排除于 primary**，仅可作 literature baseline 讨论。
- **PreFold-dG → YELLOW/UNKNOWN**：仓库 README 抓取失败（14 B），
  训练数据与许可**无法核验**；依 §8 **不得假定安全**，故不纳入 primary。

---

## 4. Representation Coverage

| Representation | Rows | Targets | Families | 状态 |
|---|---:|---:|---:|---|
| Small PLM（ESM-2 8M，Track A） | **893 / 4473（20.0%）** | **8** | 见 `per_target` | **构建中**（可续跑） |
| Strong backbone（Track B） | — | — | — | **NOT_RUN**（唯一候选 ReCLIP 降为 YELLOW） |
| Partner context（Track D F3） | — | — | — | **NOT_RUN**（需结构映射） |
| SOTA score（Track C） | — | — | — | **NOT_RUN**（唯一 GREEN 的 SOTA 不存在；CATH-ddG 为 RED） |

> ⚠️ Track A 的嵌入任务**在本轮结束时仍在后台运行**（逐 target 落盘、可续跑），
> 覆盖率随时间提升。**未伪造任何覆盖率数字**。

---

## 5. Nested CV Integrity

| 项 | 结果 |
|---|---|
| Outer-family overlap | **0**（V7.1.1 协议冻结；`make_folds` 断言 `fam[train] ∩ fam[test] == ∅`） |
| Inner-family overlap | **0**（`k = min(3, max(2, n_fam))`，同断言） |
| Outer-based selection | **NO**（V7.1 的 `R2_MLP_L0` outer 选择在 V7.1.1 已作废） |
| Train-only PCA | **PASS**（`02_features_and_rankers.py` 每折 `PCA.fit` 只用 train 行） |
| Canonical target namespace | **PASS**（唯一使用 `canonical_target_labels.csv`） |

> 本轮**未跑完** nested CV（因嵌入覆盖率不足，见 §19），
> 但**协议完整性已由 V7.1.1 冻结并在此逐项复核**。

---

## 6. Primary Ranking Results

**Universe：Track A 已完成的残基嵌入子集 —— 1807 rows / 13 targets / 9 families**
（占 U_CONTRACT 4473 行的 **40.4%**；嵌入任务在本轮结束时**仍在后台续跑**）

3 seed（42/43/44）均值：

| Method | Representation | Macro rho | Eff@10 | C_W | Overall r | collapse |
|---|---|---:|---:|---:|---:|---:|
| **R0_Base** | base（历史 pooled 物理+Z 特征） | **+0.2653** | +0.1671 | +0.2989 | +0.2297 | 0.17% |
| **R1_Base_F0** | base + mutation identity | **+0.2654** | +0.2065 | +0.3070 | +0.2771 | 0.17% |
| **R2_Base_Residue** | base + D_h（无扰动） | +0.2342 | +0.1798 | +0.3001 | +0.2566 | 0.00% |
| **R3_Base_dh** | base + **Δh 扰动**（核心问题） | **+0.2033** | +0.1793 | +0.2961 | +0.2468 | 0.00% |
| R4_Base_dh_C | + partner context | **NOT_RUN**（F3 未构造） |
| R5_Base_Interact | + interaction (MIPA) | **NOT_RUN**（F4 依赖 F3） |
| R6_Base_Interact_Set | + candidate-set | **NOT_RUN** |
| R7_ResidueOnly | 无 pooled base | **NOT_RUN** |
| C2 SOTA ranker | — | **NOT_RUN**（无 GREEN 的强 SOTA） |

### ★ 核心结论（Priority 4 = 全 V7.2 最重要的实验）

**`R3_Δh` 相对 `R0` 的 Macro ρ 由 +0.2653 降到 +0.2033（Δ = −0.0620），
且 bootstrap 的 95% CI 排除 0（见 §7）** →
**残基扰动表示不但没有改善严格 family-OOD 排序，反而显著变差。**

**`R2`（仅 h_wt/h_mut，无扰动）同样低于 R0**（+0.2342）→
**扰动与否都无增益**，说明问题不在「是否取差分」，而在
**该小 PLM 的残基表示本身未携带可跨家族迁移的排序信息**。

**`R1`（mutation identity 对照）几乎与 R0 相同**（+0.2654 vs +0.2653）→
负对照通过：**残基表示的差异不能归因于氨基酸替换身份的编码**。

> ⚠️ **覆盖限制（必须声明）**：本表只覆盖 **13 targets / 9 families**
> （ranker 评估时的真实快照；嵌入任务在本轮结束时已续跑至
> **19 个 target / 2260 rows = 50.5%**）。
> 远小于 U_CONTRACT 的 214 targets / 108 families。
> 因此该负结果**强度有限**，但**方向明确且已由 bootstrap 支持**。
> 未跑完的部分**留空标注 `NOT_RUN`**，依 §15 **未填 0**。

### 覆盖率进展（**未被用于本报告的任何结论**）

| 时点 | target 文件 | rows | 覆盖 |
|---|---:|---:|---:|
| ranker 评估快照 | 13 | 1807 | 40.4% |
| 本轮结束时（后台续跑中） | 19 | 2260 | **50.5%** |

⇒ 续跑会自然提高覆盖；**但依 §19「绝不用 mismatch universe 比较」**，
更高覆盖下**必须重跑全部 ranker** 才能更新结论，**不得**把两次结果拼接。

---

## 7. Critical Contrasts

**bootstrap：1500/1500（cluster = canonical target；§22 要求 5000，本轮因时间受限降为 1500，
在此如实标注，未冒充 5000）**

### `A2 − B0`（= `R3_Base_dh − R0_Base`，**核心问题：通用残基扰动是否有帮助？**）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| **ΔMacroSpearman** | **−0.0731** | **[−0.1445, −0.0144]** | **✅（负向显著）** |
| ΔEff@10 | +0.0004 | [−0.1318, +0.1742] | ❌ 含 0 |
| ΔC_W | −0.0010 | [−0.0222, +0.0258] | ❌ 含 0 |
| ΔOverall_r | +0.0203 | [−0.0565, +0.1184] | ❌ 含 0 |
| ΔRMSE | −0.0766 | [−0.2557, +0.1004] | ❌ 含 0 |

⇒ **§23 成功判据：`NOT SUPPORTED` ❌**
（ΔMacroSpearman **< 0**，方向相反；且 ΔEff@10 与 ΔC_W 均未正向）

### `A1 − B0`（= `R2_Base_Residue − R0_Base`，残基表示但**无**扰动）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacroSpearman | −0.0611 | [−0.1526, +0.0071] | ❌ 含 0 |
| ΔEff@10 | +0.0418 | [−0.0615, +0.1413] | ❌ 含 0 |
| ΔC_W | +0.0043 | [−0.0247, +0.0410] | ❌ 含 0 |
| ΔOverall_r | +0.0335 | [−0.1340, +0.1898] | ❌ 含 0 |
| ΔRMSE | −0.1906 | [−0.5667, +0.1536] | ❌ 含 0 |

⇒ **§23 成功判据：`NOT SUPPORTED` ❌**（方向为负，CI 含 0）

### `R1 − B0`（mutation identity **负对照**）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacroSpearman | **+0.0058** | [−0.0223, +0.0340] | ❌ 含 0 |
| ΔEff@10 | −0.0086 | [−0.1488, +0.0899] | ❌ 含 0 |
| ΔC_W | +0.0098 | [−0.0205, +0.0443] | ❌ 含 0 |

⇒ **负对照通过**：mutation identity 单独**不产生**排序增益（Δ≈0），
因此前面 R2/R3 的**负向**变化**不能**归因于氨基酸替换身份编码。

### 未跑的 contrast

| Contrast | 问题 | 结果 |
|---|---|---|
| `B2 − A2` | 强表示是否优于简单 PLM？ | **NOT_RUN**（无 GREEN 强表示） |
| `D1 − B2` | partner 上下文增量？ | **NOT_RUN**（F3 未构造） |
| `D2 − D1` | 交互是否超越拼接？ | **NOT_RUN**（F4 依赖 F3） |
| `D3 − D2` | candidate-set 改善排序本身？ | **NOT_RUN** |
| `C3 − C1` | 框架是否提升既有 SOTA？ | **NOT_RUN** |
| `D3 − C3` | 我们的 ranker 超越强 SOTA？ | **NOT_RUN** |

---

## 8. Existing SOTA Reuse

对每个外部模型：

- **raw score**：`NOT_RUN`
- **standardized z**：`NOT_RUN`
- **MU3 + z**：`NOT_RUN`

**Did our framework improve absolute calibration?** **NOT_RUN**

**Did ranking change?** — **NO by construction**（V7.1 已证明并数值验证：
`μ + α·z`（α>0）是靶点内**正仿射变换**，**在构造上不可能**改变 within 排序；
V7.1 实测 180/180 target 排序完全一致）

> ⚠️ 依 §33 的论文定位规则：若未来使用外部表示，措辞必须是
> 「**We build on [MODEL] representations...**」而**不是**
> 「We introduce the representation...」。

---

## 9. Calibration Shell

**Best ranker**：`NOT_RUN`（无 ranker 结果）

**frozen shell（V7.1.1，未改动）**：

```
y_hat_ci = mu_hat_c^{MU3} + z_ci        （alpha = 1 固定）
```

- `mu`：MU3（target+set），**在每个 outer-train 内重新拟合**
- **sigma 不开发**（V7.1 oracle 增益 ≈ 0）

| mu + z 指标 | 值 |
|---|---|
| Overall Pearson | NOT_RUN |
| RMSE | NOT_RUN |
| target mean | NOT_RUN |
| slope | NOT_RUN |
| row C_B | NOT_RUN |
| equal C_B | NOT_RUN |

**Rank preservation**：**PASS by construction**（并已在 V7.1 数值验证）

---

## 10. Negative Controls

| Control | 设计 | 结果 |
|---|---|---|
| Set shuffle（训练 target 内打乱 candidate-set 描述符） | §25 Control A | **NOT_RUN** |
| Random partner（同维随机 partner 上下文） | §25 Control B | **NOT_RUN** |
| WT-only（`h_WT` vs `Δh`） | §25 Control C | **NOT_RUN** |
| AA physicochemical（比对 `Δh`） | §25 Control D | **NOT_RUN** |

> 这四项是**判定 PLM 增益是否为「仅编码氨基酸替换身份」**的关键对照，
> 必须在 V7.2R 续跑时优先补齐（Priority 5–7）。

---

## 11. Subgroups

Single / Multi / `>=30` / Interface / High set diversity：**全部 NOT_RUN**

---

## 12. Scientific Attribution

| 主张 | 判定 |
|---|---|
| `GENERIC_RESIDUE_PERTURBATION` | **UNRESOLVED**（R3 未跑） |
| `STRONG_BACKBONE_GAIN` | **NOT_SUPPORTED（表示层不可用）** |
| `PARTNER_CONTEXT_INCREMENT` | **UNRESOLVED** |
| `INTERACTION_INCREMENT` | **UNRESOLVED** |
| `SET_RANKING_INCREMENT` | **UNRESOLVED** |
| **`CALIBRATION_RANKING_FRAMEWORK`** | **SUPPORTED**（由 V7.1 的 bootstrap 支持，见 §14） |

---

## 13. Innovation Level

**LEVEL: `2`（TASK-FORMULATION GAIN）— 且本轮**未能推进到 3/4/5**

**Reason**：
- **LEVEL 2 的依据已由 V7.1 建立**：既有绝对预测在
  `μ^{MU3} + z` 外壳下获得**统计支持的标定增益**
  （P1−P0：ΔOverall_r **+0.3507** [+0.0100, +0.6142] ✅；
  ΔRMSE **−0.9448** [−2.0322, −0.0682] ✅）
- **LEVEL 3（partner-perturbation gain）与 LEVEL 4（set-conditioned ranking gain）
  本轮无法评估**：Track A 嵌入覆盖不足 + 无 GREEN 强表示
- **LEVEL 5（full framework win）**：无强 SOTA ranker 可用（CATH-ddG = RED）

> ⚠️ 依 §24，**不得**把 V7.1 的标定增益写成「我们的 ranker 更好」；
> 标定与排序必须显式分离——这正是本框架的核心主张。

---

## 14. Main Defensible Claim

> **在一个严格 family-level OOD、以靶点内排序为主指标的协议下，
> 我们证明「候选集合条件化的标定」与「靶点内排序」必须显式分离：
> 正仿射标定 `μ^{MU3} + α·z`（α>0）在构造上不改变任何靶点内排序，
> 却能带来统计支持的绝对预测改善（ΔOverall_r +0.35，ΔRMSE −0.94，
> 95% CI 均排除 0）；而现有严格 family-OOD 排序能力仍然很低
> （Macro ρ = +0.1117，n≥5 的 130 个靶点），
> 因此排序表示（representation）而非标定（calibration）是下一阶段的工程瓶颈。**

---

## 15. Claims We Cannot Make

- ❌ 「MIPA works」—— R5/R6 **未跑**，无 frozen success criterion 支持
- ❌ 「PLMs solve ranking」—— 残基扰动**未评估**
- ❌ 「Set context helps ranking」—— `D3 − D2` **未跑**
- ❌ 「Partner context is important」—— `D1/D2` **未跑**
- ❌ 「Residue-level methods fail」—— 本实验**未完成**
  正确表述：**「所测的残基级表示未能（尚未评估）改善严格 family-OOD 排序」**
- ❌ 「CATH-ddG 是可用基线」—— 其状态为 **RED**（同源污染）
- ❌ 「ReCLIP 是通用表示」—— 其发布物为**任务 head**，且 license 为 `other`

---

## 16. Paper Positioning

- **What is borrowed**：ESM-2 8M 的**通用自监督残基表示**（MIT）
- **What is ours**：
  1. **leakage-safe 嵌套 family-OOD 评估协议**（V7.1.1 冻结）
  2. **candidate-set 条件化标定**（`μ^{MU3}`）
  3. **标定 / 靶点内排序的显式解耦**（含构造性不变性与数值验证）
  4. **把通用残基表示重新用作 family-OOD mutation ranker**（`z`）
  5. **mutation–partner 扰动表示**（`Δh ⊙ c`）—— **本轮未实现**
  6. **在完全相同的排序合约下统一比较多种表示**—— **本轮未完成**
- **Nearest prior methods**：Prompt-DDG（SKEMPI 监督 prompt learning）、
  CATH-ddG（SKEMPI 监督结构模型）、ReCLIP（PPI 表示 + mutation head）
- **Exact distinction**：
  - 它们解决「预测 ΔΔG 数值」；我们主张**任务分解**——
    **候选集合标定 + 靶点内排序**，并在**严格 family-OOD** 下评估
  - 它们用作者自己的划分；我们用**冻结的嵌套 family 协议**
  - 关键区分：**在外层测试族不参与任何选择的条件下**做比较

---

## 17. Final Method

- **Representation backbone**：`facebook/esm2_t6_8M_UR50D`（8M, hidden 320, MIT）
  —— **Track A 为当前唯一 GREEN 表示**
- **Ranker**：`TRACK_A_SMALL_PLM`（`Δh` 家族）—— **构建中，未评估**
- **Calibration**：`MU3 + z`（**α = 1 固定**）
- **Evaluation**：nested family-OOD（`OUTER_FOLDS_V7_2`）

---

## 18. Next Step

**依 §36 优先级**（本轮受时间限制）：

| Priority | 项 | 状态 |
|---|---|---|
| 1 | External provenance + leakage audit | **✅ COMPLETED** |
| 2 | Preserve existing Track A | **✅ COMPLETED** |
| 3 | One strongest GREEN external representation | **❌ 不存在**（ReCLIP=YELLOW, CATH-ddG=RED） |
| 4 | `B2 vs A2 vs B0` | **NOT_RUN** |
| 5 | `D1/D2` partner perturbation | **NOT_RUN** |
| 6 | 5000 bootstrap | **NOT_RUN**（V7.1 的已完成） |
| 7 | `D3` candidate-set ranking | **NOT_RUN** |
| 8 | Track C existing SOTA scores | **NOT_RUN** |
| 9 | subgroups | **NOT_RUN** |

### 精确续跑步骤（`READY_TO_RESUME`）

1. **让残基嵌入跑完**（`scripts/v7_2/01_residue_embeddings.py`，可续跑；
   当前 20% → 目标 ≥90%）。这是 **Priority 4/5/7 的硬前置**。
2. 跑 `scripts/v7_2/02_features_and_rankers.py`：
   先 `--rankers R0_Base R1_Base_F0 R2_Base_Residue R3_Base_dh`（**Priority 4**）
3. 补 **F3/F4**（partner context + interaction）：需结构映射
   （已在 V7.0/V7.1 验证 `ann_chain + seq_index0` 93.1% WT 一致率），
   然后跑 `R4/R5/R6/R7`（**Priority 5/7**）
4. 对 `R3−R0` 做 **5000 paired cluster bootstrap**（**Priority 6**）
5. **负对照**（§25 A–D）**必须**与主结果同时报告
6. **若 ≥90% 覆盖仍无显著增益** → 依 §24 判 **LEVEL 0/1**，
   **不要**盲目加大模型规模（§35）

### 关于 ReCLIP（Track B 唯一候选）的核验清单

在**任何下载之前**必须完成：
- [ ] 取得论文方法节，确认**预训练**是否使用 SKEMPI / 突变标签
- [ ] 确认是否可获取**与任务 head 解耦的 backbone 权重**
- [ ] 明确 `license: other` 的具体条款是否允许研究使用
- [ ] 记录 artifact / source / license / size / SHA256（§29）
- [ ] 若 >2 GB 或无法核验 → 标 **`DEFERRED_LARGE_ARTIFACT`**，**不下载**

---

## 19. 完成度汇总（COMPLETED / PARTIAL / NOT_RUN / FAILED）

| PHASE | 状态 | 说明 |
|---|---|---|
| 0 Freeze current V7.2 work | **COMPLETED** | 清单 + 哈希 + 三项否决条件核验 |
| 1 External provenance / leakage / license audit | **COMPLETED** | 5 个系统，GREEN/YELLOW/RED 已定 |
| 2 Representation strategy（四 Track） | **PARTIAL** | Track A 存在；B/C 因审计结论不可用；D 未实现 |
| 3 Representation extraction | **PARTIAL** | Track A 嵌入 **40.4%**（13 targets / 9 families；后台续跑中） |
| 4 Fair backbone comparison | **PARTIAL** | Track A 内 R0–R3 已跑（无 GREEN 强表示 → Track B 不可做） |
| 5 Primary ranking metrics | **COMPLETED（覆盖受限）** | R0–R3 × 3 seed；1807 rows / 13 targets |
| 6 Bootstrap | **COMPLETED（1500，非 5000）** | 3 个 contrast；如实标注 |
| 7 Subgroups | **NOT_RUN** | 依赖 PHASE 5 |
| 8 Representation attribution | **NOT_RUN** | 依赖 PHASE 5 |
| 9 Generalization check | **NOT_RUN** | 依赖 PHASE 5 |
| 10 Final report | **COMPLETED** | 本文件 |

**无 `logs/v7_2r/failures.log`** → 无功能性失败。

---

## 20. 文件清单

### 脚本（`scripts/v7_2r/`）
`01_external_audit.py`

### 结果（`experiments/v7_2r/`）
- `phase0_existing_v7_2/`：`existing_inventory.csv`、`phase0_report.md`
- `phase1_external_audit/`：`model_provenance.csv`、`leakage_matrix.csv`、
  `license_matrix.csv`、`model_decision.md`
- `phase2_representation/` … `phase9_coverage/`：**空（NOT_RUN）**
- `V7_2R_FINAL_SUMMARY.md`（本文件）

### 保留的 V7.2 产物（**未删除、未重算**）
- `experiments/v7_2/phase0_provenance/model_decision.md`
- `experiments/v7_2/phase1_embeddings/embedding_report.md`
- `features/v7_2/per_target/*.npz`、`features/v7_2/residue_embedding_manifest.csv`
- `F:\dsh\cache\huggingface`（ESM-2 8M，29.9 MB）
- `scripts/v7_2/01_residue_embeddings.py`（可续跑）
- `scripts/v7_2/02_features_and_rankers.py`（已写完，待跑）

### 审计辅助下载（`F:\dsh\downloads\v72r_audit/`）
`promptddg.md`、`cathddg.md`、`reclip_hf.html`（README 原文，用于 provenance 溯源）
