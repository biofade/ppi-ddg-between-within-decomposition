# V7.2R PHASE 1 — 外部模型 Provenance / Leakage / License 审计

- 时间: 2026-10-14 18:34:00
- **本 PHASE 完成前不下载任何外部 checkpoint**（§11）
- 候选系统数：**5**（≤5，符合 §6）

## 1.1 审计主表

| Model | Representation | Checkpoint | SKEMPI 重叠 | License | Status | Allowed use |
|---|---|---|---|---|---|---|
| **ESM-2 8M (small PLM)** | generic protein language model (masked LM) | YES | NO | MIT | GREEN | primary: Track A baseline + Δh perturbation |
| **Prompt-DDG** | ΔΔG prediction (supervised) | license 未声明 | **YES** | NOT DECLARED | **YELLOW** | EXPLORATORY_ONLY（附录/次要）；**不能**决定主结论 |
| **CATH-ddG (DDAffinity)** | ΔΔG prediction (supervised) | license 未声明 | **YES** | NOT DECLARED | **RED（对本项目的 primary 评估）** | **排除**于 primary；仅可作 literature baseline 讨论 |
| **ReCLIP** | contrastive pretraining of residue-level PPI representations | **license=other，不确定** | **uncertain**（mutation 分类头的训练数据未在 model card 声明） | **other** | **YELLOW（本轮由 GREEN 候选降级）** | **EXPLORATORY_ONLY**；不得进入 primary |
| **PreFold-dG** | ΔΔG prediction from folding-model intermediate representations | unknown | likely YES（同类任务） | unknown | **YELLOW / UNKNOWN** | **不纳入 primary**（文档不足，无法完成 §7 审计） |

## 1.2 GREEN / YELLOW / RED

- **GREEN**：['ESM-2 8M (small PLM)']
- **YELLOW**：['Prompt-DDG', 'ReCLIP', 'PreFold-dG']
- **RED**：['CATH-ddG (DDAffinity)']

### ESM-2 8M (small PLM) → GREEN

- **理由**：与 SKEMPI 标签完全无关的通用自监督预训练；MIT 双许可
- 训练数据：UniRef50 (sequences, ~2022)
- `can_extract_frozen_representation_without_task_checkpoint?`：**YES**

### Prompt-DDG → **YELLOW**

- **理由**：SKEMPI v2 监督训练；我们的数据集同样源自 SKEMPI v2.0 → 无法排除 outer-family 重叠（§8 `UNKNOWN_TRAINING_OVERLAP`）；且 license 未声明
- 训练数据：**SKEMPI v2**（README 明确 `data/get_skempi_v2.sh`）
- `can_extract_frozen_representation_without_task_checkpoint?`：**UNCERTAIN**

### CATH-ddG (DDAffinity) → **RED（对本项目的 primary 评估）**

- **理由**：① SKEMPI v2.0 监督；② 我们 V6.0/V6.1 已把其预测列用于同源评估链路 → 本项目的 outer-test 与其训练数据同源，**确定存在重叠风险**（§9 RED）
- 训练数据：**SKEMPI v2.0**（README 明确 `skempi_parallel.py --subset skempi_v2`）
- `can_extract_frozen_representation_without_task_checkpoint?`：**UNCERTAIN**

### ReCLIP → **YELLOW（本轮由 GREEN 候选降级）**

- **理由**：① HF 实际发布物是 **XGBoost mutation 分类头**而非通用 backbone；② 该 head 训练数据未声明 → 依 §8 标 `UNKNOWN_TRAINING_OVERLAP`；③ license = `other`（非 MIT/Apache）→ 研究使用不确定。三条独立理由使其**无法满足 §7 末条的 frozen-representation 条件**
- 训练数据：PDB 复合物结构（对比学习）；**HF 发布物为 XGBoost mutation 分类头**
- `can_extract_frozen_representation_without_task_checkpoint?`：****NO**（仓库只提供任务 head，未提供可解耦 backbone）**

### PreFold-dG → **YELLOW / UNKNOWN**

- **理由**：仓库 README 不可得；训练数据与许可无法核验 → 依 §8 `UNKNOWN_TRAINING_OVERLAP`，不得假定安全
- 训练数据：unknown
- `can_extract_frozen_representation_without_task_checkpoint?`：**UNKNOWN**

## 1.3 关键结论

**GREEN（可用于 primary）**：只有 **ESM-2 8M**（通用自监督，MIT，已下载 29.9 MB）

**ReCLIP 是 Track B 的首选候选**，但：
- README 抓取失败（14 B）→ 其**预训练是否使用 SKEMPI 标签无法确认**
- license 未声明
⇒ 依 §8，**标 `UNKNOWN_TRAINING_OVERLAP`**，
  在完成核验（HF 文件清单 + 论文方法节）之前**不下载、不纳入 primary**。

**Prompt-DDG**：README 明确使用 **SKEMPI v2**，且 license 未声明 → **YELLOW**（exploratory only）

**CATH-ddG**：SKEMPI v2.0 监督 + CATH/PPIformer split 均已使用；
且**我们 V6.0/V6.1 已把它的预测列用于同源评估链路** →
**RED（对本项目 primary 评估）**，仅可作 literature baseline 讨论。

## 1.4 依 §37 的 STOP/GO 判定

- **存在 1 个 GREEN 表示**（ESM-2 8M，且**已在 V7.2 完成部分嵌入**）
- ReCLIP 待核验；未核验前**不下载**
⇒ **结论：继续 Track A + 我们的 perturbation 模型**（§37 第一分支），
  并把 ReCLIP 的核验列为 Track B 的**前置条件**而非阻塞项。

