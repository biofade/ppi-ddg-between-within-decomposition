# V7.2R-1 PHASE 4 — 外部审计修复

- 时间: 2026-10-14 19:22:22
- 修复 V7.2R PHASE 1 的**三处错误/粗糙之处**

## 4.1 Prompt-DDG — 许可与用途拆分（**修复**）

| 对象 | 判定 |
|---|---|
| **代码许可** | **MIT** |
| **预训练 ΔΔG checkpoint** | **SKEMPI 监督** |
| **该 checkpoint 的 primary 用途** | **NOT ALLOWED** |
| **架构 / 代码复用** | **ALLOWED**（subject to attribution） |

```
Prompt-DDG:
  CODE LICENSE = MIT
  PRETRAINED DDG CHECKPOINT = SKEMPI-supervised
  PRIMARY USE OF CHECKPOINT = NOT ALLOWED
  ARCHITECTURE/CODE REUSE = ALLOWED subject to attribution
```

> V7.2R 原报告把 `license_code` 与 `license_weights` 一并写成 `NOT DECLARED`，**这是错误的**：代码为 MIT。
> 正确区分：**代码 MIT（可复用）** vs **权重 SKEMPI 监督（primary 禁用）**。

## 4.2 ReCLIP — 三个对象必须分开（**修复**）

| 对象 | 内容 | 许可 | 判定 |
|---|---|---|---|
| **A. ReCLIP 论文** | bioRxiv 2026, `10.64898/2026.06.01.729118` | — | 参考文献 |
| **B. `SiweiLab/ReCLIP` 源码** | 特征构建 pipeline + 任务脚本（**不含**训练好的 mutation head） | 见仓库 LICENSE | **架构/代码可复用候选** |
| **C. `RiverZ/reclip` HF 发布物** | `classifier_heads/mutation/reclip_mutation_xgb_*.json/.pkl` | **`other`** | **YELLOW**（除能排除训练重叠） |

```
ReCLIP CODE / ARCHITECTURE      : candidate for reuse
ReCLIP HF mutation XGBoost head : YELLOW unless training overlap excluded
```

**README 已核实的关键事实**：

- README 明言：*「The release excludes earlier binary mutation pipelines, legacy cross-attention experiments, ESM-pLM/ESum-pLM folders, and global-embedding and local automation experiment scripts.」*
  → **发布的源码不含已训练的分类头**
- README 明言：*「MINT checkpoints and trained task-specific classifier heads are external artifacts.」*
- 特征构建器：`scripts/{four_classes_mutation,ptm,peptide,clinvar,ablation}/`
- **系统要求：64-bit Linux**（Ubuntu）；参考配置 Python 3.8 / torch 1.12.1 / CUDA 11.3.1；**实测于 NVIDIA A100 40GB 与 T4**
- CPU-only **被支持**（`--device cpu` / `--xgb-device cpu`），但 README 明言「substantially slower for feature extraction」，因 ESM2 与 MINT 皆为大型 PLM
- **依赖 MINT**：`git clone https://github.com/VarunUllanat/mint.git mint`，checkpoint 放 `mint/mint.ckpt`

## 4.3 MINT 作为 foundation representation 的审计（**新增**）

| 项 | 值 |
|---|---|
| 仓库 | `https://github.com/VarunUllanat/mint` |
| checkpoint 来源 | `https://huggingface.co/varunullanat2012/mint` → `mint.ckpt` |
| **许可** | **MIT**（HF `cardData.license = mit`，tag `license:mit`） |
| 训练目标 | 掩码/对比式 PLM 预训练，**面向多聚体相互作用**（MINT = Multimeric INteraction Transformer） |
| 训练语料 | README：*「First PLM to be trained on large-scale PPI data」*；**具体 PDB/SKEMPI 组成未在 README 明确** |
| 参数规模 | 主干为 **ESM2-650M** 类（config `esm2_t33_650M_UR50D.json`） |
| **checkpoint 大小** | **3253.8 MB（≈3.25 GB）** ← **实测 HTTP HEAD** |
| embedding API | `mint.helpers.extract.MINTWrapper(cfg, ckpt, device)` → `wrapper(chains, chain_ids)` |
| 输出层级 | **默认 (N, 1280)** = 链级/对级；`sep_chains=True` → **(N, 2560)** |
| PPI 条件化机制 | 多聚体 transformer：**链间上下文注意力**（chain-aware） |

### 关键区分（§9 要求）

| 类型 | MINT 情况 |
|---|---|
| **task-label leakage**（ΔΔG/SKEMPI 标签训练） | **未见证据**（预训练为自监督 PLM；另有单独的 `bernett_mlp.pth` 任务 head，属 external artifact） |
| **foundation-pretraining family exposure**（预训练语料含我们的测试家族） | **可能**（「large-scale PPI data」组成未公开到可比对程度）→ 标 **`RESIDUAL_PRETRAINING_EXPOSURE_POSSIBLE`** |

> **论文措辞强制要求（§9）**：
> 
> > family-OOD refers to **task-specific supervised training**, **not** guaranteed absence from foundation-model pretraining.
>
> 即：我们可以声称「外层测试家族未参与**任务监督训练**」，
> **但不能**声称「基础模型预训练完全没见过这些家族」。

- 写出 `license_matrix_repaired.csv`、`model_provenance_repaired.csv`

