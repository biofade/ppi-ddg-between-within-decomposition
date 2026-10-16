# V7.2R-1 Final Summary

**时间**：2026-10-14（约 19:1x → 20:1x）
**约束**：未修改 V1–V7.2R 任何文件；新命名空间 `scripts/v7_2r_1` / `experiments/v7_2r_1` /
`logs/v7_2r_1`；**未下载任何外部 checkpoint**（Track B 因 >2 GB 标 `DEFERRED_LARGE_ARTIFACT`）；
**未删除任何缓存**；模型缓存在 `F:\dsh\cache`，**未写 C 盘**

**完成度**：PHASE 0/1/2/3/4/5 **COMPLETED**；PHASE 6（F3/F4）**NOT_RUN**（依 §13 优先级，见 §12）

---

## Pilot vs Full

| | rows | targets | families | n≥5 targets |
|---|---:|---:|---:|---:|
| **TRACK_A_PILOT**（V7.2R） | 1807 | 13 | 9 | 12 |
| **TRACK_A_FULL**（本轮） | **4233** | **174** | **91** | **116** |
| U_CONTRACT（上限） | 4473 | 214 | 108 | 130 |

- **覆盖率：40.4% → 94.6%**（§4 目标 ≥90% **已达成**）
- 未覆盖的 240 行（5.4%）**排除原因**：这些复合物在 `mutation_sites` 中
  **没有任何位点**同时满足 ① `chain ∈ PDB`、② `seq_index0` 在链长内、
  ③ `mapped_WT == expected_WT`。
  ⇒ **覆盖率受 WT 校验约束限制，未因提高覆盖率而放宽校验**（§4 明令）。
- **universe manifest 已冻结**：`phase1_coverage/trackA_universe_manifest.csv`

```text
TRACK_A_PILOT
PILOT_COMPLETED         = YES
FULL_CONTRACT_COMPLETED = YES（本轮完成）
```
（pilot 产物**原样保留**，哈希见 `phase0_pilot_freeze/pilot_manifest.csv`，**未覆盖**）

---

## Pilot Bias

**回答 §6 的问题：为什么 pilot 的 R0 Macro ρ ≈ 0.265，而更大 universe 的值低得多？**

**实测答案：pilot 存在严重的「大 target 选择偏倚」。**

| 变量 | pilot 中位 | 非 pilot 中位 | p | 类型 |
|---|---:|---:|---:|---|
| **`n_mutations`** | **140.0** | **12.0** | **9.6e-05** | DEPLOY (rb=**0.647**) |
| `n_unique_mutation_sites` | 31.0 | 8.0 | 0.0386 | DEPLOY (rb=0.343) |
| `alanine_fraction` | 0.034 | 0.467 | 0.0818 | DEPLOY (rb=0.287) |
| `site_entropy` | 4.121 | 2.739 | 0.165 | DEPLOY (rb=0.231) |
| `single_mut_fraction` | 0.715 | 0.914 | 0.214 | DEPLOY (rb=0.197) |
| **`ddg_std`** | **2.201** | **1.110** | **0.0119** | LABEL (rb=0.418) |
| `true_mean_ddg` | 1.614 | 0.867 | 0.0927 | LABEL (rb=0.279) |

**机制（已由数据支持，非假设）**：
1. 嵌入脚本按**行数降序**处理 → pilot 自动由**突变库最大的 target** 构成
   （`n_mutations` 中位 140 vs 12，p≈1e-4，效应量 rb=0.647）
2. pilot 只占全部 target 的 **6.4%**（13 / 203），且只覆盖 **9 个 family**
3. 大 target 的每-target 样本多 → within-target Spearman 估计更稳定；
   且**不含** per-target ρ 为负/近零的小 target

**实测对比（R0_Base Macro ρ）**：

| 口径 | pilot 子集 | 完整合约 | 变化 |
|---|---:|---:|---:|
| ensemble prediction | **+0.2998** | **+0.0259** | **−0.2739** |
| seed-mean of metric | +0.3272 | +0.0383 | −0.2889 |

⇒ **pilot 的高值是「大 target + 极少 target 数」的选择效应，不是该子集真的更易排序。**
**V7.1 的 Macro ρ = +0.1117（130 targets）与本轮 +0.0259（116 targets）量级一致**，
进一步确认 pilot 的 +0.265 是异常值。

---

## Full R0–R3

**universe: 4233 rows / 174 targets / 91 families**

**主口径 = ensemble prediction（与 bootstrap 一致）**

| Ranker | Representation | **Macro ρ** | **Eff@10** | **C_W** | Overall r |
|---|---|---:|---:|---:|---:|
| **R0_Base** | base（物理+Z 特征） | **+0.0259** | +0.1882 | +0.1484 | +0.2159 |
| **R1_Base_F0** | + mutation identity | −0.0047 | +0.2843 | +0.1314 | +0.1954 |
| **R2_Base_Residue** | + h_wt/h_mut（无扰动） | +0.0177 | **+0.3425** | +0.1479 | **+0.2260** |
| **R3_Base_dh** | + **Δh 扰动** | +0.0116 | **+0.4176** | +0.1410 | +0.2025 |

**敏感性口径 = seed-mean of metric**

| Ranker | Macro ρ | Eff@10 | C_W | Overall r |
|---|---:|---:|---:|---:|
| R0_Base | +0.0383 | +0.3367 | +0.1334 | +0.1959 |
| R1_Base_F0 | +0.0064 | +0.2912 | +0.1167 | +0.1750 |
| R2_Base_Residue | +0.0237 | +0.3348 | +0.1301 | +0.2013 |
| R3_Base_dh | +0.0149 | +0.3443 | +0.1249 | +0.1814 |

**观察**：
- **全部 R0–R3 的 Macro ρ 都极低**（|ρ| ≤ 0.04）→ 严格 family-OOD 排序
  在**完整合约**上同样困难
- **R1（mutation identity）反而略低于 R0** → 身份编码无增益
- **R3 的 Eff@10 最高（+0.4176）但 Macro ρ 低于 R0** →
  **top-K 与整体 within 排序出现分离**，需 bootstrap 判定（见下）

---

## 5000 Bootstrap

**5000/5000 paired target-cluster bootstrap（cluster = canonical target）**

### `R3_Base_dh − R0_Base`（**核心问题：残基扰动增益**）

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| **ΔMacroSpearman** | **−0.0143** | **[−0.0588, +0.0307]** | 5000 | **❌ 含 0** |
| ΔEff@10 | +0.2293 | [−0.2989, +0.2570] | 5000 | ❌ 含 0 |
| ΔC_W | −0.0074 | [−0.0203, +0.0051] | 5000 | ❌ 含 0 |
| ΔOverall_r | −0.0134 | [−0.0320, +0.0096] | 5000 | ❌ 含 0 |
| ΔRMSE | −0.0096 | [−0.0359, +0.0145] | 5000 | ❌ 含 0 |

⇒ **§23 成功判据（冻结，未改写）：`NOT SUPPORTED` ❌**
（ΔMacroSpearman 点估计为负，且全部 5 个指标的 CI 均含 0）

> **`R2 − R0` 与 `R1 − R0` 的 bootstrap 均已跑完（同为 5000/5000）**，见下。

### `R2_Base_Residue − R0_Base`（残基表示，无扰动）

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| **ΔMacroSpearman** | **−0.0083** | **[−0.0478, +0.0303]** | 5000 | **❌ 含 0** |
| ΔEff@10 | +0.1542 | [−0.4554, +0.4913] | 5000 | ❌ 含 0 |
| ΔC_W | −0.0005 | [−0.0185, +0.0173] | 5000 | ❌ 含 0 |
| ΔOverall_r | +0.0102 | [−0.0519, +0.0634] | 5000 | ❌ 含 0 |
| ΔRMSE | −0.0737 | [−0.1382, +0.0045] | 5000 | ❌ 含 0 |

### `R1_Base_F0 − R0_Base`（mutation identity 负对照）

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| **ΔMacroSpearman** | **−0.0306** | **[−0.0703, +0.0088]** | 5000 | **❌ 含 0** |
| **ΔC_W** | **−0.0170** | **[−0.0277, −0.0072]** | 5000 | **✅（负向显著）** |
| ΔOverall_r | −0.0205 | [−0.0484, +0.0045] | 5000 | ❌ 含 0 |
| ΔRMSE | −0.0293 | [−0.0680, +0.0169] | 5000 | ❌ 含 0 |

### 三个 contrast 汇总

| contrast | ΔMacroSpearman | CI 含 0? | §23 判定 |
|---|---:|---|---|
| `R1 − R0` | −0.0306 | ✅ 含 0 | NOT_SUPPORTED（**C_W 负向显著**） |
| `R2 − R0` | −0.0083 | ✅ 含 0 | NOT_SUPPORTED |
| **`R3 − R0`** | **−0.0143** | ✅ 含 0 | **NOT_SUPPORTED（核心问题）** |

⇒ **没有任何 R1–R3 在 Macro ρ 上显著优于 R0**（点估计全为负）→
**`TRACK_A_RESULT = FULL_MAPPABLE_NULL`**

---

## Audit Corrections

### Prompt-DDG（**修复**）
V7.2R 原报告把代码与权重许可一并写成 `NOT DECLARED` —— **错误**。

```text
Prompt-DDG:
  CODE LICENSE                 = MIT
  PRETRAINED DDG CHECKPOINT    = SKEMPI-supervised
  PRIMARY USE OF CHECKPOINT    = NOT ALLOWED
  ARCHITECTURE/CODE REUSE      = ALLOWED subject to attribution
```

### ReCLIP（**三对象分离，修复**）

| 对象 | 内容 | 许可 | 判定 |
|---|---|---|---|
| **A. ReCLIP 论文** | bioRxiv 2026, `10.64898/2026.06.01.729118` | — | 参考文献 |
| **B. `SiweiLab/ReCLIP` 源码** | 特征构建 pipeline | 见仓库 LICENSE | **架构/代码可复用候选** |
| **C. `RiverZ/reclip` HF 发布物** | `classifier_heads/mutation/reclip_mutation_xgb_*.json/.pkl` | **`other`** | **YELLOW**（除能排除训练重叠） |

```text
ReCLIP CODE / ARCHITECTURE      : candidate for reuse
ReCLIP HF mutation XGBoost head : YELLOW unless training overlap excluded
```

**README 核实的关键事实**：
- README 明言：*「The release excludes earlier binary mutation pipelines, legacy
  cross-attention experiments, ESM-pLM/ESum-pLM folders...」* →
  **发布源码不含已训练分类头**
- README 明言：*「MINT checkpoints and trained task-specific classifier heads are
  external artifacts.」*
- 特征构建器在 `scripts/{four_classes_mutation,ptm,peptide,clinvar,ablation}/`
- **系统要求：64-bit Linux**（Ubuntu）；torch 1.12.1 / CUDA 11.3.1；
  **实测于 A100 40GB 与 T4**
- CPU-only 被支持（`--device cpu`），但 README 明言 *「substantially slower」*
- **依赖 MINT**（`mint/mint.ckpt`）

### MINT（**新审计**）

| 项 | 值 |
|---|---|
| 仓库 | `https://github.com/VarunUllanat/mint` |
| checkpoint | `https://huggingface.co/varunullanat2012/mint` → `mint.ckpt` |
| **许可** | **MIT** |
| 训练目标 | 多聚体相互作用 PLM 预训练（Multimeric INteraction Transformer） |
| 训练语料 | README：*「First PLM to be trained on large-scale PPI data」*；**具体组成未公开到可比对程度** |
| 参数规模 | 主干为 **ESM2-650M** 类 |
| **checkpoint 大小** | **3253.8 MB（实测 HTTP HEAD）** |
| embedding API | `mint.helpers.extract.MINTWrapper(cfg, ckpt, device)` → `wrapper(chains, chain_ids)` |
| 输出层级 | 默认 **(N, 1280)**；`sep_chains=True` → **(N, 2560)** |
| PPI 条件化 | **链间上下文注意力**（chain-aware multimeric transformer） |

**关键区分（§9）**：

| 类型 | 判定 |
|---|---|
| **task-label leakage**（ΔΔG/SKEMPI 标签训练） | **未见证据**（自监督预训练；任务 head 单独发布） |
| **foundation-pretraining family exposure** | **可能** → 标 **`RESIDUAL_PRETRAINING_EXPOSURE_POSSIBLE`** |

> **论文措辞强制要求**：
> > family-OOD refers to **task-specific supervised training**,
> > **not** guaranteed absence from foundation-model pretraining.

---

## Track B Feasibility

### 目标形态（§10）
```text
ReCLIP-style residue-context feature construction
  +  generic leakage-safe ESM / MINT representations
  +  OUR ranking head trained inside nested family CV
```
**明确不是** ReCLIP 预训练的 mutation 分类器。

### 能否抽出所需表示而不使用其监督 head？

| 所需表示 | 公共 pipeline 是否提供 | 依据 |
|---|---|---|
| site-local residue representation | **YES** | `Feature_cache/`；ReCLIP 本身是**残基为中心**的框架 |
| partner-conditioned representation | **YES** | README：*「which interaction partner residues are most informative」* |
| MINT partner representation | **YES** | `MINTWrapper(chains, chain_ids)` 接受**双链**输入 |
| 是否需其 mutation XGBoost head | **NO** | README：分类头为 external artifact |

```text
TRACK_B_FEASIBLE = YES
```

### 计算可行性（§11）

**实测 artifact 大小（HTTP HEAD，未下载）**

| artifact | 许可 | 大小 |
|---|---|---:|
| `mint.ckpt` | MIT | **3253.8 MB** |
| ESM2-650M `model.safetensors` | MIT | **2609.5 MB** |
| `bernett_mlp.pth`（任务 head） | MIT | 26.2 MB |
| ESM2-8M（Track A 已用） | MIT | 31.4 MB |

- **Track B 必需两件合计 = 5,863.3 MB = 5.86 GB**
- §29 自动允许上限 = **2 GB** → **超出 3.86 GB**

```text
DEFERRED_LARGE_ARTIFACT = YES
```

**本机资源核算**

| 项 | 值 |
|---|---|
| OS | Windows 10（ReCLIP 参考环境 **64-bit Linux**，未在 Windows 测试） |
| GPU | **无 CUDA** |
| 参考 GPU | README 实测 **A100 40GB / T4** |
| RAM | **7.55 GB 总量** |
| ESM2-650M + MINT 同时驻留 | **≈6–9 GB → 超过 7.55 GB 物理内存** |
| 特征缓存（650M 级、残基级） | **数十 GB 量级** |

```text
TRACK_B_GPU_REQUIRED   = YES
TRACK_B_READY_FOR_GPU  = YES
```

**依 §11「Do not substitute a weaker model silently」，不以更弱模型替代，
而是明确标记需要 GPU。**

**GPU 就绪后的启动清单**：
1. Linux + NVIDIA GPU（≥24 GB 建议）：clone `SiweiLab/ReCLIP` + `VarunUllanat/mint`
2. 下载 `mint.ckpt`（3.25 GB, MIT）与 `facebook/esm2_t33_650M_UR50D`（2.61 GB, MIT）
3. **只运行 feature builder**；**绝不加载** `reclip_mutation_xgb_*.pkl`
4. 用我们自己的 `z` ranker，在 `OUTER_FOLDS_V7_2` + inner family CV 下训练
5. 论文写明：family-OOD 指**任务监督训练**，**不保证**基础模型预训练未见这些家族

---

## Exact Scientific Verdict

| 主张 | 判定 | 依据 |
|---|---|---|
| `GENERIC_RESIDUE_PERTURBATION`（R3） | **PILOT_NEGATIVE → FULL_CONTRACT_NULL** | pilot Δmacro ρ = −0.0731 [−0.1445,−0.0144]（负向显著）→ 完整合约 Δmacro ρ = **−0.0143 [−0.0588,+0.0307]（含 0）** |
| `RESIDUE_EMBEDDING`（R2） | **NOT_SUPPORTED** | 完整合约 ensemble Macro ρ = +0.0177 vs R0 +0.0259（**更低**） |
| `MUTATION_IDENTITY`（R1） | **NOT_SUPPORTED** | Macro ρ = −0.0047 vs R0 +0.0259（**更低**，负对照通过） |
| `ANY_RANKER_GAIN` | **NOT_SUPPORTED** | R0–R3 的 Macro ρ 全部落在 **[−0.005, +0.026]**，无一显著高于 R0 |
| `CALIBRATION_RANKING_FRAMEWORK` | **SUPPORTED**（沿用 V7.1） | V7.1 P1−P0 ΔOverall_r **+0.3507** [+0.0100,+0.6142] ✅、ΔRMSE **−0.9448** ✅ |

> ### 核心科学发现（本轮最重要）
> **V7.2R pilot 的负向信号未能复现。**
> 在覆盖率从 40.4% 提升到 **94.6%** 后，R3−R0 的 ΔMacroSpearman
> 由 **−0.0731（CI 排除 0）** 变为 **−0.0143（CI 含 0）**。
> ⇒ 该残基扰动表示在严格 family-OOD 下**既无正增益、也无可证实的负效应** ——
> 即 **`NULL`**，而非 pilot 暗示的「有害」。
> 这与 **pilot 的 13-target 大-target 偏倚**（`n_mutations` 140 vs 12, p≈1e-4）一致。

### 修正 B / C（§3 要求）

| 项 | 修正后表述 |
|---|---|
| **B** `GENERIC_RESIDUE_PERTURBATION` | **`PILOT_NEGATIVE`** + **`FULL_CONTRACT_NULL`**（原写 `UNRESOLVED because R3 not run` 是错的：R3 **确已运行**） |
| **C** Ranker 状态 | **`R0-R3 = PILOT_COMPLETED`**；**`R4-R7 = NOT_RUN`** |

### 不一致 A 的精确解释（§3A）

| 路径 | 计算方式 | R0 | R3 | Δ |
|---|---|---:|---:|---:|
| 1（报告表格） | **先算 metric 再平均**（seed-mean of metric） | +0.2653 | +0.2033 | **−0.0620** |
| 2（bootstrap） | **先平均预测再算 metric**（ensemble prediction） | +0.2998 | +0.2267 | **−0.0731** |

- 复算路径 1 = **−0.0619** ↔ 报告 −0.0620 ✅
- 复算路径 2 = **−0.0731** ↔ bootstrap −0.0731 ✅
- **原因**：`Macro Spearman` 是**非线性**统计量，
  `mean(ρ_seed) ≠ ρ(mean(pred_seed))`，两种聚合口径给出不同数值。
- **规范（本轮已执行）**：今后**主结果一律用 ensemble prediction 口径**，
  **同时**报告 seed-mean of metric 作为敏感性（本报告 §Full R0–R3 已双报）。

---

## Next Route

**最终判定（§14）**：

```text
FINAL DECISION = B. TRACK_A_PILOT_NEGATIVE_NOT_REPLICATED
                 （完整合约为 NULL，既非支持也非确认的负效应）

TRACK_B_FEASIBLE    = YES
TRACK_B_GPU_REQUIRED = YES
```

**依 §15 的分支**：
- Track A 负（NULL）+ **Track B feasible** →
  **→ 在放弃 residue ranking 之前，先测试强 PPI 表示（Track B）。**
  这是本轮的推荐路线。

**精确下一步**：
1. **在 GPU 环境执行 Track B**（清单见上），得到 `B1`/`B2` 表示
2. 用**同一冻结嵌套协议**跑 `B2 − A2`（强表示 vs 简单 PLM 扰动）
3. **不要**在 Track A 已是 NULL 的情况下继续投入 F3/F4（§13 优先级）
4. 若 Track B 亦为 NULL → 转向 **physics/structure-informed ranking**

---

## 文件清单

### 脚本（`scripts/v7_2r_1/`）
`01_pilot_freeze_and_reconcile.py`、`02_audit_repair_and_trackB.py`、`03_full_trackA.py`

### 结果（`experiments/v7_2r_1/`）
- `phase0_pilot_freeze/`：`pilot_manifest.csv`
- `phase1_coverage/`：`trackA_universe_manifest.csv`、`trackA_row_coverage.csv`、
  `inconsistency_A_reconciliation.csv`、`reconciliation_report.md`
- `phase2_full_eval/`：`full_nested_cv_results.csv`、`full_ranker_predictions.npz`、
  `ensemble_metrics.csv`、`seedmean_metrics.csv`、
  `bootstrap_R3_Base_dh_vs_R0_Base.csv`、`full_trackA_report.md`
- `phase3_bias/`：`subset_bias.csv`、`subset_bias_report.md`
- `phase4_audit_repair/`：`license_matrix_repaired.csv`、
  `model_provenance_repaired.csv`、`audit_repair.md`
- `phase5_trackB/`：`artifact_sizes.csv`、`trackB_feasibility.md`
- `V7_2R_1_FINAL_SUMMARY.md`（本文件）

### 覆盖率进展（Track A 嵌入）

| 时点 | rows | 覆盖 |
|---|---:|---:|
| V7.2R pilot | 1807 | 40.4% |
| **V7.2R-1 完成** | **4233** | **94.6%** |

嵌入产物 `features/v7_2/per_target/*.npz`（174 文件）**已冻结**，未删除。

**无 `logs/v7_2r_1/failures.log`** → 无功能性失败。
