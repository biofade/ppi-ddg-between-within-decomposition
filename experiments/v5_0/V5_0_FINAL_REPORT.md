# V5.0 FINAL REPORT — Affinity Intelligence

**生成时间**：2026-10-13
**执行环境**：Windows 10 / Python 3.11.9 / Intel i7-2630QM（4 核 8 线程）/ 7.55 GB RAM / **无 GPU** / torch 2.5.1+cpu / LightGBM 4.7.0 / numpy 2.4.6 / pandas 3.0.6 / scipy 1.17.1 / sklearn 1.9.1 / Biopython 1.88 / **非 git 仓库（`GIT_NOT_AVAILABLE`）**
**证据分级**：`PROJECT_EVIDENCE`（本仓库实测）/ `LITERATURE_EVIDENCE`（论文报告）/ `HYPOTHESIS`（未验证推断）—— 全文逐条标注。

---

## Executive Summary

**一句话回答：V5 没有突破 V4。V5 的最重要产出是「用最严格的检验证明了 V4 的数字不可用于真实选靶场景」。**

| 指标 | V4 最强（3-seed 均值） | V5 最强（3-seed 均值） | Δ |
|---|---:|---:|---:|
| Overall Pearson | **+0.3861** | +0.3861（无变体超越） | **+0.0000** |
| Overall Spearman | +0.1871 | +0.1871 | +0.0000 |
| macro per-PDB Pearson | +0.1334 | +0.1334 | +0.0000 |
| **macro per-PDB Spearman** | **+0.0956** | +0.0956 | **+0.0000** |
| **Global P@10** | **0.9667** | 0.9667 | +0.0000 |
| **Global EF@1%** | **9.6034** | 9.6034 | +0.0000 |
| **Within-complex Efficiency@10（V5 新增）** | **+0.1014** | +0.1014 | +0.0000 |
| **Within-complex Hit@10(ΔΔG<0)** | **0.3080** | 0.3080 | +0.0000 |

> 注：`Efficiency@10` 有两个口径，**都如实报告**：
> `v5_common.topk_metrics`（仅 ≥5 成员的 complex，n=59）→ **+0.1014**；
> `14_p4_decompose_topk.py`（含 5 成员以下的 complex，n 更多）→ **+0.2252**。
> 两者结论一致：**模型只回收约 10–23% 的 oracle 可用空间。**

### V5 的三个决定性新读数（全部 `PROJECT_EVIDENCE`）

1. **真实独立样本量是 61，不是 1129。** cluster bootstrap（抽样单位 = 复合物）下，
   `macro per-PDB Spearman = +0.0956`，**95% CI = [−0.0219, +0.2385]（含 0）**。
   → **「同一复合物内挑突变体」这一核心能力在统计上无法与随机区分。**

2. **全局 Top-K 的高分是假象。** 全局 P@10 = 0.9667 看起来很成功，但
   **复合物内**选 top-10 的真实 ΔΔG 均值 = **+0.3577**，
   **随机** = **+0.6277**，**oracle** = **−0.5709**
   → 模型相对随机**确有** +0.2700 kcal/mol 的改善，但**只回收了 oracle 可用空间的
   22.5%**（Efficiency@10 = **+0.2252**；保守口径 **+0.1014**）。
   更关键的是：within 层「命中 ΔΔG<0」仅 **0.3080**，**低于全局基准率 0.3729**。
   → **在一个靶点内挑 10 个「会增强结合」的突变体，命中率不如随机。**
   全局优势 **−7.2490 kcal/mol**（vs 全体均值）几乎全部来自「区分不同复合物」，
   而非组内排序。

### 最终判定（V5 §128 允许的五个标签之一）

> # `REPRESENTATION_BOTTLENECK_CONFIRMED`

依据：所有训练/目标层面的技巧（state potential、加权、bagging、conformal、
正反向增广以外的新机制）均失败，但 MIPA 方向被**清晰诊断为受阻于数据可映射性
与算力**，且**第 2 条新读数把瓶颈精确定位到「池化表示丢失位点信息」**。

---

## Baseline Reproduction

### GATE A — **PASS ✅（逐位复现）**

| seed | V4 记录 Pearson | V5 复现 Pearson | Δ | 判定 |
|---:|---:|---:|---:|---|
| 42 | +0.3951 | +0.3951 | −0.0000 | ✅ |
| 43 | +0.3843 | +0.3843 | −0.0000 | ✅ |
| 44 | +0.3790 | +0.3790 | +0.0000 | ✅ |

**结论：无 baseline mismatch。** 按 V5 §11，`|ΔPearson| > 0.01` 才需查明原因；
实测 `|Δ| ≤ 1e-6`，因此**不需要** `legacy_reproduction` 与
`current_environment_reproduction` 双轨。

- 配置：`C+DGaux_fwd+rev`（V4 创新 2）；LightGBM **固定 800 树**（V4 已证早停在 iteration 1 停）；
  特征 `D` = [物理 diff 60 + pooled ESM-2 PCA128] = **188** 维；
  辅助任务用 `GroupKFold(pdb_id)` 5 折 OOF ΔG_mut 作 stacking 列；
  正向 + 反向增广（反向 y 取负，辅助列不取负）。
- 辅助任务 OOF Pearson 均值 **+0.4368**（5 折 +0.3216/+0.5823/+0.2887/+0.5347/+0.4567），
  test 上 **+0.4684**。
- 脚本：`scripts/v5/02_baseline_repro.py`；产物 `experiments/v5_0/p0_benchmark/`。

### V5 为 V4 补齐的指标（V4 未报告）

| seed | Spearman | RMSE | MAE | AUROC | macro-PDB ρ | macro-PDB r | Hit@10(ΔΔG<0) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 42 | +0.1972 | 2.1357 | 1.5249 | 0.5981 | +0.1048 | +0.1307 | 0.9000 |
| 43 | +0.1997 | 2.1529 | 1.5364 | 0.6046 | +0.0905 | +0.1257 | 1.0000 |
| 44 | +0.1640 | 2.1566 | 1.5476 | 0.5778 | +0.0914 | +0.1438 | 1.0000 |
| **均值** | **+0.1871** | **2.1484** | **1.5363** | **0.5935** | **+0.0956** | **+0.1334** | 0.9667 |

---

## Data Audit

`experiments/v5_0/manifests/data_manifest.csv`：登记 **30** 个关键文件，**30/30 存在**，
全部记录 `size / SHA256 / rows / cols / mtime`。要点：

| 文件 | 大小 (B) | rows × cols | 说明 |
|---|---:|---|---|
| `data/ppb_affinity_proteinea/filtered.csv` | 5,933,651 | 8207 × 21 | 主基准（官方划分来源） |
| `features/v2.2_proteinea_features.csv` | — | 8207 × 78 | 物理 + 元数据特征 |
| `features/v2_esm2_cache.npz` | 52,481,097 | 8207 | pooled ESM-2（640×2 链，fp16） |
| `features/v2_esm2_perres_v2.npz` | 1,913,516,901 | 6885 | per-residue 原始（1.9 GB） |
| `features/v2.5_perres_cache.npz` | 582,502,134 | 6885 | per-residue 紧凑（582 MB） |
| `data/v4.0/ddg_pairs.csv` | 2,238,640 | 15052 × 17 | ΔΔG 配对（含 reverse） |
| `data/v4.0/ddg_split.json` | 82,187 | 8 keys | V4 三重隔离 split |

### 数据规模复核（**与任务 prompt 一致，但修正了 2 处**）

| 项 | prompt 值 | 实测值 | 一致? |
|---|---|---|---|
| 原始总数据 | ~8207 | **8207** | ✅ |
| WT×mutant 正向 pair | ~7526 | **7526** | ✅ |
| 映射特征后 | ~7523 | **7523**（丢 3 对） | ✅ |
| complexes | ~231 | **231** | ✅ |
| ΔΔG>0 比例 | ~74.5% | **74.5%** | ✅ |
| single-point | ~68.6% | **68.6%** | ✅ |
| alanine scanning | ~59% | **59.0%** | ✅ |
| 最大单 complex pairs | ~573 | **573** | ✅ |
| 划分 | — | train 5265 / val 1129 / test 1129；复合物 118 / 51 / **61** | 补充 |
| **test ΔΔG<0 比例** | **未提及** | **37.3%**（而 train 74.5%） | ⚠️ **新发现** |
| **train/test 均值差** | **未提及** | train +1.798 vs test **+0.523**（Δ **−1.275**） | ⚠️ 新发现 |

> ⚠️ **V5 新发现的关键协变量偏移**：训练集 74.5% 的 ΔΔG>0（削弱结合），
> 而测试集只有 **37.3%** 为 ΔΔG>0。这直接影响所有「命中 ΔΔG<0」类指标的解读
> （详见 Top-K Utility 一节），并说明**测试集不是训练集分布的随机抽样**。

---

## Leakage Audit

### 自动断言（`scripts/v5/02_baseline_repro.py`，9/9 通过）

| 检查 | 值 | 期望 | 结论 |
|---|---:|---:|---|
| `pdb_overlap_train_val` | 0 | 0 | ✅ |
| `pdb_overlap_train_test` | 0 | 0 | ✅ |
| `pdb_overlap_val_test` | 0 | 0 | ✅ |
| `family_overlap_train_val` | 0 | 0 | ✅ |
| `family_overlap_train_test` | 0 | 0 | ✅ |
| `family_overlap_val_test` | 0 | 0 | ✅ |
| `split_exclusive`（三侧互斥） | 1 | 1 | ✅ |
| `split_covers_all`（全覆盖） | 1 | 1 | ✅ |
| `pca_fit_rows_from_train_only` | 1 | 1 | ✅ |

### 其他泄漏防护（逐条对应 V5 §87–§90）

- **PCA**：`fit` 只在 train 行的（WT ∪ MUT）并集上，变换 train/val/test。已断言 test 行不在 fit 集。
- **标准化**：`standardize_train` 的 median/std 只在 train 行拟合。
- **OOF/stacking**：辅助任务用 `GroupKFold(pdb_id)`，每折断言
  `len(intersect(train_fold, oof_fold)) == 0`。
- **forward/reverse**：正反向 pair 由同一 `pair_index` 映射，**必然同侧**（split 以 pair 为单位）。
- **重复测量**：同一 `(pdb_id, mutation)` 不在跨侧（属于同一 pair 的 WT/MUT 行）。

### 主动上报的自身缺陷（**4 项**，全部不掩盖）

**缺陷 1 — P4 train OOF 的折级系统偏移**

P4 的 conformal 校准所用 train OOF 残差存在**折级系统性偏移**
（fold 0：`y.mean=+2.590` 而 `pred.mean=−1.421`，偏差 **−4.011**；
反向预测 Pearson = **+0.4110**，证明不是符号 bug 而是**折级分布漂移**）。
→ 处置：如实标注为**方法学风险**（`logs/p4_oof_diagnosis.md`），
并说明 constant-q conformal 恰好吸收该偏移（这是其设计目的），
但代价是 **q 被折级偏差主导，无法反映样本级不确定度**。

**缺陷 2 — `efficiency_at_10` 曾系统性反号（`INVALID_DUE_TO_BUG`，已修正并全部重算）**

- **发现方式**：用**合成 oracle 自检**（把真值当预测输入）发现
  「完美 oracle 得 **−0.6667**、完全反向得 **+1.0000**」—— 符号完全颠倒。
  若不做这个自检，该列会被当作「模型比随机差」而写入报告，**结论方向会反向**。
- **根因**：每个 complex 的 random 基线**只抽 1 次**，在成员数少的 complex 上方差极大，
  使分母 `v_or − v_rnd` 可正可负，符号不稳定。
- **修正**：random 基线改为 **200 次抽样取均值**；efficiency 统一定义为
  `(random − model) / (random − oracle)`。合成自检复核：
  oracle **+1.0000**、anti **−0.9816**、noisy **+0.7644** ✅
- **处置（V5 §123 全流程）**：保存证据（`logs/BUGFIX_efficiency_at_10.md`）→
  写 reproduction（合成自检表）→ 修 bug → **重跑全部下游脚本**
  （`scripts/v5/16_rerun_after_bugfix.py`）→ 旧产物重命名为 `*.INVALID_DUE_TO_BUG`
  → **不只更新最终表**。
- **影响范围（已逐一核查）**：仅 `efficiency_at_10` / `utility_*` / `regret_at_10`。
  Pearson / Spearman / macro per-PDB / RMSE / MAE / AUROC / P@k / EF@x% / Hit@k
  **不依赖该函数，因此不受影响** → **P2/P3/P4 的主要结论不变**。

---

## P0 External Benchmark

### 已从一手来源核实（`LITERATURE_EVIDENCE`）

**来源**：CATH-ddG, *Bioinformatics* 2025, **41(Suppl 1): i362–i372**, doi:10.1093/bioinformatics/btaf228
（PMID 40662779 / PMC12261453）。全文 XML 与官方仓库数据已落盘至 `F:\dsh\downloads\cathddg\`。
**官方代码**：https://github.com/ak422/CATH-ddG —— **`license` 字段为 `None`（仓库未声明许可）**，已记录。

**官方 split 实测**（从仓库 `data/SKEMPI2/skempi_v2_with_all_results.csv` 的 `cath_fold` 列读出）：
`train` 5893 / `val` **813** / `NaN` 379（共 7085）。论文所称「held-out CATH blind test set
（813 mutations / 53 PPIs）」**即 `val` 分区**。

### 官方 held-out CATH 测试集上的 SOTA（Table 2）

| 方法 | Overall PearsonR | Overall SpearmanR | **per-PPI SpearmanR** | RMSE | MAE | AUROC |
|---|---:|---:|---:|---:|---:|---:|
| FoldX | 0.4603 | 0.5246 | 0.4303 | 2.3083 | 1.6837 | 0.7541 |
| flex ddG | 0.6081 | 0.6096 | 0.4542 | 2.0643 | 1.4946 | 0.7637 |
| ProteinMPNN | 0.1837 | 0.3668 | 0.2543 | 2.5559 | 1.8612 | 0.6788 |
| ESM-IF | 0.1421 | −0.1813 | −0.1331 | 2.5738 | 1.9376 | 0.3916 |
| RDE-Network | 0.4588 | 0.4703 | 0.2880 | 2.3103 | 1.7029 | 0.7452 |
| DiffAffinity | 0.3144 | 0.3626 | 0.2494 | 2.4684 | 1.7393 | 0.6251 |
| PPIFORMER | 0.2654 | 0.3145 | 0.2164 | 2.5069 | 1.8191 | 0.6740 |
| Prompt-DDG | 0.3538 | 0.3970 | 0.3025 | 2.4320 | 1.7556 | 0.6851 |
| DDAffinity | 0.3984 | 0.4026 | 0.1956 | 2.3849 | 1.7871 | 0.7248 |
| **CATH-ddG** | **0.6150** | **0.6269** | **0.4940** | **2.0504** | **1.5047** | **0.7805** |

> 我们的对应读数：Overall Pearson **+0.3861**、**macro per-PDB Spearman +0.0956**。

### 判定：`NOT_DIRECTLY_COMPARABLE`（V5 §19）

**四条事实性差异**（不做辩解性论证，只陈述）：

1. **数据集不同**：我们 = PPB-Affinity（8207 行，五源聚合）；他们 = SKEMPI v2.0（6706 条）。
2. **划分机制不同**：我们 = 复合物级 + ≤40% 序列同一性 + 家族级；他们 = CATH 同源超家族
   hold-out + 序列同一性 **<30%** + **剔除全部 easy mutation（最大 TM-score <0.6）**。
   → **他们的约束更强**，而训练量仍有 5893 条（我们 5265 配对 / 118 复合物）。
3. **任务粒度不同**：他们逐条突变（single/multiple 分开）；我们是 WT×mutant 配对差分 + 正反向增广。
4. **表示不同**：他们用 3D 几何编码器（ProteinMPNN 主干）+ FoldX 能量项 + 微调 ESM2 + CATH 域自监督
   （`num_layers=3`、`hidden_dim=128`、`k1=20/k2=7/k3=3`、`patch_size=256`、3 个模型集成）；
   我们只用 ESM-2 **池化**嵌入 128 维 + 60 个物理特征。

**唯一允许的陈述**：双方各自最严格协议下的「组内排序能力」分别为
**+0.0956（我们，per-PDB）** 与 **0.4940（CATH-ddG，per-PPI）**，差距约 **0.40**，
与 V5「池化表示丢失位点信息」的诊断一致。

### 我们的模型跑官方 split 的可行性核查（`PROJECT_EVIDENCE`）

| 检查项 | 结果 |
|---|---|
| 官方 test（813 条）的 complex 在 PPB-Affinity 中 | **53 / 53 = 100%** |
| 官方 test 逐条突变可被 PPB-Affinity 覆盖 | **813 / 813 = 100%** |
| 官方 test 单点突变 | 543 条，**543/543 可覆盖** |
| 官方 train（5893 条）的 complex 覆盖 | 228 / 287 = 79.4% |
| 官方 train 单点突变覆盖 | 3604 / 4404 = 81.8% |
| **(PDB, 突变) 精确交集** | **0** —— PPB-Affinity 不含 CATH 的突变条目 |

**阻塞原因（精确）**：pooled ESM-2 缓存是为 PPB-Affinity **自己的**序列提取的；
跑官方 split 必须为 CATH 的 WT/mutant 序列**重新做 ESM-2 推理**。
本机 **CPU-only**，且 **`esm` / `fair-esm` 未安装、无本地权重缓存**（实测：
`esm` MISSING、`fair_esm` MISSING、`C:\Users\KELE\.cache\huggingface` 仅 1 文件 0.0 MB）。

→ **标记 `NOT_RUN_RESOURCE_BLOCKED`**（V5 §0.1 / §16 / §36）。
已完成的替代工作：官方 split 标签获取 ✅、官方数值逐表提取 ✅、sign convention 核实 ✅、
可运行流程写入 `p0_report.md` §6 与 `p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md` ✅。

**已核实的 sign convention**（V5 §17）：CATH 表 `ddG` 均值为正、正值为削弱结合；
我方 `ΔΔG = ΔG_mut − ΔG_WT`（`ΔG = 0.592·ln KD`）同样以正值为削弱结合 →
**两侧一致，无需反号**。我方实测 mean **+1.5215**、ΔΔG>0 占 74.5%。

**SOTA claim 等级（V5 §95）**：**LEVEL 0 `PROMISING_INTERNAL_RESULT`**。本报告**不主张**任何对外部 SOTA 的胜负结论。

---

## P1 MIPA

### 状态：`NOT_RUN_RESOURCE_BLOCKED`（双重阻塞）

**阻塞 A：GATE B mutation mapping 覆盖率不足**（V5 §124 `<70%` 分支）。

**GATE B 全过程**（三个版本，含一次自我推翻）：

| 版本 | 脚本 | 自报覆盖 | 独立验证 | 处置 |
|---|---|---:|---|---|
| v1 | `04_build_mutation_sites.py` | 74.06%（5574 pair） | `05` 用真实 PDB 残基号验证：E1 通过率 **7.96%**（568/7134）；99/226 pdb offset 自相矛盾 | **`INVALID_DUE_TO_BUG`** |
| v2 | `06_build_residue_map.py`（序列引导） | 61.35%（4828/7869 位点） | 残基一致率 **76.82%**（5930/7719）；链 identity=1.000 占 63.06% | 仍不足 |
| v3（严格） | `11_strict_exact_subset.py` | **13.12%**（987/7523 pair；test **189**/1129） | 位置编号 + 氨基酸身份**双重精确匹配** | `SOLVE_MAPPING_FIRST` |

**v1 的根因（已定位）**：只解出「dataset 序列内下标」，**从未确定 dataset 序列 ↔ 具体 PDB 链**；
`05` 又用「取最长链」的错误假设去验证 → 得到 7.96% 的通过率。这是一次
**由独立验证捕获的自身错误**，v1 的 `numbering_offset` / `pdb_resseq` 列已标
`INVALID_DUE_TO_BUG`，**不用于任何科学结论**（V5 §0.4 / §123）。

**前置检查的关键事实**（`03a_seq_sanity.py` / `03b_case_study.py`）：

- 配体侧 **36.4%**、受体侧 **24.1%** 的 pair **两侧序列完全相同**；
- 变化模式：`lig+rec` 39.53% / 仅 `rec` 36.37% / 仅 `lig` 24.10%；
- 多链拼接（含 `,`）：配体 30.03%、受体 14.76%；
- 朴素 1-based 位置命中率 73.37%（0-based 仅 0.43%）→ 存在**编号 offset**；
- 失配样例呈**恒定偏移**（如 1A22 全部 +10）→ 编号差异而非 construct 不一致；
- 但 v2 的 offset 分布仍有 `|offset|` 达 1745 的**假匹配**，故必须收紧到 offset=0。

**结论（严格表述）**：**不得**把该结果解读为「MIPA 假设被证伪」。
只能说「**在现有数据可映射性下，MIPA 无法被检验**」。
这本身是一条有信息量的负结论：**数据集只提供 `Ligand Sequences` / `Receptor Sequences`
（多链 `,` 拼接），没有「dataset 序列 ↔ 具体 PDB 链」的显式对应表，也没有独立 WT 序列字段** ——
这是**数据构建层面的限制**，不是方法层面的失败。

**阻塞 B：ESM-2 不可用**。MIPA 需要 `h_m^WT` / `h_m^mut` / `h_j^partner` 的**逐残基**表示；
本机无 `fair-esm`、无本地权重、无 GPU。

### 已交付的可运行流程

`experiments/v5_0/p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md`：含提取脚本骨架
（`esm2_t30_150M_UR50D`、第 30 层、只缓存 mutation site + local partner、
fp16 embedding + float32 距离 + offset 数组、每 640D block 在 train-only 上 PCA 到 32 维）、
必实现断言、P1-0…P1-6 消融矩阵、coverage fairness 三档
（`FULL baseline` / `SUBSET-matched baseline` / `MIPA subset`）、成功判据。

### 对 V5 §104 五个问题的回答

| 问题 | 回答 |
|---|---|
| **Q1** mutation-site Δh 是否比 pooled ΔZ 更有预测力？ | **无法回答**（`NOT_TESTED_RESOURCE_BLOCKED`）。阻塞于 mapping 覆盖 + ESM-2 缺失。 |
| **Q2** partner environment 是否提供独立增益？ | **无法回答**（同上）。 |
| **Q3** elementwise perturbation × partner interaction 是否关键？ | **无法回答**（同上）。 |
| **Q4** 提升主要发生在 cross-target 还是 within-target？ | **间接证据**：pooled 表示在 **within-target 上几乎无能力**（macro per-PDB ρ = +0.0956，CI 含 0；组内 Top-10 比随机更差）。这从反面支持「位点信息才是 within-target 的关键」。 |
| **Q5** single vs multi mutation 行为是否不同？ | **未做**（原计划在 MIPA 消融内）。 |

---

## P1B Inverse Folding

**状态：`NOT_RUN_RESOURCE_BLOCKED`。**

- 需 ProteinMPNN 权重做 `L_m = log P(a_mut|structure) − log P(a_WT|structure)`；
  本机 **无 GPU**，且历史上 V2.3/V2.5 的 ESM-2/per-residue 提取均在外部 GPU 平台完成。
- 另需**同一套 mutation↔残基映射**（与 P1 共用），而该映射正是 GATE B 的失败点。
- 未做任何「只搬 USP-ddG 的 FoldX / MoE / CATH encoder」的尝试（V5 §45 明确禁止）。

**按 V5 §46，成功判据待将来执行时使用**：`baseline + log-odds` 的
`ΔPearson ≥ +0.02` 或 `Δmacro per-structure Spearman ≥ +0.04` 且方向一致。

---

## P2 State Potential

**状态：DONE。判定：全部失败（GATE D = NO）。**

设计与执行（`scripts/v5/07_p2_state_potential.py`）：

- **P2-0** `pair_plain`：V4 的 `C_D_fwd+rev`（无 ΔG 辅助）
- **P2-1** `pair_dgaux`：V4 最强 `C+DGaux_fwd+rev`（主对照）
- **P2-2** `state_absolute`：训练 `h(x)→ΔG`，`ΔΔG = h(mut) − h(WT)`
- **P2-3** `state_centered`：`h(x)→ΔG_centered`（减 train 内 per-complex 均值），再差分

**State 行统计**：唯一变体行 **5108**（train 3241 / test 999，**交集 0**）；
state 模型 fit 行 = **3241**（剔除任何在 test 出现的变体）→ **无 state 泄漏**。

| seed | P2-0 pair_plain | **P2-1 pair_dgaux** | P2-2 state_absolute | P2-3 state_centered |
|---:|---:|---:|---:|---:|
| 42 | +0.3674 | **+0.4040** | +0.0797 | −0.0758 |
| 43 | +0.3434 | **+0.3740** | +0.0343 | −0.1104 |
| 44 | +0.3497 | **+0.3676** | +0.0642 | −0.0596 |
| **均值** | **+0.3535** | **+0.3819** | **+0.0594** | **−0.0819** |

**macro per-PDB Spearman**：P2-0 +0.0933 / **P2-1 +0.1199** / P2-2 −0.0055 / P2-3 −0.0298。

### GATE D 判定（3-seed 平均，vs P2-1）

| variant | ΔPearson | Δmacro-PDB ρ | 判定 |
|---|---:|---:|---|
| P2-0 pair_plain | −0.0284 | −0.0266 | 未达标准 |
| **P2-2 state_absolute** | **−0.3225** | **−0.1255** | **未达标准** |
| **P2-3 state_centered** | **−0.4639** | **−0.1498** | **未达标准** |

配对 cluster bootstrap（单位 = PDB，n_clusters = 61，same-seed）：
P2-2 vs P2-1 `ΔPearson = −0.3225`，**BH q = 0.0480**（显著更差）；
P2-3 vs P2-1 `ΔPearson = −0.4639`，**BH q = 0.0400**（显著更差）。

### 一致性指标（V5 §54）

- **state 模型的 `E_anti` = 0 由构造保证**（`ΔΔG = h(mut) − h(WT)` 天然反对称）。
- pair 模型的 `E_anti` 用「独立反向训练模型」估计（下界）。
- **这是 P2 唯一的正面结论**：state 结构确实买到了严格的数学一致性，
  但**没有买到任何预测性能**。

### 对 V5 §105 四个问题的回答

| 问题 | 回答 |
|---|---|
| **Q1** exact antisymmetry 是否带来性能提升？ | **否。** 反对称性由构造保证（E_anti = 0），但 ΔPearson **−0.3225**（BH q = 0.0480）。 |
| **Q2** transitivity 是否改善？ | **无法受益**：本数据只有 WT→mutant 边，无 mutant→mutant 实验边，无法构造真实三元组。 |
| **Q3** absolute ΔG supervision 是否真的帮助 ΔΔG？ | **间接否**：P2-1 证明「把 ΔG_mut 的 OOF 预测作为**一列特征**」有稳定正收益（+0.0284 vs P2-0），但**把 ΔG 当作要差分的 state 函数**则严重更差。→ **监督信号有用，state 差分形式有害。** |
| **Q4** complex condition offset 是否通过 state difference 被部分抵消？ | **未观察到**：P2-3（显式去 per-complex 均值）反而最差（−0.4639）。 |

**失败价值（V5 §56）**：结论是「**objective / thermodynamic consistency 已不是主要瓶颈**」，
研究资源应集中到 **representation**。这与 Executive Summary 的最终判定一致。

---

## P3 Family-balanced Training

**状态：DONE。判定：全部失败（成功标准 §66 未满足）。**

### Train 内 complex 规模（仅 train 统计，**未使用任何 test 信息**）

| 项 | 值 |
|---|---:|
| train complex 数 | **118** |
| train pairs | 5265 |
| 每 complex pairs：中位 / 均值 / 最大 | 12 / 44.6 / **573** |
| **最大 complex 占 train rows** | **10.88%** |
| top-5 complex 的 pairs | [573, 534, 358, 300, 282] |
| **top-5 合计占比** | **38.9%** |

### P3-A sample weighting（test n=1129, n_complex=61）

| seed | W0 equal | W1 1/√n | W2 1/n |
|---:|---:|---:|---:|
| 42 | **+0.3951** | +0.3342 | +0.1911 |
| 43 | **+0.3843** | +0.3270 | +0.2227 |
| 44 | **+0.3790** | +0.3659 | +0.2381 |
| **均值** | **+0.3861** | **+0.3424** | **+0.2173** |

| variant | ΔPearson vs W0 | Δmacro-PDB ρ | 判定 |
|---|---:|---:|---|
| W1 1/√n | **−0.0437** | −0.0332 | **DROP** |
| W2 1/n | **−0.1688** | −0.0673 | **DROP** |

→ **反向复杂度加权一律降低性能**，且幅度远超任何「更公平」的收益。

### P3-B complex-balanced bagging（`--skip-bagging` 后只保留 K=5 试点）

| seed | K | 聚合 | Pearson | macro-PDB ρ |
|---:|---:|---|---:|---:|
| 42 | 5 | pred_mean | +0.1525 | −0.0293 |
| 42 | 5 | rank_mean | +0.1212 | −0.0268 |

→ **比 W0 基线（+0.3951）低 0.24 以上**。按 V5 §125「不允许 endless experimentation」，
在 K=5 试点已明确为负结果后**未继续全量运行**（如实记录该决定）。

### P3-D cross-split 稳定性（V5 §65/§66，12 个随机 complex-level 划分，128 树）

| scheme | n_splits | mean | **SD** | range |
|---|---:|---:|---:|---|
| W0_equal | 12 | +0.1930 | **0.2025** | [−0.1593, +0.5059] |
| W2_inv_n | 12 | +0.1616 | **0.1753** | [−0.1457, +0.5005] |

- SD 变化：0.2025 → 0.1753（**−13.4%**）
- **GATE：要求 SD 下降 ≥20% → 未满足。**

> **最重要的观察**：**cross-split SD ≈ 0.20**，而所有加权方案带来的 mean 变化
> 都在 ±0.17 以内 → **划分随机性比任何加权策略都更重要**。
> 这与 V5 §12 的 cluster bootstrap 结论一致（真实独立样本量 = 118 个训练复合物）。

### 对 V5 §106 三个问题的回答

| 问题 | 回答 |
|---|---|
| **Q1** 大 complex 是否支配 baseline？ | **是**：top-5 complex 占 train 的 **38.9%**，最大单个占 **10.88%**。 |
| **Q2** inverse-complex weighting 是否提高 macro metrics？ | **否，反而显著降低**（W1 −0.0437、W2 −0.1688；macro-PDB ρ 同样下降）。 |
| **Q3** 是否降低 split-to-split variance？ | **略降但不达阈值**：SD **−13.4%**（0.2025 → 0.1753），要求 ≥20%。 |

---

## P4 Conformal Top-K

**状态：DONE。判定：Top-K 零改善（GATE E = NO）。**

### 方法（严格按 V5 §69–§73，**未使用 LambdaRank / xendcg**）

- train 内 complex-level GroupKFold（5 折）→ OOF 残差 `e_i = y_i − ŷ_i`
- one-sided conformal：`U(x) = ŷ(x) + q_{1−α}`，按 `U` 从小到大选；α ∈ {0.10, 0.20, 0.25}
- Mondrian：按 mutation type（丙氨酸/非丙氨酸）× mutation count（单/多）分组校准，
  组样本 < 100 时自动回退全局
- **test 只用于最终评估**；q 只来自 train OOF

### 结果：plain 与所有 conformal 变体**完全相同**

| variant | P@10 | EF@1% | Hit10(<−0.5) | sel10 mean ΔΔG |
|---|---:|---:|---:|---:|
| `plain` | 0.9667 | 9.6034 | 0.9667 | −6.7258 |
| `conformal_global_a10/20/25` | 0.9667 | 9.6034 | 0.9667 | −6.7258 |
| `mondrian_a10` | 0.9667 | 9.2700 | 0.9667 | −5.3423 |
| `mondrian_a20` | 0.9667 | 9.3034 | 0.9667 | −5.0607 |
| `mondrian_a25` | 0.9333 | 9.0032 | 0.9333 | −4.6760 |

**数学必然性**：global conformal 的 `q` 是**常数**（对所有候选相同），
因此 `U = ŷ + q` 与 `ŷ` **同序** → Top-K 与 plain **逐位相同**。
这不是实现失败，而是方法本身的必然结果。Mondrian 的 q 依赖样本（按变异类型/个数），
因此**改变了排序，但只让指标变差**（EF@1% 9.6034 → 9.0032~9.3034，
`sel10 mean ΔΔG` −6.73 → −4.68~−5.34）。

配对检验：24 个对比中 **23 个 Δ = 0.0000**，其余为负 → **全部 `DROP`**。

### Calibration coverage 诊断（V5 §79，train OOF）

| nominal | empirical coverage | q |
|---:|---:|---:|
| 90% | **0.9000** | 3.2169 |
| 80% | — | — |
| 75% | **0.7500** | 3.2169 |

→ **校准本身是正确的**（empirical ≈ nominal）。按 V5 §78：
「若 calibration 误差 coverage 正确，但 Top-K 不提高 → **问题不是 confidence estimation，
而是主模型没有足够局部 ranking signal**」。

### 对 V5 §107 四个问题的回答

| 问题 | 回答 |
|---|---|
| **Q1** prediction uncertainty 是否经过 calibration？ | **是**（empirical coverage ≈ nominal）。 |
| **Q2** 高置信候选是否真的更准？ | **未检验出增益**：constant-q 不改变排序；Mondrian 的组内 q 反而更差。 |
| **Q3** Top-10 命中是否改善？ | **否**（P@10/EF@1%/Hit10 全部不优于 plain，Δ ≤ 0）。 |
| **Q4** OOD proteins 上 uncertainty 是否更大？ | **未观察到正向关系**（见 OOD Analysis，误差与相似度非单调）。 |

---

## Statistical Analysis

### 主推断口径：cluster bootstrap（V5 §12）

**这是 V5 相对于 V4 最重要的方法学升级**：V4 使用 row-level bootstrap，
把 **1129 行当成 1129 个独立样本**；而同一复合物的多个突变体高度相关。

| 指标（seed 42） | point | 95% CI（**cluster, n_clusters = 61**） |
|---|---:|---|
| overall Pearson | +0.3951 | [+0.0730, +0.6743] |
| **macro per-PDB Spearman** | **+0.1048** | **[−0.0219, +0.2385]** |
| macro per-PDB Pearson | +0.1307 | [−0.0063, +0.2712] |

| 指标（seed 43） | point | 95% CI |
|---|---:|---|
| overall Pearson | +0.3843 | [+0.0357, +0.6690] |
| macro per-PDB Spearman | +0.0905 | [−0.0379, +0.2180] |

| 指标（seed 44） | point | 95% CI |
|---|---:|---|
| overall Pearson | +0.3790 | [+0.0346, +0.6649] |
| macro per-PDB Spearman | +0.0914 | [−0.0335, +0.2171] |

**判读**：overall Pearson 的 CI 很宽但**不含 0**（跨复合物区分确实存在）；
**macro per-PDB Spearman 的 CI 在所有三个种子上都包含 0**。

> **这是 V5 最硬的一条结论**：在正确的抽样单位下，
> 「同一复合物内排序」这一能力**统计上不可与随机区分**。

### 多重比较与 practical significance（V5 §13）

所有配对比较使用 **BH-FDR**（`q < 0.05`）+ practical threshold：
overall Pearson `Δ ≥ 0.02`；macro per-PDB Spearman `Δ ≥ 0.04`；
P@10 `Δ ≥ 0.10`；EF@1% 相对 `≥ +20%`；Hit@10 `Δ ≥ 0.10`。
不满足 practical 但统计显著的自动标 `STATISTICALLY_SIGNIFICANT_BUT_SMALL`。

**`paired_comparisons.csv` 共 24 行，决策分布**：

| decision | 行数 |
|---|---:|
| `DROP` | **24** |
| `KEEP` | **0** |

其中显著（BH q < 0.05）的只有 P2-2（q = 0.0480）与 P2-3（q = 0.0400）——
**两者都是显著更差**。

### 同 seed 配对（V5 §84）

所有比较**严格同 seed**（candidate s42 vs baseline s42）。
V4 曾出现跨 seed 配对把 BH q 从 3/6 显著变成 q = 0.0512 的教训，V5 未重犯。

---

## OOD Analysis

### Error vs nearest-train similarity（V5 §80，Figure 2）

以 pooled WT 嵌入对最近 train 变体嵌入的**最大余弦相似度**为 OOD 度量，
按四分位分箱（seed 42 基线）：

| 分位 | n | mean sim | **mean \|error\| (kcal/mol)** |
|---|---:|---:|---:|
| Q1（最 OOD） | 287 | 0.522 | **1.233** |
| Q2 | 294 | 0.682 | **2.019** |
| Q3 | 292 | 0.847 | **1.353** |
| Q4（最相似） | 256 | 0.950 | **1.481** |

**判读（如实报告，不做美化）**：
- 误差与相似度**不是单调关系**：最相似的四分位（sim 0.950）误差 **1.481**，
  反而**高于**最 OOD 的四分位（sim 0.522）的 **1.233**。
- 最高误差出现在 Q2（sim 0.682，误差 2.019）。
- → **不能用「越 OOD 越差」的简单叙事解释本模型**；本模型的误差结构
  主要由**复合物间的标签尺度差异**（test ΔΔG mean 仅 +0.523 而 train +1.798）驱动，
  而非由序列/结构相似度驱动。

> ⚠️ 该结论与「相似度度量本身」有关：这里用 pooled 嵌入余弦作为代理。
> V5 §63 已明确警告不要用 3-mer Jaccard 当精准 biological identity；
> 本处使用嵌入余弦是**可用的最好代理**，但仍是代理。

### 条件一致性分组（V5 §91，仅作 error analysis）

已知条件混杂（`PROJECT_EVIDENCE`，来自 V4 记录并已复核）：
- 温度字段覆盖仅 **14.5%**；方法字段 **30.8%**
- V4 实测：**仅 54.6% 的 ΔΔG 配对温度一致、57.2% 方法一致**；
  条件不一致的配对 `|ΔΔG|` 中位 **1.5023** vs 一致的 **1.0121**（放大 **48%**）
- **未做强对抗去条件化**（V5 §91 明确禁止），也**未根据 assay subgroup 删除任何 test 行**（V5 §92）

### 协变量偏移（V5 新量化）

| 项 | train | test | Δ |
|---|---:|---:|---:|
| ΔΔG mean | **+1.798** | **+0.523** | **−1.275** |
| ΔΔG<0 占比 | **74.5%** | **37.3%** | **−37.2 pp** |
| 用 `Z_diff` 区分 train/test 的 AUC | — | ≈ **0.699** | — |

→ 测试集与训练集**分布显著不同**（AUC 0.699 意味着有相当强的可分性）。
这是**跨靶泛化的真实难度**，也解释了为什么「命中 ΔΔG<0」类指标在 test 上数值较低。

---

## Top-K Utility

**这是 V5 最重要的诊断，也是对 V4 数字最重要的修正。**

### 全局池 vs 复合物内（V5 §76 要求同时报告；以 `plain` 部署模型为准）

| 层 | Hit@10(ΔΔG<0) | Hit@10(<−0.5) | Hit@10(<−1.0) | 选中 10 个的 mean ΔΔG |
|---|---:|---:|---:|---:|
| **Global（全部 1129 候选）** | **0.9667** | 0.9667 | 0.9667 | **−6.7258** |
| **Within-complex（真实场景）** | **0.3080** | 0.1720 | 0.1173 | **+0.3577** |
| Random 基线（within） | — | — | — | **+0.6277** |
| **Oracle 上界（within）** | — | — | — | **−0.5709** |
| 全局基准率（随机选 10 的 ΔΔG<0 期望） | **0.3729** | — | — | — |
| 全体候选的 mean ΔΔG | — | — | — | **+0.5232** |

**within 层 `P@10` 三方对照**：模型 **0.1653** / 随机 **0.1401** / oracle 0.3920。

### 四项分解读数

1. **全局层**：选中 10 个的平均 ΔΔG = **−6.7258**，比全体均值（+0.5232）低
   **7.2490** kcal/mol → 这部分优势**几乎全部来自跨复合物区分**
   （挑中的复合物整体更负）。
2. **复合物内层**：模型比随机好 **+0.2700** kcal/mol；oracle 比随机好
   **+1.1986** kcal/mol → **只回收 22.5%**。
3. **`Efficiency@10 = (random − model) / (random − oracle) = +0.2252`**
   （0 = 与随机无异；1 = 达到 oracle）。
   > 交叉核对：`v5_common.topk_metrics` 走「≥5 成员的 complex」子集（n=59）
   > 给出 **+0.1014**。两条路径口径不同（成员数下限不同），结论一致 ——
   > **模型只回收了约 10–23% 的 oracle 可用空间。**
4. **关键矛盾**：within 层 `hit10(ΔΔG<0) = 0.3080`，**低于全局基准率 0.3729**
   → **在一个靶点内挑 10 个「会增强结合」的突变体，命中率不如随机。**

> **READ CAREFULLY（避免误读）**：
> - `Efficiency@10` 为**正**（+0.2252），表示按模型排序选出的 10 个**确实优于随机**
>   （+0.3577 vs 随机 +0.6277，ΔΔG 越负越好）—— 但**只回收了 22.5% 的可用空间**。
> - 然而 `hit10(ΔΔG<0)`（0.3080）**低于**该 test 集上随机选择的期望（0.3729）。
>   两个读数不矛盾：模型倾向于选中「相对更负但仍为正值」的候选，
>   因此在「是否 <0」这个二值判据下反而不如随机。
>   这暴露了**测试集基准率大幅偏移**（test 仅 37.3% 为 ΔΔG<0，而 train 74.5%）的影响。
> - 无论取哪个读数，结论一致：**模型不具备可用的组内挑选能力**。

### Figure 4 说明

`fig4_topk_enrichment.png` 左图显示全局池曲线（model 与 oracle 几乎重合——
这正是「被跨复合物效应主导」的视觉证据）；右图显示复合物内 regret，
模型曲线**远离 0**（oracle 水平），说明组内选择能力缺失。

### 对 V5 §108/§129 的直接回答

> **V5 的核心产品目标（「实验只能做 10 个 clones 时选哪 10 个」）
> 在当前模型下不被满足。** 全局 P@10 = 0.9667 会在报告中产生
> 「模型已经能挑出好突变」的错误印象，而真正可用的指标
> （复合物内 Hit@10 = 0.3034 < 基准率 0.3729）说明**相反的事实**。

---

## Failure Analysis

### 失败清单（全部保留，不删除）

| # | 方向 | 结果 | 证据 | 标签 |
|---|---|---|---|---|
| 1 | MIPA mutation→residue mapping v1 | 覆盖 74.06% 但独立验证仅 7.96% | `manifests/mapping_pdb_validation.md` | **`INVALID_DUE_TO_BUG`** |
| 2 | MIPA v2（序列引导） | 残基一致率 76.82%，覆盖 61.35% | `manifests/gate_b_residue_map_report.md` | `SOLVE_MAPPING_FIRST` |
| 3 | MIPA v3（严格 offset=0） | 覆盖 13.12%，test 仅 189 pair | `manifests/strict_exact_report.md` | 统计功效不足 |
| 4 | P2-2 state absolute | ΔPearson **−0.3225**（BH q = 0.0480） | `p2_state_potential/p2_results.csv` | `DROP` |
| 5 | P2-3 state centered | ΔPearson **−0.4639**（BH q = 0.0400） | 同上 | `DROP` |
| 6 | P3 W1 1/√n 加权 | ΔPearson **−0.0437** | `p3_balancing/p3_results.csv` | `DROP` |
| 7 | P3 W2 1/n 加权 | ΔPearson **−0.1688** | 同上 | `DROP` |
| 8 | P3 complex bagging K=5 | ΔPearson **−0.24** 以上 | 同上 | `DROP`（K=5 后停止） |
| 9 | P3 cross-split SD 目标 | SD 仅降 **13.4%**（要求 ≥20%） | `p3_balancing/p3_stability.csv` | 未达 GATE |
| 10 | P4 global conformal | Top-K **逐位相同**（常数 q） | `p4_conformal/p4_results.csv` | 数学必然 |
| 11 | P4 Mondrian conformal | EF@1% 9.6034 → 9.0032~9.3034 | 同上 | `DROP` |
| 12 | P0 外部 benchmark 评估 | 未运行 | — | **`NOT_RUN_RESOURCE_BLOCKED`** |
| 13 | P1B inverse folding | 未运行 | — | **`NOT_RUN_RESOURCE_BLOCKED`** |
| 14 | `efficiency_at_10` 反号（自身 bug） | 合成 oracle 自检得 −0.6667 | `logs/BUGFIX_efficiency_at_10.md` | **`INVALID_DUE_TO_BUG`**（已修 + 全部重算） |

### 主动上报的自身缺陷（3 项）

1. **GATE B v1 的映射是错误的**（7.96% 通过率）。由**独立的 PDB 验证**捕获，
   而非靠alignment 自证。已标 `INVALID_DUE_TO_BUG` 并重做两版。
2. **P4 的 train OOF 存在折级系统偏移**（fold 0 偏差 −4.011；
   反向预测 Pearson = +0.4110 证明非符号 bug）。如实标为方法学风险。
3. **P3 bagging 未全量运行**（K=5 试点为负后按 V5 §125 停止）。
   这是**有意识的成本决策**，已在报告中显式记录，不隐藏。

### 失败的模式（跨版本规律）

V5 的四条 CPU 主线（P2/P3/P4 + MIPA）失败原因**高度一致**：

> **它们都在「训练策略/目标函数/校准」层面操作，而瓶颈在「输入表示」。**

- P2 改变了学习目标的形式（pair → state difference）→ 更差
- P3 改变了样本权重 → 更差
- P4 改变了选择准则（用同一批预测）→ 无法改善（因为预测本身没有组内分辨力）
- MIPA（唯一针对表示的方向）→ 被数据可映射性阻塞

---

## Final Model

### 组成（V5 §93 规则：只有独立验证有效的模块才能进入）

| 模块 | 独立验证结果 | 是否入选 |
|---|---|---|
| 特征：物理 60 + pooled ESM-2 PCA128 | 继承 V2.2–V2.3.2，V5 Gate A 逐位复现 | ✅ 入选 |
| 特征块 `D = x_mut − x_wt`（反对称） | V4 创新 1；V5 复现 `C_D_fwd+rev` | ✅ 入选 |
| 正/反向增广 | V4 创新 1；V5 Gate A 复现 | ✅ 入选 |
| ΔG_mut OOF 辅助列 | V4 创新 2；V5 Gate A **逐位复现**（Δ ≤ 1e-6） | ✅ 入选 |
| LightGBM 固定 800 树（不用早停） | V4 教训；V5 沿用 | ✅ 入选 |
| P2 state potential | ΔPearson −0.32 | ❌ 不入选 |
| P3 complex/family 加权、bagging | 全部下降 | ❌ 不入选 |
| P4 conformal 选择器 | Top-K 零改善 | ❌ 不入选 |
| MIPA | 未验证 | ❌ 不入选 |

### 最终模型 = V4 基线，无任何 V5 模块入选

```
ΔΔG_hat = LightGBM_800(
             X = [ phys_diff(60) , Z_diff(128) , ΔG_mut_OOF(1) ]     # 189 维
          )
  训练集 = train(5265) ∪ reverse(train)      # 反向 y 取负，辅助列不取负
  PCA     仅 train 拟合；标准化仅 train 拟合
  随机种子 42 / 43 / 44（同 seed 配对比较）
```

**3-seed 均值**：Pearson **+0.3861**、Spearman +0.1871、macro per-PDB ρ **+0.0956**、
RMSE 2.1484、MAE 1.5363、AUROC 0.5935。

> **诚实陈述**：V5 **没有**产出比 V4 更好的模型。
> V5 的产出是**对这个模型的适用边界给出了严格界定**。

---

## Comparison with Literature

### 唯一允许的横向陈述

| 来源 | 协议 | 指标 | 数值 | 可否比较 |
|---|---|---|---:|---|
| **CATH-ddG**（*Bioinformatics* 2025, btaf228, Table 2） | 官方 held-out CATH blind test（53 PPI / 813 mut；seq-id <30%；无 easy mutation） | **per-PPI SpearmanR** | **0.4940** | 与我方 per-structure 口径**最接近**，但数据集/任务粒度不同 → `NOT_DIRECTLY_COMPARABLE` |
| CATH-ddG | 同上 | Overall PearsonR | **0.6150** | 同上 |
| **本项目 V5** | PPB-Affinity，复合物级 + ≤40% seq-id + 家族级三重隔离（test 61 complex / 1129 pair） | **macro per-PDB Spearman** | **+0.0956** | — |
| 本项目 V5 | 同上 | Overall Pearson | **+0.3861** | — |

### 其他已核实的 SOTA 数字（`LITERATURE_EVIDENCE`）

- **PPIFORMER split**（Table 3，14 PPI / 219 mut）：CATH-ddG Pearson **0.63** / Spearman **0.56**；
  flex ddG 0.57/0.55；FoldX 0.56/0.53；PPIFORMER 0.46/0.44；RDE-Network 0.30/0.24
- **HER2 case study**（Table 4，419 变体）：CATH-ddG **0.607/0.639**；DiffAffinity 0.590/0.524；
  GearBind+Ensemble 0.579/0.608
- **SARS-CoV-2 RBD**（Table 5，285 变体）：CATH-ddG Pearson **0.579**；DDAffinity 0.456；FoldX 0.385

### 三条结构性差异（解释差距的关键，且我方**不占优**）

1. **表示**：CATH-ddG 用**逐残基几何编码器**（ProteinMPNN 主干，`num_layers=3`）+
   显式建模「突变位点 + 界面」的 **mutational microenvironment**（256 近邻残基 patch）；
   我方用**池化**嵌入（整链平均）→ **位点信息被平均掉**。
2. **能量项**：CATH-ddG 集成 FoldX 能量项（`loss_foldx` 权重 1.0）；
   我方 V2.5 已证伪「能量总和拼入 GBDT」（BH q = 0.0080 更差）—— 但**这不能推出
   「FoldX 无用」**，因为**接入方式不同**（他们按残基/位点对齐，我方是全复合物求和）。
3. **自监督**：CATH-ddG 有 CATH 域多标签自监督（`loss_cath` 权重 0.5）+ 微调 ESM2 + 3 模型集成；
   我方为纯监督 + 单模型 ×（多种子复现）。

> ⚠️ **严禁的反向推论**：不得写「我们的 0.386 低于他们的 0.615，所以我们方法更差」——
> 两者 split、数据集、任务粒度全不同（V5 §19）。同理**不得**写「我们更严格所以其实更好」。

---

## Limitations

1. **算力是硬约束**：CPU-only（4 核 i7-2630QM / 7.55 GB RAM，实测空闲常仅 1.5–1.9 GB）。
   导致：MIPA 无法提取逐残基表示；P0 无法在官方 CATH split 上评估；
   P3 bagging 只跑了 K=5 试点；per-residue CNN 类方法全部不可行。
2. **GATE B 未通过**：mutation→residue 映射覆盖率 13.12%（严格）/ 61.35%（宽松），
   均 < 70%。**这使得 MIPA 这一 V5 最高优先级假设完全未被检验。**
3. **test 规模极小**：**61 个复合物 / 1129 配对**。cluster bootstrap 的 CI 因而很宽
   （overall Pearson CI 宽达 0.60）。任何「显著」结论都需谨慎。
4. **协变量偏移未解决**：test ΔΔG mean +0.523 vs train +1.798（`Z_diff` AUC ≈ 0.699）。
   **未**采用任何 test adaptation（V5 §60 禁止）。
5. **条件混杂未解决**：温度覆盖 14.5%、方法覆盖 30.8%；仅 ~55% 配对条件一致。
6. **state potential 的 transitivity 未被真正检验**：数据只有 WT→mutant 边，
   无 mutant→mutant 实验边，无法构造真实三元组闭环。
7. **conformal 的 train OOF 有折级偏差**（fold 0 偏差 −4.011），
   使 q 被折级偏移主导，无法反映样本级不确定度。
8. **OOD 度量是代理**：用 pooled 嵌入余弦，非真实序列/结构同一性。
9. **`NOT_DIRECTLY_COMPARABLE`**：与所有文献数字均不可直接比较（数据集/划分/指标口径不同）。
10. **非 git 仓库**：无 commit hash，只能以文件 SHA256 + mtime 作为版本证据。
11. **P3-B bagging 未全量运行**（有意识的成本决策，已记录）。
12. **P4 的 α 只测 3 个值、P2 的 λ 未做 train-only 选择实验**（V5 §51 的 latent graph potential
    部分未实现，因为 P2-A 已明显失败，按 V5 §125 不再深入该线）。

---

## Reproducibility Commands

所有命令在本机（`F:\DSH_Galatea\Affinity_Intelligence`）实测可运行。

```powershell
# ---------- PHASE 0: 环境与数据清单 ----------
python -u scripts/v5/00_env_snapshot.py            # -> audit/environment.md（追加，不覆盖）
python -u scripts/v5/01_data_manifest.py           # -> manifests/data_manifest.{csv,md}

# ---------- GATE A: V4 基线复现（核心，必须先通过） ----------
python -u scripts/v5/02_baseline_repro.py --seeds 42 43 44 --trees 800
#   -> p0_benchmark/baseline_repro_{results.csv,predictions.npz,report.md}
#   期望: gateA PASS, Pearson = .3951/.3843/.3790（逐位）

# ---------- GATE B: mutation 映射（三版，含一次自我推翻） ----------
python -u scripts/v5/03a_seq_sanity.py             # 前置：WT/mut 序列一致性
python -u scripts/v5/03b_case_study.py             # 前置：annotation 可恢复性
python -u scripts/v5/04_build_mutation_sites.py    # v1（结果 INVALID_DUE_TO_BUG）
python -u scripts/v5/05_validate_mapping_pdb.py    # v1 的独立验证 -> 7.96% -> 推翻
python -u scripts/v5/06_build_residue_map.py       # v2 序列引导 -> 76.82% 残基一致
python -u scripts/v5/11_strict_exact_subset.py     # v3 严格 EXACT -> 13.12%

# ---------- P0: 外部 benchmark（数值已取，评估被算力阻塞） ----------
python -u scripts/v5/10_p0_write_benchmark.py
#   -> p0_benchmark/p0_sota_benchmark.csv (37 行), p0_source_provenance.md, p0_report.md

# ---------- P2: Thermodynamic State Potential ----------
python -u scripts/v5/07_p2_state_potential.py --seeds 42 43 44
#   -> p2_state_potential/p2_{results.csv,predictions.npz,report.md}  (~24 min)

# ---------- P3: Complex/Family Balanced Boosting ----------
python -u scripts/v5/08_p3_balancing.py --skip-bagging
#   -> p3_balancing/p3_{results.csv,stability.csv,report.md}  (~5 min)
#   注: 去掉 --skip-bagging 会跑 20 members × K∈{5,10,20}（极慢，且 K=5 已为负）

# ---------- P4: Conformal Top-K ----------
python -u scripts/v5/09_p4_conformal.py            # -> p4_conformal/*  (~4 min)
python -u scripts/v5/13_diag_p4_oof.py             # OOF 折级偏差诊断
python -u scripts/v5/14_p4_decompose_topk.py       # ** 关键: 全局 vs 复合物内 Top-K **

# ---------- 汇总与图 ----------
python -u scripts/v5/12_consolidate_results.py     # -> all_experiments.csv, paired_comparisons.csv
python -u scripts/v5/15_make_figures.py            # -> figures/fig1..fig6 (~1 min)

# ---------- 被阻塞实验的可运行流程（需 GPU + ESM-2 权重） ----------
# 见 experiments/v5_0/p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md
# 见 experiments/v5_0/p0_benchmark/p0_report.md §6
```

**外部数据溯源**（`PROJECT_EVIDENCE`，已落盘至 `F:\dsh\downloads\cathddg\`）：

| 文件 | 大小 (B) | 获取方式 |
|---|---:|---|
| `cathddg_oa.bin`（Europe PMC 全文 XML, PMC12261453） | 246,881 | Europe PMC REST |
| `skempi_v2_with_all_results.csv`（含官方 `cath_fold`） | 1,291,260 | GitHub blobs API |
| `CATH.yml` | 2,373 | GitHub contents API (base64) |
| `zero_domain_foldseek.csv` | 1,960 | GitHub contents API (base64) |
| `README.md` | 3,160 | raw.githubusercontent |

> 按 V5 §120：URL / 下载日期 / 大小已记录；**官方仓库 `license` 为 `None`**（未声明许可）已如实记录。
> 按 V5 §121：未执行任何来源不明的 shell script；未写入任何 API key。
> 注：`raw.githubusercontent.com` 对本仓库部分文件返回 0 字节，改用 GitHub API 获取（已记录）。

---

## Recommended V5.1

按「投入产出比 × 可验证性」排序。**每条都明确指出它绕开了 V5 的哪个失败原因**
（V5 §3 要求）。

### R1（最高优先）—— 修复 mutation↔residue 映射，这是解锁 MIPA 的唯一前置

- **做法**（V5.1 应严格照此，不要重复 V5 的错误假设）：
  1. **不要**从 mutant 序列反推 WT；
  2. 用 **PDB 官方 SEQRES**（而非仅 ATOM）与 dataset 序列做 global alignment；
  3. 建立 **SEQRES 下标 ↔ ATOM 残基号**双表，显式处理缺失残基、
     insertion code、altLoc、链 ID 重映射；
  4. 对 `,` 拼接的**多链分别**做子链级 alignment，**不允许跨子链共享 offset**；
  5. **必须用独立证据验证**（V5 的教训：alignment 自证会给出 7.96% 的假通过率）。
- **绕开的失败原因**：V5 GATE B 的三次失败全部源于「拿 alignment 结果自证」+
  「用最长链假设 PDB 链对应」。R1 用**官方注释 + 逐残基双表 + 独立验证**取代之。
- **预期**：把严格 EXACT 覆盖从 13% 提到 ≥90%，从而**真正检验 MIPA**。
  这是 V5 §130 最高优先级假设（`Δh_m ⊙ c_partner`）的唯一入口。
- **失败退路**：若仍 <70%，则**公开结论**为「此数据的构建体标注不足以支持位点级建模」
  —— 这本身是可发表的负面结果，并把资源转向 R2。

### R2 —— 用「已知 WT 结构 + 官方突变注释」重建位点对齐表示（不依赖数据集序列）

- **做法**：完全绕开 dataset 的 `Ligand/Receptor Sequences`，改为
  本地 3031 个 PDB 的 **SEQRES 链序列**作为 WT 的唯一真源，
  再用官方突变注释（`CHAIN_ORIGppMUT`）在该序列上定位。
- **绕开的失败原因**：V5 的映射失败根因是**依赖 dataset 序列与 PDB 的隐式对应**；
  R2 让 PDB 成为真源，dataset 只提供标签。
- **预期**：即使 R1 部分成功，R2 也能独立给出更干净的位点定义。
- **失败退路**：若 PDB 与注释的链对应仍不可靠（多链、同源二聚体），
  则改用「**仅单链 PPI + 明确链 ID**」子集，牺牲覆盖换正确性。

### R3 —— 把「划分严格性」做成显式 4 级基准（不提升任何指标，但产出最可发表资产）

- **做法**：在**同一模型**上系统跑 4 档划分（mutation-level random / complex-level /
  ≤40% seq-id / ≤30% seq-id + 家族级），每档报同一套指标 + cluster bootstrap。
- **绕开的失败原因**：V5 已实测「random vs 家族级差 +0.50」「cross-split SD ≈ 0.20」，
  但从未把这一现象**系统化为一篇方法论结果**。R3 把它变成主贡献。
- **预期**：量化「文献高数字中有多少来自划分宽松」。V5 已有 3 个独立实例
  （V2.0 +0.1799、V2.2 +0.5034、V4/V5 +0.25）。
- **失败退路**：若各档差异不如预期显著，则至少产出「划分不是唯一因素」的诚实结论。

### R4 —— 只在 within-complex 口径上评估，并把「全局 P@10」从主指标中移除

- **做法**：把 V5 新增的 `within-complex Hit@10 / Efficiency@10 / Regret` 作为
  **主评估指标**；全局 P@10 / EF@1% 降级为诊断量，并**强制同时报基准率**。
- **绕开的失败原因**：V5 发现全局 P@10 = 0.9667 会掩盖「复合物内 Hit@10 = 0.3034
  低于基准率 0.3729」这一事实。不修指标口径，后续任何改进都会被假象误导。
- **预期**：不再出现「指标很好但实际不可用」的情况。
- **失败退路**：无（这是评估纪律，不是方法创新）。

### R5（低优先，需 GPU）—— 仅当 R1 成功后，才尝试 MIPA-IF 组合

- 按 V5 §112：`MIPA`（learned local interaction perturbation）+
  ProteinMPNN（structure-conditioned residue compatibility）。
- **前置条件（硬性）**：MIPA 与 inverse folding **各自独立验证有效**后才允许组合（V5 §114）。
  V5 两者都未验证 → **V5.1 不得提前实现该组合**。

### 明确**不建议**的方向（V5 已证伪，V5.1 不要重复）

| 方向 | V5 证据 |
|---|---|
| 加权/bagging 处理 complex 不平衡 | W1 −0.0437、W2 −0.1688、bagging −0.24 |
| state-potential 差分形式 | ΔPearson −0.3225（q = 0.0480） |
| conformal 改善 Top-K | 零改变（常数 q）/ 更差（Mondrian） |
| 堆特征 / 能量总和 / 生物学代理 / 热力学分解 | 继承 V2.5：21 个特征集 **0** 显著 |
| 排序损失（LambdaRank/xendcg） | 继承 V2.0/V2.4：两次失败 |
| 单纯换更大 ESM | 继承 V3.1 记录 + AINN-P1 证据 |
| 用 NN 替换 GBDT（表格特征） | 继承 V2.3.2/V3.1/V4.0：三次证伪 |

---

## 附：产物索引

| 路径 | 内容 |
|---|---|
| `experiments/v5_0/README.md` | 执行台账（DONE/FAILED/BLOCKED/NEXT） |
| `experiments/v5_0/audit/environment.md` | 环境快照（追加式，不覆盖） |
| `experiments/v5_0/manifests/data_manifest.{csv,md}` | 30 个文件 SHA256 + rows/cols |
| `experiments/v5_0/manifests/{seq_sanity,case_study,gate_b_mapping_report,gate_b_residue_map_report,mapping_pdb_validation,strict_exact_report}.md` | GATE B 全过程（含负结果） |
| `experiments/v5_0/p0_benchmark/p0_sota_benchmark.csv` | 37 行外部 SOTA（V5 §18 列） |
| `experiments/v5_0/p0_benchmark/{p0_report,p0_source_provenance}.md` | P0 结论 + 来源溯源 |
| `experiments/v5_0/p0_benchmark/baseline_repro_*` | GATE A 产物 |
| `experiments/v5_0/p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md` | MIPA 可运行流程（被阻塞） |
| `experiments/v5_0/p2_state_potential/*` | P2 结果 + 报告 |
| `experiments/v5_0/p3_balancing/*` | P3 结果 + 稳定性 + 报告 |
| `experiments/v5_0/p4_conformal/*` | P4 结果 + 校准 + **Top-K 分解** |
| `experiments/v5_0/all_experiments.csv` | 统一实验登记（V5 §99 列） |
| `experiments/v5_0/paired_comparisons.csv` | 配对比较 + BH-FDR + 决策（V5 §100 列） |
| `experiments/v5_0/final/{results_summary.md,variant_summary.csv}` | 汇总 |
| `experiments/v5_0/figures/fig1..fig6*.png` | V5 §98 要求的 6 张图 |
| `scripts/v5/*.py` | 15 个可复现脚本 |
| `experiments/v5_0/logs/` | 全部运行日志 + `p4_oof_diagnosis.md` + `failures.log` |
| `data/v5_0/*` | mutation 映射产物（含 `INVALID_DUE_TO_BUG` 的 v1 与修正后的 v2/v3） |

### V5 对既有文件的保护

- **未修改任何 V1.1–V4.0 文件**（全部只读）。
- V5 产物全部落在 `scripts/v5/`、`experiments/v5_0/`、`data/v5_0/`。
- 外部下载全部落在 `F:\dsh\downloads\`（按磁盘纪律，未写入 C 盘）。
