# V5.0 P0 — External SOTA Benchmark Gate

生成: 2026-10-13 20:17:23

## 1. 结论（一句话）

**外部 benchmark 的数字已从一手来源（论文全文 + 官方仓库）完整提取；但我们自己的模型在该官方 split 上的运行结果为 `NOT_RUN_RESOURCE_BLOCKED`（阻塞于 CPU-only 下的 ESM-2 推理成本，非网络问题）。**

## 2. 官方 held-out CATH 测试集上的 SOTA（LITERATURE_EVIDENCE）

> 来源: CATH-ddG Table 2；813 mutations / 53 PPIs；CATH 同源超家族 hold-out，序列同一性 <30%，无 easy mutation。

| 类别 | 方法 | Overall PearsonR | Overall SpearmanR | **per-PPI SpearmanR** | RMSE | MAE | AUROC |
|---|---|---:|---:|---:|---:|---:|---:|
| 能量 | FoldX | 0.4603 | 0.5246 | 0.4303 | 2.3083 | 1.6837 | 0.7541 |
| 能量 | flex ddG | 0.6081 | 0.6096 | 0.4542 | 2.0643 | 1.4946 | 0.7637 |
| 预训练 | ProteinMPNN | 0.1837 | 0.3668 | 0.2543 | 2.5559 | 1.8612 | 0.6788 |
| 预训练 | ESM-IF | 0.1421 | -0.1813 | -0.1331 | 2.5738 | 1.9376 | 0.3916 |
| 序列+结构预训练 | RDE-Network | 0.4588 | 0.4703 | 0.2880 | 2.3103 | 1.7029 | 0.7452 |
| 序列+结构预训练 | DiffAffinity | 0.3144 | 0.3626 | 0.2494 | 2.4684 | 1.7393 | 0.6251 |
| 序列+结构预训练 | PPIFORMER | 0.2654 | 0.3145 | 0.2164 | 2.5069 | 1.8191 | 0.6740 |
| 序列+结构预训练 | Prompt-DDG | 0.3538 | 0.3970 | 0.3025 | 2.4320 | 1.7556 | 0.6851 |
| 监督学习 | DDAffinity | 0.3984 | 0.4026 | 0.1956 | 2.3849 | 1.7871 | 0.7248 |
| 监督学习（本文） | **CATH-ddG** | **0.6150** | **0.6269** | **0.4940** | 2.0504 | 1.5047 | 0.7805 |

### 与我们的数字对照

| 项 | 我们（PROJECT_EVIDENCE） | CATH-ddG（LITERATURE_EVIDENCE） |
|---|---:|---:|
| Overall Pearson | 0.3861 | 0.6150 |
| Overall Spearman | 0.1871 | 0.6269 |
| **per-structure/PPI Spearman** | **0.0956** | **0.4940** |
| test pairs / complexes | 1129 / 61 | 813 / 53 |

**判定：`NOT_DIRECTLY_COMPARABLE`（V5 §19）** —— 理由（全部为事实，非辩解）：

1. **数据集不同**：我们用 PPB-Affinity（8,207 行，聚合 PDBbind/SAbDab/SKEMPI/Affinity-Benchmark/ATLAS 五源）；CATH-ddG 用 SKEMPI v2.0（6,706 条）。
2. **划分机制不同**：我们 = 复合物级 + ≤40% 序列同一性 + 家族级三重隔离；CATH-ddG = CATH 同源超家族 hold-out + 序列同一性 <30% + **剔除全部 easy mutation**（最大 TM-score <0.6）。**他们额外施加了「无 easy mutation」这一更强约束**，同时又拥有 5.9k 训练突变（我们 5.3k 配对但只涉 118 个复合物）。
3. **任务粒度不同**：他们逐条突变（single/multiple 分开报）；我们构造的是 WT×mutant **配对差分**并做正/反向增广。
4. **特征不同**：他们用 3D 几何编码器（ProteinMPNN 主干）+ FoldX 能量项 + 微调 ESM2 + CATH 域自监督；我们只用 ESM-2 **池化**嵌入 + 60 个物理特征。

> 因此**不得**写「我们 0.386 > 他们 0.615 之外的任何比较结论」。唯一允许的陈述是：**在双方各自最严格的协议下，per-PPI/组内排序能力分别为 0.096 与 0.494，差距约 0.40**，这与 V5 §4 诊断出的「pooled 表示丢失位点信息」一致。

## 3. 更严格协议下的其他 SOTA（Table 3，PPIFORMER split）

| 方法 | PearsonR | SpearmanR |
|---|---:|---:|
| FoldX | 0.56 | 0.53 |
| flex ddG | 0.57 | 0.55 |
| GEMME | 0.41 | 0.38 |
| MSA Transformer | 0.36 | 0.31 |
| ESM-IF | 0.18 | 0.18 |
| RDE-Network | 0.30 | 0.24 |
| PPIFORMER | 0.46 | 0.44 |
| DiffAffinity | 0.41 | 0.31 |
| Prompt-DDG | 0.13 | 0.12 |
| DDAffinity | 0.38 | 0.34 |
| CATH-ddG | 0.63 | 0.56 |

- 该 split 仅 14 PPIs / 219 mutations，是最挑剔的公开协议之一。

## 4. Case study 数字（Table 4/5）

| 数据集 | 方法 | PearsonR | SpearmanR |
|---|---|---:|---:|
| HER2 (419) | FoldX | 0.246 | 0.329 |
| HER2 (419) | GearBind+P | 0.515 | 0.517 |
| HER2 (419) | GearBind+Ensemble | 0.579 | 0.608 |
| HER2 (419) | RDE-Network | 0.465 | 0.481 |
| HER2 (419) | Prompt-DDG | 0.413 | 0.432 |
| HER2 (419) | DiffAffinity | 0.590 | 0.524 |
| HER2 (419) | DDGPred | 0.387 | 0.388 |
| HER2 (419) | GearBind | 0.478 | 0.475 |
| HER2 (419) | DDAffinity | 0.386 | 0.356 |
| HER2 (419) | CATH-ddG | 0.607 | 0.639 |
| SARS-CoV-2 RBD (285) | FoldX | 0.385 | — |
| SARS-CoV-2 RBD (285) | RDE-Network | 0.403 | — |
| SARS-CoV-2 RBD (285) | DiffAffinity | 0.305 | — |
| SARS-CoV-2 RBD (285) | DDAffinity | 0.456 | — |
| SARS-CoV-2 RBD (285) | CATH-ddG | 0.579 | — |

## 5. 我们的模型能否跑在完全相同的官方 split 上？（可行性核查）

**核查结论：结构层面完全可行，但受阻于 CPU 推理成本。**（以下全部为 PROJECT_EVIDENCE，来自 `data/` 与官方表的实际比对）

| 检查项 | 结果 |
|---|---|
| 官方 test（813 条）的 complex 是否都在 PPB-Affinity 中 | **53 / 53 = 100%** |
| 官方 test（813 条）逐条突变是否可被 PPB-Affinity 覆盖 | **813 / 813 = 100%** |
| 官方 test 单点突变 | 543 条，**543/543 可覆盖** |
| 官方 train（5,893 条）的 complex 覆盖 | 228 / 287 = 79.4% |
| 官方 train 单点突变覆盖 | 3,604 / 4,404 = 81.8% |
| **(PDB, 突变) 精确交集** | **0** —— PPB-Affinity 不含 CATH 的突变条目 |

**阻塞原因（精确表述）**：我们的 pooled ESM-2 缓存是为 PPB-Affinity 自己的 8,207 行序列提取的。要跑官方 CATH split，必须为 CATH 的 WT/mutant 序列**重新做 ESM-2 推理**（需先从 PDB 重建 WT 链序列，再逐条施加突变）。本机为 **CPU-only（4 核 i7-2630QM）**，ESM-2 150M 的 CPU 推理速度不足以在合理时间内完成数千条序列的提取（历史实测: V2.3 的 8,207 行提取是在 GPU 上完成的）。

**已完成的替代工作（按 V5 §0.1 的降级顺序）**：

1. ✅ 官方 split 标签已获取并落盘（`cath_fold` 列）；
2. ✅ 官方报告数值已从论文全文逐表提取（本文件 §2–§4）；
3. ✅ sign convention 已核实一致（见 `p0_source_provenance.md`）；
4. ✅ 可运行流程已写入 `p0_report.md` §6 与 `p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md`；
5. ❌ **未运行**我们模型在官方 split 上的评估 → 标记 `NOT_RUN_RESOURCE_BLOCKED`。

## 6. 可运行流程（一旦有 GPU 即可一次跑通）

```bash
# 1) 从 PDB 重建 WT 链序列（本地 3031 个 PDB）
python scripts/v5/11_p0_build_cath_pairs.py --pdb-dir data/pdb \
       --cath-csv F:/dsh/downloads/cathddg/skempi_v2_with_all_results.csv
# 2) 为 CATH WT/mutant 序列做 ESM-2 池化推理（需 GPU）
python scripts/v5/gpu_cath_esm2.py --device cuda --out features/v5_cath_esm2.npz
# 3) 在 CATH train 上训练，在官方 val(=blind test) 上评估
python scripts/v5/12_p0_eval_cath.py --cache features/v5_cath_esm2.npz
```

## 7. SOTA claim 等级（V5 §95）

> **LEVEL 0 — `PROMISING_INTERNAL_RESULT`**

理由：仅在内部 development benchmark 上有效；外部 benchmark 评估因算力受限未执行，**不具备** Level 2/3 的资格。本报告**不主张**任何对外部 SOTA 的胜负结论。

