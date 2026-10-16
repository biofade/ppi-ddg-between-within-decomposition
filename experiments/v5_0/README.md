# V5.0 — Affinity Intelligence 执行台账

> 本文件记录每个阶段的 **DONE / FAILED / BLOCKED / NEXT**。
> 规则：负结果与失败同等保留；已确认有 bug 的结果标 `INVALID_DUE_TO_BUG`；
> 结论区分 `PROJECT_EVIDENCE` / `LITERATURE_EVIDENCE` / `HYPOTHESIS`。

## 关键规范（V5）

- **历史 PPB test 自 V5 起定义为 `DEVELOPMENT_BENCHMARK`**，不再是 blind test
  （V2→V4 已多次查看，存在 adaptive test reuse）。任何 SOTA 声明必须另依赖
  外部 CUT（CATH hold-out / 官方 SKEMI split / 新 lockbox / 未触碰外部盲集）。
- **主统计推断 = cluster bootstrap**（抽样单位 family，不足则 PDB；携带该 cluster
  全部 mutation rows，5000 次）。row-level bootstrap 仅作 legacy 次级分析。
- **禁止** test 参与 PCA / 标准化 / 特征选择 / 超参 / 早停 / conformal calibration /
  stacking。V5 所有选择只用 train（conformal 额外只用 train OOF）。
- 决策标签由代码规则生成（`v5_common.decide`），不人工挑数字。

## 阶段状态

| 阶段 | 状态 | 产物 | 备注 |
|---|---|---|---|
| PHASE 0 目录骨架 | **DONE** | `experiments/v5_0/{audit,manifests,configs,logs,p0_benchmark,p1_mipa,p1b_inverse_folding,p2_state_potential,p3_balancing,p4_conformal,final,figures}` | 未触碰任何 V1–V4 文件 |
| PHASE 0 环境快照 | **DONE** | `audit/environment.md` | Windows 10 / Python 3.11.9 / i7-2630QM 4 核 / 7.55 GB RAM / **无 GPU** / torch 2.5.1+cpu / LightGBM 4.7.0 / numpy 2.4.6 / pandas 3.0.6 / scipy 1.17.1 / sklearn 1.9.1 / Biopython 1.88；**不是 git 仓库 → `GIT_NOT_AVAILABLE`** |
| PHASE 0 数据清单 | **DONE** | `manifests/data_manifest.{csv,md}` | 登记 30 个关键文件，**30/30 存在**，全部记录 size + SHA256 + rows/cols + mtime |
| **GATE A 基线复现** | **DONE ✅ PASS** | `p0_benchmark/baseline_repro_{results.csv,predictions.npz,report.md}` | 3 seeds **逐位复现**（Δ = −0.0000/−0.0000/+0.0000）；见下 |
| GATE A 泄漏审计 | **DONE ✅ 9/9 通过** | 同上 | PDB 交集 0/0/0；family 交集 0/0/0；split 互斥且全覆盖；PCA fit 行 100% 来自 train |
| GATE B 序列一致性前置检查 | **DONE** | `manifests/seq_sanity.md`, `case_study.md` | 见下 |
| GATE B mutation 映射 v1 | **INVALID_DUE_TO_BUG** | `04_build_mutation_sites.py` 产出 | 见下方 bug 记录 |
| GATE B mutation 映射 v1 验证 | **DONE（负结果）** | `manifests/mapping_pdb_validation.md` | E1 通过率仅 **7.96%** → 推翻 v1 |
| GATE B mutation 映射 v2 | **DONE** | `manifests/gate_b_residue_map_report.md`, `data/v5_0/residue_map.csv` | 残基一致率 **76.82%**；覆盖 **61.35%** → `SOLVE_MAPPING_FIRST` |
| P0 外部 SOTA benchmark | **NOT_RUN_RESOURCE_BLOCKED** | 待定 | 需 CATH-ddG / USP-ddG 官方 split；见下方 BLOCKED 说明 |
| P1 MIPA | **IN_PROGRESS（限 EXACT 子集）** | `p1_mipa/` | 按 V5 §34 必须做 matched-baseline |
| P1B inverse folding | NOT_RUN_RESOURCE_BLOCKED | `p1b_inverse_folding/` | 需 ProteinMPNN 权重；无 GPU |
| P2 state potential | **DONE（负结果）** | `p2_state_potential/` | `GATE D = NO`：state absolute ΔPearson **−0.3225**（q=0.0480）、centered **−0.4639**（q=0.0400） |
| P3 family-balanced training | **DONE（负结果）** | `p3_balancing/` | W1 −0.0437、W2 −0.1688、bagging −0.24；cross-split SD 仅降 13.4%（要求 ≥20%） |
| P4 conformal Top-K | **DONE（零改善）** | `p4_conformal/` | constant-q 不改变排序（数学必然）；Mondrian 更差；**Top-K 分解为核心诊断** |
| Bug 修正重算链 | **DONE** | `logs/BUGFIX_efficiency_at_10.md` | `efficiency_at_10` 曾反号 → 标 `INVALID_DUE_TO_BUG` → 修 → 全部重跑 |

---

## V5.0 三个决定性结论（PROJECT_EVIDENCE）

### 1. 真实独立样本量是 **61**，不是 1129（cluster bootstrap）

| 指标（seed 42） | point | 95% CI（cluster, n_clusters = 61） |
|---|---:|---|
| overall Pearson | +0.3951 | [+0.0730, +0.6743] |
| **macro per-PDB Spearman** | **+0.1048** | **[−0.0219, +0.2385]**（**含 0**） |

三个种子的 macro per-PDB Spearman CI **全部含 0** → **组内排序能力统计上不可与随机区分**。

### 2. 全局 Top-K 的高分是**跨复合物效应 + 高基准率**的假象

| 层 | Hit@10(ΔΔG<0) | 选中 10 个的 mean ΔΔG |
|---|---:|---:|
| Global（1129 候选） | 0.9667 | −6.7258 |
| **Within-complex（真实场景）** | **0.3080** | **+0.3577** |
| Random（within） | — | +0.6277 |
| Oracle（within） | — | −0.5709 |
| 全局基准率 | **0.3729** | — |

→ `Efficiency@10 = +0.2252`（只回收 22.5% 的 oracle 空间）；
within 层 `hit10(ΔΔG<0) = 0.3080` **低于基准率 0.3729**。
**test 集仅 37.3% 为 ΔΔG<0，而 train 为 74.5%**（协变量偏移）。

### 3. 四条 CPU 主线全部为负结果

| 主线 | 结论 |
|---|---|
| P2 state potential | 反对称性由构造保证（E_anti=0），但 ΔPearson **−0.3225** |
| P3 加权 / bagging | 全部**下降**（−0.0437 ~ −0.24） |
| P4 conformal | Top-K **零改善**（常数 q 不改排序） |
| P1 MIPA | **GATE B 阻塞**（严格 EXACT 覆盖 13.12%） |

**最终判定：`REPRESENTATION_BOTTLENECK_CONFIRMED`**（V5 §128）

完整报告：`V5_0_FINAL_REPORT.md`

---

## GATE B 全过程（PROJECT_EVIDENCE）

### B-1 序列一致性前置检查（`03a_seq_sanity.py`）

| 侧 | \|Δlen\| 中位 | p95 | 序列相同比例 | identity 中位 |
|---|---:|---:|---:|---:|
| 配体 | 0 | 223 | 36.4% | 0.9975 |
| 受体 | 0 | 223 | 24.1% | 0.9924 |

- 变化模式：`lig+rec` 39.53% / 仅 `rec` 36.37% / 仅 `lig` 24.10%
- 多链拼接（含 `,`）：配体 30.03%，受体 14.76%
- 朴素 1-based 位置命中率 **73.37%**（0-based 仅 0.43%）→ 存在**编号 offset**
- 失配样例呈**恒定偏移**（如 1A22 全部 +10）→ 是编号差异而非 construct 不一致

### B-2 v1 映射（`04`）与其**独立验证失败**（`05`）

- v1 报告：主实验可用 5574 pair（**74.06%** 覆盖），offset 范围 [−425, +645]
- **独立验证（用真实 PDB 残基序号）**：E1 通过率仅 **568 / 7134 = 7.96%**；
  99 / 226 个 pdb 的 offset **自相矛盾**；11.0% 位点 |offset| > 100
- **根因**：v1 只解出「dataset 序列内下标」，**从未确定 dataset 序列 ↔ 具体 PDB 链**；
  `05` 又用「取最长链」的错误假设去验证。
- **处置**：v1 的 `numbering_offset` / `pdb_resseq` 列标 **`INVALID_DUE_TO_BUG`**，
  不用于任何科学结论（V5 §0.4 / §123）。

### B-3 v2 映射（`06`，按 V5 §25「PDB sequence ↔ dataset sequence global alignment」）

| validated_status | 位点数 | 占比 |
|---|---:|---:|
| `EXACT` | 3798 | 48.27% |
| `AA_MISMATCH` | 1737 | 22.07% |
| `LOW_IDENTITY_CHAIN` | 1102 | 14.00% |
| `ALIGNMENT_CONFIDENT` | 1030 | 13.09% |
| `INDEX_NOT_IN_ALIGNMENT` | 150 | 1.91% |
| `GAP_AT_MUTATION_SITE` | 52 | 0.66% |

- **链 identity = 1.000 的位点占 63.06%**；≥0.95 占 **81.07%** → 链对应基本找对
- **独立验证：`pdb_resseq_real` 处真实残基 == annotation WT 残基 = 5930 / 7719 = 76.82%**
  （v1 为 7.96%）
- 可用位点（EXACT + ALIGNMENT_CONFIDENT）= **4828 / 7869 = 61.35%**，覆盖 **3587** pair

### GATE B 判定

按 V5 §124 阈值（≥90% full / 70–90% subset / <70% 先解决）：

> **GATE B = `SOLVE_MAPPING_FIRST`（61.35%）**

**处置（不伪造、不强行解释）**：
1. **不做** full MIPA。
2. P1 退化为 **EXACT-only 最小消融**，且必须在**同一 intersection rows** 上
   重跑 baseline（V5 §34 coverage fairness：`FULL baseline` / `SUBSET-matched baseline` /
   `MIPA subset` 三者并列）。
3. 同时输出 `P1_GPU_RUN_INSTRUCTIONS.md`（V5 §36）：因为要真正解决映射，
   需要用 **PDB 官方 SEQRES + 完整体 WT 序列**重做，而不是从 mutant 序列反推。
4. 不把 61.35% 的映射结果当作「MIPA 被证伪」——**只能**说「在当前映射质量下
   无法得出关于 MIPA 的结论」。

---

## GATE A 结果（PROJECT_EVIDENCE）

### 复现判定

| seed | V4 记录 Pearson | V5 复现 Pearson | Δ | 判定 |
|---:|---:|---:|---:|---|
| 42 | +0.3951 | +0.3951 | −0.0000 | ✅ |
| 43 | +0.3843 | +0.3843 | −0.0000 | ✅ |
| 44 | +0.3790 | +0.3790 | +0.0000 | ✅ |

**结论：V4 最强配置在本环境逐位复现，无 baseline mismatch（无需 legacy/current 双轨）。**

### V5 补齐的指标（V4 未报告，PROJECT_EVIDENCE）

| seed | Spearman | RMSE | MAE | AUROC | **macro-PDB Spearman** | macro-PDB Pearson | Hit@10 (ΔΔG<0) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 42 | +0.1972 | 2.1357 | 1.5249 | 0.5981 | **+0.1048** | +0.1307 | 0.9000 |
| 43 | +0.1997 | 2.1529 | 1.5364 | 0.6046 | **+0.0905** | +0.1257 | 1.0000 |
| 44 | +0.1640 | 2.1566 | 1.5476 | 0.5778 | **+0.0914** | +0.1438 | 1.0000 |

**三 seed 均值：overall Pearson = +0.3861，macro-PDB Spearman = +0.0956。**

### Cluster bootstrap（V5 主推断，抽样单位 = PDB，n_clusters = **61**）

| 指标 | point | 95% CI |
|---|---:|---|
| overall Pearson (seed 42) | +0.3951 | [+0.0730, +0.6743] |
| **macro-PDB Spearman (seed 42)** | **+0.1048** | **[−0.0219, +0.2385]** |
| macro-PDB Pearson (seed 42) | +0.1307 | [−0.0063, +0.2712] |

> **关键读数**：真实的独立样本量是 **61 个复合物**，不是 1129 行。
> 在这个口径下，**within-complex 排序能力（macro-PDB Spearman ≈ +0.10）
> 的 95% CI 包含 0** —— 即「同一复合物内挑突变体」这一核心能力
> **在统计上尚不能与随机区分**。这正面验证了 V5 的出发点。

### 数据规模（复核，与 prompt 一致）

- 配对：**7523**（forward；15052 行含 reverse；映射到 v2.2 后丢弃 3 对）
- 划分：train 5265 / val 1129 / test 1129；复合物 train 118 / val 51 / **test 61**
- 特征：`D`（反对称块）**188** 维；`S`（对称块）**188** 维
- 辅助任务 `ΔG_mut` 折外 Pearson 均值 **+0.4368**（5 折），test **+0.4684**

---

## BLOCKED 记录

### P0 外部 SOTA benchmark — `NOT_RUN_RESOURCE_BLOCKED`

**阻塞原因**：需要 CATH-ddG（Bioinformatics 2025, DOI 10.1093/bioinformatics/btaf228）、
USP-ddG（Bioinformatics 2026, btag249）、PreFold-dG（btag489）的**官方 split 文件**
与官方报告数值。按 V5 §96「不要引用二手博客」，必须从论文原文/官方 repo 读取。

**已完成的替代工作**：见 `p0_benchmark/` 下的文献核对脚本与说明。
**未做**：不训练其大网络（按 V5 §16 只做「同 split 跑我们的模型」）。

---

## NEXT

1. **GATE B**：构建 mutation↔PDB 残基映射，覆盖率需 ≥90% 才做 full MIPA。
2. P1 MIPA 最小消融矩阵（P1-0…P1-6），同一 test rows 比较 + coverage fairness。
3. P2 / P3 / P4 为纯 CPU，可并行推进（但按 V5 §0.3 必须独立验证）。
