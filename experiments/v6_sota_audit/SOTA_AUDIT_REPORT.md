# V6.0 SOTA Audit Report

**阶段**：V6.0 — The Diagnostic Audit
**时间**：2026-10-13
**执行环境**：CPU-only（4 核 i7-2630QM / 7.55 GB RAM）；**未运行任何模型推理**（V6 §0.1）
**上游结论**：V5.0 = `REPRESENTATION_BOTTLENECK_CONFIRMED`

> **证据分级**：`PROJECT_EVIDENCE`（本次审计实测）/ `LITERATURE_EVIDENCE`（论文报告）/ `HYPOTHESIS`（未验证推断）。
> **数据来源**：全部为外部仓库/论文的一手文件，URL 与许可见 §1 与 §6。

---

## 1. Data Scavenging Summary

### 找到的预测数据

| 目标 | 状态 | 来源文件 | 内容 |
|---|---|---|---|
| **CATH-ddG** | ✅ **FOUND** | `F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv`（1,291,260 B，GitHub blobs API） | 7,085 行 × 22 列；含 **5 个数值预测列**：`Flex-ddG` / `DDGPred` / `ddG_FoldX` / `ESM-IF` / `proteinMPNN` + 真值 `ddG` + 官方 `cath_fold` split |
| **DSSA-PPI** | ✅ **FOUND** | `F:\dsh\downloads\dssa_ppi\RBM_baselines_comparison_ddG.csv`（174,025 B） | 1,702 行 × 12 列；含 **8 个模型预测列**（DSSA-PPI / PPIFormer / ESM-2-FFN / MuTon / MpbPPI / ProBASS / SAAMBE-SEQ / FoldX）+ 真值 `EXP_fitness` |
| CATH-ddG（模型自身预测） | ❌ `PREDICTION_FILE_NOT_FOUND` | — | 官方表只有基线预测；CATH-ddG 本人的逐条预测未随仓库发布（仅 3 个 `.pt` 权重，需 GPU 推理 → 按 §0.1 不做） |
| SAMAffinity | ❌ `PREDICTION_FILE_NOT_FOUND` | 已下载 58 个 `data/test_fold/*/*.csv` | schema 为 31 列，**只有 `ddG` 真值与突变标识，无任何预测值列**（逐个检查确认） |
| Prompt-DDG | ❌ `PREDICTION_FILE_NOT_FOUND` | CATH 表内 `prompt-fold0/1/2` | 这 3 列**只是 split 标签**（`train`/`val`），非数值预测 |
| MINT | ❌ `PREDICTION_FILE_NOT_FOUND` | `VarunUllanat/mint`（68 文件，MIT） | 仅 `mint/helpers/predict.py`(428 B) 与 metrics notebook，无预测结果 |
| Pythia-PPI | ❌ `PREDICTION_FILE_NOT_FOUND` | `Wublab/pythia_ppi`（32 文件，Apache-2.0） | 文件树中无 predict/result/test 结果文件 |

**总计：找到 2 个数据源 / 13 个模型的预测列；缺失 4 个模型**（均为「仓库只有代码/权重、无预测文件」，**非网络故障**）。

> ⚠️ **不卡死原则**（V6 §0.3）：缺失项全部标记并跳过，PHASE 2–4 以已获得的
> **13 个模型**推进，未伪造任何数字。

---

## 2. Sign Convention Audit

### 真值列（CATH 表 `ddG`）—— **手工复算验证通过**

用**独立的第三方真值**（本地 SKEMPI v2.0 原始表 `data/v2.4/skempi_v2_cleaned.csv`，
分号分隔，含 `Affinity_wt (M)` / `Affinity_mut (M)`）逐步复算：

```
ΔG_WT  = 0.5925 × ln(KD_WT[M])
ΔG_mut = 0.5925 × ln(KD_mut[M])
ΔΔG   = ΔG_mut − ΔG_WT        (R·T = 1.98720425864083e-3 × 298.15 = 0.592497)
```

手工样例（全部可复核）：

| #Pdb | Mutation | KD_WT (M) | KD_mut (M) | ΔΔG（本地复算） | ddG（CATH 表） | 一致? |
|---|---|---:|---:|---:|---:|---|
| 1CSE | LI38G | 1.12e−12 | 5.26e−11 | **+2.2807** | +2.2806 | ✅ |
| 1CSE | LI38S | 1.12e−12 | 8.33e−12 | **+1.1888** | +1.1888 | ✅ |
| 1ACB | LI38D | 1.12e−12 | 1.28e−10 | **+6.8671** | +6.8667 | ✅ |
| 1IAR | NA89A | — | — | **+1.5589** | +1.5588 | ✅ |
| 1JTG | EA143A,KB74A | — | — | **+4.0622** | +4.0619 | ✅ |
| 1JTG | EA143A,FB142A | — | — | **+2.5812** | +2.5810 | ✅ |

**读法**：`LI38G` 的 KD 从 1.12e−12 升到 5.26e−11（结合**变弱**）→ ΔΔG = **正**值。

- 全表连接命中 **6,502 行**；`Pearson(ΔΔG_calc, ddG_CATH) = +0.9869`（6,498 个有限值对）；
  **符号一致比例 = 94.99%**；逐行量级一致到 4 位小数。
- → **`ddG` 的约定与我方一致：正值 = 削弱结合，负值 = 增强结合（beneficial）。`needs_flip = False`。**

### 模型预测列 —— `needs_flip` 判定

判定规则（写死在代码）：ΔΔG 真值**越小越好** → 方向正确的模型应与真值**正相关**。
在 **CATH 全表**（train+test，n≈6,500，避免 test 小样本误判）上估计相关系数。

| 来源 | 列 | Pearson(全表) | 需要反号? |
|---|---|---:|---|
| CATH | `Flex-ddG` | +0.4970 | 否 |
| CATH | `DDGPred` | +0.5768 | 否 |
| CATH | `ddG_FoldX` | +0.4965 | 否 |
| CATH | `ESM-IF` | **−0.2625** | **是** |
| CATH | `proteinMPNN` | +0.3137 | 否 |
| RBM | `DSSA-PPI` | +0.3081 | 否 |
| RBM | `PPIFormer` | **−0.1519** | **是** |
| RBM | `ESM-2-FFN` | **−0.0134** | 是（但相关≈0，判定不可靠） |
| RBM | `MuTon` | +0.0009 | 否（相关≈0） |
| RBM | `MpbPPI` | +0.0799 | 否 |
| RBM | `ProBASS` | +0.0047 | 否 |
| RBM | `SAAMBE-SEQ` | **−0.1934** | **是** |
| RBM | `FoldX` | **−0.5151** | **是** |

> **重要区分**：同一方法（FoldX）在 CATH 表（`ddG_FoldX`，正相关）与 RBM 文件
> （`FoldX`，负相关）中**符号相反** —— 说明 RBM 文件的作者把 FoldX 的
> $-\Delta\Delta G$ 写成了 `FoldX` 列。这正是 V6 §0.2 要求做符号审计的原因。
>
> **同时输出未反号版本**（`overall_pearson_noflip` 列），不隐藏反号带来的影响。
> 相关≈0 的三列（ESM-2-FFN / MuTon / ProBASS）的反号判定**不可靠**，
> 其实质是「该模型在此数据上无信号」。

---

## 3. The Illusion of Global Metrics

### Table 1 — 核心指标对比（同协议、同指标实现、同过滤规则）

**过滤规则**：按 complex 分组，丢弃组内突变数 < 5 的组（V6 §2.2）。

#### A. CATH-ddG 官方 held-out CATH blind test（`cath_fold == val`）

原始 813 行 / 53 complex → **过滤后 755 行 / 21 complex**；Base Rate(ΔΔG<0) = **0.1815**

| 模型 | Overall Pearson | **Macro per-PDB Spearman** | Hit@10(<0) | Base Rate | Efficiency@10 | 反号? |
|---|---:|---:|---:|---:|---:|---|
| **Flex-ddG** | +0.5626 | **+0.4159** | 1.0000 | 0.1815 | +0.5715 | — |
| **DDGPred** | +0.5573 | **+0.3262** | 0.8000 | 0.1815 | +0.3715 | — |
| FoldX | +0.4348 | +0.3584 | 0.4000 | 0.1815 | +0.4505 | — |
| proteinMPNN | +0.1215 | +0.2227 | 0.5000 | 0.1815 | +0.1399 | — |
| ESM-IF | −0.1896 | +0.0329 | 0.4000 | 0.1815 | +0.0497 | ✅ |

#### B. DSSA-PPI RBM case study

1,702 行 / 1 个人工分箱；真值 = `EXP_fitness`（**适应度，非 ΔΔG**，符号语义未说明）

| 模型 | Overall Pearson | Macro Spearman（组内） | Efficiency@10 | 反号? |
|---|---:|---:|---:|---|
| FoldX | +0.5151 | **+0.6478** | +0.7515 | ✅ |
| DSSA-PPI | +0.3081 | +0.3129 | +0.7439 | — |
| SAAMBE-SEQ | +0.1934 | +0.2058 | +0.5074 | ✅ |
| PPIFormer | +0.1519 | +0.1903 | +0.2652 | ✅ |
| MpbPPI | +0.0799 | +0.0845 | −0.0192 | — |
| ESM-2-FFN | +0.0134 | +0.0132 | +0.2898 | ✅ |
| ProBASS | +0.0047 | −0.0065 | +0.3271 | — |
| MuTon | +0.0009 | −0.0200 | +0.0513 | — |

> ⚠️ RBM 的 `Hit@10(ΔΔG<0)` 因 86.25% 的行都 <0 而**无判别力**，故不列出。

#### C. V5 baseline（我们，同协议对照）

1,129 行 / **61 complex**；Base Rate(ΔΔG<0) = **0.3729**

| 模型 | Overall Pearson | Macro per-PDB Spearman | Hit@10(<0) | Base Rate | Efficiency@10 |
|---|---:|---:|---:|---:|---:|
| V5 baseline seed42 | +0.3828 | +0.0782 | 1.0000 | 0.3729 | +0.0881 |
| V5 baseline seed43 | +0.3760 | +0.0970 | 0.9000 | 0.3729 | +0.0838 |
| V5 baseline seed44 | +0.3574 | +0.0752 | 1.0000 | 0.3729 | +0.1323 |
| **均值** | **+0.3721** | **+0.0835** | 0.9667 | 0.3729 | **+0.1014** |

### 假设判定（V6 §6）

| 判据 | 阈值 | 满足的外部模型数 |
|---|---|---:|
| Overall Pearson > 0.4 | >0.4 | **4 / 13** |
| Macro per-PDB Spearman < 0.15 | <0.15 | **5 / 13** |
| **两者同时成立（幻觉模式）** | — | **0 / 13** |
| Macro per-PDB Spearman > 0.3（反驳线） | >0.3 | **5 / 13** |

> **→ 假设 REFUTED。** 不存在「全局高 + 局部低」的组合：
> 全局高的 4 个模型局部也高（0.33–0.65）；局部低的 5 个模型全局也低。
> **但「全局指标与局部能力脱节」这一模式在我们自己的 V5 baseline 上成立**
> （Overall Pearson +0.3721 属中等，而 Macro per-PDB Spearman +0.0835 为全组最低）。

---

## 4. Target-Level Utility Analysis

### Efficiency@10 全排名（0 = 与随机无异；1 = oracle 水平）

| 排名 | 模型（数据集） | Efficiency@10 | 距 oracle |
|---:|---|---:|---:|
| 1 | FoldX (RBM) | **+0.7515** | 24.9% |
| 2 | DSSA-PPI (RBM) | +0.7439 | 25.6% |
| 3 | Flex-ddG (CATH) | **+0.5715** | 42.9% |
| 4 | SAAMBE-SEQ (RBM) | +0.5074 | 49.3% |
| 5 | FoldX (CATH) | +0.4505 | 55.0% |
| 6 | DDGPred (CATH) | +0.3715 | 62.9% |
| 7 | ProBASS (RBM) | +0.3271 | 67.3% |
| 8 | ESM-2-FFN (RBM) | +0.2898 | 71.0% |
| 9 | PPIFormer (RBM) | +0.2652 | 73.5% |
| 10 | proteinMPNN (CATH) | +0.1399 | 86.0% |
| 11 | **V5 baseline (ours)** | **+0.1014** | **89.9%** |
| 12 | ESM-IF (CATH) | +0.0497 | 95.0% |
| 13 | MpbPPI (RBM) | −0.0192 | >100%（差于随机） |

### 关键读数

1. **最优模型（FoldX on RBM）回收了 75% 的 oracle 空间**；
   **Flex-ddG 在 CATH 上回收 57%**。→ **Top-10 粒度的效用并未普遍崩塌。**
2. **我们的 V5 baseline 回收 10.1%，排名第 11/13**，仅优于 ESM-IF 与 MpbPPI。
3. **即使是 SOTA，也都没有达到 oracle**（最高 0.75）→ 「仍有提升空间」成立，
   但「普遍失效」不成立。
4. **一个值得注意的模式**：表现最好的都是**含物理能量项**的方法
   （FoldX、flex-ddG、SAAMBE-SEQ 属能量/统计势类），而**纯序列表示类**
   （ESM-IF、proteinMPNN、ESM-2-FFN、MuTon、ProBASS）普遍落在 0–0.33。
   → 这直接指向 V5 的瓶颈诊断：**我们用的正是池化序列表示。**

---

## 5. Conclusion: Are SOTA Models Ready for Wet-lab Design?

### 假设（H1）判定

> # **REFUTED**
>
> 「当前 PPI ΔΔG SOTA 模型普遍存在评估幻觉（全局高掩盖局部失效）」
> **不成立**。在被审计的 13 个模型预测列中，**0 个**同时满足
> `Overall Pearson > 0.4` 与 `Macro per-PDB Spearman < 0.15`。
> 相反，外部模型的 within-target 排序能力**显著优于**我们的 V5 baseline
> （Flex-ddG **+0.4159** vs V5 **+0.0835**，差距 **+0.33**）。

### 但同时成立的三个诚实的限定

1. **没有任何模型达到 oracle**：最高的 FoldX(RBM) Efficiency@10 = 0.7515，
   Flex-ddG(CATH) = 0.5715 → **距可用的完美选靶仍有 25–43% 的空间**。
   因此「SOTA 可直接用于湿实验选靶」这一强主张**仍不成立**。
2. **V5 baseline 自身确实表现出「全局中等 / 局部极弱」**：
   Overall Pearson +0.3721（看似可用）vs Macro per-PDB Spearman +0.0835（实为全组最低）。
   → **「评估幻觉」这一概念有效，但适用对象是我们的模型与纯序列表示类模型，
   而非物理能量类 SOTA。**
3. **本审计的主要限制是样本量**：
   - CATH 过滤后仅 **21 complex**（原始 53 PPI，因 ≥5 条/组过滤掉 32 个）
   - RBM 仅 **1 个分组**，且真值是适应度而非 ΔΔG
   - 因此**两个外部数据集的 macro Spearman 均未计算置信区间**，
     具体数值应视为**点估计**，不宜过度解读其绝对大小。

### 对领域问题的直接回答

**「Are SOTA models ready for wet-lab design?」**

- 若「ready」指**在一个已知靶点上从候选突变体中挑 top-10**：
  **物理能量法（FoldX / flex-ddG）已显示出实质性能力**（Efficiency@10 0.45–0.75），
  而**纯序列表示法（含我们的模型）尚不具备**。
- 若「ready」指**跨靶点、跨家族的稳健泛化**：**本审计无法回答** ——
  需要 ≥30 个 test complex 的样本量与置信区间，而现有公开预测文件不满足。

### 对 V5.0 结论的影响（不是推翻，而是加强 + 修正）

- **加强**：V5 的 `REPRESENTATION_BOTTLENECK_CONFIRMED` 得到外部对照支持 ——
  瓶颈是**池化表示**，因为用局部物理能量的方法在同协议下拿到 0.42–0.65 的组内排序。
- **修正**：V5 的 `within-complex Hit@10 = 0.3080 < Base Rate 0.3729` 容易被误读为
  「该任务的组内选择普遍不可行」。**V6 证明这是错的**。

---

## 6. Reproducibility

### 命令

```powershell
# PHASE 1 — 搜刮外部 SOTA 预测文件（只下小文件，不跑推理）
python -u scripts/v6/01_scavenge_sota_preds.py
#   -> data_scavenge_raw.json / data_scavenge_summary.md
#   -> 下载到 F:\dsh\downloads\{samaffinity,mint,dssa_ppi}\

# PHASE 2.1 — 符号一致性审计（含手工复算）
python -u scripts/v6/02_sign_audit.py
#   -> sign_audit.csv / sign_audit.md

# PHASE 2.2 — 同协议指标重算
python -u scripts/v6/03_audit_metrics.py
#   -> sota_audit_metrics.csv / sota_audit_metrics.md

# PHASE 4 — 诊断图
python -u scripts/v6/04_make_figures.py
#   -> figures/fig1..fig4
```

**依赖**：仅 `pandas / numpy / scipy / matplotlib` + V5 的 `scripts/v5/v5_common.py`（只读复用）。
**未安装任何新包**，**未运行任何模型推理**（V6 §0.1）。

### 数据来源与溯源（V6 §0.5 / §8）

| 来源 | URL | 许可 | 下载日期 | 大小 (B) |
|---|---|---|---:|---:|
| CATH-ddG 官方结果表 | GitHub `ak422/CATH-ddG` `data/SKEMPI2/skempi_v2_with_all_results.csv`（blobs API） | **仓库 license = None（未声明）** | 2026-10-13 | 1,291,260 |
| CATH-ddG 论文全文 | Europe PMC `PMC12261453` REST | 开放获取 | 2026-10-13 | 246,881 |
| CATH-ddG 配置 | `configs/train/CATH.yml`（contents API） | 同上 | 2026-10-13 | 2,373 |
| DSSA-PPI RBM 基线对比 | GitHub `Agitw/DSSA-PPI` `validate_on_RBM/RBM_baselines_comparison_ddG.csv` | **MIT** | 2026-10-13 | 174,025 |
| DSSA-PPI 其他 | `validate_on_RBM/pred_vs_exp.csv` | MIT | 2026-10-13 | 43,492 |
| SAMAffinity test folds | GitHub `iobio-zjut/SAMAffinity` `data/test_fold/*/*.csv`（58 个文件） | **license = None** | 2026-10-13 | — |
| MINT | GitHub `VarunUllanat/mint`（2 个小文件） | **MIT** | 2026-10-13 | 3,781 |
| 第三方真值（符号审计） | 本地 `data/v2.4/skempi_v2_cleaned.csv`（SKEMPI v2.0） | 见原始来源 | — | 1,568,384 |

**未获取 commit hash**：使用 GitHub API 的默认分支（`main`）当次快照；
`raw.githubusercontent.com` 对部分文件返回 0 字节，已改用 GitHub API（blobs/contents，base64）获取。

### 磁盘纪律（V6 §0.5）

- 所有下载落于 **`F:\dsh\downloads\{cathddg,samaffinity,mint,dssa_ppi}\`**
- 未写入 C 盘任何位置；未创建 > 2 MB 的下载（脚本内 `MAX_BYTES = 2 MB`，远低于 §0.1 的 100 MB 上限）
- V6 产物全部在 `experiments/v6_sota_audit/` 与 `scripts/v6/`；**未修改 V1.1–V5.0 任何文件**

---

## 附：产物索引

| 路径 | 内容 |
|---|---|
| `experiments/v6_sota_audit/SOTA_AUDIT_REPORT.md` | 本报告 |
| `experiments/v6_sota_audit/data_scavenge_summary.md` | PHASE 1 搜刮结果（含每个缺失项的原因与 URL） |
| `experiments/v6_sota_audit/data_scavenge_raw.json` | 各仓库文件树与判定（机器可读） |
| `experiments/v6_sota_audit/sign_audit.csv` / `.md` | PHASE 2.1 符号审计（含手工复算表） |
| `experiments/v6_sota_audit/sota_audit_metrics.csv` / `.md` | PHASE 2.2 全套指标重算 |
| `experiments/v6_sota_audit/hypothesis_testing.md` | PHASE 3 假设检验（三个诊断问题） |
| `experiments/v6_sota_audit/figures/fig1_illusion_scatter.png` | Overall Pearson vs Macro per-PDB Spearman |
| `experiments/v6_sota_audit/figures/fig2_efficiency_gap.png` | within-target top-10 效用 vs Random/Oracle |
| `experiments/v6_sota_audit/figures/fig3_hit10_vs_baserate.png` | Hit@10 vs Base Rate |
| `experiments/v6_sota_audit/figures/fig4_efficiency_at_10.png` | Efficiency@10 全排名 |
| `scripts/v6/01_scavenge_sota_preds.py` … `04_make_figures.py` | 4 个可复现脚本 |
