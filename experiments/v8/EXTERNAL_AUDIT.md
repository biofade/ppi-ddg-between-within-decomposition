# V8 External Model Audit（§14/§49/§72）— PHASE 1

> **状态：PARTIAL — 两个外部模型已核实，其中 HERMES 已构建匹配宇宙**
> **不臆造**：未核实者记 `NOT_FOUND` / `NOT_CHECKED`，与「已核实不存在」严格区分。

---

## A.0 ★★ HERMES — 官方发布行级预测，单点匹配率 **98.90%**

**§9 reuse-first 的首选资源类型：官方 row-level 预测 ⇒ 无需推理、无需 checkpoint、无需 GPU。**

| 项 | 值 |
|---|---|
| official_repository | https://github.com/StatPhysBio/hermes ＋ **预测仓库** https://github.com/StatPhysBio/hermes_experiments |
| 论文 | bioRxiv [2024.07.09.602403](https://www.biorxiv.org/content/10.1101/2024.07.09.602403v3)（**preprint，未同行评审**）|
| license | **MIT**（repo）；Zenodo 权重 **CC-BY-4.0** |
| **official_predictions（行级）** | ✅ **YES** —— `skempi/results/hermes_py_000_ft_cdna117k_ddg_st/skempi_v2_cleaned_NO_1KBH-hermes_py_000_ft_cdna117k_ddg_st-use_mt_structure=0.csv`（**403,691 B / 6,993 行**）|
| 预测列 | `PDB_filename, mutant_chain, mutant, Mutation(s)_cleaned, #Pdb, log_proba_mt__minus__log_proba_wt` |
| **输出语义** | ⚠️ **不是 ΔΔG，是 log-odds**（`log P(mt) − log P(wt)`）；**越大 = 越有利**，而 V8 约定 **ΔΔG>0 = 结合变弱** ⇒ **方向相反**，须取负：`hermes_proxy = −(log_proba_mt − log_proba_wt)` |
| **GPU** | ✅ **完全不需要**（无需推理）|
| training_data | 预训练 ProteinNet CASP12 30% ID；微调 cDNA117k（稳定性）＋ **SKEMPI v2.0（结合）** ＋ RaSP ＋ Megascale |
| **§50 污染状态** | ⚠️ SKEMPI v2.0 微调 ⇒ **不可自动进 Tier 1**，定为 **TIER 2** |

### 已完成的匹配审计（`32c_hermes_final.py`，严格 1:1，4473 行无膨胀）

```
canonical 行                    = 4473   （断言校验：无行膨胀）
命中 HERMES                     = 2865   (64.05%)
仅单点突变子集  n=2897  命中     = 2865   (98.90%)   ← ★
覆盖 target                     = 195 / 214
覆盖 family                     = 101 / 108
命中子集中 n>=5 的 target        = 94                ← 可算 Macro 指标
方向校验 r(y_true, hermes_proxy) = +0.1688            ← **为正 ⇒ 符号约定正确**
```

**⇒ 结论：**
1. **符号约定经数据验证**（r 为正），`hermes_proxy` 可直接使用。
2. **单点突变匹配率 98.90%** ⇒ 可构造 **`U_HERMES_CANONICAL`**（2865 行 / 195 target / 94 个 n≥5）。
   **这是 V8 的第一个外部匹配宇宙**，使「外部模型是否复现能力分化」首次可检验。
3. **限制**：多点突变行几乎全未命中（HERMES `mutant` 列以单点为主）
   ⇒ **外部对比必须限定在单点子集**。

产出：`external/hermes_matched_final.csv`、`external/hermes_matched_meta.json`、
`external/hermes_audit.json`；下载件在 `F:\dsh\downloads\hermes\`。

---

## A.1 ★ CATH-ddG — 与 V8 任务设定**最接近**的外部工作

| 项 | 值 |
|---|---|
| 论文 | Yu, Bi, Ma, Li, Wang (2025). *CATH-ddG: Towards robust mutation effect prediction on PPI out of CATH homologous superfamily.* **Bioinformatics 41(Suppl_1): i362–i372**, ISMB/ECCB 2025. DOI [10.1093/bioinformatics/btaf228](https://doi.org/10.1093/bioinformatics/btaf228) |
| repo / license | https://github.com/ak422/CATH-ddG / **CC BY 4.0** |
| **task setting** | **"new training/test splits that mitigate data leakage according to the CATH homologous superfamily"** ⇒ 与 V8 的 task-specific family-OOD 同类 |
| **污染认知** | 摘要：*"deep learning models may **overestimate performance due to potential data leakage**"* ⇒ 与 V8 的 H4 同一关切 |
| **更严格的难度定义** | *"a **hard mutation** … maximum **TM-score < 0.6** compared to the training set"* ⇒ **比 V8 的家族划分更细**（结构级） |
| checkpoints | ✅ 3 组 Google Drive 权重：**held-out CATH**、**held-out PPIFORMER**、case study |
| 推理脚本 | ✅ `test_DDAffinity.py <config> --device <dev>`；`--output_results` 输出**行级** CSV |
| **官方行级预测文件** | ❌ **未发布**（须自行推理）|
| **cpu_inference_feasible** | ⚠️ **很可能可以**：`--device` 默认 `'cuda'` 但代码用 `torch.load(..., map_location=args.device)` 与 `recursive_to(batch, args.device)` ⇒ device-agnostic。**未实跑验证** ⇒ `LIKELY_YES_UNVERIFIED` |
| per-protein 指标 | ✅ 含 `eval_perprotein_modes` / `eval_skempi_three_modes` ⇒ 他们也做 per-protein 评价 |
| training_data | SKEMPI v2.0 + 其自有 CATH 超家族划分 |

---

## A.2 PPIformer

| 项 | 值 |
|---|---|
| 论文 / repo / license | ICLR 2024 [arXiv:2310.18515](https://arxiv.org/abs/2310.18515) / https://github.com/anton-bushuiev/PPIformer / **MIT** |
| checkpoints | ✅ Zenodo `weights.zip`：`ddg_regression/{0,1,2}.ckpt`（3 模型 ensemble）|
| 推理 | ✅ `predict_ddg(models, ppi_path, muts)`；输入 `(PDB/PPIRef 文件, 突变列表)` |
| **cpu_inference_feasible** | ✅ **README 原文即写 `device = 'cuda' if torch.cuda.is_available() else 'cpu'`** |
| training_data | 预训练 PPIRef；**微调 SKEMPI v2.0** |
| **§50** | ⚠️ SKEMPI 微调 ⇒ 不可自动进 Tier 1；定为 **TIER 2** |
| 归属注记 | 作者群即 [R1] 作者，**自采界面相似度非泄漏划分** ⇒ 属「非泄漏基准」阵营 |

---

## A.3 [R1] 泄漏量化论文（非模型，但直接约束 V8 新颖性）

**Bushuiev et al. 2024, "Revealing data leakage in protein interaction benchmarks"**
[arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)

| 划分策略 | 结构性泄漏率（iDist 实测）|
|---|---:|
| PPI-code 划分 | **86%** |
| PDB-code 划分 | **65%** |
| PDB 提交时间 | 61% |
| **序列相似度（MMseqs2, 30% identity）** | **30%** |
| Linclust | 37% |
| **SKEMPI v2.0 + PPI-code 划分** | **56%** |
| **SKEMPI v2.0 + 领域专家分组** | **0%** |

**⇒ 对 V8 的三条硬约束：**
1. 「首次发现宽松划分膨胀性能」**已被占据** ⇒ 禁用该声明。
2. V8 的 `L3_seq40` 正落在 R1 证明**仍有 30% 泄漏**的那一类 ⇒ 须承认。
3. R1 推荐的**界面 3D 相似度划分，V8 未实现** ⇒ 明确局限。

---

## A.4 其他已核实（子代理审计，摘要）

| method | license | 官方行级预测 | checkpoint | CPU | §50 污染 |
|---|---|---|---|---|---|
| **RDE-Network** | **Apache-2.0**（可商用）| ❌ 未发布（须推理）| ✅ Google Drive `DDG_RDE_Network_30k.pt` | ⚠️ `--device` 可传 cpu，未验证 | **SKEMPI v2.0 训练** |
| **Prompt-DDG** | **MIT** | ❌ 未发布 | ✅ Google Drive | ⚠️ 代码级 `cuda if available else cpu` | **SKEMPI v2.0 训练** |
| **FoldX** | 专有；学术免费但**需注册登录**，仅编译二进制，**time-limited** | ❌ | N/A | ✅ 原生 CPU-only | N/A（经验势）|
| **PreFold-dG** | 论文确认（Bioinformatics 2026, [btag489](https://academic.oup.com/bioinformatics/article-pdf/42/Supplement_2/btag489/70726266/btag489.pdf)）| `NOT_CHECKED` | `NOT_CHECKED` | `NOT_CHECKED` | `NOT_CHECKED` |
| DDGPred / ProteinMPNN / ESM-IF / BA-DDG / EBM-DDG / Light-DDG / Pythia-PPI / H3-DDG / MARS-DDG / USP-ddG / ReCLIP | — | `NOT_CHECKED` | — | — | — |

**⇒ FoldX 不可再分发**（须注册 + 限时二进制）⇒ **不适合作为可复现基准的常规项**。

---

## B. GPU 判定（§51/§52）— PHASE 13 定稿

```
GPU_JUSTIFICATION
-----------------
最强候选 = CATH-ddG
COND_1 official row-level prediction unavailable = YES
COND_2 official checkpoint exists                = YES
COND_3 official inference code exists            = YES
COND_4 CPU impractical                           = NO  (device-parametric; LIKELY_YES_UNVERIFIED)
=> GPU_REQUIRED = NO        (§52 要求四条件全部成立)

附加事实：HERMES 提供**官方行级预测** ⇒ 连推理都不需要 ⇒ GPU 更无必要。
```

**结论：不生成 Kaggle notebook。** 全部 V8 结论来自已有内部预测 + 官方外部预测 + CPU 统计。

---

## C. 对 V8 结论的影响

1. **H4 新颖性被显著削减**：CATH-ddG 与 [R1] 均已处理「泄漏导致高估」。
   **V8 仍独有**的是**三轴能力分解**在该问题上的应用 ——
   `C_B` 从 **+0.43 塌到 ≈0** 的协方差机制，以及
   「Overall Pearson 膨胀 **4.5×** 而 Efficiency@10 只膨胀 **2.5%**」的**不均匀性**。
2. **§68 `STRONGLY_SUPPORTED` 的条件现在**首次可检验**：
   已建成 `U_HERMES_CANONICAL`（2865 行 / 195 target / 94 个 n≥5）
   ⇒ 下一步用 V8 的 matched-universe 口径对 HERMES 重算三轴，
   检验**能力分化是否外部可复现**。
   **尚未执行** ⇒ 仍记 `EXTERNAL_MODEL_REPRODUCTION = NOT_RUN`。
3. **本审计本身是有价值的交付**：给出了三条可复现的外部对比路径
   （HERMES 官方预测 / CATH-ddG 权重 / RDE-Network Apache-2.0）。
