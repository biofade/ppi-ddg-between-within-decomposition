# V8 CR-PPI-ddG — FINAL PACKAGE README

**用途**：供最终科学审计、verdict 升级判断、novelty 确认与预印本准备。
**生成方式**：`python scripts/v8/45_build_final_package.py`
**未包含**（按要求）：模型权重、原始 PDB 全库、Python venv、Kaggle/Colab 环境、
167 个历史 prediction artifact、V1–V6 旧文件。

---

## 0. 一页状态

```text
UPSTREAM_HASH_CHANGE = 0            (24 项冻结，含 V7.3/V7.4/V7.4R2 制品)
GPU_USED = NO                       GPU_REQUIRED = NO
V7_4R2_NEEDED = OPTIONAL

H1 Global-vs-within discordance = WEAKLY_SUPPORTED
H2 Target-weighting sensitivity = PASS
H3 Selection alignment          = FAIL
H4 Family-OOD split inflation   = PASS
H5 Benchmark-guided remedy      = PASS

BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE          (§67 LEVEL 1)
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE             (6 独立方法族 < 8 门槛)
PUBLICATION_READINESS   = PARTIAL / CONDITIONAL
```

---

## 1. 审计者应先读的 5 个文件（按顺序）

| # | 文件 | 为什么先读 |
|---|---|---|
| 1 | `README.md`（本文件） | 状态与已知缺陷 |
| 2 | `POOL_PROVENANCE_AUDIT.md`（在 `reports/`） | **★ 最重要的方法学文件**：所有秩相关数字的池定义与可复核来源 |
| 3 | `V8_FINAL_SUMMARY.md` §22 | §85 机器可读执行块（28/28 字段） |
| 4 | `EXTERNAL_AUDIT.md` | PHASE 1 外部验证：已接入 / 已阻断 / 未核实 |
| 5 | `BENCHMARK_NOVELTY_MATRIX.md` §0 | **新颖性已被削减的部分（必须先读）** |

---

## 2. 已知缺陷与已修正项（**审计重点**）

### 2.1 已在本次审计中修正

| # | 问题 | 处置 |
|---|---|---|
| 1 | `fig1` 未标注其 160 模型池含 **51 个 V6.2 伪重复键** | **已重绘**：双面板 (a) ALL / (b) provenance-clean，图内加 CAVEAT 文字 |
| 2 | scenario 表的首选来自 **permissive split 训练**的模型 | **已补** `scen_scenario_recommendations_provenance_clean.csv`：另给 STRICT-only 首选 |
| 3 | 缺失池 provenance 审计 | **已补** `POOL_PROVENANCE_AUDIT.md` + `pool_provenance_audit.csv` |
| 4 | 缺 metric-engine 自检报告 | **已补** `bootstrap/METRIC_ENGINE_SELFTEST.md` |
| 5 | `rho(Overall r, Macro rho)` 被报出多个不同值 | **已定案**：见 §3，必须附 `pool_id` 引用 |

### 2.2 未完成 / 已知限制

| # | 项 | 状态 |
|---|---|---|
| 1 | 外部独立方法族仅 **6**（< 8 门槛）⇒ 能力分化在外部**不可检** | **NOT_TESTABLE**，有证据的阻断（见 `EXTERNAL_INFERENCE_FEASIBILITY.md`） |
| 2 | `target-cluster` bootstrap 的 **target 抽样**分量未完成 | 点估计已用规范实现逐位复核；仅有 **Fisher CI（模型池抽样）**。原因：`v8metrics` 完整调用 **14.5 s/迭代**（≈8 h / 2000 次） |
| 3 | PHASE 1 外部方法审计覆盖 **6 / 17** | 其余 11 个 `NOT_CHECKED` |
| 4 | `L2_complex` 与 `L3_seq40` 指标逐位相同 | 疑为上游同一掩码，**未核实**，故未声称两者独立 |
| 5 | 无独立外部数据集 | `EXTERNAL_DATASET_VALIDATION = NOT_RUN` |

### 2.3 本次工作中**被自己推翻**的结论（透明度记录）

| 轮次 | 曾主张 | 结果 |
|---|---|---|
| 8 | 分化由**多点突变**驱动 | **否定**（单点/多点同号，`MULTIPOINT_NOT_THE_DRIVER`） |
| 9 | V6.2 **极端点**机械制造负相关 | **否定**（极端点 0/37） |
| 10 | `+0.81 vs −0.03` 边界 | **一度撤回**（因独立实现符号相反） |
| 11 | 上述撤回 | **撤回本身是错的** —— 逐位诊断证明有 bug 的是新脚本 `36`，不是 `35`；原结论**恢复** |

**净结果**：`rho(Overall r, Macro rho)` 在 10 个池上的观测范围是 **−0.268 .. +0.700**；
**任何单值都不是 benchmark 的性质**。

---

## 3. ★ `rho(Overall r, Macro rho)` 的池依赖（审计核心）

| pool_id | n | PERMISSIVE | STRICT | **rho(r,macro)** | rho(r,eff@10) |
|---|---:|---:|---:|---:|---:|
| `P1_U_CANONICAL_all` | 160 | 7 | 139 | **+0.1837** | −0.0071 |
| `P2_U_CANONICAL_nrows4473` | 67 | 7 | 55 | **−0.2685** | +0.3029 |
| `P3_P2_minus_V6_2` | 16 | 0 | 11 | **+0.1861** | +0.1987 |
| `P4_U_CANONICAL_minus_V6_2` | 109 | 0 | 95 | **+0.4088** | +0.3085 |
| `P5_pool_composition_ALL` | 178 | 37 | 133 | **+0.2789** | +0.2722 |
| `P5a_..._PERMISSIVE` | 37 | 37 | 0 | **+0.6262** | **+0.8088** |
| `P5b_..._STRICT` | 133 | 0 | 133 | **+0.2214** | **−0.0289** |
| `P6_HERMES_variants` | 7 | 0 | 0 | **+0.6786** | +0.7143 |
| `P7_EXT_CANONICAL` | 5 | 0 | 0 | **+0.7000** | +0.7000 |
| `P8_EXT_HELDOUT_CATH` | 5 | 0 | 0 | **+0.7000** | +0.4000 |

**观测范围：−0.268 .. +0.700。** 完整定义与来源脚本见 `POOL_PROVENANCE_AUDIT.md`。

**唯一**在同一宇宙、同一行数下仅改池组成就翻号的对照：
`P2`（含 51 个 V6.2 键）= **−0.2685** → `P3`（去掉它们）= **+0.1861**。

> ⚠️ `P5a` 的 `rho(r,eff@10)=+0.8088` 与 `P5b` 的 `−0.0289` 是**本包内最值得审的对比**：
> 宽松划分训练的制度下全局指标几乎完美代理选择效用，严格制度下**完全失效**。
> 该对比的**点估计经规范实现逐位复核**（`discrepancy_diagnosis.csv`），
> 但**尚无 target-cluster 区间**。

---

## 4. 目录导览

```text
V8_FINAL_PACKAGE/
├── README.md                        本文件（状态 + 已知缺陷 + 池审计摘要）
├── V8_FINAL_SUMMARY.md              最终总报告（含 §22 §85 执行块，28/28 字段）
├── EXTERNAL_AUDIT.md                PHASE 1 外部审计
├── BENCHMARK_NOVELTY_MATRIX.md      新颖性矩阵（§0 为已核实的重叠削减）
├── model_registry.csv               211 预测键 / 12 制品 / 9 族
├── BENCHMARK_SCOREBOARD.csv         202 条模型×宇宙三轴读数
├── MANIFEST.json                    全部文件 + sha256 前 16 位
├── reports/                         21 份关键分析报告
├── external/                        PHASE 1 外部验证（匹配表 + 审计 + meta）
├── tables/                          T1..T8 + TABLE_INDEX
├── figures/                         fig1..fig7（fig1 已重绘）
├── bootstrap/                       统计验证（bootstrap / 逐宇宙 / 池审计 / 自检）
└── cards/                           7 张模型能力卡
```

---

## 5. PHASE 1 外部验证：实际接入了什么

### 5.1 成功接入（**零推理**，全部来自已发布预测）

| 方法 | 来源 | 宇宙 | 匹配率 | Overall r | Macro ρ |
|---|---|---|---:|---:|---:|
| `HERMES` py000/py050/bp000/bp050 | `hermes_experiments` 官方行级预测 | `U_HERMES_CANONICAL` | 98.9%（单点） | +0.100 ~ +0.129 | +0.063 ~ +0.081 |
| `proteinMPNN` 002/020/030 | 同上 | `U_HERMES_CANONICAL` | — | +0.158 ~ +0.184 | +0.080 ~ +0.095 |
| `DDGPred` | CATH-ddG 表 | `U_EXT_CANONICAL` | 84.6% | +0.3357 | +0.1047 |
| `ddG_FoldX` | CATH-ddG 表 | `U_EXT_CANONICAL` | 85.4% | +0.2727 | +0.1502 |
| `Flex-ddG` | CATH-ddG 表 | `U_EXT_CANONICAL` | 84.0% | +0.2619 | +0.1320 |
| `ESM-IF` | CATH-ddG 表 | `U_EXT_CANONICAL` | 85.4% | −0.1014 | −0.0346 |
| 上述 5 者 | CATH-ddG 官方 held-out | `U_EXT_HELDOUT_CATH`（813 行/53 PPI/21 n≥5） | 全部可用 | +0.14 ~ +0.61 | −0.033 ~ +0.416 |

**跨宇宙一致性**：两个 CATH-ddG 宇宙的方法排序 **ρ = +0.90**（n=5）。

### 5.2 阻断（有证据，未强行绕过）

* **PPIformer** —— 上游 `ppiformer/data/transforms.py:471` **无条件 `rel_pos.cuda()`**，
  作者以 `# TODO` 自承未实现 CPU 分支；要求 `torch>=2.7`（本机 2.5.1）；可用 RAM 1.9 GB。
* **RDE-Network / Prompt-DDG** —— 未发布逐突变预测，需自跑；Prompt-DDG 无任意
  (结构, 突变) 的评分 CLI。
* **CATH-ddG 自身** —— 88% 行属其训练集，且**未发布自身逐突变预测**（仅 per-protein 表）。
* `HER2.csv` / `S285.csv` —— **无基线预测列**。
* `prompt-fold0/1/2` 列 —— 原样为 `train`/`val` 字符串，**不是数值**。

### 5.3 每模型的审计覆盖

| 审计项 | 覆盖情况 |
|---|---|
| mapping audit | ✅ 键格式 `(PDB4, 位置, 突变)`；`discrepancy_diagnosis.csv` 记录逐位比对 |
| sign-convention audit | ✅ `hermes_proxy = −(log_proba_mt − log_proba_wt)`，方向校验 **r=+0.1688 > 0** |
| provenance / split audit | ⚠️ **部分**：HERMES/CATH-ddG 的 SKEMPI v2.0 微调已记录（TIER 2）；`Flex-ddG/DDGPred/FoldX` 的**原始训练集未核实** |

---

## 6. 场景表的 provenance-clean 版本（审计点 6）

原 `scenario_recommendations.csv` 的 A/C 场景首选取自 **`V6_2_split4`（宽松划分训练）**。

补充的 `scen_scenario_recommendations_provenance_clean.csv` 另给 **STRICT-only** 首选：

| 场景 | 主指标 | ALL 池首选（含宽松） | **STRICT-only 首选** | 是否改变 |
|---|---|---|---|---|
| A 定量 ΔΔG | `rmse` ↓ | `L1_random_A1_base_set` (PERMISSIVE) | **`F8_Full_s42`** | ✅ 改变 |
| B 文库排序 | `macro_spearman` ↑ | `R1_Base_F0_s42` (STRICT) | `R1_Base_F0_s42` | ❌ 不变 |
| C 湿实验短名单 | `efficiency@10` ↑ | `L1_random_A3_base_local` (PERMISSIVE) | **`p_P5_set_only_LSR`** | ✅ 改变 |

**⇒ 3 个场景中 2 个的推荐在 provenance-clean 后**改变**。**

---

## 7. 统计验证状态

| 项 | 状态 |
|---|---|
| `rank_agreement_bootstrap.csv` | ✅ **5000/5000 valid**，target-cluster，seed=42 |
| `metrics_decomposition_qc.csv` | ✅ 分解恒等式残差最大 **3.331e-16** |
| `METRIC_ENGINE_SELFTEST.md` | ✅ 14/14 PASS（合成数据四用例） |
| `discrepancy_diagnosis.csv` | ✅ 规范实现逐位复核（max \|Δ\| = 5.9e-3） |
| H1 子集结果（67 / 10 / 57 / 109 / 160） | ✅ `rankrev_*` 与 `POOL_PROVENANCE_AUDIT` |
| **target-cluster CI（规范实现）** | ❌ **未完成**（14.5 s/迭代 ⇒ ≈8 h）。仅有 Fisher CI（模型池抽样） |

---

## 8. 建议的审计顺序

1. 核 §2.3 的自我推翻记录是否被诚实反映在 `V8_FINAL_SUMMARY.md`；
2. 用 `POOL_PROVENANCE_AUDIT.md` 逐个复核你关心的 rho 数字；
3. 打开 `figures/fig1_*.png`，确认 CAVEAT 文字与双面板是否满足要求；
4. 打开 `bootstrap/scen_scenario_recommendations_provenance_clean.csv`，确认场景推荐；
5. 判断 `EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE` 是否应被接受为定稿结论，
   或要求补齐（需 ≥2 个额外独立外部方法）。
