# V6.1+ — Overnight Statistical & Paper-Prep Marathon

**状态**：全部 9 个 PHASE 完成（PHASE 0–8），耗时约 70 分钟（预算 5 小时，未超时）
**上游**：V6.0 = `REFUTED`（「SOTA 评估幻觉」假设被证伪）

## PHASE 清单

| PHASE | 内容 | 状态 | 关键产物 |
|---|---|---|---|
| 0 | 环境检查 | DONE | `phase0_environment.md` |
| 1 | CI + 交集 + 类别统计 | DONE | `phase1_ci/`（**类别差异 p = 0.0038**） |
| 2 | Split Inflation | DONE | `phase2_split_inflation/`（**+0.2663 Pearson**） |
| 3 | Pairwise 对比 | DONE | `phase3_pairwise/`（4/10 显著） |
| 4 | V5 vs CATH 分布 | DONE | `phase4_v5_comparison/`（p < 0.0001） |
| 5 | 论文级 Figure | DONE | `phase5_figures/`（6 张 300 dpi） |
| 6 | 论文大纲 | DONE | `phase6_outline/paper_outline.md` |
| 7 | 工具包 + 找模型 | DONE | `phase7_toolkit/` |
| 8 | 汇总报告 | DONE | `V6_1_FINAL_SUMMARY.md` |

## 三条核心结论

1. **物理/能量类显著优于纯序列表示类**：per-complex Spearman
   **+0.3871 vs +0.1278**，配对 Δ **+0.2593** [+0.1262, +0.3944]，
   Wilcoxon **p = 0.0038**，效应量 r = 0.61（大），14/21 complex 方向一致。

2. **Split Inflation = +0.2663 Pearson（相对 +46.1%）**：同一模型，仅换划分
   （mutation-random → family-level）即可让 Overall Pearson 从 **+0.8436**
   降到 **+0.5774**。且**4 档的 within-target Spearman 全在 0 附近**
   （−0.034 ~ +0.057）→ 该模型的全局能力**完全来自跨复合物区分**。

3. **我们自己的模型属于「序列表示类」档位**：V5 median per-complex Spearman
   **+0.0939**，显著低于 CATH 物理类 **+0.4179**（Mann-Whitney p < 0.0001），
   而与 CATH 纯序列类 **+0.1314** 无显著差异（p = 0.5598）。

## 执行命令（全部可复现）

```powershell
python -u scripts/v6_1/00_build_standardized_preds.py
python -u scripts/v6_1/fastmetrics.py                      # 自检
python -u scripts/v6_1/01_cluster_bootstrap_ci.py
python -u scripts/v6_1/02_testset_intersection.py
python -u scripts/v6_1/03_category_paired_test.py
python -u scripts/v6_1/04_split_inflation.py --max-minutes 85
python -u scripts/v6_1/05_pairwise_comparison.py
python -u scripts/v6_1/06_v5_cath_distribution.py
python -u scripts/v6_1/08_make_paper_figures.py
python -u scripts/v6_1/07_search_external_models.py
python -u scripts/v6_1/09_zenodo_abag_within_target.py
python -u scripts/v6_1/10_final_verify.py
```

## 约束遵守

- ✅ 未修改 V1.1–V6.0 任何文件
- ✅ 产物全部在 `v6_1` 命名空间（含 `scripts/v6_1/`）
- ✅ CPU only，未安装任何 GPU 包
- ✅ 下载单文件最大 50.5 MB（< 100 MB 上限）；全部落于 F 盘，未写 C 盘
- ✅ 未跑任何模型推理
- ✅ **未伪造任何数字**；不可用的数据源已明确记录原因
- ✅ 每完成一个 PHASE 立即更新 `PROGRESS.md`

## 关键方法学说明

- **配对检验**（PHASE 1.3）：同一批 complex 的两类模型**必须配对**
  → 用 **Wilcoxon 符号秩**，**明确不用 Mann-Whitney U**（两组不独立）
- **不配对检验**（PHASE 4）：V5 与 CATH 是不同数据集、不同 complex
  → 此处**允许**用 Mann-Whitney U
- **cluster bootstrap**：抽样单位 = complex，携带其全部突变行；
  避免把同一靶点的多个突变体误当独立样本而低估 CI
- **性能优化口径不变**：`scripts/v6_1/fastmetrics.py` 的 `macro_spearman`
  用 `np.argsort` 双排序替代 `scipy.rankdata`，**与 `v5_common` 差异 0.00e+00**（快 17×）；
  `eff10` 的 bootstrap 内部随机抽样次数由 200 降为 20（point estimate 仍用 200，已标注）
