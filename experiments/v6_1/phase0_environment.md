# V6.1+ PHASE 0 — 环境检查与工作目录

**时间**：2026-10-13
**状态**：DONE

---

## 1. 环境

| 项 | 值 |
|---|---|
| OS | Windows-10-10.0.19045-SP0 |
| Python | 3.11.9（`E:\sjtu\download\python.exe`） |
| CPU | Intel i7-2630QM，4 核 / 8 线程（**CPU-only，无 GPU**） |
| RAM | 7.55 GB 总量（`psutil` 未安装，用 OS 计数器读取；实测空闲常 1.5–2.5 GB） |
| 磁盘 F: | **59.1 GB 可用** |
| 磁盘 C: | 13 GB 可用（**按磁盘纪律，本任务不写 C 盘**） |

### 依赖

| 包 | 版本 |
|---|---|
| numpy | 2.4.6 |
| pandas | 3.0.6 |
| scipy | 1.17.1 |
| sklearn | 1.9.1 |
| lightgbm | 4.7.0 |
| matplotlib | 3.11.2 |
| torch | 2.5.1+cpu（本次不使用） |

> **未安装任何新包**（V6.1 约束 3）。

---

## 2. 关键输入文件检查

| 文件 | 状态 | 大小 (B) |
|---|---|---:|
| `experiments/v6_sota_audit/standardized_preds.csv` | ❌ **MISSING** | 0 |
| `data/v4.0/ddg_split.json` | ✅ | 82,187 |
| `experiments/v5_0/p0_benchmark/baseline_repro_predictions.npz` | ✅ | 62,577 |
| `features/v2.2_proteinea_features.csv` | ✅ | 11,852,893 |
| `experiments/v6_sota_audit/sota_audit_metrics.csv` | ✅ | 4,564 |
| `experiments/v5_0/p4_conformal/p4_predictions.npz` | ✅ | 165,225 |
| `data/v4.0/ddg_pairs.csv` | ✅ | 2,238,640 |
| `features/v2_esm2_cache.npz` | ✅ | 52,481,097 |

### ⚠️ 一处规格与实际不符（如实记录，不假设）

V6.1 PHASE 1 的规格引用 `experiments/v6_sota_audit/standardized_preds.csv`
（`dataset_name == 'CATH'`）。**该文件在 V6.0 中并未生成** —— V6.0 实际产出的是
`sota_audit_metrics.csv`（已汇总指标）与 `sign_audit.csv`（符号审计），
而**逐条 (complex, mutation, true, pred) 的长表没有落盘**。

**处置**：在 PHASE 1 中**先补建该长表**（`scripts/v6_1/00_build_standardized_preds.py`），
数据源为已核实的 `F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv`
与 `F:\dsh\downloads\dssa_ppi\RBM_baselines_comparison_ddG.csv`，
符号按 V6.0 `sign_audit.csv` 的 `needs_flip` 统一。**不修改 V6.0 任何文件。**

---

## 3. 本次运行的时间预算（V6.1 规格）

| Phase | 任务 | 预算 | 状态 |
|---|---|---|---|
| 0 | 环境检查 | 5 min | **DONE** |
| 1 | CI + 交集 + 类别统计 | 45 min | IN_PROGRESS |
| 2 | Split Inflation 量化 | 90 min | PENDING |
| 3 | Pairwise 对比 | 30 min | PENDING |
| 4 | V5 vs CATH 分布 | 30 min | PENDING |
| 5 | 论文级 Figure | 30 min | PENDING |
| 6 | 论文大纲 | 30 min | PENDING |
| 7 | 工具包 + 找模型 | 30 min | PENDING |
| 8 | 汇总报告 | 15 min | PENDING |

---

## 4. 已知约束与风险（提前声明）

1. **CATH 过滤后仅 21 个 complex**（原始 53 PPI，因「组内 ≥5 条」过滤掉 32 个）
   → 所有 CI 都会很宽，这是**数据规模硬限制**，不是实现问题。
2. **RBM 只有 1 个人工分箱** → 该数据集**无法计算 CI**，只能作趋势性支持。
3. **V5 与 CATH 的数据集不同**（PPB-Affinity 五源 vs SKEMPI v2.0），
   交集可能很小甚至为 0 → PHASE 4 只做**分布对比**，不做同集对比。
4. **无 GPU**，PHASE 2 的 4 档划分 × 3 seeds × 多指标必须在 LightGBM 固定树数下完成，
   已按 90 分钟预算规划；超时则至少完成 L1 与 L4（最极端两档）。

---

## 5. 目录结构（已创建）

```
experiments/v6_1/
  README.md
  PROGRESS.md
  phase1_ci/
  phase2_split_inflation/
  phase3_pairwise/
  phase4_v5_comparison/
  phase5_figures/
  phase6_outline/
  phase7_toolkit/
  phase8_summary/
scripts/v6_1/
logs/v6_1/
```
