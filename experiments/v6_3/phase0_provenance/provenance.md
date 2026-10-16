# V6.3 PHASE 0 — Inputs & Provenance

**时间**：2026-10-14
**约束**：未修改 V1.1–V6.2.1 任何文件（本任务**只读**上游产物）；`v6_3` 命名空间；
CPU only；**未下载任何文件**；未伪造数字

---

## 1. 冻结输入（SHA256，见 `input_manifest.csv`）

| 名称 | 路径 | 大小 (B) | SHA256（前 40） | 角色 |
|---|---|---:|---|---|
| `matched_analysis_set` | `experiments/v6_2/phase2_split4/matched_analysis_set.csv` | 754,492 | `8d209eea08cc6058fb235b059c9e3725fcbc9ba8` | V6.2 canonical 行宇宙（7523） |
| `split4_predictions` | `experiments/v6_2/phase2_split4/split4_predictions.npz` | 1,696,726 | `ba711ea1b349836ac18b44ba67f6ab2f0a0b680b` | **L4 M1/M0 预测（本审计 primary）** |
| `d_only_predictions` | `experiments/v6_2_1/phase3_d_only/d_only_predictions.npz` | 484,102 | `98b64c2b02140f232d993251c32d5a02daea7f2e` | V6.2.1 真 D-only M2 |
| `target_level_stats` | `experiments/v6_2_1/phase4_target_size/target_level_stats.csv` | 117,839 | `63c6d2a905ca2bcfca6adee8ad61b92d32b2838b` | V6.2.1 target 统计 |
| `primary_table` | `experiments/v6_2_1/phase5_corrected_stats/primary_table.csv` | 796 | `d02e6d79d9b23a6f53a367b9e1e6eb1aa847edf3` | V6.2.1 修正核心表 |
| `route_statistics` | `experiments/v6_2_1/phase5_corrected_stats/route_statistics.csv` | 1,077 | `b63cfe240567ddd56b746018a5d2d1c79e505f74` | V6.2.1 route stats |
| `ddg_pairs` | `data/v4.0/ddg_pairs.csv` | 2,238,640 | `e43b19d69a8b54efe5cc4b3e647e87657106e86c` | ΔΔG 配对（`n_mut`/`mutation`/assay/temp） |
| `ddg_split` | `data/v4.0/ddg_split.json` | 82,187 | `92d4e8d78f0778ff82f12896aeedc9607e5e8cf7` | V4 划分 + `family_of_complex` |
| `v22_features` | `features/v2.2_proteinea_features.csv` | 11,852,893 | `134a406175c975541ae40aed37c3f824ad151942` | **58 个数值结构描述符** + 序列/实验元数据 |
| `skempi_v2` | `data/v2.4/skempi_v2_cleaned.csv` | 1,568,384 | `72eaccc14c370bb61d5053cba06702e315d38195` | 第三方真值（符号审计用） |
| `esm2_cache` | `features/v2_esm2_cache.npz` | 52,481,097 | `2b6fdcc331dd5030c70ffff557d3cdbc37111a0b` | pooled ESM-2 嵌入 |
| `v5_common` | `scripts/v5/v5_common.py` | 27,748 | `91bee0814f66737b50edadd937d967644bbf3ec7` | 特征装配与超参配置源 |

## 2. 环境

| 项 | 值 |
|---|---|
| Python | 3.11.9 (`E:\sjtu\download\python.exe`) |
| CPU | Intel i7-2630QM, 4 核 / 8 线程 |
| RAM | 7.55 GB |
| GPU | **无**（CPU only，遵守约束 3） |
| 关键依赖 | numpy 2.4.6 / pandas 3.0.6 / scipy 1.17.1 / sklearn 1.9.1 |

## 3. 关键决策基线（来自 V6.2.1，本轮将检验/修正）

```
L4 family-OOD，M1 seed-ensemble:
  r_full = +0.1347
  C_B    = −0.0121     ← 本轮 PHASE 5 证明这是 weighting 产物
  C_W    = +0.1468
```

V6.2.1 分层（≥30 组 target-mean Pearson = −0.0008）是本轮的**待解释现象**。

## 4. 元数据可用性预判（决定哪些分析可做）

| 维度 | 可用性 |
|---|---|
| 结构/界面描述符 | ✅ 58 个数值列，非空率 100% |
| mutation 组成 | ✅ 可从 `ddg_pairs.mutation` 解析（本轮修复了单字母码解析 bug） |
| `n_mut`（每变体突变数） | ✅ 100% |
| assay 类型 | ✅ 11 类 |
| temperature | ⚠️ 部分（`struct_Temperature(K)` 68.5%） |
| **source / paper** | ❌ **全部 = SKEMPI v2.0（单一来源）→ leave-one-source-out `UNAVAILABLE`** |
| 实验条件（pH 等） | ❌ 不存在 → `UNAVAILABLE` |

## 5. 判定原则（V6.3 §1.1 / §4）

- **`LABEL_DERIVED`**：用到 true ΔΔG 的字段（`true_mean_ddg`、`ddg_std`、
  `beneficial_fraction`、`within_spearman` …）→ **仅可用于机制解释**，
  **不可**进入任何可部署模型
- **`DEPLOYABLE`**：预测前可知
- **⚠️ 本轮审计中发现的一处分类错误（已更正）**：PHASE 1 字典最初把
  `signed_mean_error_M0` / `pred_mean_ddg` / `pred_mean_M0` 标为 `DEPLOYABLE`，
  但它们实为**模型输出派生量**，不是可部署输入描述符。
  已在最终汇总中剔除并说明；PHASE 6 的特征白名单**从未包含**它们
  （见 `scripts/v6_3/04_set_calibration.py` 的 `TARGET_ONLY` / `SET_ONLY` 白名单）
- **只能说 association，不能说 causation**（本项目**没有** causal design）
