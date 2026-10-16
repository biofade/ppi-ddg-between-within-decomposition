# V7.2R PHASE 0 — 冻结并登记现有 V7.2 产物（TRACK_A_SMALL_PLM）

- 时间: 2026-10-14 18:33:59
- **不删除、不重算**；只在完成完整性审计后重新登记（§31）

## 0.1 现有产物清单

| 类别 | 路径 | 大小 | 标识 |
|---|---|---:|---|
| v72_report | `experiments\v7_2\phase0_provenance\model_decision.md` | 2,354 | `d398653d5dc00558` |
| v72_report | `experiments\v7_2\phase1_embeddings\embedding_report.md` | 654 | `3c19fe29f9b51ed2` |
| v72_manifest | `features\v7_2\residue_embedding_manifest.csv` | 147 | `ff2bdff5f52b5215` |
| v72_residue_embeddings | `features/v7_2/per_target/*.npz` | 1,105,582 | `10_files` |
| model_cache | `F:\dsh\cache\huggingface` | 31,385,421 | `esm2_t6_8M_UR50D` |

## 0.2 残基嵌入状态

- `per_target/*.npz`：**10** 个 target 文件
- 累计 rows：**1414** / 4473（**31.6%**）
- 模型：`facebook/esm2_t6_8M_UR50D`（7.84M 参数，hidden 320，**MIT**）
- 缓存：`F:\dsh\cache\huggingface`（**31.4 MB**，未写 C 盘）

## 0.3 完整性审计（§5 三项否决条件）

| 判据 | 结果 |
|---|---|
| target namespace 是否错误 | **NO** —— V7.2 的嵌入脚本用 `canonical_target_labels.csv` 的 `target_id_canonical` 作为 target 名；按 `row_id` 联合（§3.5） |
| 是否发生 outer 泄漏 | **NO** —— 嵌入是**与标签无关的通用 PLM 前向**，不含任何 fold 信息 |
| 测试标签是否进入特征构造 | **NO** —— F0/F1/F2 只用 mutation 字符串与WT/MUT 序列；无 ΔΔG |

⇒ **现有 V7.2 产物有效**，登记为 **`TRACK_A_SMALL_PLM`**。

### 可复用性

| 产物 | reusable | 说明 |
|---|---|---|
| 残基嵌入（per_target npz） | **YES** | 与标签无关，任何 fold 均可直接用 |
| 模型缓存（ESM-2 8M） | **YES** | MIT，F 盘 |
| F0/F1/F2 特征构造代码 | **YES** | 无标签泄漏 |
| `canonical_target_labels.csv` | **YES**（权威） | V7.1.1 冻结 |
| **R0–R7 ranker 结果** | **NOT YET** | V7.2 的 02 脚本尚未跑完（嵌入覆盖不足） |

