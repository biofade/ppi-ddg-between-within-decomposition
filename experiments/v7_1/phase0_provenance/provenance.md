
# V7.1 PHASE 0 — Provenance

**时间**：2026-10-14
**约束**：未修改 V1–V7.0 任何文件（只读上游产物）；新命名空间 `v7_1`；
CPU only；**未下载任何文件**；未使用 target-ID embedding / family lookup / 测试标签；
所有预处理 fit on train only；未跨 row universe 比较；负结果全部保留。

## 冻结输入（14 个，SHA256 见 input_manifest.csv）

| 名称 | 角色 |
|---|---|
| `v62_matched_set` | canonical 行宇宙（7523） |
| `v70_block_L` / `v70_block_S` / `v70_block_G` | V7.0 三个特征块 |
| `v70_matched_rows` | V7.0 的 4473 行索引（用于 canonical→V7.0 下标映射） |
| `v70_ablation_pred` | V7.0 A0/A1 预测（ranker R0/R1） |
| `v70_loss_pred` | V7.0 MLP L0 预测（ranker R2） |
| `v621_donly_pred` | 历史 D-only M2（ranker R3） |
| `v63_target_table` | target 级机制表（S/G 特征源） |
| `ddg_split` | family 定义（**重新加载，不重新生成**） |
| `mutation_sites` | mutation→残基映射表 |
| `v22_features` | 结构与实验描述符 |
| `esm2_cache` | pooled ESM-2（**仅链级 640 维**） |
| `v5_common` | 特征装配与超参源 |

## 关键可用性结论（决定 PHASE 3/4 的特征选择）
- `esm2_cache` **只有 pooled 链级嵌入**，**无 residue-level 嵌入**
  → §PHASE 2 的 `z` 只能用**现有 ranker 的预测分数**做靶点内标准化（可部署）
- 结构描述符（58 个数值列）与 set 描述符（25 列）齐备
- `source` 全部为 SKEMPI v2.0（单一来源）→ source 分层不可做

## 环境
Python 3.11.9 / i7-2630QM 4 核 / 7.55 GB RAM / **无 GPU** / torch 2.5.1+cpu
