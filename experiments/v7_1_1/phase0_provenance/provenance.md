
# V7.1.1 PHASE 0 — Provenance

**时间**：2026-10-14
**约束**：不修改 V1–V7.1 任何文件；新命名空间 `v7_1_1`；CPU only；
**未下载任何模型**；**未训练任何新模型**；未开始 residue-level embedding；未开始 V7.2；
模型选择不使用测试标签；旧报告文本错误**不修改**，只在 V7.1.1 记录 correction。

## 冻结输入（17 个，SHA256 见 input_manifest.csv）

| 组 | 名称 |
|---|---|
| V7.1 输出 | `v71_u_local_rows`、`v71_u_full_rows`、`v71_ranker_metrics`、`v71_z_predictions`、`v71_location_pred`、`v71_location_models`、`v71_scale_models`、`v71_lsr_models`、`v71_lsr_predictions`、`v71_oracle`、`v71_contrasts`、`v71_stratification`、`v71_final_summary` |
| V7.0 输出 | `v70_matched_rows`、`v70_block_L` |
| 上游权威 | `canonical_matched_set`、`ddg_split_family` |

## 环境
Python 3.11.9 / i7-2630QM 4 核 / 7.55 GB RAM / **无 GPU**

## 本轮不做什么
- 不训练任何模型（只做对账、资格审计、协议冻结）
- 不实现 `R_RESIDUE_*` 家族
- 不修改任何 V7.0/V7.1 文件（含其报告文本错误）
