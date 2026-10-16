# V7.0 PHASE 0 + 1 — 冻结 + 特征块构造

- 时间: 2026-10-14 14:02:40

## PHASE 0 — 输入冻结

- `matched_analysis_set`: 754,492 B  `8d209eea08cc6058…`
- `split4_predictions`: 1,696,726 B  `ba711ea1b349836a…`
- `split4_results`: 4,253 B  `52d812e77c1f0169…`
- `d_only_predictions`: 484,102 B  `98b64c2b02140f23…`
- `v63_target_table`: 285,776 B  `b2b241cd2afc0dac…`
- `v63_feature_dict`: 7,653 B  `7495a0f977daaf40…`
- `ddg_split`: 82,187 B  `92d4e8d78f0778ff…`
- `mutation_sites`: 706,966 B  `ca23d3b42eadf078…`
- `v22_features`: 11,852,893 B  `134a406175c97554…`
- `esm2_cache`: 52,481,097 B  `2b6fdcc331dd5030…`
- `v5_common`: 27,748 B  `91bee0814f66737b…`
- `standardized_preds`: 1,270,276 B  `bebfc4fad0137d0f…`

### Baseline B0（= 历史 M0，**复用不重写**）

| 项 | 值 |
|---|---|
| 定义 | pooled baseline：`[phys_diff(60), Z_diff(128)]` + fwd/rev 增广 |
| 实现 | `scripts/v6_2/02_split4_ddg.py`（V6.2 §2.0 已审计） |
| 预测来源 | `experiments/v6_2/phase2_split4/split4_predictions.npz` 的 `L4_family_M0_s{42,43,44}` |
| LightGBM | `LGB_PARAMS` + `n_estimators=300` + `random_state=seed` |
| 特征维度 | 188（60 物理 + 128 pooled ESM-2 PCA） |

