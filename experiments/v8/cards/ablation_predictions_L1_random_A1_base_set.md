# Model Capability Card

**Model**: `L1_random_A1_base_set`
**Family**: V6_2_split4
**Provenance status**: GREEN
**Evidence tier**: TIER_1_INTERNAL
**Prediction file**: `experiments/v7_0/phase2_ablation/ablation_predictions.npz`
**Join mode**: canonical_align

## Universe

Universe: `U_CANONICAL`
Rows: 4473 | Targets: 214 | n>=5: 130 | min_n: 5

## Calibration

Overall Pearson: 0.7729
Overall Spearman: 0.6898
RMSE: 1.5464
MAE: 1.1319
Target-mean Pearson: 0.7834
Target-mean Spearman: 0.7221
Attenuation slope: 0.4448
C_B row: 0.5167
C_B equal: 0.5135

## Ranking

Macro Spearman: 0.0747
Macro Kendall tau-b: 0.0472
Pairwise δ=0: 0.5232 | δ=0.2: 0.5228 | δ=0.5: 0.5358
C_W: 0.2563

## Selection

Hit@1: 0.1231 | Hit@5: 0.4692 | Hit@10: 0.7231
Precision@10: 0.7069 | Recall@10: 0.7069
Efficiency@5: 0.7471 | Efficiency@10: 0.8517
NDCG@5: 0.7131 | NDCG@10: 0.7783
EF@1%: 3.8014 | EF@5%: 2.3249

## Split sensitivity

本项目对该模型的**划分敏感性**记录于 `SPLIT_INFLATION_REPORT.md`（仅 V6.2 族有跨划分预测）。

## Target-weighting sensitivity

见 `WEIGHTING_SENSITIVITY_REPORT.md`；C_B 对加权敏感，C_W 稳定。

## Main strength / weakness

（依三轴读数人工填写；本卡只提供数据）

## Appropriate use case

Quantitative ΔΔG / Library ranking / Wet-lab shortlist / Other —— **需按三轴分别判断**。

## Caveats

* 本卡为**自动生成**，未做人工解读。
* `C_B`/`C_W` 是**协方差信号贡献**，不是独立标准化能力（§29）。
* 选择效用轴（Eff@10 等）在本 universe 上**区分度有限**（见 `DECISION_ALIGNMENT_REPORT.md`）。
