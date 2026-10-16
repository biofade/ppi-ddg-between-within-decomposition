# Model Capability Card

**Model**: `R1_Base_F0_s42`
**Family**: V7.2R1_representation
**Provenance status**: GREEN
**Evidence tier**: TIER_1_INTERNAL
**Prediction file**: `experiments/v7_2/phase5_nested_cv/ranker_predictions.npz`
**Join mode**: canonical_rowid

## Universe

Universe: `U_CANONICAL`
Rows: 1807 | Targets: 13 | n>=5: 12 | min_n: 5

## Calibration

Overall Pearson: 0.3241
Overall Spearman: 0.3072
RMSE: 2.9141
MAE: 2.1718
Target-mean Pearson: 0.7671
Target-mean Spearman: 0.3901
Attenuation slope: 0.3464
C_B row: -0.0027
C_B equal: 0.7888

## Ranking

Macro Spearman: 0.3388
Macro Kendall tau-b: 0.2740
Pairwise δ=0: 0.6374 | δ=0.2: 0.6414 | δ=0.5: 0.6502
C_W: 0.3268

## Selection

Hit@1: 0.0833 | Hit@5: 0.2500 | Hit@10: 0.2500
Precision@10: 0.3250 | Recall@10: 0.3250
Efficiency@5: 0.7000 | Efficiency@10: 0.7109
NDCG@5: 0.6818 | NDCG@10: 0.6911
EF@1%: 0.6667 | EF@5%: 2.8945

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
