# Model Capability Card

**Model**: `L1_random_M1_s42`
**Family**: V6_2_split4
**Provenance status**: YELLOW
**Evidence tier**: TIER_1_INTERNAL
**Prediction file**: `experiments/v6_2/phase2_split4/split4_predictions.npz`
**Join mode**: mas_align

## Universe

Universe: `U_MAS`
Rows: 7523 | Targets: 230 | n>=5: 172 | min_n: 5

## Calibration

Overall Pearson: 0.6990
Overall Spearman: 0.5742
RMSE: 1.7360
MAE: 1.2702
Target-mean Pearson: 0.8161
Target-mean Spearman: 0.6582
Attenuation slope: 0.4684
C_B row: 0.4275
C_B equal: 0.5588

## Ranking

Macro Spearman: 0.0876
Macro Kendall tau-b: 0.0599
Pairwise δ=0: 0.5298 | δ=0.2: 0.5359 | δ=0.5: 0.5427
C_W: 0.2715

## Selection

Hit@1: 0.1337 | Hit@5: 0.4593 | Hit@10: 0.6802
Precision@10: 0.6302 | Recall@10: 0.6302
Efficiency@5: 0.7131 | Efficiency@10: 0.8137
NDCG@5: 0.6889 | NDCG@10: 0.7513
EF@1%: 3.5630 | EF@5%: 1.9285

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
