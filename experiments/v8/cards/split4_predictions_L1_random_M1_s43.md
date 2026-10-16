# Model Capability Card

**Model**: `L1_random_M1_s43`
**Family**: V6_2_split4
**Provenance status**: YELLOW
**Evidence tier**: TIER_1_INTERNAL
**Prediction file**: `experiments/v6_2/phase2_split4/split4_predictions.npz`
**Join mode**: mas_align

## Universe

Universe: `U_MAS`
Rows: 7523 | Targets: 230 | n>=5: 172 | min_n: 5

## Calibration

Overall Pearson: 0.6987
Overall Spearman: 0.5777
RMSE: 1.7361
MAE: 1.2681
Target-mean Pearson: 0.8120
Target-mean Spearman: 0.6613
Attenuation slope: 0.4550
C_B row: 0.4293
C_B equal: 0.5401

## Ranking

Macro Spearman: 0.0765
Macro Kendall tau-b: 0.0458
Pairwise δ=0: 0.5231 | δ=0.2: 0.5262 | δ=0.5: 0.5296
C_W: 0.2694

## Selection

Hit@1: 0.1105 | Hit@5: 0.4360 | Hit@10: 0.6860
Precision@10: 0.6320 | Recall@10: 0.6320
Efficiency@5: 0.7188 | Efficiency@10: 0.8167
NDCG@5: 0.6907 | NDCG@10: 0.7517
EF@1%: 2.7316 | EF@5%: 1.9871

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
