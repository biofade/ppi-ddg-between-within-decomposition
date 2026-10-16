# HERMES 语义审计（§13.1–13.2）

* 官方来源：`https://github.com/StatPhysBio/hermes_experiments/skempi/results/hermes_py_000_ft_cdna117k_ddg_st/`
* 官方量：`log_proba_mt__minus__log_proba_wt`
* 语义：dimensionless mutation-effect score = log P(mutant) - log P(wild-type) under a structure-conditioned inverse-folding model
* 训练任务：**cDNA117k stability fine-tune (per official variant naming)**
* SKEMPI binding 标签暴露：**NOT for the cDNA117k variant used here; the model name must NOT be read as SKEMPI-supervised binding regression**
* 预测语义：score only; NOT calibrated kcal/mol ddG
* 符号约定：hermes_proxy = -(log_proba_mt - log_proba_wt); aligned so that V8 ddG>0 = weaker binding; direction verified empirically (r=+0.1688 > 0)
* **RMSE = NOT_APPLICABLE；MAE = NOT_APPLICABLE**
* 突变类型：single-point only in the primary matched set (multi-point rows do not map)

## 映射 QC（§13.2）

```
canonical rows before merge = 4473
canonical rows after  merge = 4473   (no expansion)
duplicate canonical row_id  = 0
exact match                 = 2865
ambiguous                   = 0
unmatched                   = 1608
matched single-point        = 2865
```

> **禁用**：用原始 log-odds 声称 calibrated kcal/mol 精度（§13.1）。
