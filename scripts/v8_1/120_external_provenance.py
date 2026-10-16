# -*- coding: utf-8 -*-
"""scripts/v8_1/120_external_provenance.py -- external 预测列的 provenance / 单位 / 语义。

来源核实（真实文件与代码）：
  * 数据表：F:/dsh/downloads/cathddg/skempi_v2.csv（CATH-ddG 发布的 SKEMPI 注释表）
  * 载入：scripts/v8_1/05_external_clean.py —— 直接取列名 'Flex-ddG','DDGPred',
    'ddG_FoldX','ESM-IF','proteinMPNN'，**未做任何变换**（无符号翻转、无 affine 校准）
  * 项目自身的语义审计 scripts/v8_1/04_external_audits.py 对同类打分明确写：
      "prediction_semantics": "score only; NOT calibrated kcal/mol ddG"
      "RMSE": "NOT_APPLICABLE"

本脚本量化每个外部列的原始量纲，判定 RMSE 是否可用，并输出 Supplementary 用表。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)
MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]

# 每个方法的引用（与 references.bib 对应）
CITE = {"Flex-ddG": "flexddg", "DDGPred": "ddgpred", "ddG_FoldX": "foldx5",
        "ESM-IF": "esmif", "proteinMPNN": "proteinmpnn"}
# 原文方法类型
KIND = {"Flex-ddG": "Rosetta ensemble ΔΔG predictor",
        "DDGPred": "deep-learning ΔΔG predictor",
        "ddG_FoldX": "empirical force-field ΔΔG predictor",
        "ESM-IF": "inverse-folding sequence model",
        "proteinMPNN": "inverse-folding sequence model"}

rows = []
for m in MODELS:
    s = E[[m, "ddG"]].apply(pd.to_numeric, errors="coerce").dropna()
    v = s[m]
    ddg_like = bool(v.std() > 0.5)
    rows.append(dict(
        method=m, raw_column=m, raw_n=int(len(v)),
        raw_mean=float(v.mean()), raw_sd=float(v.std()),
        raw_min=float(v.min()), raw_max=float(v.max()),
        corr_raw_with_experimental=float(np.corrcoef(v, s["ddG"])[0, 1]),
        raw_quantity=("predicted binding-affinity change (ΔΔG)" if ddg_like
                      else "dimensionless mutation-effect score"),
        units=("kcal/mol (ΔΔG-like scaling)" if ddg_like else "dimensionless (no kcal/mol scaling)"),
        transformation="none (raw column used as supplied)",
        orientation="used as supplied; higher = stronger binding",
        rmse_valid=ddg_like,
        cite=CITE[m], underlying_model=KIND[m]))
D = pd.DataFrame(rows)
D.to_csv("v8_1/EXTERNAL_PROVENANCE.csv", index=False, encoding="utf-8")

print("=" * 100)
print("外部预测列的量纲与语义")
print("=" * 100)
print(D[["method", "raw_sd", "raw_min", "raw_max", "corr_raw_with_experimental",
         "units", "rmse_valid"]].to_string(index=False))
print()
n_invalid = int((~D.rmse_valid).sum())
print(f"  RMSE 不可用的方法数 = {n_invalid} / {len(D)}")
for r in D[~D.rmse_valid].itertuples():
    print(f"    - {r.method}: raw sd = {r.raw_sd:.4f} —— 无量纲分数，RMSE(kcal/mol) 无意义")

json.dump(dict(source_table="F:/dsh/downloads/cathddg/skempi_v2.csv",
               loader="scripts/v8_1/05_external_clean.py (raw columns, no transformation)",
               semantic_precedent=("scripts/v8_1/04_external_audits.py: "
                                   "'score only; NOT calibrated kcal/mol ddG', "
                                   "RMSE = NOT_APPLICABLE"),
               columns=D.to_dict("records")),
          open("v8_1/EXTERNAL_PROVENANCE.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print()
print("wrote v8_1/EXTERNAL_PROVENANCE.{csv,json}")
