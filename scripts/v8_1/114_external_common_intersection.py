# -*- coding: utf-8 -*-
"""scripts/v8_1/114_external_common_intersection.py -- 5 方法共同交集敏感性分析。

背景：五个外部方法的预测来自同一张 SKEMPI 表（scripts/v8_1/05_external_clean.py），
只是各自缺失行不同 —— Flex-ddG 1,755 / DDGPred 1,756 / 其余 1,776。
因此「共同交集」定义明确：五个方法**都**有非缺失预测的 EXACT_UNIQUE 行。

本脚本：
  1. 在 EXACT_UNIQUE 映射上取五方法共同交集
  2. 在交集上重算 r / C_B / C_W / macro Spearman / RMSE
  3. 比较交集口径与可变覆盖口径的排名与 rank agreement
  4. 输出 Supplementary 用的小表
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

sys.path.insert(0, str(Path("scripts/v8")))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
A = R / "v8_1"

MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]

# ---------- 复刻 05_external_clean.py 的映射与 status ----------
lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
CAN = lab[["pair_row", "target_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id"})
ms = pd.read_csv(R / "data/v5_0/mutation_sites.csv", low_memory=False)
print("mutation_sites 列:", [c for c in ms.columns][:14])


def norm_key(pdb, mut):
    return (str(pdb).strip().upper() + "|" +
            str(mut).replace(" ", "").strip().upper())


# 尝试从 mutation_sites 组合 key；若结构不同则退回 05 的逻辑
E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)
print("skempi_v2 列:", [c for c in E.columns][:16])


def parse_ext(txt):
    import re
    m = re.match(r"^([A-Za-z])(\d+)([A-Za-z])$", str(txt).strip())
    return (m.group(2), m.group(3)) if m else (np.nan, np.nan)


E[["pos", "mutaa"]] = pd.DataFrame([parse_ext(x) for x in E["Mutation(s)_cleaned"]],
                                   index=E.index)
print()
print("E 行数:", len(E))
print("含 key_norm 列:", "key_norm" in E.columns)
keycols = [c for c in E.columns if "key" in c.lower() or "hash" in c.lower()]
print("key 相关列:", keycols)
print()
print(E.head(3).to_string()[:900])
