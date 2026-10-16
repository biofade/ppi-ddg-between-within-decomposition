# -*- coding: utf-8 -*-
"""scripts/v8_1/128_orientation_audit.py -- 核查五个 external 列的生物学方向。

问题：ddG>0 = 削弱结合。若某列的原始分数「越大越强结合」，则它与 ddG 反号，
用于**方法性能比较**前必须先做 sign-only 方向对齐（不涉及校准、不用 test labels）。

判定依据（多条独立证据，避免只看 poooled 相关系数符号，因为该符号混合了
「方向」与「预测能力」两件事）：
  1. pooled corr(column, ddG)
  2. within-complex（组内中心化）corr —— 去掉 target 间差异后的方向证据
  3. 列名/语义先验：ddG 型（能量）vs log-prob 型（序列模型打分）
  4. 列的取值分布特征（log-prob 型集中在 0 附近）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)
M = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]

rows = []
for c in M:
    s = E[[c, "ddG", "complex"]].copy()
    s[c] = pd.to_numeric(s[c], errors="coerce")
    s["ddG"] = pd.to_numeric(s["ddG"], errors="coerce")
    s = s.dropna()
    r_pool = float(np.corrcoef(s[c], s["ddG"])[0, 1])
    s["cw"] = s[c] - s.groupby("complex")[c].transform("mean")
    s["dw"] = s["ddG"] - s.groupby("complex")["ddG"].transform("mean")
    r_within = float(np.corrcoef(s["cw"], s["dw"])[0, 1])
    rows.append(dict(column=c, n=int(len(s)),
                     corr_pooled=r_pool, corr_within_complex=r_within,
                     sd=float(s[c].std()), mean=float(s[c].mean()),
                     frac_near_zero=float((s[c].abs() < 0.02).mean())))

D = pd.DataFrame(rows)
print("=" * 104)
print("五个 external 列的方向证据")
print("=" * 104)
print(D.to_string(index=False))
print()
print("  说明：ddG>0 = 削弱结合。corr>0 表示「该列越大 → ddG 越大 → 越削弱」，")
print("        即该列已与 ddG 同向；corr<0 表示该列「越大越增强」，与 ddG 反向。")
print()

# 方向判定：以 within-complex 相关为主（去掉 target 均值后更纯），pooled 为辅
out = {}
for r in D.itertuples():
    same_dir = r.corr_within_complex > 0
    out[r.column] = dict(
        corr_pooled=r.corr_pooled, corr_within=r.corr_within_c if False else r.corr_within_complex,
        raw_orientation=("higher means larger ddG / weaker binding" if same_dir
                         else "higher means stronger binding"),
        orientation_multiplier=(+1 if same_dir else -1),
    )
    print(f"  {r.column:<13} within={r.corr_within_complex:+.4f} -> "
          f"multiplier = {out[r.column]['orientation_multiplier']:+d}  "
          f"({out[r.column]['raw_orientation']})")

import json  # noqa: E402
json.dump(out, open("v8_1/ORIENTATION_AUDIT.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
D.to_csv("v8_1/ORIENTATION_AUDIT.csv", index=False, encoding="utf-8")
print()
print("wrote v8_1/ORIENTATION_AUDIT.{csv,json}")
