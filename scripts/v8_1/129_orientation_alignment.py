# -*- coding: utf-8 -*-
"""scripts/v8_1/129_orientation_alignment.py -- 方向对齐后重算 external panel。

方向判定（scripts/v8_1/128_orientation_audit.py，基于 within-complex 相关）：
    Flex-ddG / DDGPred / FoldX / proteinMPNN : multiplier = +1
    ESM-IF                                    : multiplier = -1（越大越强结合）

sign-only 方向对齐不涉及校准、不使用 test labels，只是把 predictor 方向统一到 y 方向。
在对齐后重算：pooled / macro / C_B / C_W / 排名 / rank agreement。
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
M = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]
MULT = {"Flex-ddG": +1, "DDGPred": +1, "ddG_FoldX": +1, "ESM-IF": -1, "proteinMPNN": +1}

PROF = R / "v8_1/external_clean_profiles.csv"
d = pd.read_csv(PROF).set_index("model")


def aligned_row(m, cols):
    """方向对齐后的读数：multiplier=-1 时 r/C_B/C_W 全部取反，macro 取反，RMSE 保持（若有意义）。"""
    o = {}
    for c in cols:
        v = d.loc[m, c]
        o[c] = (-v if MULT[m] < 0 and pd.notna(v) else v)
    return o


# ---------- 口径 A：raw（论文现用）与 B：方向对齐 ----------
rowsA, rowsB = [], []
for m in M:
    a = dict(model=m, n_rows=int(d.loc[m, "n_rows"]),
             n_ge5=int(d.loc[m, "n_targets_ge5"]))
    for c in ("overall_pearson", "macro_spearman", "C_B_row", "C_W", "rmse"):
        a[c] = float(d.loc[m, c])
    a["orientation_multiplier"] = MULT[m]
    rowsA.append(a)

    b = dict(a)
    for c in ("overall_pearson", "macro_spearman", "C_B_row", "C_W"):
        b[c] = MULT[m] * a[c]
    rowsB.append(b)

A = pd.DataFrame(rowsA).sort_values("overall_pearson", ascending=False)
B = pd.DataFrame(rowsB).sort_values("overall_pearson", ascending=False)


def rank_stats(df, tag):
    pr = df.set_index("model").overall_pearson.rank(ascending=False)
    mr = df.set_index("model").macro_spearman.rank(ascending=False)
    P = list(pr.sort_values().index); Mk = list(mr.sort_values().index)
    rho, _ = spearmanr(df.overall_pearson, df.macro_spearman)
    tau, _ = kendalltau(df.overall_pearson, df.macro_spearman)
    vals = df.overall_pearson.to_numpy(); tgt = df.macro_spearman.to_numpy()
    cnt = sum(1 for pm in itertools.permutations(range(len(df)))
              if abs(spearmanr(vals[list(pm)], tgt)[0]) >= abs(rho) - 1e-12)
    pp = cnt / 120
    swaps = [i + 1 for i in range(len(P)) if P[i] != Mk[i]]
    print(f"\n--- {tag} ---")
    print(f"  Pooled 排名: {P}")
    print(f"  Macro  排名: {Mk}")
    print(f"  位置不一致位数: {swaps}   DDGPred 双轴第一: {P[0]=='DDGPred' and Mk[0]=='DDGPred'}")
    print(f"  Spearman rho = {rho:+.4f}   Kendall tau = {tau:+.4f}   精确置换 p = {pp:.4f}")
    return dict(ranking_pooled=P, ranking_macro=Mk, mismatch_positions=swaps,
                rho=float(rho), tau=float(tau), perm_p=float(pp))


print("=" * 100)
print("口径 A：raw（论文当前使用）")
print("=" * 100)
print(A[["model", "orientation_multiplier", "n_rows", "n_ge5", "overall_pearson",
         "macro_spearman", "C_B_row", "C_W", "rmse"]].to_string(index=False))
ra = rank_stats(A, "raw")

print()
print("=" * 100)
print("口径 B：sign-only 方向对齐（ESM-IF × −1）")
print("=" * 100)
print(B[["model", "orientation_multiplier", "n_rows", "n_ge5", "overall_pearson",
         "macro_spearman", "C_B_row", "C_W", "rmse"]].to_string(index=False))
rb = rank_stats(B, "oriented")

print()
print("=" * 100)
print("结论对比")
print("=" * 100)
print(f"  rho : {ra['rho']:+.4f}  ->  {rb['rho']:+.4f}")
print(f"  tau : {ra['tau']:+.4f}  ->  {rb['tau']:+.4f}")
print(f"  p   : {ra['perm_p']:.4f}  ->  {rb['perm_p']:.4f}")
print(f"  DDGPred 双轴第一: {ra['ranking_pooled'][0]=='DDGPred' and ra['ranking_macro'][0]=='DDGPred'}"
      f"  ->  {rb['ranking_pooled'][0]=='DDGPred' and rb['ranking_macro'][0]=='DDGPred'}")

json.dump(dict(multipliers=MULT, raw=ra, oriented=rb,
               table_raw=A.to_dict("records"), table_oriented=B.to_dict("records")),
          open(R / "v8_1/EXTERNAL_ORIENTATION.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
A.to_csv(R / "v8_1/EXTERNAL_TABLE_RAW.csv", index=False, encoding="utf-8")
B.to_csv(R / "v8_1/EXTERNAL_TABLE_ORIENTED.csv", index=False, encoding="utf-8")
print("\nwrote v8_1/EXTERNAL_ORIENTATION.json + 两张表")
