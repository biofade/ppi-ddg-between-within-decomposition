# -*- coding: utf-8 -*-
"""scripts/v8_1/115_external_intersection_final.py -- 5 方法共同交集敏感性分析（最终）。

严格复刻 scripts/v8_1/05_external_clean.py 的映射与 EXACT_UNIQUE 判定，
然后在「五个方法都有非缺失预测」的共同交集上重算，与可变覆盖口径比较。
"""
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
A = R / "v8_1"
MIN_N = 5
MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]

lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
ms = pd.read_csv(R / "data/v5_0/mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y_ALL = np.asarray(z["y"], float)

C = lab.merge(ms.sort_values("pair_row").groupby("pair_row", as_index=False).agg(
    pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
    n_sites=("wt_aa", "size")), on="pair_row", how="left")
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                     .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_mut = ms.groupby("pair_row").mut_aa.first().reindex(
    C.pair_row.to_numpy()).to_numpy()
C["key_norm"] = (C.pdb4.astype(str) + "|"
                 + pd.Series(_pos).astype("float").astype("Int64").astype(str).to_numpy()
                 + "|" + pd.Series(_mut).astype(str).str.upper().to_numpy())
kc = C.key_norm.value_counts()

E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)


def parse_ext(txt):
    t = str(txt).strip().upper().replace(" ", "")
    m = re.match(r"^([A-Z]?)([A-Z])(\d+)([A-Z])$", t)
    return (None, None) if not m else (int(m.group(3)), m.group(4))


E[["pos", "mutaa"]] = pd.DataFrame([parse_ext(x) for x in E["Mutation(s)_cleaned"]],
                                   index=E.index)
E["pdb4"] = E.complex.astype(str).str.upper().str.strip().str[:4]
E["key_norm"] = (E.pdb4 + "|" + E.pos.astype("float").astype("Int64").astype(str)
                 + "|" + E.mutaa)
for c in MODELS + ["ddG"]:
    if c in E.columns:
        E[c] = pd.to_numeric(E[c], errors="coerce")
EXT = E.groupby("key_norm", as_index=False).agg({c: "mean" for c in MODELS + ["ddG"]})
C2 = C.merge(EXT, on="key_norm", how="left", suffixes=("", "_ext"))
C2["y"] = Y_ALL
C2["status"] = np.where(~C2.key_norm.isin(set(EXT.key_norm)), "UNMATCHED",
                        np.where(C2.key_norm.map(kc) > 1, "AMBIGUOUS_KEY",
                                 "EXACT_UNIQUE"))
EU = C2[C2.status == "EXACT_UNIQUE"].copy()
print(f"EXACT_UNIQUE 行数 = {len(EU)}   targets = {EU.target_id_canonical.nunique()}")


def metrics(sub, pred_col, y_col="y", tgt_col="target_id_canonical"):
    y = pd.to_numeric(sub[y_col], errors="coerce").to_numpy(float)
    p = pd.to_numeric(sub[pred_col], errors="coerce").to_numpy(float)
    t = sub[tgt_col].astype(str).to_numpy()
    m = np.isfinite(y) & np.isfinite(p)
    y, p, t = y[m], p[m], t[m]
    if len(y) < 10:
        return None
    uq, inv = np.unique(t, return_inverse=True)
    n_c = np.bincount(inv).astype(float)
    N = len(y)
    my, mp = y.mean(), p.mean()
    sd = y.std() * p.std()
    By = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bp = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_B = np.sum((n_c / N) * (By - my) * (Bp - mp)) / sd
    C_W = (np.sum((n_c / N) * np.array(
        [(y[inv == k] - By[k]).dot(p[inv == k] - Bp[k]) / n_c[k]
         for k in range(len(uq))]))) / sd
    sp = [spearmanr(y[inv == k], p[inv == k]).correlation
          for k in range(len(uq)) if n_c[k] >= MIN_N]
    sp = [x for x in sp if np.isfinite(x)]
    return dict(n_rows=int(N), n_targets=int(len(uq)),
                n_ge5=int((n_c >= MIN_N).sum()),
                pearson=float(C_B + C_W), C_B=float(C_B), C_W=float(C_W),
                macro=float(np.mean(sp)) if sp else np.nan,
                rmse=float(np.sqrt(np.mean((y - p) ** 2))))


# ---- 口径 A：可变覆盖（论文现用）----
rowsA = {}
for c in MODELS:
    sub = EU[EU[c].notna()]
    r = metrics(sub, c)
    if r:
        r["model"] = c; rowsA[c] = r

# ---- 口径 B：五方法共同交集 ----
mask = np.ones(len(EU), bool)
for c in MODELS:
    mask &= EU[c].notna().to_numpy()
EUc = EU[mask].copy()
print(f"五方法共同交集行数 = {len(EUc)}   "
      f"targets = {EUc.target_id_canonical.nunique()}   "
      f"targets(n>=5) = {int((EUc.target_id_canonical.value_counts() >= MIN_N).sum())}")
rowsB = {}
for c in MODELS:
    r = metrics(EUc, c)
    if r:
        r["model"] = c; rowsB[c] = r

A_ = pd.DataFrame(rowsA).T.reset_index(drop=True)
B_ = pd.DataFrame(rowsB).T.reset_index(drop=True)
A_ = A_.sort_values("pearson", ascending=False)
B_ = B_.sort_values("pearson", ascending=False)

print()
print("=" * 96)
print("口径 A：可变覆盖（每方法各自 uniquely-mapped subset）")
print("=" * 96)
print(A_[["model", "n_rows", "n_targets", "n_ge5", "pearson", "macro", "C_B", "C_W", "rmse"]]
      .to_string(index=False))
print()
print("=" * 96)
print("口径 B：五方法共同交集")
print("=" * 96)
print(B_[["model", "n_rows", "n_targets", "n_ge5", "pearson", "macro", "C_B", "C_W", "rmse"]]
      .to_string(index=False))


def rank_report(df, tag):
    pr = df.set_index("model").pearson.rank(ascending=False)
    mr = df.set_index("model").macro.rank(ascending=False)
    P = list(pr.sort_values().index)
    M = list(mr.sort_values().index)
    rho, _ = spearmanr(df.pearson, df.macro)
    tau, _ = kendalltau(df.pearson, df.macro)
    cnt = sum(1 for pm in itertools.permutations(range(len(df)))
              if abs(spearmanr(df.pearson.to_numpy()[list(pm)],
                               df.macro.to_numpy())[0]) >= abs(rho) - 1e-12)
    pp = cnt / 120
    swaps = [i + 1 for i in range(len(P)) if P[i] != M[i]]
    print()
    print(f"--- {tag} ---")
    print(f"  Pooled 排名: {P}")
    print(f"  Macro  排名: {M}")
    print(f"  位置交换: {swaps}  DDGPred 双轴第一: {P[0]=='DDGPred' and M[0]=='DDGPred'}")
    print(f"  Spearman rho = {rho:+.4f}   Kendall tau = {tau:+.4f}   精确置换 p = {pp:.4f}")
    return dict(ranking_pooled=P, ranking_macro=M, swap_positions=swaps,
                rho=float(rho), tau=float(tau), perm_p=float(pp))


ra = rank_report(A_, "口径 A 可变覆盖")
rb = rank_report(B_, "口径 B 共同交集")

out = dict(
    coverage_variable=rowsA, coverage_intersection=rowsB,
    rank_variable=ra, rank_intersection=rb,
    n_intersection_rows=int(len(EUc)),
    n_intersection_targets=int(EUc.target_id_canonical.nunique()),
    n_intersection_ge5=int((EUc.target_id_canonical.value_counts() >= MIN_N).sum()),
    n_exact_unique=int(len(EU)),
)
json.dump(out, open(A / "EXTERNAL_COMMON_INTERSECTION.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
A_.to_csv(A / "EXTERNAL_COVERAGE_VARIABLE.csv", index=False, encoding="utf-8")
B_.to_csv(A / "EXTERNAL_COVERAGE_INTERSECTION.csv", index=False, encoding="utf-8")
print()
print("wrote v8_1/EXTERNAL_COMMON_INTERSECTION.json + 两张 csv")
