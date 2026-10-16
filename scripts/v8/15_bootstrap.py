# -*- coding: utf-8 -*-
"""scripts/v8/15_bootstrap.py -- V8 PHASE 15：target-cluster bootstrap（§55/§56）。

优先检验 H1 的核心量：**模型排名之间的秩一致性是否显著偏离 0**。
重采样单位 = **target**（不是行），seed=42，N=5000，要求有效迭代 = 5000。

同时做**模型层面**的配对对比（§42 regret 的 bootstrap）。
CPU-first：预计算 target 索引，避免在每次迭代里 groupby。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
BS = OUT / "bootstrap"
BS.mkdir(parents=True, exist_ok=True)
N_BOOT = 5000
SEED = 42
MIN_N = 5

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
CAN = lab[["pair_row", "target_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id"})
TGT = CAN.target_id.to_numpy(str)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
           allow_pickle=True)
Y = np.asarray(z["y"], float)

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))]
uc = uc.head(40)          # 控制成本：取注册顺序前 40 个（含全部主要族）
print(f"bootstrap 模型数 = {len(uc)} | rows = {len(Y)} | targets = {len(set(TGT))}")

# ---- 载入预测矩阵 P (n_models × n_rows) ----
_cache: dict = {}
P_list, names, fams = [], [], []
for r in uc.itertuples():
    rel = r.prediction_file
    if rel not in _cache:
        _cache[rel] = np.load(R / rel, allow_pickle=True)
    zz = _cache[rel]
    if r.model_name not in zz.files:
        continue
    a = np.asarray(zz[r.model_name])
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        continue
    if a.ndim != 1 or len(a) != len(Y):
        continue
    P_list.append(a); names.append(r.model_name); fams.append(r.model_family)
P = np.vstack(P_list)
M_, N_ = P.shape
print(f"预测矩阵 P = {P.shape}")

# ---- 预计算 target 索引 ----
uq = np.unique(TGT)
tidx = [np.where(TGT == t)[0] for t in uq]
tidx = [i for i in tidx if i.size >= MIN_N]
print(f"n>={MIN_N} 的 target 数 = {len(tidx)}")

# ---- 单次评价：给定行子集，返回模型层面的指标矩阵 ----
def metrics_for(rows_sel):
    y = Y[rows_sel]
    t = TGT[rows_sel]
    uq2 = np.unique(t)
    idx2 = {u: np.where(t == u)[0] for u in uq2}
    elig = [u for u in uq2 if idx2[u].size >= MIN_N]
    out = {k: np.full(M_, np.nan) for k in
           ("overall_pearson", "macro_spearman", "efficiency@10", "ndcg@10")}
    for mi in range(M_):
        p = P[mi, rows_sel]
        m = np.isfinite(y) & np.isfinite(p)
        if m.sum() < 50:
            continue
        yy, pp, tt = y[m], p[m], t[m]
        out["overall_pearson"][mi] = M._pearson(yy, pp)
        sp, ef, nd = [], [], []
        uq3 = np.unique(tt)
        ix3 = {u: np.where(tt == u)[0] for u in uq3}
        for u in uq3:
            i = ix3[u]
            if i.size < MIN_N:
                continue
            sp.append(M._spearman(yy[i], pp[i]))
            n = i.size; ar = np.arange(n)
            ty = yy[i]; tp = pp[i]
            true_order = np.lexsort((ar, -(-ty)))
            rank = np.empty(n, int); rank[true_order] = ar
            gain = 1.0 + (n - rank)
            ideal = np.lexsort((ar, -gain))
            pred = np.lexsort((ar, tp))
            kk = max(1, min(10, n))
            disc = np.log2(np.arange(2, kk + 2))
            dcg = (gain[pred[:kk]] / disc).sum()
            idcg = (gain[ideal[:kk]] / disc).sum()
            nd.append(float(dcg / idcg) if idcg > 0 else np.nan)
            gsel = gain[pred[:kk]].sum(); gid = gain[ideal[:kk]].sum()
            ef.append(float(gsel / gid) if gid > 0 else np.nan)
        out["macro_spearman"][mi] = float(np.nanmean(sp)) if sp else np.nan
        out["efficiency@10"][mi] = float(np.nanmean(ef)) if ef else np.nan
        out["ndcg@10"][mi] = float(np.nanmean(nd)) if nd else np.nan
    return out


def rank_agreement(mt, a, b):
    v = np.isfinite(mt[a]) & np.isfinite(mt[b])
    if v.sum() < 8:
        return np.nan, np.nan
    return (float(spearmanr(mt[a][v], mt[b][v]).correlation),
            float(kendalltau(mt[a][v], mt[b][v]).correlation))


PAIRS = [("overall_pearson", "macro_spearman"),
         ("overall_pearson", "efficiency@10"),
         ("overall_pearson", "ndcg@10"),
         ("macro_spearman", "efficiency@10")]

allrows = np.arange(N_)
base = metrics_for(allrows)
point = {}
for a, b in PAIRS:
    sp, kt = rank_agreement(base, a, b)
    point[(a, b)] = (sp, kt)
    print(f"  point: {a} vs {b}: rho={sp:+.4f} tau={kt:+.4f}")

print(f"\n=== target-cluster bootstrap N={N_BOOT} ===")
rs = np.random.RandomState(SEED)
acc = {p: {"sp": [], "kt": []} for p in PAIRS}
t0 = time.time()
valid = 0
for it in range(N_BOOT):
    pick = rs.randint(0, len(tidx), len(tidx))
    rows_sel = np.concatenate([tidx[j] for j in pick])
    if np.unique(TGT[rows_sel]).size < 10:
        continue
    mt = metrics_for(rows_sel)
    ok = True
    tmp = {}
    for a, b in PAIRS:
        sp, kt = rank_agreement(mt, a, b)
        if not (np.isfinite(sp) and np.isfinite(kt)):
            ok = False; break
        tmp[(a, b)] = (sp, kt)
    if not ok:
        continue
    valid += 1
    for p in PAIRS:
        acc[p]["sp"].append(tmp[p][0]); acc[p]["kt"].append(tmp[p][1])
    if (it + 1) % 500 == 0:
        print(f"  {it+1}/{N_BOOT} valid={valid} "
              f"({(time.time()-t0)/60:.1f} min)", flush=True)

print(f"\n有效迭代 = {valid} / {N_BOOT}  ({(time.time()-t0)/60:.1f} min)")
rows = []
for a, b in PAIRS:
    sp = np.asarray(acc[(a, b)]["sp"]); kt = np.asarray(acc[(a, b)]["kt"])
    lo_s, hi_s = np.percentile(sp, [2.5, 97.5])
    lo_k, hi_k = np.percentile(kt, [2.5, 97.5])
    exc_s = bool(lo_s > 0 or hi_s < 0)
    rows.append({"pair": f"{a} vs {b}", "metric": "spearman",
                 "point": point[(a, b)][0], "ci_low": lo_s, "ci_high": hi_s,
                 "ci_excludes_zero": exc_s, "n_valid": valid})
    rows.append({"pair": f"{a} vs {b}", "metric": "kendall",
                 "point": point[(a, b)][1], "ci_low": lo_k, "ci_high": hi_k,
                 "ci_excludes_zero": bool(lo_k > 0 or hi_k < 0), "n_valid": valid})
    print(f"  {a} vs {b}: rho={point[(a,b)][0]:+.4f} [{lo_s:+.4f},{hi_s:+.4f}] "
          f"{'✅排除0' if exc_s else '❌含0'} | tau={point[(a,b)][1]:+.4f} "
          f"[{lo_k:+.4f},{hi_k:+.4f}]")
B = pd.DataFrame(rows)
B.to_csv(BS / "rank_agreement_bootstrap.csv", index=False, encoding="utf-8")
np.savez(BS / "bootstrap_samples.npz",
         **{f"{a}|{b}|{m}": np.asarray(acc[(a, b)][m])
            for a, b in PAIRS for m in ("sp", "kt")})
print(f"\n写出 rank_agreement_bootstrap.csv / bootstrap_samples.npz")
