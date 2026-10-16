# -*- coding: utf-8 -*-
"""scripts/v8/15_bootstrap_fast.py -- V8 PHASE 15：**向量化** target-cluster bootstrap（§55/§56）。

前一版按 (model × target × iteration) 三重 Python 循环，实测 ~4 小时。
本版把所有模型堆成一个矩阵，用 NumPy 一次性算完一次 bootstrap：

    P : (M 模型, N 行)          Y : (N,)
    T : (N,) target 代码        G : (N,) 组号（重采样后重排）

一次迭代的成本降到 O(M·N) 的向量运算，**无 Python 逐模型循环**。

被检验的量（H1）：模型**排名**之间的秩一致性
    ρ / τ ( Overall Pearson , Macro Spearman )
    ρ / τ ( Overall Pearson , Efficiency@10 )
    ρ / τ ( Overall Pearson , NDCG@10 )
    ρ / τ ( Macro Spearman  , Efficiency@10 )
"""

from __future__ import annotations

import argparse
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
OUT = R / "experiments" / "v8"
BS = OUT / "bootstrap"
BS.mkdir(parents=True, exist_ok=True)
SEED = 42
MIN_N = 5
K = 10


def build_data(max_models=None):
    lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    TGT = lab.target_id_canonical.to_numpy(str)
    z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
               allow_pickle=True)
    Y = np.asarray(z["y"], float)
    S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
    uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))]
    if max_models:
        uc = uc.head(max_models)
    cache, P, names = {}, [], []
    for r in uc.itertuples():
        rel = r.prediction_file
        if rel not in cache:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        zz = cache[rel]
        if r.model_name not in zz.files:
            continue
        a = np.asarray(zz[r.model_name])
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if a.ndim != 1 or len(a) != len(Y):
            continue
        P.append(a); names.append(r.model_name)
    return Y, TGT, np.vstack(P), names


def metrics_vec(Y, T, P, tcode, n_t, min_n=MIN_N, k=K):
    """给定重排后的 Y/T/P 与目标编码 tcode（0..n_t-1，已按组连续排序），
    向量化返回 (M, 4) 指标矩阵：overall_pearson, macro_spearman, eff@10, ndcg@10。"""
    M_, N_ = P.shape
    cnt = np.bincount(tcode, minlength=n_t)
    # 组边界（因为已按 tcode 稳定排序，组是连续的）
    starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    elig = cnt >= min_n
    # ---- overall pearson（逐模型，向量化） ----
    Yc = Y - Y.mean()
    Pc = P - P[:, None].mean(axis=1)[:, None] if False else P - P.mean(axis=1,
                                                                     keepdims=True)
    denom = np.sqrt((Yc ** 2).sum()) * np.sqrt((Pc ** 2).sum(axis=1))
    denom[denom == 0] = np.nan
    overall = (Pc @ Yc) / denom
    # ---- 组内排名 / gain / 预测序 ----
    macro_sp = np.zeros(M_)
    eff = np.zeros(M_)
    ndcg = np.zeros(M_)
    for g in np.where(elig)[0]:
        s = starts[g]; e = s + cnt[g]; n = cnt[g]
        yg = Y[s:e]
        pg = P[:, s:e]                       # (M, n)
        # 真值名次（utility = -y 降序）；稳定 => 先按 y 升序
        true_order = np.argsort(yg, kind="stable")     # 最小 y = 最有利
        rank = np.empty(n, int); rank[true_order] = np.arange(n)
        gain = 1.0 + (n - rank)
        ideal = np.argsort(-gain, kind="stable")
        kk = max(1, min(k, n))
        disc = np.log2(np.arange(2, kk + 2))     # 长度 = kk
        gid = (gain[ideal[:kk]] / disc).sum()
        # 预测：pg 越小越有利 -> 逐模型稳定排序
        pred = np.argsort(pg, axis=1, kind="stable")[:, :kk]     # (M, kk)
        gsel_sum = gain[pred].sum(axis=1)                        # (M,)
        dcg = (gain[pred] / disc).sum(axis=1)
        if gid > 0:
            eff += gsel_sum / gid
            ndcg += dcg / gid
        # Spearman：组内对 y 与 p 各自排名后求相关（向量化）
        ry = np.argsort(np.argsort(yg, kind="stable"), kind="stable").astype(float)
        rp = np.argsort(np.argsort(pg, axis=1, kind="stable"),
                        axis=1, kind="stable").astype(float)
        ry_c = ry - ry.mean()
        rp_c = rp - rp.mean(axis=1, keepdims=True)
        d = np.sqrt((ry_c ** 2).sum()) * np.sqrt((rp_c ** 2).sum(axis=1))
        d[d == 0] = np.nan
        macro_sp += (rp_c @ ry_c) / d
    n_elig = int(elig.sum())
    if n_elig:
        macro_sp /= n_elig; eff /= n_elig; ndcg /= n_elig
    return np.vstack([overall, macro_sp, eff, ndcg]).T


def rank_agree(mt, i, j):
    v = np.isfinite(mt[:, i]) & np.isfinite(mt[:, j])
    if v.sum() < 8:
        return np.nan, np.nan
    a, b = mt[v, i], mt[v, j]
    if a.std() == 0 or b.std() == 0:
        return np.nan, np.nan
    return (float(spearmanr(a, b).correlation),
            float(kendalltau(a, b).correlation))


PAIRS = [(0, 1, "overall_pearson vs macro_spearman"),
         (0, 2, "overall_pearson vs efficiency@10"),
         (0, 3, "overall_pearson vs ndcg@10"),
         (1, 2, "macro_spearman vs efficiency@10")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=5000)
    ap.add_argument("--max-models", type=int, default=None)
    ap.add_argument("--max-minutes", type=float, default=90.0)
    a = ap.parse_args()

    Y, TGT, P, names = build_data(a.max_models)
    M_, N_ = P.shape
    print(f"模型 = {M_} | 行 = {N_} | targets = {len(set(TGT))}")
    uq, tcode = np.unique(TGT, return_inverse=True)
    n_t = len(uq)
    print(f"target 编码完成（{n_t} 个）")

    mt0 = metrics_vec(Y, TGT, P, tcode, n_t)
    point = {}
    for i, j, lab in PAIRS:
        sp, kt = rank_agree(mt0, i, j)
        point[(i, j)] = (sp, kt)
        print(f"  point {lab}: rho={sp:+.4f} tau={kt:+.4f}")

    # 每个 target 的行索引
    tidx = [np.where(tcode == u)[0] for u in range(n_t)]
    tidx = [x for x in tidx if x.size >= MIN_N]
    print(f"n>={MIN_N} target = {len(tidx)}")

    rs = np.random.RandomState(SEED)
    acc = {p[:2]: {"sp": [], "kt": []} for p in PAIRS}
    t0 = time.time(); valid = 0
    for it in range(a.n_boot):
        pick = rs.randint(0, len(tidx), len(tidx))
        rows = np.concatenate([tidx[j] for j in pick])
        # 重排为按组连续
        order = np.argsort(np.concatenate([np.full(tidx[j].size, i)
                                           for i, j in enumerate(pick)]),
                           kind="stable")
        rows = rows[order]
        newcode = np.concatenate([np.full(tidx[j].size, i)
                                  for i, j in enumerate(pick)])[order]
        Ys = Y[rows]; Ps = P[:, rows]
        mt = metrics_vec(Ys, TGT[rows], Ps, newcode, len(pick))
        tmp = {}
        ok = True
        for i, j, lab in PAIRS:
            sp, kt = rank_agree(mt, i, j)
            if not (np.isfinite(sp) and np.isfinite(kt)):
                ok = False; break
            tmp[(i, j)] = (sp, kt)
        if not ok:
            continue
        valid += 1
        for p in PAIRS:
            acc[p[:2]]["sp"].append(tmp[p[:2]][0])
            acc[p[:2]]["kt"].append(tmp[p[:2]][1])
        if (it + 1) % 250 == 0:
            el = (time.time() - t0) / 60
            eta = el / (it + 1) * (a.n_boot - it - 1)
            print(f"  {it+1}/{a.n_boot} valid={valid} {el:.1f}min ETA {eta:.1f}min",
                  flush=True)
            # 增量落盘：即使被中断也保留部分结果（诚实标注 n_valid）
            try:
                np.savez(BS / "bootstrap_samples.partial.npz",
                         **{f"{i}_{j}_{m}": np.asarray(acc[(i, j)][m])
                            for i, j, _ in PAIRS for m in ("sp", "kt")},
                         n_valid=np.asarray([valid]))
            except Exception:
                pass
            if el > a.max_minutes:
                print("  软时限，停止"); break

    print(f"\n有效迭代 = {valid} / {a.n_boot}  ({(time.time()-t0)/60:.1f} min)")
    rows_out = []
    for i, j, lab in PAIRS:
        sp = np.asarray(acc[(i, j)]["sp"]); kt = np.asarray(acc[(i, j)]["kt"])
        def ci(v):
            return np.percentile(v, [2.5, 97.5]) if len(v) >= 20 else (np.nan, np.nan)
        ls, hs = ci(sp); lk, hk = ci(kt)
        ex_s = bool(np.isfinite(ls) and (ls > 0 or hs < 0))
        ex_k = bool(np.isfinite(lk) and (lk > 0 or hk < 0))
        rows_out.append({"pair": lab, "metric": "spearman",
                         "point": point[(i, j)][0], "ci_low": ls, "ci_high": hs,
                         "ci_excludes_zero": ex_s, "n_valid": valid})
        rows_out.append({"pair": lab, "metric": "kendall",
                         "point": point[(i, j)][1], "ci_low": lk, "ci_high": hk,
                         "ci_excludes_zero": ex_k, "n_valid": valid})
        print(f"  {lab}\n     rho={point[(i,j)][0]:+.4f} [{ls:+.4f},{hs:+.4f}] "
              f"{'✅排除0' if ex_s else '❌含0'}\n"
              f"     tau={point[(i,j)][1]:+.4f} [{lk:+.4f},{hk:+.4f}] "
              f"{'✅排除0' if ex_k else '❌含0'}")
    B = pd.DataFrame(rows_out)
    B.to_csv(BS / "rank_agreement_bootstrap.csv", index=False, encoding="utf-8")
    np.savez(BS / "bootstrap_samples.npz",
             **{f"{i}_{j}_{m}": np.asarray(acc[(i, j)][m])
                for i, j, _ in PAIRS for m in ("sp", "kt")},
             model_names=np.asarray(names))
    print("\n写出 rank_agreement_bootstrap.csv / bootstrap_samples.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
