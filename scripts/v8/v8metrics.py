# -*- coding: utf-8 -*-
"""scripts/v8/v8metrics.py -- V8 三轴能力度量引擎（CPU-first，向量化）。

**方向约定（§36）**：全项目 `ΔG = 0.592*ln(KD[M])`，`ΔΔG > 0 = 结合变弱`。
因此「有利用户的突变」= **ΔΔG 更负**。内部统一转为
`utility = -ddG`（越大越有利），所有 Top-K/NDCG/EF 基于 utility。

三轴：
  A 校准  : Pearson/Spearman/RMSE/MAE/target-mean/attenuation/C_B(row,equal)
  B 靶内排序: MacroSpearman(n>=5)/MacroKendall/Pairwise(delta)/C_W
  C 选择效用: Hit@K / Precision@K / Recall@K / Efficiency@K / NDCG@K / EF@pct

全部为确定性实现；bootstrap 时用 `prep_universe()` 预计算目标索引后复用。
"""

from __future__ import annotations

import numpy as np
from scipy.stats import kendalltau, spearmanr

PRIMARY_MIN_N = 5
K_PRIMARY = 10


# ---------------------------------------------------------------- 基础
def _pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 2:
        return np.nan
    sa, sb = a.std(), b.std()
    if sa <= 0 or sb <= 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def _spearman(a, b):
    if len(a) < 2:
        return np.nan
    r = spearmanr(a, b).correlation
    return float(r) if np.isfinite(r) else np.nan


def _kendall_tau_b(a, b):
    if len(a) < 2:
        return np.nan
    r = kendalltau(a, b, variant="b").correlation
    return float(r) if np.isfinite(r) else np.nan


def _pairwise_acc(y, p, delta=0.0):
    """同靶内两两比较；|Δy| <= delta 的真值对按预设排除（delta>0 时）。
    delta=0 时保留全部非零真值差；预测并列计 0.5。"""
    n = len(y)
    if n < 2:
        return np.nan
    iu = np.triu_indices(n, 1)
    dy = y[iu[0]] - y[iu[1]]
    dp = p[iu[0]] - p[iu[1]]
    keep = np.abs(dy) > delta if delta > 0 else (dy != 0)
    if not keep.any():
        return np.nan
    dyv = dy[keep]; dpv = dp[keep]
    agree = (np.sign(dyv) == np.sign(dpv)) & (dpv != 0)
    tie = (dpv == 0)
    return float((agree.sum() + 0.5 * tie.sum()) / dyv.size)


def prep_universe(y, p, target, family=None, min_n=PRIMARY_MIN_N):
    """预计算目标索引（bootstrap 复用，避免每次 groupby）。"""
    tgt = np.asarray(target).astype(str)
    uq = np.unique(tgt)
    idx = {t: np.where(tgt == t)[0] for t in uq}
    keep = [t for t in uq if idx[t].size >= min_n]
    fam = None
    if family is not None:
        fam = {}
        fa = np.asarray(family).astype(str)
        for t in uq:
            fam[t] = fa[idx[t][0]]
    return {"y": np.asarray(y, float), "p": np.asarray(p, float),
            "target": tgt, "uq": uq, "idx": idx,
            "eligible": keep, "family": fam, "min_n": min_n}


# ---------------------------------------------------------------- 轴 A：校准
def axis_calibration(U, weights=None):
    y, p = U["y"], U["p"]
    out = {"overall_pearson": _pearson(y, p),
           "overall_spearman": _spearman(y, p),
           "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
           "mae": float(np.mean(np.abs(y - p)))}
    # 目标均值
    uq = U["uq"]; idx = U["idx"]
    By = np.array([y[idx[t]].mean() for t in uq])
    Bp = np.array([p[idx[t]].mean() for t in uq])
    out["target_mean_pearson"] = _pearson(By, Bp)
    out["target_mean_spearman"] = _spearman(By, Bp)
    # 衰减斜率（预测目标均值 ~ 真值目标均值）
    if By.std() > 0:
        out["attenuation_slope"] = float(np.polyfit(By, Bp, 1)[0])
    else:
        out["attenuation_slope"] = np.nan
    # C_B / C_W（§28 严格分解）
    #   Cov(y,p) = Σ_t n_t (B_y,t - ȳ)(B_p,t - p̄)  +  Σ_t Σ_i W_y W_p
    #   ȳ, p̄ 与 B 必须用**同一套**（行级）整体均值；跨项 Σ n_t (B-ȳ) ≡ 0。
    n = y.size
    cnt = np.array([idx[t].size for t in uq], float)
    sd = y.std() * p.std()
    Wy = np.empty_like(y); Wp = np.empty_like(p)
    for t in uq:
        i = idx[t]
        Wy[i] = y[i] - y[i].mean(); Wp[i] = p[i] - p[i].mean()
    if sd > 0:
        ybar, pbar = y.mean(), p.mean()
        C_Bnum = float(np.sum(cnt * (By - ybar) * (Bp - pbar)))
        C_Wnum = float(np.sum(Wy * Wp))
        out["C_B_row"] = C_Bnum / (n * sd)
        out["C_W_row"] = C_Wnum / (n * sd)
        out["C_B_equal"] = float(
            np.mean((By - By.mean()) * (Bp - Bp.mean())) / sd)
    else:
        out["C_B_row"] = out["C_W_row"] = out["C_B_equal"] = np.nan
    # 分解 QC：r_full ≈ C_B + C_W（恒等式，残差应为数值 0）
    out["decomp_residual"] = float(out["overall_pearson"]
                                   - (out["C_B_row"] + out["C_W_row"]))
    return out


# ---------------------------------------------------------------- 轴 B：靶内排序
def axis_ranking(U):
    uq = U["uq"]; idx = U["idx"]; eligible = U["eligible"]
    y, p = U["y"], U["p"]
    sp, kt, pw0, pw2, pw5 = [], [], [], [], []
    for t in eligible:
        i = idx[t]
        sp.append(_spearman(y[i], p[i]))
        kt.append(_kendall_tau_b(y[i], p[i]))
        pw0.append(_pairwise_acc(y[i], p[i], 0.0))
        pw2.append(_pairwise_acc(y[i], p[i], 0.2))
        pw5.append(_pairwise_acc(y[i], p[i], 0.5))

    def _nm(v):
        v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
        return float(v.mean()) if v.size else np.nan
    cal = axis_calibration(U)
    return {"macro_spearman": _nm(sp), "macro_kendall": _nm(kt),
            "pairwise_d0": _nm(pw0), "pairwise_d02": _nm(pw2),
            "pairwise_d05": _nm(pw5), "C_W": cal["C_W_row"],
            "n_eligible_targets": len(eligible)}


# ---------------------------------------------------------------- 轴 C：选择效用
def _topk_from_ranking(util, k):
    """确定性 Top-K：按 utility 降序，**稳定排序**（tie 用原始下标）。"""
    n = len(util)
    k = min(k, n)
    order = np.lexsort((np.arange(n), -util))   # 主键 -util，次键原下标
    return order[:k]


def axis_selection(U, ks=(1, 5, 10), ef_pcts=(1, 5)):
    """轴 C：选择效用。

    方向（§36）：`ΔΔG` 越负越有利。预测 `tp` **越小** = 越有利。
    真值 gain：`g_i = 1 + (n - rank_i)`，`rank_i` 为真值 utility 降序名次
    （最优 gain 最大）。`ideal_rank` 为按 gain 降序的确定性顺序。
    """
    uq = U["uq"]; idx = U["idx"]; eligible = U["eligible"]
    y, p = U["y"], U["p"]
    ks = tuple(ks); ef_pcts = tuple(ef_pcts)
    acc = {f"hit@{k}": [] for k in ks}
    acc.update({f"precision@{k}": [] for k in ks})
    acc.update({f"recall@{k}": [] for k in ks})
    acc.update({f"efficiency@{k}": [] for k in ks})
    acc.update({f"ndcg@{k}": [] for k in ks})
    acc.update({f"ef@{pc}": [] for pc in ef_pcts})
    for t in eligible:
        i = idx[t]
        ty = y[i]; tp = p[i]
        n = len(i)
        ar = np.arange(n)
        # 真值名次（utility 降序，稳定）
        true_order = np.lexsort((ar, -(-ty)))          # utility=-ty 降序
        rank = np.empty(n, int); rank[true_order] = ar  # 0 = 最优
        gain = 1.0 + (n - rank)                        # 1..n
        ideal_order = np.lexsort((ar, -gain))          # gain 降序
        # 预测选择（tp 升序）
        pred_order = np.lexsort((ar, tp))
        for k in ks:
            kk = max(1, min(k, n))
            sel = pred_order[:kk]
            true_top = set(true_order[:kk].tolist())
            acc[f"hit@{k}"].append(1.0 if true_order[0] in set(sel.tolist()) else 0.0)
            inter = len(set(sel.tolist()) & true_top)
            acc[f"precision@{k}"].append(inter / kk)
            acc[f"recall@{k}"].append(inter / len(true_top))
            # Efficiency@K = 选中 gain 和 / 理想前 K gain 和
            gsel = gain[sel].sum(); gid = gain[ideal_order[:kk]].sum()
            acc[f"efficiency@{k}"].append(float(gsel / gid) if gid > 0 else np.nan)
            # NDCG@K = DCG@K / IDCG@K（同一 gain，log2 折扣）
            disc = np.log2(np.arange(2, kk + 2))
            dcg = (gain[sel] / disc).sum()
            idcg = (gain[ideal_order[:kk]] / disc).sum()
            acc[f"ndcg@{k}"].append(float(dcg / idcg) if idcg > 0 else np.nan)
        for pc in ef_pcts:
            kk = max(1, int(round(n * pc / 100.0)))
            sel = pred_order[:kk]
            true_top = set(true_order[:kk].tolist())
            frac = len(set(sel.tolist()) & true_top) / kk
            acc[f"ef@{pc}"].append(float(frac / (kk / n)))

    def _nm(v):
        v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
        return float(v.mean()) if v.size else np.nan
    return {k: _nm(v) for k, v in acc.items()}


# ---------------------------------------------------------------- 汇总
def compute_all(U):
    A = axis_calibration(U)
    B = axis_ranking(U)
    C = axis_selection(U)
    out = {}
    out.update(A); out.update(B); out.update(C)
    out["n_rows"] = int(U["y"].size)
    out["n_targets"] = int(len(U["uq"]))
    out["n_targets_ge5"] = int(len(U["eligible"]))
    out["min_n"] = U["min_n"]
    return out


METRIC_COLS = [
    # 校准
    "overall_pearson", "overall_spearman", "rmse", "mae",
    "target_mean_pearson", "target_mean_spearman", "attenuation_slope",
    "C_B_row", "C_W_row", "C_B_equal",
    # 排序
    "macro_spearman", "macro_kendall", "pairwise_d0", "pairwise_d02",
    "pairwise_d05", "C_W",
    # 选择
    "hit@1", "hit@5", "hit@10", "precision@5", "precision@10",
    "recall@5", "recall@10", "efficiency@5", "efficiency@10",
    "ndcg@5", "ndcg@10", "ef@1", "ef@5",
]
