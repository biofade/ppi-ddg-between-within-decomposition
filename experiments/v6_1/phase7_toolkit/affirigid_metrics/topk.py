# -*- coding: utf-8 -*-
"""affirigid_metrics.topk — Top-K 效用与富集指标.

设计原则
--------
* **ΔΔG 越负越好**（beneficial = 增强结合）。所有内部计算统一转成
  「分数越大越好」的方向，避免符号混淆。
* 所有函数只接受 numpy 数组 / pandas Series，**无隐藏全局状态**。
* 每个指标的语义在 docstring 里写死，便于跨论文对照。

指标定义
--------
设 `y` = 真实 ΔΔG，`p` = 预测 ΔΔG，`g` = 复合物分组。

* `best_true`  = 真值最负的前 `top_frac`（默认 10%）→ 「真强结合体」
* `p_at_k`     = 预测 top-k 中落在 `best_true` 内的比例（Precision@k）
* `ef_xpct`    = 预测 top-x% 命中 `best_true` 的数量 / 随机期望（富集因子）
* `hitk_*`     = 预测 top-k 中 ΔΔG < 0 / < −0.5 / < −1.0 的比例
* `efficiency_at_k` = `(Random − Model) / (Random − Oracle)`
  —— 在**每个复合物内**选 top-k，再对复合物做 macro 平均。
  0 = 与随机无异；1 = 达到 oracle；**负值 = 比随机更差**。
"""

from __future__ import annotations

import numpy as np

__all__ = ["topk_metrics", "efficiency_at_k", "per_structure_spearman",
           "per_structure_pearson"]

_SEED = 42


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    aa = a - a.mean()
    bb = b - b.mean()
    den = np.sqrt((aa * aa).sum() * (bb * bb).sum())
    return float((aa * bb).sum() / den) if den > 0 else float("nan")


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.empty(len(a))
    ra[np.argsort(a, kind="mergesort")] = np.arange(len(a))
    rb = np.empty(len(b))
    rb[np.argsort(b, kind="mergesort")] = np.arange(len(b))
    return _pearson(ra, rb)


def per_structure_spearman(y, p, groups, min_n: int = 5) -> float:
    """每个复合物内单独算 Spearman，再 macro 平均。

    **这是「靶点内排序能力」的核心指标。** `min_n` 以下的复合物被跳过
    （组内过少时 Spearman 极不稳定）。
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    g = np.asarray(groups)
    vals = []
    for c in np.unique(g):
        m = g == c
        if int(m.sum()) < min_n:
            continue
        v = _spearman(y[m], p[m])
        if v == v:
            vals.append(v)
    return float(np.mean(vals)) if vals else float("nan")


def per_structure_pearson(y, p, groups, min_n: int = 5) -> float:
    """每个复合物内单独算 Pearson，再 macro 平均。"""
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    g = np.asarray(groups)
    vals = []
    for c in np.unique(g):
        m = g == c
        if int(m.sum()) < min_n:
            continue
        v = _pearson(y[m], p[m])
        if v == v:
            vals.append(v)
    return float(np.mean(vals)) if vals else float("nan")


def efficiency_at_k(y, p, groups, k: int = 10, n_rand: int = 200,
                    seed: int = _SEED, min_group_n: int = 5) -> float:
    """`Efficiency@K = (Random − Model) / (Random − Oracle)`（复合物内，macro 平均）。

    参数
    ----
    y, p    : 真实值与预测值（ΔΔG，**越负越好**）
    groups  : 复合物 / 靶点 ID（分组键）
    k       : 每个复合物的实验预算；若组内样本少于 k，取 `min(5, n)`
    n_rand  : 随机基线抽样次数（越多越稳；200 为默认）
    min_group_n : 组内样本下限，低于该值的组被跳过

    返回
    ----
    float；0 = 与随机无异，1 = 达到 oracle，<0 = 比随机更差。
    若无可评估的组，返回 NaN。
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    g = np.asarray(groups)
    rng = np.random.RandomState(seed)
    eff = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < min_group_n:
            continue
        kk = min(k, len(m)) if len(m) >= k else min(5, len(m))
        sy, sp = y[m], p[m]
        sel = np.argpartition(sp, kk - 1)[:kk] if kk < len(m) else np.arange(len(m))
        orc = np.argpartition(sy, kk - 1)[:kk] if kk < len(m) else np.arange(len(m))
        v_sel = float(-sy[sel].mean())
        v_or = float(-sy[orc].mean())
        rnd_idx = rng.randint(0, len(m), size=(n_rand, kk))
        v_rnd = float(-sy[rnd_idx].mean())
        den = v_rnd - v_or
        if abs(den) > 1e-9:
            eff.append((v_rnd - v_sel) / den)
    return float(np.mean(eff)) if eff else float("nan")


def topk_metrics(y, p, groups=None, top_frac: float = 0.10,
                 n_rand: int = 200, seed: int = _SEED) -> dict:
    """一次算全套 Top-K / 富集指标 + 可选的 within-group Efficiency@K。

    返回 dict，键名见模块 docstring。所有「越多越好」的量已在内部统一方向。
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    n = len(y)
    out: dict[str, float] = {"n": int(n)}
    if n < 10:
        return out
    score = -y            # 越大越好
    pred = -p
    best_true = score >= np.quantile(score, 1.0 - top_frac)
    out["base_rate_best_true"] = float(best_true.mean())

    for k in (5, 10, 20):
        if k > n:
            continue
        idx = np.argsort(-pred, kind="mergesort")[:k]
        out[f"p_at_{k}"] = float(best_true[idx].mean())
        out[f"hit{k}_beneficial"] = float((y[idx] < 0).mean())
        out[f"hit{k}_ddg_m05"] = float((y[idx] < -0.5).mean())
        out[f"hit{k}_ddg_m10"] = float((y[idx] < -1.0).mean())
        out[f"sel{k}_mean_ddg"] = float(y[idx].mean())
        out[f"sel{k}_best_ddg"] = float(y[idx].min())
    out["mean_ddg_all"] = float(y.mean())
    out["oracle10_mean_ddg"] = float(np.sort(y)[:min(10, n)].mean())

    for pct, name in ((1.0, "1pct"), (5.0, "5pct"), (10.0, "10pct")):
        k = max(1, int(round(n * pct / 100.0)))
        idx = np.argsort(-pred, kind="mergesort")[:k]
        hits = int(best_true[idx].sum())
        exp = k * float(best_true.mean())
        out[f"ef_{name}"] = float(hits / exp) if exp > 0 else float("nan")
        out[f"ef_{name}_hits"] = hits
        out[f"ef_{name}_k"] = k

    if groups is not None:
        out["efficiency_at_10"] = efficiency_at_k(y, p, groups, k=10,
                                                  n_rand=n_rand, seed=seed)
        out["macro_pdb_spearman"] = per_structure_spearman(y, p, groups)
        out["macro_pdb_pearson"] = per_structure_pearson(y, p, groups)
    return out
