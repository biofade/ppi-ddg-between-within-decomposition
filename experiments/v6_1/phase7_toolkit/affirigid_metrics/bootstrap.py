# -*- coding: utf-8 -*-
"""affirigid_metrics.bootstrap — cluster bootstrap 置信区间与配对检验.

核心概念：**cluster bootstrap**
--------------------------------
在 PPI ΔΔG 任务里，同一个复合物的多个突变体**不独立**。
把 N 行当成 N 个独立样本会**低估置信区间**。
正确做法是把「复合物」当作抽样单位：每次有放回地抽复合物，
并携带该复合物的**全部**突变行。

本模块提供：
* `cluster_bootstrap_ci`   —— 单模型某指标的 95% CI
* `cluster_paired_bootstrap` —— 两模型在某指标上的配对差 CI 与 p 值
* `paired_wilcoxon`        —— 配对 Wilcoxon 符号秩（跨复合物的配对均值）
"""

from __future__ import annotations

from typing import Callable

import numpy as np

__all__ = ["cluster_bootstrap_ci", "cluster_paired_bootstrap",
           "paired_wilcoxon", "split_by_group"]

_SEED = 42


def split_by_group(groups) -> tuple[np.ndarray, list[np.ndarray]]:
    """返回 (unique_groups, [每组的行索引])。顺序按组名排序，保证确定性。"""
    g = np.asarray(groups)
    order = np.argsort(g, kind="mergesort")
    gs = g[order]
    uniq, starts = np.unique(gs, return_index=True)
    idxs = []
    for i, s in enumerate(starts):
        e = starts[i + 1] if i + 1 < len(starts) else len(gs)
        idxs.append(order[s:e])
    return uniq, idxs


def cluster_bootstrap_ci(metric: Callable, y, p, groups,
                         n_boot: int = 5000, seed: int = _SEED,
                         min_groups: int = 5) -> dict:
    """cluster bootstrap 95% CI。

    `metric` 签名必须为 `metric(y, p, groups) -> float`
    （若为二元 `metric(y, p)` 也可，会自动降级）。
    `min_groups`：某次重采样若 unique group 少于此值则跳过该次迭代。
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    uniq, cidx = split_by_group(groups)
    n_cl = len(uniq)
    pt = float(_call(metric, y, p, groups))
    if n_cl < 3:
        return {"point": pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "sd": float("nan"), "n_clusters": n_cl, "n_valid": 0,
                "n_skipped": 0}
    rs = np.random.RandomState(seed)
    vals, skipped = [], 0
    for _ in range(n_boot):
        pick = rs.randint(0, n_cl, n_cl)
        if len(np.unique(pick)) < min_groups:
            skipped += 1
            continue
        rows = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i, dtype=np.int64)
                              for i, j in enumerate(pick)])
        v = _call(metric, y[rows], p[rows], rep)
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return {"point": pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "sd": float("nan"), "n_clusters": n_cl, "n_valid": len(vals),
                "n_skipped": skipped}
    a = np.asarray(vals, float)
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": pt, "ci_low": float(lo), "ci_high": float(hi),
            "sd": float(a.std(ddof=1)), "n_clusters": n_cl,
            "n_valid": len(vals), "n_skipped": skipped,
            "p_gt_zero": float((a > 0).mean())}


def cluster_paired_bootstrap(metric: Callable, y, pa, pb, groups,
                             n_boot: int = 5000, seed: int = _SEED,
                             min_groups: int = 5) -> dict:
    """两模型在某指标上的**配对** cluster bootstrap：Δ = metric(b) − metric(a)。"""
    y = np.asarray(y, float)
    pa = np.asarray(pa, float)
    pb = np.asarray(pb, float)
    uniq, cidx = split_by_group(groups)
    n_cl = len(uniq)
    ma = float(_call(metric, y, pa, groups))
    mb = float(_call(metric, y, pb, groups))
    obs = mb - ma
    if n_cl < 3:
        return {"delta": obs, "ci_low": float("nan"), "ci_high": float("nan"),
                "p_two_sided": float("nan"), "n_clusters": n_cl, "n_valid": 0}
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        pick = rs.randint(0, n_cl, n_cl)
        if len(np.unique(pick)) < min_groups:
            continue
        rows = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i, dtype=np.int64)
                              for i, j in enumerate(pick)])
        va = _call(metric, y[rows], pa[rows], rep)
        vb = _call(metric, y[rows], pb[rows], rep)
        if va == va and vb == vb:
            d.append(vb - va)
    if len(d) < 20:
        return {"delta": obs, "ci_low": float("nan"), "ci_high": float("nan"),
                "p_two_sided": float("nan"), "n_clusters": n_cl,
                "n_valid": len(d)}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    dc = d - d.mean()
    pl = (np.sum(dc <= -abs(obs)) + 1) / (len(d) + 1)
    pr = (np.sum(dc >= abs(obs)) + 1) / (len(d) + 1)
    return {"delta": float(obs), "ci_low": float(lo), "ci_high": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl, pr))),
            "n_clusters": n_cl, "n_valid": len(d)}


def paired_wilcoxon(a, b) -> dict:
    """配对 Wilcoxon 符号秩检验（`a`、`b` 为同一批单位的配对值）。

    **用于「同一批复合物上的两类模型」**——此时不能用 Mann-Whitney U
    （两组不独立）。
    """
    from scipy.stats import wilcoxon
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = a - b
    n = len(d)
    if n < 5 or np.allclose(d, 0):
        return {"W": float("nan"), "p": float("nan"), "n": n,
                "delta_mean": float(d.mean()) if n else float("nan")}
    w, p = wilcoxon(a, b)
    return {"W": float(w), "p": float(p), "n": n,
            "delta_mean": float(d.mean()),
            "n_pos": int((d > 0).sum()), "n_neg": int((d < 0).sum())}


def _call(metric: Callable, y, p, groups):
    """兼容 3 参数与 2 参数 metric。"""
    try:
        return metric(y, p, groups)
    except TypeError:
        return metric(y, p)
