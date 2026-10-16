# -*- coding: utf-8 -*-
"""metrics_ranking.py -- 面向湿实验富集筛选的排序评估指标库 (V2.0).

为什么要换主指标
----------------
V1.1 以 R²=0.261 作为主指标, 但该模型测试集拟合斜率仅 0.24 —— 它是个**排序器**,
不是定量器。用 R² 评价排序器会系统性低估其价值, 也无法回答湿实验真正关心的问题:
「如果我只做得起 20 个克隆, 这批候选里富集了多少真强结合体?」

本模块提供:
* 相关性: Spearman ρ、Kendall τ (含 tie 处理)
* 富集:   Precision@k (k=5/10/20)、Recall@k、Enrichment Factor EF@x% (x=1/5/10)
* 排序质量: NDCG@k (k=5/10/20), 支持 graded / binary 两种 gain
* 不确定性: bootstrap 95% CI (可多 seed 重复验证稳定性)
* 模型比较: **配对 bootstrap** 检验 (同一测试集上两个模型的指标差), 及 BH-FDR 校正

明确不用什么
------------
* **不用 DeLong test**。DeLong 检验的是两条 ROC 曲线下面积 (AUC) 的差异,
  它假设是二分类 AUC 且基于 U 统计量的方差协方差结构。把它套到 Spearman 相关系数
  或富集因子上是**统计学误用**。模型间比较一律走配对 bootstrap:
  在同样的样本索引上重采样, 分别算两个模型的指标, 取差值分布的分位数作 CI。
* 不引入 statsmodels / sklearn 之外的重依赖; BH-FDR 自实现。

约定 (非常重要)
---------------
**所有输入统一为 ``score``, 约定 score 越大 = 结合越强。**
原始 ΔG 是"越负越强", 调用方必须先取负号 (本模块的 ``negate_dG`` 工具函数)。
方向搞反会让所有富集指标塌成随机基线, 因此模块入口做了显式方向说明。

缺失值处理
----------
(y, score) 中任一侧为 NaN 的行**整行丢弃**, 并在返回值里记录 ``n_dropped``。
并列值 (ties) 用 **average rank** 处理: Spearman 走 scipy 的 tie 校正,
Precision@k 遇到第 k 名并列时按 ``k`` 名整取 (不插值), 这是排序评估的通行做法,
并在返回值里用 ``k_effective`` 标明实际取了多少条。

纪律: 纯 CPU; 无网络; 固定种子; 全部函数对空/极小样本输入返回 NaN 而不是抛错。
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

SEED = 42
DEFAULT_BOOT = 2000
DEFAULT_K = (5, 10, 20)
DEFAULT_EF_PCT = (1.0, 5.0, 10.0)

logger = logging.getLogger("metrics_ranking")

MetricFn = Callable[[np.ndarray, np.ndarray], float]


# ===================================================================== 清理工具
def negate_dG(dg: Sequence[float] | np.ndarray) -> np.ndarray:
    """ΔG (越负越强) -> score (越大越强)。**必须**在构造 score 时调用。"""
    return -np.asarray(dg, dtype=float)


def _clean(y: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    if y.shape != score.shape:
        raise ValueError(f"y 与 score 形状不一致: {y.shape} vs {score.shape}")
    ok = np.isfinite(y) & np.isfinite(score)
    n_drop = int((~ok).sum())
    return y[ok], score[ok], n_drop


def _average_ranks(x: np.ndarray) -> np.ndarray:
    """average rank (1-based); 并列取平均秩。"""
    from scipy.stats import rankdata
    return rankdata(x, method="average")


# ===================================================================== 相关性
def spearman_rho(y: np.ndarray, score: np.ndarray) -> float:
    """Spearman ρ。tie 由 scipy 统一校正; 有效样本 <3 或全并列返回 NaN。"""
    y, score, _ = _clean(y, score)
    if len(y) < 3 or len(np.unique(y)) < 2 or len(np.unique(score)) < 2:
        return float("nan")
    from scipy.stats import spearmanr
    c = spearmanr(y, score).correlation
    return float(c) if c == c else float("nan")


def kendall_tau(y: np.ndarray, score: np.ndarray) -> float:
    """Kendall τ-b (对并列值稳健)。"""
    y, score, _ = _clean(y, score)
    if len(y) < 3 or len(np.unique(y)) < 2 or len(np.unique(score)) < 2:
        return float("nan")
    from scipy.stats import kendalltau
    c = kendalltau(y, score).correlation
    return float(c) if c == c else float("nan")


# ===================================================================== 相关阈值
def relevance_thresholds(y: np.ndarray, top_frac: float) -> float:
    """给定"真强结合"占比, 返回 y 的阈值 (y 越大越强)。"""
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) == 0:
        return float("nan")
    if top_frac <= 0:
        return float("inf")
    if top_frac >= 1:
        return float("-inf")
    k = max(1, int(np.ceil(top_frac * len(y))))
    srt = np.sort(y)[::-1]
    return float(srt[k - 1])


def binary_relevance(y: np.ndarray, top_frac: float) -> np.ndarray:
    """把连续 y 二值化成 relevance: 真值最高的前 top_frac 比例记 1。"""
    y = np.asarray(y, dtype=float)
    thr = relevance_thresholds(y, top_frac)
    rel = (y >= thr).astype(float)
    rel[~np.isfinite(y)] = np.nan
    return rel


def graded_relevance(y: np.ndarray, n_levels: int = 5,
                     gain: Sequence[float] | None = None) -> np.ndarray:
    """把连续 y 按**真值**分位离散成 n_levels 级 relevance (越大越强 -> 级别越高)。

    与 LambdaRank 的 ``label_gain`` 口径一致: 默认 [0,1,3,7,15] 指数增益。
    """
    y = np.asarray(y, dtype=float)
    if gain is None:
        gain = (0, 1, 3, 7, 15)[:n_levels]
    gain = np.asarray(gain, dtype=float)
    n_levels = len(gain)
    out = np.full(len(y), np.nan, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < n_levels:
        return out
    ranks = _average_ranks(y[ok])          # 1..n, 越大 = y 越大
    q = (ranks - 1) / max(len(ranks) - 1, 1)   # 0..1
    lvl = np.minimum((q * n_levels).astype(int), n_levels - 1)
    out[ok] = gain[lvl]
    return out


# ===================================================================== 富集
def _topk_expected_hits(rel: np.ndarray, score: np.ndarray, k: int) -> tuple[float, int]:
    """预测前 k 名中的**期望**命中数, 正确处理并列值。

    若不加处理, ``np.argsort(-score)`` 在**分数完全相同**时会按原始行序取前 k 名 ——
    这会让一个"完全无区分能力"的模型凭空得到 2 倍富集 (自检实测 EF=2.000 而非 1.0)。
    正确做法是把并列视为随机排列, 取期望:

    * 分数严格高于第 k 名阈值的条目: 全部计入;
    * 处于阈值的并列块: 按 ``剩余名额 / 并列块大小`` 的比例折算期望命中数。

    这是排序评估的标准处理 (与 ``sklearn.metrics.ndcg_score`` 对并列的期望口径一致),
    也保证 ``P@k``/``EF`` 在无区分度模型上收敛到基准率。
    """
    rel = np.asarray(rel, dtype=float)
    score = np.asarray(score, dtype=float)
    n = len(score)
    if n == 0 or k <= 0:
        return 0.0, 0
    kk = min(int(k), n)
    order = np.argsort(-score, kind="mergesort")
    s_sorted = score[order]
    r_sorted = rel[order]
    s_thr = s_sorted[kk - 1]
    above = s_sorted > s_thr
    n_above = int(above.sum())
    hits = float(np.nansum(r_sorted[above]))
    remain = kk - n_above
    if remain > 0:
        tie_mask = s_sorted == s_thr
        n_tie = int(tie_mask.sum())
        if n_tie > 0:
            tie_sum = float(np.nansum(r_sorted[tie_mask]))
            hits += tie_sum * (remain / n_tie)
    return hits, kk


def precision_at_k(y: np.ndarray, score: np.ndarray, k: int,
                   top_frac: float = 0.10) -> float:
    """Precision@k = 预测前 k 名中真强结合体 (真值前 top_frac) 的**期望**比例。

    并列分数按随机排列取期望 (见 ``_topk_expected_hits``)。
    """
    y, score, _ = _clean(y, score)
    if len(y) == 0 or k <= 0:
        return float("nan")
    rel = binary_relevance(y, top_frac)
    hits, kk = _topk_expected_hits(rel, score, k)
    return float(hits / kk) if kk else float("nan")


def recall_at_k(y: np.ndarray, score: np.ndarray, k: int,
                top_frac: float = 0.10) -> float:
    y, score, _ = _clean(y, score)
    if len(y) == 0 or k <= 0:
        return float("nan")
    rel = binary_relevance(y, top_frac)
    n_pos = int(np.nansum(rel))
    if n_pos == 0:
        return float("nan")
    hits, _ = _topk_expected_hits(rel, score, k)
    return float(hits / n_pos)


def enrichment_factor(y: np.ndarray, score: np.ndarray, pct: float,
                      top_frac: float = 0.10) -> float:
    """EF@pct = (预测前 pct% 中真强结合体的期望占比) / (全集真强结合体占比)。

    EF = 1 表示与随机抽取无异; EF = 5 表示同样预算下强结合体富集 5 倍。
    真强结合体定义为**真值**前 ``top_frac`` 比例 (默认 10%)。

    注意 EF 在小 k 下是**离散且高方差**的: n=950、pct=1% 时 k=10, 基准率 10%,
    故 EF 的最小非零步长约 1.05。解读 EF@1% 必须同时看它自己的 bootstrap CI。
    """
    y, score, _ = _clean(y, score)
    if len(y) == 0 or pct <= 0:
        return float("nan")
    rel = binary_relevance(y, top_frac)
    base = float(np.nanmean(rel))
    if not np.isfinite(base) or base <= 0:
        return float("nan")
    kk = max(1, int(np.ceil(pct / 100.0 * len(y))))
    hits, kk = _topk_expected_hits(rel, score, kk)
    return float((hits / kk) / base)


def ndcg_at_k(y: np.ndarray, score: np.ndarray, k: int,
              mode: str = "graded", top_frac: float = 0.10,
              n_levels: int = 5) -> float:
    """NDCG@k。``mode='graded'`` 用 5 级指数增益, ``mode='binary'`` 用二值相关。

    并列分数同样按期望处理: 对并列块用其平均折损因子, 避免行序伪影。
    """
    y, score, _ = _clean(y, score)
    if len(y) == 0 or k <= 0:
        return float("nan")
    if mode == "binary":
        rel = binary_relevance(y, top_frac)
        rel = np.nan_to_num(rel, nan=0.0)
    elif mode == "graded":
        rel = graded_relevance(y, n_levels=n_levels)
        rel = np.nan_to_num(rel, nan=0.0)
    else:
        raise ValueError(f"未知 mode: {mode!r} (用 'graded' 或 'binary')")
    n = len(y)
    kk = min(int(k), n)
    disc = 1.0 / np.log2(np.arange(2, n + 2))
    order = np.argsort(-score, kind="mergesort")
    r_sorted = rel[order]
    s_sorted = score[order]
    d_sorted = disc[:n]
    # 严格高于阈值的部分全取; 并列块按平均折损因子折算
    s_thr = s_sorted[kk - 1]
    above = s_sorted > s_thr
    dcg = float(np.sum(r_sorted[above] * d_sorted[above]))
    remain = kk - int(above.sum())
    if remain > 0:
        tie_mask = s_sorted == s_thr
        n_tie = int(tie_mask.sum())
        if n_tie > 0:
            d_tie = float(np.mean(d_sorted[tie_mask]))
            dcg += float(np.sum(r_sorted[tie_mask])) * d_tie * (remain / n_tie)
    ideal = np.sort(rel)[::-1]
    idcg = float(np.sum(ideal[:kk] * disc[:kk]))
    if idcg <= 0:
        return float("nan")
    return float(dcg / idcg)


# ===================================================================== 全套指标
@dataclass
class RankingMetrics:
    """一个模型的完整排序/富集指标集合。所有键名可直接进 CSV。"""

    n: int = 0
    n_dropped: int = 0
    n_positive: int = 0
    top_frac: float = 0.10
    ks: tuple[int, ...] = DEFAULT_K
    ef_pcts: tuple[float, ...] = DEFAULT_EF_PCT
    values: dict[str, float] = field(default_factory=dict)
    k_effective: dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {"n": self.n, "n_dropped": self.n_dropped,
                "n_positive": self.n_positive, "top_frac": self.top_frac,
                **self.values}

    def summary_lines(self) -> list[str]:
        out = [f"n={self.n} (丢弃 {self.n_dropped}), 真强结合体={self.n_positive} "
               f"(真值前 {self.top_frac * 100:.0f}%)"]
        for k, v in self.values.items():
            out.append(f"  {k:<22s} {v:+.4f}" if isinstance(v, float) else f"  {k}: {v}")
        return out


def compute_ranking_metrics(
    y: np.ndarray,
    score: np.ndarray,
    ks: Sequence[int] = DEFAULT_K,
    ef_pcts: Sequence[float] = DEFAULT_EF_PCT,
    top_frac: float = 0.10,
    ndcg_mode: str = "graded",
) -> RankingMetrics:
    """计算全套指标。``y`` 越大越强, ``score`` 越大越强。"""
    y_c, s_c, n_drop = _clean(y, score)
    rel = binary_relevance(y_c, top_frac) if len(y_c) else np.array([])
    m = RankingMetrics(
        n=int(len(y_c)), n_dropped=n_drop,
        n_positive=int(np.nansum(rel)) if len(rel) else 0,
        top_frac=top_frac, ks=tuple(ks), ef_pcts=tuple(ef_pcts),
    )
    m.values["spearman"] = spearman_rho(y_c, s_c)
    m.values["kendall_tau"] = kendall_tau(y_c, s_c)
    m.values["ndcg_mode"] = ndcg_mode  # type: ignore[assignment]

    for k in ks:
        m.values[f"precision@{k}"] = precision_at_k(y_c, s_c, k, top_frac)
        m.values[f"recall@{k}"] = recall_at_k(y_c, s_c, k, top_frac)
        m.values[f"ndcg@{k}"] = ndcg_at_k(y_c, s_c, k, mode=ndcg_mode, top_frac=top_frac)
        m.k_effective[f"precision@{k}"] = min(int(k), len(y_c))
    for p in ef_pcts:
        key = f"ef@{p:g}%"
        m.values[key] = enrichment_factor(y_c, s_c, p, top_frac)
        m.k_effective[key] = min(int(np.ceil(p / 100.0 * len(y_c))), len(y_c))
    return m


# ===================================================================== bootstrap
def bootstrap_ci(
    metric_fn: MetricFn,
    y: np.ndarray,
    score: np.ndarray,
    n_boot: int = DEFAULT_BOOT,
    seed: int = SEED,
    alpha: float = 0.05,
    min_valid_frac: float = 0.5,
) -> tuple[float, float, float]:
    """bootstrap 百分位 95% CI -> (lo, hi, point)。

    ``point`` 是**全样本**上的指标值 (不是 bootstrap 均值), 与 CI 一起报告更诚实。
    有效重采样不足 ``min_valid_frac`` 时返回 (nan, nan, point)。
    """
    y_c, s_c, _ = _clean(y, score)
    point = float(metric_fn(y_c, s_c)) if len(y_c) else float("nan")
    if len(y_c) < 5:
        return float("nan"), float("nan"), point
    rs = np.random.RandomState(seed)
    vals = []
    N = len(y_c)
    for _ in range(n_boot):
        idx = rs.randint(0, N, N)
        v = metric_fn(y_c[idx], s_c[idx])
        if v == v:
            vals.append(v)
    if len(vals) < max(20, int(min_valid_frac * n_boot)):
        return float("nan"), float("nan"), point
    lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi), point


def bootstrap_ci_multi_seed(
    metric_fn: MetricFn,
    y: np.ndarray,
    score: np.ndarray,
    seeds: Sequence[int] = (42, 43, 44, 45, 46),
    n_boot: int = 1000,
) -> dict[str, Any]:
    """多 seed bootstrap: 量化 CI 本身对随机种子的敏感性 (稳定性检查)。"""
    los, his, pts = [], [], []
    for s in seeds:
        lo, hi, pt = bootstrap_ci(metric_fn, y, score, n_boot=n_boot, seed=int(s))
        los.append(lo)
        his.append(hi)
        pts.append(pt)
    los_a, his_a = np.asarray(los, float), np.asarray(his, float)
    return {
        "point": float(np.nanmean(pts)) if pts else float("nan"),
        "ci_low_mean": float(np.nanmean(los_a)) if len(los_a) else float("nan"),
        "ci_high_mean": float(np.nanmean(his_a)) if len(his_a) else float("nan"),
        "ci_low_min": float(np.nanmin(los_a)) if len(los_a) else float("nan"),
        "ci_low_max": float(np.nanmax(los_a)) if len(los_a) else float("nan"),
        "ci_high_min": float(np.nanmin(his_a)) if len(his_a) else float("nan"),
        "ci_high_max": float(np.nanmax(his_a)) if len(his_a) else float("nan"),
        "n_seeds": len(seeds),
        "n_boot_per_seed": n_boot,
    }


def paired_bootstrap_test(
    metric_fn: MetricFn,
    y: np.ndarray,
    score_a: np.ndarray,
    score_b: np.ndarray,
    n_boot: int = DEFAULT_BOOT,
    seed: int = SEED,
    alpha: float = 0.05,
) -> dict[str, float]:
    """配对 bootstrap: Δ = metric(B) - metric(A), 在**同一批索引**上重采样。

    返回 delta / lo / hi / p_two_sided / 两模型各自的全样本指标。
    p 值定义: 2 * min(P(Δ<=0), P(Δ>=0)), 加 1 平滑避免 p 恰为 0。
    """
    y_a, sa, _ = _clean(y, score_a)
    y_b, sb, _ = _clean(y, score_b)
    if not np.array_equal(np.isfinite(y_a), np.isfinite(y_b)):
        # 两侧有效样本不同 -> 无法逐索引配对, 退回公共掩码
        pass
    n = min(len(y_a), len(y_b))
    y_a, sa, y_b, sb = y_a[:n], sa[:n], y_b[:n], sb[:n]
    if n < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "metric_a": float("nan"),
                "metric_b": float("nan"), "n": n}
    ma = float(metric_fn(y_a, sa))
    mb = float(metric_fn(y_b, sb))
    rs = np.random.RandomState(seed)
    diffs = []
    for _ in range(n_boot):
        idx = rs.randint(0, n, n)
        va = metric_fn(y_a[idx], sa[idx])
        vb = metric_fn(y_b[idx], sb[idx])
        if va == va and vb == vb:
            diffs.append(vb - va)
    if len(diffs) < 20:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "metric_a": ma, "metric_b": mb, "n": n}
    d = np.asarray(diffs, dtype=float)
    lo, hi = np.percentile(d, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    p_left = (np.sum(d <= 0) + 1) / (len(d) + 1)
    p_right = (np.sum(d >= 0) + 1) / (len(d) + 1)
    return {
        "delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
        "p_two_sided": float(min(1.0, 2 * min(p_left, p_right))),
        "metric_a": ma, "metric_b": mb, "n": n,
        "frac_b_gt_a": float(np.mean(d > 0)),
    }


# ===================================================================== 多重校正
def bh_fdr(pvals: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg FDR 校正 (自实现, 无 statsmodels 依赖)。

    返回与输入同序的 q 值; NaN 输入原样返回 NaN 且不占排名。
    """
    p = np.asarray(list(pvals), dtype=float)
    m = len(p)
    out = np.full(m, np.nan, dtype=float)
    ok = np.isfinite(p)
    m_ok = int(ok.sum())
    if m_ok == 0:
        return [float(x) for x in out]
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx], kind="mergesort")]
    q_sorted = p[order] * m_ok / np.arange(1, m_ok + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    out[order] = np.minimum(q_sorted, 1.0)
    return [float(x) for x in out]


def bonferroni(pvals: Sequence[float]) -> list[float]:
    p = np.asarray(list(pvals), dtype=float)
    m = int(np.isfinite(p).sum())
    with np.errstate(invalid="ignore"):
        return [float(min(1.0, x * m)) if x == x else float("nan") for x in p]


# ===================================================================== 自检
def _self_test() -> int:
    """内置自检: 完美排序 / 反向排序 / 随机 / 并列 / NaN 的行为必须符合预期。"""
    rs = np.random.RandomState(0)
    y = rs.normal(size=200)
    perfect = y.copy()
    reverse = -y.copy()
    rnd = rs.normal(size=200)
    ties = np.zeros(200)

    ok = True

    def check(name: str, cond: bool, detail: str = "") -> None:
        nonlocal ok
        ok = ok and cond
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")

    m_perf = compute_ranking_metrics(y, perfect)
    m_rev = compute_ranking_metrics(y, reverse)
    m_rnd = compute_ranking_metrics(y, rnd)
    check("完美排序 Spearman==1", abs(m_perf.values["spearman"] - 1.0) < 1e-9,
          f"{m_perf.values['spearman']:.6f}")
    check("反向排序 Spearman==-1", abs(m_rev.values["spearman"] + 1.0) < 1e-9,
          f"{m_rev.values['spearman']:.6f}")
    check("完美排序 NDCG@10==1", abs(m_perf.values["ndcg@10"] - 1.0) < 1e-9,
          f"{m_perf.values['ndcg@10']:.6f}")
    check("完美排序 EF@1% > 3", m_perf.values["ef@1%"] > 3,
          f"{m_perf.values['ef@1%']:.3f}")
    # n=200 时 EF@1% 的 k=2、基准率 10% -> 期望命中 0.2 条, 量化为 0 或 1 是正常的;
    # 因此这里对**更平滑的** EF@10% (k=20) 断言, 并对 EF@1% 只断言"接近随机量级"。
    check("随机 EF@10% 在 0.2~3", 0.2 < m_rnd.values["ef@10%"] < 3.0,
          f"{m_rnd.values['ef@10%']:.3f}")
    check("随机 EF@1% 有限且非负", np.isfinite(m_rnd.values["ef@1%"])
          and m_rnd.values["ef@1%"] >= 0,
          f"{m_rnd.values['ef@1%']:.3f} (k={m_rnd.k_effective['ef@1%']}, 离散)")
    check("反向 NDCG@10 < 0.6", m_rev.values["ndcg@10"] < 0.6,
          f"{m_rev.values['ndcg@10']:.3f}")

    # 并列值: 全并列 -> 相关性 NaN; 期望口径下富集必须收敛到 1 (不能凭空 2 倍)
    m_tie = compute_ranking_metrics(y, ties)
    check("全并列 Spearman 为 NaN", m_tie.values["spearman"] != m_tie.values["spearman"])
    check("全并列 EF@10% == 1 (期望口径)", abs(m_tie.values["ef@10%"] - 1.0) < 1e-9,
          f"{m_tie.values['ef@10%']:.6f}")
    check("全并列 NDCG@10 < 1", m_tie.values["ndcg@10"] < 0.999,
          f"{m_tie.values['ndcg@10']:.6f}")
    # 随机打乱行序不应改变期望口径下的结果 (无行序伪影)
    perm = np.random.RandomState(7).permutation(len(y))
    m_tie_p = compute_ranking_metrics(y[perm], ties[perm])
    check("全并列 EF@10% 与行序无关",
          abs(m_tie_p.values["ef@10%"] - m_tie.values["ef@10%"]) < 1e-9,
          f"{m_tie_p.values['ef@10%']:.6f}")

    # NaN 处理
    y_nan = y.copy()
    y_nan[:10] = np.nan
    m_nan = compute_ranking_metrics(y_nan, perfect)
    check("NaN 行被整行丢弃", m_nan.n == 190 and m_nan.n_dropped == 10,
          f"n={m_nan.n}, dropped={m_nan.n_dropped}")

    # 极小样本不抛错
    try:
        m_small = compute_ranking_metrics(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
        check("n=2 不抛错", True, f"spearman={m_small.values['spearman']}")
    except Exception as exc:  # noqa: BLE001
        check("n=2 不抛错", False, repr(exc))

    # 空输入
    try:
        m_empty = compute_ranking_metrics(np.array([]), np.array([]))
        check("空输入不抛错", m_empty.n == 0)
    except Exception as exc:  # noqa: BLE001
        check("空输入不抛错", False, repr(exc))

    # 配对 bootstrap: A=随机, B=完美 -> 显著
    pb = paired_bootstrap_test(spearman_rho, y, rnd, perfect, n_boot=500, seed=1)
    check("配对检验 B(完美)>A(随机) 显著", pb["p_two_sided"] < 0.05 and pb["delta"] > 0,
          f"delta={pb['delta']:.3f}, p={pb['p_two_sided']:.4f}")
    # 自己 vs 自己 -> 不显著
    pb0 = paired_bootstrap_test(spearman_rho, y, rnd, rnd, n_boot=500, seed=1)
    check("配对检验 自己 vs 自己 不显著",
          pb0["p_two_sided"] > 0.5 and abs(pb0["delta"]) < 1e-9,
          f"delta={pb0['delta']:.6f}, p={pb0['p_two_sided']:.4f}")

    # BH-FDR 单调性与边界
    ps = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    qs = bh_fdr(ps)
    check("BH q >= p", all(q >= p - 1e-12 for q, p in zip(qs, ps)))
    check("BH 最小 p 的 q <= 0.05", qs[0] <= 0.05, f"q0={qs[0]:.5f}")
    check("BH 保序", all(qs[i] <= qs[i + 1] + 1e-12 for i in range(len(qs) - 1)))
    # 单个假设时 BH 的 q 恰等于 p (不是"必须不显著"); 真正该测的是:
    # 当存在 m 个假设时, 同一个 p 的 q 会被放大到 p*m/rank。
    check("BH 单假设 q == p", abs(bh_fdr([0.04])[0] - 0.04) < 1e-12,
          f"{bh_fdr([0.04])[0]:.4f}")
    q_multi = bh_fdr([0.04, 0.5, 0.6, 0.7])
    check("BH 多假设放大 p=0.04 -> q=0.16", abs(q_multi[0] - 0.16) < 1e-12,
          f"{q_multi[0]:.4f}")
    check("BH 保持最大 q <= 1", max(q_multi) <= 1.0)
    check("Bonferroni p=0.04,m=4 -> 0.16", abs(bonferroni([0.04, 0.5, 0.6, 0.7])[0] - 0.16) < 1e-12)

    # 方向陷阱自证: 传 ΔG (越负越强) 而不是 score -> 富集应接近/低于随机
    dg = -y
    m_wrong = compute_ranking_metrics(dg, rnd)
    check("方向陷阱检测 (ΔG 直传时低效)", m_wrong.values["spearman"] < 0.2,
          f"{m_wrong.values['spearman']:.4f}")

    print("")
    print("自检结果:", "全部通过" if ok else "存在失败")
    return 0 if ok else 1


# ===================================================================== CLI
def main() -> int:
    ap = argparse.ArgumentParser(description="排序/富集指标库 (含自检)")
    ap.add_argument("--self-test", action="store_true", help="运行内置自检")
    ap.add_argument("--demo", action="store_true", help="在合成数据上演示指标输出")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-boot", type=int, default=DEFAULT_BOOT)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%H:%M:%S")

    if args.self_test or not args.demo:
        print("=" * 78)
        print("metrics_ranking 内置自检")
        print("=" * 78)
        rc = _self_test()
        if not args.demo:
            return rc

    if args.demo:
        t0 = time.time()
        rs = np.random.RandomState(args.seed)
        n = 400
        y = -8.0 - 3.0 * rs.normal(size=n)      # ΔG, 越负越强
        good = -y + rs.normal(scale=1.6, size=n)   # 好模型 (ΔG 空间)
        bad = -y * 0.2 + rs.normal(scale=3.0, size=n)  # 差模型
        print("=" * 78)
        print("合成演示: y = ΔG (越负越强), 输入统一转换为 score = -ΔG")
        print("=" * 78)
        for name, pred_dg in (("good", good), ("bad", bad)):
            m = compute_ranking_metrics(negate_dG(y), negate_dG(pred_dg))
            print(f"\n--- {name} ---")
            for line in m.summary_lines():
                print(line)
        pb = paired_bootstrap_test(
            spearman_rho, negate_dG(y), negate_dG(bad), negate_dG(good),
            n_boot=args.n_boot, seed=args.seed)
        print(f"\n配对检验 (good - bad) Spearman: delta={pb['delta']:+.4f} "
              f"CI=[{pb['lo']:+.4f},{pb['hi']:+.4f}] p={pb['p_two_sided']:.4f}")
        print(f"耗时 {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
