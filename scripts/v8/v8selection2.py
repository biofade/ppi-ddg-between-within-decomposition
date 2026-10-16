# -*- coding: utf-8 -*-
"""scripts/v8/v8selection2.py -- AXIS C v2：**比例型**选择效用（修复区分度地板）。

## 为什么需要 v2

`Efficiency@K = Σ gain(selected) / Σ gain(ideal top-K)` 中
`gain = 1..n` 是**名次线性**的，导致：
  * n≈20–30 时，任何还算合理的模型都能抓到前 K 中一大半 ⇒ 比值被压在高位
  * 实测 `U_CANONICAL` 上 **IQR 仅 0.0143**（p10=0.8036, p90=0.8306）
  * ⇒ AXIS C 无法区分模型，H3 无法被有意义地检验

## v2 定义（全部为**集合重叠比例**，对 n 不敏感）

方向约定：`ΔΔG` 越负越有利。真值有利序 = y 升序；预测有利序 = p 升序。

| 指标 | 定义 | 直觉 |
|---|---|---|
| `util_topq` | 真值前 q% 的集合 ∩ 预测前 q% 的集合，除以 \|真值前 q%\| | 「我最想做的 q% 里，有多少被模型排进它自己的前 q%」 |
| `util_topm` | 同上，但用**固定条数** m（而非比例），m = min(10, max(3, n//3)) | 短名单场景 |
| `util_pw` | 同靶内**有序对**一致比例（含并列 0.5） | 全序吻合度，与 Spearman 同族但直接可读 |
| `precision_at_true_topq` | 预测前 q% 中有多少落在真值前 q% | 精确率视角 |

判据：若 v2 的 IQR / 极差显著大于 `Efficiency@10`，则 AXIS C 的区分度问题被修复。
"""

from __future__ import annotations

import numpy as np


def util_pw(y, p):
    """同靶内有符号对一致比例（并列计 0.5）。"""
    n = len(y)
    if n < 2:
        return np.nan
    iu = np.triu_indices(n, 1)
    dy = y[iu[0]] - y[iu[1]]
    dp = p[iu[0]] - p[iu[1]]
    keep = dy != 0
    if not keep.any():
        return np.nan
    dyv, dpv = dy[keep], dp[keep]
    agree = (np.sign(dyv) == np.sign(dpv)) & (dpv != 0)
    return float((agree.sum() + 0.5 * (dpv == 0).sum()) / dyv.size)


def util_overlap(y, p, k):
    """真值前 k（有利）∩ 预测前 k（有利），除以 k。k 会被裁剪到 n。"""
    n = len(y)
    k = max(1, min(k, n))
    ar = np.arange(n)
    true_top = np.argsort(np.lexsort((ar, y)))[:k]      # y 升序 = 最有利
    pred_top = np.argsort(np.lexsort((ar, p)))[:k]
    return float(len(set(true_top.tolist()) & set(pred_top.tolist())) / k)


def util_precision_true_topq(y, p, q=0.1):
    """预测前 q% 中落在真值前 q% 的比例；随机基线 = q%。"""
    n = len(y)
    k = max(1, int(round(n * q)))
    ar = np.arange(n)
    true_top = set(np.argsort(np.lexsort((ar, y)))[:k].tolist())
    pred_top = np.argsort(np.lexsort((ar, p)))[:max(1, k)]
    inter = len(set(pred_top.tolist()) & true_top)
    return float(inter / max(1, len(pred_top)))


def util_overlap_norm(y, p, k):
    """把 util_overlap 归一化到「随机基线 0、完美 1」。"""
    n = len(y)
    k = max(1, min(k, n))
    base = k / n
    o = util_overlap(y, p, k)
    return float((o - base) / (1.0 - base)) if base < 1 else np.nan


def axis_selection_v2(U, qs=(0.05, 0.10, 0.20, 0.30)):
    """对每个 eligible target 求 v2 效用，再取宏平均。"""
    idx = U["idx"]; eligible = U["eligible"]
    y, p = U["y"], U["p"]
    acc = {f"util_topq@{int(q*100)}": [] for q in qs}
    acc.update({f"util_topq_norm@{int(q*100)}": [] for q in qs})
    acc.update({f"util_prec_topq@{int(q*100)}": [] for q in qs})
    acc["util_topm"] = []
    acc["util_topm_norm"] = []
    acc["util_pw"] = []
    for t in eligible:
        i = idx[t]
        ty, tp = y[i], p[i]
        n = len(i)
        for q in qs:
            k = max(1, int(round(n * q)))
            acc[f"util_topq@{int(q*100)}"].append(util_overlap(ty, tp, k))
            acc[f"util_topq_norm@{int(q*100)}"].append(util_overlap_norm(ty, tp, k))
            acc[f"util_prec_topq@{int(q*100)}"].append(
                util_precision_true_topq(ty, tp, q))
        m = min(10, max(3, n // 3))
        acc["util_topm"].append(util_overlap(ty, tp, m))
        acc["util_topm_norm"].append(util_overlap_norm(ty, tp, m))
        acc["util_pw"].append(util_pw(ty, tp))

    def _nm(v):
        v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
        return float(v.mean()) if v.size else np.nan
    return {k: _nm(v) for k, v in acc.items()}


V2_COLS = (["util_topm", "util_topm_norm", "util_pw"]
           + [f"util_topq@{q}" for q in (5, 10, 20, 30)]
           + [f"util_topq_norm@{q}" for q in (5, 10, 20, 30)]
           + [f"util_prec_topq@{q}" for q in (5, 10, 20, 30)])
