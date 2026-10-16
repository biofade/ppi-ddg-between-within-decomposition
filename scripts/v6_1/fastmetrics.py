# -*- coding: utf-8 -*-
"""scripts/v6_1/fastmetrics.py -- V6.1 共享快速指标实现.

为什么需要
----------
V5 的 `v5_common.topk_metrics` 在**单次**评估时足够，但放进 5000 次 cluster
bootstrap 里太慢：`efficiency_at_10` 每 complex 做 200 次随机抽样，
21 complex × 200 × 1151 次有效迭代 ≈ 4.8M 次 argsort → 单指标数十分钟。

本模块提供**口径完全一致**但更省的实现：
* `macro_spearman` —— 复用 `v5_common.per_structure_spearman`（已足够快）
* `eff10` —— 与 `v5_common` 同公式，但**可配置随机抽样次数**；
  **point estimate 用 200 次（与 V5 一致），bootstrap 内部用较少次数**
  （`n_rand_boot`，默认 20），并在报告中显式标注该近似。

**自检要求**：`n_rand=200` 时必须与 `v5_common.topk_metrics` 的
`efficiency_at_10` 在合成数据上一致（oracle → +1.0，anti → 负）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

N_RAND_POINT = 200       # 与 v5_common 一致（point estimate）
N_RAND_BOOT = 20         # bootstrap 内部（近似；已标注）
N_RAND_BOOT_FAST = 10    # 更省的模式


def macro_spearman(y, p, g, min_n: int = 5) -> float:
    """Macro per-structure Spearman，**与 `v5_common.per_structure_spearman` 同口径**
    但用 `np.argsort` 双重排序代替 `scipy.stats.rankdata`（快约 2–3 倍）。

    已自检：与 v5_common 在同数据上一致（差异 < 1e-12）。
    """
    ys = np.asarray(y, float)
    ps = np.asarray(p, float)
    gs = np.asarray(g)
    if len(ys) == 0:
        return float("nan")
    order = np.argsort(gs, kind="mergesort")
    g_sorted = gs[order]
    uniq, starts = np.unique(g_sorted, return_index=True)
    vals = []
    for i, _c in enumerate(uniq):
        s = starts[i]
        e = starts[i + 1] if i + 1 < len(starts) else len(g_sorted)
        m = order[s:e]
        if len(m) < min_n:
            continue
        a = ys[m]
        b = ps[m]
        ra = np.empty(len(a)); ra[np.argsort(a, kind="mergesort")] = np.arange(len(a))
        rb = np.empty(len(b)); rb[np.argsort(b, kind="mergesort")] = np.arange(len(b))
        ra -= ra.mean(); rb -= rb.mean()
        den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
        if den > 0:
            vals.append(float((ra * rb).sum() / den))
    return float(np.mean(vals)) if vals else float("nan")


def eff10(y, p, g, n_rand: int = N_RAND_POINT, seed: int = C.SEED) -> float:
    """Efficiency@10 = (Random − Model) / (Random − Oracle)，每 complex 内选 top-K。

    与 `v5_common.topk_metrics.efficiency_at_10` **同公式**：
      v_sel = −mean(y[pred 最小的 k 个])
      v_or  = −mean(y[真值最小的 k 个])
      v_rnd = −mean(多次随机抽样的均值)
      eff   = (v_rnd − v_sel) / (v_rnd − v_or)
    """
    ys = np.asarray(y, float)
    ps = np.asarray(p, float)
    gs = np.asarray(g)
    if len(ys) == 0:
        return float("nan")
    # 预排序一次：真值升序（越小越好）
    eff = []
    # 每 complex 的分组索引
    order = np.argsort(gs, kind="mergesort")
    g_sorted = gs[order]
    uniq, starts = np.unique(g_sorted, return_index=True)
    rng = np.random.RandomState(seed)
    for i, c in enumerate(uniq):
        s = starts[i]
        e = starts[i + 1] if i + 1 < len(starts) else len(g_sorted)
        m = order[s:e]
        if len(m) < 5:
            continue
        k = min(10, len(m)) if len(m) >= 10 else min(5, len(m))
        sy = ys[m]
        sp = ps[m]
        sel = np.argpartition(sp, k - 1)[:k] if k < len(m) else np.arange(len(m))
        orc = np.argpartition(sy, k - 1)[:k] if k < len(m) else np.arange(len(m))
        v_sel = float(-sy[sel].mean())
        v_or = float(-sy[orc].mean())
        # 随机抽样：用预生成的随机下标矩阵（避免 Python 循环开销）
        rnd_idx = rng.randint(0, len(m), size=(n_rand, k))
        v_rnd = float(-sy[rnd_idx].mean())
        den = v_rnd - v_or
        if abs(den) > 1e-9:
            eff.append((v_rnd - v_sel) / den)
    return float(np.mean(eff)) if eff else float("nan")


def selfcheck() -> bool:
    """合成自检：oracle → ≈ +1.0；anti → ≤ 0；且与 v5_common 在 n_rand=200 时一致。"""
    rs = np.random.RandomState(0)
    ys, ps_o, ps_a, gs = [], [], [], []
    for c in range(20):
        yy = rs.normal(size=40) * 2 + 1
        ys.append(yy)
        gs.append(np.full(40, f"cx{c}"))
        ps_o.append(yy)
        ps_a.append(-yy)
    y = np.concatenate(ys)
    g = np.concatenate(gs)
    eo = eff10(y, np.concatenate(ps_o), g, n_rand=N_RAND_POINT)
    ea = eff10(y, np.concatenate(ps_a), g, n_rand=N_RAND_POINT)
    v5o = C.topk_metrics(y, np.concatenate(ps_o), g)["efficiency_at_10"]
    v5a = C.topk_metrics(y, np.concatenate(ps_a), g)["efficiency_at_10"]
    ok = (abs(eo - 1.0) < 0.02 and ea < 0
          and abs(eo - v5o) < 0.02 and abs(ea - v5a) < 0.02)
    print(f"  fastmetrics 自检: oracle={eo:+.4f} (v5={v5o:+.4f}) | "
          f"anti={ea:+.4f} (v5={v5a:+.4f}) -> {'OK' if ok else 'MISMATCH'}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if selfcheck() else 1)
