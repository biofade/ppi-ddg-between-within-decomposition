# -*- coding: utf-8 -*-
"""scripts/v2.5/05_conservation_features.py -- V2.5 任务 3.1: 进化保守性 (ESM-2 代理).

**先说明限制（诚实优先）**
--------------------------
任务书首选方案是 **MSA + Shannon 熵**（UniRef50）。但:
* 本机**没有 MSA**（未下载 UniRef，也未跑 jackhmmer/MMseqs2）；
* 下载并比对 7614 条序列的 MSA 远超本机 CPU 与磁盘预算（>10 GB 且需数小时）。

因此本脚本实现**任务书给出的备选方案**：用 **ESM-2 表示**近似保守性。
具体用两个可用的代理量（都**不使用标签**）：

1. **表示偏离度** `dev_i = ||h_i − mean_j h_j||₂ / sqrt(D)`
   * 依据: 功能/保守位点的嵌入往往更"典型"，而异常位点偏离均值；
   * 这是**近似**，不是严格的进化熵；已在报告标注。
2. **池化前后的重构误差**: `||h_i||₂` 的分布 (作为第二代理量，互相校验)。

汇总到**界面残基**层面: 均值 / 最大 / 标准差 / 高偏离残基比例。

> ⚠️ **本文件不产出"Shannon 熵"**。任何把它称作熵的说法都是不诚实的。
> 报告中一律称为 **conservation proxy (表示偏离度)**。

输出: `features/v2.5_conservation_features.csv`、`features/v2.5_conservation_report.md`
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v2.5_conservation_features.csv"
OUT_MD = C.FEAT / "v2.5_conservation_report.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def _per_residue_dev(data, side: str, chunk: int = 512):
    """分块生成每条链的逐残基偏离度与掩码 (生成器, 内存友好)。"""
    X = data[f"{side}_topk"]; M = data[f"{side}_mask"]
    n, K, D = X.shape
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        xb = np.asarray(X[s:e], dtype=np.float32)
        mb = M[s:e]
        mu = (xb * mb[..., None]).sum(1, keepdims=True) / np.maximum(
            mb.sum(1, keepdims=True)[..., None], 1.0)
        dev = np.linalg.norm(xb - mu, axis=2) / np.sqrt(D)
        yield s, e, np.where(mb, dev, np.nan), mb
        del xb, mb


def _threshold(data, side: str, tr: np.ndarray, q: float) -> float:
    """训练集内所有有效残基偏离度的分位数 (分块累积, 不整块物化)。"""
    vals = []
    for s, e, dev, mb in _per_residue_dev(data, side):
        sel = tr[s:e]
        if sel.any():
            sub = dev[sel]
            msub = mb[sel]
            vals.append(sub[msub])
    allv = np.concatenate(vals) if vals else np.zeros(0)
    return float(np.nanpercentile(allv, q * 100)) if len(allv) else float("nan")


def _high_frac(data, side: str, thr: float) -> np.ndarray:
    """每条链中偏离度超过阈值的高偏离残基比例。"""
    out = []
    for s, e, dev, mb in _per_residue_dev(data, side):
        out.append(np.nansum((dev > thr) & mb, axis=1)
                   / np.maximum(mb.sum(1), 1))
    return np.concatenate(out) if out else np.zeros(0)


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 保守性代理特征")
    ap.add_argument("--high-frac", type=float, default=0.8,
                    help="'高偏离残基'的阈值分位 (默认 top 20%)")
    args = ap.parse_args()

    t0 = time.time()
    emit("# V2.5 任务 3.1 — 进化保守性（ESM-2 表示代理）")
    emit("")
    emit("> ⚠️ **没有 MSA，因此没有 Shannon 熵**。这里用的是**表示偏离度代理**。")
    emit("")

    data = C.load_aligned(use_perres_pooled=True)
    n = len(data["y"])
    emit(f"## 1. 数据")
    emit("")
    emit(f"- 行数 {n}; per-residue 张量 (n, K, 640) fp16")

    # 阈值从**训练集**的偏离度分布取 (不涉及标签, 但保持与划分一致的习惯)
    tr = data["tr"]

    out_rows: list[dict] = []
    for side in ("lig", "rec"):
        # **分块计算** (整块 float32 化会 OOM: 6885x32x640 float64 = 1.05 GiB)
        M = data[f"{side}_mask"]
        st = C._chunked_residue_reduce(data, side, "stats")
        dev = np.where(M, st["dev_mean"][:, None], np.nan)  # 占位, 下面重算阈值
        pre = f"cv_{side}"
        # 阈值: 训练集内所有有效残基偏离度的 80 分位
        # (需要逐残基偏离度 -> 再分块取一次分布)
        thr = _threshold(data, side, data["tr"], args.high_frac)
        high = _high_frac(data, side, thr)
        out_rows.append(pd.DataFrame({
            "pdb_ids": data["pdb_ids"],
            f"{pre}_dev_mean": st["dev_mean"], f"{pre}_dev_max": st["dev_max"],
            f"{pre}_dev_std": st["dev_std"],
            f"{pre}_high_dev_frac": high,
            f"{pre}_norm_mean": st["norm_mean"], f"{pre}_norm_std": st["norm_std"],
            f"{pre}_n_res": st["n_res"],
        }))
        emit(f"- {side}: 偏离度阈值 (训练集 {args.high_frac*100:.0f} 分位) = {thr:.4f}; "
             f"均值 {np.nanmean(st['dev_mean']):.4f}")

    df = pd.concat(out_rows, axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")

    cols = [c for c in df.columns if c.startswith("cv_")]
    emit("")
    emit(f"## 2. 特征 ({len(cols)} 个)")
    emit("")
    emit("| 特征 | 均值 | 标准差 | 范围 |")
    emit("|---|---:|---:|---|")
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce")
        emit(f"| `{c}` | {v.mean():.4f} | {v.std():.4f} | "
             f"[{np.nanmin(v):.4f}, {np.nanmax(v):.4f}] |")

    emit("")
    emit("## 3. 方法说明与局限")
    emit("")
    emit("1. **这不是 Shannon 熵**。Shannon 熵需要 MSA（同源序列比对），")
    emit("   本项目**没有 MSA**，也未下载 UniRef50（磁盘/CPU 预算不允许: "
         "7614 条序列的比对 >10 GB 且需数小时）。")
    emit("2. 本实现用 **ESM-2 逐残基表示的偏离度** `||h_i − mean_j h_j||₂/√D` 作为")
    emit("   **保守性代理**。动机: 语言模型的表示对功能位点通常更敏感；")
    emit("   但**偏离度大 ≠ 保守性高**，两者只是相关，方向甚至可能是反的。")
    emit("3. 因此这组特征的**可解释性弱**，在消融中应视为「又一个嵌入派生量」，")
    emit("   而不是「进化信息」。报告中一律称其为 conservation **proxy**。")
    emit("4. 阈值取自训练集分布（不使用标签，但保证了划分一致性）。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV.name} | {OUT_MD.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
