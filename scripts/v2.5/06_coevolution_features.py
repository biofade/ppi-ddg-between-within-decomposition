# -*- coding: utf-8 -*-
"""scripts/v2.5/06_coevolution_features.py -- V2.5 任务 3.2: 共进化耦合 (简化).

**先说明什么做不到**
--------------------
严格共进化分析（DCA / GREMLIN / EVcouplings）需要 **MSA**。
本项目没有 MSA（见 `05_conservation_features.py` 的说明）。

另一个候选是 **ESM-2 注意力矩阵** —— 但那需要重新跑一次前向并保存 attention。
本次提取的 per-residue 缓存**只存了 layer-30 的表示，没有存 attention**
（见 `logs/v2.5/perres_schema.log` 的字段清单）。**重新提取需要再次使用 Kaggle GPU。**

因此本脚本实现**可用的简化代理**，并明确标注其地位：

1. **嵌入空间内的残基间耦合**: 对每个条目，计算界面残基对的表示余弦相似度
   * `coupling_ij = cos(h_i, h_j)`
   * 汇总: 均值 / 最大 / 标准差 / 高耦合对数
   * ⚠️ 这是**表示相似度**，不是共进化信号（共进化指跨序列的协同突变）。
2. **跨链耦合**: 只在配体-受体残基对之间算（同一链内的相似度通常更高，
   因为它反映的是序列局部性，而非耦合）。

> ⚠️ 报告中一律称为 **coupling proxy（表示相似度）**，**不**称"共进化"。

输出: `features/v2.5_coevolution_features.csv`、`features/v2.5_coevolution_report.md`
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v2.5_coevolution_features.csv"
OUT_MD = C.FEAT / "v2.5_coevolution_report.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 耦合代理特征")
    ap.add_argument("--high-cos", type=float, default=0.5,
                    help="'高耦合对'的余弦阈值")
    args = ap.parse_args()

    t0 = time.time()
    emit("# V2.5 任务 3.2 — 耦合代理（表示相似度）")
    emit("")
    emit("> ⚠️ **没有 MSA -> 无法做真正的共进化 (DCA)**;")
    emit("> ⚠️ **per-residue 缓存没有存 attention** -> 也无法用 ESM-2 注意力。")
    emit("> 因此这里是**表示相似度代理**, 报告中不称其为共进化。")
    emit("")

    data = C.load_aligned(use_perres_pooled=True)
    n = len(data["y"])
    emit("## 1. 数据")
    emit("")
    emit(f"- 行数 {n}")

    # **分块计算**: 整块 float32 化会 OOM (6885x32x640 float32 = 538 MiB,
    # 而本机空闲内存常在 1-2 GB, 且同时还有其它数组)。
    Xl = data["lig_topk"]; Xr = data["rec_topk"]
    Ml = data["lig_mask"]; Mr = data["rec_mask"]
    n = Xl.shape[0]
    chunk = 256

    feats = np.full((n, 6), np.nan, dtype=np.float32)
    t1 = time.time()
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        gl_all = np.asarray(Xl[s:e], dtype=np.float32)
        gr_all = np.asarray(Xr[s:e], dtype=np.float32)
        ml = Ml[s:e]; mr = Mr[s:e]
        # 逐行 L2 归一化 (就地, 不再复制整块)
        gl_all /= np.maximum(np.linalg.norm(gl_all, axis=2, keepdims=True), 1e-8)
        gr_all /= np.maximum(np.linalg.norm(gr_all, axis=2, keepdims=True), 1e-8)
        for j in range(e - s):
            gl = gl_all[j][ml[j]]
            gr = gr_all[j][mr[j]]
            if len(gl) == 0 or len(gr) == 0:
                continue
            S = gl @ gr.T                      # (kl, kr) 余弦
            flat = S.ravel()
            feats[s + j, 0] = flat.mean()
            feats[s + j, 1] = flat.max()
            feats[s + j, 2] = flat.std()
            feats[s + j, 3] = float((flat > args.high_cos).mean())
            feats[s + j, 4] = S.max(axis=1).mean()
            feats[s + j, 5] = S.max(axis=0).mean()
            del S, flat, gl, gr
        del gl_all, gr_all
        if (s // chunk) % 5 == 0:
            emit(f"  {e}/{n} ({time.time()-t1:.0f}s)")

    df = pd.DataFrame({
        "pdb_ids": data["pdb_ids"],
        "co_cross_cos_mean": feats[:, 0],
        "co_cross_cos_max": feats[:, 1],
        "co_cross_cos_std": feats[:, 2],
        "co_cross_high_frac": feats[:, 3],
        "co_cross_ligmax_mean": feats[:, 4],
        "co_cross_recmax_mean": feats[:, 5],
    })
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")
    cols = [c for c in df.columns if c.startswith("co_")]
    emit("")
    emit(f"## 2. 特征 ({len(cols)} 个)")
    emit("")
    emit("| 特征 | 均值 | 标准差 | 范围 |")
    emit("|---|---:|---:|---|")
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce")
        emit(f"| `{c}` | {np.nanmean(v):.4f} | {np.nanstd(v):.4f} | "
             f"[{np.nanmin(v):.4f}, {np.nanmax(v):.4f}] |")

    emit("")
    emit("## 3. 局限（重要）")
    emit("")
    emit("1. **这不是共进化**。真正的共进化（协同突变）需要跨物种同源序列的 MSA。")
    emit("   这里算的是**同一结构内跨链残基表示的余弦相似度**。")
    emit("2. 相似度高**不等于**耦合强: 两个残基表示相似可能只是因为它们都是疏水/"
         "都带正电, 与协同进化无关。因此这组特征的信息量**很可能很低**, ")
    emit("   预期在消融中贡献很小 —— 这本身是一个**可检验的预测**, 见任务 5。")
    emit("3. ESM-2 **attention 未保存**: 本次 per-residue 提取只输出了 layer-30 表示。")
    emit("   若要真正的注意力耦合, 需要重跑 notebook 并保存 attention（需 GPU）。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV.name} | {OUT_MD.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
