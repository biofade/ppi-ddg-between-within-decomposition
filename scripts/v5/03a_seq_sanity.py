# -*- coding: utf-8 -*-
"""scripts/v5/03a_seq_sanity.py -- V5.0: WT/mutant 序列一致性前置检查（只读诊断）.

为什么先做这一步
----------------
GATE B 要基于「pair 内 WT row 序列 vs mutant row 序列」做 global alignment。
但如果 WT 与 mutant 行**根本不是同一 construct**（序列长度差异远大于突变数），
那么 alignment 结果没有意义的 —— 必须先量化这一点，否则整个 MIPA 建立
在错位映射上（V5 §23「不允许静默错位」）。

本脚本量化（**不做任何映射决定**，只报事实）:
1. 每 pair 两侧序列长度差的分布；
2. 用 Biopython global alignment 算两侧 identity；
3. 差异位点数 vs annotation 声明的突变数是否自洽；
4. 按「长度差 ≤ 2×突变数」与「identity 高」分桶，估计**可映射**的 pair 比例。

输出: experiments/v5_0/manifests/seq_sanity.md
"""

from __future__ import annotations

import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT = C.V5 / "manifests" / "seq_sanity.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    t0 = time.time()
    emit("# V5.0 — WT/mutant 序列一致性前置检查（GATE B 之前）")
    emit("")
    emit("> 目的: 确认 pair 内 WT row 与 mutant row 是否表示**同一 construct**。")
    emit("> 若序列长度差异远大于突变数，或 identity 很低，则 alignment 映射无效。")
    emit("")

    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score = 2.0
    al.mismatch_score = -1.0
    al.open_gap_score = -5.0
    al.extend_gap_score = -1.0

    raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    fwd = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = fwd[fwd["direction"] == "forward"].reset_index(drop=True)
    emit(f"- forward pair: **{len(fwd)}**")
    emit("")

    # ---------- 1. 同样 pdb 内 WT 行是否唯一 ----------
    emit("## 1. 每个 pdb 的 WT 行数（若 >1，说明 WT 定义本身有歧义）")
    emit("")
    wtc = fwd.groupby("pdb_id")["wt_row"].nunique()
    emit(f"- 唯一 pdb: {len(wtc)}；其中 WT 行数 >1 的 pdb: **{int((wtc > 1).sum())}**")
    emit(f"- 每 pdb 的 WT 行数分布: {wtc.value_counts().sort_index().head(8).to_dict()}")
    emit("")

    # ---------- 2. 长度差与 identity ----------
    emit("## 2. 两侧序列长度差与 global-alignment identity")
    emit("")
    recs = []
    n = len(fwd)
    for i, r in fwd.iterrows():
        try:
            wi, mi = int(r["wt_row"]), int(r["mut_row"])
        except Exception:  # noqa: BLE001
            continue
        row = {"pdb_id": r["pdb_id"], "n_ann": int(r["n_mut"]) if pd.notna(r["n_mut"]) else 0}
        for side, col in (("lig", "Ligand Sequences"), ("rec", "Receptor Sequences")):
            w = str(raw.at[wi, col]); m = str(raw.at[mi, col])
            row[f"{side}_len_wt"] = len(w)
            row[f"{side}_len_mut"] = len(m)
            row[f"{side}_dlen"] = len(m) - len(w)
            row[f"{side}_same"] = int(w == m)
            if w == m:
                row[f"{side}_ident"] = 1.0
                row[f"{side}_ndiff"] = 0
            elif len(w) > 2200 or len(m) > 2200:
                row[f"{side}_ident"] = np.nan
                row[f"{side}_ndiff"] = np.nan
            else:
                try:
                    a = al.align(w, m)[0]
                    aw, am = str(a[0]), str(a[1])
                    nd = sum(1 for x, y in zip(aw, am) if x != "-" and y != "-" and x != y)
                    Lw = len(w)
                    row[f"{side}_ident"] = (Lw - nd) / max(1, Lw)
                    row[f"{side}_ndiff"] = nd
                except Exception:  # noqa: BLE001
                    row[f"{side}_ident"] = np.nan
                    row[f"{side}_ndiff"] = np.nan
        recs.append(row)
        if (i + 1) % 1500 == 0:
            emit(f"  ... {i + 1}/{n}")
    d = pd.DataFrame(recs)
    emit("")
    emit("| 侧 | \\|Δlen\\| 中位 | \\|Δlen\\| p95 | \\|Δlen\\| 最大 | 序列相同比例 | identity 中位 | identity p05 |")
    emit("|---|---:|---:|---:|---:|---:|---:|")
    for side in ("lig", "rec"):
        a = d[f"{side}_dlen"].abs()
        idt = d[f"{side}_ident"].dropna()
        emit(f"| {side} | {a.median():.0f} | {a.quantile(0.95):.0f} | {a.max():.0f} | "
             f"{d[f'{side}_same'].mean() * 100:.1f}% | {idt.median():.4f} | "
             f"{idt.quantile(0.05):.4f} |")
    emit("")

    # ---------- 3. 差异位点数 vs annotation 突变数 ----------
    emit("## 3. 差异位点数 vs annotation 声明的突变数（自洽性）")
    emit("")
    emit("> 若 alignment 差异数远大于 annotation 声明的突变数，说明两侧不是同一 construct。")
    emit("")
    tot_ndiff = d["lig_ndiff"].fillna(0) + d["rec_ndiff"].fillna(0)
    ratio = tot_ndiff / d["n_ann"].clip(lower=1)
    emit(f"- annotation 突变数: 中位 {d['n_ann'].median():.0f}, 均值 {d['n_ann'].mean():.2f}, 最大 {d['n_ann'].max()}")
    emit(f"- alignment 总差异位点数: 中位 {tot_ndiff.median():.0f}, 均值 {tot_ndiff.mean():.1f}, 最大 {tot_ndiff.max():.0f}")
    emit(f"- 差异数 / 突变数: 中位 **{ratio.median():.2f}**, p95 {ratio.quantile(0.95):.2f}")
    emit("")
    self_consistent = (ratio <= 2.0) & (tot_ndiff >= d["n_ann"])
    emit(f"- **自洽（差异数 ≤ 2×突变数）的 pair: {int(self_consistent.sum())} / {len(d)} = "
         f"{self_consistent.mean() * 100:.2f}%**")
    emit("")

    # ---------- 4. 按 identity 分桶估计可映射比例 ----------
    emit("## 4. 可映射性估计（identity 分桶，取两侧较低者）")
    emit("")
    lo = d[["lig_ident", "rec_ident"]].min(axis=1)
    bins = [(0.99, 1.001, "≥0.99"), (0.95, 0.99, "0.95–0.99"),
            (0.90, 0.95, "0.90–0.95"), (0.70, 0.90, "0.70–0.90"),
            (0.0, 0.70, "<0.70")]
    emit("| 较低侧 identity | pair 数 | 占比 |")
    emit("|---|---:|---:|")
    for lo_, hi_, name in bins:
        m = (lo >= lo_) & (lo < hi_)
        emit(f"| {name} | {int(m.sum())} | {m.mean() * 100:.2f}% |")
    emit(f"| NaN（过长跳过） | {int(lo.isna().sum())} | {lo.isna().mean() * 100:.2f}% |")
    emit("")
    emit("## 5. 结论（事实陈述，不含映射决定）")
    emit("")
    emit(f"- 可映射候选（identity ≥0.95 且自洽）= "
         f"**{int(((lo >= 0.95) & self_consistent).sum())} / {len(d)}**")
    emit(f"- 完全不可映射（较低侧 identity <0.70）= **{int((lo < 0.70).sum())}**")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    d.to_csv(C.V5 / "manifests" / "seq_sanity_pairs.csv", index=False, encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
