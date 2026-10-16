# -*- coding: utf-8 -*-
"""scripts/v6_1/06_v5_vs_cath_distribution.py -- V6.1 PHASE 4: V5 vs CATH 分布对比.

目标
----
V5 与 CATH 是**不同数据集**（PHASE 1.2 已证交集仅 8 个 PDB），因此
**不能**做同集比较；但可以对比 **per-complex Spearman 的分布形状**：
我们的模型是否整体性地低于物理方法？

步骤
----
1. V5：从 `p4_predictions.npz` 取 61 个 complex 的 per-complex Spearman（3 seeds）
2. CATH：21 个 complex × 5 模型的 per-complex Spearman
3. 分布对比 + **Mann-Whitney U**（此处**允许**使用：V5 与 CATH 的样本
   **不是配对**的 —— 它们来自不同数据集、不同 complex，符合 U 检验的
   「独立两组」前提。这与 PHASE 1.3 的配对情形**本质不同**。）
4. 同时给出 median / IQR / n

输出
----
* `experiments/v6_1/phase4_v5_comparison/v5_cath_distribution.csv`
* `experiments/v6_1/phase4_v5_comparison/distribution_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH4 = C.ROOT / "experiments" / "v6_1" / "phase4_v5_comparison"
PH4.mkdir(parents=True, exist_ok=True)
OUT_CSV = PH4 / "v5_cath_distribution.csv"
OUT_MD = PH4 / "distribution_report.md"
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
MIN_PS_N = 5
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
    emit("# V6.1 PHASE 4 — V5 baseline vs CATH 模型的 per-complex Spearman 分布")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    emit("> **统计方法说明**：此处**允许**使用 Mann-Whitney U —— V5 与 CATH 来自")
    emit("> **不同数据集、不同 complex**，两组**不配对**，符合 U 检验的「独立两组」前提。")
    emit("> 这与 PHASE 1.3（同一批 complex 的两类模型，必须配对）**本质不同**。")
    emit("")

    rows = []
    # ---------- V5 ----------
    z = C.V5 / "p4_conformal" / "p4_predictions.npz"
    v5_vals: list[float] = []
    if z.exists():
        zz = np.load(z, allow_pickle=True)
        y = zz["test_y"]
        g = np.array([str(x) for x in zz["test_group"]])
        for key in ("plain_s42", "plain_s43", "plain_s44"):
            if key not in zz.files:
                continue
            p = zz[key]
            for cx in np.unique(g):
                m = np.where(g == cx)[0]
                if len(m) < MIN_PS_N:
                    continue
                v = C.spearman(y[m], p[m])
                if v == v:
                    v5_vals.append(float(v))
                    rows.append({"source": "V5_baseline", "model": "V5 baseline",
                                 "seed_key": key, "complex_id": str(cx),
                                 "per_complex_spearman": round(float(v), 4)})
    va = np.array(v5_vals, float)
    emit(f"## V5 baseline")
    emit("")
    emit(f"- per-complex Spearman 样本数: **{len(va)}**"
         f"（{len(set(r['complex_id'] for r in rows))} complex × 3 seeds）")
    emit(f"- median = **{np.median(va):+.4f}**，IQR = "
         f"[{np.percentile(va, 25):+.4f}, {np.percentile(va, 75):+.4f}]，"
         f"mean = {va.mean():+.4f}")
    emit("")

    # ---------- CATH ----------
    emit("## CATH 各模型")
    emit("")
    emit("| 模型 | n | median | IQR | mean | vs V5: Δmedian | Mann-Whitney p | 显著? |")
    emit("|---|---:|---:|---|---:|---:|---:|---|")
    cath_vals: dict[str, np.ndarray] = {}
    if STD.exists():
        d = pd.read_csv(STD, low_memory=False)
        ca = d[d["dataset_name"] == "CATH"]
        for mo in sorted(ca["model"].unique()):
            vals = []
            for cx, gg in ca[ca["model"] == mo].groupby("complex_id"):
                if len(gg) < MIN_PS_N:
                    continue
                v = C.spearman(gg["true_ddG"].to_numpy(float),
                               gg["pred_ddG"].to_numpy(float))
                if v == v:
                    vals.append(float(v))
                    rows.append({"source": "CATH", "model": mo, "seed_key": "n/a",
                                 "complex_id": str(cx),
                                 "per_complex_spearman": round(float(v), 4)})
            arr = np.array(vals, float)
            cath_vals[mo] = arr
            if len(arr) < 3 or len(va) < 3:
                emit(f"| {mo} | {len(arr)} | {np.median(arr):+.4f} | — | — | — | — | — |")
                continue
            try:
                _u, p = mannwhitneyu(arr, va, alternative="two-sided")
            except Exception:  # noqa: BLE001
                p = float("nan")
            dm = float(np.median(arr) - np.median(va))
            emit(f"| {mo} | {len(arr)} | {np.median(arr):+.4f} | "
                 f"[{np.percentile(arr, 25):+.4f}, {np.percentile(arr, 75):+.4f}] | "
                 f"{arr.mean():+.4f} | {dm:+.4f} | {p:.4f} | "
                 f"{'✅' if p < 0.05 else '❌'} |")
    emit("")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------- 类别合并 ----------
    emit("## 类别合并对比")
    emit("")
    phys = np.concatenate([cath_vals.get(m, np.array([])) for m in ("Flex-ddG", "FoldX")]) \
        if cath_vals else np.array([])
    seq = np.concatenate([cath_vals.get(m, np.array([])) for m in ("ESM-IF", "proteinMPNN")]) \
        if cath_vals else np.array([])
    emit("| 组 | n | median | mean | vs V5 median Δ | Mann-Whitney p |")
    emit("|---|---:|---:|---:|---:|---:|")
    for name, arr in (("V5 baseline", va), ("CATH 物理/能量类", phys),
                      ("CATH 纯序列类", seq)):
        if len(arr) < 3:
            emit(f"| {name} | {len(arr)} | — | — | — | — |")
            continue
        if len(va) >= 3 and name != "V5 baseline":
            try:
                _u, p = mannwhitneyu(arr, va, alternative="two-sided")
            except Exception:  # noqa: BLE001
                p = float("nan")
        else:
            p = float("nan")
        emit(f"| {name} | {len(arr)} | {np.median(arr):+.4f} | {arr.mean():+.4f} | "
             f"{np.median(arr) - np.median(va):+.4f} | {p:.4f} |")
    emit("")

    # ---------- 结论 ----------
    emit("## 结论")
    emit("")
    if len(phys) >= 3 and len(va) >= 3:
        try:
            _u, p_phys = mannwhitneyu(phys, va, alternative="two-sided")
        except Exception:  # noqa: BLE001
            p_phys = float("nan")
        emit(f"- V5 baseline median per-complex Spearman = **{np.median(va):+.4f}**")
        emit(f"- CATH 物理/能量类 median = **{np.median(phys):+.4f}**"
             f"（Δ = **{np.median(phys) - np.median(va):+.4f}**），"
             f"Mann-Whitney **p = {p_phys:.4f}**")
        if p_phys == p_phys and p_phys < 0.05:
            emit("")
            emit("> ### ✅ **跨数据集支持**：物理/能量类方法的 per-complex Spearman")
            emit("> 分布显著高于我们的 V5 baseline。")
            emit("> **但必须标注**：这是**跨数据集**比较（不同来源、不同 complex），")
            emit("> **不是**同集比较；V5 与 CATH 的交集仅 8 个 PDB（PHASE 1.2）。")
        else:
            emit("")
            emit("> ### ⚠️ 差异不显著（可能受样本量限制）")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
