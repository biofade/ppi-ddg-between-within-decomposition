# -*- coding: utf-8 -*-
"""scripts/v6_1/03_category_paired_test.py -- V6.1 PHASE 1.3: 类别配对检验.

问题
----
「物理/能量类方法」是否显著优于「纯序列表示类方法」？

**关键方法学要求（V6.1 明确禁止 Mann-Whitney U）**
--------------------------------------------------
42 个 per-complex 分数（21 complex × 2 类）**不是独立样本**：
同一个 complex 的物理分与序列分来自**同一批突变体**，必须**配对**。
只用 complex 作为独立单位 → 配对的 n = 21。

方法
----
1. 对每个 complex c：
   `phys[c] = mean(Flex-ddG[c], FoldX[c])`
   `seq[c]  = mean(ESM-IF[c], proteinMPNN[c])`
   （DDGPred 单列另做一组敏感性分析）
2. 得到 **21 对** (phys, seq)，配对差值 `d[c] = phys[c] − seq[c]`
3. **方法 A**：`scipy.stats.wilcoxon(d)` —— Wilcoxon 符号秩检验（配对）
4. **方法 B**：cluster bootstrap（抽 complex，5000 次）→ Δ 的 95% CI
5. 仅保留两侧都有预测的 complex（交集 complex）

**只用 complex 内 Spearman（per-complex）**，因为这是「靶点内排序能力」的
自然配对单位；不混用全局指标（那会破坏配对结构）。

输出
----
* `experiments/v6_1/phase1_ci/category_paired_report.md`
* `experiments/v6_1/phase1_ci/category_paired_values.csv`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH1 = C.ROOT / "experiments" / "v6_1" / "phase1_ci"
PH1.mkdir(parents=True, exist_ok=True)
OUT_MD = PH1 / "category_paired_report.md"
OUT_CSV = PH1 / "category_paired_values.csv"
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"

N_BOOT = 5000
SEED = 42
MIN_COMPLEX_N = 5            # per-complex Spearman 的最小突变数
PHYS = ["Flex-ddG", "FoldX"]
SEQ = ["ESM-IF", "proteinMPNN"]
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def per_cx_spearman(df: pd.DataFrame, model: str) -> dict[str, float]:
    out = {}
    g0 = df[df["model"] == model]
    for cx, gg in g0.groupby("complex_id"):
        if len(gg) < MIN_COMPLEX_N:
            continue
        v = C.spearman(gg["true_ddG"].to_numpy(float),
                       gg["pred_ddG"].to_numpy(float))
        if v == v:
            out[str(cx)] = float(v)
    return out


def main() -> int:
    t0 = time.time()
    emit("# V6.1 PHASE 1.3 — 类别配对检验（物理/能量类 vs 纯序列表示类）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 数据: CATH held-out，`{STD.name}`，per-complex Spearman"
         f"（每 complex ≥{MIN_COMPLEX_N} 条突变）")
    emit("")
    emit("> **方法学说明**：本检验**不使用 Mann-Whitney U**（V6.1 明确禁止）——")
    emit("> 42 个 per-complex 分数不是独立样本，同一 complex 的两类分数来自")
    emit("> **同一批突变体**，必须**配对**。独立单位是 **complex**（n=21）。")
    emit("")

    if not STD.exists():
        emit("!! `standardized_preds.csv` 不存在 → 无法执行")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    d = pd.read_csv(STD, low_memory=False)
    ca = d[d["dataset_name"] == "CATH"].copy()
    emit(f"- CATH: {len(ca)} 行 / {ca['complex_id'].nunique()} complex")
    emit("")

    # ---------- 逐 complex 的 per-complex Spearman ----------
    per_model: dict[str, dict[str, float]] = {}
    for mo in PHYS + SEQ + ["DDGPred"]:
        per_model[mo] = per_cx_spearman(ca, mo)
        emit(f"- `{mo}`: 可用 complex **{len(per_model[mo])}**")
    emit("")

    # 配对用的交集 complex
    common = None
    for mo in PHYS + SEQ:
        s = set(per_model[mo])
        common = s if common is None else (common & s)
    common = sorted(common or [])
    emit(f"## 配对样本：**{len(common)} 个 complex**（两类模型均有 per-complex 分数）")
    emit("")
    if len(common) < 5:
        emit("!! 配对样本过少（<5），无法做检验")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1

    rows = []
    emit("| complex | " + " | ".join(PHYS) + " | " + " | ".join(SEQ)
         + " | phys_mean | seq_mean | delta |")
    emit("|---|" + "---:|" * (len(PHYS) + len(SEQ) + 3))
    phys_v, seq_v = [], []
    for cx in common:
        pv = float(np.mean([per_model[m][cx] for m in PHYS]))
        sv = float(np.mean([per_model[m][cx] for m in SEQ]))
        phys_v.append(pv); seq_v.append(sv)
        row = {"complex_id": cx, "phys_mean": pv, "seq_mean": sv,
               "delta": pv - sv}
        for m in PHYS + SEQ:
            row[m] = per_model[m][cx]
        rows.append(row)
        emit(f"| `{cx}` | " + " | ".join(f"{per_model[m][cx]:+.3f}" for m in PHYS)
             + " | " + " | ".join(f"{per_model[m][cx]:+.3f}" for m in SEQ)
             + f" | {pv:+.4f} | {sv:+.4f} | **{pv - sv:+.4f}** |")    pv = np.asarray(phys_v, float)
    sv = np.asarray(seq_v, float)
    dv = pv - sv
    vals = pd.DataFrame(rows)
    vals.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("")

    # ---------- 汇总 ----------
    emit("## 汇总（配对，n = {}）".format(len(dv)))
    emit("")
    emit("| 量 | 值 |")
    emit("|---|---:|")
    emit(f"| phys 组均值（per-complex Spearman） | **{pv.mean():+.4f}** |")
    emit(f"| seq  组均值 | **{sv.mean():+.4f}** |")
    emit(f"| 配对差 delta 均值 | **{dv.mean():+.4f}** |")
    emit(f"| delta 中位 | {np.median(dv):+.4f} |")
    emit(f"| delta SD | {dv.std(ddof=1):.4f} |")
    emit(f"| delta > 0 的 complex 数 | **{int((dv > 0).sum())} / {len(dv)}** |")
    emit("")

    # ---------- 方法 A: Wilcoxon 符号秩 ----------
    emit("## 方法 A — Wilcoxon 符号秩检验（配对）")
    emit("")
    try:
        if np.allclose(dv, 0):
            emit("- 所有差值均为 0 → 检验无定义")
            w_stat, w_p = float("nan"), float("nan")
        else:
            w_stat, w_p = wilcoxon(pv, sv)
            emit(f"- 统计量 W = **{w_stat:.1f}**")
            emit(f"- **p = {w_p:.6f}**"
                 + ("  → **显著（p < 0.05）**" if w_p < 0.05 else "  → 不显著"))
            # 效应量 r = Z / sqrt(N)
            n = len(dv)
            try:
                from scipy.stats import norm
                # 正态近似 Z
                z = (w_stat - n * (n + 1) / 4) / np.sqrt(n * (n + 1) * (2 * n + 1) / 24)
                r_eff = abs(z) / np.sqrt(n)
                emit(f"- 效应量 r = |Z|/√N = **{r_eff:.4f}**"
                     f"（0.1 小 / 0.3 中 / 0.5 大）")
            except Exception:  # noqa: BLE001
                pass
    except Exception as exc:  # noqa: BLE001
        emit(f"- ⚠️ 检验失败: {type(exc).__name__}: {exc}")
        w_stat, w_p = float("nan"), float("nan")
    emit("")

    # ---------- 方法 B: cluster bootstrap ----------
    emit("## 方法 B — Cluster bootstrap（抽 complex，5000 次）")
    emit("")
    rs = np.random.RandomState(SEED)
    n = len(dv)
    deltas = []
    for _ in range(N_BOOT):
        idx = rs.randint(0, n, n)
        deltas.append(dv[idx].mean())
    deltas = np.asarray(deltas, float)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p_boot = float((deltas <= 0).mean())
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| Δ (phys − seq) point | **{dv.mean():+.4f}** |")
    emit(f"| 95% CI | **[{lo:+.4f}, {hi:+.4f}]** |")
    emit(f"| bootstrap SD | {deltas.std(ddof=1):.4f} |")
    emit(f"| 单侧 p（H0: Δ≤0） | **{p_boot:.6f}** |")
    emit(f"| CI 是否排除 0 | {'✅ 是' if lo > 0 else '❌ 否'} |")
    emit("")

    # ---------- 敏感性分析：DDGPred 替代 phys 组 ----------
    emit("## 敏感性分析 — 用 DDGPred 单列代表「几何/序列混合类」")
    emit("")
    dd = per_model.get("DDGPred", {})
    common2 = sorted(set(dd) & set(per_model[SEQ[0]]) & set(per_model[SEQ[1]]))
    if len(common2) >= 5:
        d2 = np.array([dd[c] - float(np.mean([per_model[m][c] for m in SEQ]))
                       for c in common2], float)
        try:
            _w2, p2 = wilcoxon([dd[c] for c in common2],
                               [float(np.mean([per_model[m][c] for m in SEQ]))
                                for c in common2])
        except Exception:  # noqa: BLE001
            p2 = float("nan")
        emit(f"- n = {len(d2)}；Δ(DDGPred − seq) 均值 = **{d2.mean():+.4f}**；"
             f"Wilcoxon p = **{p2:.6f}**")
    else:
        emit(f"- 可用配对过少（{len(common2)}）→ 跳过")
    emit("")

    # ---------- 结论 ----------
    emit("## 结论")
    emit("")
    sig = (w_p == w_p) and (w_p < 0.05)
    emit(f"- 物理/能量类（{'+'.join(PHYS)}）per-complex Spearman 均值 "
         f"**{pv.mean():+.4f}**")
    emit(f"- 纯序列表示类（{'+'.join(SEQ)}）per-complex Spearman 均值 "
         f"**{sv.mean():+.4f}**")
    emit(f"- 配对差 **{dv.mean():+.4f}**，95% CI **[lo={lo:+.4f}, hi={hi:+.4f}]**，"
         f"Wilcoxon **p = {w_p:.6f}**")
    emit("")
    if sig and lo > 0:
        emit("> ### ✅ **类别差异显著**：物理/能量类在 within-target 排序上"
             "显著优于纯序列表示类。")
        emit(f"> 方向一致性：**{int((dv > 0).sum())}/{len(dv)}** 个 complex 上"
             f"物理类更好。")
    elif lo > 0:
        emit("> ### ⚠️ **方向一致但统计不显著**：bootstrap CI 排除 0，"
             "但 Wilcoxon 未达显著（n=21 功效受限）。")
    else:
        emit("> ### ❌ **类别差异不显著**：CI 含 0。")
    emit("")
    emit(f"> ⚠️ **样本量限制**：n = **{len(dv)}** 个 complex。"
         f"V6.1 PHASE 1.3 明确要求「样本量不足时明确标注」。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    print(f"写出 {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
