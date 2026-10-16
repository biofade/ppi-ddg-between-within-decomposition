# -*- coding: utf-8 -*-
"""scripts/v6_2/04_ranking_metrics.py -- V6.2 PHASE 4: 排序指标 + δ 敏感性.

### 4.1 四指标（target-macro，**不是** pair-micro）
* Spearman ρ_s（n≥5）
* Kendall τ_b（n≥2）
* Pairwise Accuracy A_c = #concordant / #total（n≥2）
* Macro 平均

### 4.2 δ-sensitive Concordance（**预注册规则，写死**）
δ ∈ {0, 0.2, 0.5}；真实 `|y_i − y_j| ≤ δ` → **不比较**；
prediction tie → 记 0.5；某 target 在 δ 后无 eligible pair → **记 NA，不进 macro mean**；
**报告 `n_eligible_targets`**。

### 4.3 n 阈值敏感性：n≥2 / n≥3 / n≥5

### 4.4 CATH 类别比较（**修正：common eligible targets**）
eligible = **四个模型 valid_targets 的交集**；Phys = mean(FoldX, Flex-ddG)；
Learned = mean(ESM-IF, ProteinMPNN)；Wilcoxon(Phys, Learned)。

**不预设「接近 53 complex」，让数据自己说话。**

输出：`phase4_ranking_metrics/` 下 4 个文件
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH4 = C.ROOT / "experiments" / "v6_2" / "phase4_ranking_metrics"
PH4.mkdir(parents=True, exist_ok=True)
NPZ = PH2 / "split4_predictions.npz"
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
DELTAS = [0.0, 0.2, 0.5]
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def kendall_tau_b(a, b) -> float:
    """Kendall τ_b，含 tie 校正。"""
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = len(a)
    if n < 2:
        return float("nan")
    conc = disc = 0
    ta = tb = 0
    for i in range(n):
        for j in range(i + 1, n):
            da = a[i] - a[j]; db = b[i] - b[j]
            if da == 0 and db == 0:
                ta += 1; tb += 1
            elif da == 0:
                ta += 1
            elif db == 0:
                tb += 1
            elif (da > 0) == (db > 0):
                conc += 1
            else:
                disc += 1
    n0 = conc + disc
    den = np.sqrt((n0 + ta) * (n0 + tb))
    return float((conc - disc) / den) if den > 0 else float("nan")


def pairwise_acc(a, b, delta: float = 0.0):
    """δ-sensitive pairwise accuracy。

    规则（预注册，写死）：
      * 真实 |y_i − y_j| ≤ δ → **不比较**
      * prediction tie → 记 0.5
      * 无 eligible pair → 返回 (nan, 0)
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = len(a)
    if n < 2:
        return float("nan"), 0
    score = []
    for i in range(n):
        for j in range(i + 1, n):
            dy = a[i] - a[j]
            if abs(dy) <= delta:
                continue
            dp = b[i] - b[j]
            if dp == 0:
                score.append(0.5)
            elif (dy > 0) == (dp > 0):
                score.append(1.0)
            else:
                score.append(0.0)
    if not score:
        return float("nan"), 0
    return float(np.mean(score)), len(score)


def macro_metrics(y, p, g, min_n: int, delta: float = 0.0):
    """返回 (macro over targets) 的 spearman / kendall / pairwise_acc 与 eligible 数。"""
    sp, kt, ac, n_el = [], [], [], 0
    for cx in np.unique(g):
        m = np.where(g == cx)[0]
        if len(m) < min_n:
            continue
        yy, pp = y[m], p[m]
        s = C.spearman(yy, pp)
        if s == s:
            sp.append(s)
        t = kendall_tau_b(yy, pp)
        if t == t:
            kt.append(t)
        a, n_pair = pairwise_acc(yy, pp, delta)
        if a == a and n_pair > 0:
            ac.append(a); n_el += 1
    return (float(np.mean(sp)) if sp else float("nan"),
            float(np.mean(kt)) if kt else float("nan"),
            float(np.mean(ac)) if ac else float("nan"), n_el)


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 4 — 排序指标 + δ 敏感性")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- δ 预注册取值: **{DELTAS}**（规则写死，不调参）")
    emit("")

    # ============ 4.1–4.3 V6.2 split4 上的排序指标 ============
    rows = []
    if NPZ.exists():
        z = np.load(NPZ, allow_pickle=True)
        keys = list(z.files)
        # 修复：原 `k.split("_")[0]` 把 "L1_random_y" 截成 "L1"
        levels = [k[:-2] for k in sorted(keys) if k.endswith("_y")]
        emit("## 4.1–4.3 split4 上的排序指标")
        emit("")
        for min_n in (2, 3, 5):
            emit(f"### n≥{min_n}")
            emit("")
            emit("| 档 | 模型 | MacroSpearman | MacroKendall | MacroAcc | n_target |")
            emit("|---|---|---:|---:|---:|---:|")
            for lv in levels:
                yk, gk = f"{lv}_y", f"{lv}_pdb"
                if yk not in keys:
                    continue
                y = z[yk]; g = np.array([str(x) for x in z[gk]])
                for mk in ("M0", "M1"):
                    ks = [f"{lv}_{mk}_s{s}" for s in (42, 43, 44)
                          if f"{lv}_{mk}_s{s}" in keys]
                    if not ks:
                        continue
                    # seed 内先平均 macro（再对 seed 平均）
                    vals = [macro_metrics(y, z[k], g, min_n) for k in ks]
                    arr = np.array([[v[0], v[1], v[2]] for v in vals], float)
                    nel = max(v[3] for v in vals)
                    emit(f"| {lv} | {mk} | {np.nanmean(arr[:, 0]):+.4f} | "
                         f"{np.nanmean(arr[:, 1]):+.4f} | "
                         f"{np.nanmean(arr[:, 2]):+.4f} | {nel} |")
                    rows.append({"min_n": min_n, "level": lv, "model": mk,
                                 "macro_spearman": float(np.nanmean(arr[:, 0])),
                                 "macro_kendall": float(np.nanmean(arr[:, 1])),
                                 "macro_pairwise_acc": float(np.nanmean(arr[:, 2])),
                                 "n_eligible_targets": int(nel)})
            emit("")
        pd.DataFrame(rows).to_csv(PH4 / "ranking_metrics.csv", index=False,
                                  encoding="utf-8")
        # n_sensitivity.md
        n_lines = ["# n 阈值敏感性（n≥2 / n≥3 / n≥5）", "",
                   "> 注意 `n≥2` 时单靶点只有 1 个 pair，Spearman 恒为 ±1，",
                   "> 因此 **MacroSpearman(n≥2) 会被 n=2 的靶点主导**，",
                   "> 这也是为何需同时报 Kendall 与 Pairwise Accuracy。", ""]
        nd = pd.DataFrame(rows)
        for lv in nd["level"].unique():
            n_lines.append(f"## {lv}")
            n_lines.append("")
            n_lines.append("| 指标 | n≥2 | n≥3 | n≥5 |")
            n_lines.append("|---|---:|---:|---:|")
            for c, lab in (("macro_spearman", "MacroSpearman"),
                           ("macro_kendall", "MacroKendall"),
                           ("macro_pairwise_acc", "MacroAcc"),
                           ("n_eligible_targets", "n_target")):
                vals = []
                for mn in (2, 3, 5):
                    s = nd[(nd["level"] == lv) & (nd["min_n"] == mn) & (nd["model"] == "M1")]
                    vals.append(f"{s[c].iloc[0]:.4f}" if len(s) else "—")
                n_lines.append(f"| {lab} | {vals[0]} | {vals[1]} | {vals[2]} |")
            n_lines.append("")
        (PH4 / "n_sensitivity.md").write_text("\n".join(n_lines) + "\n",
                                             encoding="utf-8")
    else:
        emit(f"!! `{NPZ}` 不存在 → 跳过 4.1–4.3")
        emit("")

    # ============ 4.2 δ 敏感性（在 L1/L4 的 M0/M1 上） ============
    emit("## 4.2 δ-sensitive Concordance（δ ∈ {0, 0.2, 0.5}）")
    emit("")
    d_rows = []
    if NPZ.exists():
        emit("| 档 | 模型 | δ | MacroSpearman | MacroKendall | MacroAcc(δ) | n_eligible |")
        emit("|---|---|---:|---:|---:|---:|---:|")
        for lv in ("L1_random", "L4_family"):
            if f"{lv}_y" not in z.files:
                continue
            y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
            for mk in ("M0", "M1"):
                ks = [f"{lv}_{mk}_s{s}" for s in (42, 43, 44)
                      if f"{lv}_{mk}_s{s}" in z.files]
                if not ks:
                    continue
                for dl in DELTAS:
                    vals = [macro_metrics(y, z[k], g, 2, dl) for k in ks]
                    arr = np.array([[v[0], v[1], v[2]] for v in vals], float)
                    nel = max(v[3] for v in vals)
                    emit(f"| {lv} | {mk} | {dl} | {np.nanmean(arr[:, 0]):+.4f} | "
                         f"{np.nanmean(arr[:, 1]):+.4f} | "
                         f"{np.nanmean(arr[:, 2]):+.4f} | {nel} |")
                    d_rows.append({"level": lv, "model": mk, "delta": dl,
                                   "macro_spearman": float(np.nanmean(arr[:, 0])),
                                   "macro_kendall": float(np.nanmean(arr[:, 1])),
                                   "macro_pairwise_acc": float(np.nanmean(arr[:, 2])),
                                   "n_eligible_targets": int(nel)})
    emit("")
    pd.DataFrame(d_rows).to_csv(PH4 / "delta_sensitivity.csv", index=False,
                                encoding="utf-8")
    dl_lines = ["# δ-sensitive Concordance（预注册 δ ∈ {0, 0.2, 0.5}）", "",
                "**规则（写死，不调参）**：真实 `|y_i − y_j| ≤ δ` → 不比较；",
                "prediction tie → 记 0.5；无 eligible pair 的靶点记 NA，不进 macro mean。", ""]
    if d_rows:
        dd = pd.DataFrame(d_rows)
        dl_lines += ["| 档 | 模型 | δ | MacroAcc(δ) | n_eligible |",
                     "|---|---|---:|---:|---:|"]
        for _, r in dd.iterrows():
            dl_lines.append(f"| {r['level']} | {r['model']} | {r['delta']} | "
                            f"{r['macro_pairwise_acc']:.4f} | "
                            f"{r['n_eligible_targets']} |")
    (PH4 / "delta_sensitivity.md").write_text("\n".join(dl_lines) + "\n",
                                             encoding="utf-8")

    # ============ 4.4 CATH 类别比较（common eligible targets） ============
    emit("## 4.4 CATH 类别比较（**common eligible targets**）")
    emit("")
    emit("> **修正点**：eligible = **四个模型 valid_targets 的交集**；")
    emit("> 不让只有部分模型参与的 target 进入 paired 分析。")
    emit("> **不预设**「接近 53 complex」—— 让数据自己说话。")
    emit("")
    cat_rows = []
    cat_lines = ["# CATH 类别比较（common eligible targets）", "",
                 "eligible = intersection(FoldX, Flex-ddG, ESM-IF, ProteinMPNN) 的 valid targets。",
                 "Phys = mean(FoldX, Flex-ddG)；Learned = mean(ESM-IF, ProteinMPNN)。", ""]
    if STD.exists():
        d = pd.read_csv(STD, low_memory=False)
        ca = d[d["dataset_name"] == "CATH"].copy()
        for dl in DELTAS:
            # 每模型每 target 的 pairwise acc
            per = {}
            for mo in ("FoldX", "Flex-ddG", "ESM-IF", "proteinMPNN"):
                g0 = ca[ca["model"] == mo]
                vals = {}
                for cx, gg in g0.groupby("complex_id"):
                    if len(gg) < 5:
                        continue
                    a, n_pair = pairwise_acc(gg["true_ddG"].to_numpy(float),
                                             gg["pred_ddG"].to_numpy(float), dl)
                    if a == a and n_pair > 0:
                        vals[str(cx)] = a
                per[mo] = vals
            common = set(per["FoldX"]) & set(per["Flex-ddG"]) & \
                set(per["ESM-IF"]) & set(per["proteinMPNN"])
            common = sorted(common)
            emit(f"### δ = {dl}")
            emit("")
            emit(f"- 各模型 valid targets: FoldX {len(per['FoldX'])}, "
                 f"Flex-ddG {len(per['Flex-ddG'])}, ESM-IF {len(per['ESM-IF'])}, "
                 f"proteinMPNN {len(per['proteinMPNN'])}")
            emit(f"- **common eligible targets = {len(common)}**")
            if len(common) >= 5:
                phys = np.array([(per["FoldX"][c] + per["Flex-ddG"][c]) / 2
                                 for c in common], float)
                lrnd = np.array([(per["ESM-IF"][c] + per["proteinMPNN"][c]) / 2
                                 for c in common], float)
                from scipy.stats import wilcoxon
                try:
                    w, p = wilcoxon(phys, lrnd)
                except Exception:  # noqa: BLE001
                    w, p = float("nan"), float("nan")
                emit("")
                emit("| 组 | mean MacroAcc | median |")
                emit("|---|---:|---:|")
                emit(f"| 物理（FoldX+Flex-ddG） | {phys.mean():.4f} | "
                     f"{np.median(phys):.4f} |")
                emit(f"| 序列（ESM-IF+proteinMPNN） | {lrnd.mean():.4f} | "
                     f"{np.median(lrnd):.4f} |")
                emit(f"| **Δ（物理−序列）** | **{phys.mean() - lrnd.mean():+.4f}** | "
                     f"{np.median(phys) - np.median(lrnd):+.4f} |")
                emit("")
                emit(f"- 配对 Wilcoxon: W = {w:.1f}, **p = {p:.6f}** "
                     f"({'✅ 显著' if p == p and p < 0.05 else '❌ 不显著'})")
                cat_rows.append({"delta": dl, "n_common_targets": len(common),
                                 "phys_mean": float(phys.mean()),
                                 "learned_mean": float(lrnd.mean()),
                                 "delta_mean": float(phys.mean() - lrnd.mean()),
                                 "wilcoxon_W": float(w), "wilcoxon_p": float(p)})
                cat_lines += [f"## δ = {dl}", "",
                              f"- common eligible targets: **{len(common)}**",
                              f"- 物理 mean: **{phys.mean():.4f}** | "
                              f"序列 mean: **{lrnd.mean():.4f}** | "
                              f"Δ = **{phys.mean() - lrnd.mean():+.4f}**",
                              f"- Wilcoxon p = **{p:.6f}**", ""]
            else:
                emit("- ⚠️ common targets 不足 5 → 不适合作配对检验（如实记录）")
                cat_lines += [f"## δ = {dl}", "",
                              f"- common eligible targets 仅 **{len(common)}** → 不足，跳过", ""]
            emit("")
    pd.DataFrame(cat_rows).to_csv(PH4 / "cath_extended.csv", index=False,
                                  encoding="utf-8")
    (PH4 / "cath_extended.md").write_text("\n".join(cat_lines) + "\n",
                                          encoding="utf-8")

    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH4 / "ranking_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH4}/ 下: ranking_metrics.csv, n_sensitivity.md, "
          f"delta_sensitivity.{'{csv,md}'}, cath_extended.{'{csv,md}'}, "
          f"ranking_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
