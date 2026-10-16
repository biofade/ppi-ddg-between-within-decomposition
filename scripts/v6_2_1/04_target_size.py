# -*- coding: utf-8 -*-
"""scripts/v6_2_1/04_target_size.py -- V6.2.1 PHASE 4: Target-size stratification.

待解释的现象
------------
L4（family-level）上：
* `macro-target` Between Pearson **> 0**（V6.2: M1 = +0.5340）
* 但 **row-weighted `C_B` ≈ 0 / 负**（V6.2: M1 = −0.0112）

假设：**target mutation-count imbalance** —— 少数 mutation-rich target 主导
row-weighted 统计，而 macro-target 对每个 target 等权。

做法（预注册，不得事后挑阈值）
------------------------------
* **A. Tertiles**（按每 target 的突变数分低/中/高三档）
* **B. Fixed thresholds**：2–9 / 10–29 / ≥30
* **C. Continuous**：`|target_mean_error|` vs `log(n_mutations)` 的 Spearman

输出
----
* `phase4_target_size/target_level_stats.csv`
* `phase4_target_size/target_size_bins.csv`
* `phase4_target_size/target_size_report.md`
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
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

P3 = C.ROOT / "experiments" / "v6_2_1" / "phase3_d_only"
PH4 = C.ROOT / "experiments" / "v6_2_1" / "phase4_target_size"
PH4.mkdir(parents=True, exist_ok=True)
DONLY = P3 / "d_only_predictions.npz"
V62 = C.ROOT / "experiments" / "v6_2" / "phase2_split4" / "split4_predictions.npz"
OUT_T = PH4 / "target_level_stats.csv"
OUT_B = PH4 / "target_size_bins.csv"
OUT_MD = PH4 / "target_size_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def decomp_terms(y, p, g):
    """返回 C_B 的分子项 V_B 与总体 Pearson。"""
    y = np.asarray(y, float); p = np.asarray(p, float); g = np.asarray(g)
    order = np.argsort(g, kind="mergesort")
    gs = g[order]
    uniq, _ = np.unique(gs, return_index=True)
    codes = np.searchsorted(uniq, gs)
    cnt = np.bincount(codes, minlength=len(uniq))
    sy = np.bincount(codes, weights=y[order], minlength=len(uniq))
    sp = np.bincount(codes, weights=p[order], minlength=len(uniq))
    By = np.empty(len(y)); Bp = np.empty(len(p))
    By[order] = sy[codes] / cnt[codes]
    Bp[order] = sp[codes] / cnt[codes]
    Wy, Wp = y - By, p - Bp

    def cov(a, b):
        return float(((a - a.mean()) * (b - b.mean())).mean())
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return dict(r_full=np.nan, C_B=np.nan, C_W=np.nan, By=By, Bp=Bp)
    return dict(r_full=cov(y, p) / den, C_B=cov(By, Bp) / den,
                C_W=cov(Wy, Wp) / den, By=By, Bp=Bp)


def main() -> int:
    t0 = time.time()
    emit("# V6.2.1 PHASE 4 — Target-size Stratification")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 目标：解释 L4 上「macro-target Between Pearson > 0」")
    emit("  与「row-weighted C_B ≈ 0/负」的并存")
    emit("- **预注册分层**：A 三分位 / B 固定阈值(2–9,10–29,≥30) / C 连续相关")
    emit("- **不事后挑阈值**")
    emit("")

    if not DONLY.exists() or not V62.exists():
        emit("!! 缺预测文件")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    z2 = np.load(DONLY, allow_pickle=True)
    z6 = np.load(V62, allow_pickle=True)

    # ---------- 4.1/4.3 每 target 统计（plot-ready） ----------
    emit("## 4.3 每 target 统计（plot-ready CSV）")
    emit("")
    rows = []
    for lv in ("L1_random", "L4_family"):
        y = z6[f"{lv}_y"]
        g = np.array([str(x) for x in z6[f"{lv}_pdb"]])
        P = {}
        for mk in ("M0", "M1"):
            ks = [f"{lv}_{mk}_s{s}" for s in (42, 43, 44)
                  if f"{lv}_{mk}_s{s}" in z6.files]
            if ks:
                P[mk] = np.mean([z6[k] for k in ks], axis=0)
        ks2 = [f"{lv}_M2_s{s}" for s in (42, 43, 44)
               if f"{lv}_M2_s{s}" in z2.files]
        if ks2:
            P["M2"] = np.mean([z2[k] for k in ks2], axis=0)
        m = np.isfinite(y)
        for v in P.values():
            m &= np.isfinite(v)
        yy = y[m].astype(float); gg = g[m]
        uq = np.unique(gg)
        for cx in uq:
            sel = gg == cx
            rec = {"level": lv, "pdb_id": str(cx),
                   "n_mutations": int(sel.sum()),
                   "true_mean_ddg": float(yy[sel].mean())}
            for mk, pv in P.items():
                pp = pv[m].astype(float)
                rec[f"pred_mean_{mk}"] = float(pp[sel].mean())
                rec[f"mean_error_{mk}"] = float(pp[sel].mean() - yy[sel].mean())
                rec[f"abs_mean_error_{mk}"] = abs(rec[f"mean_error_{mk}"])
                rec[f"within_spearman_{mk}"] = (
                    FM.macro_spearman(yy[sel], pp[sel], np.array([cx] * int(sel.sum())),
                                      min_n=2))
            rows.append(rec)
    tl = pd.DataFrame(rows)
    tl.to_csv(OUT_T, index=False, encoding="utf-8")
    emit(f"- 写出 **{len(tl)}** 行（{tl.level.nunique()} 档 × "
         f"{tl[~tl.duplicated(['level','pdb_id'])].shape[0]} 个 target-档）")
    emit("")

    # ---------- 4.1 分层 ----------
    emit("## 4.1 分层（A 三分位 / B 固定阈值）")
    emit("")
    bin_rows = []
    for lv in ("L1_random", "L4_family"):
        sub = tl[tl.level == lv].copy()
        if not len(sub):
            continue
        emit(f"### {lv}")
        emit("")
        emit(f"- target 数: **{len(sub)}**；突变数 min/median/max = "
             f"**{sub.n_mutations.min()}/{int(sub.n_mutations.median())}/"
             f"{sub.n_mutations.max()}**")
        emit("")
        # A. tertiles
        try:
            sub["bin_A"] = pd.qcut(sub.n_mutations, 3,
                                   labels=["low", "medium", "high"],
                                   duplicates="drop")
        except Exception:  # noqa: BLE001
            sub["bin_A"] = "all"
        # B. fixed
        sub["bin_B"] = pd.cut(sub.n_mutations, bins=[0, 9, 29, 10 ** 9],
                              labels=["2-9", "10-29", ">=30"])
        for scheme, col in (("A_tertile", "bin_A"), ("B_fixed", "bin_B")):
            emit(f"**方案 {scheme}**")
            emit("")
            emit("| bin | n_target | n_rows | macro mean-Pearson | mean-abs-error | "
                 "row-w C_B | row-w C_W | Overall r | Macro within-ρ |")
            emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for b in sub[col].cat.categories if hasattr(sub[col], "cat") else \
                    sorted(sub[col].dropna().unique()):
                ss = sub[sub[col] == b]
                if not len(ss):
                    continue
                y = z6[f"{lv}_y"]; g = np.array([str(x) for x in z6[f"{lv}_pdb"]])
                ks = [f"{lv}_M1_s{s}" for s in (42, 43, 44)
                      if f"{lv}_M1_s{s}" in z6.files]
                P = np.mean([z6[k] for k in ks], axis=0)
                m = np.isfinite(y) & np.isfinite(P) & np.isin(g, ss.pdb_id.values)
                yy = y[m].astype(float); pp = P[m].astype(float); gg = g[m]
                if len(yy) < 20:
                    continue
                d = decomp_terms(yy, pp, gg)
                # macro mean-Pearson（跨 target 的 均值相关）
                by = np.array([yy[gg == c].mean() for c in np.unique(gg)])
                bp = np.array([pp[gg == c].mean() for c in np.unique(gg)])
                mr = C.pearson(by, bp) if len(by) >= 3 else float("nan")
                mae = float(np.abs(ss[[c for c in ss.columns
                                       if c.startswith("abs_mean_error_")][0]]).mean())
                mw = FM.macro_spearman(yy, pp, gg, min_n=5)
                emit(f"| {b} | {len(ss)} | {len(yy)} | {mr:+.4f} | {mae:.4f} | "
                     f"{d['C_B']:+.4f} | {d['C_W']:+.4f} | {d['r_full']:+.4f} | "
                     f"{mw:+.4f} |")
                bin_rows.append({"level": lv, "scheme": scheme, "bin": str(b),
                                 "n_target": len(ss), "n_rows": len(yy),
                                 "macro_mean_pearson": mr,
                                 "mean_abs_error": mae,
                                 "roww_C_B": d["C_B"], "roww_C_W": d["C_W"],
                                 "overall_r": d["r_full"],
                                 "macro_within_spearman": mw})
            emit("")
    pd.DataFrame(bin_rows).to_csv(OUT_B, index=False, encoding="utf-8")

    # ---------- 4.1C 连续分析 ----------
    emit("## 4.1C 连续分析：`|target_mean_error|` vs `log(n_mutations)`")
    emit("")
    emit("| 档 | 模型 | Spearman | p | n_target |")
    emit("|---|---|---:|---:|---:|")
    from scipy.stats import spearmanr
    cont_rows = []
    for lv in ("L1_random", "L4_family"):
        sub = tl[tl.level == lv]
        if not len(sub):
            continue
        for mk in ("M0", "M1", "M2"):
            c = f"abs_mean_error_{mk}"
            if c not in sub.columns:
                continue
            v = sub[[c, "n_mutations"]].dropna()
            if len(v) < 5:
                continue
            rho, p = spearmanr(np.log(v.n_mutations), v[c])
            emit(f"| {lv} | {mk} | {rho:+.4f} | {p:.4f} | {len(v)} |")
            cont_rows.append({"level": lv, "model": mk, "spearman": float(rho),
                              "p": float(p), "n_target": len(v)})
    emit("")

    # ---------- 判定 ----------
    emit("## 判定")
    emit("")
    verdict = "INCONCLUSIVE"
    if bin_rows:
        b = pd.DataFrame(bin_rows)
        # 异质性判据：各 bin 的 macro_mean_pearson 或 row-w C_B 差异是否显著
        l4 = b[b.level == "L4_family"]
        if len(l4) >= 2:
            spread_cb = float(l4.roww_C_B.max() - l4.roww_C_B.min())
            spread_mr = float(l4.macro_mean_pearson.max()
                              - l4.macro_mean_pearson.min())
            emit(f"- L4 各 bin 的 row-weighted `C_B` 极差: **{spread_cb:.4f}**")
            emit(f"- L4 各 bin 的 macro mean-Pearson 极差: **{spread_mr:.4f}**")
            emit(f"- L4 每 target 突变数范围: "
                 f"**{tl[tl.level=='L4_family'].n_mutations.min()}–"
                 f"{tl[tl.level=='L4_family'].n_mutations.max()}**")
            if spread_cb > 0.05 or spread_mr > 0.15:
                verdict = "TARGET_SIZE_HETEROGENEITY_PRESENT"
            else:
                verdict = "NO_CLEAR_TARGET_SIZE_EFFECT"
    emit("")
    emit(f"### 标记：**`{verdict}`**")
    emit("")
    emit("### 对「macro Between > 0 但 row-weighted C_B ≈ 0」的解释")
    emit("")
    emit("两个统计量**权重不同**，因此可以并存，**不矛盾**：")
    emit("")
    emit("| 统计量 | 权重 | 回答的问题 |")
    emit("|---|---|---|")
    emit("| `macro`-target mean Pearson | **每个 target 等权** | 「一个随机**靶点**的均值能否被预测」 |")
    emit("| row-weighted `C_B` | **按突变数加权**（`w_c = n_c/N`） | 「一个随机**突变记录**的跨靶点分量」 |")
    emit("")
    emit("⇒ 当 mutation-rich target 的均值预测**较差/反号**、而多数小 target 尚可时，")
    emit("**大 target 会主导 row-weighted 版本并把它拉向 0 或负值**，")
    emit("而等权 macro 版本仍为正。本节的分层表正是用来检验这一机制的。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_T}")
    print(f"写出 {OUT_B}")
    print(f"写出 {OUT_MD}")
    print(f"\n>>> {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
