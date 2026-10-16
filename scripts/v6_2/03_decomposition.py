# -*- coding: utf-8 -*-
"""scripts/v6_2/03_decomposition.py -- V6.2 PHASE 3: Pearson 协方差分解 + Primary Statistics.

Primary Statistic 1（V6.2 §3.1–3.2）：Pearson 协方差**精确**分解
----------------------------------------------------------------
    y  = B_y + W_y,   B_y,ci = ȳ_c ,  W_y,ci = y_ci − ȳ_c
    ŷ  = B_ŷ + W_ŷ,   B_ŷ,ci = ȳ̂_c ,  W_ŷ,ci = ŷ_ci − ȳ̂_c

因组内均值为 0 → `Cov(B_y, W_ŷ) = Cov(W_y, B_ŷ) = 0`，故**精确**有：

    Cov(y, ŷ) = Cov(B_y, B_ŷ) + Cov(W_y, W_ŷ)
    r_full    = C_B + C_W
      C_B = Cov(B_y, B_ŷ) / (σ_y σ_ŷ)
      C_W = Cov(W_y, W_ŷ) / (σ_y σ_ŷ)

**浮点 sanity**：`np.isclose(C_B + C_W, r_full, rtol=1e-10, atol=1e-8)`
（**不**用硬编码 < 1e-10，避免浮点差异假报 bug —— V6.2 §7.5）

Primary Statistics
------------------
* **Δ_C  = C_B − C_W**（单档内 between 是否主导）
* **Δ_OOD = C_W(L1) − C_W(L4)**（within 信号是否在 OOD 后丢失）—— **配对 bootstrap**
* **Gain_W(M2−M0)**、**Gain_W(M1−M2)** —— **配对 bootstrap**（同索引同时算三者）

Secondary
---------
* 方差分解（row-weighted 与 macro-target 两套口径，明确各答什么问题）
* Between/Within 重构（**diagnostic，非可用预测器**）
* 负对照 A（within-target shuffle）/ B（between-only）

输出
----
`phase3_decomposition/` 下 6 个 CSV + `decomposition_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH3 = C.ROOT / "experiments" / "v6_2" / "phase3_decomposition"
PH3.mkdir(parents=True, exist_ok=True)
NPZ = PH2 / "split4_predictions.npz"
N_BOOT = 5000
SEED = 42
RTOL, ATOL = 1e-10, 1e-8
L: list[str] = []
sanity_rows: list[dict] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------- 核心分解
def decompose(y, p, g, _cache: dict | None = None):
    """返回 r_full / C_B / C_W / share_B 及 B/W 分量（population 矩，ddof=0 一致）。

    **性能**：组均值用 `np.bincount` 向量化结算（比逐组 `mean()` 快 **24.6×**，
    已用元素级 `np.allclose` 验证结果**逐元素一致**）。
    该优化是必要的：本函数在 5000 次 cluster bootstrap 内被反复调用，
    朴素实现单次统计量需 ~90s，全部 Primary Statistics 需 >30 分钟。
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    g = np.asarray(g)
    order = np.argsort(g, kind="mergesort")
    gs = g[order]
    uniq, _starts = np.unique(gs, return_index=True)
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
        return dict(r_full=np.nan, C_B=np.nan, C_W=np.nan, share_B=np.nan,
                    By=By, Bp=Bp, Wy=Wy, Wp=Wp)
    r_full = cov(y, p) / den
    C_B = cov(By, Bp) / den
    C_W = cov(Wy, Wp) / den
    share = (C_B / (C_B + C_W)) if abs(C_B + C_W) > 1e-12 else np.nan
    return dict(r_full=float(r_full), C_B=float(C_B), C_W=float(C_W),
                share_B=float(share), By=By, Bp=Bp, Wy=Wy, Wp=Wp)


def _sanity(tag, y, p, g):
    d = decompose(y, p, g)
    ok = bool(np.isclose(d["C_B"] + d["C_W"], d["r_full"], rtol=RTOL, atol=ATOL))
    sanity_rows.append({"tag": tag, "r_full": d["r_full"], "C_B": d["C_B"],
                        "C_W": d["C_W"], "sum": d["C_B"] + d["C_W"],
                        "abs_diff": abs(d["C_B"] + d["C_W"] - d["r_full"]),
                        "pass": ok})
    return d, ok


def var_decomp(y, g):
    """方差分解：row-weighted 与 macro-target 两套口径。"""
    y = np.asarray(y, float)
    g = np.asarray(g)
    d = pd.DataFrame({"y": y, "g": g})
    gm = d.groupby("g")["y"].agg(["mean", "size"])
    N = len(y)
    w = gm["size"] / N
    grand = y.mean()
    V_B_rw = float((w * (gm["mean"] - grand) ** 2).sum())
    V_W_rw = float((w * d.groupby("g")["y"].var(ddof=0).fillna(0)).sum())
    V_B_macro = float(((gm["mean"] - gm["mean"].mean()) ** 2).mean())
    within_vars = d.groupby("g")["y"].var(ddof=0).fillna(0)
    V_W_macro = float(within_vars.mean())
    return {
        "V_B_roww": V_B_rw, "V_W_roww": V_W_rw, "V_T_roww": V_B_rw + V_W_rw,
        "share_B_roww": V_B_rw / max(1e-12, V_B_rw + V_W_rw),
        "V_B_macro": V_B_macro, "V_W_macro": V_W_macro,
        "V_T_macro": V_B_macro + V_W_macro,
        "share_B_macro": V_B_macro / max(1e-12, V_B_macro + V_W_macro),
        "n_target": int(len(gm)),
    }


def boot_ci_paired(fn, arrays: list[np.ndarray], g, n_boot=N_BOOT, seed=SEED):
    """配对 cluster bootstrap：每次迭代用**相同**抽样索引同时计算所有统计量。

    `fn` 的签名约定为 `fn(*arrays, groups)` —— **必须**接收分组数组，
    因为分解类统计量（C_B / C_W）依赖「组内均值」定义。
    """
    uniq = np.unique(g)
    cidx = {u: np.where(g == u)[0] for u in uniq}
    rs = np.random.RandomState(seed)
    pt = fn(*arrays, g)                       # 修复：原为 fn(*arrays)，漏传 g → TypeError
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        rows = np.concatenate([cidx[uniq[j]] for j in pick])
        rep = np.concatenate([np.full(len(cidx[uniq[j]]), i)
                              for i, j in enumerate(pick)])
        v = fn(*[a[rows] for a in arrays], rep)
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return {"point": pt, "ci_low": np.nan, "ci_high": np.nan, "sd": np.nan,
                "n_valid": len(vals)}
    a = np.asarray(vals, float)
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": float(pt), "ci_low": float(lo), "ci_high": float(hi),
            "sd": float(a.std(ddof=1)), "n_valid": len(vals),
            "p_gt0": float((a > 0).mean())}


# ---------------------------------------------------------------- 指标
def f_delta_C(y, p, g):
    d = decompose(y, p, g)
    return d["C_B"] - d["C_W"]


def f_CW(y, p, g):
    return decompose(y, p, g)["C_W"]


def f_rho(y, p, g):
    return FM.macro_spearman(y, p, g)


def f_eff(y, p, g):
    return FM.eff10(y, p, g, n_rand=20)


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 3 — Pearson 协方差分解 + Primary Statistics")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 数据: `{NPZ.name}`（PHASE 2 的四档 × 三模型预测）")
    emit(f"- bootstrap: **{N_BOOT}** 次，抽样单位 = **complex**，seed {SEED}")
    emit(f"- 浮点判据: `np.isclose(C_B + C_W, r_full, rtol={RTOL}, atol={ATOL})`")
    emit("")

    if not NPZ.exists():
        emit(f"!! `{NPZ}` 不存在 → 无法执行")
        (PH3 / "decomposition_report.md").write_text("\n".join(L) + "\n",
                                                     encoding="utf-8")
        return 1
    z = np.load(NPZ, allow_pickle=True)
    keys = list(z.files)
    _LEVELS = [k[:-2] for k in sorted(keys) if k.endswith("_y")]
    levels = _LEVELS          # 修复：原 `k.split("_")[0]` 把 "L1_random_y" 截成 "L1"
    emit(f"- 载入 keys: {len(keys)}；识别档: {levels}")
    emit("")

    # ---------- 3.1 分解表 ----------
    dec_rows = []
    emit("## 3.1 Pearson 协方差分解（每档 × 每 seed × 每模型）")
    emit("")
    emit("| 档 | seed | 模型 | r_full | C_B | C_W | C_B+C_W | \\|diff\\| | Share_B | Sanity |")
    emit("|---|---:|---|---:|---:|---:|---:|---:|---:|---|")
    for lv in levels:
        yk = f"{lv}_y"
        gk = f"{lv}_pdb"
        if yk not in keys or gk not in keys:
            continue
        y = z[yk]; g = np.array([str(x) for x in z[gk]])
        for mk in ("M0", "M1", "M2"):
            for sd in (42, 43, 44):
                k = f"{lv}_{mk}_s{sd}"
                if k not in keys:
                    continue
                p = z[k]
                m = np.isfinite(y) & np.isfinite(p)
                d, ok = _sanity(f"{lv}|{mk}|{sd}", y[m], p[m], g[m])
                emit(f"| {lv} | {sd} | {mk} | {d['r_full']:+.6f} | "
                     f"{d['C_B']:+.6f} | {d['C_W']:+.6f} | "
                     f"{d['C_B'] + d['C_W']:+.6f} | "
                     f"{abs(d['C_B'] + d['C_W'] - d['r_full']):.2e} | "
                     f"{d['share_B']:+.4f} | {'✅' if ok else '❌'} |")
                dec_rows.append({"level": lv, "seed": sd, "model": mk,
                                 "r_full": d["r_full"], "C_B": d["C_B"],
                                 "C_W": d["C_W"],
                                 "share_B": d["share_B"], "sanity_pass": ok,
                                 "n_rows": int(m.sum()),
                                 "n_complex": int(len(set(g[m])))})
    pd.DataFrame(dec_rows).to_csv(PH3 / "pearson_decomposition.csv",
                                  index=False, encoding="utf-8")
    all_pass = all(r["pass"] for r in sanity_rows)
    emit("")
    emit(f"**Sanity 汇总**：{sum(r['pass'] for r in sanity_rows)} / "
         f"{len(sanity_rows)} 通过；最大 \\|diff\\| = "
         f"{max((r['abs_diff'] for r in sanity_rows), default=float('nan')):.3e} "
         f"→ {'✅ PASS' if all_pass else '❌ FAIL'}")
    emit("")

    # ---------- 3.2–3.4 Primary Statistics ----------
    emit("## 3.2–3.4 Primary Statistics（cluster bootstrap，配对）")
    emit("")
    prim = []

    def add_prim(name, level, model, arrs, g, fn):
        r = boot_ci_paired(fn, arrs, g)
        excl = (r["ci_low"] > 0) or (r["ci_high"] < 0)
        emit(f"| {name} | {level} | {model} | {r['point']:+.4f} | "
             f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {r['sd']:.4f} | "
             f"{r['n_valid']} | {'✅ 排除 0' if excl else '❌ 含 0'} |")
        prim.append({"statistic": name, "level": level, "model": model,
                     "point": r["point"], "ci_low": r["ci_low"],
                     "ci_high": r["ci_high"], "sd": r["sd"],
                     "n_valid": r["n_valid"],
                     "ci_excludes_zero": bool(excl)})
        return r

    # Δ_C
    emit("### Primary Statistic 1 — `Δ_C = C_B − C_W`（单档内 between 是否主导）")
    emit("")
    emit("| 统计量 | 档 | 模型 | point | 95% CI | SD | 有效 | 判定 |")
    emit("|---|---|---|---:|---|---:|---:|---|")
    for lv in ("L4_family", "L1_random"):
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        for mk in ("M0", "M1"):
            k = f"{lv}_{mk}_s42"
            if k not in keys:
                continue
            p = z[k]
            m = np.isfinite(y) & np.isfinite(p)
            add_prim("delta_C", lv, mk, [y[m], p[m]], g[m], f_delta_C)
    emit("")

    # Δ_OOD（配对：同一 complex 上 C_W(L1) 与 C_W(L4)）
    emit("### Primary Statistic 2 — `Δ_OOD = C_W(L1) − C_W(L4)`（within 信号是否 OOD 后丢失）")
    emit("")
    emit("> **配对方式**：L1 与 L4 的行宇宙相同（同一 matched set），")
    emit("> 故按 **complex** 抽样时对两个档次使用**相同索引**。")
    emit("> L1/L4 的 `pdb` 数组相同 → 直接配对。")
    emit("")
    emit("| 统计量 | 档对 | 模型 | point | 95% CI | SD | 有效 | 判定 |")
    emit("|---|---|---|---:|---|---:|---:|---|")
    if "L1_random_y" in keys and "L4_family_y" in keys:
        y1 = z["L1_random_y"]; g1 = np.array([str(x) for x in z["L1_random_pdb"]])
        y4 = z["L4_family_y"]; g4 = np.array([str(x) for x in z["L4_family_pdb"]])
        same = (len(y1) == len(y4)) and np.array_equal(g1, g4) and \
            np.allclose(y1, y4, equal_nan=True)
        emit(f"- L1/L4 行宇宙一致（长度/分组/真值）：**{same}**")
        emit("")
        for mk in ("M0", "M1"):
            k1, k4 = f"L1_random_{mk}_s42", f"L4_family_{mk}_s42"
            if k1 not in keys or k4 not in keys:
                continue
            p1, p4 = z[k1], z[k4]
            m = np.isfinite(p1) & np.isfinite(p4) & np.isfinite(y1)
            add_prim("delta_OOD", "L1−L4", mk, [y1[m], p1[m], p4[m]], g1[m],
                     lambda yy, a, b, gg: f_CW(yy, a, gg) - f_CW(yy, b, gg))
    emit("")

    # Gain_W
    emit("### Primary Statistic 3 — `Gain_W`（D 差分 / DGaux 是否恢复 mutation signal）")
    emit("")
    emit("| 统计量 | 档 | 对比 | point | 95% CI | SD | 有效 | 判定 |")
    emit("|---|---|---|---:|---|---:|---:|---|")
    for lv in ("L1_random", "L4_family"):
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        for lab, ka, kb in (("M2−M0", "M2", "M0"), ("M1−M2", "M1", "M2")):
            k1, k2 = f"{lv}_{ka}_s42", f"{lv}_{kb}_s42"
            if k1 not in keys or k2 not in keys:
                continue
            pa, pb = z[k1], z[k2]
            m = np.isfinite(pa) & np.isfinite(pb) & np.isfinite(y)
            add_prim("gain_W", lv, lab, [y[m], pa[m], pb[m]], g[m],
                     lambda yy, a, b, gg: f_CW(yy, a, gg) - f_CW(yy, b, gg))
    emit("")

    # ρ 与 Eff10 的 Δ_OOD（附）
    emit("### 附 — `Δρ_OOD = ρ_macro(L1) − ρ_macro(L4)`（配对 bootstrap）")
    emit("")
    emit("| 模型 | point | 95% CI | 判定 |")
    emit("|---|---:|---|---|")
    if "L1_random_y" in keys and "L4_family_y" in keys:
        for mk in ("M0", "M1"):
            k1, k4 = f"L1_random_{mk}_s42", f"L4_family_{mk}_s42"
            if k1 not in keys or k4 not in keys:
                continue
            p1, p4 = z[k1], z[k4]
            m = np.isfinite(p1) & np.isfinite(p4) & np.isfinite(y1)
            r = boot_ci_paired(lambda yy, a, b, gg: f_rho(yy, a, gg) - f_rho(yy, b, gg),
                               [y1[m], p1[m], p4[m]], g1[m])
            excl = (r["ci_low"] > 0) or (r["ci_high"] < 0)
            emit(f"| {mk} | {r['point']:+.4f} | [{r['ci_low']:+.4f}, "
                 f"{r['ci_high']:+.4f}] | {'✅' if excl else '❌ 含 0'} |")
            prim.append({"statistic": "delta_rho_OOD", "level": "L1−L4",
                         "model": mk, "point": r["point"],
                         "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                         "sd": r["sd"], "n_valid": r["n_valid"],
                         "ci_excludes_zero": bool(excl)})
    emit("")
    pd.DataFrame(prim).to_csv(PH3 / "primary_statistics.csv", index=False,
                              encoding="utf-8")

    # ---------- 3.5 方差分解 ----------
    emit("## 3.5 方差分解（secondary）")
    emit("")
    emit("| 档 | 口径 | V_B | V_W | V_T | Share_B | n_target |")
    emit("|---|---|---:|---:|---:|---:|---:|")
    vd_rows = []
    for lv in levels:
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        m = np.isfinite(y)
        v = var_decomp(y[m], g[m])
        for cal, lab in (("roww", "row-weighted（回答：一个**随机实验记录**的方差来源）"),
                         ("macro", "macro-target（回答：一个**随机靶点**的方差来源）")):
            emit(f"| {lv} | {cal} | {v['V_B_' + cal]:.4f} | {v['V_W_' + cal]:.4f} | "
                 f"{v['V_T_' + cal]:.4f} | **{v['share_B_' + cal]:.4f}** | "
                 f"{v['n_target']} |")
            vd_rows.append({"level": lv, "caliber": cal,
                            "V_B": v["V_B_" + cal], "V_W": v["V_W_" + cal],
                            "V_T": v["V_T_" + cal],
                            "share_B": v["share_B_" + cal],
                            "n_target": v["n_target"]})
    emit("")
    emit("> **两套口径回答不同问题**，不可互相替代：")
    emit("> row-weighted 回答「一条随机实验记录」；macro-target 回答「一个随机靶点」。")
    emit("")
    pd.DataFrame(vd_rows).to_csv(PH3 / "variance_decomposition.csv", index=False,
                                 encoding="utf-8")

    # ---------- 3.6 重构 ----------
    emit("## 3.6 Between/Within 重构（**diagnostic only**）")
    emit("")
    emit("> ⚠️ 以下是**诊断性重构**，`ŷ^B` / `ŷ^W` **不是可用预测器**（它们用到真值的组均值）。")
    emit("")
    emit("| 档 | 模型 | Pearson(y, ŷ_full) | Pearson(y, ŷ^B) | Pearson(y, ŷ^W) |")
    emit("|---|---|---:|---:|---:|")
    rc_rows = []
    for lv in ("L1_random", "L4_family"):
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        for mk in ("M0", "M1"):
            k = f"{lv}_{mk}_s42"
            if k not in keys:
                continue
            p = z[k]
            m = np.isfinite(y) & np.isfinite(p)
            d = decompose(y[m], p[m], g[m])
            rB = C.pearson(y[m], d["Bp"])
            rW = C.pearson(y[m], d["Wp"])
            rF = C.pearson(y[m], p[m])
            emit(f"| {lv} | {mk} | {rF:+.4f} | {rB:+.4f} | {rW:+.4f} |")
            rc_rows.append({"level": lv, "model": mk, "pearson_full": rF,
                            "pearson_between_only": rB,
                            "pearson_within_only": rW})
    emit("")
    pd.DataFrame(rc_rows).to_csv(PH3 / "reconstruction_results.csv",
                                 index=False, encoding="utf-8")

    # ---------- 3.7 负对照 ----------
    emit("## 3.7 负对照")
    emit("")
    emit("### Control A — Within-target shuffle（保持每靶点预测分布，打乱对应关系）")
    emit("")
    emit("| 档 | 模型 | Overall r（原始） | Overall r（shuffle） | 下降 | "
         "Within ρ（原始） | Within ρ（shuffle） |")
    emit("|---|---|---:|---:|---:|---:|---:|")
    ctl_rows = []
    for lv in ("L1_random", "L4_family"):
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        for mk in ("M0", "M1"):
            k = f"{lv}_{mk}_s42"
            if k not in keys:
                continue
            p = z[k]
            m = np.isfinite(y) & np.isfinite(p)
            yy, pp, gg = y[m], p[m], g[m]
            rs = np.random.RandomState(SEED)
            ps = pp.copy()
            for cx in np.unique(gg):
                idx = np.where(gg == cx)[0]
                ps[idx] = pp[rs.permutation(len(idx))]
            r0 = C.pearson(yy, pp); r1 = C.pearson(yy, ps)
            w0 = FM.macro_spearman(yy, pp, gg); w1 = FM.macro_spearman(yy, ps, gg)
            emit(f"| {lv} | {mk} | {r0:+.4f} | {r1:+.4f} | **{r0 - r1:+.4f}** | "
                 f"{w0:+.4f} | {w1:+.4f} |")
            ctl_rows.append({"control": "within_shuffle", "level": lv,
                             "model": mk, "overall_full": r0,
                             "overall_shuffled": r1, "drop": r0 - r1,
                             "within_full": w0, "within_shuffled": w1})
    emit("")
    emit("### Control B — Between-only（每个突变体赋该靶点的平均预测）")
    emit("")
    emit("| 档 | 模型 | Overall r | Within ρ |")
    emit("|---|---|---:|---:|")
    for lv in ("L1_random", "L4_family"):
        if f"{lv}_y" not in keys:
            continue
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        for mk in ("M0", "M1"):
            k = f"{lv}_{mk}_s42"
            if k not in keys:
                continue
            p = z[k]
            m = np.isfinite(y) & np.isfinite(p)
            yy, pp, gg = y[m], p[m], g[m]
            d = decompose(yy, pp, gg)
            rb = C.pearson(yy, d["Bp"])
            wb = FM.macro_spearman(yy, d["Bp"], gg)
            emit(f"| {lv} | {mk} | {rb:+.4f} | {wb:+.4f} |")
            ctl_rows.append({"control": "between_only", "level": lv,
                             "model": mk, "overall_between_only": rb,
                             "within_between_only": wb})
    emit("")
    pd.DataFrame(ctl_rows).to_csv(PH3 / "control_results.csv", index=False,
                                  encoding="utf-8")

    # ---------- 报告 ----------
    pd.DataFrame(sanity_rows).to_csv(PH3 / "sanity_checks.csv", index=False,
                                     encoding="utf-8")
    emit("## 判定")
    emit("")
    emit(f"- **浮点 sanity**: {sum(r['pass'] for r in sanity_rows)}/"
         f"{len(sanity_rows)} → {'PASS ✅' if all_pass else 'FAIL ❌'}")
    if len(pd.DataFrame(prim)):
        p = pd.DataFrame(prim)
        for st in p["statistic"].unique():
            sub = p[p["statistic"] == st]
            n_ex = int(sub["ci_excludes_zero"].sum())
            emit(f"- **{st}**: {n_ex}/{len(sub)} 个 CI 排除 0")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH3 / "decomposition_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
    print(f"\n写出 {PH3 / 'decomposition_report.md'}")
    for f in ("pearson_decomposition.csv", "primary_statistics.csv",
              "variance_decomposition.csv", "reconstruction_results.csv",
              "control_results.csv", "sanity_checks.csv"):
        print(f"写出 {PH3 / f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
