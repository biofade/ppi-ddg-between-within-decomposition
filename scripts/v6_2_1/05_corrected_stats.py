# -*- coding: utf-8 -*-
"""scripts/v6_2_1/05_corrected_stats.py -- V6.2.1 PHASE 5: 修正后的核心统计.

整合四项修复：
1. **正确的** within-target shuffle（PHASE 1）
2. **真正的** M2 = D-only（PHASE 3）
3. **已审计**的树数（PHASE 2：实际统一 300）
4. **target-size** 分层（PHASE 4）

并输出 Primary route statistics（全部 paired cluster bootstrap 5000/5000）。

**关键纪律**（V6.2.1 §关键科学纪律）：不得再用 shuffle 单独推出
「global Pearson 全部来自 between-target」；必须**同时**报告
`C_B/C_W` 的 before/after。
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
PH5 = C.ROOT / "experiments" / "v6_2_1" / "phase5_corrected_stats"
PH5.mkdir(parents=True, exist_ok=True)
DONLY = P3 / "d_only_predictions.npz"
V62 = C.ROOT / "experiments" / "v6_2" / "phase2_split4" / "split4_predictions.npz"
SHUF = C.ROOT / "experiments" / "v6_2_1" / "phase1_shuffle_audit" / "shuffle_repeats.csv"
OUT_P = PH5 / "primary_table.csv"
OUT_R = PH5 / "route_statistics.csv"
OUT_S = PH5 / "corrected_shuffle.csv"
OUT_MD = PH5 / "corrected_stats_report.md"
N_BOOT = 5000
SEEDS = (42, 43, 44)
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def decompose(y, p, g):
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
        return dict(r_full=np.nan, C_B=np.nan, C_W=np.nan)
    return dict(r_full=cov(y, p) / den, C_B=cov(By, Bp) / den,
                C_W=cov(Wy, Wp) / den)


def boot_paired(fn, arrays, g, n_boot=N_BOOT, seed=42):
    uniq = np.unique(g)
    cidx = {u: np.where(g == u)[0] for u in uniq}
    rs = np.random.RandomState(seed)
    pt = fn(*arrays, g)
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        rows = np.concatenate([cidx[uniq[j]] for j in pick])
        rep = np.concatenate([np.full(len(cidx[uniq[j]]), i)
                              for i, j in enumerate(pick)])
        v = fn(*[a[rows] for a in arrays], rep)
        if v == v:
            vals.append(v)
    a = np.asarray(vals, float)
    if len(a) < 20:
        return {"point": float(pt), "ci_low": np.nan, "ci_high": np.nan,
                "n_valid": len(a)}
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": float(pt), "ci_low": float(lo), "ci_high": float(hi),
            "n_valid": len(a)}


def f_CW(y, p, g):
    return decompose(y, p, g)["C_W"]


def main() -> int:
    t0 = time.time()
    emit("# V6.2.1 PHASE 5 — 修正后的核心统计")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 整合：正确 shuffle + 真 M2(D-only) + 树数审计(300) + target-size 分层")
    emit("- 全部区间为 **paired cluster bootstrap {}/{}**（单位 = complex）".format(
        N_BOOT, N_BOOT))
    emit("")

    z6 = np.load(V62, allow_pickle=True)
    z2 = np.load(DONLY, allow_pickle=True)
    P: dict[str, dict[str, np.ndarray]] = {}
    truth: dict[str, np.ndarray] = {}
    grp: dict[str, np.ndarray] = {}
    for lv in ("L1_random", "L4_family"):
        y = z6[f"{lv}_y"]
        g = np.array([str(x) for x in z6[f"{lv}_pdb"]])
        cur = {}
        for mk in ("M0", "M1"):
            ks = [f"{lv}_{mk}_s{s}" for s in SEEDS if f"{lv}_{mk}_s{s}" in z6.files]
            if ks:
                cur[mk] = np.mean([z6[k] for k in ks], axis=0)
        ks2 = [f"{lv}_M2_s{s}" for s in SEEDS if f"{lv}_M2_s{s}" in z2.files]
        if ks2:
            cur["M2"] = np.mean([z2[k] for k in ks2], axis=0)
        m = np.isfinite(y)
        for v in cur.values():
            m &= np.isfinite(v)
        truth[lv] = y[m].astype(float)
        grp[lv] = g[m]
        P[lv] = {k: v[m].astype(float) for k, v in cur.items()}

    # ---------- 5.1 Primary table ----------
    emit("## 5.1 Primary table（seed ensemble）")
    emit("")
    emit("| Split | Model | r_full | C_B | C_W | Macro ρ | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|")
    rows = []
    for lv in ("L1_random", "L4_family"):
        yy, gg = truth[lv], grp[lv]
        for mk in ("M0", "M2", "M1"):
            if mk not in P[lv]:
                continue
            pp = P[lv][mk]
            d = decompose(yy, pp, gg)
            mr = FM.macro_spearman(yy, pp, gg, min_n=5)
            e10 = FM.eff10(yy, pp, gg, n_rand=20)
            emit(f"| {lv} | {mk} | {d['r_full']:+.4f} | {d['C_B']:+.4f} | "
                 f"{d['C_W']:+.4f} | {mr:+.4f} | {e10:+.4f} |")
            rows.append({"split": lv, "model": mk, "r_full": d["r_full"],
                         "C_B": d["C_B"], "C_W": d["C_W"],
                         "macro_spearman": mr, "eff10": e10,
                         "n_rows": len(yy), "n_target": len(set(gg))})
    pd.DataFrame(rows).to_csv(OUT_P, index=False, encoding="utf-8")
    emit("")

    # ---------- 5.2 Primary route statistics ----------
    emit("## 5.2 Primary route statistics（paired cluster bootstrap）")
    emit("")
    emit("| 统计量 | 对比 | point | 95% CI | 有效 | 排除 0? |")
    emit("|---|---|---:|---|---:|---|")
    rrows = []

    def rec(name, contrast, fn, arrays, g):
        r = boot_paired(fn, arrays, g)
        excl = (r["ci_low"] > 0) or (r["ci_high"] < 0)
        emit(f"| {name} | {contrast} | {r['point']:+.4f} | "
             f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {r['n_valid']} | "
             f"{'✅' if excl else '❌ 含 0'} |")
        rrows.append({"statistic": name, "contrast": contrast,
                      "point": r["point"], "ci_low": r["ci_low"],
                      "ci_high": r["ci_high"], "n_valid": r["n_valid"],
                      "ci_excludes_zero": bool(excl)})
        return r

    # A. OOD loss of within signal
    y1, g1 = truth["L1_random"], grp["L1_random"]
    y4, g4 = truth["L4_family"], grp["L4_family"]
    same = (len(y1) == len(y4)) and np.array_equal(g1, g4)
    emit(f"- L1/L4 row universe 相同: **{same}**")
    emit("")
    for mk in ("M0", "M2", "M1"):
        if mk not in P["L1_random"] or mk not in P["L4_family"]:
            continue
        rec("Delta_OOD_W = C_W(L1)−C_W(L4)", mk,
            lambda y, a, b, gg: f_CW(y, a, gg) - f_CW(y, b, gg),
            [y1, P["L1_random"][mk], P["L4_family"][mk]], g1)

    emit("")
    # B/C. Gain_W
    for lv in ("L1_random", "L4_family"):
        yy, gg = truth[lv], grp[lv]
        for lab, a, b in (("Gain_D = C_W(M2)−C_W(M0)", "M2", "M0"),
                          ("Gain_DGaux = C_W(M1)−C_W(M2)", "M1", "M2"),
                          ("Gain_total = C_W(M1)−C_W(M0)", "M1", "M0")):
            if a not in P[lv] or b not in P[lv]:
                continue
            rec(lab, lv,
                lambda y, x, z, gg: f_CW(y, x, gg) - f_CW(y, z, gg),
                [yy, P[lv][a], P[lv][b]], gg)
    pd.DataFrame(rrows).to_csv(OUT_R, index=False, encoding="utf-8")
    emit("")

    # ---------- 5.3 Corrected shuffle statistic ----------
    emit("## 5.3 Corrected shuffle（同时报告 C_B/C_W before/after）")
    emit("")
    emit("> **不得**用 shuffle 单独推出「global Pearson 全部来自 between-target」。")
    emit("")
    emit("| 档 | 模型 | r_full | C_B | C_W | r_shuf mean±SD | C_B_shuf | "
         "C_W_shuf | ShuffleLoss |")
    emit("|---|---|---:|---:|---:|---|---:|---:|---:|")
    srows = []
    if SHUF.exists():
        sh = pd.read_csv(SHUF)
        for _, r in sh.iterrows():
            lv, mk = str(r["level"]), str(r["model"])
            if lv not in truth or mk not in P[lv]:
                continue
            yy, gg, pp = truth[lv], grp[lv], P[lv][mk]
            d = decompose(yy, pp, gg)
            loss = d["r_full"] - r["r_shuffle_mean"]
            emit(f"| {lv} | {mk} | {d['r_full']:+.4f} | {d['C_B']:+.4f} | "
                 f"{d['C_W']:+.4f} | {r['r_shuffle_mean']:+.4f} ± "
                 f"{r['r_shuffle_sd']:.4f} | {r['C_B_shuffle_mean']:+.4f} | "
                 f"{r['C_W_shuffle_mean']:+.4f} | **{loss:+.4f}** |")
            srows.append({"level": lv, "model": mk, "r_full": d["r_full"],
                          "C_B": d["C_B"], "C_W": d["C_W"],
                          "r_shuffle_mean": r["r_shuffle_mean"],
                          "r_shuffle_sd": r["r_shuffle_sd"],
                          "C_B_shuffle_mean": r["C_B_shuffle_mean"],
                          "C_W_shuffle_mean": r["C_W_shuffle_mean"],
                          "ShuffleLoss": loss})
    pd.DataFrame(srows).to_csv(OUT_S, index=False, encoding="utf-8")
    emit("")
    emit("**读法**：`ShuffleLoss = r_full − r_shuffle_mean`")
    emit("= **组内排序的真实增量贡献**（因为 `C_B` 在 shuffle 下守恒，被减掉了）。")
    emit("")

    # ---------- 判定汇总 ----------
    emit("## 5.4 修正后的判定要点")
    emit("")
    rr = pd.DataFrame(rrows)
    emit("| 判据 | 值 | 含义 |")
    emit("|---|---:|---|")
    if len(rr):
        for lab in ("Delta_OOD_W = C_W(L1)−C_W(L4)",
                    "Gain_D = C_W(M2)−C_W(M0)",
                    "Gain_DGaux = C_W(M1)−C_W(M2)"):
            s = rr[rr.statistic == lab]
            if not len(s):
                continue
            n_ex = int(s.ci_excludes_zero.sum())
            emit(f"| {lab} | {n_ex}/{len(s)} CI 排除 0 | "
                 f"{'有统计支持的恢复' if n_ex else '**无**统计支持的恢复'} |")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_P}")
    print(f"写出 {OUT_R}")
    print(f"写出 {OUT_S}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
