# -*- coding: utf-8 -*-
"""scripts/v7_1/03_bootstrap.py -- V7.1 PHASE 7: 配对 cluster bootstrap + 分层 + weighting.

Primary comparisons（§PHASE 7）：P1−P0 / P2−P0 / P3−P0 / P3−P1 / P3−P2
指标：Overall Pearson / RMSE / macro target-mean Pearson / attenuation slope /
      C_W / Macro Spearman / Efficiency@10

含 target-size 分层（2–9 / 10–29 / ≥30）与 5 种 weighting 的 `C_B`（诊断用）。
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

V71 = C.ROOT / "experiments" / "v7_1"
PH5, PH7 = V71 / "phase5_lsr", V71 / "phase7_bootstrap"
PH7.mkdir(parents=True, exist_ok=True)
NPZ = PH5 / "lsr_predictions.npz"
OUT_C = PH7 / "contrasts.csv"
OUT_S = PH7 / "stratification.csv"
OUT_W = PH7 / "weighting.csv"
OUT_MD = PH7 / "bootstrap_report.md"
N_BOOT = 3000
MIN_PS_N = 5
L: list[str] = []

CONTRASTS = [("p_P1_location_only", "p_P0_raw_ranker", "P1-P0 location gain"),
             ("p_P2_scale_only", "p_P0_raw_ranker", "P2-P0 scale gain"),
             ("p_P3_full_LSR", "p_P0_raw_ranker", "P3-P0 full LSR"),
             ("p_P3_full_LSR", "p_P1_location_only", "P3-P1 scale incremental"),
             ("p_P3_full_LSR", "p_P2_scale_only", "P3-P2 location incremental")]


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def _parts(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    return uq, By, Bp, cnt


def m_pearson(y, p, g):
    return C.pearson(y, p)


def m_rmse(y, p, g):
    return float(np.sqrt(np.mean((y - p) ** 2)))


def m_tmean(y, p, g):
    _u, By, Bp, _c = _parts(y, p, g)
    return C.pearson(By, Bp) if len(By) >= 3 else np.nan


def m_slope(y, p, g):
    _u, By, Bp, _c = _parts(y, p, g)
    return (float(np.polyfit(By, Bp, 1)[0])
            if len(By) >= 3 and By.std() > 0 else np.nan)


def m_cw(y, p, g):
    uq, By, Bp, cnt = _parts(y, p, g)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return np.nan
    return float((cnt / cnt.sum() * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den


def m_ms(y, p, g):
    return FM.macro_spearman(y, p, g, min_n=MIN_PS_N)


def m_eff(y, p, g):
    return FM.eff10(y, p, g, n_rand=5)


FN = {"Overall_r": m_pearson, "RMSE": m_rmse, "target_mean_r": m_tmean,
      "atten_slope": m_slope, "C_W": m_cw, "Macro_spearman": m_ms, "Eff10": m_eff}


def pair_boot(fn, yy, pa, pb, gg, n_boot=N_BOOT, seed=42):
    uq = np.unique(gg)
    cidx = [np.where(gg == u)[0] for u in uq]
    rs = np.random.RandomState(seed)
    pt = fn(yy, pa, gg) - fn(yy, pb, gg)
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uq), len(uq))
        r_ = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i)
                              for i, j in enumerate(pick)])
        v = fn(yy[r_], pa[r_], rep) - fn(yy[r_], pb[r_], rep)
        if v == v:
            vals.append(v)
    a = np.asarray(vals, float)
    if len(a) < 20:
        return dict(point=pt, lo=np.nan, hi=np.nan, n=len(a), excl=False)
    lo, hi = np.percentile(a, [2.5, 97.5])
    return dict(point=float(pt), lo=float(lo), hi=float(hi), n=len(a),
                excl=bool((lo > 0) or (hi < 0)))


def main() -> int:
    t0 = time.time()
    emit("# V7.1 PHASE 7 — Bootstrap + 分层 + Weighting")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- paired cluster bootstrap **{N_BOOT}/{N_BOOT}**，单位 = complex")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    emit(f"- n = {len(y)} rows / {len(set(g))} targets")
    emit("")

    rows = []
    for a, b, lab in CONTRASTS:
        if a not in z.files or b not in z.files:
            emit(f"- {lab}: 缺预测")
            continue
        pa = z[a].astype(float); pb = z[b].astype(float)
        m = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
        yy, gg = y[m], g[m]
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for met in FN:
            r = pair_boot(FN[met], yy, pa[m], pb[m], gg)
            emit(f"| Δ{met} | {r['point']:+.4f} | [{r['lo']:+.4f}, "
                 f"{r['hi']:+.4f}] | {r['n']} | "
                 f"{'✅' if r['excl'] else '❌ 含 0'} |")
            rows.append({"contrast": lab, "metric": met, "point": r["point"],
                         "ci_low": r["lo"], "ci_high": r["hi"],
                         "n_valid": r["n"], "ci_excludes_zero": r["excl"]})
        emit("")
        pd.DataFrame(rows).to_csv(OUT_C, index=False, encoding="utf-8")

    # ---- 分层 ----
    emit("## 7.1 Target-size 分层")
    emit("")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    lf = pd.read_csv(V71 / "phase1_coverage" / "block_L_rescued.csv",
                     low_memory=False)
    sub = can[can.row.isin(set(lf.pair_row.astype(int)))]
    nper = sub.groupby("pdb_id").size()
    emit("| 组 | 模型 | n_target | n_rows | target-mean r | slope | RMSE | C_W | "
         "Macro ρ | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    srows = []
    for grp, lo, hi in (("2-9", 2, 9), ("10-29", 10, 29), (">=30", 30, 10 ** 9)):
        tgs = set(nper[(nper >= lo) & (nper <= hi)].index.astype(str))
        sel = np.isin(g, list(tgs))
        if sel.sum() < 30:
            continue
        for nm in ("p_P0_raw_ranker", "p_P3_full_LSR", "p_P1_location_only"):
            if nm not in z.files:
                continue
            p = z[nm].astype(float)
            m2 = sel & np.isfinite(y) & np.isfinite(p)
            if m2.sum() < 20:
                continue
            yy, pp, gg = y[m2], p[m2], g[m2]
            emit(f"| {grp} | {nm} | {len(set(gg))} | {int(m2.sum())} | "
                 f"{m_tmean(yy, pp, gg):+.4f} | {m_slope(yy, pp, gg):+.4f} | "
                 f"{m_rmse(yy, pp, gg):.4f} | {m_cw(yy, pp, gg):+.4f} | "
                 f"{m_ms(yy, pp, gg):+.4f} | {m_eff(yy, pp, gg):+.4f} |")
            srows.append({"grp": grp, "model": nm, "n_target": len(set(gg)),
                          "n_rows": int(m2.sum()),
                          "target_mean_r": m_tmean(yy, pp, gg),
                          "atten_slope": m_slope(yy, pp, gg),
                          "rmse": m_rmse(yy, pp, gg), "C_W": m_cw(yy, pp, gg),
                          "macro_spearman": m_ms(yy, pp, gg),
                          "eff10": m_eff(yy, pp, gg)})
    srows_df = pd.DataFrame(srows)
    srows_df.to_csv(OUT_S, index=False, encoding="utf-8")
    emit("")
    # attenuation 是否改善
    a0 = srows_df[(srows_df.grp == ">=30") & (srows_df.model == "p_P0_raw_ranker")]
    a3 = srows_df[(srows_df.grp == ">=30") & (srows_df.model == "p_P3_full_LSR")]
    att_red = "N/A"
    if len(a0) and len(a3):
        s0 = float(a0.atten_slope.iloc[0]); s3 = float(a3.atten_slope.iloc[0])
        att_red = "YES" if abs(s3 - 1) < abs(s0 - 1) else "NO"
        emit(f"- `>=30` attenuation slope: P0 = **{s0:+.4f}** → "
             f"P3 = **{s3:+.4f}**；距 1 的距离 "
             f"{abs(s0 - 1):.4f} → {abs(s3 - 1):.4f}")
        emit(f"- **ATTENUATION_REDUCED: `{att_red}`**")
    emit("")

    # ---- weighting ----
    emit("## 7.2 Weighting 敏感性（诊断用，不用于选 winner）")
    emit("")
    emit("| weighting | P0 C_B | P3 C_B |")
    emit("|---|---:|---:|")
    wrows = []
    for wname in ("row", "equal", "sqrt", "capped30", "log"):
        out = {}
        for nm in ("p_P0_raw_ranker", "p_P3_full_LSR"):
            p = z[nm].astype(float)
            m = np.isfinite(y) & np.isfinite(p)
            yy, pp, gg = y[m], p[m], g[m]
            uq, By, Bp, cnt = _parts(yy, pp, gg)
            if wname == "row":
                w = cnt / cnt.sum()
            elif wname == "equal":
                w = np.ones(len(uq)) / len(uq)
            elif wname == "sqrt":
                w = np.sqrt(cnt); w /= w.sum()
            elif wname == "capped30":
                cc = np.minimum(cnt, 30.0); w = cc / cc.sum()
            else:
                cc = np.log1p(cnt); w = cc / cc.sum()
            den = float(yy.std()) * float(pp.std())
            out[nm] = float((w * (By - yy.mean()) * (Bp - pp.mean())).sum()) / den
        emit(f"| {wname} | {out['p_P0_raw_ranker']:+.4f} | "
             f"{out['p_P3_full_LSR']:+.4f} |")
        wrows.append({"weighting": wname, "P0_C_B": out["p_P0_raw_ranker"],
                      "P3_C_B": out["p_P3_full_LSR"]})
    emit("")
    pd.DataFrame(wrows).to_csv(OUT_W, index=False, encoding="utf-8")

    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_C}\n写出 {OUT_S}\n写出 {OUT_W}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
