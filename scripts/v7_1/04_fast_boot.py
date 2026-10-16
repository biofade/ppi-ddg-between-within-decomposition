# -*- coding: utf-8 -*-
"""scripts/v7_1/04_fast_boot.py -- V7.1 PHASE 7 收尾（剩余 contrasts + 分层 + weighting）.

只算 **P3−P0 / P3−P1 / P3−P2** 的 4 个核心指标 + 分层 + weighting，
1500 次（如实标注），每步写盘。
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
PH7 = V71 / "phase7_bootstrap"
NPZ = V71 / "phase5_lsr" / "lsr_predictions.npz"
N_BOOT = 1500
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def _parts(y, p, g):
    uq = np.unique(g)
    return (uq,
            np.array([y[g == c].mean() for c in uq]),
            np.array([p[g == c].mean() for c in uq]),
            np.array([(g == c).sum() for c in uq], float))


def m_r(y, p, g):
    return C.pearson(y, p)


def m_rmse(y, p, g):
    return float(np.sqrt(np.mean((y - p) ** 2)))


def m_tm(y, p, g):
    _u, By, Bp, _c = _parts(y, p, g)
    return C.pearson(By, Bp) if len(By) >= 3 else np.nan


def m_sl(y, p, g):
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


FN = {"Overall_r": m_r, "RMSE": m_rmse, "target_mean_r": m_tm,
      "atten_slope": m_sl, "C_W": m_cw}


def boot(fn, yy, pa, pb, gg, n=N_BOOT, seed=42):
    uq = np.unique(gg)
    cidx = [np.where(gg == u)[0] for u in uq]
    rs = np.random.RandomState(seed)
    pt = fn(yy, pa, gg) - fn(yy, pb, gg)
    vals = []
    for _ in range(n):
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
    emit("# V7.1 PHASE 7（收尾）— 剩余 contrasts + 分层 + weighting")
    emit("")
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**，单位 = complex（如实标注）")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    emit(f"- n = {len(y)} / {len(set(g))} targets")
    emit("")

    # 读已有 contrasts（P1/P2 部分），追加剩余
    prev = pd.read_csv(PH7 / "contrasts.csv") if (PH7 / "contrasts.csv").exists() \
        else pd.DataFrame()
    emit(f"- 已有 contrasts 行数: {len(prev)}")
    emit("")
    rows = prev.to_dict("records") if len(prev) else []
    for a, b, lab in (("p_P3_full_LSR", "p_P0_raw_ranker", "P3-P0 full LSR"),
                      ("p_P3_full_LSR", "p_P1_location_only", "P3-P1 scale incr"),
                      ("p_P3_full_LSR", "p_P2_scale_only", "P3-P2 location incr")):
        pa = z[a].astype(float); pb = z[b].astype(float)
        m = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
        yy, gg = y[m], g[m]
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for met in FN:
            r = boot(FN[met], yy, pa[m], pb[m], gg)
            emit(f"| Δ{met} | {r['point']:+.4f} | [{r['lo']:+.4f}, "
                 f"{r['hi']:+.4f}] | {r['n']} | "
                 f"{'✅' if r['excl'] else '❌ 含 0'} |")
            rows.append({"contrast": lab, "metric": met, "point": r["point"],
                         "ci_low": r["lo"], "ci_high": r["hi"],
                         "n_valid": r["n"], "ci_excludes_zero": r["excl"]})
        emit("")
        pd.DataFrame(rows).to_csv(PH7 / "contrasts.csv", index=False,
                                  encoding="utf-8")

    # ---- 分层 ----
    emit("## Target-size 分层")
    emit("")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    lf = pd.read_csv(V71 / "phase1_coverage" / "block_L_rescued.csv",
                     low_memory=False)
    sub = can[can.row.isin(set(lf.pair_row.astype(int)))]
    nper = sub.groupby("pdb_id").size()
    emit("| 组 | 模型 | n_target | n_rows | target-mean r | slope | RMSE | C_W |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|")
    srows = []
    for grp, lo, hi in (("2-9", 2, 9), ("10-29", 10, 29), (">=30", 30, 10 ** 9)):
        tgs = set(nper[(nper >= lo) & (nper <= hi)].index.astype(str))
        sel = np.isin(g, list(tgs))
        if sel.sum() < 30:
            continue
        for nm in ("p_P0_raw_ranker", "p_P1_location_only", "p_P3_full_LSR"):
            p = z[nm].astype(float)
            m2 = sel & np.isfinite(y) & np.isfinite(p)
            if m2.sum() < 20:
                continue
            yy, pp, gg = y[m2], p[m2], g[m2]
            emit(f"| {grp} | {nm} | {len(set(gg))} | {int(m2.sum())} | "
                 f"{m_tm(yy, pp, gg):+.4f} | {m_sl(yy, pp, gg):+.4f} | "
                 f"{m_rmse(yy, pp, gg):.4f} | {m_cw(yy, pp, gg):+.4f} |")
            srows.append({"grp": grp, "model": nm, "n_target": len(set(gg)),
                          "n_rows": int(m2.sum()),
                          "target_mean_r": m_tm(yy, pp, gg),
                          "atten_slope": m_sl(yy, pp, gg),
                          "rmse": m_rmse(yy, pp, gg), "C_W": m_cw(yy, pp, gg)})
    S = pd.DataFrame(srows)
    S.to_csv(PH7 / "stratification.csv", index=False, encoding="utf-8")
    emit("")
    a0 = S[(S.grp == ">=30") & (S.model == "p_P0_raw_ranker")]
    a3 = S[(S.grp == ">=30") & (S.model == "p_P3_full_LSR")]
    if len(a0) and len(a3):
        s0, s3 = float(a0.atten_slope.iloc[0]), float(a3.atten_slope.iloc[0])
        t0v, t3v = float(a0.target_mean_r.iloc[0]), float(a3.target_mean_r.iloc[0])
        emit(f"- `>=30` slope: P0 **{s0:+.4f}** → P3 **{s3:+.4f}**"
             f"（距 1：{abs(s0 - 1):.4f} → {abs(s3 - 1):.4f}）")
        emit(f"- `>=30` target-mean r: P0 **{t0v:+.4f}** → P3 **{t3v:+.4f}**")
        emit(f"- **ATTENUATION_REDUCED: "
             f"`{'YES' if abs(s3 - 1) < abs(s0 - 1) else 'NO'}`**")
    emit("")

    # ---- weighting ----
    emit("## Weighting 敏感性（诊断用）")
    emit("")
    emit("| weighting | P0 C_B | P3 C_B |")
    emit("|---|---:|---:|")
    wrows = []
    for wn in ("row", "equal", "sqrt", "capped30", "log"):
        out = {}
        for nm in ("p_P0_raw_ranker", "p_P3_full_LSR"):
            p = z[nm].astype(float)
            m = np.isfinite(y) & np.isfinite(p)
            yy, pp, gg = y[m], p[m], g[m]
            uq, By, Bp, cnt = _parts(yy, pp, gg)
            if wn == "row":
                w = cnt / cnt.sum()
            elif wn == "equal":
                w = np.ones(len(uq)) / len(uq)
            elif wn == "sqrt":
                w = np.sqrt(cnt); w /= w.sum()
            elif wn == "capped30":
                cc = np.minimum(cnt, 30.0); w = cc / cc.sum()
            else:
                cc = np.log1p(cnt); w = cc / cc.sum()
            den = float(yy.std()) * float(pp.std())
            out[nm] = float((w * (By - yy.mean()) * (Bp - pp.mean())).sum()) / den
        emit(f"| {wn} | {out['p_P0_raw_ranker']:+.4f} | "
             f"{out['p_P3_full_LSR']:+.4f} |")
        wrows.append({"weighting": wn, "P0_C_B": out["p_P0_raw_ranker"],
                      "P3_C_B": out["p_P3_full_LSR"]})
    pd.DataFrame(wrows).to_csv(PH7 / "weighting.csv", index=False, encoding="utf-8")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH7 / "bootstrap_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n写出 contrasts.csv / stratification.csv / weighting.csv / "
          "bootstrap_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
