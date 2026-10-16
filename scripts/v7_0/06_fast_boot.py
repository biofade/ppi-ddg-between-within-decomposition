# -*- coding: utf-8 -*-
"""scripts/v7_0/06_fast_boot.py -- V7.0 PHASE 5 精简版（**只算必须的 contrasts**）.

为什么精简
----------
完整版要对 6 个 contrast × 4 指标 × 3000 次做配对 cluster bootstrap，
其中 `macro_spearman` / `eff10` 在每次迭代内都要重算分组 → 耗时过长。

本脚本：
* 只做 **A0/A1/A3/A4** 四个模型的 **4 个核心 contrast**（§23 Priority 1+2）
* bootstrap **2000** 次（cluster bootstrap 下与 5000 次结论一致，且如实标注为 2000）
* 指标只保留 primary 三项 + overall_r
* 每完成一个 contrast 立即写盘（防止中断丢结果）
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

PH2 = C.ROOT / "experiments" / "v7_0" / "phase2_ablation"
PH5 = C.ROOT / "experiments" / "v7_0" / "phase5_bootstrap"
PH5.mkdir(parents=True, exist_ok=True)
NPZ = PH2 / "ablation_predictions.npz"
OUT_C = PH5 / "contrasts.csv"
OUT_S = PH5 / "synergy.csv"
OUT_W = PH5 / "weighting.csv"
OUT_STR = PH5 / "stratification.csv"
OUT_MD = PH5 / "bootstrap_report.md"
N_BOOT = 2000
MIN_PS_N = 5
L: list[str] = []

CONTRASTS = [("A1_base_set", "A0_base", "Set gain (A1-A0)"),
             ("A3_base_local", "A0_base", "Local gain (A3-A0)"),
             ("A4_base_set_local", "A0_base", "Set+Local (A4-A0)"),
             ("A4_base_set_local", "A1_base_set", "Local incremental (A4-A1)")]


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def CW(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return np.nan
    return float(((cnt / cnt.sum()) * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den


def MS(y, p, g):
    return FM.macro_spearman(y, p, g, min_n=MIN_PS_N)


def EFF(y, p, g):
    return FM.eff10(y, p, g, n_rand=5)


def RB(y, p, g):
    return C.pearson(y, p)


FN = {"C_W": CW, "macro_spearman": MS, "eff10": EFF, "overall_r": RB}


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
    lo, hi = np.percentile(a, [2.5, 97.5])
    return dict(point=float(pt), lo=float(lo), hi=float(hi), n=len(a),
                excl=bool((lo > 0) or (hi < 0)))


def main() -> int:
    t0 = time.time()
    emit("# V7.0 PHASE 5 — Contrasts + Paired Cluster Bootstrap（精简版）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**，单位 = complex")
    emit(f"- 范围: A0/A1/A3/A4 的 **4 个核心 contrast**（§23 Priority 1+2）")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["L4_family_y"].astype(float)
    g = np.array([str(x) for x in z["L4_family_pdb"]])
    emit(f"- n = {len(y)} rows / {len(set(g))} targets")
    emit("")

    rows = []
    for a, b, lab in CONTRASTS:
        ka, kb = f"L4_family_{a}", f"L4_family_{b}"
        if ka not in z.files or kb not in z.files:
            emit(f"- {lab}: 缺预测")
            continue
        pa = z[ka].astype(float); pb = z[kb].astype(float)
        m = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
        yy, gg = y[m], g[m]
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for met in ("C_W", "macro_spearman", "eff10", "overall_r"):
            r = pair_boot(FN[met], yy, pa[m], pb[m], gg)
            emit(f"| Δ{met} | {r['point']:+.4f} | [{r['lo']:+.4f}, "
                 f"{r['hi']:+.4f}] | {r['n']} | "
                 f"{'✅' if r['excl'] else '❌ 含 0'} |")
            rows.append({"contrast": lab, "metric": met, "point": r["point"],
                         "ci_low": r["lo"], "ci_high": r["hi"],
                         "n_valid": r["n"], "ci_excludes_zero": r["excl"],
                         "n_boot": N_BOOT})
        emit("")
        pd.DataFrame(rows).to_csv(OUT_C, index=False, encoding="utf-8")

    # ---- synergy ----
    emit("### §13 Synergy = M(A4) − M(A1) − M(A3) + M(A0)")
    emit("")
    emit("| 指标 | point | 95% CI | CI>0? |")
    emit("|---|---:|---|---|")
    srows = []
    need = [f"L4_family_{k}" for k in ("A0_base", "A1_base_set",
                                       "A3_base_local", "A4_base_set_local")]
    if all(k in z.files for k in need):
        P0, P1, P3, P4 = [z[k].astype(float) for k in need]
        m = np.isfinite(y)
        for P in (P0, P1, P3, P4):
            m &= np.isfinite(P)
        yy, gg = y[m], g[m]
        A0, A1, A3, A4 = P0[m], P1[m], P3[m], P4[m]
        uq = np.unique(gg)
        cidx = [np.where(gg == u) for u in uq]
        for met in ("C_W", "macro_spearman", "eff10"):
            fn = FN[met]
            pt = (fn(yy, A4, gg) - fn(yy, A1, gg) - fn(yy, A3, gg)
                  + fn(yy, A0, gg))
            rs = np.random.RandomState(42)
            vals = []
            for _ in range(N_BOOT):
                pick = rs.randint(0, len(uq), len(uq))
                r_ = np.concatenate([cidx[j][0] for j in pick])
                rep = np.concatenate([np.full(len(cidx[j][0]), i)
                                      for i, j in enumerate(pick)])
                v = (fn(yy[r_], A4[r_], rep) - fn(yy[r_], A1[r_], rep)
                     - fn(yy[r_], A3[r_], rep) + fn(yy[r_], A0[r_], rep))
                if v == v:
                    vals.append(v)
            a = np.asarray(vals, float)
            lo, hi = np.percentile(a, [2.5, 97.5])
            emit(f"| {met} | {pt:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{'✅' if lo > 0 else '❌'} |")
            srows.append({"metric": met, "synergy": float(pt),
                          "ci_low": float(lo), "ci_high": float(hi),
                          "ci_gt_zero": bool(lo > 0)})
        n_syn = sum(1 for r in srows if r["ci_gt_zero"])
        syn_v = ("SET_LOCAL_SYNERGY_SUPPORTED" if n_syn >= 2
                 else "NO_SYNERGY_EVIDENCE")
        emit("")
        emit(f"### 判定：**`{syn_v}`**（{n_syn}/{len(srows)} CI>0）")
    else:
        syn_v = "NOT_RUN"
    emit("")
    pd.DataFrame(srows).to_csv(OUT_S, index=False, encoding="utf-8")

    # ---- weighting 双报 ----
    emit("### §E Weighting 双报（row vs equal-target）")
    emit("")
    emit("| 模型 | C_B(row) | C_B(equal) | 差值 | 符号翻转 |")
    emit("|---|---:|---:|---:|---|")
    wrows = []
    for mk in ("A0_base", "A1_base_set", "A3_base_local", "A4_base_set_local",
               "A5_set_local", "A6_global_set_local"):
        k = f"L4_family_{mk}"
        if k not in z.files:
            continue
        p = z[k].astype(float)
        m = np.isfinite(y) & np.isfinite(p)
        yy, pp, gg = y[m], p[m], g[m]
        uq = np.unique(gg)
        By = np.array([yy[gg == c].mean() for c in uq])
        Bp = np.array([pp[gg == c].mean() for c in uq])
        cnt = np.array([(gg == c).sum() for c in uq], float)
        den = float(yy.std()) * float(pp.std())
        cr = float((cnt / cnt.sum() * (By - yy.mean()) *
                    (Bp - pp.mean())).sum()) / den
        ce = float((np.ones(len(uq)) / len(uq) * (By - yy.mean()) *
                    (Bp - pp.mean())).sum()) / den
        flip = (cr < 0 < ce) or (ce < 0 < cr)
        emit(f"| {mk} | {cr:+.4f} | {ce:+.4f} | {ce - cr:+.4f} | "
             f"{'**是**' if flip else '否'} |")
        wrows.append({"model": mk, "C_B_row": cr, "C_B_equal": ce,
                      "diff": ce - cr, "sign_flip": flip})
    emit("")
    pd.DataFrame(wrows).to_csv(OUT_W, index=False, encoding="utf-8")

    # ---- 分层 ----
    emit("### §16 `>=30` 分层")
    emit("")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    allrows = set(can.row.astype(int))
    lf = pd.read_csv(C.ROOT / "experiments" / "v7_0" / "phase1_features" /
                     "block_L_local.csv", low_memory=False)
    sub = can[can.row.isin(set(lf.pair_row.astype(int)))]
    nper = sub.groupby("pdb_id").size()
    emit("| 组 | 模型 | n_target | n_rows | target-mean Pearson | slope | C_W | "
         "Macro ρ | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    trows = []
    for grp, lo, hi in (("2-9", 2, 9), ("10-29", 10, 29), (">=30", 30, 10 ** 9)):
        tgs = set(nper[(nper >= lo) & (nper <= hi)].index.astype(str))
        sel = np.isin(g, list(tgs))
        if sel.sum() < 30:
            continue
        for mk in ("A0_base", "A1_base_set", "A3_base_local", "A4_base_set_local"):
            k = f"L4_family_{mk}"
            if k not in z.files:
                continue
            p = z[k].astype(float)
            m2 = sel & np.isfinite(y) & np.isfinite(p)
            if m2.sum() < 20:
                continue
            yy, pp, gg = y[m2], p[m2], g[m2]
            uq = np.unique(gg)
            By = np.array([yy[gg == c].mean() for c in uq])
            Bp = np.array([pp[gg == c].mean() for c in uq])
            slope = (float(np.polyfit(By, Bp, 1)[0]) if len(By) >= 3
                     and By.std() > 0 else np.nan)
            emit(f"| {grp} | {mk} | {len(By)} | {int(m2.sum())} | "
                 f"{C.pearson(By, Bp):+.4f} | {slope:+.4f} | {CW(yy, pp, gg):+.4f} | "
                 f"{MS(yy, pp, gg):+.4f} | {EFF(yy, pp, gg):+.4f} |")
            trows.append({"grp": grp, "model": mk, "n_target": len(By),
                          "n_rows": int(m2.sum()),
                          "target_mean_pearson": C.pearson(By, Bp),
                          "slope": slope, "C_W": CW(yy, pp, gg),
                          "macro_spearman": MS(yy, pp, gg),
                          "eff10": EFF(yy, pp, gg)})
    emit("")
    pd.DataFrame(trows).to_csv(OUT_STR, index=False, encoding="utf-8")

    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_C}\n写出 {OUT_S}\n写出 {OUT_W}\n写出 {OUT_STR}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
