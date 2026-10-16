# -*- coding: utf-8 -*-
"""scripts/v7_0/04_bootstrap.py -- V7.0 PHASE 5: 对比 + 配对 cluster bootstrap.

§12 关键 contrasts（paired cluster bootstrap，5000/5000，单位 = complex）：
  A1−A0（Set gain） / A3−A0（Local gain） / A4−A0（Set+Local）
  A4−A1（Local incremental） / A4−A3（Set incremental） / A6−A4（Global incremental）

指标：ΔC_W / ΔMacroSpearman / ΔEfficiency@10 / ΔOverall Pearson

§13 Synergy = Metric(A4) − Metric(A1) − Metric(A3) + Metric(A0)
  对 C_W / MacroSpearman / Eff@10 做 bootstrap；CI > 0 → SET_LOCAL_SYNERGY_SUPPORTED

§16 另报 `>=30` 分层 + attenuation slope。
**所有 between 指标同时报 row-weighted 与 equal-target**（§E 强制）。
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
OUT_T = PH5 / "top_models_bootstrap.csv"
OUT_STR = PH5 / "stratification.csv"
OUT_MD = PH5 / "bootstrap_report.md"
N_BOOT = 3000
MIN_PS_N = 5
L: list[str] = []

CONTRASTS = [("A1_base_set", "A0_base", "Set gain (A1-A0)"),
             ("A3_base_local", "A0_base", "Local gain (A3-A0)"),
             ("A4_base_set_local", "A0_base", "Set+Local (A4-A0)"),
             ("A4_base_set_local", "A1_base_set", "Local incremental (A4-A1)"),
             ("A4_base_set_local", "A3_base_local", "Set incremental (A4-A3)"),
             ("A6_global_set_local", "A4_base_set_local",
              "Global incremental (A6-A4)")]
METRICS = [("C_W", None), ("macro_spearman", None), ("eff10", None),
           ("overall_r", None)]


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def _CB(y, p, g, equal=False):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    w = np.ones(len(uq)) / len(uq) if equal else cnt / cnt.sum()
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return np.nan
    return float((w * (By - y.mean()) * (Bp - p.mean())).sum()) / den


def _CW(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    w = cnt / cnt.sum()
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return np.nan
    return float((w * np.array([np.mean((y[g == c] - By[i]) *
                                       (p[g == c] - Bp[i]))
                               for i, c in enumerate(uq)])).sum()) / den


def m_CW(y, p, g):
    return _CW(y, p, g)


def m_MS(y, p, g):
    return FM.macro_spearman(y, p, g, min_n=MIN_PS_N)


def m_EFF(y, p, g):
    return FM.eff10(y, p, g, n_rand=5)


def m_R(y, p, g):
    return C.pearson(y, p)


FN = {"C_W": m_CW, "macro_spearman": m_MS, "eff10": m_EFF, "overall_r": m_R}


def boot_paired(fn, ya, pa, pb, g, n_boot=N_BOOT, seed=42):
    uq = np.unique(g)
    cidx = {u: np.where(g == u)[0] for u in uq}
    rs = np.random.RandomState(seed)
    pt = fn(ya, pa, g) - fn(ya, pb, g)
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uq), len(uq))
        rows = np.concatenate([cidx[uq[j]] for j in pick])
        rep = np.concatenate([np.full(len(cidx[uq[j]]), i)
                              for i, j in enumerate(pick)])
        v = fn(ya[rows], pa[rows], rep) - fn(ya[rows], pb[rows], rep)
        if v == v:
            vals.append(v)
    a = np.asarray(vals, float)
    if len(a) < 20:
        return dict(point=pt, lo=np.nan, hi=np.nan, n=len(a))
    lo, hi = np.percentile(a, [2.5, 97.5])
    return dict(point=float(pt), lo=float(lo), hi=float(hi), n=len(a),
                excl=bool((lo > 0) or (hi < 0)))


def main() -> int:
    t0 = time.time()
    emit("# V7.0 PHASE 5 — Contrasts + Paired Cluster Bootstrap")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- bootstrap: **{N_BOOT}** 次，单位 = **complex**")
    emit(f"- 数据: `{NPZ.name}`（PHASE 2，seed ensemble 预测）")
    emit("")
    if not NPZ.exists():
        emit("!! npz 缺失")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    z = np.load(NPZ, allow_pickle=True)
    y = z["L4_family_y"].astype(float)
    g = np.array([str(x) for x in z["L4_family_pdb"]])
    # 修复：原先 `k.split("_",1)[1]` 会得到 "family_A0_base" → models 为空。
    PRE = "L4_family_"
    models = sorted({k[len(PRE):] for k in z.files
                     if k.startswith(PRE)
                     and not k.endswith("_y") and not k.endswith("_pdb")
                     and not any(k.endswith(f"_s{s}") for s in (42, 43, 44))})
    models = [m for m in models if f"{PRE}{m}" in z.files]
    emit(f"- 可用模型: {models}")
    emit("")

    # ---------- contrasts ----------
    emit("## §12 关键 contrasts（Δ = 前者 − 后者）")
    emit("")
    rows = []
    for a, b, lab in CONTRASTS:
        ka, kb = f"L4_family_{a}", f"L4_family_{b}"
        if ka not in z.files or kb not in z.files:
            emit(f"- {lab}: 缺预测，跳过")
            continue
        pa, pb = z[ka].astype(float), z[kb].astype(float)
        m = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
        yy, gg = y[m], g[m]
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for met, _ in METRICS:
            r = boot_paired(FN[met], yy, pa[m], pb[m], gg)
            emit(f"| Δ{met} | {r['point']:+.4f} | "
                 f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['n']} | "
                 f"{'✅' if r.get('excl') else '❌ 含 0'} |")
            rows.append({"contrast": lab, "model_a": a, "model_b": b,
                         "metric": met, "point": r["point"],
                         "ci_low": r["lo"], "ci_high": r["hi"],
                         "n_valid": r["n"],
                         "ci_excludes_zero": r.get("excl", False)})
        emit("")
    pd.DataFrame(rows).to_csv(OUT_C, index=False, encoding="utf-8")

    # ---------- synergy ----------
    emit("## §13 Synergy = M(A4) − M(A1) − M(A3) + M(A0)")
    emit("")
    emit("| 指标 | point | 95% CI | 判定 |")
    emit("|---|---:|---|---|")
    need = ["L4_family_A0_base", "L4_family_A1_base_set",
            "L4_family_A3_base_local", "L4_family_A4_base_set_local"]
    srows = []
    if all(k in z.files for k in need):
        p0 = z[need[0]].astype(float); p1 = z[need[1]].astype(float)
        p3 = z[need[2]].astype(float); p4 = z[need[3]].astype(float)
        m = np.isfinite(y) & np.isfinite(p0) & np.isfinite(p1) & \
            np.isfinite(p3) & np.isfinite(p4)
        yy, gg = y[m], g[m]
        P0, P1, P3, P4 = p0[m], p1[m], p3[m], p4[m]
        for met, _ in METRICS:
            fn = FN[met]
            pt = (fn(yy, P4, gg) - fn(yy, P1, gg) - fn(yy, P3, gg)
                  + fn(yy, P0, gg))
            uq = np.unique(gg)
            cidx = {u: np.where(gg == u)[0] for u in uq}
            rs = np.random.RandomState(42)
            vals = []
            for _ in range(N_BOOT):
                pick = rs.randint(0, len(uq), len(uq))
                r_ = np.concatenate([cidx[uq[j]] for j in pick])
                rep = np.concatenate([np.full(len(cidx[uq[j]]), i)
                                      for i, j in enumerate(pick)])
                v = (fn(yy[r_], P4[r_], rep) - fn(yy[r_], P1[r_], rep)
                     - fn(yy[r_], P3[r_], rep) + fn(yy[r_], P0[r_], rep))
                if v == v:
                    vals.append(v)
            a = np.asarray(vals, float)
            lo, hi = np.percentile(a, [2.5, 97.5])
            syn = bool(lo > 0)
            emit(f"| {met} | {pt:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{'✅ CI>0' if syn else '❌'} |")
            srows.append({"metric": met, "synergy": float(pt),
                          "ci_low": float(lo), "ci_high": float(hi),
                          "ci_gt_zero": syn, "n_valid": len(a)})
        n_syn = sum(1 for r in srows if r["ci_gt_zero"])
        verdict = ("SET_LOCAL_SYNERGY_SUPPORTED" if n_syn >= 2
                   else "NO_SYNERGY_EVIDENCE")
        emit("")
        emit(f"### 判定：**`{verdict}`**（{n_syn}/{len(srows)} 个指标 CI>0）")
    else:
        verdict = "NOT_RUN"
        emit("| — | — | — | 缺 A0/A1/A3/A4 预测 |")
    emit("")
    pd.DataFrame(srows).to_csv(OUT_S, index=False, encoding="utf-8")

    # ---------- §16 分层 + attenuation ----------
    emit("## §16 `>=30` target 分层（baseline vs V7 候选）")
    emit("")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    rowmap = can.set_index("row")
    n_per_tg = rowmap.groupby("pdb_id").size()
    trows = []
    emit("| 组 | 模型 | n_target | n_rows | target-mean Pearson | slope | "
         "C_W | Macro ρ | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for grp, lo, hi in (("2-9", 2, 9), ("10-29", 10, 29), (">=30", 30, 10 ** 9)):
        tgs = set(n_per_tg[(n_per_tg >= lo) & (n_per_tg <= hi)].index.astype(str))
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
            By = np.array([yy[gg == c].mean() for c in np.unique(gg)])
            Bp = np.array([pp[gg == c].mean() for c in np.unique(gg)])
            slope = (float(np.polyfit(By, Bp, 1)[0]) if len(By) >= 3
                     and By.std() > 0 else np.nan)
            tm = C.pearson(By, Bp) if len(By) >= 3 else np.nan
            emit(f"| {grp} | {mk} | {len(By)} | {int(m2.sum())} | "
                 f"{tm:+.4f} | {slope:+.4f} | {m_CW(yy, pp, gg):+.4f} | "
                 f"{m_MS(yy, pp, gg):+.4f} | {m_EFF(yy, pp, gg):+.4f} |")
            trows.append({"grp": grp, "model": mk, "n_target": len(By),
                          "n_rows": int(m2.sum()),
                          "target_mean_pearson": tm, "slope": slope,
                          "C_W": m_CW(yy, pp, gg),
                          "macro_spearman": m_MS(yy, pp, gg),
                          "eff10": m_EFF(yy, pp, gg)})
    emit("")
    pd.DataFrame(trows).to_csv(OUT_STR, index=False, encoding="utf-8")

    # ---------- weighting 双报 ----------
    emit("## §E Weighting 双报（row-weighted 与 equal-target `C_B`）")
    emit("")
    emit("| 模型 | C_B(row) | C_B(equal) | 差值 | 符号是否翻转 |")
    emit("|---|---:|---:|---:|---|")
    wrows = []
    for mk in models:
        k = f"L4_family_{mk}"
        p = z[k].astype(float)
        m = np.isfinite(y) & np.isfinite(p)
        cr = _CB(y[m], p[m], g[m], False)
        ce = _CB(y[m], p[m], g[m], True)
        flip = (cr < 0 < ce) or (ce < 0 < cr)
        emit(f"| {mk} | {cr:+.4f} | {ce:+.4f} | {ce - cr:+.4f} | "
             f"{'**是**' if flip else '否'} |")
        wrows.append({"model": mk, "C_B_row": cr, "C_B_equal": ce,
                      "diff": ce - cr, "sign_flip": flip})
    emit("")
    pd.DataFrame(wrows).to_csv(OUT_T, index=False, encoding="utf-8")

    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_C}\n写出 {OUT_S}\n写出 {OUT_STR}\n写出 {OUT_T}")
    print(f"写出 {OUT_MD}")
    print(f"\n>>> synergy: {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
