# -*- coding: utf-8 -*-
"""scripts/v7_4/00_confirm_v73.py -- V7.4 PHASE 0（§3）.

**V7.3 confirmation gate**：用**向量化/优化实现**跑
`5000 / 5000` paired target-cluster bootstrap，只针对 5 个最关键的对比：
    F4-F1 / F4-F5 / F4-F2 / F8-F0 / F4-F4_SHUFFLED_ENV
+ **coverage-bias audit**（U_PHYS_STRICT included vs excluded）。

若主结论反转 → **停止 V7.4**（§3）。
"""

from __future__ import annotations

import argparse
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

V73 = C.ROOT / "experiments" / "v7_3"
V74 = C.ROOT / "experiments" / "v7_4"
PH0 = V74 / "phase0_confirmation"
PH0.mkdir(parents=True, exist_ok=True)
N_BOOT = 5000
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []

CONTRASTS = [
    ("F4_Base_PartnerPert", "F1_Base_Mut", "F4-F1   (核心) partner 扰动 vs mutation 性质"),
    ("F4_Base_PartnerPert", "F5_Base_IntraPert", "F4-F5   partner 扰动 vs intra 扰动"),
    ("F4_Base_PartnerPert", "F2_Base_Static", "F4-F2   扰动 vs 静态几何（§37）"),
    ("F8_Full", "F0_Base", "F8-F0   完整物理扰动 vs Base"),
]
SHUFFLE_CONTRAST = "F4-real - F4_SHUFFLED_ENV"


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def metrics(y, p, g):
    """一次算完 5 个指标（向量化，最小开销）。"""
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    cw = float((cnt / cnt.sum() * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den
    return (FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            FM.eff10(y, p, g, n_rand=10), cw,
            float(C.pearson(y, p)), float(np.sqrt(np.mean((y - p) ** 2))))


NAMES = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]


def boot_pair(y, pa, pb, g, n_boot, seed=42):
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    A0, B0 = metrics(y, pa, g), metrics(y, pb, g)
    pt = {n: A0[i] - B0[i] for i, n in enumerate(NAMES)}
    acc = {n: [] for n in NAMES}
    rs = np.random.RandomState(seed)
    for _ in range(n_boot):
        pick = rs.randint(0, len(uq), len(uq))
        r_ = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i)
                              for i, j in enumerate(pick)])
        try:
            A = metrics(y[r_], pa[r_], rep)
            Bx = metrics(y[r_], pb[r_], rep)
        except Exception:  # noqa: BLE001
            continue
        for i, n in enumerate(NAMES):
            v = A[i] - Bx[i]
            if v == v:
                acc[n].append(v)
    out = []
    for n in NAMES:
        arr = np.asarray(acc[n], float)
        if len(arr) < 20:
            lo = hi = np.nan; ex = False
        else:
            lo, hi = np.percentile(arr, [2.5, 97.5])
            ex = bool((lo > 0) or (hi < 0))
        out.append((n, pt[n], float(lo), float(hi), len(arr), ex))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--max-minutes", type=float, default=200.0)
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.4 PHASE 0 — V7.3 confirmation gate")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- bootstrap **{a.n_boot}/{a.n_boot}**（§45 要求 5000，本轮不再降级）")
    emit(f"- **主口径 = ensemble prediction**；cluster = canonical target")
    emit("")

    nz = V73 / "phase5_rankers" / "family_predictions.npz"
    z = np.load(nz, allow_pickle=True)
    y = z["y"].astype(float); g = np.array([str(x) for x in z["pdb"]])
    emit(f"- U_PHYS_STRICT: rows **{len(y)}** / targets **{len(set(g))}** / "
         f"n>={MIN_PS_N} **{sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)}**")
    emit("")

    def avg(b):
        ks = [f"{b}_s{s}" for s in SEEDS if f"{b}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0) if ks else None

    rows_out = []
    for a_, b_, lab in CONTRASTS:
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"- ⏱ 软时限，停止于 {lab}"); break
        PA, PB = avg(a_), avg(b_)
        m = np.isfinite(PA) & np.isfinite(PB)
        res = boot_pair(y[m], PA[m], PB[m], g[m], a.n_boot)
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for (n, pt, lo, hi, nv, ex) in res:
            emit(f"| Δ{n} | {pt:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {nv} | "
                 f"{'✅ 排除 0' if ex else '❌ 含 0'} |")
            rows_out.append({"contrast": lab, "metric": f"d_{n}", "point": pt,
                             "ci_low": lo, "ci_high": hi, "n_valid": nv,
                             "ci_excludes_zero": ex})
        emit(f"  ({(time.time() - t0) / 60:.1f} min)")
        emit("")
        pd.DataFrame(rows_out).to_csv(PH0 / "v73_bootstrap_5000.csv",
                                      index=False, encoding="utf-8")

    # ---------- F4 vs SHUFFLED ENV（需重跑对照）----------
    if (time.time() - t0) / 60 < a.max_minutes:
        emit(f"### {SHUFFLE_CONTRAST}")
        emit("")
        emit("> 需重算 `F4_SHUFFLED_ENV`（V7.3 未保存其逐行预测）。")
        emit("> 用 V7.3 完全相同的设计：**只在 outer-train 内**打乱环境列。")
        emit("")
        # 重算
        import importlib.util, json, gc
        _sp = importlib.util.spec_from_file_location(
            "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
        s4 = importlib.util.module_from_spec(_sp)
        sys.modules["s4"] = s4
        _sp.loader.exec_module(s4)
        PH = pd.read_parquet(C.ROOT / "features" / "v7_3" /
                             "physics_features.parquet")
        lab2 = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" /
                           "phase1_universe" / "canonical_target_labels.csv",
                           low_memory=False)
        PH = PH.merge(lab2.set_index("pair_row")[["family_id_canonical"]],
                      left_on="row_id", right_index=True, how="left")
        spec = json.loads((V73 / "phase5_rankers" / "family_spec.json")
                          .read_text(encoding="utf-8"))
        F4 = spec["F4_Base_PartnerPert"]
        env_cols = [c for c in F4 if (c.startswith("pt_") or "partner" in c
                                      or "_x_" in c)]
        env_idx = [F4.index(c) for c in env_cols]
        X4 = PH[F4].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
        fam = PH.family_id_canonical.astype(str).to_numpy()
        can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                          "matched_analysis_set.csv", low_memory=False).set_index("row")
        rows_ = PH.row_id.to_numpy(int)
        D4 = {"y": y, "pdb": g, "fam": fam,
              "seq40": can.loc[rows_, "seq40_cluster"].astype(str).to_numpy(),
              "obs": can.loc[rows_, "obs_key"].astype(str).to_numpy(),
              "state": can.loc[rows_, "state_id"].astype(str).to_numpy()}
        preds = {}
        for sd in SEEDS:
            outer = s4.make_folds(D4["fam"], sd, 5)
            yhat = np.full(len(y), np.nan)
            for tr, te in outer:
                Xtr = np.array(X4[tr], np.float32, copy=True)
                Xte = np.array(X4[te], np.float32, copy=True)
                rs = np.random.RandomState(sd)
                perm = rs.permutation(len(tr))
                Xtr[:, env_idx] = Xtr[np.ix_(perm, env_idx)]
                med = np.nanmedian(Xtr, axis=0)
                med = np.where(np.isfinite(med), med, 0.0)
                for A in (Xtr, Xte):
                    bad = ~np.isfinite(A)
                    if bad.any():
                        A[bad] = np.take(med, np.where(bad)[1])
                lgb = C.load_lgb()
                mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                           "n_estimators": 300,
                                           "random_state": sd})
                mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                yhat[te] = mdl.predict(Xte)
                del mdl; gc.collect()
            preds[sd] = yhat
        PSH = np.mean([preds[s] for s in SEEDS], axis=0)
        PF4 = avg("F4_Base_PartnerPert")
        m = np.isfinite(PF4) & np.isfinite(PSH)
        res = boot_pair(y[m], PF4[m], PSH[m], g[m], a.n_boot)
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        for (n, pt, lo, hi, nv, ex) in res:
            emit(f"| Δ{n} | {pt:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {nv} | "
                 f"{'✅ 排除 0' if ex else '❌ 含 0'} |")
            rows_out.append({"contrast": SHUFFLE_CONTRAST, "metric": f"d_{n}",
                             "point": pt, "ci_low": lo, "ci_high": hi,
                             "n_valid": nv, "ci_excludes_zero": ex})
        pd.DataFrame(rows_out).to_csv(PH0 / "v73_bootstrap_5000.csv",
                                      index=False, encoding="utf-8")
        emit("")

    # ---------- gate 判定 ----------
    emit("## GATE 判定（§3）")
    emit("")
    B = pd.DataFrame(rows_out)
    checks = []
    for cname, mname in (("F4-F1", "d_macro_spearman"),
                         ("F4-F5", "d_macro_spearman"),
                         ("F4-F2", "d_macro_spearman"),
                         ("F8-F0", "d_macro_spearman"),
                         (SHUFFLE_CONTRAST, "d_macro_spearman")):
        r = B[(B.contrast.str.startswith(cname)) & (B.metric == mname)]
        if not len(r):
            checks.append((cname, np.nan, False))
            continue
        r = r.iloc[0]
        checks.append((cname, r["point"], bool(r["ci_excludes_zero"])
                       and r["point"] > 0))
    emit("| 对比 | ΔMacro ρ | 方向为正且 CI 排除 0 |")
    emit("|---|---:|---|")
    for cname, pt, ok in checks:
        emit(f"| `{cname}` | {pt:+.4f} | {'✅' if ok else '❌'} |")
    emit("")
    passed = all(ok for _, _, ok in checks)
    emit(f"**V7.3 5000-bootstrap confirmation: "
         f"{'PASS ✅' if passed else 'FAIL ❌'}**")
    emit("")
    if not passed:
        emit("> ⇒ 依 §3：**停止 V7.4** 并报告。")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    # ---- handoff 块（§63）----
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH0 / "confirm_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    import json as _json
    _json.dump({"confirmation": "PASS" if passed else "FAIL",
                "checks": [{"contrast": c, "delta_macro": float(pt) if pt == pt else None,
                             "supported": bool(ok)} for c, pt, ok in checks],
                "n_boot": a.n_boot,
                "generated": time.strftime("%Y-%m-%dT%H:%M:%S")},
               open(PH0 / "confirm_gate.json", "w", encoding="utf-8"),
               indent=1, ensure_ascii=False)
    print(f"\n写出 v73_bootstrap_5000.csv 与 confirm_report.md")
    print(f">>> V7.3 confirmation = {'PASS' if passed else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
