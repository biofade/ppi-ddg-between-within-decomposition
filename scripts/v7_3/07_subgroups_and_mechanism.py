# -*- coding: utf-8 -*-
"""scripts/v7_3/07_subgroups_and_mechanism.py -- V7.3 PHASE 8+9.

PHASE 8（§46/§47，**只在主结果冻结之后**）
* 子群：single/multi、alanine、charge/hydrophobicity-changing、体积增减、
  界面状态（DIRECT/NEAR/DISTAL）、target 规模（2-9/10-29/>=30）
* 机制消融：逐块去掉 electrostatic / hydrophobic / packing / H-bond /
  salt-bridge / aromatic / burial / intra / partner / multi-site
  （用**重新训练**测 ΔMacro ρ，不靠 feature importance 推机制，§47）

PHASE 9（§49）：`y = mu^MU3 + z`，alpha=1，报 Pearson/RMSE/target-mean/
slope/row C_B/equal C_B，并断言**靶点内排序不变**。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

FEAT = C.ROOT / "features" / "v7_3"
V73 = C.ROOT / "experiments" / "v7_3"
PH8 = V73 / "phase8_mechanism"
PH9 = V73 / "phase9_calibration"
for _d in (PH8, PH9):
    _d.mkdir(parents=True, exist_ok=True)
TREES = 300
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    return {"macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "eff10": FM.eff10(y, p, g, n_rand=20),
            "C_W": float((cnt / cnt.sum() * np.array(
                [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den,
            "overall_r": float(C.pearson(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
            "macro_mean_pearson": float(C.pearson(By, Bp))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-minutes", type=float, default=55.0)
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.3 PHASE 8+9 — 机制/子群 + 标定")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    PH = pd.read_parquet(FEAT / "physics_features.parquet")
    FD = pd.read_csv(V73 / "phase3_features" / "feature_dict.csv")
    spec = json.loads((V73 / "phase5_rankers" / "family_spec.json")
                      .read_text(encoding="utf-8"))
    lab = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    PH = PH.merge(lab.set_index("pair_row")[["target_id_canonical",
                                            "family_id_canonical"]],
                  left_on="row_id", right_index=True, how="left")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False).set_index("row")
    d = C.build_pair_features(verbose=False)
    rows = PH.row_id.to_numpy(int)
    y = d["y"][rows].astype(float)
    g = PH.target_id_canonical.astype(str).to_numpy()
    fam = PH.family_id_canonical.astype(str).to_numpy()
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can.loc[rows, "state_id"].astype(str).to_numpy()}
    emit(f"- U_PHYS_STRICT: rows **{len(y)}** / targets **{len(set(g))}** / "
         f"families **{len(set(fam))}**")
    emit("")

    def run(cols, seeds=SEEDS):
        X = PH[cols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
        preds = {}
        for sd in seeds:
            outer = s4.make_folds(D4["fam"], sd, 5)
            yhat = np.full(len(y), np.nan)
            for tr, te in outer:
                Xtr = np.array(X[tr], np.float32, copy=True)
                Xte = np.array(X[te], np.float32, copy=True)
                med = np.nanmedian(Xtr, axis=0)
                med = np.where(np.isfinite(med), med, 0.0)
                for A in (Xtr, Xte):
                    bad = ~np.isfinite(A)
                    if bad.any():
                        A[bad] = np.take(med, np.where(bad)[1])
                lgb = C.load_lgb()
                mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                           "n_estimators": TREES,
                                           "random_state": sd})
                mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                yhat[te] = mdl.predict(Xte)
                del mdl; gc.collect()
            preds[sd] = yhat
        return preds

    F4 = spec["F4_Base_PartnerPert"]
    F8 = spec["F8_Full"]
    base = [c for c in PH.columns if c.startswith("base_")]
    P8 = run(F8)
    ens8 = np.mean([P8[s] for s in SEEDS], axis=0)
    ok8 = np.isfinite(ens8)
    emit(f"- `F8_Full` ensemble Macro ρ = "
         f"**{mset(y[ok8], ens8[ok8], g[ok8])['macro_spearman']:+.4f}**")
    emit("")

    # ================= PHASE 8.1 子群（§46）=================
    emit("## 8.1 子群分析（§46，基于 F8 冻结预测）")
    emit("")
    def _num(col, default=np.nan):
        if col in PH.columns:
            return pd.to_numeric(PH[col], errors="coerce").to_numpy(np.float64)
        return np.full(len(PH), default, np.float64)
    nsites = _num("n_mut_sites")
    nsite = nsites
    dchg = _num("d_charge__mean")
    dhph = _num("d_hydropathy__mean")
    dvol = _num("d_volume__mean")
    iface = _num("iface_state")
    tsize = pd.Series(g).map(pd.Series(g).value_counts()).to_numpy()

    SUBS = [
        ("single mutation", nsites <= 1),
        ("multi mutation", nsites >= 2),
        ("charge-changing (|Δq|>0.5)", np.abs(dchg) > 0.5),
        ("hydrophobicity-changing (|Δh|>1)", np.abs(dhph) > 1.0),
        ("large volume increase (ΔV>50)", dvol > 50),
        ("large volume decrease (ΔV<-50)", dvol < -50),
        ("DIRECT_INTERFACE", iface == 0),
        ("NEAR_INTERFACE", iface == 1),
        ("DISTAL", iface == 2),
        ("target 2-9", tsize <= 9),
        ("target 10-29", (tsize >= 10) & (tsize <= 29)),
        ("target >=30", tsize >= 30),
    ]
    emit("| 子群 | n_rows | n_target | Macro ρ | Eff@10 | C_W |")
    emit("|---|---:|---:|---:|---:|---:|")
    srows = []
    for name, m in SUBS:
        m = np.asarray(m) & ok8
        if m.sum() < 30 or len(set(g[m])) < 3:
            emit(f"| {name} | {int(m.sum())} | {len(set(g[m]))} | "
                 f"样本不足 | | |")
            continue
        mm = mset(y[m], ens8[m], g[m])
        emit(f"| {name} | {int(m.sum())} | {len(set(g[m]))} | "
             f"{mm['macro_spearman']:+.4f} | {mm['eff10']:+.4f} | "
             f"{mm['C_W']:+.4f} |")
        srows.append({"subgroup": name, "n_rows": int(m.sum()),
                      "n_target": len(set(g[m])), **mm})
    pd.DataFrame(srows).to_csv(PH8 / "subgroups.csv", index=False,
                               encoding="utf-8")
    emit("")

    # ================= PHASE 8.2 机制消融（§47）=================
    emit("## 8.2 机制消融（§47，逐块**重新训练**，非 feature importance）")
    emit("")
    BLOCKS = {
        "F_ELEC": ["F_ELEC"], "F_HYD": ["F_HYD"], "F_PACK": ["F_PACK"],
        "F_HBOND": ["F_HBOND"], "F_SALT": ["F_SALT"], "F_AROM": ["F_AROM"],
        "F_BURIAL": ["F_BURIAL", "F_IFACE"], "F_EPI_GEOM": ["F_EPI_GEOM"],
        "F_MUT": ["F_MUT"], "F_ENV_PARTNER": ["F_ENV_PARTNER"],
        "F_ENV_INTRA": ["F_ENV_INTRA"],
    }
    base_full = mset(y[ok8], ens8[ok8], g[ok8])
    emit("| 去掉的块 | ΔMacro ρ | ΔEff@10 | ΔC_W |")
    emit("|---|---:|---:|---:|")
    arows = []
    for bname, fams in BLOCKS.items():
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"| （软时限，停止） | | | |"); break
        drop = set(FD[FD.family.isin(fams)].feature.tolist())
        cols = [c for c in F8 if c not in drop]
        if len(cols) == len(F8):
            continue
        P = run(cols, seeds=[42])           # 单 seed 以控时（如实标注）
        e = P[42]
        ok = np.isfinite(e)
        mm = mset(y[ok], e[ok], g[ok])
        emit(f"| 去掉 `{bname}` | {base_full['macro_spearman'] - mm['macro_spearman']:+.4f} | "
             f"{base_full['eff10'] - mm['eff10']:+.4f} | "
             f"{base_full['C_W'] - mm['C_W']:+.4f} |")
        arows.append({"dropped": bname, "delta_macro": base_full["macro_spearman"]
                      - mm["macro_spearman"],
                      "delta_eff": base_full["eff10"] - mm["eff10"],
                      "delta_cw": base_full["C_W"] - mm["C_W"]})
        gc.collect()
    pd.DataFrame(arows).to_csv(PH8 / "mechanism_ablation.csv", index=False,
                               encoding="utf-8")
    emit("")
    emit("> ⚠️ 消融为**单 seed（42）**以控时；方向性结论，**非**精确量级（如实标注）。")
    emit("> ⚠️ 依 §47：**不从 LightGBM feature importance 推断机制**。")
    emit("")

    # ================= PHASE 9 标定（§49）=================
    emit("## 9. 标定外壳（§49）：`y = mu^MU3 + z`，alpha = 1")
    emit("")
    # MU3：target+set；用与 V7.1 相同来源（block_S 的 set 描述符）
    S = pd.read_csv(C.ROOT / "experiments" / "v7_0" / "phase1_features" /
                    "block_S_set.csv", low_memory=False)
    Scols = [c for c in S.columns if c != "pdb_id"]
    uq_t = np.unique(g)
    mu_true = np.array([y[g == c].mean() for c in uq_t])
    n_t = np.array([(g == c).sum() for c in uq_t])
    fam_t = pd.Series(fam).groupby(pd.Series(g)).first().reindex(uq_t).to_numpy()
    from sklearn.linear_model import RidgeCV
    Xt = S.set_index("pdb_id").reindex(uq_t)[Scols].apply(
        pd.to_numeric, errors="coerce").to_numpy(float)
    Xt = np.array(Xt, float, copy=True)
    med = np.nanmedian(Xt, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    bad = ~np.isfinite(Xt)
    if bad.any():
        Xt[bad] = np.take(med, np.where(bad)[1])
    sdv = Xt.std(0); sdv[sdv == 0] = 1
    Zt = (Xt - Xt.mean(0)) / sdv
    # family-aware CV 得 mu_hat
    from sklearn.model_selection import GroupKFold
    k = min(5, max(2, len(set(fam_t))))
    gkf = GroupKFold(n_splits=k)
    mu_hat = np.full(len(uq_t), np.nan)
    for tr, te in gkf.split(Zt, mu_true, groups=fam_t):
        m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(Zt[tr], mu_true[tr])
        mu_hat[te] = m.predict(Zt[te])
    mu_map = dict(zip(uq_t, mu_hat))
    # z 来自 F8 冻结预测（ensemble）→ 靶点内标准化
    z = np.zeros(len(y))
    for c in uq_t:
        m = g == c
        v = ens8[m]
        s = v.std(ddof=0)
        z[m] = 0.0 if s < 1e-6 else (v - v.mean()) / s
    yhat_cal = np.array([mu_map.get(c, np.nan) + zz for c, zz in zip(g, z)])
    okc = np.isfinite(yhat_cal)
    # 排序保持断言（§49）
    keep_rank = True
    for c in uq_t:
        m = np.where(g == c)[0]
        if len(m) < 2:
            continue
        a1 = np.argsort(np.argsort(ens8[m]))
        a2 = np.argsort(np.argsort(yhat_cal[m]))
        if not np.array_equal(a1, a2):
            keep_rank = False
            break
    By = np.array([y[g == c].mean() for c in uq_t])
    Bp = np.array([yhat_cal[g == c].mean() for c in uq_t])
    cnt = np.array([(g == c).sum() for c in uq_t], float)
    den = float(y.std()) * float(yhat_cal.std())
    cbr = float((cnt / cnt.sum() * (By - y.mean()) * (Bp - yhat_cal.mean())).sum()) / den
    cbe = float((np.ones(len(uq_t)) / len(uq_t) * (By - y.mean())
                 * (Bp - yhat_cal.mean())).sum()) / den
    sl = float(np.polyfit(By, Bp, 1)[0]) if By.std() > 0 else np.nan
    emit("| 指标 | raw z | **mu^MU3 + z** |")
    emit("|---|---:|---:|")
    emit(f"| Overall Pearson | {C.pearson(y[ok8], ens8[ok8]):+.4f} | "
         f"**{C.pearson(y[okc], yhat_cal[okc]):+.4f}** |")
    emit(f"| RMSE | {np.sqrt(np.mean((y[ok8]-ens8[ok8])**2)):.4f} | "
         f"**{np.sqrt(np.mean((y[okc]-yhat_cal[okc])**2)):.4f}** |")
    emit(f"| target-mean Pearson | {base_full['macro_mean_pearson']:+.4f} | "
         f"**{C.pearson(By, Bp):+.4f}** |")
    emit(f"| attenuation slope | — | **{sl:+.4f}** |")
    emit(f"| row C_B | — | **{cbr:+.4f}** |")
    emit(f"| equal-target C_B | — | **{cbe:+.4f}** |")
    emit("")
    emit(f"- **rank preservation：{'PASS ✅' if keep_rank else 'FAIL ❌'}**"
         f"（§49 要求 `within-target ranking unchanged`）")
    emit("")
    pd.DataFrame([{"metric": "overall_pearson_raw",
                   "value": C.pearson(y[ok8], ens8[ok8])},
                  {"metric": "overall_pearson_cal",
                   "value": C.pearson(y[okc], yhat_cal[okc])},
                  {"metric": "rmse_raw",
                   "value": float(np.sqrt(np.mean((y[ok8]-ens8[ok8])**2)))},
                  {"metric": "rmse_cal",
                   "value": float(np.sqrt(np.mean((y[okc]-yhat_cal[okc])**2)))},
                  {"metric": "target_mean_pearson", "value": float(C.pearson(By, Bp))},
                  {"metric": "atten_slope", "value": sl},
                  {"metric": "C_B_row", "value": cbr},
                  {"metric": "C_B_equal", "value": cbe},
                  {"metric": "rank_preserved", "value": float(keep_rank)}]
                 ).to_csv(PH9 / "calibration_summary.csv", index=False,
                          encoding="utf-8")
    np.savez(PH9 / "calibrated_predictions.npz",
             y=y.astype(np.float32), pdb=g.astype(str),
             z_raw=ens8.astype(np.float32), y_cal=yhat_cal.astype(np.float32))
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH8 / "phase8_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (PH9 / "phase9_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 subgroups.csv / mechanism_ablation.csv / calibration_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
