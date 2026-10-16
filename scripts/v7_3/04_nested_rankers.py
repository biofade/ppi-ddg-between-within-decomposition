# -*- coding: utf-8 -*-
"""scripts/v7_3/04_nested_rankers.py -- V7.3 PHASE 5：F0–F9 冻结嵌套 family-CV.

特征族（§26，已冻结）
-------------------
| 族 | 构成 |
|---|---|
| `F0` | Base |
| `F1` | Base + mutation property deltas |
| `F2` | Base + STATIC_LOCAL_CONTROL（BLOCK L） |
| `F3` | Base + partner static environment |
| **`F4`** | Base + mutation-conditioned **partner** 扰动 ← 核心 |
| `F5` | Base + mutation-conditioned **intra** 扰动 |
| `F6` | F4 的 partner + intra 扰动 |
| `F7` | F6 + burial/SASA 条件化 |
| `F8` | F7 + 多突变位点耦合几何 |
| `F9` | **PHYSICS-ONLY**（无 Base） |

对比（§27）：`F4−F1`（核心）、`F8−F0`、`F2−F0`、`F1−F0`、`F4 vs F5`、`F9−F0`

纪律：train-only 中位数插补 / 无 PCA（物理特征保持低维可解释，§31）；
learner 统一 LightGBM；Ridge 仅作**可解释诊断**（§28），**不用于选优**。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
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
PH5 = V73 / "phase5_rankers"
PH5.mkdir(parents=True, exist_ok=True)
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


def fam_cols(keys, FD):
    """按 family 名取物理特征列。"""
    return FD[FD.family.isin(keys)].feature.tolist()


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    from scipy.stats import kendalltau
    kv = []
    for c in uq:
        m = g == c
        if m.sum() < 2:
            continue
        t = kendalltau(y[m], p[m]).correlation
        if t == t:
            kv.append(t)
    ken = float(np.mean(kv)) if kv else np.nan
    num = d2 = 0.0
    for c in uq:
        m = np.where(g == c)[0]
        if len(m) < 2:
            continue
        yy_ = y[m]; pp_ = p[m]
        dy = yy_[:, None] - yy_[None, :]
        dp = pp_[:, None] - pp_[None, :]
        iu = np.triu_indices(len(m), 1)
        dyv = dy[iu]; dpv = dp[iu]
        keep = np.abs(dyv) > 0.5
        if not keep.any():
            continue
        a_ = dyv[keep]; b_ = dpv[keep]
        d2 += a_.size
        num += float((np.sign(a_) == np.sign(b_))[b_ != 0].sum())
        num += 0.5 * float((b_ == 0).sum())
    pw = float(num / d2) if d2 > 0 else np.nan
    return {"macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "eff10": FM.eff10(y, p, g, n_rand=20),
            "C_W": float((cnt / cnt.sum() * np.array(
                [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den,
            "overall_r": float(C.pearson(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
            "kendall": ken, "pairwise05": pw}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fams", nargs="*", default=None)
    ap.add_argument("--max-minutes", type=float, default=70.0)
    ap.add_argument("--ridge", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.3 PHASE 5 — 冻结嵌套 family-CV（F0–F9）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- learner: **LightGBM**（trees={TREES}）；Ridge = 可解释诊断（不选优）")
    emit(f"- 预处理：**train-only** 中位数插补 + 标准化；**无 PCA**（§31）")
    emit("")

    PH = pd.read_parquet(FEAT / "physics_features.parquet")
    FD = pd.read_csv(V73 / "phase3_features" / "feature_dict.csv")
    lab = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    tmap = lab.set_index("pair_row")[["target_id_canonical",
                                     "family_id_canonical"]]
    PH = PH.merge(tmap, left_on="row_id", right_index=True, how="left")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False).set_index("row")
    d = C.build_pair_features(verbose=False)

    rows = PH.row_id.to_numpy(int)
    y = d["y"][rows].astype(float)
    g = PH.target_id_canonical.astype(str).to_numpy()
    fam = PH.family_id_canonical.astype(str).to_numpy()
    emit(f"- **U_PHYS_STRICT: rows {len(y)} / targets {len(set(g))} / "
         f"families {len(set(fam))}**")
    emit(f"- n>={MIN_PS_N} target："
         f"**{sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)}**")
    emit("")

    base = [f"base_{i}" for i in range(188)]
    bl = [c for c in PH.columns if c.startswith("bl_")]
    f_mut = fam_cols(["F_MUT"], FD)
    f_envp = fam_cols(["F_ENV_PARTNER"], FD)
    f_envi = fam_cols(["F_ENV_INTRA"], FD)
    f_part = fam_cols(["F_ELEC", "F_HYD", "F_PACK", "F_HBOND", "F_SALT",
                       "F_AROM"], FD)
    # partner-only / intra-only 扰动列
    part_only = [c for c in f_part if "partner" in c or "_x_" in c]
    intra_only = [c for c in f_part if "intra" in c]
    f_burial = fam_cols(["F_BURIAL", "F_IFACE"], FD)
    f_epi = fam_cols(["F_EPI_GEOM"], FD)

    FAMILIES = {
        "F0_Base": base,
        "F1_Base_Mut": base + f_mut,
        "F2_Base_Static": base + bl,
        "F3_Base_EnvPartner": base + f_envp,
        "F4_Base_PartnerPert": base + f_mut + part_only,
        "F5_Base_IntraPert": base + f_mut + intra_only,
        "F6_Partner_Intra": base + f_mut + part_only + intra_only,
        "F7_F6_Burial": base + f_mut + part_only + intra_only + f_burial,
        "F8_Full": base + f_mut + part_only + intra_only + f_burial + f_epi,
        "F9_PhysicsOnly": f_mut + part_only + intra_only + f_burial + f_epi,
    }
    if a.fams:
        FAMILIES = {k: v for k, v in FAMILIES.items() if k in a.fams}
    emit("## 5.1 冻结特征族（§26）")
    emit("")
    emit("| 族 | n_features | 构成 |")
    emit("|---|---:|---|")
    for k, v in FAMILIES.items():
        emit(f"| `{k}` | {len(v)} | {'+'.join(v[:0]) or ''}"
             f"{'含 Base' if k != 'F9_PhysicsOnly' else '**无 Base**'} |")
    emit("")
    emit("> 出现在多个族中的列已去重（见 `family_spec.json`）")
    emit("")
    import json
    json.dump({k: v for k, v in FAMILIES.items()},
              open(PH5 / "family_spec.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)

    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can.loc[rows, "state_id"].astype(str).to_numpy()}
    res, preds = [], {}
    leak = {"outer": 0, "inner": 0}
    emit("## 5.2 嵌套 family-CV")
    emit("")
    for sd in SEEDS:
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"- ⏱ 软时限，停止于 seed {sd}"); break
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk, cols in FAMILIES.items():
            try:
                X = PH[cols].apply(pd.to_numeric, errors="coerce").to_numpy(
                    np.float32)
                yhat = np.full(len(y), np.nan)
                for tr, te in outer:
                    if set(fam[tr]) & set(fam[te]):
                        leak["outer"] += 1
                    # train-only 中位数插补（§32）
                    med = np.nanmedian(X[tr], axis=0)
                    med = np.where(np.isfinite(med), med, 0.0)
                    Xtr = np.array(X[tr], np.float32, copy=True)
                    Xte = np.array(X[te], np.float32, copy=True)
                    for A, B in ((Xtr, None), (Xte, None)):
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
                ok = np.isfinite(yhat)
                m = mset(y[ok], yhat[ok], g[ok])
                res.append({"seed": sd, "family": mk, "n_features": len(cols),
                            "n_rows": int(ok.sum()),
                            "n_target": len(set(g[ok])), **m})
                preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
                emit(f"| {mk} | s{sd} | Macro ρ={m['macro_spearman']:+.4f} | "
                     f"Eff@10={m['eff10']:+.4f} | C_W={m['C_W']:+.4f} | "
                     f"Kendall={m['kendall']:+.4f} | r={m['overall_r']:+.4f} |")
                pd.DataFrame(res).to_csv(PH5 / "family_results.csv", index=False,
                                         encoding="utf-8")
            except Exception as ex:  # noqa: BLE001
                emit(f"| {mk} | s{sd} | FAILED {type(ex).__name__}: {ex} |")
                gc.collect()
    R = pd.DataFrame(res)
    preds["y"] = y.astype(np.float32); preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str); preds["row_id"] = rows.astype(np.int64)
    np.savez(PH5 / "family_predictions.npz", **preds)
    emit("")
    emit(f"- **泄漏断言**：outer violation = **{leak['outer']}**")
    emit(f"- **train-only 插补/标准化：PASS**；**无 PCA**")
    emit("")
    if len(R):
        cols2 = ["macro_spearman", "eff10", "C_W", "kendall", "pairwise05",
                 "overall_r", "rmse"]
        sm = R.groupby("family")[cols2].mean()
        emit("### 3-seed 均值（敏感性口径）")
        emit("")
        emit("| 族 | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Overall r |")
        emit("|---|---:|---:|---:|---:|---:|---:|")
        for k, r in sm.iterrows():
            emit(f"| `{k}` | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
                 f"{r['C_W']:+.4f} | {r['kendall']:+.4f} | "
                 f"{r['pairwise05']:+.4f} | {r['overall_r']:+.4f} |")
        sm.to_csv(PH5 / "family_seedmean.csv", encoding="utf-8")
        # ensemble 口径
        emit("")
        emit("### ensemble prediction 口径（主口径）")
        emit("")
        emit("| 族 | Macro ρ | Eff@10 | C_W |")
        emit("|---|---:|---:|---:|")
        ens = {}
        for mk in FAMILIES:
            ks = [f"{mk}_s{s}" for s in SEEDS if f"{mk}_s{s}" in preds]
            if not ks:
                continue
            Pm = np.mean([preds[k].astype(float) for k in ks], axis=0)
            ok = np.isfinite(Pm)
            ens[mk] = mset(y[ok], Pm[ok], g[ok])
            emit(f"| `{mk}` | {ens[mk]['macro_spearman']:+.4f} | "
                 f"{ens[mk]['eff10']:+.4f} | {ens[mk]['C_W']:+.4f} |")
        pd.DataFrame([{"family": k, **v} for k, v in ens.items()]).to_csv(
            PH5 / "family_ensemble.csv", index=False, encoding="utf-8")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH5 / "phase5_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 family_results.csv / family_ensemble.csv / phase5_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
