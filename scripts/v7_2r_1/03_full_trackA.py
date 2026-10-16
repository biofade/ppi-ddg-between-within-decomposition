# -*- coding: utf-8 -*-
"""scripts/v7_2r_1/03_full_trackA.py -- V7.2R-1 PHASE 1(终)/2 — 完整 Track A 评估.

1. 冻结 **Track A universe manifest**（覆盖率、targets、families、排除原因）
2. 在冻结 universe 上跑 **R0–R3 × 3 seed** 的**冻结 V7.1.1 嵌套 family-CV**
3. **5000/5000** paired target-cluster bootstrap（R1−R0 / R2−R0 / R3−R0）
4. 主结果用 **ensemble prediction 口径**（与 bootstrap 一致），
   **同时**报告 seed-mean of metric 作为敏感性

输出目录：`experiments/v7_2r_1/`
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import re
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
_spec2 = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec2)
sys.modules["v7f"] = v7f
_spec2.loader.exec_module(v7f)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

FEAT = C.ROOT / "features" / "v7_2"
PER = FEAT / "per_target"
V711 = C.ROOT / "experiments" / "v7_1_1"
V70 = C.ROOT / "experiments" / "v7_0"
R1 = C.ROOT / "experiments" / "v7_2r_1"
PH1, PH2 = R1 / "phase1_coverage", R1 / "phase2_full_eval"
PH6 = C.ROOT / "experiments" / "v7_2" / "phase7_bootstrap" / "v7_2r_1"
for d in (PH1, PH2, PH6):
    d.mkdir(parents=True, exist_ok=True)

TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
EPS = 1e-6
SEEDS = [42, 43, 44]
N_BOOT = 5000
L: list[str] = []

RANKERS = {
    "R0_Base": ["base"],
    "R1_Base_F0": ["base", "F0"],
    "R2_Base_Residue": ["base", "F1"],
    "R3_Base_dh": ["base", "F2"],
}


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def std_within_target(s, g, eps=EPS):
    z = np.zeros(len(s), float); coll = np.zeros(len(s), bool)
    for c in np.unique(g):
        m = np.where(g == c)[0]; v = s[m]; sd = v.std(ddof=0)
        if sd < eps:
            coll[m] = True
        else:
            z[m] = (v - v.mean()) / sd
    return z, coll


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    return {
        "macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        "eff10": FM.eff10(y, p, g, n_rand=10),
        "C_W": float((cnt / cnt.sum() * np.array(
            [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
             for i, c in enumerate(uq)])).sum()) / den,
        "overall_r": float(C.pearson(y, p)),
        "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=float, default=0.0)
    ap.add_argument("--bootstrap", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.2R-1 — 完整 Track A 评估（冻结嵌套 family-CV）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- learner: **LightGBM**（trees={TREES}）；PCA={PCA_NEW}（**train-only**）")
    emit(f"- 协议: V7.1.1 冻结 `OUTER_FOLDS_V7_2` + inner family CV")
    emit(f"- 主结果口径: **ensemble prediction**；敏感性: seed-mean of metric")
    emit("")

    # ---------- 载入嵌入 + 冻结 universe ----------
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]

    hw, hm = {}, {}
    for p in sorted(PER.glob("*.npz")):
        try:
            z = np.load(p)
            hw[p.stem] = z["h_wt"]; hm[p.stem] = z["h_mut"]
        except Exception:  # noqa: BLE001
            continue
    hit = np.zeros(len(lab), bool)
    pos = pd.Series(np.arange(len(lab)), index=lab.pair_row.to_numpy(int))
    for t in hw:
        z = np.load(PER / f"{t}.npz")
        ii = pos.reindex(z["row_id"]).to_numpy()
        ok = np.isfinite(ii)
        if ok.any():
            hit[ii[ok].astype(int)] = True
    cov = hit.mean()
    emit("## PHASE 1 — Track A universe manifest（**冻结**）")
    emit("")
    sub = lab[hit]
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit(f"| U_CONTRACT rows | {len(lab)} |")
    emit(f"| **Track A covered rows** | **{int(hit.sum())}** |")
    emit(f"| **coverage** | **{cov * 100:.1f}%** |")
    emit(f"| targets | **{sub.target_id_canonical.nunique()}** / "
         f"{lab.target_id_canonical.nunique()} |")
    emit(f"| families | **{sub.family_id_canonical.nunique()}** / "
         f"{lab.family_id_canonical.nunique()} |")
    n5 = (sub.groupby("target_id_canonical").size() >= MIN_PS_N).sum()
    emit(f"| n>={MIN_PS_N} targets（Macro ρ 资格） | **{n5}** |")
    emit("")
    miss = lab[~hit]
    if len(miss):
        emit(f"- **未覆盖 {len(miss)} 行**，涉及 "
             f"{miss.target_id_canonical.nunique()} 个 target。")
        emit(f"- **排除原因**：这些复合物在 `mutation_sites` 中**没有任何位点**"
             f"同时满足 `chain ∈ PDB`、`seq_index0` 在链长内、"
             f"且 `mapped_WT == expected_WT`（§19：绝不接受 WT mismatch）。")
        emit(f"  即覆盖率受**WS 校验约束**限制，**未因提高覆盖率而放宽校验**。")
    emit("")
    rows = sub.pair_row.to_numpy(int)
    # 组装
    idx = {int(r): i for i, r in enumerate(rows)}
    HW = np.zeros((len(rows), 320), np.float32)
    HM = np.zeros((len(rows), 320), np.float32)
    for t, arr in hw.items():
        z = np.load(PER / f"{t}.npz")
        for i, r in enumerate(z["row_id"]):
            j = idx.get(int(r))
            if j is None:
                continue
            HW[j] = arr[i].astype(np.float32)
            HM[j] = hm[t][i].astype(np.float32)
    y = d["y"][rows].astype(float)
    g = can_i.loc[rows, "pdb_id"].astype(str).to_numpy()
    fam = can_i.loc[rows, "family_label"].astype(str).to_numpy()
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    # F0
    fwd = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = fwd[fwd.direction == "forward"].reset_index(drop=True)
    f0 = []
    for r in rows:
        s = str(fwd["mutation"].iloc[r]).upper(); muts = []
        for tok in re.split(r"[;,]", s):
            tok = tok.strip()
            if not tok:
                continue
            body = tok.split("_", 1)[-1] if "_" in tok else tok
            for m in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])", body):
                muts.append((m.group(1), m.group(3)))
        if not muts:
            muts = [("X", "X")]
        vs = [v7f.mut_identity(a, b) for a, b in muts]
        ks = list(vs[0].keys())
        f0.append([float(np.mean([v[k] for v in vs])) for k in ks])
    F0 = np.array(f0, np.float32)
    F1 = np.hstack([HW, HM])
    dH = HM - HW
    F2 = np.hstack([dH, np.linalg.norm(dH, axis=1, keepdims=True),
                    (np.sum(HW * HM, axis=1, keepdims=True) /
                     (np.linalg.norm(HW, axis=1, keepdims=True) *
                      np.linalg.norm(HM, axis=1, keepdims=True) + 1e-9))])
    FIT = {"base": base, "F0": F0, "F1": F1, "F2": F2}
    emit(f"- 特征块维度：base={base.shape[1]}, F0={F0.shape[1]}, "
         f"F1={F1.shape[1]}, F2={F2.shape[1]}")
    emit(f"- **universe: rows {len(y)} / targets {len(set(g))} / "
         f"families {len(set(fam))}**")
    emit("")
    pd.DataFrame([{"rows": int(hit.sum()), "targets": sub.target_id_canonical.nunique(),
                   "families": sub.family_id_canonical.nunique(),
                   "coverage": float(cov), "n_ge5": int(n5),
                   "contract_rows": len(lab)}]).to_csv(
        PH1 / "trackA_universe_manifest.csv", index=False, encoding="utf-8")
    lab.assign(covered=hit).to_csv(PH1 / "trackA_row_coverage.csv", index=False,
                                   encoding="utf-8")

    if cov < args.min_cov:
        emit(f"!! 覆盖率 {cov * 100:.1f}% < 要求 {args.min_cov * 100:.0f}% → 停止")
        (PH2 / "full_trackA_report.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
        return 1

    # ---------- PHASE 2：嵌套 CV ----------
    emit("## PHASE 2 — 完整 R0–R3 嵌套 family-CV")
    emit("")
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rows, "state_id"].astype(str).to_numpy()}
    res, preds = [], {}
    for sd in SEEDS:
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk, blocks in RANKERS.items():
            try:
                yhat = np.full(len(y), np.nan)
                for tr, te in outer:
                    assert not (set(fam[tr]) & set(fam[te])), "outer 泄漏"
                    mtr, mte = [], []
                    for b in blocks:
                        X = FIT[b]
                        if b == "base":
                            mtr.append(X[tr]); mte.append(X[te]); continue
                        from sklearn.decomposition import PCA
                        mu = np.nanmean(X[tr], axis=0)
                        sv = np.nanstd(X[tr], axis=0); sv[sv == 0] = 1
                        Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
                        k = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                        pc = PCA(n_components=k, random_state=sd).fit(Ztr)
                        mtr.append(pc.transform(Ztr).astype(np.float32))
                        mte.append(pc.transform(Zte).astype(np.float32))
                    Xtr = np.hstack(mtr).astype(np.float32)
                    Xte = np.hstack(mte).astype(np.float32)
                    lgb = C.load_lgb()
                    mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                               "n_estimators": TREES,
                                               "random_state": sd})
                    mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                    yhat[te] = mdl.predict(Xte)
                    del mdl; gc.collect()
                ok = np.isfinite(yhat)
                m = mset(y[ok], yhat[ok], g[ok])
                res.append({"seed": sd, "ranker": mk, "n_rows": int(ok.sum()),
                            "n_target": len(set(g[ok])), **m})
                preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
                emit(f"| {mk} | s{sd} | Macro ρ={m['macro_spearman']:+.4f} | "
                     f"Eff@10={m['eff10']:+.4f} | C_W={m['C_W']:+.4f} | "
                     f"r={m['overall_r']:+.4f} |")
            except Exception as ex:  # noqa: BLE001
                emit(f"| {mk} | s{sd} | FAILED {type(ex).__name__}: {ex} |")
                gc.collect()
    R = pd.DataFrame(res)
    R.to_csv(PH2 / "full_nested_cv_results.csv", index=False, encoding="utf-8")
    preds["y"] = y.astype(np.float32); preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str); preds["row_id"] = rows.astype(np.int64)
    np.savez(PH2 / "full_ranker_predictions.npz", **preds)
    emit("")

    # ---------- 两种口径 ----------
    emit("## 2.1 两种聚合口径（**必须同时报告**）")
    emit("")
    sm = R.groupby("ranker")[["macro_spearman", "eff10", "C_W", "overall_r"]].mean()
    ens = {}
    for mk in RANKERS:
        ks = [f"{mk}_s{s}" for s in SEEDS if f"{mk}_s{s}" in preds]
        if not ks:
            continue
        P = np.mean([preds[k].astype(float) for k in ks], axis=0)
        o = np.isfinite(P)
        ens[mk] = mset(y[o], P[o], g[o])
    emit("| Ranker | 口径 | Macro ρ | Eff@10 | C_W | Overall r |")
    emit("|---|---|---:|---:|---:|---:|")
    for mk in RANKERS:
        if mk in sm.index:
            s = sm.loc[mk]
            emit(f"| {mk} | seed-mean of metric | {s['macro_spearman']:+.4f} | "
                 f"{s['eff10']:+.4f} | {s['C_W']:+.4f} | {s['overall_r']:+.4f} |")
        if mk in ens:
            e = ens[mk]
            emit(f"| {mk} | **ensemble prediction** | {e['macro_spearman']:+.4f} | "
                 f"{e['eff10']:+.4f} | {e['C_W']:+.4f} | {e['overall_r']:+.4f} |")
    emit("")
    pd.DataFrame([{"ranker": k, "pathway": "ensemble", **v} for k, v in ens.items()]).to_csv(
        PH2 / "ensemble_metrics.csv", index=False, encoding="utf-8")
    sm.to_csv(PH2 / "seedmean_metrics.csv", encoding="utf-8")

    # ---------- PHASE：bootstrap ----------
    if args.bootstrap:
        emit("## PHASE 3 — 5000/5000 paired target-cluster bootstrap")
        emit("")
        uq = np.unique(g)
        cidx = [np.where(g == c)[0] for c in uq]
        for a, b in (("R3_Base_dh", "R0_Base"), ("R2_Base_Residue", "R0_Base"),
                     ("R1_Base_F0", "R0_Base")):
            ka = [f"{a}_s{s}" for s in SEEDS if f"{a}_s{s}" in preds]
            kb = [f"{b}_s{s}" for s in SEEDS if f"{b}_s{s}" in preds]
            if not ka or not kb:
                emit(f"### {a} − {b} → NOT_RUN"); emit(""); continue
            PA = np.mean([preds[k].astype(float) for k in ka], axis=0)
            PB = np.mean([preds[k].astype(float) for k in kb], axis=0)
            rs = np.random.RandomState(42)
            pt = {k: mset(y, PA, g)[k] - mset(y, PB, g)[k]
                  for k in ("macro_spearman", "eff10", "C_W", "overall_r", "rmse")}
            acc = {k: [] for k in pt}
            for _ in range(N_BOOT):
                pick = rs.randint(0, len(uq), len(uq))
                r_ = np.concatenate([cidx[j] for j in pick])
                rep = np.concatenate([np.full(len(cidx[j]), i)
                                      for i, j in enumerate(pick)])
                try:
                    ma = mset(y[r_], PA[r_], rep); mb = mset(y[r_], PB[r_], rep)
                except Exception:  # noqa: BLE001
                    continue
                for k in pt:
                    v = ma[k] - mb[k]
                    if v == v:
                        acc[k].append(v)
            emit("")
            emit(f"### {a} − {b}")
            emit("")
            emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
            emit("|---|---:|---|---:|---|")
            rr = []
            for k in pt:
                arr = np.asarray(acc[k], float)
                if len(arr) < 20:
                    lo = hi = np.nan; ex = False
                else:
                    lo, hi = np.percentile(arr, [2.5, 97.5])
                    ex = bool((lo > 0) or (hi < 0))
                emit(f"| Δ{k} | {pt[k]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                     f"{len(arr)} | {'✅' if ex else '❌ 含 0'} |")
                rr.append({"contrast": f"{a}-{b}", "metric": f"d_{k}",
                           "point": pt[k], "ci_low": float(lo),
                           "ci_high": float(hi), "n_valid": len(arr),
                           "ci_excludes_zero": ex})
            ok = (rr[0].get("ci_excludes_zero") or False)
            dmap = {r["metric"]: r for r in rr}
            sup = (dmap["d_macro_spearman"]["point"] > 0
                   and dmap["d_macro_spearman"]["ci_excludes_zero"]
                   and (dmap["d_eff10"]["point"] > 0 or dmap["d_C_W"]["point"] > 0)
                   and dmap["d_rmse"]["ci_high"] < 0.5)
            emit("")
            emit(f"- **§23 成功判据（冻结）：{'SUPPORTED ✅' if sup else 'NOT SUPPORTED ❌'}**")
            emit("")
            pd.DataFrame(rr).to_csv(PH2 / f"bootstrap_{a}_vs_{b}.csv", index=False,
                                    encoding="utf-8")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH2 / "full_trackA_report.md").write_text("\n".join(L) + "\n",
                                               encoding="utf-8")
    print(f"\n写出 PHASE1/2 文件（universe rows={int(hit.sum())}, "
          f"coverage={cov * 100:.1f}%）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
