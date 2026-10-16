# -*- coding: utf-8 -*-
"""scripts/v8/38_boundary_bootstrap_clean.py -- 用**规范实现**（v8metrics）做边界 bootstrap。

## 背景

`36_boundary_bootstrap.py` 用内联一/二阶量算 `overall_pearson`，实测**有 bug**
（给 PERMISSIVE 组 −0.5579，而规范实现 `v8metrics` 给 **+0.8088**，逐位复核见
`37_diagnose_discrepancy.py`：pool 文件与 v8metrics 重算 max |Δ| = 5.9e-3）。

⇒ `36` 的 bootstrap 结果**作废**。本脚本用 `v8metrics.compute_all` 重做。

## 设计

* 模型池：`pool_composition_models.csv` 的 PERMISSIVE(37) + STRICT(133)
* 每次迭代：重采样 target -> 重建行子集 -> 对每个模型调用 `compute_all`
* 为控制成本，默认 2000 次（并**如实报告 n_valid**）；`--n-boot` 可调
* 检验量：Δ = ρ_PERMISSIVE(Overall r, Eff@10) − ρ_STRICT(...)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
EXT = OUT / "external"
SEED = 42
MIN_N = 5
K = 10


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--max-minutes", type=float, default=40.0)
    a = ap.parse_args()

    H = pd.read_csv(EXT / "hermes_matched_final.csv", low_memory=False)
    z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
                allow_pickle=True)
    Y = np.asarray(z["y"], float)
    lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    H = H.sort_values("pair_row").reset_index(drop=True)
    TGT = H.target_id_canonical.astype(str).to_numpy()
    FAM = H.family_id_canonical.astype(str).to_numpy()
    assert len(H) == len(Y) == len(lab)

    POOL = pd.read_csv(EXT / "pool_composition_models.csv", low_memory=False)
    POOL = POOL[POOL.provenance.isin(["PERMISSIVE", "STRICT"])].dropna(
        subset=["overall_pearson", "efficiency@10"]).reset_index(drop=True)
    print(f"模型池 = {len(POOL)} "
          f"(PERMISSIVE {int((POOL.provenance=='PERMISSIVE').sum())} / "
          f"STRICT {int((POOL.provenance=='STRICT').sum())})")

    # ---------- 重建预测矩阵（与 35 完全相同的制品列表与对齐逻辑）----------
    NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
               "z_raw", "y_cal", "family", "seed"}
    def is_label(k):
        kl = str(k).lower()
        return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

    ART = [
        ("experiments/v6_2/phase2_split4/split4_predictions.npz", "mas"),
        ("experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas"),
        ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
        ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
        ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
        ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
        ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
        ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
        ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
        ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    ]
    ROWS_ALL = lab.pair_row.to_numpy(int)
    pos_all = pd.Series(np.arange(len(ROWS_ALL)), index=ROWS_ALL)
    mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    mas_row = mas.row.to_numpy(int)
    pos_of_mas = pos_all.reindex(mas_row)

    want = set(POOL.model_name)
    P, NM, GR = [], [], []
    for rel, kind in ART:
        p = R / rel
        if not p.exists():
            continue
        zz = np.load(p, allow_pickle=True)
        rid = None
        for k in ("row_id", "pair_row", "rows", "row"):
            if k in zz.files:
                rid = np.asarray(zz[k]).ravel().astype(int); break
        for key in zz.files:
            if key not in want or key in NONPRED or is_label(key):
                continue
            av = np.asarray(zz[key])
            try:
                av = av.astype(float)
            except (ValueError, TypeError):
                continue
            if av.ndim != 1:
                continue
            if kind == "align":
                if len(av) != len(ROWS_ALL):
                    continue
                vec = av
            elif kind == "rowid":
                if rid is None or len(av) != len(rid):
                    continue
                m = pos_all.reindex(rid)
                if m.isna().any():
                    continue
                vec = np.full(len(ROWS_ALL), np.nan)
                vec[m.to_numpy(int)] = av
            else:
                if len(av) != len(mas_row):
                    continue
                vec = np.full(len(ROWS_ALL), np.nan)
                ok = pos_of_mas.notna().to_numpy()
                vec[pos_of_mas[ok].to_numpy(int)] = av[ok]
            P.append(vec); NM.append(key)
            GR.append(0 if POOL.loc[POOL.model_name == key,
                                    "provenance"].iloc[0] == "PERMISSIVE" else 1)
    P = np.vstack(P); NM = np.asarray(NM); GR = np.asarray(GR)
    print(f"重建预测矩阵 P = {P.shape} "
          f"(PERMISSIVE {int((GR==0).sum())} / STRICT {int((GR==1).sum())})")

    uq, inv = np.unique(TGT, return_inverse=True)
    cnt = np.bincount(inv)
    tidx = [np.where(inv == u)[0] for u in range(len(uq))]
    elig = [j for j in range(len(uq)) if cnt[j] >= MIN_N]
    print(f"n>={MIN_N} 的 target = {len(elig)}")

    def two_metrics(rows_sel):
        yv = Y[rows_sel]; tv = TGT[rows_sel]; fv = FAM[rows_sel]
        ov = np.full(P.shape[0], np.nan); ef = np.full(P.shape[0], np.nan)
        U = M.prep_universe(yv, np.zeros_like(yv), tv, fv)
        idx = U["idx"]; el = U["eligible"]
        for i in range(P.shape[0]):
            pi = P[i, rows_sel]
            m = np.isfinite(yv) & np.isfinite(pi)
            if m.sum() < 50:
                continue
            U2 = M.prep_universe(yv[m], pi[m], tv[m], fv[m])
            ov[i] = M.axis_calibration(U2)["overall_pearson"]
            ef[i] = (M.axis_selection(U2, ks=(10,), ef_pcts=()))["efficiency@10"]
        return ov, ef

    def rho_of(ov, ef, mask):
        v = np.isfinite(ov) & np.isfinite(ef) & mask
        if v.sum() < 6:
            return np.nan
        return float(spearmanr(ov[v], ef[v]).correlation)

    mperm = GR == 0; mstr = GR == 1
    t0 = time.time()
    ov0, ef0 = two_metrics(np.arange(len(Y)))
    p0 = rho_of(ov0, ef0, mperm); s0 = rho_of(ov0, ef0, mstr)
    print(f"\n点估计（规范实现）: PERMISSIVE rho={p0:+.4f} | "
          f"STRICT rho={s0:+.4f} | Δ={p0-s0:+.4f}   "
          f"[{time.time()-t0:.1f}s/次迭代]")
    print(f"  对照 pool 文件: PERMISSIVE +0.8088 / STRICT -0.0289")

    rs = np.random.RandomState(SEED)
    ps, ss, ds = [], [], []
    t0 = time.time(); valid = 0
    for it in range(a.n_boot):
        pick = rs.randint(0, len(elig), len(elig))
        rows_sel = np.concatenate([tidx[elig[j]] for j in pick])
        ov, ef = two_metrics(rows_sel)
        rp = rho_of(ov, ef, mperm); rq = rho_of(ov, ef, mstr)
        if not (np.isfinite(rp) and np.isfinite(rq)):
            continue
        valid += 1
        ps.append(rp); ss.append(rq); ds.append(rp - rq)
        if (it + 1) % 100 == 0:
            el = (time.time() - t0) / 60
            print(f"  {it+1}/{a.n_boot} valid={valid} {el:.1f}min "
                  f"ETA {el/(it+1)*(a.n_boot-it-1):.1f}min", flush=True)
            if el > a.max_minutes:
                print("  软时限，停止"); break

    ps = np.asarray(ps); ss = np.asarray(ss); ds = np.asarray(ds)
    print(f"\n有效迭代 = {valid}/{a.n_boot} ({(time.time()-t0)/60:.1f} min)")
    out = []
    for nm, v, pt in (("PERMISSIVE", ps, p0), ("STRICT", ss, s0),
                      ("DELTA", ds, p0 - s0)):
        if len(v) < 20:
            print(f"  {nm}: 样本过少"); continue
        lo, hi = np.percentile(v, [2.5, 97.5])
        ex = bool(lo > 0 or hi < 0)
        print(f"  {nm:<11} point={pt:+.4f}  95%CI=[{lo:+.4f},{hi:+.4f}]  "
              f"{'排除0' if ex else '含0'}")
        out.append({"quantity": nm, "point": float(pt), "ci_low": float(lo),
                    "ci_high": float(hi), "excludes_zero": ex})
    pd.DataFrame(out).to_csv(EXT / "boundary_bootstrap_clean.csv", index=False,
                             encoding="utf-8")
    pl, ph = np.percentile(ps, [2.5, 97.5]); sl, sh = np.percentile(ss, [2.5, 97.5])
    overlap = not (ph < sl or sh < pl)
    meta = {"implementation": "v8metrics.axis_calibration / axis_selection (canonical)",
            "n_boot_requested": a.n_boot, "n_valid": valid, "seed": SEED,
            "cluster_unit": "target", "K": K, "min_n": MIN_N,
            "point_PERMISSIVE": float(p0), "point_STRICT": float(s0),
            "ci_overlap": bool(overlap),
            "verdict": ("BOUNDARY_ROBUST" if not overlap else "BOUNDARY_NOT_ROBUST"),
            "note": ("取代 36_boundary_bootstrap.py（其内联 overall_pearson 有 bug）")}
    (EXT / "boundary_bootstrap_clean_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n两组 CI {'重叠' if overlap else '**不重叠**'}")
    print("写出 boundary_bootstrap_clean.csv / _meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
