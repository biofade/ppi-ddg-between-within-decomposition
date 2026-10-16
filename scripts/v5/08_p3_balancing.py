# -*- coding: utf-8 -*-
"""scripts/v5/08_p3_balancing.py -- V5.0 P3: Complex / Family Balanced Boosting.

动机（V5 §57）
-------------
不同 complex 的 pair 数差异极大（最多 ~573）。row-wise LightGBM 隐含认为
573 rows ≈ 573 个独立观测，这**不成立**（真实独立样本量 ≈ complex 数）。

实验（独立于 P1/P2，按 V5 §0.3）
-------------------------------
* **W0** equal row weight（基线）
* **W1** `w_i = 1/√n_complex(i)`
* **W2** `w_i = 1/n_complex(i)`
* **B-K** complex-balanced bagging：每个 member **先抽 complex**，再在 complex
  内最多抽 K mutations；K ∈ {5, 10, 20}；20 members；报 rank-mean 与 pred-mean。

评价重点（V5 §65/§66）
----------------------
不只 mean，还要看 **cross-split SD**：在 N 个随机 complex-level 划分上分别评估，
比较 W0 与 W2 的 SD。若 SD 下降 ≥20% 即值得保留（即使 mean 不变）。

禁止（V5 §60）
--------------
weight 的任何参数都不得基于 test label distribution。

输出
----
* `experiments/v5_0/p3_balancing/p3_results.csv`
* `experiments/v5_0/p3_balancing/p3_stability.csv`
* `experiments/v5_0/p3_balancing/p3_report.md`
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

OUTDIR = C.V5 / "p3_balancing"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "p3_results.csv"
OUT_STAB = OUTDIR / "p3_stability.csv"
OUT_PRED = OUTDIR / "p3_predictions.npz"
OUT_MD = OUTDIR / "p3_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 P3 balancing")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    ap.add_argument("--trees", type=int, default=800)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--n-members", type=int, default=20)
    ap.add_argument("--n-splits", type=int, default=12,
                    help="cross-split SD 用的随机 complex 级划分数")
    ap.add_argument("--stab-trees", type=int, default=128)
    ap.add_argument("--skip-bagging", action="store_true",
                    help="跳过 complex bagging（K=5 试点已观测为负结果且极耗 CPU）")
    ap.add_argument("--bag-ks", type=int, nargs="*", default=[5, 10, 20],
                    help="bagging 的 K 取值")
    ap.add_argument("--bag-seeds", type=int, nargs="*", default=[42],
                    help="bagging 只在这些 seed 上跑（默认仅 42，控制成本）")
    args = ap.parse_args()

    lgb = C.load_lgb()
    from sklearn.model_selection import GroupKFold
    t0 = time.time()

    emit("# V5.0 P3 — Complex / Family Balanced Boosting")
    emit("")
    emit("> 动机（V5 §57）: 单 complex 最多 ~573 pairs，row-wise 训练隐含"
         "「573 rows ≈ 573 独立观测」，这与 cluster bootstrap 的结论矛盾。")
    emit("")

    d = C.build_pair_features(verbose=True)
    blk, y, groups = d["X_blocks"], d["y"], d["groups"]
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    tr_i, va_i, te_i = np.where(mtr)[0], np.where(mva)[0], np.where(mte)[0]
    df = d["df"]

    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    Sx = np.hstack([blk["phys_mut"] + blk["phys_wt"],
                    blk["Z_mut"] + blk["Z_wt"]]).astype(np.float32)
    Dn, _ = C.standardize_train(Dx, mtr)
    Sn, _ = C.standardize_train(Sx, mtr)
    dg_mut = df["dg_mut"].to_numpy(float)

    # ---- 复现 P2-1（OOF ΔG_mut 辅助列）作为统一底座 ----
    X_aux = np.hstack([Sn, Dn])
    oof = np.full(len(y), np.nan)
    for ia, ib in GroupKFold(n_splits=args.folds).split(tr_i, groups=groups[tr_i]):
        a, b = tr_i[ia], tr_i[ib]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
        m.fit(X_aux[a], dg_mut[a]); oof[b] = m.predict(X_aux[b])
    m_aux = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
    m_aux.fit(X_aux[mtr], dg_mut[mtr])
    aux_te = m_aux.predict(X_aux[te_i])

    Daux = np.hstack([Dn, oof[:, None]]).astype(np.float32)
    Daux_te = np.hstack([Dn[te_i], aux_te[:, None]]).astype(np.float32)
    Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)

    # ---- complex 规模（只用 train 统计，绝不用 test）----
    tr_complex = groups[mtr]
    ncomp_map = pd.Series(tr_complex).value_counts().to_dict()
    n_comp_train = len(ncomp_map)
    emit("")
    emit("## 0. Train 内 complex 规模分布（仅 train，未使用 test 信息）")
    emit("")
    ncs = np.array([ncomp_map[c] for c in tr_complex])
    emit(f"- train complex 数: **{n_comp_train}**；pairs: {int(mtr.sum())}")
    emit(f"- 每 complex pairs: 中位 {int(np.median(list(ncomp_map.values())))}，"
         f"均值 {np.mean(list(ncomp_map.values())):.1f}，最大 {max(ncomp_map.values())}")
    emit(f"- 最大 complex 占比: **{max(ncomp_map.values()) / mtr.sum() * 100:.2f}%** "
         f"of train rows")
    top5 = sorted(ncomp_map.values(), reverse=True)[:5]
    emit(f"- top-5 complex 的 pairs: {top5} → 合计占 "
         f"**{sum(top5) / mtr.sum() * 100:.1f}%**")
    emit("")

    Xc = Daux[mtr]
    Xr = Daux_rev[mtr]
    Xall = np.vstack([Xc, Xr])
    yall = np.r_[y[mtr], -y[mtr]]
    gall = np.r_[tr_complex, tr_complex]

    def w_of(scheme: str) -> np.ndarray:
        base = np.array([1.0 / ncomp_map[c] for c in tr_complex])
        if scheme == "W0":
            b = np.ones(len(tr_complex))
        elif scheme == "W1":
            b = 1.0 / np.sqrt(np.array([ncomp_map[c] for c in tr_complex], float))
        else:
            b = base
        return np.r_[b, b]

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}
    variants = ["W0_equal", "W1_inv_sqrt", "W2_inv_n"]
    emit("## 1. 加权方案（test n=%d, n_complex=%d）" % (int(mte.sum()),
                                                       len(set(groups[mte]))))
    emit("")
    emit("| seed | variant | Pearson | Spearman | macro-PDB ρ | macro-PDB r | RMSE | AUROC |")
    emit("|---:|---|---:|---:|---:|---:|---:|---:|")
    for sd in args.seeds:
        for vname, scheme in zip(variants, ("W0", "W1", "W2")):
            w = w_of(scheme)
            m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                     "random_state": sd})
            m.fit(Xall, yall, sample_weight=w)
            pt = m.predict(Daux_te)
            mm = C.all_metrics(y[te_i], pt, groups[te_i])
            preds[f"{vname}_s{sd}"] = pt.astype(np.float32)
            rows.append({"experiment_id": f"V5_P3_{vname}_s{sd}", "family": "P3_balancing",
                         "variant": vname, "seed": sd, "train_n": int(mtr.sum()),
                         "test_n": int(mte.sum()),
                         "train_complex_n": n_comp_train,
                         "test_complex_n": len(set(groups[mte])),
                         "status": "OK", "notes": f"sample_weight={scheme}", **mm})
            emit(f"| {sd} | {vname} | {mm['overall_pearson']:+.4f} | "
                 f"{mm['overall_spearman']:+.4f} | {mm['macro_pdb_spearman']:+.4f} | "
                 f"{mm['macro_pdb_pearson']:+.4f} | {mm['rmse']:.4f} | "
                 f"{mm['auroc']:.4f} |")

    # ---- P3-B: complex-balanced bagging ----
    emit("")
    emit("## 2. Complex-balanced bagging（先抽 complex，再抽 ≤K mutations）")
    emit("")
    emit("| seed | K | Pearson | Spearman | macro-PDB ρ | 聚合 |")
    emit("|---:|---:|---:|---:|---:|---|")
    tr_pid = np.array([str(x) for x in groups[mtr]])
    uniq_tr = np.unique(tr_pid)
    if args.skip_bagging:
        emit("> **complex-balanced bagging 已跳过**（`--skip-bagging`）："
             "K=5 试点已观测为显著负结果，且每 member 需完整训练，"
             "按 V5 §125「不允许 endless experimentation」不再全量运行。")
    for sd in ([] if args.skip_bagging else args.bag_seeds):
        for K in args.bag_ks:
            rs = np.random.RandomState(sd * 100 + K)
            mems = []
            for mem in range(args.n_members):
                pick = rs.choice(uniq_tr, len(uniq_tr), replace=True)
                keep = []
                for c in pick:
                    idx = np.where(tr_pid == c)[0]
                    if len(idx) > K:
                        idx = rs.choice(idx, K, replace=False)
                    keep.append(idx)
                sel = np.concatenate(keep)
                Xb = Daux[mtr][sel]
                yb = y[mtr][sel]
                # 每个 member 内部再做 fwd/rev 增广
                Xb2 = np.vstack([Xb, np.hstack([-Xb[:, :Dn.shape[1]], Xb[:, Dn.shape[1]:]])])
                yb2 = np.r_[yb, -yb]
                mm_ = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                           "random_state": sd + mem})
                mm_.fit(Xb2, yb2)
                mems.append(mm_.predict(Daux_te))
            Pm = np.vstack(mems)
            from scipy.stats import rankdata
            for agg, pv in (("pred_mean", Pm.mean(axis=0)),
                            ("rank_mean", np.mean([rankdata(r) for r in Pm], axis=0))):
                mmx = C.all_metrics(y[te_i], pv, groups[te_i])
                key = f"B{K}_{agg}"
                preds[f"{key}_s{sd}"] = pv.astype(np.float32)
                rows.append({"experiment_id": f"V5_P3_{key}_s{sd}",
                             "family": "P3_balancing", "variant": key, "seed": sd,
                             "train_n": int(mtr.sum()), "test_n": int(mte.sum()),
                             "train_complex_n": n_comp_train,
                             "test_complex_n": len(set(groups[mte])),
                             "status": "OK",
                             "notes": f"complex bagging K={K} members={args.n_members} agg={agg}",
                             **mmx})
                if agg == "pred_mean":
                    emit(f"| {sd} | {K} | {mmx['overall_pearson']:+.4f} | "
                         f"{mmx['overall_spearman']:+.4f} | "
                         f"{mmx['macro_pdb_spearman']:+.4f} | pred_mean |")
                else:
                    emit(f"| {sd} | {K} | {mmx['overall_pearson']:+.4f} | "
                         f"{mmx['overall_spearman']:+.4f} | "
                         f"{mmx['macro_pdb_spearman']:+.4f} | rank_mean |")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te_i].astype(np.float32),
             test_group=groups[te_i], **preds)

    # ---- P3-D: cross-split SD（W0 vs W2）----
    emit("")
    emit("## 3. Cross-split 稳定性（V5 §65/§66）")
    emit("")
    emit(f"> 在 **{args.n_splits}** 个随机 complex-level 划分上分别训练与评估"
         f"（{args.stab_trees} 树，控制成本），比较 W0 与 W2 的 Pearson **SD**。")
    emit(f"> 注意：这是**稳定性诊断**，不是主评估口径；test 不参与任何参数选择。")
    emit("")
    all_cx = np.unique(groups)
    stab_rows = []
    for scheme, vname in (("W0", "W0_equal"), ("W2", "W2_inv_n")):
        vals = []
        for si in range(args.n_splits):
            rs = np.random.RandomState(1000 + si)
            perm = rs.permutation(all_cx)
            n_te = max(1, int(round(len(all_cx) * 0.15)))
            te_cx = set(perm[:n_te])
            tmask = np.array([g in te_cx for g in groups])
            trmask = ~tmask
            if tmask.sum() < 50 or trmask.sum() < 100:
                continue
            Dn2, _ = C.standardize_train(Dx, trmask)
            Xt = np.hstack([Dn2, oof[:, None]])[tmask]
            Xa = np.hstack([Dn2, oof[:, None]])[trmask]
            Xr2 = np.hstack([-Dn2, oof[:, None]])[trmask]
            # complex 权重（仅 train）
            cmp_ = pd.Series(groups[trmask]).value_counts().to_dict()
            if scheme == "W0":
                b = np.ones(int(trmask.sum()))
            else:
                b = 1.0 / np.array([cmp_[g] for g in groups[trmask]], float)
            m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.stab_trees,
                                     "random_state": 7})
            m.fit(np.vstack([Xa, Xr2]), np.r_[y[trmask], -y[trmask]],
                  sample_weight=np.r_[b, b])
            pv = m.predict(Xt)
            vals.append(C.pearson(y[tmask], pv))
        arr = np.array([v for v in vals if v == v])
        stab_rows.append({"scheme": vname, "n_splits": len(arr),
                          "mean": float(arr.mean()), "sd": float(arr.std(ddof=1)),
                          "min": float(arr.min()), "max": float(arr.max()),
                          "p25": float(np.percentile(arr, 25)),
                          "p75": float(np.percentile(arr, 75))})
        emit(f"- **{vname}**: n_splits={len(arr)}, mean={arr.mean():+.4f}, "
             f"**SD={arr.std(ddof=1):.4f}**, range [{arr.min():+.4f}, {arr.max():+.4f}]")
    st = pd.DataFrame(stab_rows)
    st.to_csv(OUT_STAB, index=False, encoding="utf-8")
    if len(st) == 2:
        sd0, sd2 = float(st.iloc[0]["sd"]), float(st.iloc[1]["sd"])
        drop = (sd0 - sd2) / sd0 * 100 if sd0 > 0 else float("nan")
        emit("")
        emit(f"- SD 变化: {sd0:.4f} → {sd2:.4f}（**{drop:+.1f}%**）")
        emit(f"- GATE（V5 §66 要求 SD 下降 ≥20%）: "
             f"{'**满足**' if drop >= 20 else '**未满足**'}")

    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s；峰值 RAM {C.peak_ram_mb():.0f} MB")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}\n写出 {OUT_STAB}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
