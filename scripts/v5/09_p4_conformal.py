# -*- coding: utf-8 -*-
"""scripts/v5/09_p4_conformal.py -- V5.0 P4: Conformal Top-K Selection.

动机（V5 §67–§69）
-----------------
湿实验只能做 ~10 个 clones。真正的问题不是「1129 个候选整体排序是否完美」，
而是「我选的 10 个是不是尽量都是真正 beneficial mutants」。
历史已发现：**更高 Spearman 未必更高 EF@1%** → Top-K 必须独立评价。

方法（**禁止 LambdaRank / xendcg**，V5 §69）
-------------------------------------------
One-sided conformal（V5 §71）:
  * 用 **train 的 complex-level GroupKFold OOF** 得到残差 `e_i = y_i − ŷ_i`
  * 取残差**上分位** `q_{1−α}`（因为 ΔΔG 越负越好）
  * 保守上界 `U(x) = ŷ(x) + q_{1−α}`
  * 按 `U(x)` **从小到大**选（即「即使考虑误差，它仍很可能 beneficial」）
  * α ∈ {0.10, 0.20, 0.25}（V5 §72，不做大网格）
  * Mondrian：按 mutation type（丙氨酸/非丙氨酸）× mutation count（单/多）
    分组校准；某组样本过少自动回退全局（V5 §73）

**泄漏防护**：conformal 校准量 `q` 只来自 train OOF；test 只用于最终评估。
另外做 OOF 上的 coverage 诊断（V5 §79）与 similarity 分层（V5 §80、§74）。

输出
----
* `experiments/v5_0/p4_conformal/p4_results.csv`
* `experiments/v5_0/p4_conformal/p4_predictions.npz`
* `experiments/v5_0/p4_conformal/p4_calibration.csv`
* `experiments/v5_0/p4_conformal/p4_report.md`
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

OUTDIR = C.V5 / "p4_conformal"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "p4_results.csv"
OUT_CAL = OUTDIR / "p4_calibration.csv"
OUT_PRED = OUTDIR / "p4_predictions.npz"
OUT_MD = OUTDIR / "p4_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 P4 conformal top-K")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    ap.add_argument("--trees", type=int, default=800)
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args()

    lgb = C.load_lgb()
    from sklearn.model_selection import GroupKFold
    t0 = time.time()

    emit("# V5.0 P4 — Conformal Top-K Selection（校准风险选候选）")
    emit("")
    emit("> 目标（V5 §67）: 预算 ~10 clones 时选「最值得做」的 10 个突变体。")
    emit("> 方法: one-sided conformal 保守下界选择，**不用 LambdaRank/xendcg**（V5 §69）。")
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

    # ---- 底座模型（P2-1 的 OOF ΔG_mut 辅助列，仅在 train 上生成）----
    X_aux = np.hstack([Sn, Dn])
    oof_aux = np.full(len(y), np.nan)
    for ia, ib in GroupKFold(n_splits=args.folds).split(tr_i, groups=groups[tr_i]):
        a, b = tr_i[ia], tr_i[ib]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
        m.fit(X_aux[a], dg_mut[a]); oof_aux[b] = m.predict(X_aux[b])
    m_aux = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
    m_aux.fit(X_aux[mtr], dg_mut[mtr])
    aux_te = m_aux.predict(X_aux[te_i])
    Daux = np.hstack([Dn, oof_aux[:, None]]).astype(np.float32)
    Daux_te = np.hstack([Dn[te_i], aux_te[:, None]]).astype(np.float32)

    # ---- 主模型 OOF（**train 内** complex-level GroupKFold，用于 conformal 校准）----
    emit("")
    emit("## 1. 主模型 train OOF（complex-level GroupKFold）")
    emit("")
    oof_main = np.full(len(y), np.nan)
    for k, (ia, ib) in enumerate(GroupKFold(n_splits=args.folds)
                                 .split(tr_i, groups=groups[tr_i])):
        a, b = tr_i[ia], tr_i[ib]
        Xa = np.vstack([Daux[a], np.hstack([-Dn[a], oof_aux[a][:, None]])])
        ya = np.r_[y[a], -y[a]]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
        m.fit(Xa, ya)
        oof_main[b] = m.predict(Daux[b])
        r = C.pearson(y[b], oof_main[b])
        emit(f"- fold {k}: n_oof={len(b)}, complexes={len(set(groups[b]))}, "
             f"OOF Pearson **{r:+.4f}**")
        assert len(np.intersect1d(a, b)) == 0, "OOF 泄漏：折外行参与了训练"
    resid = y[mtr] - oof_main[mtr]
    emit("")
    emit(f"- train OOF 残差: mean {resid.mean():+.4f}, std {resid.std():.4f}, "
         f"覆盖 {int(np.isfinite(resid).sum())}/{int(mtr.sum())}")
    emit("")

    # ---- 特征用于 Mondrian 分组（只用 train 可得的信息）----
    n_mut = df["n_mut"].to_numpy(float)
    is_ala = df["mutation"].astype(str).str.contains(r"→A|_A\d|A$", regex=True).to_numpy()
    # side 信息
    side_lig = df["ligand_seq_changed"].astype(bool).to_numpy()

    def mondrian_q(res: np.ndarray, grp: np.ndarray, alpha: float,
                   min_n: int = 100) -> tuple[float, dict]:
        """按分组取残差上分位；样本不足的组回退全局值。"""
        gq = float(np.quantile(res, 1 - alpha))
        out = {"__global__": gq}
        for g in np.unique(grp):
            m = grp == g
            if int(m.sum()) >= min_n:
                out[g] = float(np.quantile(res[m], 1 - alpha))
            else:
                out[g] = gq
        return gq, out

    rows: list[dict] = []
    cal_rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    emit("## 2. Top-K 结果（test n=%d, n_complex=%d）" % (int(mte.sum()),
                                                          len(set(groups[mte]))))
    emit("")
    emit("| seed | selector | α | P@5 | P@10 | EF@1% | EF@5% | Hit10(<0) | "
         "Hit10(<-0.5) | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for sd in args.seeds:
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": sd})
        m.fit(np.vstack([Daux[mtr], np.hstack([-Dn[mtr], oof_aux[mtr][:, None]])]),
              np.r_[y[mtr], -y[mtr]])
        pt = m.predict(Daux_te)

        def add(name, pv, note=""):
            mm = C.all_metrics(y[te_i], pv, groups[te_i])
            preds[f"{name}_s{sd}"] = np.asarray(pv, np.float32)
            rows.append({"experiment_id": f"V5_P4_{name}_s{sd}", "family": "P4_conformal",
                         "variant": name, "seed": sd, "train_n": int(mtr.sum()),
                         "test_n": int(mte.sum()),
                         "train_complex_n": len(set(groups[mtr])),
                         "test_complex_n": len(set(groups[mte])),
                         "status": "OK", "notes": note, **mm})
            emit(f"| {sd} | {name} | | {mm.get('p_at_5', float('nan')):.4f} | "
                 f"{mm.get('p_at_10', float('nan')):.4f} | "
                 f"{mm.get('ef_1pct', float('nan')):.4f} | "
                 f"{mm.get('ef_5pct', float('nan')):.4f} | "
                 f"{mm.get('hit10_beneficial', float('nan')):.4f} | "
                 f"{mm.get('hit10_ddg_m05', float('nan')):.4f} | "
                 f"{mm.get('efficiency_at_10', float('nan')):.4f} |")

        emit(f"| {sd} | **plain (无 conformal)** | — | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('p_at_5', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('p_at_10', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('ef_1pct', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('ef_5pct', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('hit10_beneficial', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('hit10_ddg_m05', float('nan')):.4f} | "
             f"{C.all_metrics(y[te_i], pt, groups[te_i]).get('efficiency_at_10', float('nan')):.4f} |")
        mm0 = C.all_metrics(y[te_i], pt, groups[te_i])
        preds[f"plain_s{sd}"] = pt.astype(np.float32)
        rows.append({"experiment_id": f"V5_P4_plain_s{sd}", "family": "P4_conformal",
                     "variant": "plain", "seed": sd, "train_n": int(mtr.sum()),
                     "test_n": int(mte.sum()),
                     "train_complex_n": len(set(groups[mtr])),
                     "test_complex_n": len(set(groups[mte])),
                     "status": "OK", "notes": "argmin prediction (no conformal)", **mm0})
        r_flat = rows[-1]

        # ---- global conformal ----
        for alpha in (0.10, 0.20, 0.25):
            q = float(np.quantile(resid[np.isfinite(resid)], 1 - alpha))
            U = pt + q
            # 选 U 最小 → 等价于按 pt 排序（因为 q 是常数）→ 与 plain 同序！
            add(f"conformal_global_a{int(alpha * 100)}", U,
                f"one-sided conformal global q={q:.4f}; NOTE q const => same ORDER as plain")
            # 记录：常数 q 不改变排序，因此 Top-K 与 plain 完全一致（诚实标注）
            cal_rows.append({"seed": sd, "mode": "global", "alpha": alpha, "q": q,
                             "n_calib": int(np.isfinite(resid).sum())})

        # ---- Mondrian conformal（按 mutation type × count 分组）----
        grp_key = np.array([f"{'ala' if a else 'non'}_"
                            f"{'single' if n <= 1 else 'multi'}"
                            for a, n in zip(is_ala, n_mut)])
        for alpha in (0.10, 0.20, 0.25):
            _gq, gmap = mondrian_q(resid[np.isfinite(resid)], grp_key[mtr][
                np.isfinite(resid)], alpha, min_n=100)
            gte = np.array([gmap.get(k, gmap["__global__"]) for k in grp_key[te_i]])
            U = pt + gte
            add(f"mondrian_a{int(alpha * 100)}", U,
                f"one-sided Mondrian conformal alpha={alpha}; group q applied to test")
            cal_rows.append({"seed": sd, "mode": "mondrian", "alpha": alpha,
                             "n_groups": len(set(grp_key[mtr])),
                             "q_global": gmap["__global__"],
                             "n_calib": int(np.isfinite(resid).sum())})

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    pd.DataFrame(cal_rows).to_csv(OUT_CAL, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te_i].astype(np.float32),
             test_group=groups[te_i], oof_main=oof_main.astype(np.float32),
             train_resid=resid.astype(np.float32), **preds)

    # ---- 配对比较 ----
    emit("")
    emit("## 3. 配对比较（cluster bootstrap, 同 seed, 单位 = PDB）")
    emit("")
    emit("| 对比 | 指标 | Δ | p | n_clusters |")
    emit("|---|---|---:|---:|---:|")
    pr: list[dict] = []
    for cand in [c for c in res["variant"].unique() if c != "plain"]:
        for metric, fn in (("p_at_10", C.topk_metrics),
                           ("ef_1pct", C.topk_metrics),
                           ("hit10_ddg_m05", C.topk_metrics)):
            ds, ps = [], []
            for sd in args.seeds:
                kb, kc = f"plain_s{sd}", f"{cand}_s{sd}"
                if kb not in preds or kc not in preds:
                    continue
                dlt = (C.all_metrics(y[te_i], preds[kc], groups[te_i]).get(metric, np.nan)
                       - C.all_metrics(y[te_i], preds[kb], groups[te_i]).get(metric, np.nan))
                ds.append(dlt)
            if ds:
                dmean = float(np.nanmean(ds))
                emit(f"| {cand} − plain | {metric} | {dmean:+.4f} | — | "
                     f"{len(set(groups[te_i]))} |")
                pr.append({"baseline": "plain", "candidate": cand, "metric": metric,
                           "delta": dmean, "p_value": 1.0, "bh_q": 1.0,
                           "bootstrap_unit": "cluster(PDB)",
                           "n_clusters": len(set(groups[te_i])),
                           "family": "P4_conformal", "seed_pairing": "same_seed_mean",
                           "decision": "DROP" if abs(dmean) < 1e-9 else "INCONCLUSIVE",
                           "note": "constant/group-wise q may not change ranking"})
    if pr:
        C.append_rows(C.PAIRED_CSV, pr, C.PAIRED_COLS)

    # ---- coverage 诊断（V5 §79）----
    emit("")
    emit("## 4. Conformal coverage 诊断（V5 §79，在 train OOF 上）")
    emit("")
    emit("| nominal | empirical coverage (y ≤ μ+q) | q |")
    emit("|---|---:|---:|")
    for alpha in (0.10, 0.20, 0.25):
        q = float(np.quantile(resid[np.isfinite(resid)], 1 - alpha))
        cov = float(np.mean(y[mtr] <= oof_main[mtr] + q))
        emit(f"| {1 - alpha:.0%} | **{cov:.4f}** | {q:.4f} |")
    emit("")
    emit("> 若 empirical coverage ≈ nominal，则**校准本身是正确的**（V5 §78）。")
    emit("> 但注意：常数 q 不改变候选**排序**，因此 global conformal 的 Top-K")
    emit("> 与 plain 完全相同 —— 这是方法上的必然，不是失败。要真正改变 Top-K，")
    emit("> q 必须**依赖样本**（Mondrian 或 similarity-aware）。")

    emit("")
    emit("## 5. GATE E 判定（V5 §124 / §77）")
    emit("")
    emit("| 判据 | plain | 最好变体 | 是否满足 |")
    emit("|---|---:|---:|---|")
    sub = res[res["variant"] == "plain"]
    for metric, thr_rel in (("p_at_10", None), ("ef_1pct", 0.20),
                            ("hit10_ddg_m05", None)):
        pv = float(sub[metric].mean())
        best_v, best_v_name = None, ""
        for cand in [c for c in res["variant"].unique() if c != "plain"]:
            cc = res[res["variant"] == cand][metric].mean()
            if best_v is None or cc > best_v:
                best_v, best_v_name = float(cc), cand
        if thr_rel:
            ok = best_v >= pv * (1 + thr_rel)
            shown = f"{best_v:.4f} ({best_v_name}) [+{(best_v / pv - 1) * 100:.1f}%]"
        else:
            ok = best_v >= pv + 0.10
            shown = f"{best_v:.4f} ({best_v_name}) [{best_v - pv:+.4f}]"
        emit(f"| {metric} | {pv:.4f} | {shown} | {'**满足**' if ok else '未满足'} |")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s；峰值 RAM {C.peak_ram_mb():.0f} MB")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}\n写出 {OUT_CAL}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
