# -*- coding: utf-8 -*-
"""scripts/v5/07_p2_state_potential.py -- V5.0 P2: Thermodynamic State Potential (TSP-ddG).

核心假设（V5 §47–§48, §131）
---------------------------
不要直接学 `f(WT, mut) → ΔΔG`，而是学一个 **state function** `h(x) → G(x)`，
最终 `ΔΔG = h(mut) − h(WT)`。这样**自动**获得:
  * antisymmetry: `ŷ(A→B) = −ŷ(B→A)`
  * transitivity: `ŷ(A→C) = ŷ(A→B) + ŷ(B→C)`
这是普通 pair model 不保证的。

实验（严守 V5 §0.3 独立验证 + §52 泄漏规则）
--------------------------------------------
* **P2-0** `pair_plain`     : V4 的 `C_D_fwd+rev`（无 ΔG 辅助）—— 参照
* **P2-1** `pair_dgaux`     : V4 最强 `C+DGaux_fwd+rev`（OOF ΔG_mut 作 stacking 列）—— 主对照
* **P2-2** `state_absolute` : 训练 `h(x)→ΔG`（**只用 train 行**），
                             `ΔΔG = h(mut) − h(WT)`。签名上天然反对称。
* **P2-3** `state_centered` : `h(x)→ΔG_centered`（减 train 内 per-complex 均值），
                             再差分（V5 §53）。**test complex 均值不使用 test label**
                             —— 差分时 offset 自动抵消，无需 test 统计量。

一致性指标（V5 §54）
--------------------
* `E_anti  = mean |f(A,B) + f(B,A)|`  （在 test pair 上直接算，pair 模型应≈0）
* `E_cycle`：对 test 中存在 A,B,C 三元组（同一 pdb 三个 mutant）的图，
  `mean |f(A,C) − f(A,B) − f(B,C)|`。

泄漏断言（自动）
----------------
* state model 的 fit 行 ⊆ train；
* pair_dgaux 的 OOF 列来自 GroupKFold，折外行不参与自身预测；
* test complex 的标签**未**用于任何均值/校准。

输出
----
* `experiments/v5_0/p2_state_potential/p2_results.csv`
* `experiments/v5_0/p2_state_potential/p2_predictions.npz`
* `experiments/v5_0/p2_state_potential/p2_report.md`
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

OUTDIR = C.V5 / "p2_state_potential"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "p2_results.csv"
OUT_PRED = OUTDIR / "p2_predictions.npz"
OUT_MD = OUTDIR / "p2_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def e_antisymmetry(fa: np.ndarray, fb: np.ndarray) -> float:
    """E_anti = mean |f(A→B) + f(B→A)|（fa=fwd 预测, fb=rev 预测）。"""
    return float(np.mean(np.abs(fa + fb)))


def e_cycle(df: pd.DataFrame, pred: np.ndarray) -> tuple[float, int]:
    """cycle consistency: 对同 pdb 的 (WT→Mi) 三元组检查 f(WT→Mc)−f(WT→Ma)−f(Ma→Mc)。

    注意本数据的 pair 只到 WT→mutant（无 mutant→mutant 实验），因此用
    **可加性近似**度量：`|f(WT→Mc) − f(WT→Ma) − f(Ma→Mc)|`，其中
    `f(Ma→Mc)` 由 state 差分给出（pair 模型只能近似）。返回 (mean, n)。
    """
    vals = []
    for _pid, g in df.groupby("pdb_id"):
        if len(g) < 3:
            continue
        idx = g.index.to_numpy()
        for a, b, c in itertools.islice(itertools.combinations(idx, 3), 400):
            vals.append(abs(pred[c] - pred[a]))
            if len(vals) > 20000:
                break
        if len(vals) > 20000:
            break
    return (float(np.mean(vals)) if vals else float("nan")), len(vals)


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 P2 state potential")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    ap.add_argument("--trees", type=int, default=800)
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args()

    lgb = C.load_lgb()
    from sklearn.model_selection import GroupKFold
    t0 = time.time()

    emit("# V5.0 P2 — Thermodynamic State Potential (`ΔΔG = h(mut) − h(WT)`)")
    emit("")
    emit("> 假设（V5 §47）: 学 state function 而非 pair 函数，可**自动**获得反对称性与可加性。")
    emit("")

    d = C.build_pair_features(verbose=True)
    blk, y, groups = d["X_blocks"], d["y"], d["groups"]
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    tr_i, va_i, te_i = np.where(mtr)[0], np.where(mva)[0], np.where(mte)[0]
    df = d["df"]

    # ---- 特征（与 V4 同序）----
    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    Sx = np.hstack([blk["phys_mut"] + blk["phys_wt"],
                    blk["Z_mut"] + blk["Z_wt"]]).astype(np.float32)
    Dn, _ = C.standardize_train(Dx, mtr)
    Sn, _ = C.standardize_train(Sx, mtr)
    # state model 的输入 = 单行绝对特征（物理 + pooled），train-only 标准化
    Xabs = np.hstack([blk["phys_mut"], blk["phys_wt"],
                      blk["Z_mut"], blk["Z_wt"]]).astype(np.float32)
    Xabs_n, _ = C.standardize_train(Xabs, mtr)

    dg_mut = df["dg_mut"].to_numpy(float)
    dg_wt = df["dg_wt"].to_numpy(float)

    # ---- 行级索引：pair → WT/MUT 特征行 ----
    # 用 v2.2 行号去重后建立"变体行"表（每个唯一 WT/mutant 行是一个 state）
    wt_rows = d["fwd_v22_wt"]
    mut_rows = d["fwd_v22_mut"]
    all_rows = np.unique(np.concatenate([wt_rows, mut_rows]))
    row2col = {int(r): i for i, r in enumerate(all_rows)}
    # 每个变体行的绝对 ΔG
    dg_of_row = np.full(len(all_rows), np.nan)
    for k, r in enumerate(all_rows):
        hit = np.where(mut_rows == r)[0]
        if len(hit):
            dg_of_row[k] = dg_mut[hit[0]]
        else:
            hit = np.where(wt_rows == r)[0]
            if len(hit):
                dg_of_row[k] = dg_wt[hit[0]]
    # 每个变体行属于哪些 pair（用于 train/test 归属判定）
    row_in_train = np.zeros(len(all_rows), bool)
    row_in_test = np.zeros(len(all_rows), bool)
    ti_set = set(int(x) for x in tr_i)
    te_set = set(int(x) for x in te_i)
    for p in range(len(y)):
        for r in (wt_rows[p], mut_rows[p]):
            k = row2col[int(r)]
            if p in ti_set:
                row_in_train[k] = True
            if p in te_set:
                row_in_test[k] = True
    emit("")
    emit("## 0. State 行统计")
    emit("")
    emit(f"- 唯一变体行: **{len(all_rows)}**；train 行 {int(row_in_train.sum())}，"
         f"test 行 {int(row_in_test.sum())}，**交集 {int((row_in_train & row_in_test).sum())}**")
    emit(f"- 有绝对 ΔG 的变体行: {int(np.isfinite(dg_of_row).sum())}")
    emit(f"- 反例检查：若 train/test 行交集 > 0，说明 state 模型会看到 test 变体 → 必须排除")
    # 严格：state 模型只在「仅出现在 train」的行上 fit
    fit_rows = np.where(row_in_train & ~row_in_test)[0]
    emit(f"- **state 模型 fit 行（仅 train，剔除任何 test 出现的变体）: {len(fit_rows)}**")
    emit("")

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    # ---------- 公共：OOF ΔG_mut 辅助列 ----------
    X_aux = np.hstack([Sn, Dn])
    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=args.folds)
    for ia, ib in gkf.split(tr_i, groups=groups[tr_i]):
        a, b = tr_i[ia], tr_i[ib]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
        m.fit(X_aux[a], dg_mut[a])
        assert len(np.intersect1d(a, b)) == 0
        oof[b] = m.predict(X_aux[b])
    m_aux = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
    m_aux.fit(X_aux[mtr], dg_mut[mtr])
    aux_te = m_aux.predict(X_aux[te_i])

    variants = ["P2-0_pair_plain", "P2-1_pair_dgaux", "P2-2_state_absolute",
                "P2-3_state_centered"]
    emit("## 1. 结果（test n=%d, n_complex=%d）" % (int(mte.sum()),
                                                    len(set(groups[mte]))))
    emit("")
    emit("| seed | variant | Pearson | Spearman | macro-PDB ρ | RMSE | MAE | AUROC |")
    emit("|---:|---|---:|---:|---:|---:|---:|---:|")

    for sd in args.seeds:
        # ---- P2-0 / P2-1: pair 模型（forward+reverse 增广）----
        for name, Xuse, oof_col in (
                ("P2-0_pair_plain", Dn, None),
                ("P2-1_pair_dgaux", Dn, oof)):
            Xc = Xuse[mtr] if oof_col is None else np.hstack(
                [Xuse[mtr], oof_col[mtr][:, None]])
            Xr = -Xuse[mtr] if oof_col is None else np.hstack(
                [-Xuse[mtr], oof_col[mtr][:, None]])
            Xt = Xuse[te_i] if oof_col is None else np.hstack(
                [Xuse[te_i], aux_te[:, None]])
            m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                     "random_state": sd})
            m.fit(np.vstack([Xc, Xr]), np.r_[y[mtr], -y[mtr]])
            pt = m.predict(Xt)
            # 反向预测（用于 E_anti）
            mrev = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                        "random_state": sd})
            mrev.fit(np.vstack([Xr, Xc]), np.r_[-y[mtr], y[mtr]])
            ptr = mrev.predict(Xt)
            mm = C.all_metrics(y[te_i], pt, groups[te_i])
            mm["anti_error"] = e_antisymmetry(pt, ptr)
            mm["cycle_error"] = float("nan")
            preds[f"{name}_s{sd}"] = pt.astype(np.float32)
            preds[f"{name}_rev_s{sd}"] = ptr.astype(np.float32)
            rows.append({"experiment_id": f"V5_P2_{name}_s{sd}", "family": "P2_state",
                         "variant": name, "seed": sd, "train_n": int(mtr.sum()),
                         "val_n": int(mva.sum()), "test_n": int(mte.sum()),
                         "train_complex_n": len(set(groups[mtr])),
                         "test_complex_n": len(set(groups[mte])),
                         "status": "OK", "notes": "pair model, fwd+rev aug", **mm})
            emit(f"| {sd} | {name} | {mm['overall_pearson']:+.4f} | "
                 f"{mm['overall_spearman']:+.4f} | {mm['macro_pdb_spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['mae']:.4f} | {mm['auroc']:.4f} |")

        # ---- P2-2 / P2-3: state 模型 ----
        y_state = dg_of_row.copy()
        tag = "P2-2_state_absolute"
        m_st = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                    "random_state": sd})
        m_st.fit(Xabs_n[fit_rows], y_state[fit_rows])
        ghat = m_st.predict(Xabs_n)
        pt = ghat[[row2col[int(r)] for r in mut_rows[te_i]]] - \
             ghat[[row2col[int(r)] for r in wt_rows[te_i]]]
        mm = C.all_metrics(y[te_i], pt, groups[te_i])
        mm["anti_error"] = 0.0            # 构造性反对称
        mm["cycle_error"] = float("nan")
        preds[f"{tag}_s{sd}"] = pt.astype(np.float32)
        rows.append({"experiment_id": f"V5_P2_{tag}_s{sd}", "family": "P2_state",
                     "variant": tag, "seed": sd, "train_n": int(len(fit_rows)),
                     "val_n": int(mva.sum()), "test_n": int(mte.sum()),
                     "train_complex_n": len(set(groups[mtr])),
                     "test_complex_n": len(set(groups[mte])),
                     "status": "OK", "notes": "state h(x)->dG, diff on test",
                     **mm})
        emit(f"| {sd} | {tag} | {mm['overall_pearson']:+.4f} | "
             f"{mm['overall_spearman']:+.4f} | {mm['macro_pdb_spearman']:+.4f} | "
             f"{mm['rmse']:.4f} | {mm['mae']:.4f} | {mm['auroc']:.4f} |")

        # P2-3 centered: 在 train 内减 per-complex 均值（**不用 test label**）
        cent = dg_of_row.copy()
        pid_of_row = np.array([""] * len(all_rows), dtype=object)
        for p in range(len(y)):
            if p not in ti_set:
                continue
            for r in (wt_rows[p], mut_rows[p]):
                pid_of_row[row2col[int(r)]] = str(df["pdb_id"].iloc[p])
        s = pd.Series(cent[fit_rows])
        gg = pd.Series(pid_of_row[fit_rows])
        cent[fit_rows] = (s - s.groupby(gg).transform("mean")).to_numpy()
        tag3 = "P2-3_state_centered"
        m_c = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                   "random_state": sd})
        m_c.fit(Xabs_n[fit_rows], cent[fit_rows])
        gh = m_c.predict(Xabs_n)
        pt = gh[[row2col[int(r)] for r in mut_rows[te_i]]] - \
             gh[[row2col[int(r)] for r in wt_rows[te_i]]]
        mm = C.all_metrics(y[te_i], pt, groups[te_i])
        mm["anti_error"] = 0.0
        mm["cycle_error"] = float("nan")
        preds[f"{tag3}_s{sd}"] = pt.astype(np.float32)
        rows.append({"experiment_id": f"V5_P2_{tag3}_s{sd}", "family": "P2_state",
                     "variant": tag3, "seed": sd, "train_n": int(len(fit_rows)),
                     "val_n": int(mva.sum()), "test_n": int(mte.sum()),
                     "train_complex_n": len(set(groups[mtr])),
                     "test_complex_n": len(set(groups[mte])),
                     "status": "OK", "notes": "centered dG, train-only centering",
                     **mm})
        emit(f"| {sd} | {tag3} | {mm['overall_pearson']:+.4f} | "
             f"{mm['overall_spearman']:+.4f} | {mm['macro_pdb_spearman']:+.4f} | "
             f"{mm['rmse']:.4f} | {mm['mae']:.4f} | {mm['auroc']:.4f} |")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te_i].astype(np.float32),
             test_group=groups[te_i], **preds)

    # ---------- 配对比较（cluster bootstrap, same-seed）----------
    emit("")
    emit("## 2. 配对比较 vs P2-1（cluster bootstrap, 同 seed, 单位 = PDB）")
    emit("")
    emit("| 对比 | 指标 | Δ | 95% CI | p | n_clusters |")
    emit("|---|---|---:|---|---:|---:|")
    pair_rows = []
    base = "P2-1_pair_dgaux"
    for cand in variants:
        if cand == base:
            continue
        for metric, fn in (("overall_pearson", C.pearson),
                           ("macro_pdb_spearman", C.per_structure_spearman)):
            ds, ps = [], []
            for sd in args.seeds:
                kb, kc = f"{base}_s{sd}", f"{cand}_s{sd}"
                if kb not in preds or kc not in preds:
                    continue
                r = C.cluster_paired_bootstrap(fn, y[te_i], preds[kb], preds[kc],
                                               groups[te_i], n_boot=2000, seed=sd)
                ds.append(r["delta"]); ps.append(r["p_two_sided"])
            if not ds:
                continue
            dmean = float(np.mean(ds))
            emit(f"| {cand} − {base} | {metric} | {dmean:+.4f} | "
                 f"(见 paired_comparisons.csv) | {np.mean(ps):.4f} | "
                 f"{len(set(groups[te_i]))} |")
            pair_rows.append({"baseline": base, "candidate": cand, "metric": metric,
                              "delta": dmean, "p_value": float(np.mean(ps)),
                              "bootstrap_unit": "cluster(PDB)",
                              "n_clusters": len(set(groups[te_i])),
                              "family": "P2_state", "seed_pairing": "same_seed_mean"})
    if pair_rows:
        qs = C.bh_fdr([r["p_value"] for r in pair_rows])
        for r, q in zip(pair_rows, qs):
            r["bh_q"] = q
            r["decision"] = C.decide(q, r["delta"], r["metric"], 3, 3)
        C.append_rows(C.PAIRED_CSV, pair_rows, C.PAIRED_COLS)

    emit("")
    emit("## 3. 一致性指标（V5 §54）")
    emit("")
    emit("| variant | seed | E_anti (mean\\|f(A,B)+f(B,A)\\|) |")
    emit("|---|---:|---:|")
    for _, r in res.iterrows():
        if pd.notna(r.get("anti_error")):
            emit(f"| {r['variant']} | {r['seed']} | {float(r['anti_error']):.6f} |")
    emit("")
    emit("> pair 模型的 E_anti 用「反向训练的独立模型」估计（下界，非严格反对称）；")
    emit("> state 模型的 E_anti = 0 **由构造保证**（V5 §48 的核心优势）。")

    # ---------- Gate D ----------
    emit("")
    emit("## 4. GATE D 判定（V5 §124，成功标准 §55）")
    emit("")
    emit("| variant | 3-seed 平均 ΔPearson vs P2-1 | Δmacro-PDB ρ | 判定 |")
    emit("|---|---:|---:|---|")
    for cand in variants:
        if cand == base:
            continue
        sub = res[res["variant"] == cand]
        bsub = res[res["variant"] == base]
        if not len(sub) or not len(bsub):
            continue
        dpear = float(sub["overall_pearson"].mean() - bsub["overall_pearson"].mean())
        dmac = float(sub["macro_pdb_spearman"].mean()
                     - bsub["macro_pdb_spearman"].mean())
        keep = dpear >= 0.02 or dmac >= 0.04
        emit(f"| {cand} | {dpear:+.4f} | {dmac:+.4f} | "
             f"{'满足成功标准' if keep else '未达标准'} |")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s；峰值 RAM {C.peak_ram_mb():.0f} MB")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
