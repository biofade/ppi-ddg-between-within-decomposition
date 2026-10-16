# -*- coding: utf-8 -*-
"""scripts/v6_3/04_set_calibration.py -- V6.3 PHASE 6: Candidate-Set Calibration Test.

科学问题
--------
target mean ΔΔG（`mean_i ddG_ci`）能否仅由

* **Model T**：target-only descriptors（结构/界面/链长）

预测？还是必须加入

* **Model TS**：target + **mutation-set composition**（部署可用）

才能明显改善？

**严格性**
* 单位 = target（N ≤ 203）
* **family-level CV**（相同 folds 比较 T 与 TS）—— 绝不用 random target CV 作为 primary
* **禁止**任何 label-derived 变量进入模型（ddg mean/std/beneficial fraction/rank）
* 单独报告 `>=30` mutation targets
* paired target bootstrap 给出 ΔRMSE / ΔMAE / ΔPearson 的 95% CI
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH1 = C.ROOT / "experiments" / "v6_3" / "phase1_target_table"
PH6 = C.ROOT / "experiments" / "v6_3" / "phase6_set_calibration"
PH6.mkdir(parents=True, exist_ok=True)
O_CMP = PH6 / "target_only_vs_set.csv"
O_PRED = PH6 / "predictions.csv"
O_BOOT = PH6 / "bootstrap.csv"
O_MD = PH6 / "set_calibration_report.md"
N_BOOT = 5000
L: list[str] = []

# ---- 特征白名单（严格可部署；无 label-derived） ----
TARGET_ONLY = ["struct_n_iface_res_lig", "struct_n_iface_res_rec",
               "struct_n_contact_pairs", "struct_n_contact_res_pairs",
               "struct_dsasa_total", "struct_dsasa_lig", "struct_dsasa_rec",
               "struct_interface_residue_ratio", "struct_contact_density",
               "struct_sasa_complex", "struct_n_hbond", "struct_n_salt_bridge",
               "struct_n_hydrophobic_contact", "struct_n_chainpairs_with_contact",
               "protein_length", "lig_seq_len", "rec_seq_len",
               "struct_frac_charged_iface", "struct_frac_polar_iface",
               "struct_frac_hydrophobic_iface", "struct_frac_helix_iface",
               "struct_frac_sheet_iface", "struct_chainpair_contact_density"]
SET_ONLY = ["log_n_mutations", "n_unique_mutation_sites", "mutations_per_site",
            "max_mutations_per_site", "site_entropy", "mutation_pair_entropy",
            "WT_residue_entropy", "Mut_residue_entropy", "alanine_fraction",
            "glycine_fraction", "charged_to_neutral_fraction",
            "hydrophobic_to_polar_fraction", "conservative_fraction",
            "single_mut_fraction", "mean_mutations_per_variant",
            "max_mutations_per_variant", "site_count_gini",
            "mutation_count_gini", "site_coverage_all",
            "site_coverage_interface"]


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def cv_predict(X, y, fam, seed=42):
    """family-level CV（GroupKFold by family）→ 折外预测。"""
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import GroupKFold
    n_g = len(set(fam))
    k = min(5, max(2, n_g))
    gkf = GroupKFold(n_splits=k)
    pred = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups=fam):
        mdl = RidgeCV(alphas=np.logspace(-2, 3, 20))
        mdl.fit(X[tr], y[tr])
        pred[te] = mdl.predict(X[te])
    return pred


def metrics(y, p):
    from sklearn.metrics import mean_absolute_error
    from scipy.stats import spearmanr
    if len(y) < 3:
        return dict(pearson=np.nan, spearman=np.nan, rmse=np.nan, mae=np.nan)
    r = C.pearson(y, p) if len(y) >= 3 else np.nan
    sp = spearmanr(y, p).statistic
    return dict(pearson=float(r), spearman=float(sp) if sp == sp else np.nan,
                rmse=float(np.sqrt(np.mean((y - p) ** 2))),
                mae=float(mean_absolute_error(y, p)))


def main() -> int:
    t0 = time.time()
    emit("# V6.3 PHASE 6 — Candidate-Set Calibration Test")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 单位 = **target**；**family-level CV**（T 与 TS 用相同 folds）")
    emit("- **禁止** label-derived 变量进入模型")
    emit("")
    T = pd.read_csv(PH1 / "target_mechanism_table.csv", low_memory=False)
    y = pd.to_numeric(T["true_mean_ddg"], errors="coerce").to_numpy(float)
    fam = T["family"].astype(str).to_numpy()
    keep_cols = [c for c in (TARGET_ONLY + SET_ONLY)
                 if c in T.columns and T[c].notna().sum() >= 150]
    Tfeat = [c for c in TARGET_ONLY if c in keep_cols]
    Sfeat = [c for c in SET_ONLY if c in keep_cols]
    emit(f"- 可用 target-only 特征: **{len(Tfeat)}**")
    emit(f"- 可用 set 特征: **{len(Sfeat)}**")
    emit("")

    # 标准化（全样本标准化，仅用于数值稳定；CV 内不做泄漏因为只用 train 拟合模型系数）
    def build(cols):
        X = T[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        mu = np.nanmean(X, axis=0); sd = np.nanstd(X, axis=0)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd
        return X

    XT = build(Tfeat)
    XS = build(Sfeat)
    okT = np.isfinite(y) & np.isfinite(XT).all(axis=1)
    okS = np.isfinite(y) & np.isfinite(XS).all(axis=1)
    # 为公平比较，用两者都可用的交集
    ok = okT & okS
    emit(f"- 两者都完整的 target: **{int(ok.sum())}** / {len(T)}")
    emit("")

    yv, famv = y[ok], fam[ok]
    XTv, XSv = XT[ok], XS[ok]
    XTSv = np.hstack([XTv, XSv])
    grp = pd.cut(T.loc[ok, "n_mutations"], bins=[0, 9, 29, 10 ** 9],
                 labels=["2-9", "10-29", ">=30"]).astype(str).to_numpy()

    predT = cv_predict(XTv, yv, famv)
    predTS = cv_predict(XTSv, yv, famv)

    emit("## 6.5 结果（family-level CV）")
    emit("")
    emit("| 模型 | 特征数 | n | Pearson | Spearman | RMSE | MAE |")
    emit("|---|---:|---:|---:|---:|---:|---:|")
    cmp_rows = []
    for name, pred, nf in (("T_target_only", predT, len(Tfeat)),
                           ("TS_target+set", predTS, len(Tfeat) + len(Sfeat))):
        mt = metrics(yv, pred)
        emit(f"| {name} | {nf} | {len(yv)} | {mt['pearson']:+.4f} | "
             f"{mt['spearman']:+.4f} | {mt['rmse']:.4f} | {mt['mae']:.4f} |")
        cmp_rows.append({"model": name, "n_features": nf, "subset": "ALL",
                         "n": len(yv), **mt})
    emit("")

    # >=30 单独报告
    emit("## 6.5b `>=30` mutation targets 单独报告")
    emit("")
    emit("| 模型 | n | Pearson | Spearman | RMSE | MAE |")
    emit("|---|---:|---:|---:|---:|---:|")
    for sub in ("2-9", "10-29", ">=30"):
        sel = grp == sub
        if sel.sum() < 4:
            continue
        for name, pred in (("T_target_only", predT), ("TS_target+set", predTS)):
            mt = metrics(yv[sel], pred[sel])
            emit(f"| {name} ({sub}) | {int(sel.sum())} | {mt['pearson']:+.4f} | "
                 f"{mt['spearman']:+.4f} | {mt['rmse']:.4f} | {mt['mae']:.4f} |")
            cmp_rows.append({"model": name, "n_features": np.nan, "subset": sub,
                             "n": int(sel.sum()), **mt})
    emit("")
    pd.DataFrame(cmp_rows).to_csv(O_CMP, index=False, encoding="utf-8")
    pd.DataFrame({"pdb_id": T.loc[ok, "pdb_id"].to_numpy(), "family": famv,
                  "n_mutations": T.loc[ok, "n_mutations"].to_numpy(),
                  "grp": grp, "true_mean_ddg": yv, "pred_T": predT,
                  "pred_TS": predTS}).to_csv(O_PRED, index=False,
                                             encoding="utf-8")

    # ---------- 6.6 paired bootstrap ----------
    emit("## 6.6 Paired target bootstrap（单位 = target，5000 次）")
    emit("")
    emit("| 统计量 | 全样本 | `>=30` 子集 |")
    emit("|---|---:|---:|")
    rs = np.random.RandomState(42)
    n = len(yv)
    dR, dM, dP = [], [], []
    for _ in range(N_BOOT):
        idx = rs.randint(0, n, n)
        a = metrics(yv[idx], predT[idx]); b = metrics(yv[idx], predTS[idx])
        if a["rmse"] == a["rmse"] and b["rmse"] == b["rmse"]:
            dR.append(b["rmse"] - a["rmse"])
            dM.append(b["mae"] - a["mae"])
            dP.append(b["pearson"] - a["pearson"])
    dR = np.asarray(dR); dM = np.asarray(dM); dP = np.asarray(dP)
    boot_rows = []
    for lab, arr in (("Delta_RMSE (TS-T)", dR), ("Delta_MAE (TS-T)", dM),
                     ("Delta_Pearson (TS-T)", dP)):
        lo, hi = np.percentile(arr, [2.5, 97.5])
        emit(f"| {lab} | {arr.mean():+.4f} [{lo:+.4f}, {hi:+.4f}] | — |")
        boot_rows.append({"statistic": lab, "subset": "ALL",
                          "point": float(arr.mean()), "ci_low": float(lo),
                          "ci_high": float(hi),
                          "ci_excludes_zero": bool((lo > 0) or (hi < 0))})
    # >=30 子集
    sel = grp == ">=30"
    if sel.sum() >= 6:
        ys, psT, psTS = yv[sel], predT[sel], predTS[sel]
        ns = len(ys)
        dRs = []
        for _ in range(N_BOOT):
            idx = rs.randint(0, ns, ns)
            a = metrics(ys[idx], psT[idx]); b = metrics(ys[idx], psTS[idx])
            if a["rmse"] == a["rmse"] and b["rmse"] == b["rmse"]:
                dRs.append(b["rmse"] - a["rmse"])
        dRs = np.asarray(dRs)
        if len(dRs) > 20:
            lo, hi = np.percentile(dRs, [2.5, 97.5])
            emit(f"| Delta_RMSE (TS-T) | — | {dRs.mean():+.4f} "
                 f"[{lo:+.4f}, {hi:+.4f}] |")
            boot_rows.append({"statistic": "Delta_RMSE (TS-T)", "subset": ">=30",
                              "point": float(dRs.mean()), "ci_low": float(lo),
                              "ci_high": float(hi),
                              "ci_excludes_zero": bool((lo > 0) or (hi < 0))})
    pd.DataFrame(boot_rows).to_csv(O_BOOT, index=False, encoding="utf-8")
    emit("")

    # ---------- 判定 ----------
    emit("## 判定")
    emit("")
    aT = metrics(yv, predT); aTS = metrics(yv, predTS)
    d_rmse = aTS["rmse"] - aT["rmse"]
    d_pear = aTS["pearson"] - aT["pearson"]
    # 用 bootstrap CI 判据
    rmse_row = [b for b in boot_rows
                if b["statistic"] == "Delta_RMSE (TS-T)" and b["subset"] == "ALL"]
    sig = bool(rmse_row and rmse_row[0]["ci_excludes_zero"]
               and rmse_row[0]["point"] < 0)
    emit(f"- T  (target-only):  Pearson {aT['pearson']:+.4f}, RMSE {aT['rmse']:.4f}")
    emit(f"- TS (target+set)  : Pearson {aTS['pearson']:+.4f}, RMSE {aTS['rmse']:.4f}")
    emit(f"- ΔRMSE(TS−T) = **{d_rmse:+.4f}**；ΔPearson = {d_pear:+.4f}")
    emit("")
    if aT["pearson"] != aT["pearson"] or aTS["pearson"] != aTS["pearson"]:
        verdict = "INCONCLUSIVE"
    elif sig and d_pear > 0.05:
        verdict = "SET_CONTEXT_SUPPORTED"
    else:
        verdict = "SET_CONTEXT_NOT_SUPPORTED"
    emit(f"> ### **`{verdict}`**")
    emit("")
    emit("> ⚠️ 注意：无论 T 还是 TS，若 family-CV Pearson 都很低（≈0），")
    emit("> 则说明 **target-level mean 在 family OOD 下不可由这些描述符预测**，")
    emit("> 这才是对 V7 `G(x_c, M_c)` 设计的真正约束。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    O_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {O_CMP}")
    print(f"写出 {O_PRED}")
    print(f"写出 {O_BOOT}")
    print(f"写出 {O_MD}")
    print(f"\n>>> {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
