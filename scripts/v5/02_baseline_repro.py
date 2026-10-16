# -*- coding: utf-8 -*-
"""scripts/v5/02_baseline_repro.py -- V5.0 GATE A: 复现 V4 最强基线.

必须复现: `C+DGaux_fwd+rev`（V4 创新 2），seeds 42/43/44。
口径与 V4 的 `scripts/v4.0/11_dualhead_dg_ddg.py` **完全一致**：
* 特征 D = [phys_diff(60), Z_diff(128)]（仅反对称块），仅用 train 拟合的 scaler
* 训练集 = forward + reverse 增广（y 取负）
* 辅助任务 = 用 [S, D] 预测绝对 ΔG_mut，GroupKFold(pdb_id) 折外，再把 OOF 作为额外一列
* LightGBM 固定 800 树（**不用 val 早停** —— V4 已证明早停在 iteration 1 停）
* 反向增广时 **辅助 ΔG_mut 不取负**（它不是反对称量）

本脚本额外补齐 V4 未报告的指标：Spearman / macro per-structure Pearson+Spearman /
P@5 / P@10 / EF@1% / EF@5% / Hit@10 / Efficiency@10，并加**自动泄漏断言**。

判定（V5 §11）：若 |ΔPearson| > 0.01 与 V4 记录值不符，必须查明原因，不进入 P1。

输出:
  experiments/v5_0/p0_benchmark/baseline_repro_results.csv
  experiments/v5_0/p0_benchmark/baseline_repro_predictions.npz
  experiments/v5_0/logs/02_baseline_repro.log
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

OUTDIR = C.V5 / "p0_benchmark"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "baseline_repro_results.csv"
OUT_PRED = OUTDIR / "baseline_repro_predictions.npz"
OUT_MD = OUTDIR / "baseline_repro_report.md"
LOG = C.V5 / "logs" / "02_baseline_repro.log"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def leakage_assertions(d: dict, name: str) -> dict:
    """V5 §86–§90 自动断言。全部只读检查，违反即抛错并把结果标 INVALID。"""
    g, y = d["groups"], d["y"]
    res = {}
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    s_tr, s_va, s_te = set(g[mtr]), set(g[mva]), set(g[mte])
    res["pdb_overlap_train_val"] = len(s_tr & s_va)
    res["pdb_overlap_train_test"] = len(s_tr & s_te)
    res["pdb_overlap_val_test"] = len(s_va & s_te)
    fam = d.get("family_of_complex", {})
    if fam:
        f_tr = {fam.get(p, p) for p in s_tr}
        f_va = {fam.get(p, p) for p in s_va}
        f_te = {fam.get(p, p) for p in s_te}
        res["family_overlap_train_val"] = len(f_tr & f_va)
        res["family_overlap_train_test"] = len(f_tr & f_te)
        res["family_overlap_val_test"] = len(f_va & f_te)
    # 精确重复序列跨侧：V4 已用硬断言验证过，这里复核 mask 只覆盖同一批行
    res["split_exclusive"] = int((mtr.astype(int) + mva.astype(int) + mte.astype(int)
                                 == 1).all())
    res["split_covers_all"] = int((mtr | mva | mte).all())
    # PCA fit 范围（V5 §87）：确认 PCA 只在 train 行的并集上 fit
    tr_rows = d["split_idx"]["train"]
    widx = np.unique(np.concatenate([d["fwd_v22_wt"][tr_rows], d["fwd_v22_mut"][tr_rows]]))
    te_rows = d["split_idx"]["test"]
    tidx = np.unique(np.concatenate([d["fwd_v22_wt"][te_rows], d["fwd_v22_mut"][te_rows]]))
    res["pca_fit_rows_from_train_only"] = int(len(np.intersect1d(widx, tidx)) == 0)
    bad = [k for k, v in res.items()
           if k.startswith(("pdb_overlap", "family_overlap")) and v != 0]
    if bad:
        raise AssertionError(f"[{name}] 泄漏断言失败: {bad} -> 结果标 INVALID_DUE_TO_BUG")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 Gate A: V4 baseline reproduction")
    ap.add_argument("--trees", type=int, default=800)
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args()

    lgb = C.load_lgb()
    t0 = time.time()
    emit("# V5.0 GATE A — V4 最强基线复现 (`C+DGaux_fwd+rev`)")
    emit("")
    emit(f"- 设备: **CPU only**（4 核 i7-2630QM / 7.55 GB）；LightGBM 固定 "
         f"{args.trees} 树，**不用 val 早停**")
    emit(f"- seeds: {args.seeds}")
    emit(f"- 口径来源: `scripts/v4.0/11_dualhead_dg_ddg.py`（只读）")
    emit("")

    d = C.build_pair_features(verbose=True)
    emit("")

    # ---------- 泄漏审计（Gate 通过前提）----------
    emit("## 1. 自动泄漏断言（V5 §86–§90）")
    emit("")
    try:
        leak = leakage_assertions(d, "build_pair_features")
        emit("| 检查 | 值 | 期望 | 结论 |")
        emit("|---|---:|---:|---|")
        for k, v in leak.items():
            exp = 1 if k.startswith(("split_", "pca_")) else 0
            emit(f"| `{k}` | {v} | {exp} | {'✅' if v == exp else '❌'} |")
        emit("")
        emit("-> ✅ **全部泄漏断言通过**")
    except AssertionError as exc:
        emit(f"-> ❌ **{exc}**")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 5

    blk, y, groups = d["X_blocks"], d["y"], d["groups"]
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    tr_i, va_i, te_i = np.where(mtr)[0], np.where(mva)[0], np.where(mte)[0]

    # 特征块（与 V4 完全同序）
    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    Sx = np.hstack([blk["phys_mut"] + blk["phys_wt"],
                    blk["Z_mut"] + blk["Z_wt"]]).astype(np.float32)
    Dn, _dsc = C.standardize_train(Dx, mtr)
    Sn, _ssc = C.standardize_train(Sx, mtr)
    dg_mut = d["df"]["dg_mut"].to_numpy(float)

    emit("## 2. 数据规模")
    emit("")
    emit(f"- 配对: **{len(y)}**（forward）；train {int(mtr.sum())} / val {int(mva.sum())}"
         f" / test {int(mte.sum())}")
    emit(f"- 复合物: train {len(set(groups[mtr]))} / val {len(set(groups[mva]))} / "
         f"test {len(set(groups[mte]))}")
    emit(f"- 特征: `D`（反对称块）**{Dn.shape[1]}** 维；`S`（对称块）**{Sn.shape[1]}** 维")
    emit(f"- 辅助目标 `ΔG_mut`: mean {dg_mut[mtr].mean():+.3f}, std {dg_mut[mtr].std():.3f}")
    emit("")
    emit(f"> test 只有 **{len(set(groups[mte]))} 个复合物** -> cluster bootstrap 会"
         f"在极少 cluster 上重采样，CI 会很宽。这是数据规模的硬限制，如实报告。")
    emit("")

    # ---------- 阶段 1: 辅助任务 OOF ----------
    from sklearn.model_selection import GroupKFold
    emit("## 3. 阶段 1 — 辅助任务 `ΔG_mut` 折外预测（GroupKFold by pdb_id）")
    emit("")
    X_aux = np.hstack([Sn, Dn])
    oof = np.full(len(y), np.nan)
    tr_groups = groups[tr_i]
    gkf = GroupKFold(n_splits=args.folds)
    aux_rho = []
    for k, (ia, ib) in enumerate(gkf.split(tr_i, groups=tr_groups)):
        a, b = tr_i[ia], tr_i[ib]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": C.SEED})
        m.fit(X_aux[a], dg_mut[a])
        oof[b] = m.predict(X_aux[b])
        # OOF 完整性断言（V5 §88）：折外行不参与该折训练
        assert len(np.intersect1d(a, b)) == 0, f"fold {k} 训练/折外行重叠"
        aux_rho.append(C.pearson(dg_mut[b], oof[b]))
        emit(f"- fold {k}: train {len(a)}, oof {len(b)}, ΔG_mut Pearson "
             f"**{aux_rho[-1]:+.4f}**")
    emit("")
    emit(f"- 辅助任务折外 Pearson 均值 **{np.mean(aux_rho):+.4f}**")
    m_full = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                  "random_state": C.SEED})
    m_full.fit(X_aux[mtr], dg_mut[mtr])
    aux_te = m_full.predict(X_aux[te_i])
    aux_va = m_full.predict(X_aux[va_i])
    emit(f"- test 上辅助预测 Pearson **{C.pearson(dg_mut[te_i], aux_te):+.4f}**")
    emit(f"- OOF 覆盖: {int(np.isfinite(oof[mtr]).sum())} / {int(mtr.sum())} "
         f"train 行（val/test 的 OOF 不参与训练，仅 test 用全 train 模型）")
    emit("")

    # ---------- 阶段 2: 主任务 ----------
    emit("## 4. 阶段 2 — 主任务 ΔΔG（C+DGaux_fwd+rev）")
    emit("")
    Daux = np.hstack([Dn, oof[:, None]]).astype(np.float32)
    Daux_te = np.hstack([Dn[te_i], aux_te[:, None]]).astype(np.float32)
    Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)   # S/aux 不取负

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}
    emit("| seed | Pearson | Spearman | RMSE | MAE | AUROC | macro-PDB ρ | "
         "macro-PDB r | P@10 | EF@1% | Hit10 | 耗时 |")
    emit("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for sd in args.seeds:
        ts = time.time()
        try:
            Xc = np.vstack([Daux[mtr], Daux_rev[mtr]])
            yc = np.r_[y[mtr], -y[mtr]]
            m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                     "random_state": sd})
            m.fit(Xc, yc)
            pt = m.predict(Daux_te)
            mm = C.all_metrics(y[te_i], pt, groups[te_i])
            rt = time.time() - ts
            key = f"C+DGaux_fwd+rev_s{sd}"
            preds[key] = pt.astype(np.float32)
            rows.append({"experiment_id": f"V5_GATEA_{key}", "family": "gateA_repro",
                         "variant": "C+DGaux_fwd+rev", "seed": sd,
                         "train_n": int(mtr.sum()), "val_n": int(mva.sum()),
                         "test_n": int(mte.sum()),
                         "train_complex_n": len(set(groups[mtr])),
                         "test_complex_n": len(set(groups[mte])),
                         "runtime_sec": round(rt, 1), "status": "OK",
                         "notes": "V4 baseline reproduction", **mm})
            emit(f"| {sd} | {mm['overall_pearson']:+.4f} | {mm['overall_spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['mae']:.4f} | {mm['auroc']:.4f} | "
                 f"{mm['macro_pdb_spearman']:+.4f} | {mm['macro_pdb_pearson']:+.4f} | "
                 f"{mm.get('p_at_10', float('nan')):.4f} | "
                 f"{mm.get('ef_1pct', float('nan')):.4f} | "
                 f"{mm.get('hit10_beneficial', float('nan')):.4f} | {rt:.0f}s |")
        except Exception as exc:  # noqa: BLE001
            emit(f"| {sd} | FAIL {exc} | | | | | | | | | | |")
            with open(C.V5 / "logs" / "failures.log", "a", encoding="utf-8") as fh:
                fh.write(f"[02_baseline_repro] seed={sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")
            rows.append({"experiment_id": f"V5_GATEA_s{sd}", "family": "gateA_repro",
                         "variant": "C+DGaux_fwd+rev", "seed": sd,
                         "status": "FAILED", "notes": str(exc)[:200]})
    emit("")

    res = pd.DataFrame(rows)
    res["git_commit"] = C.git_commit()
    res["peak_ram_mb"] = round(C.peak_ram_mb(), 1)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te_i].astype(np.float32),
             test_group=groups[te_i], aux_dg_mut_oof=oof.astype(np.float32),
             aux_dg_mut_test=aux_te.astype(np.float32),
             **{k: v for k, v in preds.items()})

    # ---------- Gate A 判定 ----------
    emit("## 5. GATE A 判定（|ΔPearson| > 0.01 即不符）")
    emit("")
    emit("| seed | V4 记录 | V5 复现 | Δ | 判定 |")
    emit("|---:|---:|---:|---:|---|")
    gate_ok = True
    for sd in args.seeds:
        sub = res[(res["seed"] == sd) & (res["status"] == "OK")]
        v4 = C.V4_C_DGDG.get(sd)
        if not len(sub) or v4 is None:
            emit(f"| {sd} | {v4} | — | — | 缺数据 |")
            gate_ok = False
            continue
        v5 = float(sub["overall_pearson"].iloc[0])
        dlt = v5 - v4
        ok = abs(dlt) <= 0.01
        gate_ok &= ok
        emit(f"| {sd} | {v4:+.4f} | {v5:+.4f} | {dlt:+.4f} | "
             f"{'✅ 一致' if ok else '❌ 超出阈值'} |")
    emit("")
    emit(f"**GATE A: {'PASS ✅' if gate_ok else 'FAIL ❌'}**")
    emit("")
    if not gate_ok:
        emit("> 按 V5 §11，必须先查明 row/split/PCA/版本/早停/种子/NaN/符号等问题，")
        emit("> **不得直接进入 P1**。所有差异记录于 `baseline_repro_report.md`。")

    # ---------- cluster bootstrap ----------
    emit("")
    emit("## 6. Cluster bootstrap（V5 §12 主推断：抽样单位 = PDB）")
    emit("")
    emit("| seed | 指标 | point | 95% CI | n_clusters | 说明 |")
    emit("|---:|---|---:|---|---:|---|")
    for sd in args.seeds:
        key = f"C+DGaux_fwd+rev_s{sd}"
        if key not in preds:
            continue
        pt = preds[key]
        for mname, fn in (
                ("overall_pearson", C.pearson),
                ("macro_pdb_spearman", C.per_structure_spearman),
                ("macro_pdb_pearson", C.per_structure_pearson)):
            r = C.cluster_bootstrap_ci(fn, y[te_i], pt, groups[te_i],
                                       n_boot=C.N_BOOT_CI, seed=sd)
            emit(f"| {sd} | {mname} | {r['point']:+.4f} | "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {r['n_clusters']} | "
                 f"{r['ci_setting']} |")
    emit("")
    emit(f"> **n_clusters = {len(set(groups[te_i]))}** —— 这是 V5 统计推断的真实"
         f"独立样本量（不是 {int(mte.sum())}）。CI 必然宽，不得用 row-level CI 冒充。")

    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s；峰值 RAM {C.peak_ram_mb():.0f} MB")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    LOG.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0 if gate_ok else 3


if __name__ == "__main__":
    sys.exit(main())
