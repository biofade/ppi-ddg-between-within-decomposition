# -*- coding: utf-8 -*-
"""scripts/v4.0/11_dualhead_dg_ddg.py -- V4.0 任务 4 创新 2: ΔG + ΔΔG 双任务.

文献依据
--------
* **Twin Peaks** [p7]：`L = MSE(ΔG) + λ·MSE(ΔΔG)` 双头联合；消融显示关闭双头
  后 ΔΔG Pearson **0.485→0.558** 而 ΔG **0.638→0.532** [p11] —— 两任务间存在
  **可调 trade-off**。最优 λ = **0.718** [p9]。
* **Pythia-PPI** [p8][p18]：`L = α·L_affinity + (1−α)·L_stability`，**α=0.8** 最优；
  加辅助任务使 per-structure Pearson **0.4462→0.4778**。
* **SELFprot** [p5]：去掉 binding 任务 → kcat RMSD **+154%**；去掉 kinetics →
  Ki RMSD **+106%** —— 任务间迁移的定量证据。

与创新 1 的关系（复合，不重复）
-------------------------------
创新 1 得到 `C_D_fwd+rev`（仅反对称块 + 前向/反向增广）。本创新在它之上**再加一个
ΔG 辅助任务**，检验"多任务"是否在 ΔΔG 上额外有效。因此两者是**叠加关系**：
`C` ⊂ `C + ΔG_aux`。

实现方式（为什么用 stacking 而不是双头 NN）
-------------------------------------------
任务书建议"双头"。但本机实测（见 `10_dualhead_symmetric.py` 的说明）：**MLP 在这批
表格特征上根本学不动**（三种容量配置全部第 0 轮停止改善，test Pearson ≤ +0.10）。
V3.1 与 V2.3.2 也是同一结论。因此忠实保留**双任务语义**，改用 CPU 上真正有效的
**两阶段 stacking**：

1. **阶段 1（辅助任务）**：用 `[S, D]` 预测**绝对** `ΔG_mut`。
   恒等式 `ΔG_mut = ((x_mut+x_wt) + (x_mut−x_wt))/2 = (S + D)/2`，
   故 `[S, D]` 是预测 ΔG_mut 的**充分且正确**的输入。
   用 **GroupKFold（按 pdb_id 分组）** 生成**折外**预测，避免泄漏。
2. **阶段 2（主任务）**：把折外 `ΔG_mut` 预测作为**额外一列**加到创新 1 的 `D` 上，
   预测 ΔΔG。

> 诚实标注：这是"**辅助任务提供额外监督信号**"的 stacking 实现，
> **不是**神经网络双头。报告与结论中一律如此表述，不冒充双头。

评估: Pearson / Spearman / RMSE / MAE / per-structure Pearson / AUROC；
与创新 1 的 `C_D_fwd+rev` **同种子配对** bootstrap + BH-FDR。
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v40_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.EXP / "innovation2_results.csv"
OUT_PRED = C.EXP / "innovation2_predictions.npz"
OUT_PAIR = C.EXP / "innovation2_pairs.csv"
OUT_MD = C.FEAT / "v4.0_innovation2_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 创新 2: ΔG+ΔΔG 双任务")
    ap.add_argument("--trees", type=int, default=800)
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    ap.add_argument("--folds", type=int, default=5, help="GroupKFold 折数")
    args = ap.parse_args()

    lgb = C.load_lgb()
    from sklearn.model_selection import GroupKFold
    t0 = time.time()

    emit("# V4.0 任务 4 创新 2 — ΔG + ΔΔG 双任务（Twin Peaks / Pythia-PPI 思想）")
    emit("")
    emit("> 文献依据: Twin Peaks 双头 `L=MSE(ΔG)+λ·MSE(ΔΔG)`，λ=0.718 最优，")
    emit("> 消融显示 ΔΔG 0.485→0.558 / ΔG 0.638→0.532 的 trade-off [Twin Peaks p7,p9,p11]；")
    emit("> Pythia-PPI α=0.8，per-structure Pearson 0.4462→0.4778 [p8,p18]；")
    emit("> SELFprot 任务迁移定量证据 [p5]。")
    emit("")
    emit("> **实现方式诚实标注**: 任务书建议 NN 双头；但本机实测 MLP 在这批表格特征上")
    emit("> **学不动**（见 `10_dualhead_symmetric.py` 的实测数字：三种配置全部第 0 轮停止，")
    emit("> test Pearson ≤ +0.10，RMSE ≈ 常数基线）。因此保留**双任务语义**，")
    emit("> 改用 CPU 上有效的**两阶段 stacking**（辅助任务折外预测作为额外特征）。")
    emit("> 这**不是** NN 双头，结论中不冒充。")
    emit("")

    d = C.build_pair_features(verbose=True)
    blk, y, groups = d["X_blocks"], d["y"], d["groups"]
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    # 对称/反对称分解
    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    Sx = np.hstack([blk["phys_mut"] + blk["phys_wt"],
                    blk["Z_mut"] + blk["Z_wt"]]).astype(np.float32)
    Dn, dsc = C.standardize_train(Dx, mtr)
    Sn, ssc = C.standardize_train(Sx, mtr)
    dg_mut = d["df"]["dg_mut"].to_numpy(float)
    emit("")
    emit(f"- `D`（反对称块）{Dn.shape[1]} 维; `S`（对称块）{Sn.shape[1]} 维")
    emit(f"- 辅助目标 `ΔG_mut`: mean {dg_mut[mtr].mean():+.3f}, "
         f"std {dg_mut[mtr].std():.3f}")
    emit(f"- 恒等式 `ΔG_mut = (S + D)/2` -> 用 `[S, D]` 预测它是**正确**输入")
    emit("")

    # ---------- 阶段 1: 辅助任务折外预测 ----------
    emit("## 1. 阶段 1 — 辅助任务 `ΔG_mut` 的折外预测（GroupKFold by pdb_id）")
    emit("")
    X_aux = np.hstack([Sn, Dn])
    oof = np.full(len(y), np.nan)
    tr_idx = np.where(mtr)[0]
    gkf = GroupKFold(n_splits=args.folds)
    tr_groups = groups[tr_idx]
    aux_scores = []
    for k, (i_a, i_b) in enumerate(gkf.split(tr_idx, groups=tr_groups)):
        a, b = tr_idx[i_a], tr_idx[i_b]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees})
        m.fit(X_aux[a], dg_mut[a])
        oof[b] = m.predict(X_aux[b])
        aux_scores.append(C.pearson(dg_mut[b], oof[b]))
        emit(f"- fold {k}: n_train {len(a)}, n_val {len(b)}, "
             f"ΔG_mut Pearson **{aux_scores[-1]:+.4f}**")
    emit("")
    emit(f"- 辅助任务折外 Pearson: 均值 **{np.mean(aux_scores):+.4f}** "
         f"（范围 [{min(aux_scores):+.4f}, {max(aux_scores):+.4f}]）")
    emit("> 说明辅助任务**本身是可学的**（ΔG 绝对量比 ΔΔG 好预测，与 V2.3.2 的 "
         "cross-target 0.3890 同量级）。")
    emit("")
    # test 的辅助预测: 用全 train 训练
    m_full = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees})
    m_full.fit(X_aux[mtr], dg_mut[mtr])
    aux_te = m_full.predict(X_aux[mte])
    aux_va = m_full.predict(X_aux[mva])
    emit(f"- test 上辅助预测的 Pearson **{C.pearson(dg_mut[mte], aux_te):+.4f}**")
    emit("")

    # ---------- 阶段 2: 主任务 ----------
    emit("## 2. 阶段 2 — 主任务 ΔΔG")
    emit("")
    emit("| 变体 | seed | Pearson | 95% CI | Spearman | RMSE | MAE | per-struct ρ (n) | AUROC |")
    emit("|---|---:|---:|---|---:|---:|---:|---|---:|")
    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    def fit_pred(Xtr, ytr, seed, Xte):
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,
                                 "random_state": seed})
        m.fit(Xtr, ytr)
        return m.predict(Xte)

    for sd in args.seeds:
        # 参照: 创新 1 的 C = 仅 D + fwd/rev 增广
        Xc = np.vstack([Dn[mtr], -Dn[mtr]])
        yc = np.r_[y[mtr], -y[mtr]]
        try:
            p = fit_pred(Xc, yc, sd, Dn[mte])
            mm = C.all_metrics(y[mte], p, groups[mte])
            lo, hi, _ = C.bootstrap_ci(C.pearson, y[mte], p)
            mm["pearson_ci_low"], mm["pearson_ci_high"] = lo, hi
            results.append({"model": "C_D_fwd+rev", "seed": sd, **mm})
            preds[f"C_D_fwd+rev_s{sd}"] = p.astype(np.float32)
            emit(f"| `C_D_fwd+rev` | {sd} | **{mm['pearson']:+.4f}** | "
                 f"[{lo:+.3f}, {hi:+.3f}] | {mm['spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['mae']:.4f} | "
                 f"{mm['per_structure_pearson']:+.4f} ({mm['per_structure_n']}) | "
                 f"{mm['auroc']:.4f} |")
        except Exception as exc:  # noqa: BLE001
            emit(f"| C_D_fwd+rev | {sd} | FAIL {exc} | | | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[11_dualhead_dg_ddg] C {sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")

    # 主变体: D + 辅助折外预测 (+ fwd/rev 增广)
    Daux = np.hstack([Dn, oof[:, None]]).astype(np.float32)
    Daux_te = np.hstack([Dn[mte], aux_te[:, None]]).astype(np.float32)
    # 反向增广: D 取负, 辅助 ΔG_mut 不变(它不是反对称量!)
    Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)
    for sd in args.seeds:
        try:
            Xc2 = np.vstack([Daux[mtr], Daux_rev[mtr]])
            yc2 = np.r_[y[mtr], -y[mtr]]
            p = fit_pred(Xc2, yc2, sd, Daux_te)
            mm = C.all_metrics(y[mte], p, groups[mte])
            lo, hi, _ = C.bootstrap_ci(C.pearson, y[mte], p)
            mm["pearson_ci_low"], mm["pearson_ci_high"] = lo, hi
            results.append({"model": "C+DGaux_fwd+rev", "seed": sd, **mm})
            preds[f"C+DGaux_fwd+rev_s{sd}"] = p.astype(np.float32)
            emit(f"| `C+DGaux_fwd+rev` | {sd} | **{mm['pearson']:+.4f}** | "
                 f"[{lo:+.3f}, {hi:+.3f}] | {mm['spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['mae']:.4f} | "
                 f"{mm['per_structure_pearson']:+.4f} ({mm['per_structure_n']}) | "
                 f"{mm['auroc']:.4f} |")
        except Exception as exc:  # noqa: BLE001
            emit(f"| C+DGaux_fwd+rev | {sd} | FAIL {exc} | | | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[11_dualhead_dg_ddg] C+aux {sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[mte].astype(np.float32), test_group=groups[mte],
             aux_dg_mut_oof=oof.astype(np.float32),
             aux_dg_mut_test=aux_te.astype(np.float32),
             **{k: v for k, v in preds.items()})

    # ---------- 配对检验（同种子）----------
    emit("")
    emit("## 3. 配对 bootstrap vs 创新 1 的 `C_D_fwd+rev`（同种子，BH-FDR）")
    emit("")
    emit("> 这是**复合对比**: 检验\"在创新 1 之上再加 ΔG 辅助任务\"是否额外有效。")
    emit("")
    emit("| 对比 | Δ Pearson | 95% CI | p | BH q | 判定 |")
    emit("|---|---:|---|---:|---:|---|")
    prows, pvals = [], []
    for sd in args.seeds:
        a = f"C_D_fwd+rev_s{sd}"
        k = f"C+DGaux_fwd+rev_s{sd}"
        if a in preds and k in preds:
            rr = C.paired_bootstrap(C.pearson, y[mte], preds[a], preds[k])
            prows.append({"comparison": f"{k} − {a}", **rr})
            pvals.append(rr["p_two_sided"])
    qs = C.bh_fdr(pvals)
    win = 0
    for r, q in zip(prows, qs):
        r["bh_q"] = q
        ok = (q < 0.05 and r["delta"] > 0.02)
        if ok:
            win += 1
        v = ("**胜出**" if ok else
             "更优(未达阈值)" if (q < 0.05 and r["delta"] > 0) else
             "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
        emit(f"| {r['comparison']} | {r['delta']:+.4f} | "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
             f"{q:.4f} | {v} |")
    if prows:
        pd.DataFrame(prows).to_csv(OUT_PAIR, index=False, encoding="utf-8")
    emit("")
    emit(f"**胜出计数（BH q<0.05 且 ΔPearson>0.02）: {win} / {len(prows)}**")
    if prows:
        dd = np.array([r["delta"] for r in prows])
        emit(f"- Δ Pearson: 均值 {dd.mean():+.4f}, 范围 [{dd.min():+.4f}, {dd.max():+.4f}]")

    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
