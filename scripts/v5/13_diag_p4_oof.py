# -*- coding: utf-8 -*-
"""scripts/v5/13_diag_p4_oof.py -- V5.0: 诊断 P4 的 OOF Pearson 异常（只读诊断）.

现象
----
P4 报告 fold 0 的 OOF Pearson = **−0.4234**，但 Gate A 复现的辅助任务 OOF
（同 Fold 数、同 GroupKFold、同特征）在 V4 记录中为 +0.3216（fold 0）。
同一 fold 却符号相反 → 必须先定位原因，**不得**在其上做 conformal 校准
（V5 §123 Bug protocol：先查明再继续）。

检查项
------
1. Gate A 的 OOF 是否真的复现 +0.3216（读 predictions.npz）；
2. P4 与 Gate A 的 OOF 构造差异：特征是否一致、aux 列是否参与了自身预测；
3. 判断哪种是 bug。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def main() -> int:
    emit("# V5.0 — P4 OOF 异常诊断")
    emit("")

    # ---------- 1. Gate A 的 OOF ----------
    g = np.load(C.V5 / "p0_benchmark" / "baseline_repro_predictions.npz",
                allow_pickle=True)
    emit(f"- Gate A npz 数组: {list(g.files)}")
    d = C.build_pair_features(verbose=False)
    y, groups, mtr = d["y"], d["groups"], d["mask_train"]
    tr_i = np.where(mtr)[0]
    df = d["df"]
    dg_mut = df["dg_mut"].to_numpy(float)

    if "aux_dg_mut_oof" in g.files:
        oof_g = g["aux_dg_mut_oof"]
        emit(f"- Gate A OOF 向量长度: {len(oof_g)}（当前特征表 {len(y)}）")
        # Gate A 的 oof 是按**当前** train 行索引存放的
        if len(oof_g) == len(y):
            r_all = C.pearson(dg_mut[mtr], oof_g[mtr])
            cov = np.isfinite(oof_g[mtr]).sum()
            emit(f"- Gate A OOF(train 行) Pearson = **{r_all:+.4f}**"
                 f"（覆盖 {cov}/{int(mtr.sum())}）")
        emit("")
        emit("### Gate A OOF 逐折（按 pdb_id GroupKFold 复算折归属）")
        emit("")
        from sklearn.model_selection import GroupKFold
        for k, (ia, ib) in enumerate(GroupKFold(n_splits=5)
                                     .split(tr_i, groups=groups[tr_i])):
            b = tr_i[ib]
            r = C.pearson(dg_mut[b], oof_g[b])
            emit(f"- fold {k}: n={len(b)}, complexes={len(set(groups[b]))}, "
                 f"Pearson **{r:+.4f}**")
    emit("")

    # ---------- 2. 检查 P4 的 OOF 是否把 aux 列泄漏进自身预测 ----------
    emit("## 关键检查：P4 的 OOF 是否用 aux 列预测了自己")
    emit("")
    emit("P4 的 fold 内训练矩阵构造：")
    emit("")
    emit("```python")
    emit("Xa = np.vstack([Daux[a], np.hstack([-Dn[a], oof_aux[a][:, None]])])")
    emit("```")
    emit("")
    emit("其中 `Daux[a] = [Dn[a], oof_aux[a]]`。因此 **训练特征里包含 `oof_aux`**，")
    emit("而 `oof_aux` 本身是 GroupKFold 的**折外**对该行 ΔG_mut 的预测。")
    emit("")
    emit("对**折外行 b** 而言，`oof_aux[b]` 来自「未包含 b 的模型」→ 对 b 是合法的；")
    emit("但 fold 内**训练行 a** 用的 `oof_aux[a]` 也来自「未包含 a 的模型」→ 也合法。")
    emit("**因此构成上没有泄漏。** 问题在别处。")
    emit("")

    # ---------- 3. 直接复现 P4 的 fold 0 ----------
    emit("## 直接复现 P4 fold 0（与 Gate A 同一 fold 划分与种子）")
    emit("")
    lgb = C.load_lgb()
    from sklearn.model_selection import GroupKFold
    blk = d["X_blocks"]
    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    Sx = np.hstack([blk["phys_mut"] + blk["phys_wt"],
                    blk["Z_mut"] + blk["Z_wt"]]).astype(np.float32)
    Dn, _ = C.standardize_train(Dx, mtr)
    Sn, _ = C.standardize_train(Sx, mtr)
    X_aux = np.hstack([Sn, Dn])

    # 复现 Gate A 的 aux OOF
    oof_aux = np.full(len(y), np.nan)
    for k, (ia, ib) in enumerate(GroupKFold(n_splits=5)
                                 .split(tr_i, groups=groups[tr_i])):
        a, b = tr_i[ia], tr_i[ib]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": 800,
                                 "random_state": C.SEED})
        m.fit(X_aux[a], dg_mut[a])
        oof_aux[b] = m.predict(X_aux[b])
    r_aux = C.pearson(dg_mut[mtr], oof_aux[mtr])
    emit(f"- 本次复算 aux OOF(train) Pearson = **{r_aux:+.4f}**")
    emit("")

    Daux = np.hstack([Dn, oof_aux[:, None]]).astype(np.float32)
    emit("### P4 式的第 1 折主模型（含 aux 列）")
    emit("")
    for k, (ia, ib) in enumerate(GroupKFold(n_splits=5)
                                 .split(tr_i, groups=groups[tr_i])):
        a, b = tr_i[ia], tr_i[ib]
        Xa = np.vstack([Daux[a], np.hstack([-Dn[a], oof_aux[a][:, None]])])
        ya = np.r_[y[a], -y[a]]
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": 800,
                                 "random_state": C.SEED})
        m.fit(Xa, ya)
        pb = m.predict(Daux[b])
        pb_rev = m.predict(np.hstack([-Dn[b], oof_aux[b][:, None]]))
        emit(f"- fold {k}: n={len(b)}, complexes={len(set(groups[b]))}, "
             f"Pearson(正向预测) **{C.pearson(y[b], pb):+.4f}**, "
             f"Pearson(反向预测) **{C.pearson(y[b], pb_rev):+.4f}**")
        if k == 0:
            emit(f"    fold0 y 统计: mean {y[b].mean():+.3f}, std {y[b].std():.3f}, "
                 f"ΔΔG>0 占 {(y[b] > 0).mean() * 100:.1f}%")
            emit(f"    fold0 pred 统计: mean {pb.mean():+.3f}, std {pb.std():.3f}")
    emit("")
    emit("## 判读")
    emit("")
    emit("若各折「正向预测」Pearson 的正负号在不同折间翻转，说明 **OOF 逐折 Pearson")
    emit("在 n≈1000 且组内标签高度偏斜时方差极大**，不能作为校准残差的可靠性指标；")
    emit("此时 conformal 的分位数估计会被离群折拉偏。这不属于「代码 bug」，")
    emit("但属于**方法学风险**，必须在报告中如实标注。")
    emit("")
    OUT = C.V5 / "logs" / "p4_oof_diagnosis.md"
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
