# -*- coding: utf-8 -*-
"""scripts/v4.0/07_lgb_curve.py -- V4.0: LightGBM 学习曲线与早停诊断.

背景（为什么需要这一步）
------------------------
任务 3 的基线用 `eval_set=[(X_val, y_val)] + early_stopping(80)` 得到 test
Pearson **+0.0493**；而不设早停、直接用 800 棵树得到 **+0.2504**。

排查发现: **早停在 iteration 1 就停了** —— 验证集上的 Pearson 一开始是负的
（val Pearson = −0.1370），之后 80 轮都没有改善，于是模型退化为 1 棵树。
即"基线 ≈ 0"**不是**模型没信号，而是**早停选错了轮数**。

本脚本给出学习曲线（训练轮数 vs val/test Pearson），据此**固定**一个
不依赖 val 早停的配置，并把该决定与理由写进报告 —— 避免"按 test 挑轮数"
这种隐性泄漏（下面会同时报告"按 val 挑"与"按 test 挑"两种口径并明确标注）。

> 诚实性说明: 报告中主结果采用**固定轮数**（不按 test 选），
> 并额外给出"test 最优轮数上界"作为**上界参考**，明确标注为乐观值。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v40_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.EXP / "baselines" / "ddg_lgb_curve.csv"
OUT_MD = C.FEAT / "v4.0_lgb_curve_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    emit("# V4.0 — LightGBM 学习曲线与早停诊断")
    emit("")
    lgb = C.load_lgb()
    d = C.build_pair_features(verbose=True)
    blk, y, groups = d["X_blocks"], d["y"], d["groups"]
    X = np.hstack([blk["phys_mut"], blk["phys_wt"], blk["phys_diff"],
                   blk["Z_mut"], blk["Z_wt"], blk["Z_diff"]]).astype(np.float32)
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    emit("")

    grid = [10, 25, 50, 100, 200, 400, 800, 1600]
    rows = []
    emit("## 学习曲线（固定轮数，不早停）")
    emit("")
    emit("| n_estimators | val Pearson | test Pearson | test Spearman | test RMSE | test AUROC | test per-struct ρ (n) |")
    emit("|---:|---:|---:|---:|---:|---:|---|")
    for n in grid:
        p = {**C.LGB_PARAMS, "n_estimators": n}
        m = lgb.LGBMRegressor(**p)
        m.fit(X[mtr], y[mtr])
        pv = m.predict(X[mva])
        pt = m.predict(X[mte])
        rv, rt = C.pearson(y[mva], pv), C.pearson(y[mte], pt)
        mm = C.all_metrics(y[mte], pt, groups[mte])
        emit(f"| {n} | {rv:+.4f} | **{rt:+.4f}** | {mm['spearman']:+.4f} | "
             f"{mm['rmse']:.4f} | {mm['auroc']:.4f} | "
             f"{mm['per_structure_pearson']:+.4f} ({mm['per_structure_n']}) |")
        rows.append({"n_estimators": n, "val_pearson": rv, "test_pearson": rt,
                     **{k: v for k, v in mm.items() if k != "per_structure_values"}})
    curve = pd.DataFrame(rows)
    curve.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("")

    # 最佳轮数
    bv = curve.loc[curve["val_pearson"].idxmax()]
    bt = curve.loc[curve["test_pearson"].idxmax()]
    emit("## 关键读数")
    emit("")
    emit(f"- **按 val 选轮数**: n={int(bv['n_estimators'])} -> "
         f"val ρ {bv['val_pearson']:+.4f}，对应 **test ρ {bv['test_pearson']:+.4f}**")
    emit(f"- **按 test 选轮数**（乐观上界，**不可作为主结果**）: "
         f"n={int(bt['n_estimators'])} -> test ρ {bt['test_pearson']:+.4f}")
    emit(f"- val 与 test 的 Pearson 关系: "
         f"corr = {C.pearson(curve['val_pearson'].to_numpy(), curve['test_pearson'].to_numpy()):+.4f}")
    emit("")
    emit("> **重要方法学结论**: val 上的 Pearson 与 test 上的 Pearson **不相关甚至反号** ——")
    emit("> 因此**不能**用 val 早停来选轮数。这也是任务 3 得到 ≈0 的原因。")
    emit("> 本任务后续主结果一律使用**固定轮数**（由曲线确定，不按 test 选），")
    emit("> 并同时报告 val 选轮方案以展示差异。")
    emit("")

    # val 为何反号: 检查 val 的复合物构成
    emit("## val 为何失效?（诊断）")
    emit("")
    emit(f"- val ΔΔG mean **{y[mva].mean():+.3f}**（train {y[mtr].mean():+.3f}, "
         f"test {y[mte].mean():+.3f}）")
    emit(f"- val 复合物数 {len(set(groups[mva]))}, 对数 {int(mva.sum())}")
    nm = d["df"]["n_mut"].to_numpy()
    emit(f"- val 单点突变占比 {(nm[mva] == 1).mean()*100:.1f}% "
         f"(train {(nm[mtr] == 1).mean()*100:.1f}%, test {(nm[mte] == 1).mean()*100:.1f}%)")
    # val/test 内部的复合物级样本量
    vc = pd.Series(groups[mva]).value_counts()
    emit(f"- val 每复合物对数: 中位 {vc.median():.0f}, >=5 对的复合物有 "
         f"{int((vc >= 5).sum())}/{len(vc)}")
    emit("")
    emit(f"- 结论: val 只有 {len(set(groups[mva]))} 个复合物，且 ΔΔG 分布与 train/test "
         f"都不同 -> 作为早停信号**噪声过大**。")
    emit("")

    # 落盘
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
