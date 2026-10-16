# -*- coding: utf-8 -*-
"""scripts/v4.0/05_diagnose_ddg.py -- V4.0 诊断: ΔΔG 信号是否存在于现有特征中?

为什么必须做这一步
------------------
任务 3 的四个基线全部 ≈ 0（Pearson 0.035–0.049，BH q=0.72，与常数预测无区别）。
在宣布"ΔΔG 难做"或"任务定义有问题"**之前**，必须先排除以下**我自己的**可能错误：

1. **特征构造错误** —— 差分块方向搞反？WT/mutant 行映射错位？
2. **嵌入本身对单点突变不敏感** —— pooled ESM-2 只换一个残基时几乎不变
   （SAMAffinity [p16] 明确指出这是 pooling 的固有缺陷）。
3. **标签噪声** —— ΔΔG 的量级（std 2.29）是否远大于可用信号？
4. **划分太严** —— 家族级隔离是否把同类界面全部排除，导致泛化不可能？

诊断设计（全部只用训练标签，不碰 test 做选择）
---------------------------------------------
* **D1 常数基线**：预测 train 均值，算 test Pearson/RMSE（下限参照）。
* **D2 嵌入差分的模**：`||Z_mut − Z_wt||₂` 与 |ΔΔG| / ΔΔG 的相关性。
  若这个简单量与 ΔΔG 都不相关，说明 pooled 嵌入对突变**不敏感**。
* **D3 物理特征差分**：`phys_diff` 与 ΔΔG 的相关（取最大 |r| 的列）。
* **D4 同复合物内 ΔΔG 是否可排序**（train 内部，留一复合物）：
  检验"同一复合物内不同突变体的 ΔΔG"与特征的关系 —— 这是 ΔΔG 的本义。
* **D5 标签噪声上界**：同一复合物内**同一 ΔΔG 值**重复出现次数（重复测量）。
* **D6 划分严格度**：test 的复合物在 train 中是否有**任何**同源（用序列相似度看），
  以及 test 的 ΔΔG 分布与 train 的差异（协变量偏移）。
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

OUT_MD = C.FEAT / "v4.0_ddg_diagnosis.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    emit("# V4.0 诊断 — ΔΔG 信号是否存在于现有特征中?")
    emit("")
    emit("> 触发原因: 任务 3 的四个基线全部 ≈ 0（Pearson 0.035–0.049，BH q=0.72）。")
    emit("> 在归因于\"ΔΔG 难做\"之前，先排除**我自己的特征构造错误**。")
    emit("")

    d = C.build_pair_features(verbose=True)
    b, y = d["X_blocks"], d["y"]
    groups = d["groups"]
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    emit("")

    # ---------- D1 常数基线 ----------
    emit("## D1 常数基线（下限参照）")
    emit("")
    mu = float(y[mtr].mean())
    const = np.full(int(mte.sum()), mu)
    m = C.all_metrics(y[mte], const, groups[mte])
    emit(f"- 预测值恒为 train 均值 {mu:+.4f}")
    emit(f"- test Pearson = **{m['pearson']:+.4f}**（常数预测的 Pearson 定义为 0/NaN）")
    emit(f"- test RMSE = **{m['rmse']:.4f}**, MAE = {m['mae']:.4f}")
    emit(f"- test AUROC = **{m['auroc']:.4f}**（常数应为 0.5）")
    emit("")
    emit(f"> 对比: B1 的 RMSE 2.60 vs 常数基线 RMSE {m['rmse']:.2f} -> "
         f"{'B1 略优' if 2.60 < m['rmse'] else '**B1 并不优于常数预测**'}")
    emit("")

    # ---------- D2 嵌入差分的模 ----------
    emit("## D2 pooled 嵌入对突变是否敏感?（核心诊断）")
    emit("")
    zd = np.linalg.norm(b["Z_diff"], axis=1)
    emit(f"- `||Z_mut − Z_wt||₂`: 中位 {np.median(zd):.5f}, 均值 {zd.mean():.5f}, "
         f"最大 {zd.max():.5f}")
    emit(f"- `Z` 本身的模: 中位 {np.median(np.linalg.norm(b['Z_mut'],axis=1)):.4f}")
    rel = zd / (np.linalg.norm(b["Z_mut"], axis=1) + 1e-12)
    emit(f"- **相对变化** `||ΔZ||/||Z||`: 中位 {np.median(rel)*100:.3f}% "
         f"（均值 {rel.mean()*100:.3f}%）")
    emit("")
    emit(f"- 单点突变对的 `||ΔZ||` 中位: "
         f"{np.median(zd[d['df']['n_mut'].to_numpy() == 1]):.5f}")
    emit(f"- 多点突变对的 `||ΔZ||` 中位: "
         f"{np.median(zd[d['df']['n_mut'].to_numpy() > 1]):.5f}")
    emit("")
    emit("| 关系 | Pearson | Spearman |")
    emit("|---|---:|---:|")
    for lab, yy in (("ΔΔG", y), ("|ΔΔG|", np.abs(y))):
        emit(f"| `\\|\\|ΔZ\\|\\|` vs {lab} | {C.pearson(zd, yy):+.4f} | "
             f"{C.spearman(zd, yy):+.4f} |")
    emit("")
    emit(f"- `\\|\\|ΔZ\\|\\|` vs ΔΔG 用 train 拟合的线性回归在 test 上的 Pearson:")
    # 单特征线性回归（train 拟合）
    zt, zs = zd[mtr], zd[mte]
    A = np.vstack([zt, np.ones_like(zt)]).T
    coef = np.linalg.lstsq(A, y[mtr], rcond=None)[0]
    pred = coef[0] * zs + coef[1]
    emit(f"  **{C.pearson(y[mte], pred):+.4f}**（斜率 {coef[0]:+.2f}, 截距 {coef[1]:+.2f}）")
    emit("")

    # ---------- D3 物理特征差分 ----------
    emit("## D3 物理特征差分与 ΔΔG 的相关（取最强列）")
    emit("")
    pd_ = b["phys_diff"]
    rs = np.array([C.pearson(pd_[:, i], y) if np.std(pd_[:, i]) > 0 else 0.0
                   for i in range(pd_.shape[1])])
    order = np.argsort(-np.abs(rs))[:10]
    phys_names, _ = C.load_phys_names()
    emit("| 物理差分特征 | Pearson vs ΔΔG |")
    emit("|---|---:|")
    for i in order:
        nm = phys_names[i] if i < len(phys_names) else f"col{i}"
        emit(f"| {nm} | {rs[i]:+.4f} |")
    emit("")
    emit(f"- 最大 |r| = **{np.abs(rs).max():.4f}**（{len(rs)} 个差分特征）")
    emit("")

    # ---------- D4 同复合物内可排序性 ----------
    emit("## D4 同一复合物内 ΔΔG 是否可排序?（ΔΔG 的本义）")
    emit("")
    # 只在 train 内，按复合物算组内 |ΔΔG| 的离散度
    tr_df = d["df"][mtr].copy()
    tr_df["_y"] = y[mtr]
    tr_df["_zd"] = zd[mtr]
    g = tr_df.groupby("pdb_id")["_y"]
    sz = g.size()
    sd = g.std()
    big = sz[sz >= 5].index
    emit(f"- train 内 >=5 对的复合物: **{len(big)}**")
    if len(big):
        emit(f"- 这些复合物内 ΔΔG 的 std 中位: **{sd[big].median():.4f}**")
        emit(f"- 组内 `||ΔZ||` 与 ΔΔG 的平均组内 Spearman:")
        vals = []
        for pid in big:
            s = tr_df[tr_df["pdb_id"] == pid]
            v = C.spearman(s["_zd"].to_numpy(), s["_y"].to_numpy())
            if v == v:
                vals.append(v)
        emit(f"  **{np.mean(vals) if vals else float('nan'):+.4f}** "
             f"（{len(vals)} 个复合物）")
        emit("")
        emit("> 若组内也 ≈ 0，说明 pooled 嵌入**连同一复合物内**的相对大小都排不出来 ——")
        emit("> 这是 pooling 丢失位点信息的直接证据（与 SAMAffinity [p16] 的诊断一致）。")
    emit("")

    # ---------- D5 标签重复/噪声 ----------
    emit("## D5 标签噪声与重复")
    emit("")
    tr = d["df"][mtr]
    dup = tr.groupby(["pdb_id", "mutation"]).size()
    emit(f"- train 中 `(pdb_id, mutation)` 组合数 {len(dup)}；重复(>1)的组合 "
         f"**{int((dup > 1).sum())}**")
    if (dup > 1).any():
        sub = tr[tr.set_index(["pdb_id", "mutation"]).index.isin(dup[dup > 1].index)]
        gg = sub.groupby(["pdb_id", "mutation"])["ddg"].std().dropna()
        emit(f"- 重复组合的 ΔΔG std 中位: **{gg.median():.4f}** kcal/mol"
             f"（最大 {gg.max():.4f}）")
        emit("> 这就是**同一突变被多次测量**的离散度 = 标签噪声的直接估计。")
    emit(f"- test ΔΔG std = {y[mte].std():.4f}")
    emit("")

    # ---------- D6 划分严格度 / 协变量偏移 ----------
    emit("## D6 划分严格度与协变量偏移")
    emit("")
    emit(f"- train ΔΔG: mean {y[mtr].mean():+.3f}, std {y[mtr].std():.3f}")
    emit(f"- test  ΔΔG: mean {y[mte].mean():+.3f}, std {y[mte].std():.3f}")
    emit(f"- **均值差 {y[mte].mean() - y[mtr].mean():+.3f} kcal/mol**"
         f"（若常数预测，RMSE 至少是这个量级）")
    emit("")
    emit(f"- test 的 `n_mut` 分布 vs train:")
    nm = d["df"]["n_mut"].to_numpy()
    emit(f"  - train 单点占比 {(nm[mtr] == 1).mean()*100:.1f}%, "
         f"test {(nm[mte] == 1).mean()*100:.1f}%")
    emit("")
    # 特征分布的偏移: 用 train 拟合的简单分类器区分 train/test（AUC 越高偏移越大）
    emit("- 协变量偏移检测（train vs test 可区分性）:")
    try:
        lgb = C.load_lgb()
        Xa = np.vstack([b["Z_diff"][mtr], b["Z_diff"][mte]])
        ya = np.r_[np.zeros(int(mtr.sum())), np.ones(int(mte.sum()))]
        clf = lgb.LGBMClassifier(n_estimators=100, num_leaves=15, n_jobs=2,
                                 random_state=C.SEED, verbose=-1)
        clf.fit(Xa, ya)
        from sklearn.model_selection import cross_val_predict
        pr = cross_val_predict(clf, Xa, ya, cv=3, method="predict_proba")[:, 1]
        emit(f"  - 用 `Z_diff` 区分 train/test 的 AUC ≈ **{C.auroc(ya, pr):.4f}**"
             f"（0.5=无偏移，>0.7=明显偏移）")
    except Exception as exc:  # noqa: BLE001
        emit(f"  - 偏移检测失败: {exc}")
    emit("")

    # ---------- 总结 ----------
    emit("## 总结与归因")
    emit("")
    emit("（见各节数字；结论在最终报告中给出，此处只列事实。）")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
