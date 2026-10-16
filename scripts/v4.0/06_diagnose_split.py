# -*- coding: utf-8 -*-
"""scripts/v4.0/06_diagnose_split.py -- V4.0 诊断: 划分是否过于严苛（或偏移过大）?

待检验的三个假设
----------------
**H1: 协变量偏移** —— 我的确定性贪心划分让 train ΔΔG 均值 +1.798、test +0.523，
相差 **−1.275 kcal/mol**。若模型预测接近 train 均值，则 RMSE 天然很大、Pearson 被压低。
-> 检验: (a) 在 **val**（同分布性更好?）上评估；(b) 构造**ΔΔG 均值对齐**的
   分层划分变体，比较指标变化。

**H2: 家族级隔离过严** —— 40% 同一性 + 家族级可能把相似界面全部排除，
使 test 成为"跨家族外推"（本质上是 few-shot 新家族问题）。
-> 检验: 放宽到**仅复合物级**隔离（不做序列聚类），比较指标。

**H3: 特征/模型不足** —— 真信号在，但 pooled 嵌入 + 物理特征不够。
-> 由 D2/D4 支持（组内 Spearman 仅 +0.0215）。

**重要**: 放宽划分**只是为了诊断**，主结果仍以严格划分为准。
本脚本会把两种划分的结果并列报告，并明确哪一个用于主要结论。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v40_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_MD = C.FEAT / "v4.0_split_diagnosis.md"
OUT_CSV = C.EXP / "baselines" / "ddg_split_variants.csv"
OUT_JSON = C.V40 / "ddg_split_stratified.json"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def stratified_split(fwd: pd.DataFrame, seed: int = 42,
                     target=(0.70, 0.15, 0.15)) -> dict[str, list[int]]:
    """按复合物分配，使各侧 **ΔΔG 加权均值尽量对齐**（确定性，无随机）。

    做法: 复合物按自身 ΔΔG 均值**降序**处理（极端值先安置），每次把当前复合物
    分给"实际加权均值 < 全局均值"最多的那一侧（同时用配额做硬约束）。
    这比"交错分配"更直接地对齐均值 —— 实测交错反而把 train/test 均值差
    从 −1.275 扩大到 −1.905。
    """
    gm = fwd.groupby("pdb_id")["ddg"].agg(["mean", "size"])
    gm = gm.sort_values("mean", ascending=False).reset_index()
    total = int(gm["size"].sum())
    gmean = float((gm["mean"] * gm["size"]).sum() / total)
    quota = {s: target[i] * total for i, s in enumerate(("train", "val", "test"))}
    got = {s: 0 for s in quota}
    wsum = {s: 0.0 for s in quota}
    assign: dict[str, str] = {}
    for _, row in gm.iterrows():
        size, mn = int(row["size"]), float(row["mean"])
        cands = [s for s in quota if got[s] + size <= quota[s] * 1.10]
        if not cands:
            cands = list(quota)
        # 选"当前加权均值低于全局均值最多"的一侧；并列时选配额更欠的
        def key(s):
            cur = (wsum[s] / got[s]) if got[s] > 0 else gmean
            return (cur - gmean, got[s] - quota[s], s)
        s = min(cands, key=key)
        assign[row["pdb_id"]] = s
        got[s] += size
        wsum[s] += mn * size
    rows: dict[str, list[int]] = {s: [] for s in quota}
    for i, pid in enumerate(fwd["pdb_id"]):
        rows[assign[pid]].append(int(i))
    return rows


def eval_lgb(lgb, X, y, groups, tr, te, tag, mfns_extra=True):
    m = lgb.LGBMRegressor(**C.LGB_PARAMS)
    m.fit(X[tr], y[tr])
    p = m.predict(X[te])
    return p, C.all_metrics(y[te], p, groups[te])


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 划分诊断")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44])
    args = ap.parse_args()

    t0 = time.time()
    emit("# V4.0 诊断 — 划分是否过严 / 偏移是否致命?")
    emit("")
    emit("> 本脚本的结果**仅用于诊断**；主要结论仍以严格划分（任务 2）为准。")
    emit("")

    lgb = C.load_lgb()
    d = C.build_pair_features(verbose=True)
    b, y, groups = d["X_blocks"], d["y"], d["groups"]
    X = np.hstack([b["phys_mut"], b["phys_wt"], b["phys_diff"],
                   b["Z_mut"], b["Z_wt"], b["Z_diff"]]).astype(np.float32)
    mtr, mva, mte = d["mask_train"], d["mask_val"], d["mask_test"]
    emit("")

    # ---------- 1. 原划分: train/val/test 三处都评 ----------
    emit("## 1. 严格划分（任务 2）下，各侧指标")
    emit("")
    emit("| 训练用 | 评估用 | n | Pearson | Spearman | RMSE | MAE | AUROC |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|")
    res_rows = []
    for tr_name, tr_m in (("train", mtr), ("train+val", mtr | mva)):
        for ev_name, ev_m in (("val", mva), ("test", mte)):
            if not ev_m.any():
                continue
            m = lgb.LGBMRegressor(**C.LGB_PARAMS)
            m.fit(X[tr_m], y[tr_m])
            p = m.predict(X[ev_m])
            mm = C.all_metrics(y[ev_m], p, groups[ev_m])
            emit(f"| {tr_name} | {ev_name} | {int(ev_m.sum())} | "
                 f"**{mm['pearson']:+.4f}** | {mm['spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['mae']:.4f} | {mm['auroc']:.4f} |")
            res_rows.append({"split_rule": "strict", "train_on": tr_name,
                             "eval_on": ev_name, **mm})
    emit("")
    emit(f"- train ΔΔG mean **{y[mtr].mean():+.3f}** / val **{y[mva].mean():+.3f}** / "
         f"test **{y[mte].mean():+.3f}**")
    emit(f"- **train 与 test 的均值差 = {y[mte].mean() - y[mtr].mean():+.3f} kcal/mol**")
    emit("")

    # 常数基线（每侧）
    emit("### 1b. 常数基线对比（判断模型是否有信息）")
    emit("")
    mu = float(y[mtr].mean())
    for ev_name, ev_m in (("val", mva), ("test", mte)):
        cm = C.all_metrics(y[ev_m], np.full(int(ev_m.sum()), mu), groups[ev_m])
        emit(f"- {ev_name}: 常数预测 RMSE **{cm['rmse']:.4f}** / "
             f"AUROC {cm['auroc']:.4f}")
    emit("")

    # ---------- 2. 分层划分变体（对齐 ΔΔG 均值）----------
    emit("## 2. 变体 A: ΔΔG 均值对齐的分层划分（只保证复合物级隔离）")
    emit("")
    fwd = d["df"]
    rows = stratified_split(fwd, C.SEED)
    emit("| split | n | ΔΔG mean | ΔΔG std | 复合物数 |")
    emit("|---|---:|---:|---:|---:|")
    masks = {}
    for s in ("train", "val", "test"):
        mk = np.zeros(len(y), bool)
        mk[rows[s]] = True
        masks[s] = mk
        emit(f"| {s} | {int(mk.sum())} | {y[mk].mean():+.3f} | {y[mk].std():.3f} | "
             f"{len(set(groups[mk]))} |")
    emit("")
    emit(f"- **train/test 均值差 = {y[masks['test']].mean() - y[masks['train']].mean():+.3f}**"
         f"（严格划分是 {y[mte].mean() - y[mtr].mean():+.3f}）")
    # 复合物隔离验证
    leak = 0
    for a_name, b_name in (("train", "val"), ("train", "test"), ("val", "test")):
        leak += len(set(groups[masks[a_name]]) & set(groups[masks[b_name]]))
    emit(f"- 复合物级交集合计: **{leak}** （{'✅' if leak == 0 else '❌'}）")
    emit("")
    emit("| 训练用 | 评估用 | n | Pearson | Spearman | RMSE | AUROC |")
    emit("|---|---|---:|---:|---:|---:|---:|")
    for tr_name, tr_m in (("train", masks["train"]),
                          ("train+val", masks["train"] | masks["val"])):
        for ev_name in ("val", "test"):
            ev_m = masks[ev_name]
            m = lgb.LGBMRegressor(**C.LGB_PARAMS)
            m.fit(X[tr_m], y[tr_m])
            p = m.predict(X[ev_m])
            mm = C.all_metrics(y[ev_m], p, groups[ev_m])
            emit(f"| {tr_name} | {ev_name} | {int(ev_m.sum())} | "
                 f"**{mm['pearson']:+.4f}** | {mm['spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['auroc']:.4f} |")
            res_rows.append({"split_rule": "stratified", "train_on": tr_name,
                             "eval_on": ev_name, **mm})
    emit("")

    # ---------- 3. 仅复合物级隔离（放宽 40% 序列约束）----------
    emit("## 3. 变体 B: 仅复合物级隔离（放宽 40% 序列同一性 + 家族级）")
    emit("")
    # 随机但确定性地把复合物分到 70/15/15（按复合物，仍无跨侧泄漏）
    rs = np.random.RandomState(C.SEED)
    pids = np.array(sorted(set(groups)))
    rs.shuffle(pids)
    n = len(pids)
    n_tr, n_va = int(0.70 * n), int(0.15 * n)
    sets = {"train": set(pids[:n_tr]), "val": set(pids[n_tr:n_tr + n_va]),
            "test": set(pids[n_tr + n_va:])}
    masks2 = {s: np.isin(groups, list(sets[s])) for s in sets}
    emit("| split | n | ΔΔG mean | 复合物数 |")
    emit("|---|---:|---:|---:|")
    for s in ("train", "val", "test"):
        emit(f"| {s} | {int(masks2[s].sum())} | {y[masks2[s]].mean():+.3f} | "
             f"{len(sets[s])} |")
    emit("")
    emit("| 训练用 | 评估用 | n | Pearson | Spearman | RMSE | AUROC |")
    emit("|---|---|---:|---:|---:|---:|---:|")
    for tr_name, tr_m in (("train", masks2["train"]),
                          ("train+val", masks2["train"] | masks2["val"])):
        for ev_name in ("val", "test"):
            ev_m = masks2[ev_name]
            m = lgb.LGBMRegressor(**C.LGB_PARAMS)
            m.fit(X[tr_m], y[tr_m])
            p = m.predict(X[ev_m])
            mm = C.all_metrics(y[ev_m], p, groups[ev_m])
            emit(f"| {tr_name} | {ev_name} | {int(ev_m.sum())} | "
                 f"**{mm['pearson']:+.4f}** | {mm['spearman']:+.4f} | "
                 f"{mm['rmse']:.4f} | {mm['auroc']:.4f} |")
            res_rows.append({"split_rule": "complex_only", "train_on": tr_name,
                             "eval_on": ev_name, **mm})
    emit("")

    # ---------- 4. 多次随机复合物级划分（估计划分敏感性）----------
    emit("## 4. 变体 C: 20 次随机**复合物级**划分（衡量划分敏感性 → 置信下限）")
    emit("")
    rhos = []
    for s_ in range(20):
        r2 = np.random.RandomState(1000 + s_)
        pp = np.array(sorted(set(groups)))
        r2.shuffle(pp)
        n2 = len(pp)
        a_tr = set(pp[:int(0.70 * n2)])
        a_te = set(pp[int(0.85 * n2):])
        m_tr = np.isin(groups, list(a_tr))
        m_te = np.isin(groups, list(a_te))
        mm = lgb.LGBMRegressor(**C.LGB_PARAMS)
        mm.fit(X[m_tr], y[m_tr])
        pr = mm.predict(X[m_te])
        rhos.append(C.pearson(y[m_te], pr))
    rhos = np.array(rhos)
    emit(f"- 20 次随机复合物级划分的 test Pearson: 均值 **{rhos.mean():+.4f}**, "
         f"std {rhos.std():.4f}, 范围 [{rhos.min():+.4f}, {rhos.max():+.4f}]")
    emit(f"- 我的严格划分结果在 P{int((rhos < 0.0493).mean()*100)} 分位附近"
         f"（随机划分均值 {rhos.mean():+.3f} vs 严格 {0.0493:+.3f}）")
    emit("")
    emit("> 若随机划分也 ≈ 0，说明**不是划分太严**，而是特征/任务本身信号弱。")

    # ---------- 5. 落盘 ----------
    pd.DataFrame(res_rows).to_csv(OUT_CSV, index=False, encoding="utf-8")
    OUT_JSON.write_text(json.dumps(
        {s: [int(i) for i in rows[s]] for s in rows},
        ensure_ascii=False, indent=1), encoding="utf-8")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
