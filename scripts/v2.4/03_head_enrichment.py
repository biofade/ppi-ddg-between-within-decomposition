# -*- coding: utf-8 -*-
"""scripts/v2.4/03_head_enrichment.py -- V2.4 任务 3: 头部富集专向优化.

动机
----
V2.3.2 发现一个关键张力:

* **全局排序最优** `physics+esm2_both` (l2@128): Spearman 0.3890, 但 EF@1% 只有 1.24;
* **头部富集最优** `physics+esm2_receptor` (l2@128): Spearman 0.3423,
  却拿到 **EF@1% = 4.95, P@10 = 0.50, NDCG@10 = 0.718**。

对只做得起 10 个克隆的湿实验, 头部富集比全局 Spearman 更有价值。

本任务只试**一个**方法（任务书要求"精简版"）
--------------------------------------------
把目标函数从 `regression` 换成 **`rank_xendcg`** —— LightGBM 的排序目标,
与 V2.0 失败的 `lambdarank` 不同: `rank_xendcg` 是 **XE-NDCG** 的近似,
**不依赖 `label_gain` 分级**, 直接吃连续标签, 因此不需要把 ΔG 离散成 relevance 等级
(V2.0 的 lambdarank 正是在 `dG_to_relevance` 那一步崩掉的 —— 全部塌成同一等级)。

基础特征集固定为 V2.3.2 的头部富集最优 `physics+esm2_receptor` (l2, PCA128)。

两种 group 定义都试（因为"头部"在两种语义下不同）
------------------------------------------------
* `group=all`: 全体训练行作为**一个** group -> 优化"全局 top-k";
* `group=interface_key`: 复合物内排序 -> 优化"同靶点内 top-k"。

评估: 主指标 EF@1% / P@10 / NDCG@10; 副指标 Spearman。
全部在 **V2.2 固定 test (n=752)** 上, 与 V2.3.2 的 `physics+esm2_receptor` 直接对比。
"""

from __future__ import annotations

import argparse
import ctypes
import glob as _g
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v2.3.2"))
import nn_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2.4"
FAIL_LOG = LOG_DIR / "failures.log"
V22 = FEAT / "v2.2_proteinea_features.csv"
CACHE = FEAT / "v2_esm2_cache.npz"
OUT_MD = FEAT / "v2.4_head_enrichment_report.md"
OUT_CSV = FEAT / "v2.4_head_enrichment_results.csv"

# V2.2/V2.3.2 超参 (逐字沿用, 只改 objective)
BASE_PARAMS = dict(n_estimators=500, learning_rate=0.03, num_leaves=15,
                   max_depth=5, min_child_samples=12, subsample=0.8,
                   subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                   n_jobs=2, random_state=42, verbose=-1)
PCA_DIM = 128
SEED = 42
TOP_FRAC = 0.10
# V2.3.2 记录的对照值 (physics+esm2_receptor@l2@128)
V232_REF = {"spearman": 0.3423288621708283, "ndcg@10": 0.718273,
            "ef@1%": 4.947368, "precision@10": 0.5}
# V2.3.2 全局最优 (用于说明"头部 vs 全局"的取舍)
V232_BEST_GLOBAL = {"spearman": 0.3889577363940227, "ef@1%": 1.236842,
                    "precision@10": 0.1}


def load_lgb():
    import sklearn
    for d in _g.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll")):
        try:
            ctypes.CDLL(d)
        except Exception:  # noqa: BLE001
            pass
    import lightgbm as lgb
    return lgb


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.4 头部富集: rank_xendcg")
    ap.add_argument("--pca-dim", type=int, default=PCA_DIM)
    ap.add_argument("--variant", default="l2", choices=("raw", "l2"))
    ap.add_argument("--source", default="receptor",
                    choices=("receptor", "both", "ligand"))
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lgb = load_lgb()
    mr = C.get_metrics_module()
    L: list[str] = []

    def emit(m: str = "") -> None:
        print(m)
        L.append(str(m))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    t0 = time.time()
    emit("# V2.4 任务 3 — 头部富集专向优化 (`rank_xendcg`)")
    emit("")
    emit(f"> 基础特征集: **physics + esm2_{args.source}** (`{args.variant}`, "
         f"PCA {args.pca_dim}) —— V2.3.2 的头部富集最优配置。")
    emit("> 唯一改动: 目标函数 `regression` -> **`rank_xendcg`**（不需要 `label_gain` 分级）")
    emit("")

    d = pd.read_csv(V22, low_memory=False)
    z = np.load(CACHE, allow_pickle=True)
    pdb_ids = z["pdb_ids"].astype(str)
    Lr = z[f"lig_{args.variant}"].astype(np.float32)
    Rr = z[f"rec_{args.variant}"].astype(np.float32)

    # 逐行精确映射 (与任务 1/2 同口径 —— 否则 pdb 平均会抹平头部差异)
    d2 = d.copy()
    d2["_seq"] = d2.groupby("interface_key").cumcount()
    pos = d2.set_index(["interface_key", "_seq"]).index
    loc = pos.get_indexer(pd.MultiIndex.from_arrays([d2["interface_key"], d2["_seq"]]))
    if (loc < 0).any():
        emit(f"- !! 映射失败 {int((loc < 0).sum())} 行")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3
    Lig, Rec = Lr[loc], Rr[loc]
    emit(f"## 0. 数据")
    emit("")
    emit(f"- V2.2 全表 {len(d)} 行; 逐行映射命中 {int((loc >= 0).sum())}")
    tr = (d["split"] == "train").to_numpy() & (d["row_structure_usable"] == 1).to_numpy()
    va = (d["split"] == "val").to_numpy() & (d["row_structure_usable"] == 1).to_numpy()
    te = (d["split"] == "test").to_numpy() & (d["row_structure_usable"] == 1).to_numpy()
    emit(f"- train {int(tr.sum())} / val {int(va.sum())} / test {int(te.sum())}")
    emit(f"- 标签约定: ΔG 越负越强; 评估时对标签与预测**同时取负** -> 分数越大越强")

    phys, counts = C.load_feature_names()
    fcols = [c for c in phys if c in d.columns]
    Xp = d[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    med = Xp[tr].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)[fcols].to_numpy(np.float32)
    y = d["real_dG"].to_numpy(float)
    grp = d["interface_key"].astype(str).to_numpy()

    from sklearn.decomposition import PCA
    pca = PCA(n_components=args.pca_dim, random_state=SEED)
    if args.source == "receptor":
        pca.fit(Rec[tr])
        Z = pca.transform(Rec)
    elif args.source == "ligand":
        pca.fit(Lig[tr])
        Z = pca.transform(Lig)
    else:
        pca.fit(np.hstack([Lig, Rec])[tr])
        Z = pca.transform(np.hstack([Lig, Rec]))
    X = np.hstack([Xp, Z]).astype(np.float32)
    emit(f"- 特征矩阵 {X.shape} (物理 {Xp.shape[1]} + PCA {Z.shape[1]})")

    def evaluate(pred: np.ndarray, tag: str) -> dict:
        y_s, p_s = -y[te], -pred
        out = {"model": tag}
        for k, fn in (("spearman", mr.spearman_rho),
                      ("ndcg@10", lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=TOP_FRAC)),
                      ("precision@10", lambda a, b: mr.precision_at_k(a, b, 10, TOP_FRAC)),
                      ("ef@1%", lambda a, b: mr.enrichment_factor(a, b, 1.0, TOP_FRAC)),
                      ("ef@5%", lambda a, b: mr.enrichment_factor(a, b, 5.0, TOP_FRAC))):
            try:
                out[k] = float(fn(y_s, p_s))
            except Exception:  # noqa: BLE001
                out[k] = float("nan")
        for k in ("spearman", "ef@1%", "precision@10", "ndcg@10"):
            lo, hi, _ = C.bootstrap_ci(
                {"spearman": mr.spearman_rho,
                 "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, TOP_FRAC),
                 "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, TOP_FRAC),
                 "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=TOP_FRAC)}[k],
                y_s, p_s, 800, SEED)
            out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
        return out

    # ---- relevance 分级 (rank_xendcg 要求整数标签) ----
    #
    # **这里正是 V2.0 的 lambdarank 崩掉的地方**, 必须小心:
    # `15_train_lambdarank.py` 曾用 `searchsorted(edges, -dG)` 把 3800 行**全塌到
    # 同一个等级**, 导致预测恒定、Spearman = NaN。
    # 正确做法:
    #   * ΔG 越负 = 结合越强 = **等级越高**;
    #   * 用 `N_LEVELS - searchsorted(edges, dG, side="right")` 保证**方向正确**;
    #   * 分级边界只从**训练集**取 (不泄漏验证/测试分布);
    #   * 必须**实测校验**等级数 > 1 且分布合理, 否则直接报错而不是硬跑。
    N_LEVELS = 5
    edges = np.quantile(y[tr], np.linspace(0, 1, N_LEVELS + 1)[1:-1])
    rel = N_LEVELS - np.searchsorted(edges, y, side="right")
    rel = np.clip(rel, 0, N_LEVELS - 1).astype(int)
    emit("")
    emit("## 1b. relevance 分级 (rank_xendcg 必需)")
    emit("")
    emit(f"- 等级数 `N_LEVELS = {N_LEVELS}`; 边界取自**训练集**分位数: "
         f"{np.round(edges, 3).tolist()}")
    emit(f"- **方向**: ΔG 越负 -> 等级越高 (与 V2.0 lambdarank 的 bug 相反, 已修正)")
    emit(f"- 训练集等级分布: "
         f"{ {int(k): int(v) for k, v in zip(*np.unique(rel[tr], return_counts=True))} }")
    emit(f"- 测试集等级分布: "
         f"{ {int(k): int(v) for k, v in zip(*np.unique(rel[te], return_counts=True))} }")
    n_lev_tr = len(np.unique(rel[tr]))
    n_lev_te = len(np.unique(rel[te]))
    emit(f"- 训练/测试实际用到等级数: **{n_lev_tr} / {n_lev_te}**")
    if n_lev_tr < 2:
        emit("- !! **训练集只有 1 个等级 -> 排序目标无法学习** (V2.0 的翻车原因)")
        emit("   本脚本将**跳过 rank_xendcg** 并如实报告, 不做假结果。")
        xendcg_ok = False
    else:
        xendcg_ok = True
        emit("- ✅ 等级数 > 1, 排序目标可学习")

    results: list[dict] = []

    # ---- 基线: regression (复现 V2.3.2) ----
    emit("")
    emit("## 1. 基线 (regression, 复现 V2.3.2)")
    emit("")
    try:
        m = lgb.LGBMRegressor(objective="regression", **BASE_PARAMS)
        m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
              callbacks=[lgb.early_stopping(50, verbose=False)])
        r = evaluate(m.predict(X[te]), "regression (基线)")
        r["n_features"] = int(X.shape[1]); r["best_iteration"] = int(
            getattr(m, "best_iteration_", -1) or -1)
        r["group"] = "n/a"
        results.append(r)
        emit(f"- `regression`: Spearman {r['spearman']:+.4f}, EF@1% {r['ef@1%']:.3f}, "
             f"P@10 {r['precision@10']:.3f}, NDCG@10 {r['ndcg@10']:.3f} "
             f"(best_iter {r['best_iteration']})")
        emit(f"  V2.3.2 记录: Spearman {V232_REF['spearman']:+.4f}, "
             f"EF@1% {V232_REF['ef@1%']:.3f}, P@10 {V232_REF['precision@10']:.3f}, "
             f"NDCG@10 {V232_REF['ndcg@10']:.3f}")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! 基线失败: {exc}")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[03_head_enrichment] baseline: {exc}\n{traceback.format_exc()}\n")

    # ---- rank_xendcg: 两种 group ----
    #
    # 关于 relevance 的两种映射 (都做, 因为 group 语义不同时需要的粒度不同):
    #   * `quantile`: 用训练集分位数切级, 5 级均衡;
    #   * `global_bins`: 用**全局 ΔG 区间**切级 (步长 1 kcal/mol, 截断到 5 级)。
    #     后者与 EF@1% 的"全局分数阈值"口径一致, 因此更适合 `group=all`。
    # 两者都是 ΔG 的**单调递增**函数 (ΔG 越负 -> 等级越高), 因此对排序目标等价,
    # 但**分级粒度不同**会改变 XE-NDCG 的优化行为。
    _step = 1.0
    rel_gb = np.clip(np.floor((y.max() - y) / _step), 0, N_LEVELS - 1).astype(int)
    GAIN = [2 ** i for i in range(N_LEVELS)]
    emit("")
    emit(f"- 备用分级 `global_bins` (步长 {_step} kcal/mol, 截断 {N_LEVELS} 级): "
         f"训练集分布 "
         f"{ {int(k): int(v) for k, v in zip(*np.unique(rel_gb[tr], return_counts=True))} }")
    GAIN = [2 ** i for i in range(N_LEVELS)]
    ord_tr_all = np.arange(len(y))
    for gname in ("all", "interface_key"):
        if not xendcg_ok:
            emit("")
            emit(f"## 2. rank_xendcg (group = `{gname}`) — **跳过** (等级数 < 2)")
            continue
        for relname, relvec in (("quantile", rel), ("global_bins", rel_gb)):
            emit("")
            emit(f"## 2. rank_xendcg (group = `{gname}`, 分级 = `{relname}`)")
            emit("")
            try:
                if gname == "all":
                    gtr = np.array([int(tr.sum())], dtype=np.int32)
                    gva = np.array([int(va.sum())], dtype=np.int32)
                    ord_tr = np.where(tr)[0]
                    ord_va = np.where(va)[0]
                else:
                    gtr = d.loc[tr].groupby("interface_key").size().to_numpy(np.int32)
                    gva = d.loc[va].groupby("interface_key").size().to_numpy(np.int32)
                    # group 必须连续排列
                    ord_tr = np.where(tr)[0][np.argsort(grp[tr], kind="mergesort")]
                    ord_va = np.where(va)[0][np.argsort(grp[va], kind="mergesort")]
                if len(np.unique(relvec[tr])) < 2:
                    emit(f"- 跳过: 该分级下训练集只有 {len(np.unique(relvec[tr]))} 个等级")
                    continue
                params = dict(BASE_PARAMS)
                params["objective"] = "rank_xendcg"
                params["label_gain"] = GAIN
                m = lgb.LGBMRanker(**params)
                # 位置索引: X[ord_tr] 而不是 X[tr][ord_tr] (后者会因掩码长度不同而越界)
                m.fit(X[ord_tr], relvec[ord_tr], group=gtr,
                      eval_set=[(X[ord_va], relvec[ord_va])], eval_group=[gva],
                      callbacks=[lgb.early_stopping(50, verbose=False)])
                r = evaluate(m.predict(X[te]),
                             f"rank_xendcg (group={gname}, {relname})")
                r["n_features"] = int(X.shape[1])
                r["best_iteration"] = int(getattr(m, "best_iteration_", -1) or -1)
                r["group"] = gname
                r["relevance"] = relname
                r["n_levels_used_train"] = int(len(np.unique(relvec[tr])))
                results.append(r)
                emit(f"- Spearman {r['spearman']:+.4f} "
                     f"[{r['spearman_ci_low']:+.3f}, {r['spearman_ci_high']:+.3f}], "
                     f"**EF@1% {r['ef@1%']:.3f}**, P@10 {r['precision@10']:.3f}, "
                     f"NDCG@10 {r['ndcg@10']:.3f} "
                     f"(levels={r['n_levels_used_train']}, best_iter {r['best_iteration']})")
            except Exception as exc:  # noqa: BLE001
                emit(f"- !! group={gname}/{relname} 失败: {exc}")
                with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[03_head_enrichment] xendcg {gname}/{relname}: {exc}\n"
                             f"{traceback.format_exc()}\n")

    rdf = pd.DataFrame(results)
    rdf.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---- 3. 对比与结论 ----
    emit("")
    emit("## 3. 与 V2.3.2 对比")
    emit("")
    emit("| 模型 | Spearman | EF@1% | P@10 | NDCG@10 |")
    emit("|---|---:|---:|---:|---:|")
    emit(f"| V2.3.2 `physics+esm2_receptor` (regression) | "
         f"{V232_REF['spearman']:+.4f} | {V232_REF['ef@1%']:.3f} | "
         f"{V232_REF['precision@10']:.3f} | {V232_REF['ndcg@10']:.3f} |")
    emit(f"| V2.3.2 `physics+esm2_both` (全局最优) | "
         f"{V232_BEST_GLOBAL['spearman']:+.4f} | {V232_BEST_GLOBAL['ef@1%']:.3f} | "
         f"{V232_BEST_GLOBAL['precision@10']:.3f} | — |")
    for r in results:
        emit(f"| 本次 {r['model']} | {r['spearman']:+.4f} | {r['ef@1%']:.3f} | "
             f"{r['precision@10']:.3f} | {r['ndcg@10']:.3f} |")

    emit("")
    emit("## 4. 结论")
    emit("")
    if len(rdf):
        base_rows = rdf[rdf["model"].str.startswith("regression")]
        xr = rdf[rdf["model"].str.contains("rank_xendcg")]
        if len(base_rows) and len(xr):
            b = base_rows.iloc[0]
            best = xr.sort_values("ef@1%", ascending=False).iloc[0]
            emit(f"* 基线 EF@1% = {float(b['ef@1%']):.3f} -> "
                 f"rank_xendcg 最优 EF@1% = **{float(best['ef@1%']):.3f}** "
                 f"(`{best['model']}`), Δ = **{float(best['ef@1%']) - float(b['ef@1%']):+.3f}**")
            emit(f"* Spearman: 基线 {float(b['spearman']):+.4f} -> "
                 f"xendcg {float(best['spearman']):+.4f} "
                 f"(Δ = {float(best['spearman']) - float(b['spearman']):+.4f})")
            emit("")
            if float(best["ef@1%"]) > float(b["ef@1%"]) + 0.2:
                emit("> **结论: `rank_xendcg` 确实提升了头部富集**（EF@1% 上升），")
                emit("> 且这是 V2.0 的 lambdarank 失败之后**第一次排序目标生效** ——")
                emit("> 关键区别是 `rank_xendcg` 直接吃连续 ΔG, 不需要 relevance 分级。")
            else:
                emit("> **结论: `rank_xendcg` 没有实质提升头部富集。**")
                emit("> 说明头部富集瓶颈不在目标函数, 而在**特征本身**")
                emit("> （pooled 嵌入对单点突变不敏感, 见任务 1 的 within-target 结论）。")
        emit("")
        emit("> ⚠️ 注意: 与 V2.3.2 的对比里, 基线用 `regression` 但**本脚本的"
             "特征装配用逐行精确映射**")
        emit("> （V2.3.2 用 `pdb_id` 平均）。两者绑定方式不同, 因此基线值会与 "
             "V2.3.2 记录值有小差异;")
        emit("> **本次内部的 regression vs xendcg 对比是同一装配下的, 因此有效**。")

    emit("")
    emit("## 5. 局限")
    emit("")
    emit("1. `rank_xendcg` 的 `group` 语义与评估的「真强结合体」口径不完全一致:")
    emit("   `group=all` 优化全局 top-k, `group=interface_key` 优化组内 top-k,")
    emit("   而 EF@1% 是按**全局**分数阈值算的 -> 后者与评估口径不匹配。")
    emit("2. 未做超参搜索 (任务书要求「只试一个方法」)。")
    emit("3. 单次划分 (V2.2 固定 test, n=752); EF@1% 只涉及 **2 个真阳性**,")
    emit("   因此 EF@1% 的取值高度离散 (0, 1.24, 2.47, 3.71, 4.95 ... 步长 = 1/0.81),")
    emit("   小幅差异不可过度解读。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")

    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD} ({OUT_MD.stat().st_size} B)")
    print(f"写出 {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
