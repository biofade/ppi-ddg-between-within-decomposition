# -*- coding: utf-8 -*-
"""scripts/v2.5/10_ablation.py -- V2.5 任务 5: 综合消融 (cross-target, 统一口径).

**口径统一声明**
----------------
* 全部在 **per-residue 缓存覆盖的同一批 6885 行**上, 使用**固定划分**
  (test 757 / val 965 / train 5163) —— 这就是 Proteinea cross-target 划分;
* **不用 within-target**（V2.4 已证伪: 中位仅 +0.1043, 比 cross-target 低 0.28）;
* 由于本批训练集 (5163) **小于** V2.2 全量 (6482), V2.3.2 的 0.3890 **不是**
  同一训练集下的数字 -> 报告里同时给出「本批重跑的 pooled 基线」与
  「V2.3.2 记录值」两个参照, **不混用**。

特征集 (任务书 5.1)
-------------------
1. `physics_only`
2. `physics+esm2_both`            (pooled, 本批重跑 = 主对照)
3. `physics+esm2_receptor`        (pooled 头部最优)
4. `energy_only`
5. `physics+energy`
6. `physics+energy+esm2_both`
7. `physics+energy+bio+esm2_both` (全量)
8. `physics+energy+bio+esm2_both+perres`
9. `+leaky_dH` 对照 (量化"把实测 ΔG 当特征"的泄漏幅度, 诚实展示)

输出: `features/v2.5_ablation_report.md`、`v2.5_ablation_results.csv`、
`models/v2.5_best_model.pkl`
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_MD = C.FEAT / "v2.5_ablation_report.md"
OUT_CSV = C.FEAT / "v2.5_ablation_results.csv"
OUT_PAIRS = C.FEAT / "v2.5_ablation_pairs.csv"
OUT_PRED = C.FEAT / "v2.5_ablation_predictions.npz"
BEST = C.MODEL_DIR / "v2.5_best_model.pkl"

PCA_DIMS = (32, 64, 128)
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 综合消融 (cross-target)")
    ap.add_argument("--pca-dims", type=int, nargs="*", default=list(PCA_DIMS))
    ap.add_argument("--n-boot", type=int, default=C.N_BOOT_PAIR)
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()

    lgb = C.load_lgb()
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V2.5 任务 5 — 综合消融 (cross-target, 统一口径)")
    emit("")
    emit(f"> 配对 bootstrap {args.n_boot} 次 + BH-FDR; 固定划分 (Proteinea cross-target)")
    emit("")

    data = C.load_aligned(use_perres_pooled=True)
    y = data["y"]; tr, va, te = data["tr"], data["va"], data["te"]
    Xp = data["Xp"]
    dv = data["df_v22"]
    emit("## 1. 数据")
    emit("")
    emit(f"- {len(y)} 行; train {int(tr.sum())} / val {int(va.sum())} / test {int(te.sum())}")
    emit(f"- 物理 {Xp.shape[1]} 维 (V2.2 口径)")

    # ---------- 载入外部特征表, 按 interface_key 对齐 ----------
    def load_ik(path: Path, prefix: str):
        """载入一张外部特征表并按行对齐。

        支持两种键:
        * `interface_key` (能量/热力学表, 3197 行 -> 需 reindex 展开到 6885 行);
        * **行序一致** (保守性/共进化表由 per-residue 缓存直接生成, 6885 行,
          键是 `pdb_ids`)。这类表**必须**逐行同序, 因此直接校验长度后使用。
        """
        if not path.exists():
            emit(f"- !! 缺少 `{path.name}` -> 该特征组跳过")
            return None, [], None
        t = pd.read_csv(path, low_memory=False)
        cols = [c for c in t.columns
                if c not in ("interface_key", "pdb_ids", "pdb_id")
                and t[c].dtype.kind in "fiu"
                and not c.startswith("leaky_")]
        leak = [c for c in t.columns if c.startswith("leaky_")]
        if "interface_key" in t.columns:
            # 外部表按 interface_key 唯一化后, 再 reindex 到本批 6885 行。
            # 注意: 本批的 `dv["interface_key"]` **有重复** (同一链组合下多个突变体),
            # 因此必须先去重源表的索引, 否则 reindex 会报
            # "cannot reindex on an axis with duplicate labels"。
            m = t.drop_duplicates("interface_key", keep="first").set_index("interface_key")
            sub = m.reindex(dv["interface_key"].to_numpy())
            raw = sub[cols]
        else:
            if len(t) != len(dv):
                emit(f"- !! `{path.name}` 行数 {len(t)} != {len(dv)} 且无 "
                     f"`interface_key` -> **无法对齐, 跳过**")
                return None, [], None
            raw = t[cols]
            # 逐行一致性校验 (pdb_ids 必须同序)
            if "pdb_ids" in t.columns:
                same = bool((t["pdb_ids"].astype(str).to_numpy()
                             == data["pdb_ids"]).all())
                if not same:
                    emit(f"- !! `{path.name}` 的 pdb_ids 与缓存**不同序** -> 跳过")
                    return None, [], None
        mat = raw.apply(pd.to_numeric, errors="coerce")
        med = mat[tr].median(numeric_only=True)
        mat = mat.fillna(med).fillna(0.0)
        emit(f"- `{path.name}`: {len(cols)} 个可用特征 "
             f"({len(leak)} 个 leaky 已排除), 缺失率 "
             f"{(raw.isna().any(axis=1)).mean()*100:.1f}%")
        lk = None
        if leak:
            lk = t[leak].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(np.float32)
            if len(lk) != len(dv):
                lk = None            # 无法对齐就不用
        return mat.to_numpy(np.float32), cols, lk

    emit("")
    emit("## 2. 特征块")
    emit("")
    Xe, e_cols, Xleak = load_ik(C.FEAT / "v2.5_energy_features_norm.csv", "energy")
    Xe_raw, e_raw_cols, _ = load_ik(C.FEAT / "v2.5_energy_features.csv", "energy_raw")
    Xc, c_cols, _ = load_ik(C.FEAT / "v2.5_conservation_features.csv", "conserv")
    Xv, v_cols, _ = load_ik(C.FEAT / "v2.5_coevolution_features.csv", "coevol")
    Xt, t_cols, _ = load_ik(C.FEAT / "v2.5_thermo_features.csv", "thermo")

    bio_parts = [x for x in (Xc, Xv, Xt) if x is not None]
    Xbio = np.hstack(bio_parts) if bio_parts else None
    bio_cols = (c_cols or []) + (v_cols or []) + (t_cols or [])
    emit(f"- 生物学/热力学块合计: {Xbio.shape[1] if Xbio is not None else 0} 维")

    # pooled 嵌入 (本批)
    from sklearn.decomposition import PCA
    P = data["pooled"]
    pooled_blocks: dict[str, np.ndarray] = {}
    for dim in args.pca_dims:
        for src, mat in (("both", np.hstack([P["lig"], P["rec"]])),
                         ("receptor", P["rec"])):
            p = PCA(n_components=min(dim, mat.shape[1], int(tr.sum()) - 1),
                    random_state=args.seed)
            p.fit(mat[tr])
            pooled_blocks[f"{src}|{dim}"] = p.transform(mat).astype(np.float32)

    # per-residue mean-pool 块 (任务 1.3 的表示)
    Lm = C.mean_pool_perres(data, "lig"); Rm = C.mean_pool_perres(data, "rec")
    perres_mp = np.hstack([Lm, Rm])
    perres_blocks: dict[str, np.ndarray] = {}
    for dim in args.pca_dims:
        p = PCA(n_components=min(dim, perres_mp.shape[1], int(tr.sum()) - 1),
                random_state=args.seed)
        p.fit(perres_mp[tr])
        perres_blocks[f"both|{dim}"] = p.transform(perres_mp).astype(np.float32)

    # ---------- 定义特征集 ----------
    FS: dict[str, tuple[np.ndarray | None, str]] = {
        "physics_only": (Xp, "V2.2 物理特征 (60)"),
        "energy_raw_only": (Xe_raw, "只用**原始**能量总和 (受界面规模主导)"),
        "energy_norm_only": (Xe, "只用**归一化**能量特征 (主推)"),
        "physics+energy_raw": (np.hstack([Xp, Xe_raw]) if Xe_raw is not None else None,
                               "物理 + 原始能量"),
        "physics+energy": (np.hstack([Xp, Xe]) if Xe is not None else None,
                           "物理 + 归一化能量"),
        "physics+energy_both": (
            np.hstack([Xp, Xe, Xe_raw]) if (Xe is not None and Xe_raw is not None) else None,
            "物理 + 归一化能量 + 原始能量"),
    }
    for dim in args.pca_dims:
        FS[f"physics+esm2_both@{dim}"] = (
            np.hstack([Xp, pooled_blocks[f"both|{dim}"]]), "物理 + pooled(both)")
        FS[f"physics+esm2_receptor@{dim}"] = (
            np.hstack([Xp, pooled_blocks[f"receptor|{dim}"]]), "物理 + pooled(受体)")
    for dim in args.pca_dims:
        if Xe is not None:
            FS[f"physics+energy+esm2_both@{dim}"] = (
                np.hstack([Xp, Xe, pooled_blocks[f"both|{dim}"]]),
                "物理 + 能量 + pooled")
            if Xbio is not None:
                FS[f"physics+energy+bio+esm2_both@{dim}"] = (
                    np.hstack([Xp, Xe, Xbio, pooled_blocks[f"both|{dim}"]]),
                    "物理 + 能量 + 生物学 + pooled (全量)")
                FS[f"physics+energy+bio+esm2_both+perres@{dim}"] = (
                    np.hstack([Xp, Xe, Xbio, pooled_blocks[f"both|{dim}"],
                               perres_blocks[f"both|{dim}"]]),
                    "全量 + per-residue mean-pool")
    # leaky 对照
    if Xleak is not None and Xe is not None:
        FS["LEAKY_physics+energy+bio+esm2_both+dH@128"] = (
            np.hstack([Xp, Xe, Xbio if Xbio is not None else Xe,
                       pooled_blocks["both|128"], Xleak]),
            "**故意加入 leaky_dH (含实测ΔG) 以量化泄漏**")

    emit("")
    emit(f"- 共 **{len(FS)}** 个特征集")

    # ---------- 训练 ----------
    emit("")
    emit("## 3. 结果 (固定 test, n=%d)" % int(te.sum()))
    emit("")
    emit("| 特征集 | 维度 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---:|---|---:|---:|---:|")
    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}
    models: dict[str, object] = {}
    for tag, (X, note) in FS.items():
        if X is None:
            continue
        try:
            m = lgb.LGBMRegressor(**C.LGB_PARAMS)
            m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
                  callbacks=[lgb.early_stopping(50, verbose=False)])
            pred = m.predict(X[te])
            y_s, p_s = -y[te], -pred
            r = {"feature_set": tag, "note": note, "n_features": int(X.shape[1]),
                 "best_iteration": int(getattr(m, "best_iteration_", -1) or -1)}
            for k, fn in mfns.items():
                try:
                    r[k] = float(fn(y_s, p_s))
                except Exception:  # noqa: BLE001
                    r[k] = float("nan")
            for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
                lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, C.N_BOOT_CI, args.seed)
                r[f"{k}_ci_low"], r[f"{k}_ci_high"] = lo, hi
            rows.append(r); preds[tag] = pred.astype(np.float32); models[tag] = m
            emit(f"| {tag} | {r['n_features']} | {r['spearman']:+.4f} | "
                 f"[{r['spearman_ci_low']:+.3f}, {r['spearman_ci_high']:+.3f}] | "
                 f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} |")
        except Exception as exc:  # noqa: BLE001
            emit(f"| {tag} | — | FAIL: {exc} | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[10_ablation] {tag}: {exc}\n{traceback.format_exc()}\n")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te].astype(np.float32),
             test_pdb_id=data["pdb_ids"][te], **preds)

    # ---------- 配对检验 vs 主对照 ----------
    base_tag = next((t for t in preds if t.startswith("physics+esm2_both@128")), None)
    if base_tag is None:
        base_tag = next((t for t in preds if t.startswith("physics+esm2_both")),
                        "physics_only")
    emit("")
    emit(f"## 4. 配对 bootstrap vs 主对照 `{base_tag}`")
    emit("")
    emit("| 特征集 | Δ Spearman | 95% CI | p | BH q | 判定 |")
    emit("|---|---:|---|---:|---:|---|")
    prows, pvals = [], []
    for t in preds:
        if t == base_tag:
            continue
        r = C.paired_bootstrap(mr.spearman_rho, -y[te], -preds[base_tag], -preds[t],
                               args.n_boot, args.seed)
        prows.append({"model": t, "vs": base_tag, **r})
        pvals.append(r["p_two_sided"])
    qs = C.bh_fdr(pvals)
    for r, q in zip(prows, qs):
        r["bh_q"] = q
        v = ("更优" if (q < 0.05 and r["delta"] > 0) else
             "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
        if r["model"].startswith("LEAKY"):
            v = "**泄漏对照**"
        emit(f"| {r['model']} | {r['delta']:+.4f} | "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
             f"{q:.4f} | {v} |")
    pd.DataFrame(prows).to_csv(OUT_PAIRS, index=False, encoding="utf-8")

    # ---------- 特征重要性 ----------
    emit("")
    emit("## 5. 特征重要性 (全量配置 top 15)")
    emit("")
    full_tags = [t for t in models if t.startswith("physics+energy+bio+esm2_both@")]
    if full_tags:
        ft = max(full_tags, key=lambda t: float(
            res.loc[res["feature_set"] == t, "spearman"].iloc[0]))
        m = models[ft]
        cols = ([f"phys:{c}" for c in data["Xp_cols"]]
                + (e_cols or []) + bio_cols
                + [f"pooled_both:{i}" for i in range(pooled_blocks["both|128"].shape[1])])
        imp = np.asarray(m.feature_importances_, float)
        if len(cols) == len(imp):
            top = np.argsort(-imp)[:15]
            emit(f"- 配置 `{ft}` (n_features={len(imp)})")
            emit("")
            emit("| # | 特征 | 重要性 |")
            emit("|---:|---|---:|")
            for k, i in enumerate(top, 1):
                emit(f"| {k} | `{cols[i]}` | {imp[i]:.0f} |")
        else:
            emit(f"(列名数 {len(cols)} != 重要性数 {len(imp)}, 跳过命名)")

    # ---------- 保存最佳 ----------
    if len(res):
        best = res.sort_values("spearman", ascending=False).iloc[0]
        emit("")
        emit(f"## 6. 最优配置")
        emit("")
        emit(f"- `{best['feature_set']}` -> Spearman **{best['spearman']:+.4f}** "
             f"[{best['spearman_ci_low']:+.3f}, {best['spearman_ci_high']:+.3f}], "
             f"NDCG@10 {best['ndcg@10']:.4f}, EF@1% {best['ef@1%']:.4f}, "
             f"P@10 {best['precision@10']:.4f}")
        if not str(best["feature_set"]).startswith("LEAKY"):
            try:
                with open(BEST, "wb") as fh:
                    pickle.dump({"model": models[best["feature_set"]],
                                 "feature_set": best["feature_set"],
                                 "metrics": best.to_dict(),
                                 "note": "V2.5 best (cross-target, perres-subset rows)"},
                                fh)
                emit(f"- 已写出 `{BEST.name}`")
            except Exception as exc:  # noqa: BLE001
                emit(f"- !! 保存失败: {exc}")

    # ---------- 泄漏量化 ----------
    lk = [r for r in rows if str(r["feature_set"]).startswith("LEAKY")]
    if lk:
        emit("")
        emit("## 7. 泄漏量化（诚实展示）")
        emit("")
        non_leak = [r for r in rows if not str(r["feature_set"]).startswith("LEAKY")]
        if non_leak:
            best_ok = max(non_leak, key=lambda r: r["spearman"])
            emit(f"- 不含 leaky 特征的最优: `{best_ok['feature_set']}` "
                 f"{best_ok['spearman']:+.4f}")
            emit(f"- 含 `leaky_dH_from_dG` 的对照: {lk[0]['spearman']:+.4f} "
                 f"(Δ = **{lk[0]['spearman'] - best_ok['spearman']:+.4f}**)")
            emit("- -> 这就是把**实测 ΔG 的线性组合**当特征带来的虚高幅度。")
            emit("  再次说明 `leaky_` 特征**必须排除**在正式结果之外。")

    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
