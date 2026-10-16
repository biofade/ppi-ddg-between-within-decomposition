# -*- coding: utf-8 -*-
"""scripts/v2.5/02_perres_meanpool.py -- V2.5 任务 1.2/1.3: pooled 基线重跑 + per-residue mean-pool.

**为什么必须重跑 pooled 基线**
------------------------------
per-residue 缓存的 `train` 只有 **5163** 行, 而 V2.2 可用行是 **6482** 行
(`logs/v2.5/join_check.log`)。训练集不同 -> 直接拿 V2.3.2 的 0.3890 与
per-residue 结果比是**不公平**的。

因此本脚本在**同一批 6885 行、同一固定划分**上重跑三个 LightGBM 变体:

* `pooled_both`       = `physics + esm2_both`(链级 pooled, l2)  <- V2.3.2 口径的对照重跑
* `pooled_receptor`   = `physics + esm2_receptor`(链级 pooled, l2) <- V2.3.2 头部最优重跑
* `perres_meanpool`   = `physics + 界面残基 meanpool(lig+rec)`  <- 任务 1.3
* `perres_meanpool_only` = 只用界面残基 meanpool (无物理)

PCA (pooled 与 perres 都做) **只在训练集 fit**。
所有模型在**同一个 test (n=757)** 上评估 -> 可逐行配对检验。
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
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v2.5_meanpool_results.csv"
OUT_PRED = C.FEAT / "v2.5_meanpool_predictions.npz"
OUT_MD = C.FEAT / "v2.5_meanpool_report.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 pooled 基线重跑 + perres mean-pool")
    ap.add_argument("--pca-dims", type=int, nargs="*", default=[32, 64, 128])
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()

    lgb = C.load_lgb()
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V2.5 任务 1.2/1.3 — pooled 基线重跑 + per-residue mean-pool")
    emit("")
    emit("> 目的: 在**同一批行 / 同一固定划分**上建立可比基线, 并测试"
         "「界面残基 mean-pool」是否优于「整链 mean-pool」。")
    emit("")

    data = C.load_aligned(use_perres_pooled=True)
    y = data["y"]
    tr, va, te = data["tr"], data["va"], data["te"]
    Xp = data["Xp"]
    emit("## 1. 数据（同一批行）")
    emit("")
    emit(f"- 行数 **{len(y)}**; splits: "
         f"train {int(tr.sum())} / val {int(va.sum())} / test {int(te.sum())}")
    emit(f"- 物理特征 {Xp.shape[1]} 维")
    emit(f"- **注意**: 本批的 train 是 {int(tr.sum())} 行, 小于 V2.2 全量 "
         f"{int((pd.read_csv(C.V22, low_memory=False)['split'] == 'train').sum())} 行")
    emit(f"  -> 因此 pooled 基线**必须重跑**, 不能直接引用 V2.3.2 的 {C.V232['physics+esm2_both@l2@128']:.4f}")

    # 一致性复核: 紧凑缓存的 pooled 是否等于原 pooled 缓存的对应行
    emit("")
    emit("### 1b. pooled 嵌入一致性复核")
    try:
        z = np.load(C.PERRES_CACHE, allow_pickle=True)
        pr_row = np.asarray(z["pooled_row"])
        po = C.load_pooled_direct()
        # 原缓存 lig_l2 是 float32; 紧凑缓存存了 fp16 -> 比较时容忍 fp16 精度
        a = np.abs(po["lig"][pr_row] - data["pooled"]["lig"]).max()
        b = np.abs(po["rec"][pr_row] - data["pooled"]["rec"]).max()
        emit(f"- 原始 pooled vs 紧凑缓存 pooled 的最大绝对差: "
             f"lig {a:.2e}, rec {b:.2e} (fp16 往返, 应 < 1e-3)")
        assert a < 1e-3 and b < 1e-3, "pooled 嵌入不一致"
        emit("- ✅ 一致 (仅 fp16 往返误差)")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! 复核失败: {exc}")
        with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[02_perres_meanpool] pooled check: {exc}\n{traceback.format_exc()}\n")

    from sklearn.decomposition import PCA

    # ---- 构造各特征块 ----
    emit("")
    emit("## 2. 特征块")
    emit("")
    blocks: dict[str, np.ndarray] = {}
    for dim in args.pca_dims:
        for variant, key in (("l2", "pooled"),):
            Lg = data["pooled"]["lig"]
            Rc = data["pooled"]["rec"]
            for src, mat in (("ligand", Lg), ("receptor", Rc),
                             ("both", np.hstack([Lg, Rc]))):
                p = PCA(n_components=min(dim, mat.shape[1], int(tr.sum()) - 1),
                        random_state=args.seed)
                p.fit(mat[tr])                       # 只在训练集 fit
                blocks[f"pooled|{src}|{dim}"] = p.transform(mat).astype(np.float32)
        for src in ("ligand", "receptor", "both"):
            Lm = C.mean_pool_perres(data, "lig")
            Rm = C.mean_pool_perres(data, "rec")
            mat = {"ligand": Lm, "receptor": Rm,
                   "both": np.hstack([Lm, Rm])}[src]
            p = PCA(n_components=min(dim, mat.shape[1], int(tr.sum()) - 1),
                    random_state=args.seed)
            p.fit(mat[tr])
            blocks[f"perres|{src}|{dim}"] = p.transform(mat).astype(np.float32)
    emit(f"- 生成 {len(blocks)} 个 PCA 块 (pooled / perres × ligand/receptor/both "
         f"× {args.pca_dims})")

    # ---- 训练与评估 ----
    emit("")
    emit("## 3. 训练与评估（固定 test, n=%d）" % int(te.sum()))
    emit("")
    emit("| 变体 | 特征数 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---:|---|---:|---:|---:|")

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    def run(tag: str, X: np.ndarray, note: str) -> None:
        try:
            m = lgb.LGBMRegressor(**C.LGB_PARAMS)
            m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
                  callbacks=[lgb.early_stopping(50, verbose=False)])
            pred = m.predict(X[te])
            y_s, p_s = -y[te], -pred
            row = {"variant": tag, "note": note, "n_features": int(X.shape[1]),
                   "best_iteration": int(getattr(m, "best_iteration_", -1) or -1)}
            for k, fn in mfns.items():
                try:
                    row[k] = float(fn(y_s, p_s))
                except Exception:  # noqa: BLE001
                    row[k] = float("nan")
            for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
                lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, C.N_BOOT_CI, args.seed)
                row[f"{k}_ci_low"], row[f"{k}_ci_high"] = lo, hi
            results.append(row)
            preds[tag] = pred
            emit(f"| {tag} | {row['n_features']} | {row['spearman']:+.4f} | "
                 f"[{row['spearman_ci_low']:+.3f}, {row['spearman_ci_high']:+.3f}] | "
                 f"{row['ndcg@10']:.4f} | {row['ef@1%']:.4f} | "
                 f"{row['precision@10']:.4f} |")
        except Exception as exc:  # noqa: BLE001
            emit(f"| {tag} | — | FAIL: {exc} | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[02_perres_meanpool] {tag}: {exc}\n{traceback.format_exc()}\n")

    # 物理单独 (V2.2 基线重跑)
    run("physics_only", Xp, "V2.2 物理特征, 本批行重跑")
    # pooled
    best_pooled_tag, best_pooled_rho = None, -9.0
    for dim in args.pca_dims:
        run(f"pooled_both@{dim}", np.hstack([Xp, blocks[f"pooled|both|{dim}"]]),
            "链级 pooled (l2) + 物理")
        run(f"pooled_receptor@{dim}",
            np.hstack([Xp, blocks[f"pooled|receptor|{dim}"]]),
            "链级 pooled 受体侧 + 物理 (V2.3.2 头部最优口径)")
    # per-residue mean-pool
    for dim in args.pca_dims:
        run(f"perres_meanpool_both@{dim}",
            np.hstack([Xp, blocks[f"perres|both|{dim}"]]),
            "**界面残基** mean-pool (lig+rec) + 物理")
        run(f"perres_meanpool_only@{dim}", blocks[f"perres|both|{dim}"],
            "只用界面残基 mean-pool (无物理)")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te].astype(np.float32),
             test_pdb_id=data["pdb_ids"][te], **{k: v.astype(np.float32)
                                                 for k, v in preds.items()})

    # ---- 关键对比: perres meanpool vs pooled (配对) ----
    emit("")
    emit("## 4. 关键配对检验: per-residue mean-pool vs 链级 pooled")
    emit("")
    base_tag = None
    for pref in ("pooled_both@128", "pooled_both@64", "pooled_both@32"):
        if pref in preds:
            base_tag = pref
            break
    emit(f"- 基准: `{base_tag}`")
    emit("")
    emit("| 对比 | Δ Spearman | 95% CI | p | BH q | 判定 |")
    emit("|---|---:|---|---:|---:|---|")
    comps = [t for t in preds if t.startswith("perres_meanpool")]
    prows, pvals = [], []
    for t in comps:
        r = C.paired_bootstrap(mr.spearman_rho, -y[te], -preds[base_tag], -preds[t],
                               C.N_BOOT_PAIR, args.seed)
        prows.append({"comparison": f"{t} − {base_tag}", **r})
        pvals.append(r["p_two_sided"])
    qs = C.bh_fdr(pvals)
    for r, q in zip(prows, qs):
        verd = "更优" if (q < 0.05 and r["delta"] > 0) else (
            "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
        emit(f"| {r['comparison']} | {r['delta']:+.4f} | "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
             f"{q:.4f} | **{verd}** |")
        r["bh_q"] = q
    pd.DataFrame(prows).to_csv(C.FEAT / "v2.5_meanpool_pairs.csv", index=False,
                               encoding="utf-8")

    # ---- 结论 ----
    emit("")
    emit("## 5. 结论")
    emit("")
    if len(res):
        bp = res[res["variant"].str.startswith("pooled_both")]
        pp = res[res["variant"].str.startswith("perres_meanpool_both")]
        if len(bp) and len(pp):
            emit(f"- 链级 pooled 最优: `{bp.loc[bp['spearman'].idxmax(), 'variant']}` "
                 f"Spearman = {bp['spearman'].max():+.4f}")
            emit(f"- 界面残基 mean-pool 最优: "
                 f"`{pp.loc[pp['spearman'].idxmax(), 'variant']}` "
                 f"Spearman = {pp['spearman'].max():+.4f}")
            d = float(pp["spearman"].max() - bp["spearman"].max())
            emit(f"- Δ = **{d:+.4f}**")
            emit("")
            emit("> 解读: 「界面残基 mean-pool」与「整链 mean-pool」的唯一区别是")
            emit("> **平均的范围**(界面残基 vs 整条链)。若两者接近, 说明")
            emit("> **平均操作本身就是信息瓶颈**, 而不只是「平均了哪些残基」。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    print(f"写出 {OUT_CSV} | {OUT_PRED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
