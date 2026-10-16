# -*- coding: utf-8 -*-
"""scripts/v2.4/04_binding_control.py -- V2.4 关键对照: 嵌入绑定方式对指标的影响.

为什么必须做这个对照
--------------------
任务 1/2 里 `row_exact` (逐行精确映射) 给出的数字**明显高于** V2.3.2 主结果:

| 设置 | 数字 |
|---|---|
| V2.3.2 主结果 (`pdb_id` 平均, 固定 V2.2 test) | +0.3890 |
| 任务 2 `row_exact` (GroupKFold, 全体) | +0.6734 |
| 任务 2 `row_exact` (GroupKFold, SAbDab) | +0.7993 |

这个差距太大, 必须先查清是**方法差异**还是**乐观偏差**, 否则所有后续结论都建在沙上。

已知事实 (来自诊断与实测):
* `pdb_id` 平均会把同一结构的多个突变体的序列嵌入**平均成一个** ->
  4987 个 SKEMPI 行只剩 197 个唯一嵌入; 同一复合物内**所有行特征相同**;
* `row_exact` 让每行拿到**自己那条突变序列**的嵌入 (1808 个唯一嵌入);
* V2.2 的划分在 `pdb_id` 层面是**完全不重叠的** (train/test 交集 = 0)
  -> 因此 `row_exact` 下, 测试复合物的**自身序列**只有在训练时出现过才会被"认出来"。
  同一 PDB 的其它突变体若在训练集里, 其序列与测试行**高度相似** (单点突变),
  模型可以借它**推断出这是哪个复合物**, 从而回忆起该复合物的亲和力量级 ->
  **这是 target-identification, 不是 binding-affinity 泛化**。

本脚本做三个对照 (全部在同一固定划分上, 只改一处)
-------------------------------------------------
* **A. 固定 V2.2 test + `pdb_id` 平均** = V2.3.2 口径 (should reproduce 0.3890)
* **B. 固定 V2.2 test + `row_exact`**     = 只改绑定方式
* **C. 固定 V2.2 test + `row_exact` 但把测试复合物在训练集中的
  同 pdb 行全部删除** = 严格 target-disjoint 的 row_exact

C 是关键: 若 B 明显高于 C, 则证明 **B 的好成绩来自同 pdb 的训练行 (泄漏)**;
若 C ≈ A, 则 row_exact 本身没有额外收益。
"""

from __future__ import annotations

import argparse
import ctypes
import gc
import glob as _g
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
OUT_MD = FEAT / "v2.4_binding_control_report.md"
OUT_CSV = FEAT / "v2.4_binding_control_results.csv"

PARAMS = dict(objective="regression", n_estimators=500, learning_rate=0.03,
              num_leaves=15, max_depth=5, min_child_samples=12, subsample=0.8,
              subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
              n_jobs=2, random_state=42, verbose=-1)
PCA_DIM = 128
SEED = 42
V232_RHO = 0.3889577363940227


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
    ap = argparse.ArgumentParser(description="V2.4 嵌入绑定方式对照")
    ap.add_argument("--pca-dim", type=int, default=PCA_DIM)
    ap.add_argument("--variant", default="l2", choices=("raw", "l2"))
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
    emit("# V2.4 关键对照 — 嵌入绑定方式 (`pdb_id` 平均 vs 逐行精确映射)")
    emit("")
    emit("> 三个设置**只改一处**, 全部在同一固定划分 (V2.2 的 train/val/test) 上评估。")
    emit("")

    d = pd.read_csv(V22, low_memory=False)
    z = np.load(CACHE, allow_pickle=True)
    pdb_ids = z["pdb_ids"].astype(str)
    Lr = z[f"lig_{args.variant}"].astype(np.float32)
    Rr = z[f"rec_{args.variant}"].astype(np.float32)

    emit("## 0. 前提校验")
    emit("")
    emit(f"- V2.2 全表 {len(d)} 行 | 缓存 {len(pdb_ids)} 行 | "
         f"行序一致: **{bool((d['pdb_id'].astype(str).to_numpy() == pdb_ids).all())}**")
    tr = (d["split"] == "train").to_numpy()
    va = (d["split"] == "val").to_numpy()
    te = (d["split"] == "test").to_numpy()
    emit(f"- 划分: train {int(tr.sum())} / val {int(va.sum())} / test {int(te.sum())}")
    pdb_tr = set(d.loc[tr, "pdb_id"].astype(str))
    pdb_te = set(d.loc[te, "pdb_id"].astype(str))
    emit(f"- **train/test 的 pdb_id 交集: {len(pdb_tr & pdb_te)}** "
         f"(0 = 靶点层面完全不重叠)")
    emit(f"- 但**同一 pdb 的其它突变体**可能在训练集里: "
         f"test 中有 {len(pdb_te & pdb_tr)} 个 pdb 同时出现在 train "
         f"(结构层面重叠度 = {len(pdb_te & pdb_tr)}/{len(pdb_te)})")

    phys, counts = C.load_feature_names()
    fcols = [c for c in phys if c in d.columns]
    Xp = d[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    med = Xp[tr].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)[fcols].to_numpy(np.float32)
    y = d["real_dG"].to_numpy(float)
    usable = d["row_structure_usable"].to_numpy() == 1

    # ---- 绑定 A: pdb_id 平均 ----
    emb = pd.DataFrame({"pdb_id": pdb_ids, "lig": list(Lr), "rec": list(Rr)})
    agg = emb.groupby("pdb_id", sort=False).agg(
        lig=("lig", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
        rec=("rec", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
    ).reset_index()
    key = d[["pdb_id"]].merge(agg, on="pdb_id", how="left")
    LigA = np.stack(key["lig"].to_numpy()).astype(np.float32)
    RecA = np.stack(key["rec"].to_numpy()).astype(np.float32)

    # ---- 绑定 B: 逐行精确映射 ----
    d2 = d.copy()
    d2["_seq"] = d2.groupby("interface_key").cumcount()
    pos = d2.set_index(["interface_key", "_seq"]).index
    loc = pos.get_indexer(pd.MultiIndex.from_arrays([d2["interface_key"], d2["_seq"]]))
    assert (loc >= 0).all(), "逐行映射失败"
    LigB, RecB = Lr[loc], Rr[loc]

    emit("")
    emit("| 绑定 | 全局唯一嵌入行数 | test 内唯一嵌入行数 |")
    emit("|---|---:|---:|")
    for nm, LL in (("A. `pdb_id` 平均", LigA), ("B. `row_exact`", LigB)):
        emit(f"| {nm} | {np.unique(LL, axis=0).shape[0]} / {len(LL)} | "
             f"{np.unique(LL[te], axis=0).shape[0]} / {int(te.sum())} |")

    from sklearn.decomposition import PCA
    results: list[dict] = []

    def run(tag: str, Lig: np.ndarray, Rec: np.ndarray, train_mask: np.ndarray,
            note: str, feat: str = "physics+esm2_both") -> None:
        try:
            pca = PCA(n_components=args.pca_dim, random_state=SEED)
            pca.fit(np.hstack([Lig, Rec])[train_mask])
            Z = pca.transform(np.hstack([Lig, Rec]))
            X = (Xp if feat == "physics_only" else
                 (Z.astype(np.float32) if feat == "esm2_both"
                  else np.hstack([Xp, Z]).astype(np.float32)))
            ev = te & usable
            mdl = lgb.LGBMRegressor(**PARAMS)
            mdl.fit(X[train_mask & usable], y[train_mask & usable])
            pred = mdl.predict(X[ev])
            rho = mr.spearman_rho(-y[ev], -pred)
            lo, hi, _ = C.bootstrap_ci(mr.spearman_rho, -y[ev], -pred, 1000, SEED)
            results.append({"setting": tag, "note": note, "feature_set": feat,
                            "n_train": int((train_mask & usable).sum()),
                            "n_test": int(ev.sum()), "spearman": float(rho),
                            "ci_low": lo, "ci_high": hi,
                            "n_unique_emb_train": int(np.unique(
                                Lig[train_mask & usable], axis=0).shape[0])})
            emit(f"- **{tag}** ({note}): Spearman = {rho:+.4f} "
                 f"[{lo:+.4f}, {hi:+.4f}] (train {int((train_mask & usable).sum())} 行)")
        except Exception as exc:  # noqa: BLE001
            emit(f"- !! {tag} 失败: {exc}")
            with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[04_binding_control] {tag}: {exc}\n{traceback.format_exc()}\n")

    emit("")
    emit("## 1. 对照结果 (同一固定 V2.2 test 集)")
    emit("")
    run("A", LigA, RecA, tr, "pdb_id 平均 = V2.3.2 口径")
    run("B", LigB, RecB, tr, "row_exact 全训练集")
    # C: 把 test 的 pdb 在训练集中**所有**行删掉 (严格靶点不重叠)
    te_pdbs = set(d.loc[te, "pdb_id"].astype(str))
    tr_no_te_pdb = tr & ~d["pdb_id"].astype(str).isin(te_pdbs).to_numpy()
    emit("")
    emit(f"- 设置 C 的训练集规模: {int(tr_no_te_pdb.sum())} 行 "
         f"(从 {int(tr.sum())} 行删除 {int((tr & ~tr_no_te_pdb).sum())} 行 "
         f"= 与 test 同 pdb 的训练行)")
    run("C", LigB, RecB, tr_no_te_pdb, "row_exact + 严格靶点不重叠 (删同 pdb 训练行)")
    # D: 只用 physics (对照物理特征的上限)
    run("D", LigB, RecB, tr, "physics_only (物理特征, 无嵌入)", feat="physics_only")
    # E: 只用 esm2 (无物理)
    run("E", LigB, RecB, tr, "esm2_both 单独 (无物理特征)", feat="esm2_both")

    rdf = pd.DataFrame(results)
    rdf.to_csv(OUT_CSV, index=False, encoding="utf-8")

    emit("")
    emit("## 2. 结论")
    emit("")
    if len(rdf):
        def g(tag):
            r = rdf[rdf["setting"] == tag]
            return float(r.iloc[0]["spearman"]) if len(r) else float("nan")
        A, B, Cc, D, E = g("A"), g("B"), g("C"), g("D"), g("E")
        emit(f"- V2.3.2 记录值: {V232_RHO:+.4f} | 本脚本设置 A: {A:+.4f} "
             f"(差 {A - V232_RHO:+.4f})")
        emit(f"- **绑定方式的影响 (A -> B)**: {A:+.4f} -> {B:+.4f} "
             f"= **{B - A:+.4f}**")
        emit(f"- **泄漏的影响 (B -> C, 删掉同 pdb 训练行)**: {B:+.4f} -> {Cc:+.4f} "
             f"= **{Cc - B:+.4f}**")
        emit("")
        if B - A > 0.05 and Cc - B < -0.05:
            emit("> **判定: `row_exact` 的高分主要来自 target-identification 泄漏。**")
            emit("> 同一 PDB 的其它突变体在训练集时, 模型借其序列识别出测试复合物,")
            emit("> 从而回忆起该复合物的亲和力量级 —— 这不是 binding-affinity 泛化。")
            emit("> 因此**任务 1/2 中 `row_exact` 的数字必须标注为乐观上界**,")
            emit("> 跨来源比较只有在**同一绑定方式**下才有效。")
        elif abs(Cc - A) < 0.05:
            emit("> **判定: 严格靶点不重叠后, row_exact ≈ pdb_id 平均** ->")
            emit("> 绑定方式的表面增益确实来自泄漏; 两者在无泄漏设置下等价。")
        else:
            emit("> 结果较为微妙, 需结合逐设置 CI 解读 (见上表), 不做强断言。")
        emit("")
        emit(f"- 物理特征单独 (D): {D:+.4f}; ESM-2 单独 (E): {E:+.4f}")

    emit("")
    emit("## 3. 局限")
    emit("")
    emit("1. 设置 C 删除同 pdb 训练行会**同时减少训练量**, 因此 C 偏低也可能是"
         "训练数据变少所致 —— 与任务 2 同样的**混合效应**, 无法完全分离。")
    emit("2. `interface_key` 级不重叠 (A/B) 与 `pdb_id` 级不重叠 (C) 是两种不同严格度;")
    emit("   真实部署更接近 C (遇到全新复合物)。")
    emit("3. 单一划分 (V2.2 的固定 test, n=752), 未做多次划分。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")

    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD} ({OUT_MD.stat().st_size} B)")
    print(f"写出 {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
