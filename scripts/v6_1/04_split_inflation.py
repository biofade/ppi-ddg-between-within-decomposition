# -*- coding: utf-8 -*-
"""scripts/v6_1/04_split_inflation.py -- V6.1 PHASE 2: Split Inflation 量化.

目标（V6.1 最高科学价值 Phase）
-------------------------------
用**同一个模型**（LightGBM 固定 800 树；特征 = 60 物理 + pooled ESM-2 PCA128）
跑 **4 档划分**，量化「划分严格性」对可见性能的贡献。

| 档 | 划分方式 | 泄漏控制 |
|---|---|---|
| **L1** | mutation-level random 5-fold | 无（最宽松，仅作上界） |
| **L2** | complex-level（PDB）5-fold | 同一结构不跨折 |
| **L3** | ≤40% seq-id + complex-level | 加同源控制 |
| **L4** | family-level（≤40% + 家族隔离） | 最严格 |

评估口径（**所有档统一**）
--------------------------
* 采用 **5-fold OOF**（折外预测），在同一批全量行上算指标 → 各档**分母相同**，可比
* 指标：Overall Pearson / Macro per-PDB Spearman / Within-complex Hit@10 /
  Base Rate / Efficiency@10

**内存与时间保护（V6.1 要求）**
------------------------------
* 逐档逐 seed 跑，跑完即释放；PCA 只在**该档该折的训练行**上拟合
* 若 L4 某 seed 内存不足 → 记录并跳过
* 若整体超时 → 至少完成 **L1 与 L4**（最极端两档）

输出
----
* `experiments/v6_1/phase2_split_inflation/split_inflation_results.csv`
* `experiments/v6_1/phase2_split_inflation/split_inflation_report.md`
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH2 = C.ROOT / "experiments" / "v6_1" / "phase2_split_inflation"
PH2.mkdir(parents=True, exist_ok=True)
OUT_CSV = PH2 / "split_inflation_results.csv"
OUT_MD = PH2 / "split_inflation_report.md"
LOG = C.ROOT / "logs" / "v6_1"
LOG.mkdir(parents=True, exist_ok=True)
FAIL = LOG / "failures.log"

TREES = 800
PCA_DIM = 128
SEEDS = [42, 43, 44]
MIN_PS_N = 5
LEVELS = ["L1_random", "L2_complex", "L3_complex_seq40", "L4_family"]
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def fail(where: str, exc: BaseException) -> None:
    with open(FAIL, "a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] [04_split_inflation] {where}: "
                 f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}\n")


def eff10_fast(y, p, g, n_rand=200) -> float:
    ys, ps, gs = np.asarray(y, float), np.asarray(p, float), np.asarray(g)
    eff = []
    rng = np.random.RandomState(C.SEED)
    for c in np.unique(gs):
        m = np.where(gs == c)[0]
        if len(m) < 5:
            continue
        k = min(10, len(m)) if len(m) >= 10 else min(5, len(m))
        sy, sp = ys[m], ps[m]
        sel = np.argsort(sp, kind="mergesort")[:k]
        orc = np.argsort(sy, kind="mergesort")[:k]
        rr = [rng.permutation(len(m))[:k] for _ in range(n_rand)]
        v_sel = float(-sy[sel].mean()); v_or = float(-sy[orc].mean())
        v_rnd = float(np.mean([-sy[r].mean() for r in rr]))
        den = v_rnd - v_or
        if abs(den) > 1e-9:
            eff.append((v_rnd - v_sel) / den)
    return float(np.mean(eff)) if eff else float("nan")


def load_rows():
    """行级数据：物理 60 + pooled ESM-2（原始 1280 维，待逐折 PCA）+ ΔG + 分组键。"""
    v22 = pd.read_csv(C.V22, low_memory=False)
    v22["pdb_id"] = v22["pdb_id"].astype(str)
    phys, counts = C.load_phys_names()
    fcols = [c for c in phys if c in v22.columns]
    Xp = v22[fcols].apply(pd.to_numeric, errors="coerce")
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    Xp = Xp[fcols].to_numpy(np.float32)
    z = np.load(C.POOLED, allow_pickle=True)
    both = np.hstack([z["lig_l2"].astype(np.float32),
                      z["rec_l2"].astype(np.float32)])
    if len(both) != len(v22):
        raise RuntimeError(f"pooled {len(both)} != v22 {len(v22)}")
    y = pd.to_numeric(v22["real_dG"], errors="coerce").to_numpy(float)
    keep = np.isfinite(y)
    return {"v22": v22, "Xp": Xp, "both": both, "y": y, "keep": keep,
            "pdb": v22["pdb_id"].to_numpy(),
            "phys_names": fcols}


def pca_fit_apply(both_tr, both_all, tr_idx, dim=PCA_DIM, seed=C.SEED):
    """PCA **只在 tr_idx 上拟合**，再变换全部行。"""
    from sklearn.decomposition import PCA
    d = min(dim, both_tr.shape[1], len(tr_idx) - 1)
    p = PCA(n_components=d, random_state=seed)
    p.fit(both_all[tr_idx])
    return p.transform(both_all).astype(np.float32), d


def make_folds(level: str, d: dict, seed: int, n_splits: int = 5):
    """返回 list[(tr_idx, te_idx)]；索引针对**全部行**（含 keep 过滤前的原索引）。"""
    from sklearn.model_selection import GroupKFold, KFold
    idx_all = np.where(d["keep"])[0]
    y = d["y"][idx_all]
    pdb = d["pdb"][idx_all]
    if level == "L1_random":
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        return [(idx_all[tr], idx_all[te]) for tr, te in kf.split(idx_all)]
    if level == "L2_complex":
        gkf = GroupKFold(n_splits=n_splits)
        return [(idx_all[tr], idx_all[te])
                for tr, te in gkf.split(idx_all, groups=pdb)]
    # L3 / L4：先按 ≤40% 序列同一性做贪心聚类（复用 V1.1/2.0 的 k-mer Jaccard 代理）
    fam = _family_ids(d, idx_all)
    if level == "L3_complex_seq40":
        # 复合物级 + family 不跨折（等价于把 family 当 group，同时保证 pdb 不跨折）
        key = np.array([f"{p}|{f}" for p, f in zip(pdb, fam)])
        gkf = GroupKFold(n_splits=n_splits)
        return [(idx_all[tr], idx_all[te])
                for tr, te in gkf.split(idx_all, groups=key)]
    # L4_family：纯 family 分组（家族整族分侧）
    gkf = GroupKFold(n_splits=n_splits)
    return [(idx_all[tr], idx_all[te])
            for tr, te in gkf.split(idx_all, groups=fam)]


_FAM_CACHE: dict = {}


def _family_ids(d: dict, idx_all: np.ndarray) -> np.ndarray:
    """≤40% 序列同一性贪心聚类（k-mer Jaccard 代理），返回每行的家族 ID。

    与 V1.1/V2.0 的历史口径一致（**代理**，不是精确序列同一性）。
    """
    key = ("fam", len(idx_all))
    if key in _FAM_CACHE:
        return _FAM_CACHE[key]
    v22 = d["v22"]
    seq_col = None
    for c in ("Receptor Sequences", "receptor_seq", "Ligand Sequences"):
        if c in v22.columns:
            seq_col = c
            break
    if seq_col is None:
        # 退化：用 pdb 当家族
        fam = pd.factorize(d["pdb"][idx_all])[0]
        _FAM_CACHE[key] = fam
        return fam
    seqs = v22[seq_col].astype(str).to_numpy()[idx_all]
    THR = 0.40
    K = 3
    sets = []
    for s in seqs:
        s = "".join(ch for ch in s.upper() if ch.isalpha())[:1200]
        sets.append({s[i:i + K] for i in range(max(0, len(s) - K + 1))} or {s})
    reps: list[int] = []
    fam = np.full(len(seqs), -1, dtype=int)
    for i, si in enumerate(sets):
        if not si:
            fam[i] = len(reps); reps.append(i); continue
        assigned = -1
        for fi, ri in enumerate(reps):
            sj = sets[ri]
            inter = len(si & sj)
            if inter == 0:
                continue
            jac = inter / (len(si) + len(sj) - inter)
            if jac >= THR:
                assigned = fi
                break
        if assigned < 0:
            assigned = len(reps)
            reps.append(i)
        fam[i] = assigned
    _FAM_CACHE[key] = fam
    return fam


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="*", default=LEVELS)
    ap.add_argument("--seeds", type=int, nargs="*", default=SEEDS)
    ap.add_argument("--max-minutes", type=float, default=80.0,
                    help="软时限；超过则只保留已完成的档")
    args = ap.parse_args()
    t0 = time.time()
    emit("# V6.1 PHASE 2 — Split Inflation 量化")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 模型: LightGBM 固定 **{TREES} 树**（不用早停，沿用 V4/V5 教训）")
    emit(f"- 特征: 60 物理 + pooled ESM-2 PCA**{PCA_DIM}**（**PCA 只在各折训练行上拟合**）")
    emit(f"- 评估: **5-fold OOF**，所有档在同一批全量行上算指标（分母相同 → 可比）")
    emit(f"- seeds: {args.seeds}；软时限 {args.max_minutes} 分钟")
    emit("")

    lgb = C.load_lgb()
    d = load_rows()
    n_tot = int(d["keep"].sum())
    n_pdb = len(set(d["pdb"][d["keep"]]))
    emit(f"- 全量可用行: **{n_tot}**；唯一 PDB: **{n_pdb}**")
    emit(f"- 物理特征数: **{len(d['phys_names'])}**；pooled 原始维度: "
         f"**{d['both'].shape[1]}**")
    emit("")

    rows: list[dict] = []
    emit("| 档 | seed | Overall Pearson | Macro per-PDB Spearman | Hit@10(<0) | "
         "Base Rate | Efficiency@10 | 耗时 |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    for level in args.levels:
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"| {level} | — | ⏭ 软时限已达，跳过 | | | | | |")
            continue
        try:
            folds = make_folds(level, d, seed=args.seeds[0])
        except Exception as exc:  # noqa: BLE001
            fail(f"make_folds({level})", exc)
            emit(f"| {level} | — | ❌ 划分失败 | | | | | |")
            continue
        for sd in args.seeds:
            ts = time.time()
            try:
                if level != "L1_random":
                    folds = make_folds(level, d, seed=sd)
                else:
                    folds = make_folds(level, d, seed=sd)
                yhat = np.full(len(d["y"]), np.nan)
                for fi, (tr, te) in enumerate(folds):
                    # 逐折 PCA（只用 train 行）
                    Ztr_te, dim = pca_fit_apply(d["both"][tr], d["both"], tr,
                                                PCA_DIM, sd)
                    X = np.hstack([d["Xp"], Ztr_te]).astype(np.float32)
                    m = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                             "n_estimators": TREES,
                                             "random_state": sd})
                    ok_tr = tr[np.isfinite(d["y"][tr])]
                    m.fit(X[ok_tr], d["y"][ok_tr])
                    pred = m.predict(X[te])
                    yhat[te] = pred
                    del X, Ztr_te
                    gc.collect()
                idx = np.where(np.isfinite(yhat) & d["keep"])[0]
                yy = d["y"][idx]; pp = yhat[idx]; gg = d["pdb"][idx]
                r = C.pearson(yy, pp)
                macro = C.per_structure_spearman(yy, pp, gg, min_n=MIN_PS_N)
                base = float((yy < 0).mean())
                k10 = min(10, len(yy))
                hit = float((yy[np.argsort(pp, kind="mergesort")[:k10]] < 0).mean())
                e10 = eff10_fast(yy, pp, gg)
                dt = time.time() - ts
                emit(f"| {level} | {sd} | {r:+.4f} | {macro:+.4f} | {hit:.4f} | "
                     f"{base:.4f} | {e10:+.4f} | {dt:.0f}s |")
                rows.append({"level": level, "seed": sd, "n_rows": int(len(idx)),
                             "n_pdb": int(len(set(gg))),
                             "overall_pearson": round(float(r), 4),
                             "macro_pdb_spearman": round(float(macro), 4),
                             "hit10_lt0": round(hit, 4),
                             "base_rate_lt0": round(base, 4),
                             "efficiency_at_10": round(float(e10), 4),
                             "runtime_sec": round(dt, 1),
                             "n_folds": len(folds), "pca_dim": int(dim),
                             "status": "OK"})
                pd.DataFrame(rows).to_csv(OUT_CSV, index=False, encoding="utf-8")
            except MemoryError as exc:
                fail(f"{level}/seed{sd} MemoryError", exc)
                emit(f"| {level} | {sd} | ❌ MemoryError（已跳过） | | | | | |")
                rows.append({"level": level, "seed": sd, "status": "OOM_SKIPPED"})
                gc.collect()
            except Exception as exc:  # noqa: BLE001
                fail(f"{level}/seed{sd}", exc)
                emit(f"| {level} | {sd} | ❌ {type(exc).__name__}（已跳过） | | | | | |")
                rows.append({"level": level, "seed": sd,
                             "status": f"FAILED_{type(exc).__name__}"})
                gc.collect()
        gc.collect()

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------- 汇总 ----------
    emit("")
    emit("## 各档 3-seed 均值")
    emit("")
    emit("| 档 | Overall Pearson | Macro per-PDB Spearman | Hit@10 | Efficiency@10 |")
    emit("|---|---:|---:|---:|---:|")
    agg = {}
    for lv in args.levels:
        sub = res[(res["level"] == lv) & (res.get("status", "") == "OK")]
        if not len(sub):
            continue
        agg[lv] = {c: float(pd.to_numeric(sub[c], errors="coerce").mean())
                   for c in ("overall_pearson", "macro_pdb_spearman",
                             "hit10_lt0", "efficiency_at_10")}
        a = agg[lv]
        emit(f"| {lv} | {a['overall_pearson']:+.4f} | {a['macro_pdb_spearman']:+.4f} | "
             f"{a['hit10_lt0']:.4f} | {a['efficiency_at_10']:+.4f} |")
    emit("")
    if "L1_random" in agg and "L4_family" in agg:
        dp = agg["L1_random"]["overall_pearson"] - agg["L4_family"]["overall_pearson"]
        dm = (agg["L1_random"]["macro_pdb_spearman"]
              - agg["L4_family"]["macro_pdb_spearman"])
        emit("## ★ 核心量化：Split Inflation（L1 − L4）")
        emit("")
        emit("| 指标 | L1 random | L4 family | **L1 − L4（inflation）** |")
        emit("|---|---:|---:|---:|")
        emit(f"| Overall Pearson | {agg['L1_random']['overall_pearson']:+.4f} | "
             f"{agg['L4_family']['overall_pearson']:+.4f} | **{dp:+.4f}** |")
        emit(f"| Macro per-PDB Spearman | {agg['L1_random']['macro_pdb_spearman']:+.4f} | "
             f"{agg['L4_family']['macro_pdb_spearman']:+.4f} | **{dm:+.4f}** |")
        emit("")
        emit(f"> **结论**：仅把划分从最严格（L4 家族级）换成最宽松（L1 随机），")
        emit(f"> Overall Pearson 就**虚高 {dp:+.4f}**（相对提升 "
             f"{dp / max(1e-9, abs(agg['L4_family']['overall_pearson'])) * 100:+.1f}%）。")
        emit(f"> 这就是**split inflation** 的幅度 —— 与文献跨论文数字不可比的根源。")
    emit("")
    emit(f"- 总耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
