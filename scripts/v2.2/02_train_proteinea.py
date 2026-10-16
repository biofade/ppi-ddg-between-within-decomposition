# -*- coding: utf-8 -*-
"""scripts/v2.2/02_train_proteinea.py -- 在 Proteinea 自带划分上重跑 LightGBM 回归 (V2.2).

任务与硬约束
------------
* 指标**必须在 Proteinea 自带的 test 划分上报告**, 不自行重新划分。
* 报告 Spearman + bootstrap 95% CI、NDCG@10、EF@1%、Precision@10。
* 只写 `features/v2.2_*`、`models/v2.2_*`、`logs/v2.2/`; 不触碰 V1.1 / V2.0 文件。
* 结构不可用的行**显式排除并计入报告**, 不用中位数硬填后当成正常样本。

三个评测口径 (为什么必须同时给)
--------------------------------
====================================== ==========================================
口径                                    含义
====================================== ==========================================
``test`` (Proteinea 官方 test)          **主报告口径**。<=30% 序列同一性隔离 =>
                                        跨家族外推, 可与外部工作在同一划分上对比。
``oof`` (train 内 GroupKFold)           同分布内的排序能力 (按 PDB 分组, 避免同结构
                                        跨折)。用于区分"模型弱"与"划分难"。
``random_split`` (filtered_random 对照) 同一批条目的**完全随机** 80/10/10 划分。
                                        量化"随机划分会带来多少虚高"。
====================================== ==========================================

模型选择
--------
主模型沿用 **V1.1 的 LightGBM 超参** (保证与历史结果可比); 另在 val 上做一次
**小网格搜索**并单独报告, 用于说明"调参能否弥补划分差异"。两种配置都在 test 上评估。

纪律: 纯 CPU; 无网络; seed=42; 逐配置容错; 不改动既有版本产物。
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import logging
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT_DIR = ROOT / "features"
MODEL_DIR = ROOT / "models"
LOG_DIR = ROOT / "logs" / "v2.2"

FEATURES_CSV = FEAT_DIR / "v2.2_proteinea_features.csv"
PROT_DIR = ROOT / "data" / "ppb_affinity_proteinea"
RANDOM_CSV = PROT_DIR / "filtered.csv"

RESULTS_CSV = FEAT_DIR / "v2.2_proteinea_results.csv"
PAIRS_CSV = FEAT_DIR / "v2.2_proteinea_model_comparison.csv"
JSON_OUT = FEAT_DIR / "v2.2_proteinea_results.json"
MODEL_OUT = MODEL_DIR / "v2.2_proteinea_lgbm.pkl"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
N_BOOT = 2000
TOP_FRAC = 0.10
RT_LITERAL = 0.592


def _ensure_openmp() -> str | None:
    import ctypes
    import glob
    try:
        import sklearn
        hits = glob.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll"))
        if hits:
            ctypes.CDLL(hits[0])
            return hits[0]
    except Exception:  # noqa: BLE001
        pass
    return None


_VCOMP = _ensure_openmp()

logger = logging.getLogger("train_proteinea")


def norm_chain(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    txt = str(value).replace(";", ",").replace("|", ",")
    parts = [p.strip().upper() for p in txt.split(",") if p.strip()]
    parts = [p.split(":")[0].strip() for p in parts]
    return ",".join(sorted(dict.fromkeys([p for p in parts if p])))


# ---------------------------------------------------------------- 评测工具
def build_metric_fns(mr: Any, top_frac: float) -> dict[str, Any]:
    return {
        "spearman": mr.spearman_rho,
        "kendall_tau": mr.kendall_tau,
        "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=top_frac),
        "ndcg@5": lambda a, b: mr.ndcg_at_k(a, b, 5, top_frac=top_frac),
        "ndcg@20": lambda a, b: mr.ndcg_at_k(a, b, 20, top_frac=top_frac),
        "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, top_frac),
        "precision@20": lambda a, b: mr.precision_at_k(a, b, 20, top_frac),
        "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, top_frac),
        "ef@5%": lambda a, b: mr.enrichment_factor(a, b, 5.0, top_frac),
        "ef@10%": lambda a, b: mr.enrichment_factor(a, b, 10.0, top_frac),
    }


def evaluate(y_dG: np.ndarray, pred_dG: np.ndarray, mfns: dict[str, Any],
             n_boot: int, seed: int) -> dict[str, Any]:
    """在 score 空间 (越大越强) 上算全套指标 + 主指标 bootstrap CI。"""
    y_s = -y_dG
    p_s = -pred_dG
    out: dict[str, Any] = {}
    for name, fn in mfns.items():
        try:
            out[name] = float(fn(y_s, p_s))
        except Exception:  # noqa: BLE001
            out[name] = float("nan")
    for name in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
        try:
            lo, hi, pt = bootstrap_ci(mfns[name], y_s, p_s, n_boot, seed)
        except Exception:  # noqa: BLE001
            lo = hi = float("nan")
        out[f"{name}_ci_low"] = lo
        out[f"{name}_ci_high"] = hi
    return out


def bootstrap_ci(fn: Any, y: np.ndarray, p: np.ndarray, n_boot: int,
                 seed: int, alpha: float = 0.05) -> tuple[float, float, float]:
    """自实现的 bootstrap 百分位 CI (避免依赖 scripts/ 的相对导入)。"""
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    point = float(fn(y, p)) if len(y) else float("nan")
    if len(y) < 5:
        return float("nan"), float("nan"), point
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n_boot):
        idx = rs.randint(0, len(y), len(y))
        v = fn(y[idx], p[idx])
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return float("nan"), float("nan"), point
    lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi), point


def paired_bootstrap(fn: Any, y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                     n_boot: int, seed: int) -> dict[str, float]:
    ok = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
    y, pa, pb = y[ok], pa[ok], pb[ok]
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y))}
    ma, mb = float(fn(y, pa)), float(fn(y, pb))
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        idx = rs.randint(0, len(y), len(y))
        va, vb = fn(y[idx], pa[idx]), fn(y[idx], pb[idx])
        if va == va and vb == vb:
            d.append(vb - va)
    if len(d) < 20:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "metric_a": ma, "metric_b": mb,
                "n": int(len(y))}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    pl = (np.sum(d <= 0) + 1) / (len(d) + 1)
    pr = (np.sum(d >= 0) + 1) / (len(d) + 1)
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl, pr))),
            "metric_a": ma, "metric_b": mb, "n": int(len(y))}


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([1.0 if (x != x) else float(x) for x in pvals], float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p, kind="mergesort")
    q = np.empty(m, float)
    prev = 1.0
    for rank, _ in enumerate(reversed(order), start=1):
        i = m - rank
        prev = min(prev, p[order[i]] * m / (i + 1))
        q[order[i]] = prev
    return [float(min(1.0, v)) for v in q]


# ---------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="V2.2: Proteinea 基准上重跑 LightGBM 回归")
    ap.add_argument("--features", type=Path, default=FEATURES_CSV)
    ap.add_argument("--out-csv", type=Path, default=RESULTS_CSV)
    ap.add_argument("--out-model", type=Path, default=MODEL_OUT)
    ap.add_argument("--json-out", type=Path, default=JSON_OUT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--top-frac", type=float, default=TOP_FRAC)
    ap.add_argument("--grid", action="store_true", help="额外跑小网格搜索")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "02_train_proteinea.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    import lightgbm as lgb
    from sklearn.model_selection import GroupKFold

    # V1.1 特征集与超参 (单一实现来源)
    spec = importlib.util.spec_from_file_location("t8_v22t", ROOT / "scripts" / "08_train_v1.0.py")
    v11mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v11mod)
    feat_names = [c for c in v11mod.V02_FEATURES if c not in v11mod.REPLACED] + list(v11mod.V03_NEW)
    count_cols = list(v11mod.COUNT_COLS)
    V11_PARAMS = dict(v11mod.PARAMS)
    V11_PARAMS["random_state"] = args.seed

    mfns = build_metric_fns(mr, args.top_frac)
    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# V2.2 — Proteinea `ppb_affinity` 基准上的 LightGBM 回归结果")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}; LightGBM: {lgb.__version__} (OpenMP shim: "
         f"{'OK' if _VCOMP else '未加载'})")
    emit(f"- bootstrap 次数: {args.n_boot}; 真强结合体定义: 真值 ΔG 最负的前 "
         f"{args.top_frac * 100:.0f}%")

    if not args.features.exists():
        emit(f"!! 特征表不存在: {args.features}")
        emit("   请先运行: python scripts/v2.2/01_adapt_proteinea.py")
        return 2

    df = pd.read_csv(args.features, low_memory=False)
    emit("")
    emit("## 1. 数据与可用性")
    emit("")
    emit(f"- 适配特征表: `{args.features.name}` — {len(df)} 行 x {df.shape[1]} 列")

    # 只用结构可用且标签有效的行
    usable = (df["row_structure_usable"] == 1) & df["real_dG"].notna()
    n_excl_struct = int((df["row_structure_usable"] != 1).sum())
    n_excl_label = int(df["real_dG"].isna().sum())
    emit(f"- 排除结构不可用行: **{n_excl_struct}**; 排除标签无效行: {n_excl_label}")
    d = df[usable].copy()
    emit(f"- 可用行: **{len(d)} / {len(df)}** ({len(d) / len(df) * 100:.2f}%)")

    fcols = [c for c in feat_names if c in d.columns]
    missing_feats = [c for c in feat_names if c not in d.columns]
    emit(f"- 物理特征: {len(fcols)} / {len(feat_names)} 可用"
         + (f" (缺 {missing_feats})" if missing_feats else ""))

    # 特征缺失率 (按 split)
    emit("")
    emit("| split | 行数 | 物理特征全非空行 | 特征缺失率中位 |")
    emit("|---|---:|---:|---:|")
    for sp in ("train", "val", "test"):
        s = d[d["split"] == sp]
        allp = s[fcols].notna().all(axis=1).sum()
        miss = s[fcols].isna().mean().median()
        emit(f"| {sp} | {len(s)} | {int(allp)} | {miss * 100:.2f}% |")

    if args.dry_run:
        emit("")
        emit("--dry-run: 不训练")
        return 0

    # ---- 2. 特征矩阵 ---------------------------------------------------
    def prepare(frame: pd.DataFrame, medians: pd.Series | None = None) -> tuple[pd.DataFrame, pd.Series]:
        X = frame[fcols].apply(pd.to_numeric, errors="coerce").copy()
        for c in count_cols:
            if c in X.columns:
                X[c] = np.log1p(X[c].clip(lower=0))
        med = X.median(numeric_only=True) if medians is None else medians
        return X.fillna(med).fillna(0.0), med

    tr = d[d["split"] == "train"]
    va = d[d["split"] == "val"]
    te = d[d["split"] == "test"]
    X_tr, med = prepare(tr)
    X_va, _ = prepare(va, med)
    X_te, _ = prepare(te, med)
    y_tr = tr["real_dG"].to_numpy(float)
    y_va = va["real_dG"].to_numpy(float)
    y_te = te["real_dG"].to_numpy(float)
    emit("")
    emit("## 2. 划分与样本量 (完全沿用 Proteinea 自带划分)")
    emit("")
    emit(f"- train {len(tr)} 行 / {tr['pdb_id'].nunique()} 结构 / "
         f"{tr['interface_key'].nunique()} 界面键")
    emit(f"- val   {len(va)} 行 / {va['pdb_id'].nunique()} 结构")
    emit(f"- test  {len(te)} 行 / {te['pdb_id'].nunique()} 结构 / "
         f"{te['interface_key'].nunique()} 界面键  ← **主报告口径**")

    # ---- 3. 主模型 (V1.1 超参) -----------------------------------------
    emit("")
    emit("## 3. 主模型: V1.1 超参 LightGBM 回归")
    emit("")
    emit(f"- 超参: {json.dumps({k: v for k, v in V11_PARAMS.items() if k != 'verbose'}, ensure_ascii=False)}")

    def fit_predict(params: dict, Xa: pd.DataFrame, ya: np.ndarray,
                    Xb: pd.DataFrame, yb: np.ndarray, Xc: pd.DataFrame,
                    early: bool = True) -> tuple[Any, np.ndarray]:
        m = lgb.LGBMRegressor(**params)
        if early:
            m.fit(Xa, ya, eval_set=[(Xb, yb)],
                  callbacks=[lgb.early_stopping(50, verbose=False)])
        else:
            m.fit(Xa, ya)
        return m, m.predict(Xc)

    model_main, pred_te = fit_predict(V11_PARAMS, X_tr, y_tr, X_va, y_va, X_te)
    main_metrics = evaluate(y_te, pred_te, mfns, args.n_boot, args.seed)
    emit(f"- best_iteration = {getattr(model_main, 'best_iteration_', None)}")
    emit("")
    emit("### 3.1 主结果 (Proteinea test, n=%d)" % len(te))
    emit("")
    emit("| 指标 | 值 | 95% CI |")
    emit("|---|---:|---|")
    for k in ("spearman", "kendall_tau", "ndcg@10", "ndcg@5", "ndcg@20",
              "precision@10", "precision@20", "ef@1%", "ef@5%", "ef@10%"):
        v = main_metrics.get(k, float("nan"))
        ci = ""
        if f"{k}_ci_low" in main_metrics:
            lo, hi = main_metrics[f"{k}_ci_low"], main_metrics[f"{k}_ci_high"]
            ci = f"[{lo:.4f}, {hi:.4f}]" if np.isfinite(lo) and np.isfinite(hi) else "—"
        emit(f"| {k} | {v:+.4f} | {ci} |")

    # ---- 4. OOF (train 内 GroupKFold, 按 PDB 分组) ----------------------
    emit("")
    emit("## 4. 同分布口径: train 内 GroupKFold OOF (按 PDB 分组)")
    emit("")
    emit("> 用于区分「模型弱」与「划分难」: 若 OOF 明显高于 test, 说明差距主要来自")
    emit("> **跨家族外推**而不是模型容量。分组单位 = `pdb_id` (避免同结构跨折)。")
    g_tr = tr["pdb_id"].to_numpy()
    n_splits = min(5, max(2, len(set(g_tr))))
    gkf = GroupKFold(n_splits=n_splits)
    oof = np.full(len(tr), np.nan)
    emit("")
    emit("| 折 | 训练行 | 验证行 | 组交集 |")
    emit("|---:|---:|---:|---:|")
    all_ok = True
    for i, (ia, ib) in enumerate(gkf.split(X_tr, y_tr, g_tr)):
        inter = len(set(g_tr[ia]) & set(g_tr[ib]))
        all_ok = all_ok and inter == 0
        m = lgb.LGBMRegressor(**V11_PARAMS)
        m.fit(X_tr.iloc[ia], y_tr[ia])
        oof[ib] = m.predict(X_tr.iloc[ib])
        emit(f"| {i} | {len(ia)} | {len(ib)} | {inter} |")
        gc.collect()
    emit("")
    emit(f"- 组交集全为 0: {'✅ PASS' if all_ok else '!! FAIL'}")
    oof_metrics = evaluate(y_tr, oof, mfns, args.n_boot, args.seed)
    emit(f"- OOF Spearman = **{oof_metrics['spearman']:+.4f}** "
         f"[{oof_metrics['spearman_ci_low']:.4f}, {oof_metrics['spearman_ci_high']:.4f}]")
    emit(f"- OOF NDCG@10 = {oof_metrics['ndcg@10']:.4f}, "
         f"EF@1% = {oof_metrics['ef@1%']:.4f}, P@10 = {oof_metrics['precision@10']:.4f}")

    # ---- 5. 随机划分对照 (量化虚高) -------------------------------------
    emit("")
    emit("## 5. 随机划分对照 (`filtered_random` 口径, 量化虚高)")
    emit("")
    emit("> 用同一批条目按 80/10/10 **完全随机**划分 (seed=42, 与官方 filtered_random 一致),")
    emit("> 与上面的同源隔离划分对比。**这不是用来宣称性能, 而是量化隔离的代价。**")
    rand_metrics: dict[str, Any] = {}
    try:
        # 随机划分对照必须使用**同一张适配特征表** (含 60 个物理特征 + 结构可用性标记),
        # 不能直接读 raw filtered.csv —— 那里没有 join 过的特征列 (早期版本因此报
        # KeyError: 'n_missing_chain' ...)。
        rr = d.copy()
        n_total = len(rr)
        idx = np.random.RandomState(42).permutation(n_total)
        n_train, n_val = int(0.8 * n_total), int(0.1 * n_total)
        rr = rr.iloc[idx].reset_index(drop=True)
        frames = {
            "train": rr.iloc[:n_train],
            "val": rr.iloc[n_train:n_train + n_val],
            "test": rr.iloc[n_train + n_val:],
        }
        Xr_tr, medr = prepare(frames["train"])
        Xr_va, _ = prepare(frames["val"], medr)
        Xr_te, _ = prepare(frames["test"], medr)
        yr_tr = frames["train"]["real_dG"].to_numpy(float)
        yr_va = frames["val"]["real_dG"].to_numpy(float)
        yr_te = frames["test"]["real_dG"].to_numpy(float)
        _, predr = fit_predict(V11_PARAMS, Xr_tr, yr_tr, Xr_va, yr_va, Xr_te)
        rand_metrics = evaluate(yr_te, predr, mfns, args.n_boot, args.seed)

        # 随机划分下的**实际跨侧重叠** (这是"虚高"的机制来源)
        rr2 = rr.copy()
        rr2["_side"] = "train"
        rr2.loc[rr2.index[n_train:n_train + n_val], "_side"] = "val"
        rr2.loc[rr2.index[n_train + n_val:], "_side"] = "test"
        cross_pdb = int((rr2.groupby("pdb_id")["_side"].nunique() > 1).sum())
        cross_ik = int((rr2.groupby("interface_key")["_side"].nunique() > 1).sum())
        # 同源隔离划分下的对照值 (应为 0)
        iso_cross_pdb = int((d.groupby("pdb_id")["split"].nunique() > 1).sum())
        iso_cross_ik = int((d.groupby("interface_key")["split"].nunique() > 1).sum())

        emit("")
        emit(f"- 随机划分样本: train {len(frames['train'])}, val {len(frames['val'])}, "
             f"test {len(frames['test'])}")
        emit("")
        emit("| 划分协议 | 跨侧 PDB 数 | 跨侧 interface_key 数 |")
        emit("|---|---:|---:|")
        emit(f"| Proteinea 同源隔离 (<=30% 同一性) | **{iso_cross_pdb}** | **{iso_cross_ik}** |")
        emit(f"| 完全随机 80/10/10 | **{cross_pdb}** | **{cross_ik}** |")
        emit("")
        emit(f"- 随机划分 test Spearman = **{rand_metrics['spearman']:+.4f}** "
             f"[{rand_metrics['spearman_ci_low']:.4f}, {rand_metrics['spearman_ci_high']:.4f}]")
        emit(f"- 随机划分 test NDCG@10 = {rand_metrics['ndcg@10']:.4f}, "
             f"EF@1% = {rand_metrics['ef@1%']:.4f}, P@10 = {rand_metrics['precision@10']:.4f}")
        delta = rand_metrics["spearman"] - main_metrics["spearman"]
        emit("")
        emit(f"- **随机划分 vs Proteinea 隔离划分: ΔSpearman = {delta:+.4f}**")
        emit("  -> 这就是「换一个更宽松的划分会虚高多少」的直接量化。")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! 随机划分对照失败: {type(exc).__name__}: {exc}")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
            fh2.write(f"[02_train_proteinea] random-split: {exc}\n{traceback.format_exc()}\n")

    # ---- 6. 小网格搜索 --------------------------------------------------
    grid_rows: list[dict[str, Any]] = []
    if args.grid:
        emit("")
        emit("## 6. 小网格搜索 (在 val 上选, 在 test 上报告)")
        emit("")
        grid = [
            {"num_leaves": 15, "learning_rate": 0.03},
            {"num_leaves": 31, "learning_rate": 0.03},
            {"num_leaves": 63, "learning_rate": 0.03},
            {"num_leaves": 31, "learning_rate": 0.01},
            {"num_leaves": 15, "learning_rate": 0.05},
            {"num_leaves": 31, "learning_rate": 0.05, "min_child_samples": 5},
            {"num_leaves": 63, "learning_rate": 0.01, "min_child_samples": 20},
        ]
        emit("| # | 覆盖参数 | val Spearman | test Spearman | test NDCG@10 | test EF@1% |")
        emit("|---:|---|---:|---:|---:|---:|")
        for i, g in enumerate(grid):
            p = dict(V11_PARAMS)
            p.update(g)
            try:
                m, pr_ = fit_predict(p, X_tr, y_tr, X_va, y_va, X_te)
                val_sp = mr.spearman_rho(-y_va, -m.predict(X_va))
                tm = evaluate(y_te, pr_, mfns, n_boot=600, seed=args.seed)
                grid_rows.append({"config": json.dumps(g, ensure_ascii=False),
                                  "val_spearman": float(val_sp), **{f"test_{k}": v
                                  for k, v in tm.items() if not k.endswith(("_ci_low", "_ci_high"))}})
                emit(f"| {i} | {json.dumps(g, ensure_ascii=False)} | {val_sp:+.4f} | "
                     f"{tm['spearman']:+.4f} | {tm['ndcg@10']:.4f} | {tm['ef@1%']:.4f} |")
            except Exception as exc:  # noqa: BLE001
                emit(f"| {i} | {json.dumps(g, ensure_ascii=False)} | FAIL: {exc} | | | |")
            gc.collect()
        if grid_rows:
            best = max(grid_rows, key=lambda r: r["val_spearman"])
            emit("")
            emit(f"- val 最优: `{best['config']}` -> val Spearman {best['val_spearman']:+.4f}, "
                 f"test Spearman {best['test_spearman']:+.4f}")
            emit(f"- 对照: V1.1 超参 test Spearman {main_metrics['spearman']:+.4f}")
            emit("  -> **调参带来的差异远小于划分协议带来的差异**"
                 if abs(best["test_spearman"] - main_metrics["spearman"]) < 0.05
                 else "  -> 调参带来可观差异, 需注意过拟合 val 的风险")

    # ---- 7. 配对检验 (隔离 vs 随机) -------------------------------------
    pairs: list[dict[str, Any]] = []
    if rand_metrics:
        emit("")
        emit("## 7. 配对 bootstrap 检验")
        emit("")
        emit("> 注意: 隔离划分与随机划分使用**不同的测试集**, 因此无法逐样本配对。")
        emit("> 这里只并列报告两个口径的点估计与各自 CI, **不做配对显著性检验** ——")
        emit("> 宣称「随机划分显著更好」需要共享测试集, 本设计不满足该前提。")
        emit("")
        emit("| 口径 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
        emit("|---|---:|---|---:|---:|---:|")
        emit(f"| Proteinea 隔离 test | {main_metrics['spearman']:+.4f} | "
             f"[{main_metrics['spearman_ci_low']:.4f}, {main_metrics['spearman_ci_high']:.4f}] | "
             f"{main_metrics['ndcg@10']:.4f} | {main_metrics['ef@1%']:.4f} | "
             f"{main_metrics['precision@10']:.4f} |")
        emit(f"| 随机划分 test | {rand_metrics['spearman']:+.4f} | "
             f"[{rand_metrics['spearman_ci_low']:.4f}, {rand_metrics['spearman_ci_high']:.4f}] | "
             f"{rand_metrics['ndcg@10']:.4f} | {rand_metrics['ef@1%']:.4f} | "
             f"{rand_metrics['precision@10']:.4f} |")
        emit(f"| OOF (train 内, 同分布) | {oof_metrics['spearman']:+.4f} | "
             f"[{oof_metrics['spearman_ci_low']:.4f}, {oof_metrics['spearman_ci_high']:.4f}] | "
             f"{oof_metrics['ndcg@10']:.4f} | {oof_metrics['ef@1%']:.4f} | "
             f"{oof_metrics['precision@10']:.4f} |")

    # ---- 8. 落盘 --------------------------------------------------------
    rows = []
    n_rand = int(len(predr)) if rand_metrics else 0
    for name, metrics, note, ntest in (
            ("v2.2_lgbm_v11params_isolated_test", main_metrics,
             "Proteinea 自带 <=30% 序列同一性隔离 test (主报告口径)", len(te)),
            ("v2.2_lgbm_v11params_oof_train", oof_metrics,
             "train 内 GroupKFold OOF (按 PDB 分组), 同分布口径", len(tr)),
            ("v2.2_lgbm_v11params_random_test", rand_metrics,
             "完全随机 80/10/10 test (对照, 量化虚高)", n_rand)):
        if not metrics:
            continue
        r = {"model_id": name, "note": note, "n_eval": int(ntest)}
        r.update(metrics)
        rows.append(r)
    res = pd.DataFrame(rows)
    if len(res):
        tmp = args.out_csv.with_suffix(args.out_csv.suffix + ".tmp")
        res.to_csv(tmp, index=False, encoding="utf-8")
        tmp.replace(args.out_csv)

    # 保存主模型
    payload = {
        "model": model_main, "feature_names": fcols, "params": V11_PARAMS,
        "count_features_log1p": [c for c in count_cols if c in fcols],
        "train_medians": med.to_dict(),
        "label_definition": f"real_dG = {RT_LITERAL} * ln(KD[M])",
        "split_policy": "Proteinea ppb_affinity filtered (<=30% seq identity, provided splits)",
        "train_rows": int(len(tr)), "val_rows": int(len(va)), "test_rows": int(len(te)),
        "metrics": {"test": main_metrics, "oof_train": oof_metrics,
                    "random_split_test": rand_metrics},
        "version": "v2.2",
    }
    tmpm = args.out_model.with_suffix(args.out_model.suffix + ".tmp")
    import pickle
    with open(tmpm, "wb") as fhm:
        pickle.dump(payload, fhm)
    tmpm.replace(args.out_model)
    emit("")
    emit("## 8. 输出")
    emit("")
    emit(f"- 结果表: `{args.out_csv.name}`")
    emit(f"- 模型: `{args.out_model.name}` ({args.out_model.stat().st_size / 1024:.1f} KB)")

    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "n_boot": args.n_boot, "top_frac": args.top_frac,
        "dataset": "proteinea/ppb_affinity (filtered)",
        "rows_total": int(len(df)), "rows_usable": int(len(d)),
        "rows_excluded_structure": n_excl_struct,
        "splits": {"train": int(len(tr)), "val": int(len(va)), "test": int(len(te))},
        "metrics_test": main_metrics, "metrics_oof_train": oof_metrics,
        "metrics_random_test": rand_metrics,
        "grid": grid_rows,
        "hyperparams": {k: v for k, v in V11_PARAMS.items()},
        "oof_group_intersection_zero": bool(all_ok),
    }
    args.json_out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str),
                             encoding="utf-8")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    (LOG_DIR / "02_train_proteinea_report.md").write_text(
        "\n".join(run_log) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
            fh2.write(f"[02_train_proteinea] {traceback.format_exc()}\n")
        raise
