# -*- coding: utf-8 -*-
"""scripts/v2.5/v25_common.py -- V2.5 共享机器 (数据装配 / 指标 / 检验 / 物理特征).

设计要点
--------
* **评估口径统一 cross-target**: 使用 per-residue 缓存自带的 `splits`
  (与 V2.2 Proteinea 固定划分逐行一致, 已用 split 一致率 100% 与 md5 双重确认)。
* **pooled 基线必须在同一批行上重跑**: per-residue 缓存的 train 只有 5163 行,
  而 V2.2 可用行是 6482 —— 训练集不同则不可比。因此本模块提供
  `load_aligned()` 返回**同一批 6885 行**上的 pooled 与 per-residue 特征。
* 物理特征列与超参**沿用 V2.2/V2.3.2**(从 `08_train_v1.0.py` 载入, 单一实现来源)。
* CPU only; 大数组用 fp16 存储、按需转 fp32。
"""

from __future__ import annotations

import ctypes
import glob as _glob
import importlib.util
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
MODEL_DIR = ROOT / "models"
LOG_DIR = ROOT / "logs" / "v2.5"
LOG_DIR.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "failures.log"

PERRES_CACHE = FEAT / "v2.5_perres_cache.npz"
PERRES_RAW = FEAT / "v2_esm2_perres_v2.npz"
POOLED = FEAT / "v2_esm2_cache.npz"
V22 = FEAT / "v2.2_proteinea_features.csv"

SEED = 42
TOP_FRAC = 0.10
N_BOOT_PAIR = 5000     # 配对检验次数 (任务书要求)
N_BOOT_CI = 800        # 每个模型自身 CI 次数 (算力预算)

# V2.2 / V2.3.2 逐字超参
LGB_PARAMS = dict(objective="regression", n_estimators=500, learning_rate=0.03,
                  num_leaves=15, max_depth=5, min_child_samples=12, subsample=0.8,
                  subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                  n_jobs=2, random_state=SEED, verbose=-1)

# V2.3.2 关键参照值 (cross-target, Proteinea 固定 test n=752)
V232 = {
    "physics_only": 0.23747781344414348,
    "physics+esm2_both@l2@128": 0.3889577363940227,
    "physics+esm2_receptor@l2@128": 0.3423288621708283,
    "ensemble_rankmean": 0.4370,
}

# V2.4 within-target 参照 (非主指标)
V24_WITHIN = 0.1043


def load_lgb():
    """LightGBM 需要 OpenMP shim (否则 FileNotFoundError: lib_lightgbm.dll)。"""
    import sklearn
    for d in _glob.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll")):
        try:
            ctypes.CDLL(d)
        except Exception:  # noqa: BLE001
            pass
    import lightgbm as lgb
    return lgb


def get_metrics_module():
    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    return mr


def load_phys_names() -> tuple[list[str], list[str]]:
    spec = importlib.util.spec_from_file_location("t8_v25",
                                                  ROOT / "scripts" / "08_train_v1.0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    phys = ([c for c in mod.V02_FEATURES if c not in mod.REPLACED] + list(mod.V03_NEW))
    return phys, list(mod.COUNT_COLS)


def build_metric_fns(mr: Any) -> dict[str, Callable]:
    return {
        "spearman": mr.spearman_rho,
        "kendall_tau": mr.kendall_tau,
        "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=TOP_FRAC),
        "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, TOP_FRAC),
        "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, TOP_FRAC),
        "ef@5%": lambda a, b: mr.enrichment_factor(a, b, 5.0, TOP_FRAC),
    }


def bootstrap_ci(fn: Callable, y: np.ndarray, p: np.ndarray, n_boot: int,
                 seed: int = SEED) -> tuple[float, float, float]:
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    pt = float(fn(y, p)) if len(y) else float("nan")
    if len(y) < 5:
        return float("nan"), float("nan"), pt
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n_boot):
        idx = rs.randint(0, len(y), len(y))
        v = fn(y[idx], p[idx])
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return float("nan"), float("nan"), pt
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), pt


def paired_bootstrap(fn: Callable, y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                     n_boot: int = N_BOOT_PAIR, seed: int = SEED) -> dict[str, float]:
    """配对 bootstrap: Δ = metric(B) − metric(A)。"""
    ok = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
    y, pa, pb = y[ok], pa[ok], pb[ok]
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y)), "n_distinct_deltas": 0}
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
                "n": int(len(y)), "n_distinct_deltas": 0}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    pl = (np.sum(d <= 0) + 1) / (len(d) + 1)
    pr = (np.sum(d >= 0) + 1) / (len(d) + 1)
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl, pr))),
            "metric_a": ma, "metric_b": mb, "n": int(len(y)),
            "n_distinct_deltas": int(len(np.unique(np.round(d, 10))))}


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


def load_aligned(use_perres_pooled: bool = True) -> dict[str, Any]:
    """载入**同一批行**上的物理 / pooled / per-residue 特征。

    返回 dict:
      df_v22   : 对齐后的 V2.2 行 (6885 行, 含 real_dG / split / interface_key)
      Xp       : 物理特征 (n, F) float32
      Xp_cols  : 物理特征列名
      pooled   : {'lig': (n,640), 'rec': (n,640)} float32  (来自紧凑缓存, fp16->fp32)
      lig_topk / rec_topk / lig_mask / rec_mask : per-residue 下采样张量
      y        : ΔG (n,)
      tr/va/te : 布尔掩码
      splits   : 字符串数组
    """
    z = np.load(PERRES_CACHE, allow_pickle=True)
    v22_row = np.asarray(z["v22_row"])
    splits = np.asarray(z["splits"]).astype(str)
    d = pd.read_csv(V22, low_memory=False)
    dv = d.iloc[v22_row].reset_index(drop=True)
    # 三重一致性复核 (每次载入都查, 防止缓存被换掉)
    if not bool((dv["split"].astype(str).to_numpy() == splits).all()):
        raise RuntimeError("缓存与 V2.2 的 split 不一致 -> 拒绝使用")
    if not bool((dv["pdb_id"].astype(str).to_numpy()
                 == np.asarray(z["pdb_ids"]).astype(str)).all()):
        raise RuntimeError("缓存与 V2.2 的 pdb_id 不一致 -> 拒绝使用")

    phys, counts = load_phys_names()
    fcols = [c for c in phys if c in dv.columns]
    Xp = dv[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    tr = (splits == "train")
    med = Xp[tr].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)[fcols].to_numpy(np.float32)

    out: dict[str, Any] = {
        "df_v22": dv, "Xp": Xp, "Xp_cols": fcols,
        "y": dv["real_dG"].to_numpy(float), "splits": splits,
        "tr": tr, "va": (splits == "val"), "te": (splits == "test"),
        "pdb_ids": np.asarray(z["pdb_ids"]).astype(str),
        "lig_topk": z["lig_topk"], "rec_topk": z["rec_topk"],
        "lig_mask": np.asarray(z["lig_mask"]),
        "rec_mask": np.asarray(z["rec_mask"]),
    }
    if use_perres_pooled:
        out["pooled"] = {
            "lig": np.asarray(z["lig_pooled"], dtype=np.float32),
            "rec": np.asarray(z["rec_pooled"], dtype=np.float32),
        }
    return out


def load_pooled_direct() -> dict[str, Any]:
    """直接从 pooled ESM-2 缓存载入 (用于核对紧凑缓存的 pooled 是否一致)。"""
    z = np.load(POOLED, allow_pickle=True)
    return {"lig": z["lig_l2"].astype(np.float32),
            "rec": z["rec_l2"].astype(np.float32),
            "lig_md5": z["lig_md5"].astype(str), "rec_md5": z["rec_md5"].astype(str)}


def mean_pool_perres(data: dict[str, Any], side: str) -> np.ndarray:
    """对下采样后的 per-residue 张量做均值池化 (仅对有效残基求均值)。

    这是"per-residue mean-pool 对照"的核心操作: 与链级 pooled 的区别在于
    它平均的是**界面残基**, 而不是整条链。

    内存: **分块计算** —— 整块 float32 转换在 6885x32x640 时会申请 1.05 GiB,
    在本机 (7.5 GB 且空闲 <2 GB) 会 OOM。
    """
    return _chunked_residue_reduce(data, side, "mean")


def _chunked_residue_reduce(data: dict[str, Any], side: str,
                            how: str, chunk: int = 512,
                            extra: dict[str, np.ndarray] | None = None
                            ) -> np.ndarray | tuple:
    """分块遍历 per-residue 张量, 避免一次性 float32 化。

    ``how``: "mean" -> 均值池化; "stats" -> (mean, std, max, n) 逐条统计。
    """
    X = data[f"{side}_topk"]                      # (n, K, D) fp16 (npz 数组)
    M = data[f"{side}_mask"]                      # (n, K) bool
    n, K, D = X.shape
    outs: dict[str, list] = {}
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        xb = np.asarray(X[s:e], dtype=np.float32)  # 只转这一块
        mb = M[s:e]
        if how == "mean":
            w = mb[..., None].astype(np.float32)
            outs.setdefault("mean", []).append(
                (xb * w).sum(1) / np.maximum(w.sum(1), 1.0))
        elif how == "stats":
            mu = (xb * mb[..., None]).sum(1, keepdims=True) / np.maximum(
                mb.sum(1, keepdims=True)[..., None], 1.0)
            dev = np.linalg.norm(xb - mu, axis=2) / np.sqrt(D)
            dev = np.where(mb, dev, np.nan)
            nrm = np.linalg.norm(xb, axis=2) / np.sqrt(D)
            nrm = np.where(mb, nrm, np.nan)
            with np.errstate(all="ignore"):
                outs.setdefault("dev_mean", []).append(np.nanmean(dev, axis=1))
                outs.setdefault("dev_max", []).append(np.nanmax(dev, axis=1))
                outs.setdefault("dev_std", []).append(np.nanstd(dev, axis=1))
                outs.setdefault("norm_mean", []).append(np.nanmean(nrm, axis=1))
                outs.setdefault("norm_std", []).append(np.nanstd(nrm, axis=1))
            outs.setdefault("n_res", []).append(mb.sum(1).astype(float))
        del xb, mb
    res = {k: np.concatenate(v) for k, v in outs.items()}
    return res if how == "stats" else res["mean"]
