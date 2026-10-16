# -*- coding: utf-8 -*-
"""scripts/v3.1/v31_common.py -- V3.1 共享机器（数据 / 指标 / 检验）.

**口径统一声明（硬约束 12）**
----------------------------
* 评估一律用 **Proteinea 固定划分**（train 6485 / val 965 / test 757，`row_structure_usable==1`
  过滤后为 6482 / 965 / 752）。
* 该划分已由 `01_verify_split.py` + `02_exact_seq_overlap.py` **独立复核**：
  * 结构级：train ∩ test 的 pdb_id = **0**；
  * 序列级：**两侧都 exact-match 的 test 行 = 0 / 757**；
  * （`01` 报的 22.93% >0.30 Jaccard 是 3-mer **组成相似**，不构成泄漏证据 ——
    3-mer Jaccard 会高估真实同一性。）
* 参照值：V2.3.2 的 `physics+esm2_both@l2@128` = **0.3889577363940227**（Spearman）。

CPU only；LightGBM 需先加载 OpenMP shim。
"""

from __future__ import annotations

import ctypes
import glob as _glob
import importlib.util
import pickle
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
MODEL_DIR = ROOT / "models"
LOG_DIR = ROOT / "logs" / "v3.1"
LOG_DIR.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "failures.log"

V22 = FEAT / "v2.2_proteinea_features.csv"
POOLED = FEAT / "v2_esm2_cache.npz"

SEED = 42
TOP_FRAC = 0.10
N_BOOT_PAIR = 5000
N_BOOT_CI = 800
PCA_DIM = 128

# V2.3.2 参照（同一固定划分上的历史最优）
REF_RHO = 0.3889577363940227
REF_LABEL = "V2.3.2 physics+esm2_both@l2@128"

LGB_PARAMS = dict(objective="regression", n_estimators=500, learning_rate=0.03,
                  num_leaves=15, max_depth=5, min_child_samples=12, subsample=0.8,
                  subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                  n_jobs=2, random_state=SEED, verbose=-1)


def load_lgb():
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
    spec = importlib.util.spec_from_file_location("t8_v31",
                                                  ROOT / "scripts" / "08_train_v1.0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    phys = ([c for c in mod.V02_FEATURES if c not in mod.REPLACED] + list(mod.V03_NEW))
    return phys, list(mod.COUNT_COLS)


def build_metric_fns(mr: Any) -> dict:
    return {
        "spearman": mr.spearman_rho,
        "kendall_tau": mr.kendall_tau,
        "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=TOP_FRAC),
        "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, TOP_FRAC),
        "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, TOP_FRAC),
        "ef@5%": lambda a, b: mr.enrichment_factor(a, b, 5.0, TOP_FRAC),
    }


def bootstrap_ci(fn, y: np.ndarray, p: np.ndarray, n_boot: int = N_BOOT_CI,
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


def paired_bootstrap(fn, y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                     n_boot: int = N_BOOT_PAIR, seed: int = SEED) -> dict:
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


def load_data(pca_dim: int = PCA_DIM, use_pooled: bool = True) -> dict[str, Any]:
    """载入固定划分的数据 + 物理特征 + pooled ESM-2 的 PCA 表示。

    返回 dict 含: y, tr/va/te 掩码, Xp(物理), Z(pooled PCA), groups(interface_key),
    labels(split 字符串), mutation 等元数据。
    """
    d = pd.read_csv(V22, low_memory=False)
    d = d[d["row_structure_usable"] == 1].reset_index(drop=True)
    y = d["real_dG"].to_numpy(float)
    tr = (d["split"] == "train").to_numpy()
    va = (d["split"] == "val").to_numpy()
    te = (d["split"] == "test").to_numpy()

    phys, counts = load_phys_names()
    fcols = [c for c in phys if c in d.columns]
    Xp = d[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    med = Xp[tr].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)[fcols].to_numpy(np.float32)

    out: dict[str, Any] = {
        "df": d, "y": y, "tr": tr, "va": va, "te": te,
        "Xp": Xp, "Xp_cols": fcols,
        "groups": d["interface_key"].astype(str).to_numpy(),
        "pdb": d["pdb_id"].astype(str).to_numpy(),
    }
    if use_pooled:
        z = np.load(POOLED, allow_pickle=True)
        # 缓存与 V2.2 全表同序（V2.3.2/V2.5 均已用 md5 逐行确认）
        if len(z["pdb_ids"]) != len(pd.read_csv(V22, low_memory=False)):
            raise RuntimeError("pooled 缓存与 V2.2 全表行数不一致")
        Lig = z["lig_l2"].astype(np.float32)
        Rec = z["rec_l2"].astype(np.float32)
        both = np.hstack([Lig, Rec])
        # 只保留可用行 -> 需要用未过滤表的行索引
        d_all = pd.read_csv(V22, low_memory=False)
        keep = (d_all["row_structure_usable"] == 1).to_numpy()
        both = both[keep]
        from sklearn.decomposition import PCA
        p = PCA(n_components=min(pca_dim, both.shape[1], int(tr.sum()) - 1),
                random_state=SEED)
        p.fit(both[tr])                      # 只在训练集 fit
        out["Z"] = p.transform(both).astype(np.float32)
        out["Lig"] = Lig[keep]
        out["Rec"] = Rec[keep]
        out["evr"] = float(np.sum(p.explained_variance_ratio_))
    return out


def ref_features(data: dict) -> np.ndarray:
    """V2.3.2 参照配置的特征矩阵（物理 + pooled PCA128）。"""
    return np.hstack([data["Xp"], data["Z"]]).astype(np.float32)


def fit_eval_lgb(lgb, X: np.ndarray, data: dict, tag: str, mfns: dict,
                 seed: int = SEED) -> dict:
    """在固定划分上训练 + 评估。返回指标 dict（含 test 预测）。"""
    y, tr, va, te = data["y"], data["tr"], data["va"], data["te"]
    m = lgb.LGBMRegressor(**{**LGB_PARAMS, "random_state": seed})
    m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
          callbacks=[lgb.early_stopping(50, verbose=False)])
    pred = m.predict(X[te])
    y_s, p_s = -y[te], -pred
    out = {"tag": tag, "n_features": int(X.shape[1]),
           "best_iteration": int(getattr(m, "best_iteration_", -1) or -1)}
    for k, fn in mfns.items():
        try:
            out[k] = float(fn(y_s, p_s))
        except Exception:  # noqa: BLE001
            out[k] = float("nan")
    for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
        lo, hi, _ = bootstrap_ci(mfns[k], y_s, p_s, N_BOOT_CI, seed)
        out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
    out["_pred"] = pred
    out["_model"] = m
    return out
