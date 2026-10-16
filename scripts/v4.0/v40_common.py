# -*- coding: utf-8 -*-
"""scripts/v4.0/v40_common.py -- V4.0 共享机器（配对数据 / 特征 / 指标 / 检验）.

口径声明
--------
* 任务: 预测 **ΔΔG = ΔG_mutant − ΔG_wild-type**（kcal/mol）
* 划分: `data/v4.0/ddg_split.json`（**重建**的复合物级 + 40% 序列同一性 + 家族级隔离）
  > 继承的 Proteinea split 在 ΔΔG 任务上退化（test 仅 27 对 / 10 复合物），已如实记录。
* 指标: Pearson（主，与文献可比）/ Spearman / RMSE / MAE /
  **per-structure Pearson**（SAMAffinity / Pythia-PPI 口径）/ AUROC（ΔΔG>0 vs <0）
* 检验: 配对 bootstrap 5000 次 + BH-FDR

特征（全部复用已有缓存，不重新提取）
------------------------------------
| 名称 | 维度 | 来源 |
|---|---|---|
| `phys` | 60 | `v2.2_proteinea_features.csv` 物理特征 |
| `Z` | 128 | `v2_esm2_cache.npz` pooled ESM-2 (640×2) → PCA128（**仅 train 拟合**）|
| `perres` | 128 | `v2.5_perres_cache.npz` per-residue mean-pool（可选）|

配对级特征构造: `[x_mut, x_wt, x_mut − x_wt]`（差分是 ΔΔG 的核心信号）。
"""

from __future__ import annotations

import ctypes
import glob as _glob
import json
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
DATA = ROOT / "data"
V40 = DATA / "v4.0"
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v4.0"
EXP = ROOT / "experiments" / "v4.0"
for _p in (V40, LOG_DIR, EXP, EXP / "baselines"):
    _p.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "failures.log"

RAW = DATA / "ppb_affinity_proteinea" / "filtered.csv"
V22 = FEAT / "v2.2_proteinea_features.csv"
POOLED = FEAT / "v2_esm2_cache.npz"
PERRES = FEAT / "v2.5_perres_cache.npz"
PAIRS = V40 / "ddg_pairs.csv"
SPLIT = V40 / "ddg_split.json"

SEED = 42
N_BOOT_PAIR = 5000
N_BOOT_CI = 800
PCA_DIM = 128
MIN_PS_N = 5          # per-structure Pearson 的最小对/结构（Pythia-PPI 用 10，此处放宽并报告）

LGB_PARAMS = dict(objective="regression", n_estimators=800, learning_rate=0.03,
                  num_leaves=15, max_depth=5, min_child_samples=20, subsample=0.8,
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


def load_phys_names() -> tuple[list[str], list[str]]:
    """复用 V1.0 的特征列定义（保证与 V2.2 表一致）。"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("t8_v40",
                                                  ROOT / "scripts" / "08_train_v1.0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    phys = ([c for c in mod.V02_FEATURES if c not in mod.REPLACED] + list(mod.V03_NEW))
    return phys, list(mod.COUNT_COLS)


# ---------------------------------------------------------------- metrics
def pearson(y: np.ndarray, p: np.ndarray) -> float:
    y, p = np.asarray(y, float), np.asarray(p, float)
    if len(y) < 3:
        return float("nan")
    ys, ps = y - y.mean(), p - p.mean()
    den = np.sqrt((ys ** 2).sum() * (ps ** 2).sum())
    return float((ys * ps).sum() / den) if den > 0 else float("nan")


def spearman(y: np.ndarray, p: np.ndarray) -> float:
    from scipy.stats import rankdata
    if len(y) < 3:
        return float("nan")
    return pearson(rankdata(y), rankdata(p))


def rmse(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(y, float) - np.asarray(p, float)) ** 2)))


def mae(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean(np.abs(np.asarray(y, float) - np.asarray(p, float))))


def auroc(y: np.ndarray, p: np.ndarray) -> float:
    """ΔΔG>0 (削弱结合) 作为正类的 AUROC。用秩统计实现，避免 sklearn 依赖歧义。"""
    y, p = np.asarray(y, float), np.asarray(p, float)
    lab = (y > 0).astype(int)
    n_pos, n_neg = int(lab.sum()), int((1 - lab).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    from scipy.stats import rankdata
    r = rankdata(p)
    return float((r[lab == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def per_structure_pearson(y: np.ndarray, p: np.ndarray, groups: np.ndarray,
                          min_n: int = MIN_PS_N) -> dict:
    """每个结构内单独算 Pearson 再平均（SAMAffinity / Pythia-PPI 口径）。"""
    vals, used, skipped = [], 0, 0
    for g in np.unique(groups):
        m = groups == g
        if int(m.sum()) < min_n:
            skipped += 1
            continue
        v = pearson(y[m], p[m])
        if v == v:
            vals.append(v)
            used += 1
        else:
            skipped += 1
    return {"per_structure_pearson": float(np.mean(vals)) if vals else float("nan"),
            "per_structure_n": used, "per_structure_skipped": skipped,
            "per_structure_values": vals}


def all_metrics(y: np.ndarray, p: np.ndarray, groups: np.ndarray | None = None,
                min_n: int = MIN_PS_N) -> dict:
    out = {"pearson": pearson(y, p), "spearman": spearman(y, p),
           "rmse": rmse(y, p), "mae": mae(y, p), "auroc": auroc(y, p), "n": int(len(y))}
    if groups is not None:
        ps = per_structure_pearson(y, p, groups, min_n)
        out.update({k: v for k, v in ps.items() if k != "per_structure_values"})
    return out


def bootstrap_ci(fn, y, p, n_boot: int = N_BOOT_CI, seed: int = SEED):
    y, p = np.asarray(y, float), np.asarray(p, float)
    pt = float(fn(y, p))
    if len(y) < 5:
        return float("nan"), float("nan"), pt
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n_boot):
        i = rs.randint(0, len(y), len(y))
        v = fn(y[i], p[i])
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return float("nan"), float("nan"), pt
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), pt


def paired_bootstrap(fn, y, pa, pb, n_boot: int = N_BOOT_PAIR, seed: int = SEED) -> dict:
    """配对 bootstrap: Δ = fn(y, pb) − fn(y, pa)。"""
    y, pa, pb = (np.asarray(x, float) for x in (y, pa, pb))
    y = y.astype(float)
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y))}
    ma, mb = float(fn(y, pa)), float(fn(y, pb))
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        i = rs.randint(0, len(y), len(y))
        va, vb = fn(y[i], pa[i]), fn(y[i], pb[i])
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
            "metric_a": ma, "metric_b": mb, "n": int(len(y)),
            "n_distinct_deltas": int(len(np.unique(np.round(d, 10))))}


def paired_permutation(fn, y, pa, pb, n_perm: int = 5000, seed: int = SEED) -> float:
    """配对置换检验（更保守的 p 值）：随机交换每个样本的 a/b 标签。"""
    y, pa, pb = (np.asarray(x, float) for x in (y, pa, pb))
    obs = fn(y, pb) - fn(y, pa)
    rs = np.random.RandomState(seed)
    cnt = 0
    for _ in range(n_perm):
        sw = rs.rand(len(y)) < 0.5
        a2 = np.where(sw, pb, pa)
        b2 = np.where(sw, pa, pb)
        if (fn(y, b2) - fn(y, a2)) >= obs:
            cnt += 1
    return float((cnt + 1) / (n_perm + 1))


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([1.0 if (x != x) else float(x) for x in pvals], float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p, kind="mergesort")
    q = np.empty(m, float)
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    return [float(min(1.0, v)) for v in q]


# ---------------------------------------------------------------- data
def standardize_train(X: np.ndarray, tr: np.ndarray) -> tuple[np.ndarray, dict]:
    """**仅在训练行上**拟合 median + std，防止信息从 val/test 泄漏。

    必要性（实测）: 物理特征 std 可达 293、绝对值到 1417，而 PCA 嵌入 std 仅
    0.04–0.05 —— 相差 4 个数量级。不标准化时 MLP 直接梯度爆炸输出 NaN。
    """
    X = np.asarray(X, np.float32)
    med = np.nanmedian(X[tr], axis=0)
    med = np.where(np.isfinite(med), med, 0.0).astype(np.float32)
    Xf = np.where(np.isnan(X), med[None, :], X)
    sd = Xf[tr].std(axis=0)
    sd = np.where((sd > 0) & np.isfinite(sd), sd, 1.0).astype(np.float32)
    return ((Xf - med[None, :]) / sd[None, :]).astype(np.float32), {
        "median": med, "std": sd}


def standardize_apply(X: np.ndarray, scaler: dict) -> np.ndarray:
    """用已拟合的 scaler 变换（必须先做 median 填充，否则 NaN 会残留）。"""
    X = np.asarray(X, np.float32)
    med, sd = scaler["median"], scaler["std"]
    Xf = np.where(np.isnan(X), med[None, :], X)
    return ((Xf - med[None, :]) / sd[None, :]).astype(np.float32)


def build_pair_features(use_perres: bool = False, verbose: bool = True) -> dict[str, Any]:
    """装配配对级特征。返回 dict 含 X 各块、y、groups、split 索引。"""
    log = (lambda m: print(m)) if verbose else (lambda m: None)

    pairs = pd.read_csv(PAIRS, low_memory=False)
    fwd = pairs[pairs["direction"] == "forward"].reset_index(drop=True)
    split = json.loads(SPLIT.read_text(encoding="utf-8"))
    # 01 写入的 pair_index 指向 fwd 的行号
    assert set(split["pair_index"]) == {"train", "val", "test"}
    n_fwd = len(fwd)
    for s, idx in split["pair_index"].items():
        assert max(idx) < n_fwd, f"{s} 索引越界"

    v22 = pd.read_csv(V22, low_memory=False)
    v22["pdb_id"] = v22["pdb_id"].astype(str)
    raw = pd.read_csv(RAW, low_memory=False)
    raw["PDB"] = raw["PDB"].astype(str)

    # --- 把 ddg_pairs 的 raw 行号映射到 v22 的行号 ---
    # ddg_pairs 的 wt_row/mut_row 是 filtered.csv 的原始行号
    raw_key = raw["PDB"].astype(str) + "|" + raw["KD(M)"].astype(str) + "|" + \
        raw["Mutations"].fillna("WT").astype(str)
    v22_key = v22["pdb_id"].astype(str) + "|" + v22["KD(M)"].astype(str) + "|" + \
        v22["Mutations"].fillna("WT").astype(str)
    key2v22: dict[str, int] = {}
    for i, k in enumerate(v22_key):
        key2v22.setdefault(k, i)
    wt_v22 = np.array([key2v22.get(raw_key.iloc[r], -1) for r in fwd["wt_row"]])
    mut_v22 = np.array([key2v22.get(raw_key.iloc[r], -1) for r in fwd["mut_row"]])
    n_map = int(((wt_v22 >= 0) & (mut_v22 >= 0)).sum())
    log(f"- ddg_pairs -> v2.2 行映射: **{n_map} / {n_fwd}** "
        f"({n_map / n_fwd * 100:.2f}%)")
    keep = (wt_v22 >= 0) & (mut_v22 >= 0)
    if not keep.all():
        log(f"  !! 丢弃 {int((~keep).sum())} 对（v2.2 中找不到对应行），如实记录")
    fwd = fwd[keep].reset_index(drop=True)
    wt_v22, mut_v22 = wt_v22[keep], mut_v22[keep]
    # split 索引也要同步收缩
    old_ok = np.where(keep)[0]
    remap = {int(o): int(n) for n, o in enumerate(old_ok)}
    split_idx = {s: [remap[i] for i in idx if i in remap]
                 for s, idx in split["pair_index"].items()}

    # --- 物理特征 ---
    phys, counts = load_phys_names()
    fcols = [c for c in phys if c in v22.columns]
    Xp_all = v22[fcols].apply(pd.to_numeric, errors="coerce")
    for c in counts:
        if c in Xp_all.columns:
            Xp_all[c] = np.log1p(Xp_all[c].clip(lower=0))
    Xp_all = Xp_all[fcols].to_numpy(np.float32)

    # --- pooled ESM-2 -> PCA128（仅 train 拟合）---
    z = np.load(POOLED, allow_pickle=True)
    Lig = z["lig_l2"].astype(np.float32)
    Rec = z["rec_l2"].astype(np.float32)
    both = np.hstack([Lig, Rec])
    if len(both) != len(v22):
        raise RuntimeError(f"pooled 缓存 {len(both)} != v2.2 表 {len(v22)}")

    # --- 组装 ---
    tr_rows = split_idx["train"]
    widx = np.unique(np.concatenate([wt_v22[tr_rows], mut_v22[tr_rows]]))
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(PCA_DIM, both.shape[1], len(widx) - 1),
              random_state=SEED)
    pca.fit(both[widx])
    Z_all = pca.transform(both).astype(np.float32)
    log(f"- PCA{PCA_DIM} 解释方差 {float(np.sum(pca.explained_variance_ratio_)):.4f}"
        f"（仅在 train 行的并集上拟合）")

    P_wt, P_mut = Xp_all[wt_v22], Xp_all[mut_v22]
    Z_wt, Z_mut = Z_all[wt_v22], Z_all[mut_v22]

    blocks: dict[str, np.ndarray] = {
        "phys_mut": P_mut, "phys_wt": P_wt, "phys_diff": P_mut - P_wt,
        "Z_mut": Z_mut, "Z_wt": Z_wt, "Z_diff": Z_mut - Z_wt,
    }
    if use_perres:
        pr = np.load(PERRES, allow_pickle=True)
        # perres 缓存按 v22 的 usable 行索引（v22_row）
        # 这里保守处理: 用 npz 里的 v22_row 建立映射
        raise NotImplementedError("per-residue 配对特征需单独实现（见 04_baselines 的 B4）")

    y = fwd["ddg"].to_numpy(float)
    groups = fwd["pdb_id"].astype(str).to_numpy()
    out = {"X_blocks": blocks, "y": y, "groups": groups, "df": fwd,
           "split_idx": split_idx, "fwd_v22_wt": wt_v22, "fwd_v22_mut": mut_v22,
           "Z_all": Z_all, "Xp_all": Xp_all, "pca": pca,
           "n_usable": n_fwd, "n_mapped": n_map}
    for s in ("train", "val", "test"):
        out[f"mask_{s}"] = np.zeros(len(y), bool)
        out[f"mask_{s}"][split_idx[s]] = True
    return out
