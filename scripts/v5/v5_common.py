# -*- coding: utf-8 -*-
"""scripts/v5/v5_common.py -- V5.0 共享机器.

设计原则（对应 V5 总 prompt）
-----------------------------
* **只读** V1.1–V4.0 的任何文件；V5 新产物全部写 `experiments/v5_0/` 或 `scripts/v5/`。
* **不使用 test 做任何选择**：PCA / 标准化 / 超参 / 早停 / conformal 全部只用 train
  （conformal 额外只用 train 的 OOF）。
* **cluster bootstrap** 为主推断（V5 §12）：抽样单位 = family（不足则 PDB），
  携带该 cluster 的全部 mutation rows。row-level bootstrap 仅作 legacy 次级分析。
* **BH-FDR** 多重比较 + practical significance 阈值（V5 §13）。
* **决策标签**自动生成（V5 §101）：KEEP / DROP / INCONCLUSIVE / INVALID / RESOURCE_BLOCKED。

数据口径（已完成 PHASE 0 审计）
-------------------------------
* 配对表: `data/v4.0/ddg_pairs.csv`（15052 行 = 7526 forward + 7526 reverse）
* 划分: `data/v4.0/ddg_split.json`（三重隔离；train 5268 / val 1129 / test 1129）
* 特征: 564 维 = [物理60 + pooled ESM-2 PCA128] × 3 块 (x_mut / x_wt / x_mut−x_wt)
* 特征装配与 V4 的 `scripts/v4.0/v40_common.build_pair_features` **逐位一致**
  （V5 直接复用该函数的装配逻辑，但**重新实现**以避免 V5 依赖 V4 代码路径。
  若两者不一致，baseline reproduction gate 会报错。）
"""

from __future__ import annotations

import ctypes
import glob as _glob
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
V5 = ROOT / "experiments" / "v5_0"
DATA = ROOT / "data"
FEAT = ROOT / "features"
for _p in (V5, V5 / "logs", V5 / "manifests", V5 / "configs", V5 / "figures",
           V5 / "audit"):
    _p.mkdir(parents=True, exist_ok=True)

V22 = FEAT / "v2.2_proteinea_features.csv"
POOLED = FEAT / "v2_esm2_cache.npz"
PAIRS = DATA / "v4.0" / "ddg_pairs.csv"
SPLIT = DATA / "v4.0" / "ddg_split.json"

SEED = 42
SEEDS = (42, 43, 44)
PCA_DIM = 128
MIN_PS_N = 5                 # per-structure 指标的最小对/结构（与 V4 一致以便可比）
N_BOOT_CLUSTER = 5000        # V5 主推断
N_BOOT_ROW = 5000            # legacy 次级
N_BOOT_CI = 2000

# V4 记录值（来自 experiments/v4.0/innovation*_results.csv，已由 PHASE 0 哈希固定）
V4_C_DGDG = {42: 0.395133, 43: 0.384308, 44: 0.379019}
V4_C_SYM = {42: 0.367427, 43: 0.343410, 44: 0.349652}
V4_A_FULL = {42: 0.250352, 43: 0.168290, 44: 0.260540}

LGB_PARAMS = dict(objective="regression", n_estimators=800, learning_rate=0.03,
                  num_leaves=15, max_depth=5, min_child_samples=20, subsample=0.8,
                  subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                  n_jobs=2, verbose=-1)

# practical significance 阈值（V5 §13/§102/§103）
PRACTICAL_DPEARSON = 0.02
PRACTICAL_DMACRO_SPEARMAN = 0.04
PRACTICAL_DP_AT_10 = 0.10
PRACTICAL_EF1_REL = 0.20
PRACTICAL_DHIT10 = 0.10


def log(msg: str = "") -> None:
    print(msg)
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha256_file(p: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_lgb():
    """LightGBM 4.7.0 需要 OpenMP 运行库先加载（否则 FileNotFoundError）。"""
    import sklearn
    for d in _glob.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll")):
        try:
            ctypes.CDLL(d)
        except Exception:  # noqa: BLE001
            pass
    import lightgbm as lgb
    return lgb


def load_phys_names() -> tuple[list[str], list[str]]:
    """复用 V1.0 的特征列定义（与 V2.2 表口径一致）。"""
    spec = importlib.util.spec_from_file_location(
        "t8_v5", ROOT / "scripts" / "08_train_v1.0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    phys = ([c for c in mod.V02_FEATURES if c not in mod.REPLACED]
            + list(mod.V03_NEW))
    return phys, list(mod.COUNT_COLS)


# ============================================================ 数据装配
def build_pair_features(pca_dim: int = PCA_DIM, verbose: bool = True) -> dict[str, Any]:
    """装配 V4-口径的配对级特征（与 V4 的 v40_common.build_pair_features 同逻辑）。

    返回 dict: X_blocks（8 块）、y、groups(pdb_id)、df、split_idx、mask_*。
    """
    logf = log if verbose else (lambda *_a, **_k: None)

    pairs = pd.read_csv(PAIRS, low_memory=False)
    fwd = pairs[pairs["direction"] == "forward"].reset_index(drop=True)
    split = json.loads(SPLIT.read_text(encoding="utf-8"))
    n_fwd = len(fwd)
    for s, idx in split["pair_index"].items():
        assert max(idx) < n_fwd, f"{s} 索引越界"

    v22 = pd.read_csv(V22, low_memory=False)
    v22["pdb_id"] = v22["pdb_id"].astype(str)
    raw = pd.read_csv(DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    raw["PDB"] = raw["PDB"].astype(str)

    raw_key = (raw["PDB"].astype(str) + "|" + raw["KD(M)"].astype(str) + "|"
               + raw["Mutations"].fillna("WT").astype(str))
    v22_key = (v22["pdb_id"].astype(str) + "|" + v22["KD(M)"].astype(str) + "|"
               + v22["Mutations"].fillna("WT").astype(str))
    key2v22: dict[str, int] = {}
    for i, k in enumerate(v22_key):
        key2v22.setdefault(k, i)
    wt_v22 = np.array([key2v22.get(raw_key.iloc[r], -1) for r in fwd["wt_row"]])
    mut_v22 = np.array([key2v22.get(raw_key.iloc[r], -1) for r in fwd["mut_row"]])
    keep = (wt_v22 >= 0) & (mut_v22 >= 0)
    n_map = int(keep.sum())
    logf(f"- ddg_pairs -> v2.2 行映射: **{n_map} / {n_fwd}** "
         f"({n_map / n_fwd * 100:.2f}%)")
    if not keep.all():
        logf(f"  !! 丢弃 {int((~keep).sum())} 对（v2.2 无对应行），如实记录")
    fwd = fwd[keep].reset_index(drop=True)
    wt_v22, mut_v22 = wt_v22[keep], mut_v22[keep]
    old_ok = np.where(keep)[0]
    remap = {int(o): int(n) for n, o in enumerate(old_ok)}
    split_idx = {s: [remap[i] for i in idx if i in remap]
                 for s, idx in split["pair_index"].items()}

    phys, counts = load_phys_names()
    fcols = [c for c in phys if c in v22.columns]
    Xp_all = v22[fcols].apply(pd.to_numeric, errors="coerce")
    for c in counts:
        if c in Xp_all.columns:
            Xp_all[c] = np.log1p(Xp_all[c].clip(lower=0))
    Xp_all = Xp_all[fcols].to_numpy(np.float32)

    z = np.load(POOLED, allow_pickle=True)
    both = np.hstack([z["lig_l2"].astype(np.float32),
                      z["rec_l2"].astype(np.float32)])
    if len(both) != len(v22):
        raise RuntimeError(f"pooled 缓存 {len(both)} != v2.2 表 {len(v22)}")

    tr_rows = split_idx["train"]
    widx = np.unique(np.concatenate([wt_v22[tr_rows], mut_v22[tr_rows]]))
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(pca_dim, both.shape[1], len(widx) - 1),
              random_state=SEED)
    pca.fit(both[widx])                     # **只用 train 行**
    Z_all = pca.transform(both).astype(np.float32)
    logf(f"- PCA{pca_dim} 解释方差 {float(np.sum(pca.explained_variance_ratio_)):.4f}"
         f"（仅在 train 行的并集上拟合）")

    P_wt, P_mut = Xp_all[wt_v22], Xp_all[mut_v22]
    Z_wt, Z_mut = Z_all[wt_v22], Z_all[mut_v22]
    blocks = {"phys_mut": P_mut, "phys_wt": P_wt, "phys_diff": P_mut - P_wt,
              "Z_mut": Z_mut, "Z_wt": Z_wt, "Z_diff": Z_mut - Z_wt}
    y = fwd["ddg"].to_numpy(float)
    out = {"X_blocks": blocks, "y": y,
           "groups": fwd["pdb_id"].astype(str).to_numpy(), "df": fwd,
           "split_idx": split_idx, "fwd_v22_wt": wt_v22, "fwd_v22_mut": mut_v22,
           "Z_all": Z_all, "Xp_all": Xp_all, "pca": pca,
           "n_usable": n_fwd, "n_mapped": n_map,
           "family_of_complex": split.get("family_of_complex", {})}
    for s in ("train", "val", "test"):
        m = np.zeros(len(y), bool)
        m[split_idx[s]] = True
        out[f"mask_{s}"] = m
    # 特征矩阵（与 V4 同序）
    out["X"] = np.hstack([blocks["phys_mut"], blocks["phys_wt"], blocks["phys_diff"],
                          blocks["Z_mut"], blocks["Z_wt"], blocks["Z_diff"]
                          ]).astype(np.float32)
    return out


def standardize_train(X: np.ndarray, tr: np.ndarray) -> tuple[np.ndarray, dict]:
    """仅在训练行上拟合 median+std（防 val/test 泄漏）。"""
    X = np.asarray(X, np.float32)
    med = np.nanmedian(X[tr], axis=0)
    med = np.where(np.isfinite(med), med, 0.0).astype(np.float32)
    Xf = np.where(np.isnan(X), med[None, :], X)
    sd = Xf[tr].std(axis=0)
    sd = np.where((sd > 0) & np.isfinite(sd), sd, 1.0).astype(np.float32)
    return ((Xf - med[None, :]) / sd[None, :]).astype(np.float32), {"median": med, "std": sd}


# ============================================================ 指标
def pearson(y, p) -> float:
    y, p = np.asarray(y, float), np.asarray(p, float)
    if len(y) < 3:
        return float("nan")
    ys, ps = y - y.mean(), p - p.mean()
    den = np.sqrt((ys ** 2).sum() * (ps ** 2).sum())
    return float((ys * ps).sum() / den) if den > 0 else float("nan")


def spearman(y, p) -> float:
    from scipy.stats import rankdata
    if len(y) < 3:
        return float("nan")
    return pearson(rankdata(y), rankdata(p))


def rmse(y, p) -> float:
    return float(np.sqrt(np.mean((np.asarray(y, float) - np.asarray(p, float)) ** 2)))


def mae(y, p) -> float:
    return float(np.mean(np.abs(np.asarray(y, float) - np.asarray(p, float))))


def auroc(y, p) -> float:
    """ΔΔG>0（削弱结合）为正类的 AUROC。"""
    y, p = np.asarray(y, float), np.asarray(p, float)
    lab = (y > 0).astype(int)
    n_pos, n_neg = int(lab.sum()), int((1 - lab).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    from scipy.stats import rankdata
    r = rankdata(p)
    return float((r[lab == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def per_structure(y, p, groups, min_n: int = MIN_PS_N) -> dict:
    """每结构内单独算 Pearson/Spearman 再 macro 平均。"""
    pv, sv, used, skipped = [], [], 0, 0
    for g in np.unique(groups):
        m = groups == g
        if int(m.sum()) < min_n:
            skipped += 1
            continue
        a, b = pearson(y[m], p[m]), spearman(y[m], p[m])
        if a == a:
            pv.append(a)
        if b == b:
            sv.append(b)
        used += 1
    return {"macro_pdb_pearson": float(np.mean(pv)) if pv else float("nan"),
            "macro_pdb_spearman": float(np.mean(sv)) if sv else float("nan"),
            "per_structure_n": used, "per_structure_skipped": skipped}


def topk_metrics(y, p, groups=None, top_frac: float = 0.10) -> dict:
    """Top-K 富集指标。**ΔΔG 越负越好**，故统一的 score 空间为 `s = −ΔΔG`。

    指标定义（**与 V2.5/V2.3.2 的 P@k 口径对齐**，便于跨版本比较）:
      * `best_true` = 真值最负的前 `top_frac`（默认 10%）→ 「真强结合体」
      * `p_at_k` = 预测 top-k 中落在 `best_true` 内的比例
        → **这就是 V2.5/V2.3.2 报的 Precision@k**（V4 记录 `C+DGaux` 的 P@10 ≈ 0.20）
      * `ef_xpct` = 预测 top-x% 中命中 `best_true` 的数量 / 随机期望
      * `hit k` = 预测 top-k 中 **ΔΔG<0 / <−0.5 / <−1.0** 的比例
        （「有用突变」口径，与 `p_at_k` 不同；比 base rate 更有意义）
    """
    y, p = np.asarray(y, float), np.asarray(p, float)
    n = len(y)
    out: dict[str, float] = {}
    if n < 10:
        return out
    ys, ps = -y, -p                                  # 越大越好
    best_true = ys >= np.quantile(ys, 1.0 - top_frac)
    n_best = int(best_true.sum())
    out["base_rate_best_true"] = float(n_best / n)   # 随机选中的期望比例

    def _top_idx(k):
        return np.argsort(-ps, kind="mergesort")[:k]

    for k in (5, 10):
        if k > n:
            continue
        idx = _top_idx(k)
        out[f"p_at_{k}"] = float(best_true[idx].mean())        # ← 对齐 V2/V4 口径
        out[f"hit{k}_beneficial"] = float((y[idx] < 0).mean())
        out[f"hit{k}_ddg_m05"] = float((y[idx] < -0.5).mean())
        out[f"hit{k}_ddg_m10"] = float((y[idx] < -1.0).mean())
        out[f"sel{k}_mean_ddg"] = float(y[idx].mean())         # 越负越好
        out[f"sel{k}_best_ddg"] = float(y[idx].min())
        out[f"mean_ddg_all"] = float(y.mean())
        out[f"oracle{k}_mean_ddg"] = float(np.sort(y)[:k].mean())  # 全池 oracle
    for pct, name in ((1.0, "1pct"), (5.0, "5pct")):
        k = max(1, int(round(n * pct / 100.0)))
        idx = _top_idx(k)
        hits = int(best_true[idx].sum())
        exp = k * (n_best / n) if n else 0.0
        out[f"ef_{name}"] = float(hits / exp) if exp > 0 else float("nan")
        out[f"ef_{name}_hits"] = hits
        out[f"ef_{name}_k"] = k
    # Target-level utility（V5 §108/§76）：**每个 complex 内** 选 top-K 再 macro 平均
    #
    # ⚠️ BUG 修正记录（V5 §123 Bug protocol，2026-10-13）
    # ---------------------------------------------------
    # 原实现用 `orc = np.argsort(sub_y)[:k]` —— **升序**排列 sub_y，
    # 取前 k 个得到的是 ΔΔG **最大（最差）** 的候选，而非 oracle（最负最优）。
    # 后果: `efficiency_at_10` **被系统性反号**。合成自检证据（y = -3..8, k=10）:
    #   完美 oracle 模型 → efficiency_at_10 = **−0.6667**（应为 +1.0）
    #   完全反向模型   → efficiency_at_10 = **+1.0000**（应为 ≤0）
    # 处置: 该列在所有已产出 CSV 中标 `INVALID_DUE_TO_BUG` 并**全部重算**。
    # 修正: oracle 必须取 sub_y **最小**的 k 个 → `np.argsort(sub_y)[:k]`
    #       （升序取前 k = 最小 k 个，正确）；原 bug 实为 sel/orc 语义颠倒的
    #       组合效应，现统一为「越小越好」的显式字典序比较。
    if groups is not None:
        eff, acc_sel, acc_or, acc_rnd, n_cx = [], [], [], [], 0
        rs = np.random.RandomState(SEED)
        n_rand = 200          # 随机基线用多次抽样（单次抽样在小组内方差太大）
        for g in np.unique(groups):
            m = np.where(groups == g)[0]
            if len(m) < 5:
                continue
            k = min(10, len(m)) if len(m) >= 10 else min(5, len(m))
            sub_y, sub_p = y[m], p[m]
            # 预测最优 = 预测值最小（ΔΔG 越负越好）
            sel = np.argsort(sub_p, kind="mergesort")[:k]
            # oracle = 真值最小
            orc = np.argsort(sub_y, kind="mergesort")[:k]
            # 随机基线 = 多次随机抽样的平均（降低小组方差）
            rr = [rs.permutation(len(m))[:k] for _ in range(n_rand)]
            v_sel = float(-sub_y[sel].mean())
            v_or = float(-sub_y[orc].mean())
            v_rnd = float(np.mean([-sub_y[r].mean() for r in rr]))
            # Efficiency@K = 相对随机回收了多少 oracle 空间
            # 用 (random − model)/(random − oracle)：对合理模型恒 ≥ 0，
            # oracle → +1，random → 0（与 V5 §108 的语义一致）
            den = v_rnd - v_or
            eff.append((v_rnd - v_sel) / den if abs(den) > 1e-9 else float("nan"))
            acc_sel.append(v_sel); acc_or.append(v_or); acc_rnd.append(v_rnd)
            n_cx += 1
        eff = [e for e in eff if e == e]
        if eff:
            out["efficiency_at_10"] = float(np.mean(eff))
            out["efficiency_at_10_median"] = float(np.median(eff))
            out["n_complex_eff"] = n_cx
        if acc_sel:
            out["utility_model_mean"] = float(np.mean(acc_sel))
            out["utility_oracle_mean"] = float(np.mean(acc_or))
            out["utility_random_mean"] = float(np.mean(acc_rnd))
            out["regret_at_10"] = float(np.mean(acc_or) - np.mean(acc_sel))
    return out


def per_structure_spearman(y, p, groups, min_n: int = MIN_PS_N) -> float:
    """签名兼容 cluster bootstrap（3 参数）：macro per-structure Spearman。"""
    return per_structure(y, p, groups, min_n)["macro_pdb_spearman"]


def per_structure_pearson(y, p, groups, min_n: int = MIN_PS_N) -> float:
    """签名兼容 cluster bootstrap（3 参数）：macro per-structure Pearson。"""
    return per_structure(y, p, groups, min_n)["macro_pdb_pearson"]


def all_metrics(y, p, groups=None, with_topk: bool = True) -> dict:
    out = {"n": int(len(y)), "overall_pearson": pearson(y, p),
           "overall_spearman": spearman(y, p), "rmse": rmse(y, p),
           "mae": mae(y, p), "auroc": auroc(y, p)}
    if groups is not None:
        out.update(per_structure(y, p, groups))
    if with_topk:
        out.update(topk_metrics(y, p, groups))
    return out


# ============================================================ 统计推断
def _cluster_index(groups: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    """返回 (unique_clusters, [每个 cluster 的行索引数组])。"""
    order = np.argsort(groups, kind="mergesort")
    gs = groups[order]
    uniq, starts = np.unique(gs, return_index=True)
    idxs = []
    for i, s in enumerate(starts):
        e = starts[i + 1] if i + 1 < len(starts) else len(gs)
        idxs.append(order[s:e])
    return uniq, idxs


def cluster_bootstrap_ci(fn, y, p, groups, n_boot: int = N_BOOT_CI,
                         seed: int = SEED) -> dict:
    """**cluster bootstrap**：每次抽 cluster（family/PDB），携带其全部 rows。

    这是 V5 主推断口径（V5 §12）：同一 PDB 的 mutation 相关，row-level bootstrap
    会低估 CI。

    `fn` 的签名必须能接受 `groups`（第三参数），以便 per-structure 类指标在
    重采样后按 **重采样后的 cluster 标签** 重新计算（否则 index 会错位）。
    为兼容只接受 (y, p) 的 fn，这里做一次 arity 探测。
    """
    import inspect
    y, p = np.asarray(y, float), np.asarray(p, float)
    groups = np.asarray(groups)
    uniq, cidx = _cluster_index(groups)
    n_single = int(sum(1 for ix in cidx if len(ix) == 1))
    try:
        nparams = len(inspect.signature(fn).parameters)
    except (TypeError, ValueError):
        nparams = 2
    takes_groups = nparams >= 3

    def _eval(yy, pp, gg):
        return fn(yy, pp, gg) if takes_groups else fn(yy, pp)

    pt = float(_eval(y, p, groups))
    if len(uniq) < 5:
        return {"point": pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "n_clusters": int(len(uniq)), "n_singleton_clusters": n_single,
                "ci_setting": "insufficient_clusters"}
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        rows = np.concatenate([cidx[j] for j in pick])
        # 重采样后 cluster 标签必须重编号（同一 cluster 被抽多次时视为多个副本）
        rep = np.concatenate([np.full(len(cidx[j]), i, dtype=np.int64)
                              for i, j in enumerate(pick)])
        v = _eval(y[rows], p[rows], rep)
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return {"point": pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "n_clusters": int(len(uniq)), "n_singleton_clusters": n_single,
                "ci_setting": "degenerate"}
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"point": pt, "ci_low": float(lo), "ci_high": float(hi),
            "n_clusters": int(len(uniq)), "n_singleton_clusters": n_single,
            "ci_setting": "cluster_bootstrap_effective"}


def cluster_paired_bootstrap(fn, y, pa, pb, groups,
                             n_boot: int = N_BOOT_CLUSTER,
                             seed: int = SEED) -> dict:
    """cluster 配对 bootstrap: Δ = fn(y, pb) − fn(y, pa)。

    每次重采样 cluster，保证 a/b 在同一样本上比较（配对）。
    """
    import inspect
    y, pa, pb = (np.asarray(x, float) for x in (y, pa, pb))
    groups = np.asarray(groups)
    uniq, cidx = _cluster_index(groups)
    n_cl = len(uniq)
    try:
        nparams = len(inspect.signature(fn).parameters)
    except (TypeError, ValueError):
        nparams = 2
    takes_groups = nparams >= 3

    def _eval(yy, pp, gg):
        return fn(yy, pp, gg) if takes_groups else fn(yy, pp)

    ma, mb = float(_eval(y, pa, groups)), float(_eval(y, pb, groups))
    d_pt = mb - ma
    if n_cl < 5:
        return {"delta": d_pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "p_two_sided": float("nan"), "n_clusters": n_cl,
                "bootstrap_unit": "cluster", "insufficient_clusters": True}
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        pick = rs.randint(0, n_cl, n_cl)
        rows = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i, dtype=np.int64)
                              for i, j in enumerate(pick)])
        va, vb = _eval(y[rows], pa[rows], rep), _eval(y[rows], pb[rows], rep)
        if va == va and vb == vb:
            d.append(vb - va)
    if len(d) < 20:
        return {"delta": d_pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "p_two_sided": float("nan"), "n_clusters": n_cl,
                "bootstrap_unit": "cluster", "degenerate": True}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    pl = (np.sum(d <= 0) + 1) / (len(d) + 1)
    pr = (np.sum(d >= 0) + 1) / (len(d) + 1)
    dc = d - d.mean()
    pl0 = (np.sum(dc <= -abs(d_pt)) + 1) / (len(d) + 1)
    pr0 = (np.sum(dc >= abs(d_pt)) + 1) / (len(d) + 1)
    return {"delta": float(d_pt), "ci_low": float(lo), "ci_high": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl0, pr0))),
            "p_sign": float(min(1.0, 2 * min(pl, pr))),
            "n_clusters": n_cl, "bootstrap_unit": "cluster",
            "n_distinct_deltas": int(len(np.unique(np.round(d, 10))))}


def row_paired_bootstrap(fn, y, pa, pb, n_boot: int = N_BOOT_ROW,
                         seed: int = SEED) -> dict:
    """legacy row-level 配对 bootstrap（**次级分析**，与 V4 历史可比）。"""
    y, pa, pb = (np.asarray(x, float) for x in (y, pa, pb))
    ma, mb = float(fn(y, pa)), float(fn(y, pb))
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        i = rs.randint(0, len(y), len(y))
        va, vb = fn(y[i], pa[i]), fn(y[i], pb[i])
        if va == va and vb == vb:
            d.append(vb - va)
    if len(d) < 20:
        return {"delta": mb - ma, "ci_low": float("nan"), "ci_high": float("nan"),
                "p_two_sided": float("nan"), "bootstrap_unit": "row"}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    dc = d - d.mean()
    pl0 = (np.sum(dc <= -abs(mb - ma)) + 1) / (len(d) + 1)
    pr0 = (np.sum(dc >= abs(mb - ma)) + 1) / (len(d) + 1)
    return {"delta": float(mb - ma), "ci_low": float(lo), "ci_high": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl0, pr0))),
            "bootstrap_unit": "row"}


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


def decide(q: float, delta: float, metric: str, direction_ok: int,
           n_seeds: int = 3) -> str:
    """V5 §101/§102/§103 决策标签（规则写死在代码里，不人工挑数字）。"""
    if not np.isfinite(q) or not np.isfinite(delta):
        return "RESOURCE_BLOCKED"
    if metric in ("overall_pearson",):
        thr = PRACTICAL_DPEARSON
    elif metric in ("macro_pdb_spearman", "macro_pdb_pearson"):
        thr = PRACTICAL_DMACRO_SPEARMAN
    elif metric == "p_at_10":
        thr = PRACTICAL_DP_AT_10
    elif metric == "hit10_beneficial":
        thr = PRACTICAL_DHIT10
    else:
        thr = PRACTICAL_DPEARSON
    if q < 0.05 and delta >= thr and direction_ok * 2 >= n_seeds:
        return "KEEP"
    if q < 0.05 and 0 < delta < thr:
        return "STATISTICALLY_SIGNIFICANT_BUT_SMALL"
    if q >= 0.05 and delta > 0:
        return "INCONCLUSIVE"
    if delta <= 0:
        return "DROP"
    return "INCONCLUSIVE"


# ============================================================ 实验登记
ALL_EXP_CSV = V5 / "all_experiments.csv"
PAIRED_CSV = V5 / "paired_comparisons.csv"

ALL_EXP_COLS = ["experiment_id", "family", "variant", "seed", "train_n", "val_n",
                "test_n", "train_complex_n", "test_complex_n", "overall_pearson",
                "overall_spearman", "rmse", "mae", "auroc", "macro_pdb_pearson",
                "macro_pdb_spearman", "p_at_5", "p_at_10", "ef_1pct", "ef_5pct",
                "hit10_beneficial", "hit10_ddg_m05", "hit10_ddg_m10",
                "efficiency_at_10", "anti_error", "cycle_error", "runtime_sec",
                "peak_ram_mb", "status", "notes", "git_commit"]

PAIRED_COLS = ["baseline", "candidate", "metric", "delta", "ci_low", "ci_high",
               "p_value", "bh_q", "effect_size", "bootstrap_unit", "n_clusters",
               "practical_threshold", "decision", "family", "seed_pairing"]


def append_rows(csv_path: Path, rows: list[dict], cols: list[str]) -> None:
    """追加（不覆盖已有行）。缺列填 ''。"""
    if not rows:
        return
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df = df[cols]
    if csv_path.exists():
        old = pd.read_csv(csv_path, low_memory=False)
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(csv_path, index=False, encoding="utf-8")


def git_commit() -> str:
    if not (ROOT / ".git").exists():
        return "GIT_NOT_AVAILABLE"
    import subprocess
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                              capture_output=True, text=True,
                              timeout=20).stdout.strip() or "GIT_UNKNOWN"
    except Exception:  # noqa: BLE001
        return "GIT_UNKNOWN"


def peak_ram_mb() -> float:
    """尽力报告峰值 RSS；Windows 下用 ctypes 读 PeakWorkingSetSize。"""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return float(pmc.PeakWorkingSetSize) / (1024 * 1024)
    except Exception:  # noqa: BLE001
        return float("nan")
