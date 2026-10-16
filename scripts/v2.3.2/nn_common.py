# -*- coding: utf-8 -*-
"""scripts/v2.3.2/nn_common.py -- V2.3.2 任务 3 的共享 NN 训练/评估机器.

为什么单独成模块
----------------
任务书要求 `nn_mlp.py` 与 `nn_two_tower.py` **两个**脚本。两者除了模型结构以外
完全同构 (数据、PCA、训练循环、早停、四个种子、指标、配对检验)。
把共享部分放这里, 两个脚本各自只定义**自己的模型** —— 避免复制粘贴导致的漂移。

严格约束
--------
* **CPU only**: 全程不出现 `.cuda()` / `torch.cuda`; 设备固定为 ``cpu``。
* 早停准则: **验证集 Spearman**, 连续 20 epoch 不提升则停 (任务书要求)。
* 种子: 42 / 43 / 44 / 45 共四次, 分别训练与评估。
* 指标口径与 LightGBM 侧**完全一致** (同一 `metrics_ranking.py`)。
"""

from __future__ import annotations

import gc
import importlib.util
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT_DIR = ROOT / "features"
MODEL_DIR = ROOT / "models"
LOG_DIR = ROOT / "logs" / "v2.3.2"

V22_FEATURES = FEAT_DIR / "v2.2_proteinea_features.csv"
ESM2_CACHE = FEAT_DIR / "v2_esm2_cache.npz"
PREDS_NPZ = FEAT_DIR / "v2.3.2_test_predictions.npz"    # 任务 2 的产物 (含 LightGBM 预测)
FAIL_LOG = LOG_DIR / "failures.log"

SEEDS = (42, 43, 44, 45)
PCA_DIMS = (32, 64, 128)
TOP_FRAC = 0.10
MAX_EPOCH = 300
BATCH = 64
LR = 1e-3
WD = 1e-4
DROPOUT = 0.3
PATIENCE = 20
EXPECTED_SPLITS = {"train": 6482, "val": 965, "test": 752}

# 物理特征列 (与 V2.2 一致, 从 08 载入以保证单一实现来源)
PHYS_FEATS: list[str] = []
COUNT_COLS: list[str] = []


def load_feature_names() -> tuple[list[str], list[str]]:
    global PHYS_FEATS, COUNT_COLS
    spec = importlib.util.spec_from_file_location("t8_v232nn",
                                                  ROOT / "scripts" / "08_train_v1.0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    PHYS_FEATS = ([c for c in mod.V02_FEATURES if c not in mod.REPLACED]
                  + list(mod.V03_NEW))
    COUNT_COLS = list(mod.COUNT_COLS)
    return PHYS_FEATS, COUNT_COLS


def get_metrics_module():
    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    return mr


# ---------------------------------------------------------------- 指标 / 检验
def build_metric_fns(mr: Any, top_frac: float) -> dict[str, Callable]:
    return {
        "spearman": mr.spearman_rho,
        "kendall_tau": mr.kendall_tau,
        "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=top_frac),
        "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, top_frac),
        "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, top_frac),
        "ef@5%": lambda a, b: mr.enrichment_factor(a, b, 5.0, top_frac),
    }


def bootstrap_ci(fn: Callable, y: np.ndarray, p: np.ndarray, n_boot: int,
                 seed: int) -> tuple[float, float, float]:
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
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), point


def paired_bootstrap(fn: Callable, y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                     n_boot: int, seed: int) -> dict[str, float]:
    """配对 bootstrap: Δ = metric(B) - metric(A)。与 LightGBM 侧实现同源。"""
    ok = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
    y, pa, pb = y[ok], pa[ok], pb[ok]
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y)),
                "n_distinct_deltas": 0}
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


# ---------------------------------------------------------------- 数据装配
def load_data(esm2_cache: Path = ESM2_CACHE,
              v22_features: Path = V22_FEATURES) -> dict[str, Any]:
    """载入 V2.2 特征表 + ESM-2 缓存, 做**与 pdb_id join 相同的**装配。

    为与任务 2 的 LightGBM 口径严格一致, 这里复用同一套规则:
    * join 键 = `pdb_id`; 缓存侧一对多取**平均**;
    * 只用 `row_structure_usable == 1` 且标签有效的行;
    * 物理特征用 V2.2 模型的训练集中位数兜底 + COUNT_COLS log1p。
    """
    load_feature_names()
    d = pd.read_csv(v22_features, low_memory=False)
    z = np.load(esm2_cache, allow_pickle=True)
    pdb_ids = z["pdb_ids"].astype(str)
    lig_raw, lig_l2 = z["lig_raw"].astype(np.float32), z["lig_l2"].astype(np.float32)
    rec_raw, rec_l2 = z["rec_raw"].astype(np.float32), z["rec_l2"].astype(np.float32)

    usable = (d["row_structure_usable"] == 1) & d["real_dG"].notna()
    d = d[usable].reset_index(drop=True)
    split_counts = d["split"].value_counts().to_dict()
    for sp, exp in EXPECTED_SPLITS.items():
        if int(split_counts.get(sp, 0)) != exp:
            raise RuntimeError(f"划分与 V2.2 不一致: {sp}={split_counts.get(sp)} != {exp}")

    # 缓存按 pdb_id 取平均
    emb = pd.DataFrame({"pdb_id": pdb_ids,
                        "lig_raw": list(lig_raw), "lig_l2": list(lig_l2),
                        "rec_raw": list(rec_raw), "rec_l2": list(rec_l2)})
    agg = emb.groupby("pdb_id", sort=False).agg(
        lig_raw=("lig_raw", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
        lig_l2=("lig_l2", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
        rec_raw=("rec_raw", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
        rec_l2=("rec_l2", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
    ).reset_index()
    m = d.merge(agg, on="pdb_id", how="left", validate="many_to_one")
    hit = m["lig_raw"].notna().to_numpy()
    if hit.mean() < 0.95:
        raise RuntimeError(f"join 命中率过低: {hit.mean():.4f}")
    m = m[hit].reset_index(drop=True)

    # 物理特征
    fcols = [c for c in PHYS_FEATS if c in m.columns]
    Xp = m[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in COUNT_COLS:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    tr_m = (m["split"] == "train").to_numpy()
    med = Xp[tr_m].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)[fcols]

    E = {
        "lig_raw": np.stack(m["lig_raw"].to_numpy()),
        "lig_l2": np.stack(m["lig_l2"].to_numpy()),
        "rec_raw": np.stack(m["rec_raw"].to_numpy()),
        "rec_l2": np.stack(m["rec_l2"].to_numpy()),
    }
    return {
        "df": m, "Xp": Xp, "E": E, "fcols": fcols,
        "y": m["real_dG"].to_numpy(float),
        "tr_m": tr_m,
        "va_m": (m["split"] == "val").to_numpy(),
        "te_m": (m["split"] == "test").to_numpy(),
        "split_counts": split_counts,
    }


def pca_blocks(data: dict[str, Any], dims=PCA_DIMS, seed: int = 42,
               include_raw640: bool = True) -> dict[str, np.ndarray]:
    """对 ligand/receptor/both 做 PCA, **只在训练集 fit**。返回 {'{variant}|{src}|{d}': Z}。

    ``include_raw640=True`` 时额外提供 ``'{variant}|{src}|raw640'`` = **不做 PCA 的原始嵌入**。
    这是任务书双塔规格 `Linear(640,256)` 字面要求的输入 (640 = 原始 ESM-2 维度)。
    注意: 原始 640 维在 n≈6500 训练行上**过参数化风险高** (远多于样本),
    因此同时跑 PCA 版与 raw640 版, 让数据说话, 而不是替用户选。
    """
    from sklearn.decomposition import PCA
    tr = data["tr_m"]
    out: dict[str, np.ndarray] = {}
    for variant in ("raw", "l2"):
        L = data["E"][f"lig_{variant}"]
        R = data["E"][f"rec_{variant}"]
        for src, mat in (("ligand", L), ("receptor", R), ("both", np.hstack([L, R]))):
            for dim in dims:
                nmax = min(dim, mat.shape[1], int(tr.sum()) - 1)
                p = PCA(n_components=nmax, random_state=seed)
                p.fit(mat[tr])
                out[f"{variant}|{src}|{dim}"] = p.transform(mat).astype(np.float32)
            if include_raw640:
                out[f"{variant}|{src}|raw640"] = np.ascontiguousarray(
                    mat.astype(np.float32))
    return out


# ---------------------------------------------------------------- NN 机器
def train_nn(build_model: Callable[[dict[str, int]], "Any"],
             inputs: dict[str, np.ndarray],
             n_in: dict[str, int],
             y_tr: np.ndarray, y_va: np.ndarray,
             X_tr_parts: dict[str, np.ndarray], X_va_parts: dict[str, np.ndarray],
             X_te_parts: dict[str, np.ndarray],
             seed: int) -> dict[str, Any]:
    """训练一个 NN 并返回验证/测试预测。**纯 CPU**。

    ``build_model(n_in)`` 返回一个 ``nn.Module``, 其 ``forward(batch_dict)``
    接收 ``{input_name: tensor}``, 输出 ``(B, 1)``。
    """
    import torch
    import torch.nn as nn

    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(max(1, min(4, (torch.get_num_threads() or 2))))
    device = torch.device("cpu")          # 硬约束: CPU only, 禁止 cuda

    model = build_model(n_in).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WD)
    lossf = nn.MSELoss()

    def to_tensor(parts: dict[str, np.ndarray]) -> dict[str, "torch.Tensor"]:
        return {k: torch.from_numpy(np.asarray(v, dtype=np.float32)) for k, v in parts.items()}

    tr_t = to_tensor(X_tr_parts)
    va_t = to_tensor(X_va_parts)
    te_t = to_tensor(X_te_parts)
    ytr_t = torch.from_numpy(y_tr.astype(np.float32)).view(-1, 1)

    n_tr = len(y_tr)
    best_rho, best_state, best_epoch, bad = -9.0, None, -1, 0
    mr = get_metrics_module()
    hist: list[dict[str, float]] = []
    t0 = time.time()

    for epoch in range(MAX_EPOCH):
        model.train()
        perm = np.random.RandomState(seed * 1000 + epoch).permutation(n_tr)
        for s in range(0, n_tr, BATCH):
            idx = perm[s:s + BATCH]
            batch = {k: v[idx] for k, v in tr_t.items()}
            opt.zero_grad()
            out = model(batch)
            loss = lossf(out, ytr_t[idx])
            loss.backward()
            opt.step()
        # 早停用验证集 Spearman
        model.eval()
        with torch.no_grad():
            pv = model(va_t).cpu().numpy().ravel()
        rho = mr.spearman_rho(-y_va, -pv)
        hist.append({"epoch": epoch, "val_spearman": float(rho) if rho == rho else float("nan")})
        if rho == rho and rho > best_rho + 1e-6:
            best_rho, best_epoch, bad = float(rho), epoch, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        pv = model(va_t).cpu().numpy().ravel()
        pt = model(te_t).cpu().numpy().ravel()
    return {"val_pred": pv, "test_pred": pt, "best_val_spearman": best_rho,
            "best_epoch": best_epoch, "epochs_ran": len(hist),
            "train_seconds": time.time() - t0, "history": hist}
