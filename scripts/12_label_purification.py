# -*- coding: utf-8 -*-
"""12_label_purification.py -- 阶段 1: 标签净化与测量异质性消融.

科学问题
--------
V1.1 在 40% 家族隔离测试集上 R²=0.261、拟合斜率仅 0.24。学长建议的第一优先级是
**零算力的"标签净化实验"**: 先确认这个偏低的 R² 有多少来自「测量异质性」
(Kd/Ki/IC50/Km 混装、温度混杂、方法混杂), 而不是模型容量不足。

数据现实 (必须先说清楚, 见 AUTO_ASSUMPTIONS.md A3/A4)
----------------------------------------------------
* 本项目原始表 ``PPB-Affinity.xlsx`` **只有 KD(M) 一列亲和力数值**, 不存在
  Ki / IC50 / Km 列 (19 列已全量枚举)。因此提示词要求的
  ``Kd_Ki_only`` / ``Kd_only`` / ``Km_aux`` 三个变体在本数据上**退化为 no-op**
  (与 ``full`` 逐行相同), 本脚本会如实打印 ``identical_to_full`` 标记,
  不伪造"净化增益"。
* 全部 4762 条记录统一按 ``ΔG = 0.592 × ln(Kd[M])`` 换算 (298 K 口径)。
* 真正可做的净化维度是 **测量方法 (assay_method)** 与 **实验温度 (temperature_K)**。
  这两个字段在 ``01_clean_data.py`` 的聚合里被 ``"first"`` 丢弃了,
  已由 ``12a_recover_assay_meta.py`` 从源表重建 (``data/v2_assay_meta.csv``)。

做的事
------
1. 合并 ``features/features_v1.1.csv`` (4750 行, 含 split_family) 与恢复的测量元数据。
2. 构造 8 个数据集变体 + 2 个特征改动对照 + 1 个标签口径敏感性对照。
3. 在**完全相同的家族隔离划分**下训练同一个 LightGBM 回归器 (与 V1.1 同超参),
   报告 Spearman / RMSE / MAE / R² / bootstrap 95% CI / 测试集行数 / 家族数。
4. 用**配对 bootstrap** 检验每个变体相对 ``full`` 的 Spearman 差异, 并做 BH-FDR 校正。
5. 明确回答: 净化后 Spearman/R² 是否提升? 提升多少? 若没提升, 为什么?

关键方法论约束
--------------
* 变体只在**训练集**上过滤, 测试集永远保持完整的家族隔离全集 —— 这样所有变体的
  指标都在同一个测试集上可比, 配对 bootstrap 才成立。若变体只保留 train 侧,
  报告里标注 ``test_full=True``; 若变体同时过滤 test 侧, 标注 ``test_full=False``
  且该行指标**不可与其他行直接比大小** (样本集不同)。
* 净化后样本量骤降时必须报告**统计功效上限** (以 n 和可分辨 Δρ 说明), 否则
  "没提升" 会被误读成 "净化无用"。

纪律: 纯 CPU; 无网络; 固定 seed=42; 逐变体 try/except 不中断全流程;
      失败写入 logs/v2/failures.log; 不改动任何 V1.1 产物。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import logging
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT_DIR = ROOT / "features"
DATA_DIR = ROOT / "data"
LOG_DIR = ROOT / "logs" / "v2"

FEATURES_CSV = FEAT_DIR / "features_v1.1.csv"
ASSAY_META_CSV = DATA_DIR / "v2_assay_meta.csv"
RESULT_CSV = FEAT_DIR / "v2_label_purification_results.csv"
REPORT_MD = FEAT_DIR / "v2_label_purification_report.md"
JSON_OUT = FEAT_DIR / "v2_label_purification_results.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
N_BOOT = 2000
N_PERM = 15
PERM_EVAL_ROUNDS = 30          # 仅对 full / purified 跑更稳的置换基线
RT_KCAL = 1.98720425864083e-3  # kcal/(mol*K)
REF_TEMP_K = 298.15            # 0.592 kcal/mol 对应的参考温度

# 统计功效地板: 低于此规模的变体只作参考, 不足以支撑结论
MIN_TEST_ROWS = 60
MIN_TEST_FAMILIES = 30
MIN_MAJORITY_SHARE = 0.80      # 单一方法的纯度下限

TEMP_BIN = (293.0, 303.0)      # 提示词建议的单一温度区间

logger = logging.getLogger("label_purification")


def _fmt(v: Any, d: int = 4) -> str:
    """安全格式化数值: None / NaN 一律显示为 ``—``, 绝不抛错。"""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{f:+.{d}f}" if np.isfinite(f) else "—"


# --------------------------------------------------------------------- 复用 V1.1 特征集定义
def _load_v11_defs() -> dict:
    """从 08_train_v1.0.py 载入特征集/超参, 保证与 V1.1 **单一实现来源**, 避免口径漂移。"""
    path = ROOT / "scripts" / "08_train_v1.0.py"
    spec = importlib.util.spec_from_file_location("t8_purif", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {
        "V02_FEATURES": mod.V02_FEATURES,
        "V03_NEW": mod.V03_NEW,
        "REPLACED": mod.REPLACED,
        "COUNT_COLS": mod.COUNT_COLS,
        "PARAMS": dict(mod.PARAMS),
        "N_SPLITS": mod.N_SPLITS,
        "struct_group_key": mod.struct_group_key,
    }


# --------------------------------------------------------------------- 指标
def spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr
    if len(a) < 3:
        return float("nan")
    c = spearmanr(a, b).correlation
    return float(c) if c == c else float("nan")


def regression_metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    from sklearn.metrics import mean_absolute_error, mean_squared_error
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    if len(y) < 3:
        return {k: float("nan") for k in ("spearman", "pearson", "rmse", "mae", "r2", "n")}
    rmse = float(np.sqrt(mean_squared_error(y, p)))
    var = float(np.var(y, ddof=0))
    from scipy.stats import pearsonr
    pe = pearsonr(y, p)[0]
    return {
        "spearman": spearman(y, p),
        "pearson": float(pe) if pe == pe else float("nan"),
        "rmse": rmse,
        "mae": float(mean_absolute_error(y, p)),
        "r2": float(1 - rmse ** 2 / var) if var > 0 else float("nan"),
        "n": int(len(y)),
    }


def bootstrap_spearman_ci(y: np.ndarray, p: np.ndarray, n: int = N_BOOT,
                          seed: int = SEED) -> tuple[float, float]:
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    if len(y) < 5:
        return float("nan"), float("nan")
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n):
        idx = rs.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) < 3:
            continue
        c = spearman(y[idx], p[idx])
        if c == c:
            vals.append(c)
    if len(vals) < 20:
        return float("nan"), float("nan")
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def paired_bootstrap_spearman(y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                              n: int = N_BOOT, seed: int = SEED) -> dict[str, float]:
    """配对 bootstrap: 同一批样本上重采样, 统计 ρ(B) - ρ(A) 的分布。"""
    ok = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
    y, pa, pb = y[ok], pa[ok], pb[ok]
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y))}
    rs = np.random.RandomState(seed)
    diffs = []
    for _ in range(n):
        idx = rs.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) < 3:
            continue
        ca, cb = spearman(y[idx], pa[idx]), spearman(y[idx], pb[idx])
        if ca == ca and cb == cb:
            diffs.append(cb - ca)
    if len(diffs) < 20:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y))}
    d = np.asarray(diffs)
    lo, hi = np.percentile(d, [2.5, 97.5])
    # 双侧 p 值: 2 * min(P(d<=0), P(d>=0)), 加 1 平滑避免 p=0
    p_left = (np.sum(d <= 0) + 1) / (len(d) + 1)
    p_right = (np.sum(d >= 0) + 1) / (len(d) + 1)
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(p_left, p_right))), "n": int(len(y))}


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR 校正 (自实现, 不引入 statsmodels 依赖)。"""
    p = np.asarray([1.0 if (x != x) else float(x) for x in pvals], dtype=float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p, kind="mergesort")
    q = np.empty(m, dtype=float)
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank  # 0-based position in sorted order
        val = p[order[i]] * m / (i + 1)
        prev = min(prev, val)
        q[order[i]] = prev
    return [float(min(1.0, v)) for v in q]


# --------------------------------------------------------------------- 变体定义
def build_variants(df: pd.DataFrame) -> list[dict]:
    """构造数据集变体。每个变体 = 训练侧过滤条件 (测试侧默认不过滤)。"""
    method = df["method_majority"].fillna("").to_numpy()
    tmed = df["temperature_median_k"].to_numpy(dtype=float)
    tused = df["temperature_k_used"].to_numpy(dtype=float)
    t_known = df["temperature_n"].fillna(0).to_numpy() > 0
    lo, hi = TEMP_BIN

    def m(x: str) -> np.ndarray:
        return method == x

    variants: list[dict] = [
        {
            "name": "full",
            "desc": "全部数据 (V1.1 口径)",
            "train_mask": np.ones(len(df), dtype=bool),
            "label_col": "real_dG",
        },
        {
            "name": "Kd_only",
            "desc": "只保留 Kd 标签 (源表仅 KD(M), 故与 full 逐行相同)",
            "train_mask": np.ones(len(df), dtype=bool),
            "label_col": "real_dG",
            "expect_identical": "full",
        },
        {
            "name": "Kd_Ki_only",
            "desc": "只保留 Kd/Ki (源表无 Ki 列, 故与 full 逐行相同)",
            "train_mask": np.ones(len(df), dtype=bool),
            "label_col": "real_dG",
            "expect_identical": "full",
        },
        {
            "name": "single_method_SPR",
            "desc": "只保留主流单一方法 SPR (train 侧过滤)",
            "train_mask": m("SPR"),
            "label_col": "real_dG",
        },
        {
            "name": "single_method_ITC",
            "desc": "只保留主流单一方法 ITC (train 侧过滤)",
            "train_mask": m("ITC"),
            "label_col": "real_dG",
        },
        {
            "name": "single_temp_bin",
            "desc": f"只保留单温度区间 {lo:.0f}-{hi:.0f} K (有实验温度者)",
            "train_mask": t_known & (tmed >= lo) & (tmed <= hi),
            "label_col": "real_dG",
        },
        {
            "name": "single_temp_bin_wide",
            "desc": "只保留 285-310 K 宽温度区间 (有实验温度者)",
            "train_mask": t_known & (tmed >= 285) & (tmed <= 310),
            "label_col": "real_dG",
        },
        {
            "name": "purified",
            "desc": f"SPR + 单温度区间 {lo:.0f}-{hi:.0f} K",
            "train_mask": m("SPR") & t_known & (tmed >= lo) & (tmed <= hi),
            "label_col": "real_dG",
        },
        {
            "name": "purified_wide",
            "desc": "单一方法 (SPR) + 有实验温度 (不限区间)",
            "train_mask": m("SPR") & t_known,
            "label_col": "real_dG",
        },
        {
            "name": "method_consistent",
            "desc": "剔除同界面方法冲突条目 (measurement heterogeneity 清洗)",
            "train_mask": df["method_conflict"].fillna(0).to_numpy() == 0,
            "label_col": "real_dG",
        },
        {
            "name": "Kd_at_T",
            "desc": "标签口径敏感性: ΔG = RT·ln(Kd) 用**实验温度**重算 (非主口径)",
            "train_mask": np.ones(len(df), dtype=bool),
            "label_col": "dG_at_T",
        },
        {
            "name": "paired_method_only",
            "desc": "只保留有方法标注的条目 (方法作特征可用)",
            "train_mask": method != "",
            "label_col": "real_dG",
        },
    ]
    return variants


def feature_columns(v11: dict, df: pd.DataFrame, add_meta: bool) -> list[str]:
    base = [c for c in v11["V02_FEATURES"] if c in df.columns and c not in v11["REPLACED"]]
    base += [c for c in v11["V03_NEW"] if c in df.columns]
    if add_meta:
        for c in ("temperature_k_used", "method_is_known", "temp_is_known",
                  "temperature_spread_k"):
            if c in df.columns:
                base.append(c)
    return base


# --------------------------------------------------------------------- 训练
def train_eval(
    df: pd.DataFrame,
    y: np.ndarray,
    feat_cols: list[str],
    variant_mask: np.ndarray,
    v11: dict,
    seed: int = SEED,
    perm_rounds: int = N_PERM,
) -> dict:
    """在给定**变体过滤掩码**上训练 LightGBM 回归, 在完整测试集上评估。

    ``tr_mask`` 实际为 ``variant_mask & (split_family == 'train')`` —— 变体过滤只作用于
    训练侧, 测试集恒定, 这样所有变体在同一测试集上可比, 配对 bootstrap 才成立。
    """
    import lightgbm as lgb
    from sklearn.model_selection import GroupKFold

    te_mask = df["split_family"].to_numpy() == "test"
    tr_mask = variant_mask & (~te_mask)
    X = df[feat_cols].apply(pd.to_numeric, errors="coerce")
    med = X[tr_mask].median(numeric_only=True)
    X = X.fillna(med).fillna(0.0)

    X_tr, y_tr = X[tr_mask], y[tr_mask]
    X_te, y_te = X[te_mask], y[te_mask]
    groups_tr = v11["struct_group_key"](df[tr_mask])

    n_splits = min(v11["N_SPLITS"], max(2, len(set(groups_tr))))
    gkf = GroupKFold(n_splits=n_splits)
    oof = np.full(len(X_tr), np.nan)
    fold_sp = []
    for i_tr, i_va in gkf.split(X_tr, y_tr, groups_tr):
        mdl = lgb.LGBMRegressor(**v11["PARAMS"])
        mdl.fit(
            X_tr.iloc[i_tr], y_tr[i_tr],
            eval_set=[(X_tr.iloc[i_va], y_tr[i_va])],
            callbacks=[lgb.early_stopping(50, verbose=False)],
        )
        oof[i_va] = mdl.predict(X_tr.iloc[i_va])
        sp = spearman(y_tr[i_va], oof[i_va])
        if sp == sp:
            fold_sp.append(sp)
        gc.collect()

    final = lgb.LGBMRegressor(**v11["PARAMS"])
    final.fit(X_tr, y_tr)
    pred = final.predict(X_te)

    tm = regression_metrics(y_te, pred)
    lo, hi = bootstrap_spearman_ci(y_te, pred)
    oof_m = regression_metrics(y_tr, oof)

    # 置换基线: 打乱训练标签, 量化"无信号"参考线。
    # 性能取舍: 用**单次按组划分** (GroupShuffleSplit) 而不是 5 折 CV ——
    # 置换基线的唯一用途是给出"完全没有信号时 Spearman 长什么样"的参照,
    # 单折已经足够, 而 5 折会让这一段的拟合次数变成 5 倍 (实测 12 个变体总耗时
    # 因此从 ~22 min 降到 ~6 min), 对结论没有任何影响。
    from sklearn.model_selection import GroupShuffleSplit
    rng = np.random.RandomState(seed)
    perm = []
    if perm_rounds > 0 and len(set(groups_tr)) >= 5:
        gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed)
        try:
            i_perm_tr, i_perm_va = next(gss.split(X_tr, y_tr, groups_tr))
        except StopIteration:
            i_perm_tr, i_perm_va = np.arange(len(X_tr)), np.arange(len(X_tr))
        for _ in range(perm_rounds):
            yp = rng.permutation(y_tr)
            mm = lgb.LGBMRegressor(**v11["PARAMS"])
            mm.fit(X_tr.iloc[i_perm_tr], yp[i_perm_tr])
            pr = mm.predict(X_tr.iloc[i_perm_va])
            c = spearman(yp[i_perm_va], pr)
            if c == c:
                perm.append(c)

    return {
        "n_train": int(tr_mask.sum()),
        "n_test": int(te_mask.sum()),
        "n_features": int(X_tr.shape[1]),
        "n_families_train": int(df.loc[tr_mask, "family_id"].nunique()),
        "n_families_test": int(df.loc[te_mask, "family_id"].nunique()),
        "n_structs_train": int(df.loc[tr_mask, "pdb_id"].nunique()),
        "cv_spearman_mean": float(np.nanmean(fold_sp)) if fold_sp else float("nan"),
        "cv_spearman_std": float(np.nanstd(fold_sp)) if fold_sp else float("nan"),
        "oof_spearman": oof_m["spearman"],
        "test": tm,
        "test_ci": [lo, hi],
        "perm_mean": float(np.nanmean(perm)) if perm else float("nan"),
        "perm_std": float(np.nanstd(perm)) if perm else float("nan"),
        "y_test": y_te,
        "pred_test": pred,
    }


# --------------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 1: 标签净化与测量异质性消融")
    ap.add_argument("--features", type=Path, default=FEATURES_CSV)
    ap.add_argument("--assay-meta", type=Path, default=ASSAY_META_CSV)
    ap.add_argument("--out-csv", type=Path, default=RESULT_CSV)
    ap.add_argument("--out-md", type=Path, default=REPORT_MD)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--perm-rounds", type=int, default=N_PERM)
    ap.add_argument("--perm-only", action="store_true",
                    help="只保留 full / purified 主对照, 快速冒烟")
    ap.add_argument("--dry-run", action="store_true", help="只打印变体规模, 不训练")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    (ROOT / "figures").mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "12_label_purification.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))
        # 显式 flush: 本脚本耗时长, 若 stdout 被重定向到文件会走块缓冲,
        # 导致外部看不到任何进度 (实测: 跑 18 分钟日志里一行变体结果都没有)。
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# Affinity Intelligence V2.0 — 阶段 1: 标签净化与测量异质性消融")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed} (LightGBM random_state 同为 {args.seed})")
    emit(f"- bootstrap 次数: {args.n_boot}; 置换轮数: {args.perm_rounds}")

    if not args.features.exists():
        emit(f"!! 特征表不存在: {args.features}")
        return 2
    if not args.assay_meta.exists():
        emit(f"!! 测量元数据不存在: {args.assay_meta} (先运行 12a_recover_assay_meta.py)")
        return 2

    v11 = _load_v11_defs()
    v11["PARAMS"]["random_state"] = args.seed

    df = pd.read_csv(args.features, low_memory=False)
    meta = pd.read_csv(args.assay_meta)
    emit("")
    emit("## 数据")
    emit("")
    emit(f"- 特征表: {args.features.name} — {len(df)} 行 x {df.shape[1]} 列")
    emit(f"- 测量元数据: {args.assay_meta.name} — {len(meta)} 行")
    df = df.merge(meta, on="interface_key", how="left", validate="one_to_one")
    emit(f"- 合并后: {len(df)} 行")

    n_meta_missing = int(df["method_majority"].isna().sum())
    if n_meta_missing:
        emit(f"- !! {n_meta_missing} 行未匹配到测量元数据 (按「无方法信息」处理)")
    df["method_majority"] = df["method_majority"].fillna("")
    df["method_is_known"] = (df["method_majority"] != "").astype(int)
    df["temp_is_known"] = (df["temperature_n"].fillna(0) > 0).astype(int)

    # 标签口径: 主口径 real_dG (298 K, 与 V1.1 一致); 敏感性口径用实验温度重算
    temp_used = df["temperature_k_used"].to_numpy(dtype=float)
    df["dG_at_T"] = RT_KCAL * temp_used * np.log(df["kd_m"].to_numpy(dtype=float)
                                                 if "kd_m" in df.columns
                                                 else np.exp(df["real_dG"].to_numpy() / 0.592))
    df["dG_at_T"] = np.where(np.isfinite(df["dG_at_T"]), df["dG_at_T"], df["real_dG"])

    if "family_id" not in df.columns or "split_family" not in df.columns:
        emit("!! 缺少 family_id / split_family, 无法执行家族隔离评估")
        return 2

    # ---- 红线复核 (与 11_train_v1.1.py 同口径) ---------------------------
    cross_fam = int((df.groupby("family_id")["split_family"].nunique() > 1).sum())
    cross_pdb = int((df.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    emit(f"- 红线复核: 跨侧家族={cross_fam}, 跨侧 PDB={cross_pdb} (必须为 0)")
    if cross_fam or cross_pdb:
        emit("!! 划分泄漏, 拒绝继续")
        return 3

    te_mask_full = df["split_family"].to_numpy() == "test"
    y_full = df["real_dG"].to_numpy(dtype=float)
    emit(f"- 测试集 (固定, 所有变体共用): {int(te_mask_full.sum())} 行 / "
         f"{df.loc[te_mask_full, 'family_id'].nunique()} 家族 / "
         f"{df.loc[te_mask_full, 'pdb_id'].nunique()} 结构")
    emit(f"- 训练集: {int((~te_mask_full).sum())} 行 / "
         f"{df.loc[~te_mask_full, 'family_id'].nunique()} 家族")

    # ---- 标签异质性画像 ------------------------------------------------
    emit("")
    emit("## 标签与测量异质性画像")
    emit("")
    emit("| 维度 | 取值 | 行数 | 占比 | train | test |")
    emit("|---|---|---:|---:|---:|---:|")
    for label, col in (("测量方法", "method_majority"), ("来源数据集", "source_dataset")):
        vc = df[col].replace("", "(缺失)").value_counts()
        for k, v in vc.head(8).items():
            ntr = int(((df[col].replace("", "(缺失)") == k) & (~te_mask_full)).sum())
            nte = int(((df[col].replace("", "(缺失)") == k) & te_mask_full).sum())
            emit(f"| {label} | {k} | {v} | {v / len(df) * 100:.1f}% | {ntr} | {nte} |")
    tb = pd.cut(df["temperature_median_k"], bins=[0, 290, 295, 299, 303, 310, 400],
                labels=["<290", "290-295", "295-299", "299-303", "303-310", ">310"])
    for k, v in tb.value_counts().sort_index().items():
        emit(f"| 温度区间 | {k} | {v} | {v / len(df) * 100:.1f}% | "
             f"{int(((tb == k) & (~te_mask_full)).sum())} | "
             f"{int(((tb == k) & te_mask_full).sum())} |")
    emit("")
    emit(f"- 标签口径: ΔG = 0.592·ln(Kd[M]) 固定 {REF_TEMP_K:.2f} K (V1.1 主口径)")
    emit(f"- ΔG 范围: {y_full.min():.2f} ~ {y_full.max():.2f} kcal/mol "
         f"(mean {y_full.mean():.2f}, std {y_full.std():.2f})")
    dgt = df["dG_at_T"].to_numpy(dtype=float)
    emit(f"- 实验温度重算口径 ΔG: {dgt.min():.2f} ~ {dgt.max():.2f} "
         f"(mean {dgt.mean():.2f}); 与主口径逐行 Spearman="
         f"{spearman(y_full, dgt):+.4f}")
    max_shift = float(np.max(np.abs(dgt - y_full)))
    emit(f"- 两口径最大逐行差异: {max_shift:.3f} kcal/mol")

    # ---- 变体规模预检 ------------------------------------------------
    variants = build_variants(df)
    if args.perm_only:
        keep_names = {"full", "purified", "single_method_SPR"}
        variants = [v for v in variants if v["name"] in keep_names]

    emit("")
    emit("## 变体规模预检 (训练侧过滤, 测试集恒定)")
    emit("")
    emit("| 变体 | 说明 | train 行 | train 家族 | test 行 | test 家族 | 功效 |")
    emit("|---|---|---:|---:|---:|---:|---|")
    feasible: list[dict] = []
    for v in variants:
        tm = v["train_mask"] & (~te_mask_full)
        nte = int(te_mask_full.sum())
        nte_fam = int(df.loc[te_mask_full, "family_id"].nunique())
        ok = nte >= MIN_TEST_ROWS and nte_fam >= MIN_TEST_FAMILIES and int(tm.sum()) >= 150
        v["feasible"] = bool(ok)
        emit(f"| {v['name']} | {v['desc']} | {int(tm.sum())} | "
             f"{df.loc[tm, 'family_id'].nunique()} | {nte} | {nte_fam} | "
             f"{'OK' if ok else '规模不足(参考)'} |")
    emit("")
    emit(f"- 功效地板: 测试侧 >= {MIN_TEST_ROWS} 行 且 >= {MIN_TEST_FAMILIES} 家族, "
         f"训练侧 >= 150 行。低于地板的变体仍会跑, 但结论不作数。")

    if args.dry_run:
        emit("")
        emit("--dry-run: 不做训练")
        return 0

    # ---- 逐变体训练 ---------------------------------------------------
    feat_base = feature_columns(v11, df, add_meta=False)
    feat_withmeta = feature_columns(v11, df, add_meta=True)

    emit("")
    emit("## 训练与评估 (同一家族隔离划分, 同一 LightGBM 超参)")
    emit("")
    emit(f"- 物理特征 {len(feat_base)} 个; 加入测量元数据特征后 {len(feat_withmeta)} 个")
    emit(f"- LightGBM 超参: {json.dumps({k: v for k, v in v11['PARAMS'].items() if k != 'verbose'}, ensure_ascii=False)}")
    emit("")

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}
    y_ref: np.ndarray | None = None
    base_pred: np.ndarray | None = None

    jobs: list[tuple[str, str, np.ndarray, list[str], str]] = []
    for v in variants:
        jobs.append((v["name"], v["desc"], v["train_mask"],
                     feat_base, v.get("label_col", "real_dG")))
    # 附: full 变体 + 测量元数据特征 (检验"把方法/温度当特征"是否有效)
    jobs.append(("full_plus_meta", "full + 测量元数据特征 (温度/方法可见性)",
                 np.ones(len(df), dtype=bool), feat_withmeta, "real_dG"))
    jobs.append(("purified_plus_meta", "purified + 测量元数据特征",
                 [v for v in variants if v["name"] == "purified"][0]["train_mask"],
                 feat_withmeta, "real_dG"))

    for name, desc, mask, cols, label_col in jobs:
        try:
            y = df[label_col].to_numpy(dtype=float)
            # 置换基线只对主对照跑更多轮 (置换是 O(轮数 x 5 折) 拟合, 最耗时)
            pr = PERM_EVAL_ROUNDS if name in ("full", "purified") else args.perm_rounds
            res = train_eval(df, y, cols, mask, v11, seed=args.seed,
                             perm_rounds=pr)
            tm = res["test"]
            rows.append({
                "variant": name, "description": desc, "label_col": label_col,
                "n_features": res["n_features"],
                "n_train": res["n_train"], "n_test": res["n_test"],
                "n_families_train": res["n_families_train"],
                "n_families_test": res["n_families_test"],
                "n_structs_train": res["n_structs_train"],
                "test_spearman": tm["spearman"], "test_ci_low": res["test_ci"][0],
                "test_ci_high": res["test_ci"][1],
                "test_rmse": tm["rmse"], "test_mae": tm["mae"], "test_r2": tm["r2"],
                "cv_spearman_mean": res["cv_spearman_mean"],
                "cv_spearman_std": res["cv_spearman_std"],
                "oof_spearman": res["oof_spearman"],
                "perm_mean": res["perm_mean"], "perm_std": res["perm_std"],
                "test_full": True,
                "feasible": bool(
                    res["n_test"] >= MIN_TEST_ROWS
                    and res["n_families_test"] >= MIN_TEST_FAMILIES
                    and res["n_train"] >= 150
                ),
            })
            preds[name] = res["pred_test"]
            if name == "full":
                base_pred = res["pred_test"]
                y_ref = res["y_test"]
            emit(f"  完成 {name:<22s} Spearman={tm['spearman']:+.3f} "
                 f"[{res['test_ci'][0]:+.3f},{res['test_ci'][1]:+.3f}] "
                 f"RMSE={tm['rmse']:.3f} R2={tm['r2']:+.3f} (train {res['n_train']})")
        except Exception as exc:  # noqa: BLE001
            logger.error("变体 %s 失败: %s", name, exc)
            with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[12_label_purification] {name}: {exc}\n{traceback.format_exc()}\n")
        gc.collect()

    if not rows:
        emit("!! 全部变体失败")
        return 4

    res_df = pd.DataFrame(rows)

    # ---- 配对 bootstrap vs full ---------------------------------------
    emit("")
    emit("## 配对 bootstrap 检验 (vs full, 同一测试集) 与 BH-FDR 校正")
    emit("")
    if base_pred is None or y_ref is None:
        emit("!! 缺少 full 基线预测, 跳过配对检验")
    else:
        # 全部结果**按 variant 名字做键**, 绝不依赖列表位置。
        # 早期版本用 ``zip(res_df.variant, pairs)`` 位置配对, 而 pairs 的构造顺序
        # (跳过 "full") 与 res_df 的行序不同 -> 出现「Δ 点估计落在自己的 CI 之外」
        # 这种自相矛盾的表 (实测 14 行里 9 行违反), 且 BH q 与显著性标记全错。
        pair_by_variant: dict[str, dict] = {}
        pval_by_variant: dict[str, float] = {}
        for r in res_df.itertuples():
            if r.variant == "full":
                continue
            pb = preds.get(r.variant)
            if pb is None:
                continue
            d = paired_bootstrap_spearman(y_ref, base_pred, pb, n=args.n_boot, seed=args.seed)
            pair_by_variant[r.variant] = d
            pval_by_variant[r.variant] = d["p_two_sided"]

        # BH-FDR 只作用于真正算出来的 p 值, 并按 variant 名字回填
        names_with_p = [k for k, v in pval_by_variant.items() if v == v]
        qvals = bh_fdr([pval_by_variant[k] for k in names_with_p])
        q_map = dict(zip(names_with_p, qvals))

        full_sp_ref = float(res_df.loc[res_df.variant == "full", "test_spearman"].iloc[0])
        res_df["delta_spearman_vs_full"] = res_df["test_spearman"] - full_sp_ref
        res_df["paired_p"] = res_df["variant"].map(
            {k: v["p_two_sided"] for k, v in pair_by_variant.items()})
        res_df["paired_ci_low"] = res_df["variant"].map(
            {k: v["lo"] for k, v in pair_by_variant.items()})
        res_df["paired_ci_high"] = res_df["variant"].map(
            {k: v["hi"] for k, v in pair_by_variant.items()})
        res_df["paired_delta_boot"] = res_df["variant"].map(
            {k: v["delta"] for k, v in pair_by_variant.items()})
        res_df["bh_q"] = res_df["variant"].map(q_map)

        emit("| 变体 | Δρ vs full | 配对 Δ 均值 | 95% CI (bootstrap) | 配对 p | BH q | 显著(q<0.05) |")
        emit("|---|---:|---:|---|---:|---:|---|")
        for r in res_df.to_dict("records"):
            if r["variant"] == "full":
                continue
            pv = r.get("paired_p")
            qv = r.get("bh_q")
            clo, chi = r.get("paired_ci_low"), r.get("paired_ci_high")
            dlt = r.get("delta_spearman_vs_full")
            dboot = r.get("paired_delta_boot")
            sig = "是" if (qv is not None and qv == qv and qv < 0.05) else "否"
            emit(f"| {r['variant']} | {_fmt(dlt)} | {_fmt(dboot)} | "
                 f"[{_fmt(clo)}, {_fmt(chi)}] | {_fmt(pv)} | {_fmt(qv)} | {sig} |")

        # 自检: 点估计必须落在自己的 bootstrap CI 内 (否则配对逻辑有错)
        bad = []
        for r in res_df.to_dict("records"):
            if r["variant"] == "full":
                continue
            dlt, clo, chi = (r.get("delta_spearman_vs_full"),
                             r.get("paired_ci_low"), r.get("paired_ci_high"))
            if None in (dlt, clo, chi) or not (clo == clo and chi == chi):
                continue
            if not (clo - 1e-9 <= dlt <= chi + 1e-9):
                bad.append(r["variant"])
        emit("")
        if bad:
            emit(f"- !! **自检失败**: 以下变体的 Δ 点估计落在自身 CI 之外, 配对逻辑有误: {bad}")
        else:
            emit("- ✅ 自检通过: 所有变体的 Δ 点估计都落在其 bootstrap 95% CI 内。")

        # 逐样本预测相关: 变体是否只是把 full 的预测重排
        emit("")
        emit("逐样本预测 Spearman (变体 vs full, 测试集):")
        for name in res_df["variant"]:
            pb = preds.get(name)
            if pb is not None and name != "full":
                emit(f"  {name:<22s} {spearman(base_pred, pb):+.4f}")

    # ---- 主结果表 -----------------------------------------------------
    emit("")
    emit("## 主结果表 (测试集完整, 同一家族隔离划分)")
    emit("")
    emit("| 变体 | 特征数 | train 行 | 家族(train/test) | Spearman | 95% CI | RMSE | MAE | R² | 功效 |")
    emit("|---|---:|---:|---|---:|---|---:|---:|---:|---|")
    for r in res_df.sort_values("test_spearman", ascending=False).itertuples():
        emit(f"| {r.variant} | {r.n_features} | {r.n_train} | "
             f"{r.n_families_train}/{r.n_families_test} | {r.test_spearman:+.3f} | "
             f"[{r.test_ci_low:+.3f},{r.test_ci_high:+.3f}] | {r.test_rmse:.3f} | "
             f"{r.test_mae:.3f} | {r.test_r2:+.3f} | {'OK' if r.feasible else '不足'} |")

    # ---- 退化变体核对 -------------------------------------------------
    emit("")
    emit("## 退化变体核对 (证明 Kd/Ki/Km 维度在本数据上不可分)")
    emit("")
    full_sp = float(res_df.loc[res_df.variant == "full", "test_spearman"].iloc[0])
    for nm in ("Kd_only", "Kd_Ki_only"):
        sub = res_df[res_df.variant == nm]
        if len(sub):
            same = abs(float(sub.test_spearman.iloc[0]) - full_sp) < 1e-9
            emit(f"- {nm}: Spearman={float(sub.test_spearman.iloc[0]):+.6f} "
                 f"vs full {full_sp:+.6f} -> "
                 f"{'逐行相同 (no-op 已证实)' if same else '存在差异(需排查)'}")
    emit("- 结论: 源表 `PPB-Affinity.xlsx` 只有 `KD(M)` 一列亲和力数值, 无 Ki/IC50/Km 列;")
    emit("  因此提示词设想的 `Kd_Ki_only` / `Km_aux` 变体在本数据上**不可能产生差异**。")
    emit("  这是数据现实, 不是实验失败。若要真正做 Kd/Ki/Km 标签异质性研究,")
    emit("  必须接入 SABIO-RK (Km/Ki/kcat) 或 BindingDB (Ki/IC50) —— 见阶段 4。")

    # ---- 直接回答 -----------------------------------------------------
    emit("")
    emit("## 结论")
    emit("")
    best = res_df[res_df.feasible].sort_values("test_spearman", ascending=False)
    if len(best):
        b = best.iloc[0]
        emit(f"- 可行变体中最优: **{b['variant']}** — Spearman={b['test_spearman']:+.3f} "
             f"(full={full_sp:+.3f}, Δ={b['test_spearman'] - full_sp:+.3f})")
    pur = res_df[res_df.variant == "purified"]
    if len(pur):
        p = pur.iloc[0]
        emit(f"- 强净化 `purified` (SPR+单温区): Spearman={p['test_spearman']:+.3f} "
             f"(Δ vs full {p['test_spearman'] - full_sp:+.3f}), "
             f"训练行 {int(p['n_train'])}, RMSE={p['test_rmse']:.3f}, R²={p['test_r2']:+.3f}")
        emit(f"- `purified` 与 full 的配对 p={p['paired_p']:.4f}, BH q={p['bh_q']:.4f}")
    emit("")

    # 天花板分析: 同一 target 内部多突变体之间的可预测性
    emit("### 测量异质性对性能天花板的影响")
    emit("")
    within = []
    for fam, grp in df[te_mask_full].groupby("family_id"):
        if len(grp) >= 3 and grp["real_dG"].std(ddof=0) > 0:
            within.append(len(grp))
    emit(f"- 测试集中同家族多条目家族数 (n>=3): {len(within)}, 覆盖 {sum(within)} 行")
    emit(f"- 测试集单条目家族数: {int((df[te_mask_full].groupby('family_id').size() == 1).sum())}")
    emit("- 解读: 家族隔离测试集里绝大多数家族只有 1 条的样本, 预测是**跨靶外推**;")
    emit("  在该设定下 R² 天然被压低, 与标签净化程度只有间接关系。")
    emit("")
    emit("### 诚实回答: 净化后 Spearman/R² 是否提升?")
    emit("")
    deltas = []
    for nm in ("single_method_SPR", "single_temp_bin", "single_temp_bin_wide",
               "purified", "purified_wide", "method_consistent", "full_plus_meta",
               "purified_plus_meta"):
        sub = res_df[res_df.variant == nm]
        if len(sub):
            d = float(sub.test_spearman.iloc[0]) - full_sp
            deltas.append((nm, d, bool(sub.feasible.iloc[0])))
    for nm, d, feas in deltas:
        emit(f"- {nm}: ΔSpearman = {d:+.4f} ({'功效 OK' if feas else '样本规模不足, 不可作结论'})")
    emit("")
    ok_deltas = [d for _, d, f in deltas if f and abs(d) > 1e-9]
    if not ok_deltas:
        emit("- **净化未带来可用的粒度提升**: 所有净化变体的 ΔSpearman 都在 0 附近或样本不足。")
    else:
        emit(f"- 净化带来最大 ΔSpearman = {max(ok_deltas, key=abs):+.4f} (仅统计功效达标变体)。")
    emit("- 无论结论方向, 都以上表实测数字为准, 不做任何外推。")
    emit("")
    emit("### 残留局限")
    emit("")
    emit("1. 测量元数据 (方法/温度) 覆盖率低: 方法 30.8%, 温度 14.5% (源表恢复后);")
    emit("   这是 PPB-Affinity 原始表的固有缺失, 无法通过清洗补回。")
    emit("2. 强净化变体 (single_method_ITC / purified) 训练行数偏少, 统计功效受限;")
    emit("   已用 95% CI 与配对检验量化不确定性, 未把不显著读成无效。")
    emit("3. `Kd_at_T` 只是标签口径敏感性对照, 不替换 V1.1 主口径 (保持可比)。")
    emit("4. 所有变体共用同一测试集, 因此可逐行配对比较; 但测试集本身是跨靶外推集。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    # ---- 落盘 ---------------------------------------------------------
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out_csv.with_suffix(args.out_csv.suffix + ".tmp")
    res_df.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out_csv)

    payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "seed": args.seed,
        "n_boot": args.n_boot,
        "n_perm": args.perm_rounds,
        "features_csv": str(args.features),
        "assay_meta_csv": str(args.assay_meta),
        "label_definitions": {
            "real_dG": "0.592 * ln(Kd[M])  固定 298.15 K 口径 (V1.1 主口径)",
            "dG_at_T": "R*T*ln(Kd[M]), R=1.98720425864083e-3 kcal/mol/K, T=实验温度(缺失 298.15)",
        },
        "test_set_fixed": True,
        "min_test_rows": MIN_TEST_ROWS,
        "min_test_families": MIN_TEST_FAMILIES,
        "results": rows,
        "bh_fdr_note": "BH 校正作用于所有与 full 的配对检验 p 值",
    }
    with open(JSON_OUT, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    with open(args.out_md, "w", encoding="utf-8") as fh:
        fh.write("\n".join(run_log) + "\n")

    emit("")
    emit(f"结果表: {args.out_csv}")
    emit(f"报告:   {args.out_md}")
    emit(f"JSON:   {JSON_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
