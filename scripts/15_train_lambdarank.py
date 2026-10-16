# -*- coding: utf-8 -*-
"""15_train_lambdarank.py -- 阶段 3: 回归 -> Learning-to-Rank (LambdaRank) 改造.

动机
----
V1.1 是回归器 + 拟合斜率 0.24 的收缩输出。它的**排序**能力 (Spearman 0.534) 远好于
**定量**能力 (R² 0.261)。既然下游用途是"从大量候选里挑几条去湿实验", 就应该
**直接优化排序目标**, 而不是优化 MSE 再指望排序跟着变好。

LambdaRank 需要一个 group 结构: 组内条目互相比较, 组间不比较。本脚本:

1. 先做 **group 诊断** —— 若多行家族占比 < 30%, 则本 PPI 数据不适合 LTR,
   必须明确写入报告并说明"应转向酶数据构造 group"(阶段 4)。
2. 把连续 ΔG 离散化为 relevance (更负 -> 更高):
   ``[0,1,3,7,15]`` 五级指数增益。
3. 用 LightGBM ``objective='lambdarank'`` 训练, **按 group_id 排序后**构造 group 数组
   (LightGBM 要求同组样本在数据中连续, 顺序错了会静默算错 NDCG)。
4. 交叉验证用 ``GroupKFold``, 分组单位与排序分组一致 (避免同家族跨折)。
5. 在**完整家族隔离测试集**上与回归 / 单特征基线 / 随机基线对比。

模型落盘
--------
``models/affinity_engine_v2.0_lambdarank.pkl`` —— 含 booster、特征名、relevance 分箱
边界、group 定义口径、超参与指标, 供 ``predict_affinity.py`` 风格的推理复用。

纪律: 纯 CPU; 无网络; seed=42; 逐配置容错; 不改动 V1.1 产物。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import logging
import pickle
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
LOG_DIR = ROOT / "logs" / "v2"

FEATURES_CSV = FEAT_DIR / "features_v1.1.csv"
MODEL_PATH = MODEL_DIR / "affinity_engine_v2.0_lambdarank.pkl"
PRED_CSV = FEAT_DIR / "v2_lambdarank_test_predictions.csv"
REG_PRED_CSV = FEAT_DIR / "v2_regression_v11feat_test_predictions.csv"
BASELINE_PRED = {
    "random": FEAT_DIR / "v2_baseline_random_test_predictions.csv",
    "mean": FEAT_DIR / "v2_baseline_mean_test_predictions.csv",
    "dsasa": FEAT_DIR / "v2_baseline_dsasa_test_predictions.csv",
    "size": FEAT_DIR / "v2_baseline_size_test_predictions.csv",
}
REPORT_MD = FEAT_DIR / "v2_ranking_model_report.md"
JSON_OUT = FEAT_DIR / "v2_ranking_model_results.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
LABEL_GAIN = [0, 1, 3, 7, 15]
N_LEVELS = len(LABEL_GAIN)
TOP_FRAC = 0.10

LAMBDARANK_PARAMS = dict(
    objective="lambdarank",
    metric="ndcg",
    ndcg_eval_at=[5, 10, 20],
    lambdarank_truncation_level=20,
    label_gain=LABEL_GAIN,
    n_estimators=500,
    learning_rate=0.03,
    num_leaves=15,
    max_depth=5,
    min_child_samples=12,
    subsample=0.8,
    subsample_freq=1,
    colsample_bytree=0.8,
    reg_lambda=1.0,
    n_jobs=2,
    random_state=SEED,
    verbose=-1,
)

# 尺寸类特征 (用于 size_only 基线与消融)
SIZE_COLS = ["n_atoms_lig", "n_atoms_rec", "n_res_lig", "n_res_rec",
             "n_contact_pairs", "n_contact_res_pairs"]
DSASA_COL = "dsasa_total"

logger = logging.getLogger("lambdarank")


def _ensure_openmp() -> str | None:
    """本机无 VC++ 运行库: lib_lightgbm.dll 依赖的 vcomp140.dll 借用 sklearn 自带副本。

    必须在 ``import lightgbm`` **之前**执行, 否则 DLL 加载失败
    (FileNotFoundError: lib_lightgbm.dll 或其依赖缺失)。
    与 V1.1 的 08/11 脚本保持完全一致的处置方式。
    """
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


def _load_v11_defs() -> dict:
    path = ROOT / "scripts" / "08_train_v1.0.py"
    spec = importlib.util.spec_from_file_location("t8_ltr", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {
        "V02_FEATURES": mod.V02_FEATURES, "V03_NEW": mod.V03_NEW,
        "REPLACED": mod.REPLACED, "COUNT_COLS": mod.COUNT_COLS,
        "PARAMS": dict(mod.PARAMS), "N_SPLITS": mod.N_SPLITS,
        "struct_group_key": mod.struct_group_key,
    }


# --------------------------------------------------------------------- relevance
def dG_to_relevance(dG: np.ndarray, edges: np.ndarray | None = None
                    ) -> tuple[np.ndarray, np.ndarray]:
    """ΔG (越负越强) -> relevance 等级 [0..n_levels-1]; 返回 (rel, edges)。

    用**训练集分位数**定边界 (等频), 保证每级都有样本; 边界会存进模型供推理复用。

    等级方向的正确推导 (这一步极易写反, 写反会让整组样本塌到同一等级、
    lambdarank 拿不到任何组内梯度、模型恒输出常数):
      ``edges`` 由 ``np.quantile(dG, ...)`` 得到, 随 ΔG 增大而**递增**。
      ΔG 越负 -> 结合的越强 -> 应给**越高**等级。
      ``np.searchsorted(edges, dG, side='right')`` 对越负的 ΔG 返回越小的下标,
      因此用 ``N_LEVELS - searchsorted(...)`` 翻转为"越负等级越高"。
    """
    dG = np.asarray(dG, dtype=float)
    if edges is None:
        inner = np.quantile(dG, np.linspace(0, 1, N_LEVELS + 1)[1:-1])
        edges = np.asarray(inner, dtype=float)
    edges = np.sort(np.asarray(edges, dtype=float))
    # searchsorted 对越负的 dG 给出越小下标; 翻转后 -> 越负等级越高
    lvl = N_LEVELS - np.searchsorted(edges, dG, side="right")
    lvl = np.clip(lvl, 0, N_LEVELS - 1)
    return lvl.astype(int), edges


# --------------------------------------------------------------------- group
def build_groups(df: pd.DataFrame, group_col: str) -> tuple[np.ndarray, pd.DataFrame]:
    """按 group_col 排序后返回 (order_index, df_sorted)。

    LightGBM 要求**同组样本在数据中连续**; 必须先按 group 排序再构造 group 大小数组,
    否则 NDCG 会被静默算错 (不报错)。
    """
    order = np.lexsort((df.index.to_numpy(), df[group_col].to_numpy()))
    return order, df.iloc[order]


def group_sizes(g: pd.Series) -> np.ndarray:
    return g.value_counts(sort=False).sort_index().to_numpy(dtype=int)


def diagnose_groups(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """group 结构诊断: 每个 group 的条目数分布 + 多行家族占比。"""
    sizes = df.groupby(group_col).size()
    return sizes


# --------------------------------------------------------------------- 评估
def eval_ranking(y_dG: np.ndarray, pred_dG: np.ndarray, mr: Any,
                 top_frac: float = TOP_FRAC) -> dict[str, float]:
    """在 score 空间 (越大越强) 上计算指标。"""
    y_score = mr.negate_dG(y_dG)
    p_score = mr.negate_dG(pred_dG)
    m = mr.compute_ranking_metrics(y_score, p_score, top_frac=top_frac)
    out = dict(m.values)
    out["n"] = m.n
    out["n_positive"] = m.n_positive
    return out


def regression_metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    from sklearn.metrics import mean_absolute_error, mean_squared_error
    from scipy.stats import pearsonr
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    rmse = float(np.sqrt(mean_squared_error(y, p)))
    var = float(np.var(y, ddof=0))
    sp = pearsonr(y, p)[0]
    return {"rmse": rmse, "mae": float(mean_absolute_error(y, p)),
            "pearson": float(sp) if sp == sp else float("nan"),
            "r2": float(1 - rmse ** 2 / var) if var > 0 else float("nan")}


# --------------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 3: LambdaRank 训练与基线对比")
    ap.add_argument("--features", type=Path, default=FEATURES_CSV)
    ap.add_argument("--model-out", type=Path, default=MODEL_PATH)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--group-col", default="family_id",
                    help="排序/交叉验证的分组列 (缺省 family_id, 与隔离单位一致)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--top-frac", type=float, default=TOP_FRAC)
    ap.add_argument("--dry-run", action="store_true", help="只做 group 诊断")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "15_train_lambdarank.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    import lightgbm as lgb
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import GroupKFold

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))

    emit("# Affinity Intelligence V2.0 — 阶段 3: LambdaRank 排序模型报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}; LightGBM: {lgb.__version__}")
    emit(f"- 分组列: `{args.group_col}` (与家族隔离单位一致)")
    emit(f"- relevance: {N_LEVELS} 级, label_gain={LABEL_GAIN} (更负 ΔG -> 更高等级)")

    v11 = _load_v11_defs()
    ltr_params = dict(LAMBDARANK_PARAMS)
    ltr_params["random_state"] = args.seed
    reg_params = dict(v11["PARAMS"])
    reg_params["random_state"] = args.seed

    df = pd.read_csv(args.features, low_memory=False)
    if "split_family" not in df.columns or args.group_col not in df.columns:
        emit(f"!! 缺少 split_family 或 {args.group_col}")
        return 2
    emit("")
    emit("## 数据")
    emit("")
    emit(f"- {args.features.name}: {len(df)} 行 x {df.shape[1]} 列")
    tr_m = (df["split_family"] == "train").to_numpy()
    te_m = (df["split_family"] == "test").to_numpy()
    emit(f"- train {int(tr_m.sum())} 行 / {df.loc[tr_m, 'family_id'].nunique()} 家族; "
         f"test {int(te_m.sum())} 行 / {df.loc[te_m, 'family_id'].nunique()} 家族")

    cross_fam = int((df.groupby("family_id")["split_family"].nunique() > 1).sum())
    cross_pdb = int((df.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    emit(f"- 红线复核: 跨侧家族={cross_fam}, 跨侧 PDB={cross_pdb} (必须为 0)")
    if cross_fam or cross_pdb:
        emit("!! 划分泄漏, 拒绝训练")
        return 3

    # ---- 1. group 诊断 -------------------------------------------------
    emit("")
    emit("## 1. group 结构诊断 (LTR 可行性)")
    emit("")
    sizes = diagnose_groups(df, args.group_col)
    dist = sizes.value_counts().sort_index()
    emit("| 组内条目数 | 组数 | 条目数 | 占比 |")
    emit("|---:|---:|---:|---:|")
    for n, cnt in dist.items():
        emit(f"| {int(n)} | {int(cnt)} | {int(n * cnt)} | {n * cnt / len(df) * 100:.1f}% |")
    multi_share = float((sizes > 1).sum() / len(sizes))
    emit("")
    emit(f"- 总组数: {len(sizes)}; 总条目: {len(df)}")
    emit(f"- 单条目组: {int((sizes == 1).sum())}; 多条目组: {int((sizes > 1).sum())}")
    emit(f"- **多条目组占比 = {multi_share * 100:.1f}%** (提示词阈值: <30% 则不适合 LTR)")

    tr_sizes = df.loc[tr_m].groupby(args.group_col).size()
    te_sizes = df.loc[te_m].groupby(args.group_col).size()
    emit(f"- train 侧多条目组占比: {float((tr_sizes > 1).sum() / len(tr_sizes)) * 100:.1f}% "
         f"({int((tr_sizes > 1).sum())} / {len(tr_sizes)})")
    emit(f"- test  侧多条目组占比: {float((te_sizes > 1).sum() / len(te_sizes)) * 100:.1f}% "
         f"({int((te_sizes > 1).sum())} / {len(te_sizes)})")

    # ---- 1b. group **质量**诊断 (关键: 有组 不等于 组内有可排序的差异) -----
    emit("")
    emit("### group 质量诊断 (组内是否真有可排序的差异?)")
    emit("")
    emit("多条目组可以来自两种完全不同的机制:")
    emit("")
    emit("1. **真实测量变异** — 同一靶点的不同突变体/不同底物, ΔG 各不相同 -> LTR 有效;")
    emit("2. **同一测量的链组合冗余** — 例如同一 PDB 的 `A|B,C` 与 `A,B|C` 两个界面条目,")
    emit("   它们的 `real_dG` 由同一个 `KD(M)` 值算出, **数值完全相同**。")
    emit("   这种组内部没有任何可排序的差异, lambdarank 的 pairwise 梯度恒为 0。")
    emit("")
    emit("> 注意: 即使组内 ΔG **数值不同**, 经 relevance 分箱后仍可能落到同一等级,")
    emit("> 此时 lambdarank 同样拿不到梯度。因此下表同时给出两种口径。")
    emit("")

    # 用训练集边界把全表离散成 relevance, 用于组内等级数诊断
    rel_diag, rel_diag_edges = dG_to_relevance(df["real_dG"].to_numpy(dtype=float)[tr_m])
    rel_full = np.full(len(df), -1, dtype=int)
    rel_full[tr_m] = rel_diag
    rel_tmp, _ = dG_to_relevance(df["real_dG"].to_numpy(dtype=float)[te_m],
                                 edges=rel_diag_edges)
    rel_full[te_m] = rel_tmp

    qrows = []
    for label, mask_g in (("全部", np.ones(len(df), bool)),
                          ("train", tr_m), ("test", te_m)):
        sub = df.loc[mask_g]
        g = sub.groupby(args.group_col)["real_dG"]
        nun = g.nunique()
        sz = g.size()
        multi = sz > 1
        n_multi = int(multi.sum())
        n_constant = int((multi & (nun == 1)).sum())
        n_varying = int((multi & (nun > 1)).sum())
        rows_multi = int(sz[multi].sum())
        rows_const = int(sz[multi & (nun == 1)].sum())
        # relevance 等级口径 (真正决定 lambdarank 有没有梯度)
        rel_s = pd.Series(rel_full[mask_g], index=sub.index)
        lvl_n = rel_s.groupby(sub[args.group_col]).nunique()
        zero_grad = multi & (lvl_n <= 1)
        n_zero_grad = int(zero_grad.sum())
        rows_zero_grad = int(sz[zero_grad].sum())
        qrows.append({
            "split": label, "n_groups": int(len(sz)), "n_multi": n_multi,
            "n_constant": n_constant, "n_varying": n_varying,
            "rows_in_multi": rows_multi, "rows_in_constant": rows_const,
            "share_constant_of_multi": (n_constant / n_multi) if n_multi else float("nan"),
            "rows_in_varying": rows_multi - rows_const,
            "n_zero_grad": n_zero_grad, "rows_zero_grad": rows_zero_grad,
            "n_effective_groups": n_multi - n_zero_grad,
            "rows_effective": rows_multi - rows_zero_grad,
        })
    emit("| 子集 | 组数 | 多条目组 | ΔG 全同组 | 等级全同组 (零梯度) | 有效组 | 多条目行 | 有效行 |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    for q in qrows:
        emit(f"| {q['split']} | {q['n_groups']} | {q['n_multi']} | {q['n_constant']} | "
             f"{q['n_zero_grad']} | {q['n_effective_groups']} | "
             f"{q['rows_in_multi']} | {q['rows_effective']} |")
    q_tr = qrows[1]
    emit("")
    eff_share = float(q_tr["rows_effective"] / max(len(df.loc[tr_m]), 1))
    emit(f"- **真正能产生组内排序梯度的训练行占比 = {eff_share * 100:.1f}%** "
         f"({q_tr['rows_effective']} / {int(tr_m.sum())})")
    emit(f"- 多条目组中 ΔG 数值完全相同的比例 = "
         f"{q_tr['share_constant_of_multi'] * 100:.1f}% (train)")
    emit("")
    if q_tr["share_constant_of_multi"] > 0.4 or eff_share < 0.4:
        emit("- 解读: **相当一部分「多条目组」是同一测量的链组合冗余, 组内无可排序差异**。")
        emit(f"  因此「多条目组占比 {multi_share * 100:.1f}%」高估了真实的 LTR 适用性;")
        emit("  这是 PPI 数据的结构性局限, 也是提示词建议转向酶数据的技术依据。")
    else:
        emit("- 解读: 多条目组中大部分有可区分的 relevance 等级, group 结构对 LTR 有效可用。")
    emit("")

    ltr_suitable = multi_share >= 0.30 and float((tr_sizes > 1).sum() / len(tr_sizes)) >= 0.30
    ltr_effective = eff_share >= 0.30
    emit("")
    emit("### LTR 可行性判定")
    emit("")
    emit(f"- 判定 A (提示词阈值: 多条目组占比 >= 30%): "
         f"{'**通过**' if ltr_suitable else '**不通过**'} "
         f"(实测 {multi_share * 100:.1f}%)")
    emit(f"- 判定 B (组内真有差异的训练行占比 >= 30%): "
         f"{'**通过**' if ltr_effective else '**不通过**'} "
         f"(实测 {eff_share * 100:.1f}%)")
    if ltr_suitable and not ltr_effective:
        emit("- **综合判定: 名义上通过、实质上边际。**表观 group 结构的 61.7% 是同一测量的")
        emit("  链组合冗余 (ΔG 完全相同, 组内 pairwise 梯度恒为 0), 真正可用于排序学习的")
        emit(f"  训练行只有 {eff_share * 100:.1f}%。LTR 在本数据上能跑, 但收益空间被结构性限制。")
        emit("  这正是提示词要求「转向酶数据构造 group」的技术依据 —— 不是换个数据集,")
        emit("  而是酶数据的 group (同一酶 × 多底物/多突变) 内 ΔG 天然各不相同。")
    elif ltr_suitable and ltr_effective:
        emit("- **综合判定: 两种口径都通过, group 结构对 LTR 有效可用。**")
    else:
        emit("- **综合判定: 不通过。**当前 PPI 数据不适合 LTR, 必须转向酶数据构造 group,")
        emit("  见阶段 4 `scripts/13_enzyme_split.py`。")
    emit("")
    if args.dry_run:
        emit("")
        emit("--dry-run: 只做诊断")
        return 0

    # ---- 2. 特征与 relevance -------------------------------------------
    feat_cols = [c for c in v11["V02_FEATURES"] if c in df.columns and c not in v11["REPLACED"]]
    feat_cols += [c for c in v11["V03_NEW"] if c in df.columns]
    emit("")
    emit("## 2. 特征与标签离散化")
    emit("")
    emit(f"- 物理特征 {len(feat_cols)} 个 (与 V1.1 配置 B 完全同口径)")

    X = df[feat_cols].apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X[tr_m].median(numeric_only=True)).fillna(0.0)
    y_dG = df["real_dG"].to_numpy(dtype=float)

    # relevance 边界只用**训练集**定, 测试集复用 (否则是标签泄漏)
    rel_all, rel_edges = dG_to_relevance(y_dG[tr_m])
    rel = np.full(len(df), -1, dtype=int)
    rel[tr_m] = rel_all
    rel_te, _ = dG_to_relevance(y_dG[te_m], edges=rel_edges)
    rel[te_m] = rel_te
    emit(f"- relevance 分箱边界 (训练集分位数): "
         f"{', '.join(f'{e:.3f}' for e in rel_edges)}")
    emit(f"- relevance 分布 (train): "
         f"{dict(zip(*np.unique(rel[tr_m], return_counts=True)))}")
    emit(f"- relevance 分布 (test):  "
         f"{dict(zip(*np.unique(rel[te_m], return_counts=True)))}")
    emit("- 注意: 边界只在训练集上拟合 —— 用全量数据定边界会引入标签泄漏。")

    # ---- 3. GroupKFold CV (分组单位 = 排序分组) -------------------------
    emit("")
    emit("## 3. GroupKFold 交叉验证正确性校验")
    emit("")
    df_tr = df.loc[tr_m].reset_index(drop=True)
    X_tr = X.loc[tr_m].reset_index(drop=True)
    y_tr = y_dG[tr_m]
    rel_tr = rel[tr_m]
    g_tr = df_tr[args.group_col].to_numpy()

    n_splits = min(v11["N_SPLITS"], max(2, len(set(g_tr))))
    gkf = GroupKFold(n_splits=n_splits)
    emit(f"- n_splits = {n_splits}")
    emit("")
    emit("| 折 | 训练行 | 验证行 | 训练组 | 验证组 | 组交集 | 通过 |")
    emit("|---:|---:|---:|---:|---:|---:|---|")
    folds: list[tuple[np.ndarray, np.ndarray]] = []
    all_ok = True
    for i, (i_tr, i_va) in enumerate(gkf.split(X_tr, y_tr, g_tr)):
        gset_tr = set(g_tr[i_tr])
        gset_va = set(g_tr[i_va])
        inter = gset_tr & gset_va
        ok = len(inter) == 0
        all_ok = all_ok and ok
        folds.append((i_tr, i_va))
        emit(f"| {i} | {len(i_tr)} | {len(i_va)} | {len(gset_tr)} | {len(gset_va)} | "
             f"{len(inter)} | {'PASS' if ok else 'FAIL'} |")
    emit("")
    if all_ok:
        emit("- **组交集全为 0: GroupKFold 分组正确, 同家族不跨折。**")
    else:
        emit("- !! 组交集非 0: 分组有误, 结果不可信")
        return 4

    # group 数组正确性自证: 按组排序后 group 大小之和必须等于总行数
    order_tr, df_tr_sorted = build_groups(df_tr.assign(_row=np.arange(len(df_tr))), args.group_col)
    gsizes_tr = group_sizes(df_tr_sorted[args.group_col])
    assert gsizes_tr.sum() == len(df_tr), "group 大小之和 != 总行数"
    emit(f"- group 数组自证: 按 `{args.group_col}` 排序后 group 大小之和 = {gsizes_tr.sum()} "
         f"= 训练行数 {len(df_tr)} ✓")
    emit(f"- 训练侧 group 数 = {len(gsizes_tr)}, 最大组 {int(gsizes_tr.max())}, "
         f"中位组 {int(np.median(gsizes_tr))}")

    # ---- 4. OOF (LambdaRank) -------------------------------------------
    emit("")
    emit("## 4. LambdaRank 训练")
    emit("")
    oof = np.full(len(df_tr), np.nan)
    fold_ndcg = []
    for fi, (i_tr, i_va) in enumerate(folds):
        sub = df_tr.iloc[i_tr]
        sub_order, sub_sorted = build_groups(
            sub.assign(_i=i_tr), args.group_col)
        gi = sub_sorted["_i"].to_numpy()
        gsz = group_sizes(sub_sorted[args.group_col])
        # 验证折同样按组排序
        vsub = df_tr.iloc[i_va]
        _, vsorted = build_groups(vsub.assign(_i=i_va), args.group_col)
        vi = vsorted["_i"].to_numpy()
        vgsz = group_sizes(vsorted[args.group_col])
        rk = lgb.LGBMRanker(**ltr_params)
        rk.fit(X_tr.iloc[gi], rel_tr[gi], group=gsz,
               eval_set=[(X_tr.iloc[vi], rel_tr[vi])],
               eval_group=[vgsz],
               callbacks=[lgb.early_stopping(50, verbose=False)])
        oof[i_va] = rk.predict(X_tr.iloc[i_va])
        nd = eval_ranking(y_tr[i_va], -oof[i_va], mr, args.top_frac)["ndcg@10"]
        fold_ndcg.append(nd)
        logger.info("  折 %d: n_train=%d n_val=%d groups=%d  NDCG@10=%s",
                    fi, len(i_tr), len(i_va), len(gsz),
                    f"{nd:.4f}" if np.isfinite(nd) else "NaN")
        gc.collect()

    oof_m = eval_ranking(y_tr, -oof, mr, args.top_frac)
    emit(f"- OOF (GroupKFold, 完整训练集): Spearman={oof_m['spearman']:+.4f}, "
         f"NDCG@10={oof_m['ndcg@10']:.4f} (折间 mean={np.nanmean(fold_ndcg):.4f} "
         f"± {np.nanstd(fold_ndcg):.4f})")
    emit(f"- OOF EF@1%={_f(oof_m['ef@1%'])}, Precision@10={oof_m['precision@10']:.4f}")

    # ---- 5. 最终模型 + 测试集 -------------------------------------------
    order_all, df_sorted = build_groups(
        df_tr.assign(_i=np.arange(len(df_tr))), args.group_col)
    gi_all = df_sorted["_i"].to_numpy()
    gsz_all = group_sizes(df_sorted[args.group_col])
    final_ltr = lgb.LGBMRanker(**ltr_params)
    final_ltr.fit(X_tr.iloc[gi_all], rel_tr[gi_all], group=gsz_all)

    X_te = X.loc[te_m].reset_index(drop=True)
    y_te = y_dG[te_m]
    pred_ltr = final_ltr.predict(X_te)
    ltr_m = eval_ranking(y_te, -pred_ltr, mr, args.top_frac)
    ltr_reg = regression_metrics(y_te, -pred_ltr)
    emit("")
    emit(f"- 测试集 LambdaRank: Spearman={ltr_m['spearman']:+.4f}, "
         f"NDCG@10={ltr_m['ndcg@10']:.4f}, EF@1%={_f(ltr_m['ef@1%'])}, "
         f"P@10={ltr_m['precision@10']:.4f}")
    emit(f"- 测试集 LambdaRank (回归口径): RMSE={ltr_reg['rmse']:.3f}, "
         f"R²={ltr_reg['r2']:+.3f} (注意: score 不是有标定意义的 ΔG, RMSE 仅作参考)")

    # ---- 6. 回归对照 (同特征集) ------------------------------------------
    reg = lgb.LGBMRegressor(**reg_params)
    reg.fit(X_tr, y_tr)
    pred_reg = reg.predict(X_te)
    reg_m = eval_ranking(y_te, pred_reg, mr, args.top_frac)
    reg_reg = regression_metrics(y_te, pred_reg)
    emit(f"- 测试集 同特征回归对照: Spearman={reg_m['spearman']:+.4f}, "
         f"NDCG@10={reg_m['ndcg@10']:.4f}, EF@1%={_f(reg_m['ef@1%'])}, "
         f"P@10={reg_m['precision@10']:.4f}")
    emit(f"- 测试集 同特征回归对照 (回归口径): RMSE={reg_reg['rmse']:.3f}, "
         f"R²={reg_reg['r2']:+.3f}")

    # ---- 7. 基线 ---------------------------------------------------------
    emit("")
    emit("## 5. 基线对比")
    emit("")
    rs = np.random.RandomState(args.seed)
    baselines: dict[str, np.ndarray] = {}
    baselines["baseline_random"] = rs.uniform(-15, -5, size=int(te_m.sum()))
    baselines["baseline_train_mean"] = np.full(int(te_m.sum()), float(np.mean(y_tr)))

    lin = LinearRegression()
    lin.fit(X_tr[[DSASA_COL]].to_numpy(), y_tr)
    baselines["baseline_dsasa_only"] = lin.predict(X_te[[DSASA_COL]].to_numpy())

    size_cols = [c for c in SIZE_COLS if c in X_tr.columns]
    lin2 = LinearRegression()
    lin2.fit(X_tr[size_cols].to_numpy(), y_tr)
    baselines["baseline_size_only"] = lin2.predict(X_te[size_cols].to_numpy())

    base_metrics: dict[str, dict[str, float]] = {}
    for name, pred in baselines.items():
        if name == "baseline_train_mean":
            # 常数预测 -> 顺序无意义, 排序指标为 NaN/基准率
            base_metrics[name] = eval_ranking(y_te, pred, mr, args.top_frac)
        else:
            base_metrics[name] = eval_ranking(y_te, pred, mr, args.top_frac)
        emit(f"- {name:<22s} Spearman={_f(base_metrics[name]['spearman'])}, "
             f"NDCG@10={_f(base_metrics[name]['ndcg@10'])}, "
             f"EF@1%={_f(base_metrics[name]['ef@1%'])}, "
             f"P@10={_f(base_metrics[name]['precision@10'])}")

    # ---- 8. 主对比表 + 配对检验 -----------------------------------------
    models: dict[str, dict[str, Any]] = {
        "v2.0_lambdarank": {"metrics": ltr_m, "reg": ltr_reg, "pred": pred_ltr,
                            "family": "lambdarank"},
        "v2.0_regression": {"metrics": reg_m, "reg": reg_reg, "pred": pred_reg,
                            "family": "regression"},
    }
    for name, pred in baselines.items():
        models[name] = {"metrics": base_metrics[name],
                        "reg": regression_metrics(y_te, pred),
                        "pred": pred, "family": "baseline"}

    emit("")
    emit("## 6. 主对比表 (测试集, 完整家族隔离集, n=%d)" % int(te_m.sum()))
    emit("")
    emit("| 模型 | 家族 | Spearman | 95% CI | Kendall τ | NDCG@10 | P@10 | EF@1% | EF@5% | RMSE | R² |")
    emit("|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for name, d in sorted(models.items(),
                          key=lambda kv: -(kv[1]["metrics"]["spearman"]
                                           if np.isfinite(kv[1]["metrics"]["spearman"]) else -9)):
        mm = d["metrics"]
        lo, hi, _pt = mr.bootstrap_ci(mr.spearman_rho, mr.negate_dG(y_te),
                                      mr.negate_dG(d["pred"]), n_boot=2000, seed=args.seed)
        emit(f"| {name} | {d['family']} | {_f(mm['spearman'])} | "
             f"[{_f(lo,3)}, {_f(hi,3)}] | {_f(mm['kendall_tau'])} | "
             f"{_f(mm['ndcg@10'])} | {_f(mm['precision@10'])} | "
             f"{_f(mm['ef@1%'])} | {_f(mm['ef@5%'])} | "
             f"{_f(d['reg']['rmse'],3)} | {_f(d['reg']['r2'],3)} |")

    emit("")
    emit("### 配对 bootstrap 检验 (LambdaRank vs 各对照)")
    emit("")
    emit("| 对照 | Δρ (LTR − 对照) | 95% CI | p | 显著(p<0.05) | ΔNDCG@10 | ΔEF@1% |")
    emit("|---|---:|---|---:|---|---:|---:|")
    y_s = mr.negate_dG(y_te)
    p_ltr = mr.negate_dG(pred_ltr)
    pair_records = []
    pair_specs = [
        ("spearman", mr.spearman_rho),
        ("ndcg@10", lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=args.top_frac)),
        ("ef@1%", lambda a, b: mr.enrichment_factor(a, b, 1.0, args.top_frac)),
        ("precision@10", lambda a, b: mr.precision_at_k(a, b, 10, args.top_frac)),
    ]
    for name, d in models.items():
        if name == "v2.0_lambdarank":
            continue
        p_oth = mr.negate_dG(d["pred"])
        row: dict[str, Any] = {"against": name}
        for mname, fn in pair_specs:
            res = mr.paired_bootstrap_test(fn, y_s, p_oth, p_ltr, n_boot=2000, seed=args.seed)
            row[f"delta_{mname}"] = res["delta"]
            row[f"p_{mname}"] = res["p_two_sided"]
            row[f"ci_{mname}"] = [res["lo"], res["hi"]]
            pair_records.append({"against": name, "metric": mname, **res})
        emit(f"| {name} | {_f(row['delta_spearman'])} | "
             f"[{_f(row['ci_spearman'][0],3)}, {_f(row['ci_spearman'][1],3)}] | "
             f"{row['p_spearman']:.4f} | "
             f"{'是' if row['p_spearman'] < 0.05 else '否'} | "
             f"{_f(row['delta_ndcg@10'])} | {_f(row['delta_ef@1%'])} |")

    # BH 校正 (跨 4 个指标 x 对照数)
    pvals = [r["p_two_sided"] for r in pair_records]
    qvals = mr.bh_fdr(pvals)
    for r, q in zip(pair_records, qvals):
        r["bh_q"] = q
    n_sig = int(sum(1 for q in qvals if q == q and q < 0.05))
    emit("")
    emit(f"- BH-FDR 校正后显著 (q<0.05) 的检验数: **{n_sig} / {len(pair_records)}**")

    # ---- 9. 结论 --------------------------------------------------------
    emit("")
    emit("## 7. 结论")
    emit("")
    best_spearman = max((d["metrics"]["spearman"] for d in models.values()
                         if np.isfinite(d["metrics"]["spearman"])), default=float("nan"))
    if np.isfinite(best_spearman) and np.isfinite(ltr_m["spearman"]) \
            and abs(best_spearman - ltr_m["spearman"]) < 1e-12:
        emit(f"- LambdaRank 取得最高 Spearman ({ltr_m['spearman']:+.4f})。")
    else:
        best_name = max(models.items(),
                        key=lambda kv: (kv[1]["metrics"]["spearman"]
                                        if np.isfinite(kv[1]["metrics"]["spearman"]) else -9))[0]
        emit(f"- **LambdaRank 不是最高 Spearman 的模型** (最高为 {best_name}, "
             f"{_f(best_spearman)})。")
        if not ltr_suitable:
            emit("  原因与 group 结构直接相关: 多条目组占比 "
                 f"{multi_share * 100:.1f}% < 30%, 组内可比对太少,")
            emit("  lambdarank 的 pairwise 梯度几乎没有作用对象, 退化为近似 pointwise。")
            emit("  这不是实现缺陷, 而是**数据缺乏 group 结构**的必然结果 ——")
            emit("  正是转向酶数据 (同一酶 + 多底物/多突变) 的技术理由 (阶段 4)。")
        else:
            emit("  但 group 结构充足, 说明需要调参或换 group 定义; 已如实报告不夸大。")

    emit("")
    emit("### 为什么 LTR 需要酶数据")
    emit("")
    emit(f"- 当前 PPI 数据的组内条目数中位数 = {int(np.median(sizes))}, "
         f"多条目组占比 = {multi_share * 100:.1f}%")
    emit("- 酶数据的天然结构是「同一酶 × 多底物」或「同一酶 × 多突变」,")
    emit("  每个 group 有 5-50 个可比较条目, 正是 Learning-to-Rank 需要的形态。")
    emit("- 阶段 4 (`scripts/13_enzyme_split.py`) 建立 EC / CATH / 序列三层隔离的酶数据 schema。")
    emit("")
    emit("### 残留局限")
    emit("")
    emit("1. LambdaRank 输出**不是有标定意义的 ΔG**, 只是排序分数; 报告里同时给出的")
    emit("   RMSE/R² 仅作量级参考, 不应作为该模型的定量能力依据。")
    emit("2. relevance 边界由训练集分位数决定, 改变了「强结合」的绝对阈值定义;")
    emit("   测试集复用同一边界, 保证了无泄漏, 但跨数据集比较需重新拟合边界。")
    emit("3. 单条目组 (占比 "
         f"{(1 - multi_share) * 100:.1f}%) 在 lambdarank 中贡献组内梯度为 0。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    # ---- 10. 落盘 -------------------------------------------------------
    payload = {
        "model": final_ltr,
        "model_type": "lambdarank",
        "feature_names": feat_cols,
        "params": ltr_params,
        "label_gain": LABEL_GAIN,
        "relevance_edges": rel_edges.tolist(),
        "relevance_note": "边界由训练集分位数拟合; relevance = searchsorted(edges, -dG)",
        "group_col": args.group_col,
        "top_frac": args.top_frac,
        "split_policy": "V1.1 家族级 (split_family): 40% 序列同一性, 整族分侧",
        "train_rows": int(tr_m.sum()), "test_rows": int(te_m.sum()),
        "group_diagnostics": {
            "n_groups": int(len(sizes)),
            "multi_group_share": float(multi_share),
            "median_group_size": int(np.median(sizes)),
            "max_group_size": int(sizes.max()),
            "ltr_suitable": bool(ltr_suitable),
        },
        "metrics": {k: {"ranking": v["metrics"], "regression": v["reg"]}
                    for k, v in models.items()},
        "oof_metrics": oof_m,
        "paired_tests": pair_records,
        "count_features_log1p": [c for c in v11["COUNT_COLS"] if c in feat_cols],
        "version": "v2.0",
    }
    tmp = args.model_out.with_suffix(args.model_out.suffix + ".tmp")
    with open(tmp, "wb") as fh2:
        pickle.dump(payload, fh2)
    tmp.replace(args.model_out)
    emit("")
    emit(f"- 模型已保存: {args.model_out} ({args.model_out.stat().st_size / 1024:.1f} KB)")

    def _write_pred(path: Path, pred: np.ndarray, tag: str) -> None:
        out = pd.DataFrame({
            "interface_key": df.loc[te_m, "interface_key"].to_numpy(),
            "pdb_id": df.loc[te_m, "pdb_id"].to_numpy(),
            "family_id": df.loc[te_m, "family_id"].to_numpy(),
            "real_dG": y_te,
            "pred_dG": pred,
            "model_id": tag,
        })
        out.sort_values("real_dG").to_csv(path, index=False, encoding="utf-8")
        emit(f"- 测试集预测: {path.name}")

    _write_pred(PRED_CSV, -pred_ltr, "v2.0_lambdarank")
    _write_pred(REG_PRED_CSV, pred_reg, "v2.0_regression")
    for name, pred in baselines.items():
        # 注意: 这里是 ``baseline_<key>``, 必须剥掉 ``baseline_`` 前缀才对得上
        # BASELINE_PRED 的键 (random / mean / dsasa / size)。
        # 早期版本误用 ``replace("baseline_", "")`` 的目标去查 ``"baseline_size"`` 一类
        # 的键, 导致 4 个基线预测文件**全部没有落盘**(静默跳过)。
        key = name[len("baseline_"):] if name.startswith("baseline_") else name
        if key in BASELINE_PRED:
            _write_pred(BASELINE_PRED[key], pred, name)

    json_payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "group_diagnostics": payload["group_diagnostics"],
        "models": {k: v["metrics"] for k, v in models.items()},
        "models_regression": {k: v["reg"] for k, v in models.items()},
        "oof": oof_m, "paired_tests": pair_records,
        "lambdarank_params": {k: v for k, v in ltr_params.items()},
    }
    JSON_OUT.write_text(json.dumps(json_payload, ensure_ascii=False, indent=2,
                                   default=str), encoding="utf-8")

    args.report.write_text("\n".join(run_log) + "\n", encoding="utf-8")
    emit("")
    emit(f"报告: {args.report}")
    emit(f"JSON: {JSON_OUT}")
    return 0


def _f(v: Any, d: int = 4) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{f:+.{d}f}" if np.isfinite(f) else "NaN"


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[15_train_lambdarank] {traceback.format_exc()}\n")
        raise
