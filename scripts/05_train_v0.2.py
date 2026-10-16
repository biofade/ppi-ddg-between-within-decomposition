# -*- coding: utf-8 -*-
"""
05_train_v0.2.py -- Affinity Intelligence V0.2 / 阶段三：重新训练与对比评估

输入 : features/features_mvp.csv  (V0.2 特征表: 归一化特征 + 标签派生量 + 两套划分)
输出 : models/affinity_engine_v0.2.pkl
       features/v0.2_vs_v0.1.png
       features/v0.2_feature_importance.csv
       features/v0.2_report.txt
       features/v0.2_test_predictions.csv

对比矩阵 (2 套划分 x 3 个特征配置):
  划分:  V0.1 = PDB 级结构身份 (含同源泄漏) | V0.2 = 40% 同一性家族级 (已清除泄漏)
  配置:  A 仅原始物理特征 | B 原始 + 体积/规模归一化特征 | C = B + 标签派生量 (泄漏对照)

  关键读法:
    * V0.1 x A  = 复现 V0.1 报告的旧数字 (同样是 PDB 级划分)。
    * V0.2 x B  = V0.2 主模型, **唯一可对外的真实数字**。
    * 两份差异同时来自"划分变严"与"特征变化", 报告中分别拆开:
        - 划分影响 = (V0.2 x A) - (V0.1 x A)
        - 特征影响 = (V0.2 x B) - (V0.2 x A)
    * C 配置用 real_dG 派生量入模, 用于量化"目标泄漏"的代价, 不可采信为模型能力。
纪律: n_jobs 强制 2; 纯 CPU; 无任何深度学习库。
"""

import gc
import os
import pickle
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


def _setup_cjk_font():
    """matplotlib 默认字体无 CJK 字形, 中文会渲染成豆腐块; 显式挑一个系统字体。"""
    prefer = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC",
              "WenQuanYi Zen Hei", "SimSun", "DengXian"]
    try:
        available = {f.name for f in font_manager.fontManager.ttflist}
        for name in prefer:
            if name in available:
                matplotlib.rcParams["font.sans-serif"] = [name] + \
                    matplotlib.rcParams.get("font.sans-serif", [])
                matplotlib.rcParams["axes.unicode_minus"] = False
                return name
    except Exception:                                           # noqa: BLE001
        pass
    return None


_CJK_FONT = _setup_cjk_font()


def _ensure_openmp():
    """本机无 VC++ 运行库, 借用 scikit-learn wheel 自带的 vcomp140.dll。"""
    import ctypes
    import glob
    try:
        import sklearn
        hits = glob.glob(os.path.join(os.path.dirname(sklearn.__file__),
                                      ".libs", "vcomp140.dll"))
        if hits:
            ctypes.CDLL(hits[0])
            return hits[0]
    except Exception:                                           # noqa: BLE001
        pass
    return None


_VCOMP = _ensure_openmp()

import lightgbm as lgb
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from scipy.stats import pearsonr, spearmanr

# ---------------------------------------------------------------- 路径 / 常量
ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
FEAT_DIR = os.path.join(ROOT, "features")
MODEL_DIR = os.path.join(ROOT, "models")
IN_CSV = os.path.join(FEAT_DIR, "features_mvp.csv")
MODEL_PATH = os.path.join(MODEL_DIR, "affinity_engine_v0.2.pkl")
FIG_PATH = os.path.join(FEAT_DIR, "v0.2_vs_v0.1.png")
IMP_PATH = os.path.join(FEAT_DIR, "v0.2_feature_importance.csv")
REPORT_PATH = os.path.join(FEAT_DIR, "v0.2_report.txt")
PRED_PATH = os.path.join(FEAT_DIR, "v0.2_test_predictions.csv")

N_JOBS = 2
N_SPLITS = 5
SEED = 42
TARGET = "real_dG"

META_COLS = ["pdb_id", "interface_key", "ligand_chains", "receptor_chains",
             "split", "split_v0.1", "family_id", "subgroup", "source_dataset",
             "feature_version"]

COUNT_COLS = ["n_contact_pairs", "n_contact_res_pairs", "n_iface_res_lig",
              "n_iface_res_rec", "n_hbond", "n_salt_bridge",
              "n_hydrophobic_contact", "n_chainpair_max", "n_chainpair_min",
              "n_chainpair_contacts", "n_atoms_lig", "n_atoms_rec",
              "n_res_lig", "n_res_rec", "sasa_complex", "sasa_lig_mono",
              "sasa_rec_mono", "sasa_lig_in_complex", "sasa_rec_in_complex",
              "dsasa_lig", "dsasa_rec", "dsasa_total"]

# V0.1 的原始物理特征集 (配置 A 的可比基线)
BASE_FEATURES = [
    "n_missing_chain", "chain_map_ok", "n_contact_pairs", "n_atoms_lig", "n_atoms_rec",
    "n_res_lig", "n_res_rec", "n_contact_res_pairs", "n_iface_res_lig",
    "n_iface_res_rec", "iface_res_frac_lig", "iface_res_frac_rec",
    "min_interatomic_dist", "n_hbond", "n_salt_bridge", "n_hydrophobic_contact",
    "kd_hydro_lig", "kd_hydro_rec", "kd_hydro_all", "frac_charged_iface",
    "frac_polar_iface", "frac_hydrophobic_iface", "frac_helix_iface",
    "frac_sheet_iface", "n_chainpair_max", "n_chainpair_min",
    "n_chainpair_contacts", "contact_density", "sasa_complex", "sasa_lig_mono",
    "sasa_rec_mono", "sasa_lig_in_complex", "sasa_rec_in_complex",
    "dsasa_lig", "dsasa_rec", "dsasa_total", "dsasa_frac_lig", "dsasa_frac_rec",
]

PARAMS = dict(
    objective="regression", n_estimators=500, learning_rate=0.03,
    num_leaves=15, max_depth=5, min_child_samples=12,
    subsample=0.8, subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
    n_jobs=N_JOBS, random_state=SEED, verbose=-1,
)


# ---------------------------------------------------------------- 工具
def leaky_cols(df):
    return [c for c in df.columns if c.endswith("__leaky")]


def normalized_cols(df):
    """新增的归一化特征 = 既非基础特征、也非元数据/标签/泄漏列的数字列。

    同时排除全空列 (如 PRODIGY 不可用时的 predicted_dG / prodigy_ic) ——
    否则它们会被中位数填充成常数列, 虚增特征数且无信息量。
    """
    base = set(BASE_FEATURES) | set(META_COLS) | {TARGET} | set(leaky_cols(df))
    out = []
    for c in df.columns:
        if c in base or not pd.api.types.is_numeric_dtype(df[c]):
            continue
        if df[c].isna().all():
            continue
        out.append(c)
    return out


def dead_cols(df, cols):
    """全空列清单 (会被剔除, 仅用于报告透明)。"""
    return [c for c in cols if df[c].isna().all()]


def prepare(df, cols):
    """取列 -> 丢弃全空列 -> log1p 压缩计数型 -> 中位数兜底缺失。"""
    keep = [c for c in cols if not df[c].isna().all()]
    X = df[keep].copy()
    for c in COUNT_COLS:
        if c in X.columns:
            X[c] = np.log1p(pd.to_numeric(X[c], errors="coerce").clip(lower=0))
    X = X.apply(pd.to_numeric, errors="coerce")
    med = X.median(numeric_only=True)
    return X.fillna(med).fillna(0.0)


def struct_group_key(df):
    """结构规模代理分组, 用于 GroupKFold 降低乐观偏差。"""
    import math
    keys = []
    for rec in df.to_dict("records"):
        a = int(rec.get("n_atoms_lig", 0) or 0)
        b = int(rec.get("n_atoms_rec", 0) or 0)
        ra = int(rec.get("n_res_lig", 0) or 0)
        rb = int(rec.get("n_res_rec", 0) or 0)
        bundle = sorted([math.log1p(a), math.log1p(b)])
        keys.append("%d_%d_%d_%d" % (round(bundle[0]), round(bundle[1]),
                                     round(math.log1p(ra)), round(math.log1p(rb))))
    return np.asarray(keys)


def metrics(y_true, y_pred):
    sp = spearmanr(y_true, y_pred).correlation
    pe = pearsonr(y_true, y_pred)[0]
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    var = float(np.var(y_true, ddof=0))
    return {"spearman": float(sp) if sp == sp else float("nan"),
            "pearson": float(pe) if pe == pe else float("nan"),
            "rmse": rmse, "mae": float(mean_absolute_error(y_true, y_pred)),
            "r2": float(1 - rmse ** 2 / var) if var > 0 else float("nan")}


def bootstrap_ci(y_true, y_pred, n=1000, seed=SEED):
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n):
        idx = rs.randint(0, len(y_true), len(y_true))
        c = spearmanr(y_true[idx], y_pred[idx]).correlation
        if c == c:
            vals.append(c)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def run_config(tag, cols, X, y, tr_mask, te_mask, groups_tr, log):
    """某配置 x 某划分: GroupKFold CV -> 全训练集重训 -> 测试集评估。"""
    X_tr, y_tr = X[tr_mask], y[tr_mask]
    X_te, y_te = X[te_mask], y[te_mask]
    n_splits = min(N_SPLITS, max(2, len(set(groups_tr))))
    gkf = GroupKFold(n_splits=n_splits)
    oof = np.full(len(X_tr), np.nan)
    fold_sp = []
    for i_tr, i_va in gkf.split(X_tr, y_tr, groups_tr):
        m = lgb.LGBMRegressor(**PARAMS)
        m.fit(X_tr.iloc[i_tr], y_tr[i_tr], eval_set=[(X_tr.iloc[i_va], y_tr[i_va])],
              callbacks=[lgb.early_stopping(50, verbose=False)])
        oof[i_va] = m.predict(X_tr.iloc[i_va])
        fold_sp.append(spearmanr(y_tr[i_va], oof[i_va]).correlation)
        gc.collect()
    oof_m = metrics(y_tr, oof)

    final = lgb.LGBMRegressor(**PARAMS)
    final.fit(X_tr, y_tr)
    pred = final.predict(X_te)
    tm = metrics(y_te, pred)
    lo, hi = bootstrap_ci(y_te, pred)

    # 置换基线: 打乱训练标签后重跑同一 CV, 给出 OOF Spearman 的随机水平
    rng = np.random.RandomState(SEED)
    perm = []
    for _ in range(12):
        yp = rng.permutation(y_tr)
        o = np.full(len(X_tr), np.nan)
        for i_tr, i_va in gkf.split(X_tr, yp, groups_tr):
            mm = lgb.LGBMRegressor(**PARAMS)
            mm.fit(X_tr.iloc[i_tr], yp[i_tr])
            o[i_va] = mm.predict(X_tr.iloc[i_va])
        c = spearmanr(yp, o).correlation
        if c == c:
            perm.append(c)

    log("  %-26s 特征=%3d | CV Spearman=%+.3f±%.3f | OOF=%+.3f | 置换基线=%+.3f±%.3f | "
        "测试 Spearman=%+.3f [%+.3f,%+.3f] RMSE=%.3f R2=%.3f"
        % (tag, X_tr.shape[1], np.nanmean(fold_sp), np.nanstd(fold_sp),
           oof_m["spearman"], np.nanmean(perm), np.nanstd(perm),
           tm["spearman"], lo, hi, tm["rmse"], tm["r2"]))
    return {"tag": tag, "n_features": int(X_tr.shape[1]), "cols": list(cols),
            "cv_spearman": float(np.nanmean(fold_sp)),
            "cv_spearman_std": float(np.nanstd(fold_sp)),
            "oof": oof_m, "test": tm, "test_ci": [lo, hi],
            "perm_mean": float(np.nanmean(perm)), "perm_std": float(np.nanstd(perm)),
            "n_train": int(tr_mask.sum()), "n_test": int(te_mask.sum())}, pred, final


# ---------------------------------------------------------------- 主流程
def main():
    t0 = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)
    fh_out = open(REPORT_PATH, "w", encoding="utf-8")
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(str(msg))
        fh_out.write(str(msg) + "\n")
        fh_out.flush()

    log("=" * 96)
    log("Affinity Intelligence V0.2 报告 — 真家族划分 + 无标签特征归一化")
    log("=" * 96)
    log("LightGBM %s (OpenMP 预加载: %s)" % (lgb.__version__, _VCOMP or "失败"))
    log("绘图字体: %s" % (_CJK_FONT or "未找到 CJK 字体 (图中中文可能显示为方框)"))

    df = pd.read_csv(IN_CSV)
    if "split" not in df.columns or "split_v0.1" not in df.columns:
        raise RuntimeError("特征表需同时含 split(V0.2) 与 split_v0.1 列")
    y = df[TARGET].to_numpy(dtype=float)
    tr2 = (df["split"] == "train").to_numpy()
    te2 = (df["split"] == "test").to_numpy()
    tr1 = (df["split_v0.1"] == "train").to_numpy()
    te1 = (df["split_v0.1"] == "test").to_numpy()

    lk = leaky_cols(df)
    nz = normalized_cols(df)
    base_cols = [c for c in BASE_FEATURES if c in df.columns]

    log("输入: %s" % IN_CSV)
    log("样本: %d 行" % len(df))
    log("  V0.1 划分 (PDB 级, 含同源泄漏): train=%d / test=%d"
        % (tr1.sum(), te1.sum()))
    log("  V0.2 划分 (40%% 同一性家族级)  : train=%d / test=%d"
        % (tr2.sum(), te2.sum()))
    log("特征分组: 原始物理 %d | 新增无标签归一化 %d | 标签派生(仅诊断) %d"
        % (len(base_cols), len(nz), len(lk)))
    log("  新增归一化特征: %s" % ", ".join(nz))
    log("  标签派生量 (分子含 real_dG, **不入模**): %s" % ", ".join(lk))
    dead = [c for c in df.columns if c not in META_COLS + [TARGET]
            and c not in base_cols + nz + lk and df[c].isna().all()]
    log("  已剔除的全空列: %s" % (dead or "无"))

    groups_all = struct_group_key(df)
    all_cols = base_cols + nz

    log("")
    log("[1] 对比矩阵 (同一数据, 两种划分 x 三个特征配置)")
    log("    —— 配置: A=原始物理特征 | B=A+无标签归一化 | C=B+标签派生量(泄漏对照)")
    results, preds, models = {}, {}, {}
    for split_name, tr_m, te_m, gtr in (
            ("V0.1", tr1, te1, groups_all[tr1]),
            ("V0.2", tr2, te2, groups_all[tr2])):
        for cfg_name, cols in (("A", base_cols), ("B", all_cols),
                               ("C", all_cols + lk)):
            tag = "%s x %s" % (split_name, cfg_name)
            X = prepare(df, cols)
            res, pred, mdl = run_config(tag, cols, X, y, tr_m, te_m, gtr, log)
            res.update({"split": split_name, "cfg": cfg_name,
                        "y_test": y[te_m],
                        "pdb_test": df.loc[te_m, "pdb_id"].to_numpy()})
            results[tag] = res
            preds[tag] = pred
            models[tag] = mdl
            gc.collect()
        log("")

    # ---- 影响拆解 -------------------------------------------------------
    a1, b1, c1 = (results["V0.1 x A"], results["V0.1 x B"], results["V0.1 x C"])
    a2, b2, c2 = (results["V0.2 x A"], results["V0.2 x B"], results["V0.2 x C"])

    log("[2] 影响拆解 (测试集 Spearman)")
    log("    V0.1 x A (旧基线, PDB 级 + 无归一化)      = %+.3f  <- V0.1 报告数字" % a1["test"]["spearman"])
    log("    V0.2 x A (家族级 + 无归一化)             = %+.3f" % a2["test"]["spearman"])
    log("    V0.2 x B (家族级 + 归一化, **主模型**)   = %+.3f" % b2["test"]["spearman"])
    log("    ---")
    log("    ≥ 划分变严的代价 (V0.2xA - V0.1xA)       = %+.3f"
        % (a2["test"]["spearman"] - a1["test"]["spearman"]))
    log("    ≥ 归一化特征的增益 (V0.2xB - V0.2xA)     = %+.3f"
        % (b2["test"]["spearman"] - a2["test"]["spearman"]))
    log("    ≥ 净变化 (V0.2xB - V0.1xA)               = %+.3f"
        % (b2["test"]["spearman"] - a1["test"]["spearman"]))
    log("    ≥ 目标泄漏的虚假增益 (V0.2xC - V0.2xB)   = %+.3f  (不可采信)"
        % (c2["test"]["spearman"] - b2["test"]["spearman"]))
    log("")
    log("    同口径 RMSE / R2:")
    for r in (a1, a2, b2, c2):
        log("      %-12s RMSE=%.3f  R2=%+.3f" % (r["tag"], r["test"]["rmse"], r["test"]["r2"]))

    # ---- 特征重要性 (主模型) -------------------------------------------
    imp = pd.DataFrame({
        "feature": all_cols,
        "importance_gain": models["V0.2 x B"].booster_.feature_importance("gain"),
        "importance_split": models["V0.2 x B"].booster_.feature_importance("split"),
    }).sort_values("importance_gain", ascending=False).reset_index(drop=True)
    imp["rank"] = np.arange(1, len(imp) + 1)
    imp["is_normalized"] = imp["feature"].isin(nz)
    imp["was_in_v01"] = imp["feature"].isin(BASE_FEATURES)
    imp.to_csv(IMP_PATH, index=False, encoding="utf-8")

    log("")
    log("[3] V0.2 主模型特征重要性 Top 12 (gain):")
    for r in imp.head(12).to_dict("records"):
        log("    %2d. %-26s gain=%9.1f split=%4d%s"
            % (r["rank"], r["feature"], r["importance_gain"],
               r["importance_split"], "  [新增归一化]" if r["is_normalized"] else ""))

    size_like = ["n_contact_pairs", "n_contact_res_pairs", "sasa_complex",
                 "dsasa_total", "n_atoms_lig", "n_atoms_rec", "n_iface_res_lig",
                 "n_iface_res_rec", "n_res_lig", "n_res_rec", "n_chainpair_max",
                 "n_chainpair_min", "contact_density"]
    top10 = list(imp.head(10)["feature"])
    log("")
    log("    规模类特征进入 Top10: %s" % (sorted(set(top10) & set(size_like)) or "无"))
    log("    被挤出 Top10 的规模类特征: %s"
        % (sorted(set(size_like) - set(top10)) or "无"))
    log("    Top10 中的归一化特征: %s" % (sorted(set(top10) & set(nz)) or "无"))
    log("    归一化特征排名区间: %s"
        % (sorted(imp.loc[imp["is_normalized"], "rank"].tolist())[:12]))

    # ---- 保存 -----------------------------------------------------------
    payload = {
        "model": models["V0.2 x B"],
        "feature_names": all_cols,
        "params": PARAMS,
        "target": TARGET,
        "label_transform": "dG = 0.592*ln(Kd[M]) [kcal/mol]",
        "excluded_leaky_features": lk,
        "normalized_features": nz,
        "split_policy": "V0.2 家族级: 40% 序列同一性单链接聚类, 整族分侧",
        "train_rows": int(tr2.sum()), "test_rows": int(te2.sum()),
        "metrics": {k: v["test"] for k, v in results.items()},
        "count_features_log1p": [c for c in COUNT_COLS if c in all_cols],
    }
    with open(MODEL_PATH, "wb") as fh:
        pickle.dump(payload, fh)
    log("")
    log("[4] 模型已保存: %s (%.1f KB)" % (MODEL_PATH, os.path.getsize(MODEL_PATH) / 1024))

    pred_df = pd.DataFrame({
        "pdb_id": df.loc[te2, "pdb_id"].to_numpy(),
        "interface_key": df.loc[te2, "interface_key"].to_numpy(),
        "real_dG": y[te2], "pred_dG": preds["V0.2 x B"],
        "source_dataset": df.loc[te2, "source_dataset"].to_numpy(),
    })
    pred_df["residual"] = pred_df["pred_dG"] - pred_df["real_dG"]
    pred_df.sort_values("real_dG").to_csv(PRED_PATH, index=False, encoding="utf-8")

    # ---- 对比图 ---------------------------------------------------------
    fig = plt.figure(figsize=(16.5, 10.4))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.32, wspace=0.26)

    def scatter(ax, res, color, title):
        yt, pt = res["y_test"], preds[res["tag"]]
        ax.scatter(yt, pt, s=28, alpha=0.7, color=color,
                   edgecolors="white", linewidths=0.5)
        lim = [min(yt.min(), pt.min()) - 1.5, max(yt.max(), pt.max()) + 1.5]
        ax.plot(lim, lim, "--", color="#a0aec0", lw=1.0)
        z = np.polyfit(yt, pt, 1)
        xs = np.linspace(lim[0], lim[1], 40)
        ax.plot(xs, np.polyval(z, xs), "-", color="#c53030", lw=1.4,
                label="slope=%.2f" % z[0])
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("real_dG (kcal/mol)", fontsize=8.5)
        ax.set_ylabel("predicted_dG (kcal/mol)", fontsize=8.5)
        ax.set_title("%s\nSpearman=%+.3f  RMSE=%.2f  R2=%+.2f  (n=%d)"
                     % (title, res["test"]["spearman"], res["test"]["rmse"],
                        res["test"]["r2"], res["n_test"]), fontsize=9.5)
        ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.25)

    scatter(fig.add_subplot(grid[0, 0]), a1, "#718096",
            "V0.1 划分 x A 原始特征 (旧基线)")
    scatter(fig.add_subplot(grid[0, 1]), b2, "#2b6cb0",
            "V0.2 主模型: 家族划分 x B 归一化")
    scatter(fig.add_subplot(grid[0, 2]), c2, "#c53030",
            "V0.2 划分 x C 含标签泄漏 (不可采信)")

    ax = fig.add_subplot(grid[1, 0])
    order = ["V0.1 x A", "V0.1 x B", "V0.1 x C", "V0.2 x A", "V0.2 x B", "V0.2 x C"]
    vals = [results[k]["test"]["spearman"] for k in order]
    cols = ["#a0aec0", "#90cdf4", "#feb2b2", "#718096", "#2b6cb0", "#c53030"]
    bars = ax.bar(np.arange(len(order)), vals, color=cols, alpha=0.92, width=0.62)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.012,
                "%+.3f" % v, ha="center", fontsize=8.5)
    ax.axhline(0, color="#4a5568", lw=0.8)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([k.replace(" ", "\n") for k in order], fontsize=8)
    ax.set_ylabel("test Spearman"); ax.set_title("六种组合的核心指标", fontsize=10)
    ax.grid(alpha=0.25, axis="y")

    ax = fig.add_subplot(grid[1, 1])
    labels = ["A 原始\n特征", "B +归一化\n(主模型)", "C +标签泄漏\n(对照)"]
    r2 = [a2["test"]["rmse"], b2["test"]["rmse"], c2["test"]["rmse"]]
    ax.bar(np.arange(3), r2, color=["#718096", "#2b6cb0", "#c53030"],
           alpha=0.92, width=0.55)
    for i, v in enumerate(r2):
        ax.text(i, v + 0.03, "%.3f" % v, ha="center", fontsize=9)
    ax.set_xticks(np.arange(3)); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("test RMSE (kcal/mol)")
    ax.set_title("V0.2 划分下三配置的 RMSE", fontsize=10)
    ax.grid(alpha=0.25, axis="y")

    ax = fig.add_subplot(grid[1, 2])
    top = imp.head(14).iloc[::-1]
    ax.barh(top["feature"], top["importance_gain"],
            color=["#2f855a" if v else "#718096" for v in top["is_normalized"]],
            alpha=0.9)
    ax.set_xlabel("importance (gain)")
    ax.set_title("V0.2 主模型特征重要性 Top14\n(绿 = 新增无标签归一化特征)", fontsize=10)
    ax.tick_params(axis="y", labelsize=7.5); ax.grid(alpha=0.25, axis="x")

    fig.suptitle("Affinity Intelligence: V0.1 (PDB 级划分) vs V0.2 (40% 同一性家族划分 + 特征归一化)",
                 fontsize=12.5, y=0.985)
    plt.savefig(FIG_PATH, dpi=145, bbox_inches="tight")
    plt.close()
    log("    对比图: %s" % FIG_PATH)
    log("    特征重要性: %s" % IMP_PATH)
    log("    测试集预测: %s" % PRED_PATH)

    # ---- 结论 -----------------------------------------------------------
    m = b2["test"]
    log("")
    log("[5] 结论")
    log("    * V0.2 主模型 (V0.2 家族划分 + 无标签归一化):")
    log("      Spearman=%+.3f (95%% CI [%+.3f, %+.3f]), RMSE=%.3f kcal/mol, MAE=%.3f, R2=%.3f"
        % (m["spearman"], b2["test_ci"][0], b2["test_ci"][1], m["rmse"],
           m["mae"], m["r2"]))
    log("    * 相对 V0.1 旧基线: Spearman %+.3f -> %+.3f (Δ%+.3f); RMSE %.3f -> %.3f"
        % (a1["test"]["spearman"], m["spearman"],
           m["spearman"] - a1["test"]["spearman"],
           a1["test"]["rmse"], m["rmse"]))
    log("    * 划分变严贡献 ΔSpearman=%+.3f, 归一化特征贡献 Δ=%+.3f"
        % (a2["test"]["spearman"] - a1["test"]["spearman"],
           b2["test"]["spearman"] - a2["test"]["spearman"]))
    log("    * 目标泄漏 (C) 会把 Spearman 抬到 %+.3f (虚高 %+.3f) —— 已从主模型剔除"
        % (c2["test"]["spearman"], c2["test"]["spearman"] - b2["test"]["spearman"]))
    log("    * 置换检验对照 (主模型): OOF=%+.3f vs 随机基线 %+.3f±%.3f (超 %.1fσ)"
        % (b2["oof"]["spearman"], b2["perm_mean"], b2["perm_std"],
           (b2["oof"]["spearman"] - b2["perm_mean"]) / max(b2["perm_std"], 1e-6)))
    log("    * 残留局限: 40%% 同一性仍是同源代理 (远缘同折叠未覆盖); 训练 %d 条;"
        % int(tr2.sum()))
    log("                PRODIGY 不可用, predicted_dG 缺失; 测试集仅 %d 条, 置信区间宽。"
        % int(te2.sum()))
    log("")
    log("耗时: %.1fs" % (time.time() - t0))

    fh_out.close()
    print("报告: %s" % REPORT_PATH)
    print("")
    print("阶段三完成，等待下一步指令")


if __name__ == "__main__":
    main()
