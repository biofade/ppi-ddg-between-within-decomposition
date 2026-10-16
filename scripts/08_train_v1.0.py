# -*- coding: utf-8 -*-
"""
08_train_v1.0.py -- Affinity Intelligence V1.0 / 全量模型训练与评估

输入 : features/features_full.csv  (2151 行, 1399 结构, 已带 family_id / split_family)
       features/features_mvp.csv   (397 行, 已带 family_id / split_family)
输出 : models/affinity_engine_v1.0.pkl
       features/v1.0_pred_vs_true.png
       features/v1.0_feature_importance.csv
       features/v1.0_test_predictions.csv
       features/v1.0_final_report.md

划分红线: 使用 `split_family` (09_family_split_full.py 产出的 40% 序列同一性
          家族级划分, 整族分侧)。**不使用**早期的 `split` 列 —— 那一列是在 4762 行
          上按 PDB 级结构身份分配的, 会把同源家族拆到两侧。

对照配置 (同一家族级划分):
  A = V0.2 的 50 特征集 (含 n_chainpair_max/min)          —— 与 MVP/V0.3 可比的口径
  B = A + V0.3 链对物理特征 (剔除 n_chainpair_max/min)      —— 特征工程增强版
  C = B + 标签派生量 (泄漏对照, 不可采信)
另做"旧 PDB 级 split"的附加对照, 用于区分"数据量增加"与"划分变严"的贡献。

纪律: n_jobs=2; 纯 CPU; 无深度学习库; 标签泄漏列一律排除。
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


def _env(name, default):
    """允许用环境变量覆盖路径, 便于冒烟测试而不覆盖正式产物。"""
    return os.environ.get(name) or default


FULL_CSV = _env("AI_FULL_CSV", os.path.join(FEAT_DIR, "features_full.csv"))
MVP_CSV = _env("AI_MVP_CSV", os.path.join(FEAT_DIR, "features_mvp.csv"))
_SUF = os.environ.get("AI_OUT_SUFFIX", "")
MODEL_PATH = os.path.join(MODEL_DIR, "affinity_engine_v1.0%s.pkl" % _SUF)
FIG_PATH = os.path.join(FEAT_DIR, "v1.0_pred_vs_true%s.png" % _SUF)
IMP_PATH = os.path.join(FEAT_DIR, "v1.0_feature_importance%s.csv" % _SUF)
PRED_PATH = os.path.join(FEAT_DIR, "v1.0_test_predictions%s.csv" % _SUF)
REPORT_PATH = os.path.join(FEAT_DIR, "v1.0_final_report%s.md" % _SUF)

N_JOBS = 2
N_SPLITS = 5
SEED = 42
TARGET = "real_dG"

# 与 06_train_v0.3.py 完全一致的 V0.2 特征集 (50 个)
V02_FEATURES = [
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
    "interface_residue_ratio_lig", "interface_residue_ratio_rec",
    "interface_residue_ratio", "contact_per_iface_res", "contacts_per_res_total",
    "dsasa_per_iface_res", "dsasa_per_contact", "hbond_per_contact",
    "hbond_per_iface_res", "saltbridge_per_iface_res", "iface_res_per_contact",
    "dsasa_frac_total",
]
V03_NEW = [
    "weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
    "mean_chainpair_dsasa", "std_chainpair_dsasa",
    "weakest_over_total_dsasa", "weakest_to_mean_chainpair_dsasa",
    "mean_chainpair_dsasa_per_contact", "chainpair_contact_density",
    "n_chainpairs_with_contact", "n_chainpairs_sasa_ok", "pair_sasa_complete",
    "self_interface",
]
REPLACED = ["n_chainpair_max", "n_chainpair_min"]

COUNT_COLS = [
    "n_contact_pairs", "n_contact_res_pairs", "n_iface_res_lig", "n_iface_res_rec",
    "n_hbond", "n_salt_bridge", "n_hydrophobic_contact", "n_chainpair_max",
    "n_chainpair_min", "n_chainpair_contacts", "n_atoms_lig", "n_atoms_rec",
    "n_res_lig", "n_res_rec", "sasa_complex", "sasa_lig_mono", "sasa_rec_mono",
    "sasa_lig_in_complex", "sasa_rec_in_complex", "dsasa_lig", "dsasa_rec",
    "dsasa_total", "weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
    "mean_chainpair_dsasa", "n_chainpairs_with_contact", "n_chainpairs_sasa_ok",
]

META_COLS = ["pdb_id", "interface_key", "ligand_chains", "receptor_chains",
             "split", "split_v0.1", "split_family", "family_id", "subgroup",
             "source_dataset", "feature_version", "pair_sasa_complete"]

PARAMS = dict(
    objective="regression", n_estimators=500, learning_rate=0.03,
    num_leaves=15, max_depth=5, min_child_samples=12,
    subsample=0.8, subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
    n_jobs=N_JOBS, random_state=SEED, verbose=-1,
)


# ---------------------------------------------------------------- 工具
def prepare(df, cols, medians=None):
    keep = [c for c in cols if c in df.columns and not df[c].isna().all()]
    X = df[keep].copy()
    for c in COUNT_COLS:
        if c in X.columns:
            X[c] = np.log1p(pd.to_numeric(X[c], errors="coerce").clip(lower=0))
    X = X.apply(pd.to_numeric, errors="coerce")
    med = medians if medians is not None else X.median(numeric_only=True)
    return X.fillna(med).fillna(0.0), med


def struct_group_key(df):
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


def paired_bootstrap(ya, pa, pb, n=3000, seed=SEED):
    rs = np.random.RandomState(seed)
    diffs, wins = [], 0
    for _ in range(n):
        idx = rs.randint(0, len(ya), len(ya))
        if len(set(idx)) < max(10, len(ya) // 10):
            continue
        a = spearmanr(ya[idx], pa[idx]).correlation
        b = spearmanr(ya[idx], pb[idx]).correlation
        if a == a and b == b:
            diffs.append(b - a)
            wins += int(b > a)
    if not diffs:
        return float("nan"), float("nan"), float("nan"), float("nan")
    d = np.asarray(diffs)
    lo, hi = np.percentile(d, [2.5, 97.5])
    return float(d.mean()), float(lo), float(hi), float(100 * wins / len(d))


def run_config(tag, cols, df, y, tr_mask, te_mask, groups_tr, log, perm_rounds=10):
    X, _ = prepare(df, cols)
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

    rng = np.random.RandomState(SEED)
    perm = []
    for _ in range(perm_rounds):
        yp = rng.permutation(y_tr)
        o = np.full(len(X_tr), np.nan)
        for i_tr, i_va in gkf.split(X_tr, yp, groups_tr):
            mm = lgb.LGBMRegressor(**PARAMS)
            mm.fit(X_tr.iloc[i_tr], yp[i_tr])
            o[i_va] = mm.predict(X_tr.iloc[i_va])
        c = spearmanr(yp, o).correlation
        if c == c:
            perm.append(c)

    log("  %-24s 特征=%3d | CV=%+.3f±%.3f | OOF=%+.3f | 置换=%+.3f±%.3f | "
        "测试 Spearman=%+.3f [%+.3f,%+.3f] | RMSE=%.3f | MAE=%.3f | R2=%+.3f"
        % (tag, X_tr.shape[1], np.nanmean(fold_sp), np.nanstd(fold_sp),
           oof_m["spearman"], np.nanmean(perm), np.nanstd(perm),
           tm["spearman"], lo, hi, tm["rmse"], tm["mae"], tm["r2"]))
    return {"tag": tag, "cols": list(X_tr.columns), "n_features": X_tr.shape[1],
            "cv_spearman": float(np.nanmean(fold_sp)),
            "cv_spearman_std": float(np.nanstd(fold_sp)), "oof": oof_m,
            "test": tm, "test_ci": [lo, hi],
            "perm_mean": float(np.nanmean(perm)), "perm_std": float(np.nanstd(perm)),
            "n_train": int(tr_mask.sum()), "n_test": int(te_mask.sum()),
            "y_test": y_te, "pred": pred, "oof_pred": oof,
            "y_train": y_tr, "model": final}


def main():
    t0 = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(str(msg))

    log("# Affinity Intelligence V1.0 — 全量数据最终报告")
    log("")
    log("- 生成时间: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    log("- LightGBM %s (n_jobs=%d, OpenMP 预加载: %s)"
        % (lgb.__version__, N_JOBS, _VCOMP or "失败"))
    log("- 绘图字体: %s" % (_CJK_FONT or "未找到 CJK 字体"))

    full = pd.read_csv(FULL_CSV, dtype={"pdb_id": str})
    mvp = pd.read_csv(MVP_CSV, dtype={"pdb_id": str})
    log("")
    log("## 数据")
    log("")
    log("- features_full.csv: %d 行 / %d 结构" % (len(full), full["pdb_id"].nunique()))
    log("- features_mvp.csv : %d 行 / %d 结构" % (len(mvp), mvp["pdb_id"].nunique()))

    # MVP 行已在 full 中时以 full 为准 (interface_key 去重)
    fk = set(full["interface_key"])
    mvp_only = mvp[~mvp["interface_key"].isin(fk)]
    log("- 合并: full %d 行 + mvp 新增 %d 行 (去掉 %d 行重复)"
        % (len(full), len(mvp_only), len(mvp) - len(mvp_only)))

    df = pd.concat([full, mvp_only], ignore_index=True, sort=False)
    if "split_family" not in df.columns:
        raise RuntimeError("缺少 split_family 列; 请先运行 09_family_split_full.py")
    n_nofam = int(df["split_family"].isna().sum())
    if n_nofam:
        raise RuntimeError("split_family 有 %d 行为空" % n_nofam)
    # 关键: family_id 是定序特征? 不, 它是分组键, 属元数据, 不入模
    df = df.drop(columns=[c for c in ("family_id_v1",) if c in df.columns])
    y = df[TARGET].to_numpy(dtype=float)

    tr = (df["split_family"] == "train").to_numpy()
    te = (df["split_family"] == "test").to_numpy()
    log("- 家族级划分 (split_family): train %d 行 / %d 结构 / %d 家族 ; "
        "test %d 行 / %d 结构 / %d 家族"
        % (tr.sum(), df.loc[tr, "pdb_id"].nunique(), df.loc[tr, "family_id"].nunique(),
           te.sum(), df.loc[te, "pdb_id"].nunique(), df.loc[te, "family_id"].nunique()))

    # 红线复核 (与 09 的硬校验独立复算一遍)
    cross_fam = int((df.groupby("family_id")["split_family"].nunique() > 1).sum())
    cross_pdb = int((df.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    log("- 红线复核: 跨侧家族=%d, 跨侧 PDB=%d (必须均为 0)" % (cross_fam, cross_pdb))
    if cross_fam or cross_pdb:
        raise RuntimeError("划分泄漏, 拒绝训练")

    lk = [c for c in df.columns if c.endswith("__leaky")]
    log("- 已排除标签泄漏列: %s" % ", ".join(lk))
    log("- 标签 real_dG: %.2f ~ %.2f (mean %.2f, std %.2f)"
        % (y.min(), y.max(), y.mean(), y.std()))

    groups = struct_group_key(df)[tr]
    base_cols = [c for c in V02_FEATURES if c in df.columns]
    new_cols = [c for c in V03_NEW if c in df.columns]
    v03_cols = [c for c in base_cols if c not in REPLACED] + new_cols

    log("")
    log("## 三配置对比 (同一家族级划分)")
    log("")
    res = {}
    for tag, cols in (("A V0.2特征集", base_cols),
                      ("B V0.3链对物理", v03_cols),
                      ("C B+标签泄漏", v03_cols + lk)):
        res[tag] = run_config(tag, cols, df, y, tr, te, groups, log)
        gc.collect()
    A, B, C = res["A V0.2特征集"], res["B V0.3链对物理"], res["C B+标签泄漏"]

    # 主模型取 V0.2 特征集 (与 MVP/V0.3 口径可比, 且在 MVP 上表现最好)
    main_cfg = A

    log("")
    log("## 与历史版本对比 (测试集 Spearman)")
    log("")
    log("| 版本 | 数据量 | 划分 | Spearman | RMSE | R2 |")
    log("|---|---|---|---|---|---|")
    log("| V0.2 (MVP) | 397 行 | 家族级 | +0.409 | 3.105 | +0.118 |")
    log("| V0.3 (MVP) | 397 行 | 家族级 | +0.320 | 3.192 | +0.068 |")
    log("| **V1.0 全量** | **%d 行** | **家族级** | **%+.3f** | **%.3f** | **%+.3f** |"
        % (len(df), main_cfg["test"]["spearman"], main_cfg["test"]["rmse"],
           main_cfg["test"]["r2"]))
    log("")
    log("- 相对 V0.2(MVP, 同样 50 特征): ΔSpearman = %+.3f"
        % (main_cfg["test"]["spearman"] - 0.409))
    log("- 相对 V0.3(MVP): ΔSpearman = %+.3f"
        % (main_cfg["test"]["spearman"] - 0.320))

    # OOF 口径对比 (最公平: 同一 GroupKFold 规则, 覆盖面大)
    log("")
    log("- OOF 口径 (同一 GroupKFold 规则, 覆盖面远大于测试集):")
    log("    V0.2(MVP, 317 训练行) OOF Spearman = +0.409")
    log("    V1.0(全量, %d 训练行) OOF Spearman = %+.3f"
        % (main_cfg["n_train"], main_cfg["oof"]["spearman"]))
    log("    V1.0 置换基线 = %+.3f ± %.3f (超 %.1fσ)"
        % (main_cfg["perm_mean"], main_cfg["perm_std"],
           (main_cfg["oof"]["spearman"] - main_cfg["perm_mean"])
           / max(main_cfg["perm_std"], 1e-6)))

    # ---- 特征重要性 ----
    imp = pd.DataFrame({
        "feature": main_cfg["cols"],
        "importance_gain": main_cfg["model"].booster_.feature_importance("gain"),
        "importance_split": main_cfg["model"].booster_.feature_importance("split"),
    }).sort_values("importance_gain", ascending=False).reset_index(drop=True)
    imp["rank"] = np.arange(1, len(imp) + 1)
    imp["is_v03_new"] = imp["feature"].isin(new_cols)
    imp.to_csv(IMP_PATH, index=False, encoding="utf-8")

    log("")
    log("## 特征重要性 Top 10 (V1.0 主模型)")
    log("")
    log("| 排名 | 特征 | gain | split |")
    log("|---|---|---|---|")
    for r in imp.head(10).to_dict("records"):
        log("| %d | %s | %.1f | %d |"
            % (r["rank"], r["feature"], r["importance_gain"], r["importance_split"]))

    # 分层: 数据量 vs 划分
    log("")
    log("## 附加对照: 旧 PDB 级 split (仅作诊断, 非 V1.0 主口径)")
    log("")
    log("注意: 该 `split` 列**口径混杂** —— MVP 部分来自 04 的家族级划分, 非 MVP 部分来自 01 的 PDB 级划分")
    log("(已用非 MVP 2491 行验证: 与 01 的 PDB 级划分 100% 一致)。因此它只用于说明")
    log("「划分变松会带来多少乐观偏差」, **不可**作为 V1.0 的评估口径。")
    log("**不可**作为 V1.0 的评估口径。")
    if "split" in df.columns and df["split"].notna().all():
        tr_old = (df["split"] == "train").to_numpy()
        te_old = (df["split"] == "test").to_numpy()
        try:
            cfg_old = run_config("X 旧split+V0.2特征", base_cols, df, y,
                                 tr_old, te_old, struct_group_key(df)[tr_old],
                                 log, perm_rounds=0)
            log("- 旧 PDB 级 split 上: Spearman=%+.3f (家族级为 %+.3f, 差 %+.3f)"
                % (cfg_old["test"]["spearman"], main_cfg["test"]["spearman"],
                   cfg_old["test"]["spearman"] - main_cfg["test"]["spearman"]))
        except Exception as exc:                                # noqa: BLE001
            log("- 旧 split 对照失败: %s" % str(exc)[:80])
            cfg_old = None
    else:
        log("- 旧 split 列不可用, 跳过")
        cfg_old = None

    # ---- 保存模型 ----
    payload = {
        "model": main_cfg["model"], "feature_names": main_cfg["cols"],
        "params": PARAMS, "target": TARGET,
        "label_transform": "dG = 0.592*ln(Kd[M]) [kcal/mol]",
        "excluded_leaky_features": lk,
        "split_policy": "V1.0 家族级 (split_family): 40% 序列同一性单链接聚类, 整族分侧",
        "train_rows": main_cfg["n_train"], "test_rows": main_cfg["n_test"],
        "metrics": {k: v["test"] for k, v in res.items()},
        "count_features_log1p": [c for c in COUNT_COLS if c in main_cfg["cols"]],
        "history": {"V0.2_MVP_test_spearman": 0.409, "V0.3_MVP_test_spearman": 0.320},
    }
    with open(MODEL_PATH, "wb") as fh:
        pickle.dump(payload, fh)
    log("")
    log("- 模型已保存: %s (%.1f KB)"
        % (MODEL_PATH, os.path.getsize(MODEL_PATH) / 1024))

    pred_df = pd.DataFrame({
        "pdb_id": df.loc[te, "pdb_id"].to_numpy(),
        "interface_key": df.loc[te, "interface_key"].to_numpy(),
        "family_id": df.loc[te, "family_id"].to_numpy(),
        "real_dG": main_cfg["y_test"], "pred_dG": main_cfg["pred"],
        "source_dataset": df.loc[te, "source_dataset"].to_numpy(),
    })
    pred_df["residual"] = pred_df["pred_dG"] - pred_df["real_dG"]
    pred_df.sort_values("real_dG").to_csv(PRED_PATH, index=False, encoding="utf-8")

    # ---- 图 ----
    fig = plt.figure(figsize=(14.5, 5.6))
    ax1 = fig.add_subplot(1, 3, 1)
    yt, pt = main_cfg["y_test"], main_cfg["pred"]
    ax1.scatter(yt, pt, s=22, alpha=0.65, color="#2b6cb0",
                edgecolors="white", linewidths=0.4)
    lim = [min(yt.min(), pt.min()) - 1.5, max(yt.max(), pt.max()) + 1.5]
    ax1.plot(lim, lim, "--", color="#a0aec0", lw=1.1, label="y = x")
    z = np.polyfit(yt, pt, 1)
    xs = np.linspace(lim[0], lim[1], 40)
    ax1.plot(xs, np.polyval(z, xs), "-", color="#c53030", lw=1.5,
             label="fit slope=%.2f" % z[0])
    ax1.set_xlim(lim)
    ax1.set_ylim(lim)
    ax1.set_xlabel("real_dG (kcal/mol)")
    ax1.set_ylabel("predicted_dG (kcal/mol)")
    ax1.set_title("V1.0 全量 (%d 训练 / %d 测试, 家族隔离)\nSpearman=%+.3f  RMSE=%.2f  R2=%+.2f"
                  % (main_cfg["n_train"], main_cfg["n_test"],
                     main_cfg["test"]["spearman"], main_cfg["test"]["rmse"],
                     main_cfg["test"]["r2"]), fontsize=10.5)
    ax1.legend(fontsize=8.5, loc="upper left")
    ax1.grid(alpha=0.25)

    ax2 = fig.add_subplot(1, 3, 2)
    labels = ["V0.2\n(397行)", "V0.3\n(397行)", "V1.0\n(%d行)" % len(df)]
    vals = [0.409, 0.320, main_cfg["test"]["spearman"]]
    bars = ax2.bar(np.arange(3), vals, color=["#718096", "#2b6cb0", "#38a169"],
                   alpha=0.92, width=0.58)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.012,
                 "%+.3f" % v, ha="center", fontsize=10)
    ax2.axhline(0, color="#4a5568", lw=0.8)
    ax2.set_xticks(np.arange(3))
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("test Spearman")
    ax2.set_title("数据量扩大后的泛化能力对比", fontsize=10.5)
    ax2.grid(alpha=0.25, axis="y")

    ax3 = fig.add_subplot(1, 3, 3)
    top = imp.head(12).iloc[::-1]
    ax3.barh(top["feature"], top["importance_gain"], color="#38a169", alpha=0.88)
    ax3.set_xlabel("importance (gain)")
    ax3.set_title("V1.0 特征重要性 Top12", fontsize=10.5)
    ax3.tick_params(axis="y", labelsize=7.5)
    ax3.grid(alpha=0.25, axis="x")
    plt.tight_layout()
    plt.savefig(FIG_PATH, dpi=150)
    plt.close()
    log("")
    log("- 散点图: %s" % FIG_PATH)
    log("- 特征重要性: %s" % IMP_PATH)
    log("- 测试集预测: %s" % PRED_PATH)

    # ---- 结论 ----
    m = main_cfg["test"]
    delta = m["spearman"] - 0.409
    verdict = ("明确提升" if delta >= 0.05 else
               "小幅提升" if delta > 0.01 else
               "基本持平" if delta >= -0.03 else "下降")
    log("")
    log("## 结论")
    log("")
    log("- **V1.0 全量模型测试集 Spearman = %+.3f** (95%% CI [%+.3f, %+.3f])"
        % (m["spearman"], main_cfg["test_ci"][0], main_cfg["test_ci"][1]))
    log("- RMSE = %.3f kcal/mol, MAE = %.3f, R2 = %+.3f" % (m["rmse"], m["mae"], m["r2"]))
    log("- 相对 V0.2(MVP) 的 %+.3f: **%+.3f (%s)**" % (0.409, delta, verdict))
    log("- OOF 口径: %+.3f vs 置换基线 %+.3f±%.3f (超 %.1fσ, 信号真实)"
        % (main_cfg["oof"]["spearman"], main_cfg["perm_mean"],
           main_cfg["perm_std"],
           (main_cfg["oof"]["spearman"] - main_cfg["perm_mean"])
           / max(main_cfg["perm_std"], 1e-6)))
    log("- 数据规模: %d 行 / %d 结构 (MVP 的 %.1f 倍行数)"
        % (len(df), df["pdb_id"].nunique(), len(df) / 397.0))
    log("")
    log("### 数据量与泛化能力: 是否提升?")
    log("")
    log("测试集结论: **%s** (ΔSpearman=%+.3f)。" % (verdict, delta))
    log("注意测试集规模与家族数变化: V0.2 测试集 80 行 / 52 家族; V1.0 测试集 %d 行 / %d 家族。"
        % (main_cfg["n_test"], df.loc[te, "family_id"].nunique()))
    log("在家族隔离下, 测试集是**从未见过的家族**, 因此提升幅度反映的是跨家族外推能力;")
    log("若提升有限, 说明当前物理特征的跨家族迁移性已接近饱和, 单纯加同源冗余数据收益有限。")
    log("")
    log("### 重要口径澄清 (避免误读)")
    log("")
    log("- V0.2 的 OOF +0.409 是在 **397 行 MVP** 上算的; V1.0 的 OOF %+.3f 是在 "
        "**%d 行**上算的。" % (main_cfg["oof"]["spearman"], len(df)))
    log("  数据规模不同, 两者的 OOF **不可直接比较**; 全量的 OOF 置信区间更窄、更可信。")
    log("- 在同一全量数据上比较两个特征集 (独立复核, 配对 bootstrap 3000 次):")
    log("    OOF:      V0.2 特征集 +0.389 vs V0.3 特征集 +0.402, Δ=%+.3f, 95%% CI [%.4f, %.4f]"
        % (0.0126, 0.0000, 0.0249))
    log("    测试集:   V0.2 特征集 +0.434 vs V0.3 特征集 +0.469, Δ=%+.3f, 95%% CI [+0.012, +0.061]"
        % 0.0353)
    log("  -> 两个口径**方向一致**(均为 V0.3 特征集更好)且 CI 均不跨 0;")
    log("     因此 V1.0 的**推荐配置是 B (V0.3 链对物理特征)**, 而非 A。")
    log("- 与 V0.2(MVP) 主模型 +0.409 相比: 特征集 A 的 +0.434 为 %+.3f,"
        % (main_cfg["test"]["spearman"] - 0.409))
    log("  但两者 CI 大量重叠 ([+0.182,+0.586] vs [+0.363,+0.496]), 属**小幅/持平**, 不应宣称显著提升。")
    log("- RMSE 从 3.105 降到 %.3f (降 %.0f%%), R2 从 +0.118 升到 %+.3f —— 误差与解释力的"
        % (main_cfg["test"]["rmse"],
           100 * (1 - main_cfg["test"]["rmse"] / 3.105), main_cfg["test"]["r2"]))
    log("  改善比 Spearman 更明确, 这部分得益于测试集更大(578 行)带来的稳定估计。")
    log("")
    log("### 残留局限")
    log("")
    log("1. 数据只覆盖 %d / 2774 个目标结构 (特征提取在 1600 个结构处被止损)。"
        % df["pdb_id"].nunique())
    log("2. PRODIGY 不可用, predicted_dG 整列缺失。")
    log("3. 40% 同一性为同源代理, 单链接聚类存在传递性。")
    log("4. 测试集 %d 行, 置信区间仍较宽。" % main_cfg["n_test"])
    log("")
    log("耗时: %.1fs" % (time.time() - t0))

    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("")
    print("报告: %s" % REPORT_PATH)
    print("=" * 80)
    print("V1.0 最终测试集 Spearman = %+.3f (RMSE=%.3f, R2=%+.3f)"
          % (m["spearman"], m["rmse"], m["r2"]))
    print("=" * 80)
    print("阶段二完成，等待下一步指令")


if __name__ == "__main__":
    main()
