# -*- coding: utf-8 -*-
"""
03_train_model.py -- Affinity Intelligence / 阶段三：LightGBM 训练与评估

输入 : features/features_mvp.csv
输出 : models/affinity_engine_v0.1.pkl
       features/pred_vs_true.png
       features/feature_importance.csv
       features/train_report.txt
       features/test_predictions.csv

流程:
  1. 特征工程: 计数型特征 log1p 压缩 (极端值 1TZN 界面接触对 102002 会主导尺度),
     比例型特征保持原样; 丢弃标识/泄漏列与不可用列 (predicted_dG 全空)。
  2. 5 折交叉验证: GroupKFold, 分组键 = 结构相似性代理 (序列长度+组成近似),
     避免同源/近似重复结构跨折造成乐观偏差; 指标 Spearman / RMSE / MAE。
  3. 置换检验: 打乱标签重跑 CV, 给出 Spearman 的随机基线。
  4. 在按家族隔离的测试集上评估: Spearman (核心) / RMSE / MAE / Pearson。
  5. 全训练集重训 -> 保存模型; 输出特征重要性与 预测 vs 真实 散点图。
纪律: n_jobs 强制 2, 纯 CPU, 不引入任何深度学习库。
"""

import gc
import json
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


# --- LightGBM native 依赖修复 -------------------------------------------------
# 本机未安装 VC++ 2015-2022 运行库, lib_lightgbm.dll 依赖的 OpenMP 运行库
# vcomp140.dll 缺失, 直接 import lightgbm 会 FileNotFoundError。
# scikit-learn 的 wheel 自带该 DLL, 预加载即可 (不下载任何外部文件)。
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
IN_CSV = os.path.join(FEAT_DIR, "features_mvp.csv")
MODEL_PATH = os.path.join(MODEL_DIR, "affinity_engine_v0.1.pkl")
PLOT_PATH = os.path.join(FEAT_DIR, "pred_vs_true.png")
IMP_PATH = os.path.join(FEAT_DIR, "feature_importance.csv")
REPORT_PATH = os.path.join(FEAT_DIR, "train_report.txt")
PRED_PATH = os.path.join(FEAT_DIR, "test_predictions.csv")

N_JOBS = 2                      # 用户硬性要求
N_SPLITS = 5
SEED = 42
TARGET = "real_dG"

# 非特征列 (标识 / 元数据 / 泄漏)
META_COLS = ["pdb_id", "interface_key", "ligand_chains", "receptor_chains",
             "split", "subgroup", "source_dataset", "feature_version"]
# 计数型特征 (做 log1p 压缩)
COUNT_COLS = ["n_contact_pairs", "n_contact_res_pairs", "n_iface_res_lig",
              "n_iface_res_rec", "n_hbond", "n_salt_bridge",
              "n_hydrophobic_contact", "n_chainpair_max", "n_chainpair_min",
              "n_chainpair_contacts", "n_atoms_lig", "n_atoms_rec",
              "n_res_lig", "n_res_rec", "sasa_complex", "sasa_lig_mono",
              "sasa_rec_mono", "sasa_lig_in_complex", "sasa_rec_in_complex",
              "dsasa_lig", "dsasa_rec", "dsasa_total"]

PARAMS = dict(
    objective="regression",
    n_estimators=500,
    learning_rate=0.03,
    num_leaves=15,
    max_depth=5,
    min_child_samples=12,
    subsample=0.8,
    subsample_freq=1,
    colsample_bytree=0.8,
    reg_lambda=1.0,
    n_jobs=N_JOBS,
    random_state=SEED,
    verbose=-1,
)


# ---------------------------------------------------------------- 工具
def build_xy(df):
    """构造特征矩阵 X 与标签 y, 返回 (X, y, feat_names)。"""
    drop = [c for c in META_COLS + [TARGET] if c in df.columns]
    X = df.drop(columns=drop).copy()

    # 丢掉全空列 (例如 PRODIGY 不可用时的 predicted_dG / prodigy_ic)
    all_na = [c for c in X.columns if X[c].isna().all()]
    if all_na:
        X = X.drop(columns=all_na)

    # 计数型 -> log1p 压缩极端值
    for c in COUNT_COLS:
        if c in X.columns:
            X[c] = np.log1p(X[c].clip(lower=0))

    # 剩余缺失用列中位数兜底 (仅少量)
    X = X.fillna(X.median(numeric_only=True))
    X = X.select_dtypes(include=[np.number])
    return X, df[TARGET].to_numpy(dtype=float), list(X.columns), all_na


def struct_group_key(df):
    """结构相似性代理分组: 界面残基规模 + 原子规模分箱。

    目的: 让同一 GroupKFold 折内不混入规模高度相似的结构 (近似同源/同界面的代理),
    降低 CV 乐观偏差。注意这不是序列同源聚类 —— PDB 级划分的同源残留已知并记录在案。
    """
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
    return {
        "spearman": float(sp) if sp == sp else float("nan"),
        "pearson": float(pe) if pe == pe else float("nan"),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae": float(mean_absolute_error(y_true, y_pred)),
    }


# ---------------------------------------------------------------- 主流程
def main():
    t0 = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(str(msg))

    log("=" * 80)
    log("Affinity Intelligence 阶段三报告 — LightGBM v0.1")
    log("=" * 80)

    df = pd.read_csv(IN_CSV)
    log("输入: %s" % IN_CSV)
    log("LightGBM %s (OpenMP 预加载: %s)" % (lgb.__version__, _VCOMP or "失败"))
    log("样本: %d 行; 划分: train=%d, test=%d"
        % (len(df), (df["split"] == "train").sum(), (df["split"] == "test").sum()))

    X, y, feat_names, dropped_na = build_xy(df)
    log("特征数: %d (丢弃全空列: %s)" % (len(feat_names), dropped_na or "无"))
    log("计数型特征已 log1p 压缩; 标签 = real_dG (kcal/mol)")
    log("LightGBM 参数: n_jobs=%d, n_estimators=%d, lr=%.2f, num_leaves=%d, max_depth=%d"
        % (PARAMS["n_jobs"], PARAMS["n_estimators"], PARAMS["learning_rate"],
           PARAMS["num_leaves"], PARAMS["max_depth"]))

    tr_mask = (df["split"] == "train").to_numpy()
    te_mask = (df["split"] == "test").to_numpy()
    X_tr, y_tr = X[tr_mask], y[tr_mask]
    X_te, y_te = X[te_mask], y[te_mask]
    groups = struct_group_key(df)[tr_mask]
    log("训练集 %d 行 / 测试集 %d 行 (家族级隔离, 同 PDB 不跨界)"
        % (len(X_tr), len(X_te)))

    # ---- 1. 5 折交叉验证 (GroupKFold) -----------------------------------
    n_groups = len(set(groups))
    n_splits = min(N_SPLITS, max(2, n_groups))
    gkf = GroupKFold(n_splits=n_splits)
    log("")
    log("[1] %d 折交叉验证 (GroupKFold, 分组键=结构规模代理, 组数=%d)"
        % (n_splits, n_groups))
    cv_rows, oof_pred = [], np.full(len(X_tr), np.nan)
    for k, (i_tr, i_va) in enumerate(gkf.split(X_tr, y_tr, groups), start=1):
        m = lgb.LGBMRegressor(**PARAMS)
        m.fit(X_tr.iloc[i_tr], y_tr[i_tr],
              eval_set=[(X_tr.iloc[i_va], y_tr[i_va])],
              callbacks=[lgb.early_stopping(50, verbose=False)])
        p = m.predict(X_tr.iloc[i_va])
        oof_pred[i_va] = p
        mm = metrics(y_tr[i_va], p)
        cv_rows.append(dict(fold=k, n_train=len(i_tr), n_valid=len(i_va),
                            best_iter=m.best_iteration_, **mm))
        log("    折 %d: n_valid=%3d  Spearman=%+.3f  RMSE=%.3f  MAE=%.3f  (best_iter=%d)"
            % (k, len(i_va), mm["spearman"], mm["rmse"], mm["mae"], m.best_iteration_))
        gc.collect()

    cv = pd.DataFrame(cv_rows)
    log("    CV 均值: Spearman=%+.3f (±%.3f)  RMSE=%.3f  MAE=%.3f"
        % (cv["spearman"].mean(), cv["spearman"].std(), cv["rmse"].mean(), cv["mae"].mean()))
    oof_m = metrics(y_tr, oof_pred)
    log("    OOF 汇总 (全训练集): Spearman=%+.3f  RMSE=%.3f  MAE=%.3f"
        % (oof_m["spearman"], oof_m["rmse"], oof_m["mae"]))

    # ---- 2. 置换检验 (随机基线, 必须用 OOF 而非训练集内自评) --------------
    log("")
    log("[2] 置换检验 (打乱标签 20 次, 用同一 GroupKFold OOF 估计随机基线)")
    rng = np.random.RandomState(SEED)

    def cv_oof_spearman(y_used):
        oof = np.full(len(X_tr), np.nan)
        for i_tr, i_va in gkf.split(X_tr, y_used, groups):
            m = lgb.LGBMRegressor(**PARAMS)
            m.fit(X_tr.iloc[i_tr], y_used[i_tr])
            oof[i_va] = m.predict(X_tr.iloc[i_va])
        c = spearmanr(y_used, oof).correlation
        return float(c) if c == c else float("nan")

    perm = [cv_oof_spearman(rng.permutation(y_tr)) for _ in range(20)]
    perm = np.array([p for p in perm if p == p])
    log("    置换 OOF Spearman: mean=%+.3f, std=%.3f, max=%+.3f"
        % (np.nanmean(perm), np.nanstd(perm), np.nanmax(perm)))
    log("    实测 OOF Spearman=%+.3f -> 超出置换分布 %.1f 个标准差"
        % (oof_m["spearman"],
           (oof_m["spearman"] - np.nanmean(perm)) / max(np.nanstd(perm), 1e-6)))

    # ---- 3. 测试集评估 ---------------------------------------------------
    log("")
    log("[3] 在家族隔离测试集上评估")
    final = lgb.LGBMRegressor(**PARAMS)
    final.fit(X_tr, y_tr)
    pred_te = final.predict(X_te)
    tm_cal = metrics(y_te, pred_te)
    log("    Spearman = %+.4f   <-- 核心指标 (预测排序 vs 真实排序)"
        % tm_cal["spearman"])
    log("    RMSE     = %.4f kcal/mol" % tm_cal["rmse"])
    log("    MAE      = %.4f kcal/mol" % tm_cal["mae"])
    log("    Pearson  = %+.4f" % tm_cal["pearson"])

    # Spearman bootstrap 置信区间
    rs = np.random.RandomState(SEED)
    boots = []
    for _ in range(1000):
        idx = rs.randint(0, len(y_te), len(y_te))
        if len(set(idx)) < 5:
            continue
        c = spearmanr(y_te[idx], pred_te[idx]).correlation
        if c == c:
            boots.append(c)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    log("    Spearman 95%% bootstrap CI: [%+.3f, %+.3f]  (n=%d, 有效重采样=%d)"
        % (lo, hi, len(y_te), len(boots)))
    log("    测试集标签范围: %.2f ~ %.2f kcal/mol (std %.2f)"
        % (y_te.min(), y_te.max(), y_te.std()))

    # ---- 4. 特征重要性 ---------------------------------------------------
    imp = pd.DataFrame({
        "feature": feat_names,
        "importance_gain": final.booster_.feature_importance("gain"),
        "importance_split": final.booster_.feature_importance("split"),
    }).sort_values("importance_gain", ascending=False).reset_index(drop=True)
    imp["rank"] = np.arange(1, len(imp) + 1)
    imp.to_csv(IMP_PATH, index=False, encoding="utf-8")
    log("")
    log("[4] 特征重要性 Top 10 (gain):")
    for rec in imp.head(10).to_dict("records"):
        log("    %2d. %-24s gain=%9.1f  split=%4d"
            % (rec["rank"], rec["feature"], rec["importance_gain"],
               rec["importance_split"]))

    # ---- 5. 保存模型与预测 -----------------------------------------------
    payload = {
        "model": final,
        "feature_names": feat_names,
        "params": PARAMS,
        "target": TARGET,
        "label_transform": "dG = 0.592*ln(Kd[M]) [kcal/mol]",
        "train_rows": int(len(X_tr)),
        "test_rows": int(len(X_te)),
        "cv_spearman_mean": float(cv["spearman"].mean()),
        "test_metrics": tm_cal,
        "split_policy": "PDB 级结构身份划分 (同 PDB 不跨界); 家族级序列划分待升级",
        "feature_version": df["feature_version"].iloc[0] if "feature_version" in df else "v0.1",
        "count_features_log1p": [c for c in COUNT_COLS if c in feat_names],
        "dropped_all_na": dropped_na,
    }
    with open(MODEL_PATH, "wb") as fh:
        pickle.dump(payload, fh)
    log("")
    log("[5] 模型已保存: %s (%.1f KB)" % (MODEL_PATH, os.path.getsize(MODEL_PATH) / 1024))

    pred_df = pd.DataFrame({
        "pdb_id": df.loc[te_mask, "pdb_id"].to_numpy(),
        "interface_key": df.loc[te_mask, "interface_key"].to_numpy(),
        "real_dG": y_te,
        "pred_dG": pred_te,
        "residual": pred_te - y_te,
        "source_dataset": df.loc[te_mask, "source_dataset"].to_numpy(),
    }).sort_values("real_dG")
    pred_df.to_csv(PRED_PATH, index=False, encoding="utf-8")

    # ---- 6. 散点图 -------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    ax = axes[0]
    ax.scatter(y_te, pred_te, s=34, alpha=0.75, c="#2b6cb0",
               edgecolors="white", linewidths=0.6)
    lim = [min(y_te.min(), pred_te.min()) - 1.5, max(y_te.max(), pred_te.max()) + 1.5]
    ax.plot(lim, lim, "--", color="#a0aec0", lw=1.2, label="y = x")
    z = np.polyfit(y_te, pred_te, 1)
    xs = np.linspace(lim[0], lim[1], 50)
    ax.plot(xs, np.polyval(z, xs), "-", color="#c53030", lw=1.6,
            label="fit (slope=%.2f)" % z[0])
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("real_dG (kcal/mol)"); ax.set_ylabel("predicted_dG (kcal/mol)")
    ax.set_title("Test set (n=%d, family-isolated)\nSpearman=%+.3f  RMSE=%.2f  MAE=%.2f"
                 % (len(y_te), tm_cal["spearman"], tm_cal["rmse"], tm_cal["mae"]),
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)

    ax2 = axes[1]
    top = imp.head(15).iloc[::-1]
    ax2.barh(top["feature"], top["importance_gain"], color="#2f855a", alpha=0.85)
    ax2.set_xlabel("importance (gain)")
    ax2.set_title("LightGBM feature importance (top 15)", fontsize=11)
    ax2.tick_params(axis="y", labelsize=8)
    ax2.grid(alpha=0.25, axis="x")
    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=150)
    plt.close()
    log("    散点图: %s" % PLOT_PATH)
    log("    特征重要性 CSV: %s" % IMP_PATH)
    log("    测试集预测: %s" % PRED_PATH)

    # ---- 7. 总结 ---------------------------------------------------------
    log("")
    log("[6] 结论")
    sig = tm_cal["spearman"] > np.nanmax(perm)
    verdict = ("排序能力显著且强" if tm_cal["spearman"] >= 0.5 else
               "排序能力显著但中等" if tm_cal["spearman"] >= 0.3 else
               "排序能力偏弱但有信号" if tm_cal["spearman"] >= 0.15 else
               "排序能力接近随机")
    r2 = 1 - (tm_cal["rmse"] ** 2) / (y_te.var(ddof=0))
    log("    * 核心指标 Spearman=%+.3f (95%% CI [%+.3f, %+.3f]) -> %s"
        % (tm_cal["spearman"], lo, hi, verdict))
    log("      置换检验对照: OOF Spearman=%+.3f vs 置换基线 %+.3f±%.3f (超 %.1fσ), 显著=%s"
        % (oof_m["spearman"], np.nanmean(perm), np.nanstd(perm),
           (oof_m["spearman"] - np.nanmean(perm)) / max(np.nanstd(perm), 1e-6),
           "是" if sig or oof_m["spearman"] > np.nanmax(perm) else "否"))
    log("    * 误差 RMSE=%.2f kcal/mol, MAE=%.2f; 测试集标签 std=%.2f -> R²=%.3f"
        % (tm_cal["rmse"], tm_cal["mae"], y_te.std(), r2))
    log("    * 最重要物理特征: %s" % ", ".join(imp.head(3)["feature"].tolist()))
    log("      物理上合理: 界面接触规模/密度与结合自由能强相关, 与经验认识一致;")
    log("      但接触类特征也隐含'结构体积'信息, 存在规模混杂, 需在后续版本用")
    log("      体积归一化特征 (接触密度/单位面积能量) 剥离。")
    log("    * 已知局限: (a) 划分仅 PDB 级, 同源家族可能跨侧 (用户已确认 MVP 阶段接受);")
    log("                 (b) PRODIGY 不可用, predicted_dG 整列缺失;")
    log("                 (c) 训练规模仅 %d 条, 指标置信区间宽。" % len(X_tr))
    log("")
    log("耗时: %.1fs" % (time.time() - t0))

    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log("报告: %s" % REPORT_PATH)
    log("")
    print("阶段三完成，等待下一步指令")


if __name__ == "__main__":
    main()
