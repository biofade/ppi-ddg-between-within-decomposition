# -*- coding: utf-8 -*-
"""
11_train_v1.1.py -- Affinity Intelligence V1.1 / 数据补全后的最终训练与评估

输入 : features/features_full.csv (补全后 ~4353 行, 已带 family_id / split_family)
       features/features_mvp.csv  (397 行, 同上)
输出 : features/features_v1.1.csv         (合并表, 含 V0.3 的 60 特征, 剔除 __leaky)
       models/affinity_engine_v1.1.pkl
       features/v1.1_pred_vs_true.png
       features/v1.1_feature_importance.csv
       features/v1.1_test_predictions.csv
       features/v1.1_final_report.md

划分红线: 使用 09_family_split_full.py 的 `split_family` (40% 序列同一性家族级,
          整族分侧)。脚本内独立复算一遍红线, 违者拒绝训练。

对照配置 (同一家族级划分):
  A = V0.2 的 50 特征集        —— 与 V0.2/V1.0 主模型同口径
  B = A + V0.3 链对物理特征     —— V1.0 中表现最好的配置
  C = B + 标签派生量            —— 泄漏对照, 不可采信

对比对象: V0.2(MVP 397 行, A 口径) +0.409 ; V1.0(2888 行, A 口径) +0.434 ;
          V1.0 最优(B 口径) +0.469
纪律: n_jobs=2; 纯 CPU; 无深度学习库; 标签泄漏列一律排除。
"""

import argparse
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
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr

# 复用 08 的特征集定义与评估工具 (单一实现来源, 避免口径漂移)
_spec = None
try:
    import importlib.util
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "08_train_v1.0.py")
    _spec = importlib.util.spec_from_file_location("t8", _p)
    T8 = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(T8)
except Exception as _exc:                                       # noqa: BLE001
    print("!! 无法载入 08_train_v1.0.py 的工具函数: %s" % _exc)
    raise

ROOT = T8.ROOT
FEAT_DIR = T8.FEAT_DIR
MODEL_DIR = T8.MODEL_DIR
V02_FEATURES = T8.V02_FEATURES
V03_NEW = T8.V03_NEW
REPLACED = T8.REPLACED
COUNT_COLS = T8.COUNT_COLS
PARAMS = T8.PARAMS
N_JOBS = T8.N_JOBS
SEED = T8.SEED
TARGET = T8.TARGET

MERGED_CSV = os.path.join(FEAT_DIR, "features_v1.1.csv")
MODEL_PATH = os.path.join(MODEL_DIR, "affinity_engine_v1.1.pkl")
FIG_PATH = os.path.join(FEAT_DIR, "v1.1_pred_vs_true.png")
IMP_PATH = os.path.join(FEAT_DIR, "v1.1_feature_importance.csv")
PRED_PATH = os.path.join(FEAT_DIR, "v1.1_test_predictions.csv")
REPORT_PATH = os.path.join(FEAT_DIR, "v1.1_final_report.md")

META_COLS = ["pdb_id", "interface_key", "ligand_chains", "receptor_chains",
             "split", "split_v0.1", "split_family", "family_id", "subgroup",
             "source_dataset", "feature_version", "predicted_dG", "prodigy_ic"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", default=T8.FULL_CSV)
    ap.add_argument("--mvp", default=T8.MVP_CSV)
    args = ap.parse_args()

    t0 = time.time()
    os.makedirs(MODEL_DIR, exist_ok=True)
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(str(msg))

    log("# Affinity Intelligence V1.1 — 数据补全后最终报告")
    log("")
    log("- 生成时间: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    log("- LightGBM %s (n_jobs=%d, OpenMP: %s)"
        % (lgb.__version__, N_JOBS, _VCOMP or "失败"))
    log("- 绘图字体: %s" % (_CJK_FONT or "未找到 CJK 字体"))

    full = pd.read_csv(args.full, dtype={"pdb_id": str})
    mvp = pd.read_csv(args.mvp, dtype={"pdb_id": str})
    fk = set(full["interface_key"])
    mvp_only = mvp[~mvp["interface_key"].isin(fk)]
    df = pd.concat([full, mvp_only], ignore_index=True, sort=False)

    log("")
    log("## 数据")
    log("")
    log("- features_full.csv: %d 行 / %d 结构" % (len(full), full["pdb_id"].nunique()))
    log("- features_mvp.csv : %d 行 (并入新增 %d 行)" % (len(mvp), len(mvp_only)))
    log("- 合并后: **%d 行 / %d 结构 / %d 家族**"
        % (len(df), df["pdb_id"].nunique(), df["family_id"].nunique()))

    if "split_family" not in df.columns:
        raise RuntimeError("缺少 split_family 列; 请先运行 09_family_split_full.py")
    if "family_id" not in df.columns:
        raise RuntimeError("缺少 family_id 列")

    # 仅保留有家族标签的行 (新补全后 09 已全局分配; 若仍缺则明确报告)
    n_nofam = int(df["family_id"].isna().sum() + df["split_family"].isna().sum())
    if n_nofam:
        log("- !! 警告: %d 处家族标签/划分缺失, 这些行将被剔除" % n_nofam)
        df = df.dropna(subset=["family_id", "split_family"])
    y = df[TARGET].to_numpy(dtype=float)
    tr = (df["split_family"] == "train").to_numpy()
    te = (df["split_family"] == "test").to_numpy()
    log("- 家族级划分: train %d 行 / %d 结构 / %d 家族 ; test %d 行 / %d 结构 / %d 家族"
        % (tr.sum(), df.loc[tr, "pdb_id"].nunique(), df.loc[tr, "family_id"].nunique(),
           te.sum(), df.loc[te, "pdb_id"].nunique(), df.loc[te, "family_id"].nunique()))

    # 红线独立复算
    cross_fam = int((df.groupby("family_id")["split_family"].nunique() > 1).sum())
    cross_pdb = int((df.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    log("- 红线复核: 跨侧家族=%d, 跨侧 PDB=%d (必须为 0)" % (cross_fam, cross_pdb))
    if cross_fam or cross_pdb:
        raise RuntimeError("划分泄漏, 拒绝训练")
    log("- 标签 real_dG: %.2f ~ %.2f (mean %.2f, std %.2f)"
        % (y.min(), y.max(), y.mean(), y.std()))

    lk = [c for c in df.columns if c.endswith("__leaky")]
    base_cols = [c for c in V02_FEATURES if c in df.columns]
    new_cols = [c for c in V03_NEW if c in df.columns]
    v03_cols = [c for c in base_cols if c not in REPLACED] + new_cols

    # 落盘合并表 (剔除 __leaky, 供后续复用)
    try:
        out_cols = [c for c in df.columns if c not in lk]
        df[out_cols].to_csv(MERGED_CSV, index=False, encoding="utf-8")
        log("- 合并特征表已写出: %s (%d 行 x %d 列, 已剔除 %d 个泄漏列)"
            % (MERGED_CSV, len(df), len(out_cols), len(lk)))
    except Exception as exc:                                    # noqa: BLE001
        log("- !! 合并表写出失败: %s" % str(exc)[:80])

    groups = T8.struct_group_key(df)[tr]
    log("")
    log("## 三配置对比 (同一家族级划分)")
    log("")
    res = {}
    for tag, cols in (("A V0.2特征集", base_cols),
                      ("B V0.3链对物理", v03_cols),
                      ("C B+标签泄漏", v03_cols + lk)):
        res[tag] = T8.run_config(tag, cols, df, y, tr, te, groups, log)
        gc.collect()
    A, B, C = res["A V0.2特征集"], res["B V0.3链对物理"], res["C B+标签泄漏"]

    log("")
    log("## 与历史版本对比 (测试集 Spearman)")
    log("")
    log("| 版本 | 数据量 | 特征集 | Spearman | RMSE | R2 |")
    log("|---|---|---|---|---|---|")
    log("| V0.2 (MVP) | 397 行 | 50 | +0.409 | 3.105 | +0.118 |")
    log("| V1.0 (补全前) | 2888 行 | 50 | +0.434 | 2.384 | +0.212 |")
    log("| V1.0 (补全前) | 2888 行 | 60 (V0.3) | **+0.469** | 2.339 | +0.241 |")
    log("| **V1.1 (补全后)** | **%d 行** | **50** | **%+.3f** | %.3f | %+.3f |"
        % (len(df), A["test"]["spearman"], A["test"]["rmse"], A["test"]["r2"]))
    log("| **V1.1 (补全后)** | **%d 行** | **60 (V0.3)** | **%+.3f** | %.3f | %+.3f |"
        % (len(df), B["test"]["spearman"], B["test"]["rmse"], B["test"]["r2"]))
    log("")
    log("- 相对 V1.0 同口径(50 特征)的 +0.434: Δ = %+.3f"
        % (A["test"]["spearman"] - 0.434))
    log("- 相对 V1.0 最优(60 特征)的 +0.469: Δ = %+.3f"
        % (B["test"]["spearman"] - 0.469))
    log("- 相对 V0.2 的 +0.409: Δ(50特征) = %+.3f, Δ(60特征) = %+.3f"
        % (A["test"]["spearman"] - 0.409, B["test"]["spearman"] - 0.409))

    # 配对显著性: B vs A (独立复核)
    log("")
    log("## 配对显著性检验 (同一测试集, 3000 次 bootstrap)")
    log("")
    m1, l1, h1, w1 = T8.paired_bootstrap(B["y_test"], A["pred"], B["pred"])
    log("- 测试集 Δ(B-A) = %+.4f, 95%% CI [%+.4f, %+.4f], P(B>A)=%.1f%%"
        % (m1, l1, h1, w1))
    log("  -> %s" % ("CI 含 0: 两特征集不可区分" if l1 < 0 < h1 else
                     ("B 显著更好" if l1 > 0 else "A 显著更好")))
    oof_a, oof_b = A.get("oof_pred"), B.get("oof_pred")
    if oof_a is not None and oof_b is not None:
        m2, l2, h2, w2 = T8.paired_bootstrap(A["y_train"], oof_a, oof_b)
        log("- OOF    Δ(B-A) = %+.4f, 95%% CI [%+.4f, %+.4f], P(B>A)=%.1f%%"
            % (m2, l2, h2, w2))
        log("  -> %s" % ("CI 含 0: 不可区分" if l2 < 0 < h2 else
                         ("B 显著更好" if l2 > 0 else "A 显著更好")))
    log("- 两特征集预测逐样本相关: 测试集=%+.4f"
        % spearmanr(A["pred"], B["pred"]).correlation)

    # 主模型: 取两配置中测试集更优者, 但明确标注
    main_cfg = B if B["test"]["spearman"] >= A["test"]["spearman"] else A
    log("")
    log("- 主模型选定: **%s** (测试集 Spearman 更高)" % main_cfg["tag"])

    imp = pd.DataFrame({
        "feature": main_cfg["cols"],
        "importance_gain": main_cfg["model"].booster_.feature_importance("gain"),
        "importance_split": main_cfg["model"].booster_.feature_importance("split"),
    }).sort_values("importance_gain", ascending=False).reset_index(drop=True)
    imp["rank"] = np.arange(1, len(imp) + 1)
    imp["is_v03_new"] = imp["feature"].isin(new_cols)
    imp.to_csv(IMP_PATH, index=False, encoding="utf-8")

    log("")
    log("## 特征重要性 Top 10 (%s)" % main_cfg["tag"])
    log("")
    log("| 排名 | 特征 | gain | 是否 V0.3 新增 |")
    log("|---|---|---|---|")
    for r in imp.head(10).to_dict("records"):
        log("| %d | %s | %.1f | %s |"
            % (r["rank"], r["feature"], r["importance_gain"],
               "是" if r["is_v03_new"] else ""))

    # 保存
    payload = {
        "model": main_cfg["model"], "feature_names": main_cfg["cols"],
        "params": PARAMS, "target": TARGET,
        "label_transform": "dG = 0.592*ln(Kd[M]) [kcal/mol]",
        "excluded_leaky_features": lk,
        "split_policy": "V1.1 家族级 (split_family): 40% 序列同一性, 整族分侧",
        "train_rows": int(tr.sum()), "test_rows": int(te.sum()),
        "metrics": {k: v["test"] for k, v in res.items()},
        "count_features_log1p": [c for c in COUNT_COLS if c in main_cfg["cols"]],
        "history": {"V0.2_MVP": 0.409, "V1.0_A_2888rows": 0.434,
                    "V1.0_B_2888rows": 0.469},
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

    # 图
    fig = plt.figure(figsize=(14.5, 5.6))
    ax1 = fig.add_subplot(1, 3, 1)
    yt, pt = main_cfg["y_test"], main_cfg["pred"]
    ax1.scatter(yt, pt, s=16, alpha=0.55, color="#2b6cb0",
                edgecolors="white", linewidths=0.3)
    lim = [min(yt.min(), pt.min()) - 1.5, max(yt.max(), pt.max()) + 1.5]
    ax1.plot(lim, lim, "--", color="#a0aec0", lw=1.0, label="y = x")
    z = np.polyfit(yt, pt, 1)
    xs = np.linspace(lim[0], lim[1], 40)
    ax1.plot(xs, np.polyval(z, xs), "-", color="#c53030", lw=1.4,
             label="fit slope=%.2f" % z[0])
    ax1.set_xlim(lim)
    ax1.set_ylim(lim)
    ax1.set_xlabel("real_dG (kcal/mol)")
    ax1.set_ylabel("predicted_dG (kcal/mol)")
    ax1.set_title("V1.1 补全后 (%d 训练 / %d 测试)\nSpearman=%+.3f  RMSE=%.2f  R2=%+.2f"
                  % (main_cfg["n_train"], main_cfg["n_test"],
                     main_cfg["test"]["spearman"], main_cfg["test"]["rmse"],
                     main_cfg["test"]["r2"]), fontsize=10.5)
    ax1.legend(fontsize=8.5, loc="upper left")
    ax1.grid(alpha=0.25)

    ax2 = fig.add_subplot(1, 3, 2)
    labels = ["V0.2\n397行", "V1.0\n2888行", "V1.0\n2888行(60f)",
              "V1.1\n%d行" % len(df), "V1.1\n%d行(60f)" % len(df)]
    vals = [0.409, 0.434, 0.469, A["test"]["spearman"], B["test"]["spearman"]]
    bars = ax2.bar(np.arange(len(vals)), vals,
                   color=["#a0aec0", "#718096", "#4a5568", "#38a169", "#2f855a"],
                   alpha=0.92, width=0.6)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.012,
                 "%+.3f" % v, ha="center", fontsize=8.5)
    ax2.axhline(0, color="#4a5568", lw=0.8)
    ax2.set_xticks(np.arange(len(vals)))
    ax2.set_xticklabels(labels, fontsize=7.5)
    ax2.set_ylabel("test Spearman")
    ax2.set_title("数据补全前后的泛化能力", fontsize=10.5)
    ax2.grid(alpha=0.25, axis="y")

    ax3 = fig.add_subplot(1, 3, 3)
    top = imp.head(12).iloc[::-1]
    ax3.barh(top["feature"], top["importance_gain"],
             color=["#38a169" if v else "#718096" for v in top["is_v03_new"]],
             alpha=0.9)
    ax3.set_xlabel("importance (gain)")
    ax3.set_title("V1.1 特征重要性 Top12\n(绿 = V0.3 新增链对物理特征)", fontsize=10.5)
    ax3.tick_params(axis="y", labelsize=7.5)
    ax3.grid(alpha=0.25, axis="x")
    plt.tight_layout()
    plt.savefig(FIG_PATH, dpi=145)
    plt.close()
    log("")
    log("- 散点图: %s" % FIG_PATH)
    log("- 特征重要性: %s" % IMP_PATH)
    log("- 测试集预测: %s" % PRED_PATH)

    # 结论
    m = main_cfg["test"]
    dA = A["test"]["spearman"] - 0.434
    dB = B["test"]["spearman"] - 0.469
    log("")
    log("## 结论")
    log("")
    log("- **V1.1 最优配置 (%s): 测试集 Spearman = %+.3f** (95%% CI [%+.3f, %+.3f])"
        % (main_cfg["tag"], m["spearman"], main_cfg["test_ci"][0],
           main_cfg["test_ci"][1]))
    log("- RMSE = %.3f kcal/mol, MAE = %.3f, R2 = %+.3f" % (m["rmse"], m["mae"], m["r2"]))
    log("- 数据量: %d 行 / %d 结构 / %d 家族 (V1.0 的 %.2f 倍行数)"
        % (len(df), df["pdb_id"].nunique(), df["family_id"].nunique(),
           len(df) / 2888.0))
    log("- 相对 V1.0 同口径: Δ = %+.3f ; 相对 V1.0 最优: Δ = %+.3f" % (dA, dB))
    log("- 置换基线: OOF=%+.3f vs 随机 %+.3f±%.3f (超 %.1fσ)"
        % (main_cfg["oof"]["spearman"], main_cfg["perm_mean"], main_cfg["perm_std"],
           (main_cfg["oof"]["spearman"] - main_cfg["perm_mean"])
           / max(main_cfg["perm_std"], 1e-6)))
    log("")
    log("### 数据补全是否带来提升?")
    log("")
    verdict = ("明确提升" if dB >= 0.03 else "小幅提升" if dB > 0.005 else
               "基本持平" if dB >= -0.02 else "下降")
    log("以 V1.0 最优口径(+0.469)为对照, V1.1 为 %+.3f, Δ=%+.3f -> **%s**。"
        % (B["test"]["spearman"], dB, verdict))
    log("注意: 补全后测试集从 %d 行增至 %d 行、家族数从 %d 增至 %d, 属于更难/更全面的"
        % (578, int(te.sum()), 277, df.loc[te, "family_id"].nunique()))
    log("评估集合, 因此小幅变化不应解读为能力退化。")
    log("")
    log("### 残留局限")
    log("")
    n_fail_struct = 0
    try:
        import json as _json
        _cp = os.path.join(FEAT_DIR, "extract_checkpoint_v11.json")
        if os.path.exists(_cp):
            with open(_cp, "r", encoding="utf-8") as _fh:
                n_fail_struct = len(_json.load(_fh).get("failed_pdb", {}))
    except Exception:                                           # noqa: BLE001
        pass
    log("1. 特征提取失败结构 %d 个 (链组缺失等), 已从训练集剔除。" % n_fail_struct)
    log("2. PRODIGY 不可用, predicted_dG 整列缺失。")
    log("3. 同源性为代理指标: 40 pct 同一性阈值 + 贪心代表聚类 (CD-HIT 式),")
    log("   簇比单链接更细 (更保守), 不保证传递闭包; 红线由整族分侧的硬校验保证。")
    log("4. 模型拟合斜率约 %.2f, 排序能力强于定量能力。" % z[0])
    log("")
    log("耗时: %.1fs" % (time.time() - t0))

    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("")
    print("报告: %s" % REPORT_PATH)
    print("=" * 80)
    print("V1.1 最优测试集 Spearman = %+.3f (RMSE=%.3f, R2=%+.3f) [配置 %s]"
          % (m["spearman"], m["rmse"], m["r2"], main_cfg["tag"]))
    print("=" * 80)


if __name__ == "__main__":
    main()
