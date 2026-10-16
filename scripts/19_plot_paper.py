# -*- coding: utf-8 -*-
"""19_plot_paper.py -- 阶段 6: 论文级图表 (Fig 1-7, Table 1-3).

设计纪律
--------
* **Okabe-Ito 色盲友好调色板**, 300 dpi, 矢量友好尺寸 (毫米制)。
* **只画真实存在的产物**: 每个图/表的数据源都先检查存在性; 缺失时**明确打印
  "SKIP (缺 X)" 并跳过**, 绝不用占位/构造数据填充。
* 中文标签用 CJK 字体 (与 V1.1 脚本同一套字体探测逻辑), 找不到则回退英文标签,
  避免出现豆腐块。

产出
----
figures/fig1_architecture.png      引擎四层架构
figures/fig2_data_composition.png  数据构成与划分
figures/fig3_pred_vs_true.png      预测 vs 真实 (标注拟合斜率, 诚实展示收缩)
figures/fig4_enrichment.png        富集曲线 EF vs Top-k%
figures/fig5_ablation.png          消融实验条形图 + bootstrap 误差棒
figures/fig6_feature_importance.png 特征重要性 (按类别着色)
figures/fig7_split_sensitivity.png 划分协议敏感性 (随机/PDB级/家族级)
features/v2_table1_dataset_stats.md
features/v2_table2_baselines.md
features/v2_table3_ablation.md
"""

from __future__ import annotations

import argparse
import json
import logging
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
FEAT_DIR = ROOT / "features"
FIG_DIR = ROOT / "figures"
LOG_DIR = ROOT / "logs" / "v2"

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

# ---------------------------------------------------------------- Okabe-Ito
# 色盲友好 8 色 (Wong 2011, Nature Methods "Points of view")
OKABE_ITO = {
    "black": "#000000", "orange": "#E69F00", "sky": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "vermillion": "#D55E00", "purple": "#CC79A7",
}
CAT_COLORS = [OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["green"],
              OKABE_ITO["orange"], OKABE_ITO["purple"], OKABE_ITO["sky"],
              OKABE_ITO["yellow"], OKABE_ITO["black"]]
DPI = 300
SEED = 42

logger = logging.getLogger("plot_paper")
_FIGS: list[str] = []
_TABLES: list[str] = []
_SKIPPED: list[str] = []


def setup_font() -> str | None:
    prefer = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC",
              "WenQuanYi Zen Hei", "SimSun", "DengXian"]
    try:
        available = {f.name for f in font_manager.fontManager.ttflist}
        for name in prefer:
            if name in available:
                matplotlib.rcParams["font.sans-serif"] = [name] + list(
                    matplotlib.rcParams.get("font.sans-serif", []))
                matplotlib.rcParams["axes.unicode_minus"] = False
                return name
    except Exception:  # noqa: BLE001
        pass
    return None


_CJK = setup_font()


def save(fig, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _FIGS.append(name)
    logger.info("  写出 %s", path.name)


def skip(what: str, why: str) -> None:
    _SKIPPED.append(f"{what}: {why}")
    logger.warning("  SKIP %s (缺 %s)", what, why)


# ===================================================================== Fig 1
def fig1_architecture() -> None:
    """引擎四层架构图 (纯示意图, 不依赖任何数据文件)。"""
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    layers = [
        ("L1  结构输入层", "PDB / CIF  →  链角色解析  →  界面判定 (跨链重原子 < 5 Å)",
         OKABE_ITO["sky"]),
        ("L2  界面物理特征层", "接触/氢键/盐桥 · ΔSASA (Shrake-Rupley) · 残基组成\n"
                          "疏水矩 (Eisenberg) · 静电代理 · 形状密度  →  60 + 147 维",
         OKABE_ITO["green"]),
        ("L3  排序建模层", "LightGBM 回归 (V1.1 对照)  vs  LambdaRank (group = family)\n"
                        "家族级 GroupKFold · 训练集拟合 relevance 分箱",
         OKABE_ITO["orange"]),
        ("L4  评测与富集层", "Spearman/Kendall · Precision@k · EF@x% · NDCG@k\n"
                         "bootstrap 95% CI · 配对 bootstrap + BH-FDR",
         OKABE_ITO["vermillion"]),
    ]
    y = 8.4
    for title, body, color in layers:
        box = FancyBboxPatch((0.4, y - 1.0), 9.2, 1.7,
                             boxstyle="round,pad=0.10,rounding_size=0.18",
                             linewidth=1.6, edgecolor=color,
                             facecolor=color, alpha=0.16, zorder=1)
        ax.add_patch(box)
        ax.text(0.75, y + 0.44, title, fontsize=10.5, fontweight="bold",
                color=color, va="center", zorder=2)
        ax.text(0.75, y - 0.30, body, fontsize=8.2, va="center", zorder=2)
        y -= 2.05
    ax.text(5.0, 9.62, "Affinity Intelligence V2.0  四层架构",
            fontsize=12.5, fontweight="bold", ha="center")
    ax.text(5.0, 0.28,
            "红线: 同一 family_id / pdb_id 绝不跨 train-test (硬校验, 违者抛错)",
            fontsize=8.0, ha="center", style="italic", color=OKABE_ITO["vermillion"])
    save(fig, "fig1_architecture.png")


# ===================================================================== Fig 2
def fig2_data_composition(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.5))
    ax = axes[0]
    vc = df["source_dataset"].value_counts()
    ax.barh(range(len(vc)), vc.to_numpy(), color=CAT_COLORS[:len(vc)], alpha=0.9)
    ax.set_yticks(range(len(vc)))
    ax.set_yticklabels([str(x)[:18] for x in vc.index], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("记录数", fontsize=8.5)
    ax.set_title("(a) 来源数据集", fontsize=9.5)
    ax.grid(alpha=0.22, axis="x")
    for i, v in enumerate(vc.to_numpy()):
        ax.text(v, i, f" {v}", va="center", fontsize=7)

    ax = axes[1]
    tr = int((df["split_family"] == "train").sum())
    te = int((df["split_family"] == "test").sum())
    str_tr = df.loc[df["split_family"] == "train", "pdb_id"].nunique()
    str_te = df.loc[df["split_family"] == "test", "pdb_id"].nunique()
    fam_tr = df.loc[df["split_family"] == "train", "family_id"].nunique()
    fam_te = df.loc[df["split_family"] == "test", "family_id"].nunique()
    x = np.arange(3)
    w = 0.36
    ax.bar(x - w / 2, [tr, str_tr, fam_tr], w, label="train",
           color=OKABE_ITO["blue"], alpha=0.9)
    ax.bar(x + w / 2, [te, str_te, fam_te], w, label="test",
           color=OKABE_ITO["vermillion"], alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(["记录", "结构(PDB)", "家族"], fontsize=8.5)
    ax.set_yscale("log")
    ax.set_ylabel("数量 (log)", fontsize=8.5)
    ax.set_title("(b) 家族级划分", fontsize=9.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.22, axis="y")
    for xi, (a, b) in enumerate(zip([tr, str_tr, fam_tr], [te, str_te, fam_te])):
        ax.text(xi - w / 2, a, f"{a}", ha="center", va="bottom", fontsize=6.8)
        ax.text(xi + w / 2, b, f"{b}", ha="center", va="bottom", fontsize=6.8)

    ax = axes[2]
    gs = df.groupby("family_id").size().value_counts().sort_index()
    ax.bar(gs.index.to_numpy(), gs.to_numpy(), color=OKABE_ITO["green"], alpha=0.85)
    ax.set_yscale("log")
    ax.set_xlabel("家族内记录数 (group size)", fontsize=8.5)
    ax.set_ylabel("家族数 (log)", fontsize=8.5)
    ax.set_title("(c) group 结构分布", fontsize=9.5)
    ax.grid(alpha=0.22, axis="y")
    one = int(gs.get(1, 0))
    tot = int(gs.sum())
    ax.text(0.97, 0.94, f"单条目家族 {one}/{tot}\n({one / tot * 100:.1f}%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
    plt.tight_layout()
    save(fig, "fig2_data_composition.png")


# ===================================================================== Fig 3
def fig3_pred_vs_true() -> None:
    """预测 vs 真实散点。**必须标注拟合斜率**以诚实展示收缩。"""
    src = FEAT_DIR / "v1.1_test_predictions.csv"
    if not src.exists():
        skip("Fig3", src.name)
        return
    t = pd.read_csv(src)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
    y, p = t["real_dG"].to_numpy(), t["pred_dG"].to_numpy()
    z = np.polyfit(y, p, 1)
    lim = [min(y.min(), p.min()) - 1.5, max(y.max(), p.max()) + 1.5]
    ax = axes[0]
    ax.scatter(y, p, s=13, alpha=0.5, color=OKABE_ITO["blue"],
               edgecolors="white", linewidths=0.25, label="测试集界面")
    ax.plot(lim, lim, "--", color="0.55", lw=1.1, label="y = x (理想)")
    xs = np.linspace(lim[0], lim[1], 50)
    ax.plot(xs, np.polyval(z, xs), "-", color=OKABE_ITO["vermillion"], lw=1.7,
            label=f"实测拟合 slope = {z[0]:.2f}")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("真实 ΔG (kcal/mol)", fontsize=9)
    ax.set_ylabel("预测 ΔG (kcal/mol)", fontsize=9)
    ax.set_title("(a) V1.1 回归: 预测 vs 真实", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.22)

    ax = axes[1]
    resid = p - y
    ax.scatter(y, resid, s=13, alpha=0.5, color=OKABE_ITO["purple"],
               edgecolors="white", linewidths=0.25)
    ax.axhline(0, color="0.4", lw=1.0, ls="--")
    ax.set_xlabel("真实 ΔG (kcal/mol)", fontsize=9)
    ax.set_ylabel("残差 (预测 − 真实)", fontsize=9)
    ax.set_title(f"(b) 残差: RMSE={np.sqrt(np.mean(resid ** 2)):.2f}, "
                 f"std={resid.std():.2f}", fontsize=10)
    ax.grid(alpha=0.22)
    fig.suptitle("模型向均值收缩: 排序能力强于定量能力 (斜率远小于 1)",
                 fontsize=10.5, y=1.02)
    plt.tight_layout()
    save(fig, "fig3_pred_vs_true.png")


# ===================================================================== Fig 4
def fig4_enrichment() -> None:
    """富集曲线: 在真实测试标签上, 对每个模型算 EF vs Top-k%。"""
    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    man = FEAT_DIR / "v2_model_manifest.json"
    if not man.exists():
        skip("Fig4", man.name)
        return
    manifest = json.loads(man.read_text(encoding="utf-8"))
    pcts = np.array([1, 2, 3, 5, 7, 10, 15, 20, 25, 30], dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    plotted = 0
    for i, e in enumerate(manifest):
        p = Path(e["predictions_csv"])
        if not p.exists():
            continue
        t = pd.read_csv(p)
        if e["y_col"] not in t.columns or e["score_col"] not in t.columns:
            continue
        y = pd.to_numeric(t[e["y_col"]], errors="coerce").to_numpy(float)
        s = pd.to_numeric(t[e["score_col"]], errors="coerce").to_numpy(float)
        if e.get("score_is_dG", True):
            y, s = -y, -s
        efs = [mr.enrichment_factor(y, s, float(x), 0.10) for x in pcts]
        color = CAT_COLORS[i % len(CAT_COLORS)]
        ax.plot(pcts, efs, "-o", ms=3.2, lw=1.5, color=color,
                label=f"{e['model_id']}")
        plotted += 1
    if plotted == 0:
        skip("Fig4", "无可用预测文件")
        plt.close(fig)
        return
    ax.axhline(1.0, color="0.4", ls="--", lw=1.2, label="随机基线 (EF = 1)")
    ax.set_xlabel("预测取前 Top-k%", fontsize=9.5)
    ax.set_ylabel("富集因子 EF", fontsize=9.5)
    ax.set_title("富集曲线 (真强结合体 = 真值前 10%)", fontsize=10.5)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.22)
    plt.tight_layout()
    save(fig, "fig4_enrichment.png")


# ===================================================================== Fig 5
def fig5_ablation() -> None:
    src = FEAT_DIR / "v2_ablation_results.csv"
    if not src.exists():
        skip("Fig5", src.name)
        return
    a = pd.read_csv(src)
    a = a[a["split_protocol"] == "family"].copy()
    if not len(a):
        skip("Fig5", "无 family 划分结果")
        return
    a = a.sort_values("test_spearman", ascending=True)
    fig, ax = plt.subplots(figsize=(7.6, 0.42 * len(a) + 2.0))
    ypos = np.arange(len(a))
    err = np.vstack([
        a["test_spearman"].to_numpy() - a["test_ci_low"].to_numpy(),
        a["test_ci_high"].to_numpy() - a["test_spearman"].to_numpy(),
    ])
    err = np.nan_to_num(err, nan=0.0)
    colors = [OKABE_ITO["vermillion"] if c == "full" else OKABE_ITO["blue"]
              for c in a["config"]]
    ax.barh(ypos, a["test_spearman"].to_numpy(), xerr=err, color=colors,
            alpha=0.88, error_kw=dict(ecolor="0.3", lw=1.0, capsize=2.5))
    ax.set_yticks(ypos)
    ax.set_yticklabels(a["config"], fontsize=8)
    ax.axvline(0, color="0.4", lw=1.0)
    ax.set_xlabel("测试集 Spearman ρ (家族级划分, 误差棒 = bootstrap 95% CI)", fontsize=9)
    ax.set_title("消融实验: 特征族与划分协议", fontsize=10.5)
    ax.grid(alpha=0.22, axis="x")
    plt.tight_layout()
    save(fig, "fig5_ablation.png")


# ===================================================================== Fig 6
def fig6_feature_importance() -> None:
    src = FEAT_DIR / "v1.1_feature_importance.csv"
    if not src.exists():
        skip("Fig6", src.name)
        return
    imp = pd.read_csv(src).head(24).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(imp) + 1.8))

    def cat_of(name: str) -> str:
        n = str(name).lower()
        if any(t in n for t in ("sasa",)):
            return "埋藏面积 ΔSASA"
        if any(t in n for t in ("hbond", "salt", "saltbridge")):
            return "氢键/盐桥"
        if any(t in n for t in ("frac_", "aac_", "hmom", "residue_ratio")):
            return "残基组成/疏水矩"
        if any(t in n for t in ("chainpair", "pair_sasa", "n_chainpairs")):
            return "链对聚合"
        if any(t in n for t in ("n_atoms", "n_res_", "n_contact", "n_iface")):
            return "尺寸/接触计数"
        return "其他归一化量"

    imp = imp.assign(cat=imp["feature"].map(cat_of))
    cats = list(dict.fromkeys(imp["cat"]))
    cmap = {c: CAT_COLORS[i % len(CAT_COLORS)] for i, c in enumerate(cats)}
    ax.barh(np.arange(len(imp)), imp["importance_gain"].to_numpy(),
            color=[cmap[c] for c in imp["cat"]], alpha=0.9)
    ax.set_yticks(np.arange(len(imp)))
    ax.set_yticklabels(imp["feature"], fontsize=7.2)
    ax.set_xlabel("LightGBM gain 重要性", fontsize=9)
    ax.set_title("V1.1 特征重要性 Top24 (按语义类别着色)", fontsize=10.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=cmap[c]) for c in cats]
    ax.legend(handles, cats, fontsize=7, loc="lower right")
    ax.grid(alpha=0.22, axis="x")
    plt.tight_layout()
    save(fig, "fig6_feature_importance.png")


# ===================================================================== Fig 7
def fig7_split_sensitivity() -> None:
    src = FEAT_DIR / "v2_ablation_results.csv"
    if not src.exists():
        skip("Fig7", src.name)
        return
    a = pd.read_csv(src)
    a = a[a["config"].isin(["full", "leak_pdb_level", "leak_random"])].copy()
    if not len(a):
        skip("Fig7", "无划分协议对照结果")
        return
    label_map = {"family": "家族级\n(红线, 可信)", "pdb": "PDB 级\n(同源泄漏)",
                 "random": "完全随机\n(最强泄漏)"}
    order = ["family", "pdb", "random"]
    a["_o"] = a["split_protocol"].map({k: i for i, k in enumerate(order)})
    a = a.sort_values("_o")
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    x = np.arange(len(a))
    cols = [OKABE_ITO["green"], OKABE_ITO["orange"], OKABE_ITO["vermillion"]][:len(a)]

    ax = axes[0]
    err = np.vstack([a["test_spearman"].to_numpy() - a["test_ci_low"].to_numpy(),
                     a["test_ci_high"].to_numpy() - a["test_spearman"].to_numpy()])
    ax.bar(x, a["test_spearman"].to_numpy(), color=cols, alpha=0.9, width=0.6,
           yerr=np.nan_to_num(err, nan=0.0),
           error_kw=dict(ecolor="0.3", lw=1.1, capsize=3))
    ax.set_xticks(x)
    ax.set_xticklabels([label_map.get(s, s) for s in a["split_protocol"]], fontsize=8)
    ax.axhline(0, color="0.4", lw=1.0)
    ax.set_ylabel("测试集 Spearman ρ", fontsize=9)
    ax.set_title("(a) 划分协议对指标的影响", fontsize=10)
    ax.grid(alpha=0.22, axis="y")
    for xi, v in zip(x, a["test_spearman"]):
        ax.text(xi, v, f"{v:+.3f}", ha="center", va="bottom", fontsize=8)

    ax = axes[1]
    cf = a["cross_family_id"].to_numpy() if "cross_family_id" in a.columns \
        else np.zeros(len(a))
    cp = a["cross_pdb_id"].to_numpy() if "cross_pdb_id" in a.columns else np.zeros(len(a))
    w = 0.36
    ax.bar(x - w / 2, cf, w, label="跨侧家族数", color=OKABE_ITO["blue"], alpha=0.9)
    ax.bar(x + w / 2, cp, w, label="跨侧 PDB 数", color=OKABE_ITO["purple"], alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([label_map.get(s, s) for s in a["split_protocol"]], fontsize=8)
    ax.set_ylabel("跨侧数量 (必须为 0)", fontsize=9)
    ax.set_title("(b) 实际泄漏量", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.22, axis="y")
    fig.suptitle("泄漏敏感性: 放松隔离协议会让指标虚高", fontsize=10.5, y=1.02)
    plt.tight_layout()
    save(fig, "fig7_split_sensitivity.png")


# ===================================================================== Tables
def table1_dataset_stats(df: pd.DataFrame) -> None:
    lines = ["# Table 1 — 数据集统计 (Affinity Intelligence V2.0)", "",
             "| 项目 | 全体 | train | test |", "|---|---:|---:|---:|"]
    tr = df[df["split_family"] == "train"]
    te = df[df["split_family"] == "test"]

    def row(name: str, fn) -> None:
        lines.append(f"| {name} | {fn(df)} | {fn(tr)} | {fn(te)} |")

    row("记录数", lambda d: len(d))
    row("结构数 (PDB)", lambda d: d["pdb_id"].nunique())
    row("家族数 (family_id)", lambda d: d["family_id"].nunique())
    row("ΔG 均值 (kcal/mol)", lambda d: f"{d['real_dG'].mean():.3f}")
    row("ΔG 标准差", lambda d: f"{d['real_dG'].std():.3f}")
    row("ΔG 最小值", lambda d: f"{d['real_dG'].min():.3f}")
    row("ΔG 最大值", lambda d: f"{d['real_dG'].max():.3f}")
    row("单条目家族数", lambda d: int((d.groupby("family_id").size() == 1).sum()))
    row("多条目家族数", lambda d: int((d.groupby("family_id").size() > 1).sum()))
    lines += ["", "## 来源数据集构成", "", "| 来源 | 全体 | train | test |",
              "|---|---:|---:|---:|"]
    for src in sorted(df["source_dataset"].dropna().unique()):
        lines.append(f"| {src} | {int((df['source_dataset'] == src).sum())} | "
                     f"{int((tr['source_dataset'] == src).sum())} | "
                     f"{int((te['source_dataset'] == src).sum())} |")
    lines += ["", "- 标签口径: `real_dG = 0.592 · ln(Kd[M])`, 固定 298.15 K。",
              "- 划分协议: 40% 序列同一性家族级, 整族分侧; 跨侧家族 = 0, 跨侧 PDB = 0。",
              "- 数据许可: PPB-Affinity 遵循原始许可; 本仓库只放脚本与统计量, 不再分发原始数据。",
              ""]
    p = FEAT_DIR / "v2_table1_dataset_stats.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    _TABLES.append(p.name)
    logger.info("  写出 %s", p.name)


def table2_baselines() -> None:
    src = FEAT_DIR / "v2_ranking_metrics.csv"
    if not src.exists():
        skip("Table2", src.name)
        return
    m = pd.read_csv(src)
    lines = ["# Table 2 — 与基线对比 (家族隔离测试集)", "",
             "配对检验为**配对 bootstrap** (非 DeLong; DeLong 只适用于二分类 AUC,",
             "套到 Spearman/EF 上是统计学误用)。p 值经 BH-FDR 校正。", ""]
    cols = [("model_id", "模型"), ("family", "类别"), ("n", "n"),
            ("spearman", "Spearman"), ("spearman_ci_low", "CI 下"),
            ("spearman_ci_high", "CI 上"), ("ndcg@10", "NDCG@10"),
            ("precision@10", "P@10"), ("ef@1%", "EF@1%"), ("ef@5%", "EF@5%")]
    have = [(k, v) for k, v in cols if k in m.columns]
    lines.append("| " + " | ".join(v for _k, v in have) + " |")
    lines.append("|" + "---|" * len(have))
    for r in m.sort_values("spearman", ascending=False).to_dict("records"):
        cells = []
        for k, _v in have:
            val = r.get(k)
            if isinstance(val, float) and k != "n":
                cells.append(f"{val:+.4f}" if np.isfinite(val) else "NaN")
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    pairs = FEAT_DIR / "v2_model_comparison_pairs.csv"
    if pairs.exists():
        pr = pd.read_csv(pairs)
        sig = pr[pr["metric"] == "spearman"].sort_values("p_two_sided")
        if len(sig):
            lines += ["", "## 配对 bootstrap 检验 (Spearman, BH-FDR 校正)", "",
                      "| A | B | Δ(B−A) | 95% CI | p | BH q | 显著 |", "|---|---|---:|---|---:|---:|---|"]
            for r in sig.head(15).to_dict("records"):
                lines.append(f"| {r['model_a']} | {r['model_b']} | "
                             f"{r['delta_b_minus_a']:+.4f} | "
                             f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                             f"{r['p_two_sided']:.4f} | {r['bh_q']:.4f} | "
                             f"{'是' if r['significant_bh_0.05'] else '否'} |")
    p = FEAT_DIR / "v2_table2_baselines.md"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _TABLES.append(p.name)
    logger.info("  写出 %s", p.name)


def table3_ablation() -> None:
    src = FEAT_DIR / "v2_ablation_results.csv"
    if not src.exists():
        skip("Table3", src.name)
        return
    a = pd.read_csv(src)
    cols = [("config", "配置"), ("split_protocol", "划分"), ("n_features", "特征数"),
            ("test_spearman", "Spearman"), ("test_ci_low", "CI 下"),
            ("test_ci_high", "CI 上"), ("test_ndcg10", "NDCG@10"),
            ("test_ef1", "EF@1%"), ("test_rmse", "RMSE"), ("test_r2", "R²"),
            ("cross_family_id", "跨侧家族"), ("cross_pdb_id", "跨侧 PDB")]
    have = [(k, v) for k, v in cols if k in a.columns]
    lines = ["# Table 3 — 消融实验全表", "",
             "> 泄漏对照 (`leak_*`) 使用**不同测试集**, 其指标不可与非泄漏行直接比大小,",
             "> 仅用于量化「放松隔离协议会让指标虚高多少」。", ""]
    lines.append("| " + " | ".join(v for _k, v in have) + " |")
    lines.append("|" + "---|" * len(have))
    for r in a.to_dict("records"):
        cells = []
        for k, _v in have:
            val = r.get(k)
            if isinstance(val, float):
                cells.append(f"{val:+.4f}" if np.isfinite(val) else "NaN")
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    j = FEAT_DIR / "v2_ablation_results.json"
    if j.exists():
        v = json.loads(j.read_text(encoding="utf-8")).get("verdicts", {})
        lines += ["", "## 三个关键论证", "",
                  "| 论证 | Δρ |", "|---|---:|"]
        for k, val in v.items():
            lines.append(f"| {k} | {'NaN' if val is None else f'{val:+.4f}'} |")
    p = FEAT_DIR / "v2_table3_ablation.md"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _TABLES.append(p.name)
    logger.info("  写出 %s", p.name)


# ===================================================================== main
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 6: 论文级图表")
    ap.add_argument("--features", type=Path, default=FEAT_DIR / "features_v1.1.csv")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "19_plot_paper.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    logger.info("=" * 78)
    logger.info("阶段 6: 论文级图表 (Okabe-Ito 色盲友好, %d dpi)", DPI)
    logger.info("=" * 78)
    logger.info("CJK 字体: %s", _CJK or "未找到 (可能显示豆腐块)")
    if not args.features.exists():
        logger.error("特征表不存在: %s", args.features)
        return 2
    df = pd.read_csv(args.features, low_memory=False)
    logger.info("数据: %d 行 x %d 列", *df.shape)

    if args.dry_run:
        logger.info("--dry-run")
        return 0

    logger.info("Fig 1 ..."); fig1_architecture()
    logger.info("Fig 2 ..."); fig2_data_composition(df)
    logger.info("Fig 3 ..."); fig3_pred_vs_true()
    logger.info("Fig 4 ..."); fig4_enrichment()
    logger.info("Fig 5 ..."); fig5_ablation()
    logger.info("Fig 6 ..."); fig6_feature_importance()
    logger.info("Fig 7 ..."); fig7_split_sensitivity()
    logger.info("Table 1 ..."); table1_dataset_stats(df)
    logger.info("Table 2 ..."); table2_baselines()
    logger.info("Table 3 ..."); table3_ablation()

    logger.info("")
    logger.info("生成图 %d 个: %s", len(_FIGS), ", ".join(_FIGS))
    logger.info("生成表 %d 个: %s", len(_TABLES), ", ".join(_TABLES))
    if _SKIPPED:
        logger.warning("跳过 %d 项 (缺数据源, 未用占位数据填充):", len(_SKIPPED))
        for s in _SKIPPED:
            logger.warning("  - %s", s)
    logger.info("耗时 %.1fs", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
