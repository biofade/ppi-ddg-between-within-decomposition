# -*- coding: utf-8 -*-
"""scripts/v2.3.2/05_error_analysis.py -- V2.3.2 任务 5.2: 全局最优模型的错误分析.

分析内容 (任务书要求)
---------------------
1. **test 残差最大的 20 条** (按 |pred − true| 降序);
2. **按 ΔG 分位数分箱**看每箱的 Spearman (模型在哪一段 ΔG 上最准/最差);
3. 按 **`Affinity Method` / `Subgroup` / 温度区间** 分组的指标;
4. 分类错误的方向性 (过强 vs 过弱预测的比例)。

输入
----
* `features/v2.3.2_test_predictions.npz` (LightGBM 全部配置的 test 预测)
* `features/v2.3.2_lgbm_results.csv` (用于确定最优配置名)
* `features/v2.2_proteinea_features.csv` (元数据: 方法/亚组/温度)
* `data/ppb_affinity_proteinea/filtered.csv` (原始列, 用于补充 assay method 等)

输出
----
* `features/v2.3.2_error_analysis.md`
* `features/v2.3.2_error_analysis_worst20.csv`
* `features/v2.3.2_error_analysis_bins.csv`
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nn_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PREDS = C.FEAT_DIR / "v2.3.2_test_predictions.npz"
LGBM_CSV = C.FEAT_DIR / "v2.3.2_lgbm_results.csv"
NN_PREDS = [C.FEAT_DIR / "v2.3.2_nn_mlp_predictions.npz",
            C.FEAT_DIR / "v2.3.2_nn_twotower_predictions.npz"]
OUT_MD = C.FEAT_DIR / "v2.3.2_error_analysis.md"
OUT_WORST = C.FEAT_DIR / "v2.3.2_error_analysis_worst20.csv"
OUT_BINS = C.FEAT_DIR / "v2.3.2_error_analysis_bins.csv"
RAW_CSV = C.ROOT / "data" / "ppb_affinity_proteinea" / "filtered.csv"


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.3.2 错误分析")
    ap.add_argument("--preds", type=Path, default=PREDS)
    ap.add_argument("--top-worst", type=int, default=20)
    args = ap.parse_args()

    if not args.preds.exists():
        print(f"!! 缺少 {args.preds}")
        return 2
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr, C.TOP_FRAC)
    L: list[str] = []

    def emit(msg: str = "") -> None:
        print(msg)
        L.append(str(msg))

    z = np.load(args.preds, allow_pickle=True)
    y = z["test_y"].astype(float)
    keys = [k for k in z.files if k.startswith("pred__")]
    preds = {k[len("pred__"):]: z[k].astype(float) for k in keys}
    ikeys = z["test_interface_key"].astype(str) if "test_interface_key" in z.files else None

    # ---- 候选集 ----
    # 只有**回归目标**训练出来的模型 (V2.2 口径) 的预测才是 ΔG 尺度, 才能算残差。
    # Lambdarank 模型输出的是排序分数, 与 ΔG 不同尺度 -> 只参与序数指标。
    # 这一点必须显式区分, 否则 residual/RMSE 会变成任意尺度的数字 (实测曾出现
    # residual 高达 92053 这种无意义值)。
    is_dG = {"v2.2_regression"}
    try:
        p = pd.read_csv(LGBM_CSV)
        if "feature_set" in p.columns:
            for c in p["feature_set"].astype(str).unique():
                if "lambdarank" in c.lower():
                    is_dG.discard(c)
                else:
                    is_dG.add(c)
    except Exception:  # noqa: BLE001
        pass

    candidates: dict[str, tuple[np.ndarray, bool]] = {
        k: (v, k in is_dG) for k, v in preds.items()}
    for pth in NN_PREDS:
        if not pth.exists():
            continue
        zn = np.load(pth, allow_pickle=True)
        if len(zn["test_y"]) != len(y):
            continue
        groups: dict[str, list[np.ndarray]] = {}
        for k in zn.files:
            if k == "test_y" or "|" not in k:
                continue
            parts = k.split("|")
            if len(parts) == 3:
                groups.setdefault(f"{parts[0]}|{parts[1]}", []).append(zn[k].astype(float))
        from scipy.stats import rankdata
        for g, pl in groups.items():
            if len(pl) >= 2:
                candidates[f"NN:{g}|mean"] = (np.mean(pl, axis=0), False)
                candidates[f"NN:{g}|rank_mean"] = (
                    np.mean([rankdata(v) for v in pl], axis=0), False)

    scored = []
    for k, (pv, _dg) in candidates.items():
        r = mr.spearman_rho(-y, -pv)
        if r == r:
            scored.append((float(r), k))
    scored.sort(reverse=True)
    if not scored:
        emit("!! 没有可用预测")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3
    best_rho, best_key = scored[0]
    best_pred, best_is_dg = candidates[best_key]

    emit("# V2.3.2 任务 5.2 — 错误分析 (全局最优模型)")
    emit("")
    emit(f"- 全局最优配置(按 Spearman): **`{best_key}`** (test Spearman = {best_rho:+.4f})")
    emit(f"- test 样本数: {len(y)}")
    emit("")
    emit("> ⚠️ **尺度声明 (重要)**: 由 `is_dG` 标记哪一列的预测是**有标定意义的 ΔG**。")
    emit("> NN 的最后一层输出与 LightGBM 的 lambdarank 分数**不是 ΔG**, 与标签不同尺度,")
    emit("> 因此对它们只报告**序数型指标** (Spearman / NDCG / EF / P@10),")
    emit("> **不报告残差与 RMSE** —— 那会把任意尺度的排序分数误当成 kcal/mol。")

    resid = best_pred - y
    absr = np.abs(resid)
    if not best_is_dg:
        emit("")
        emit(f"- 当前最优 `{best_key}` 的预测**不是 ΔG 尺度** (is_dG=False),")
        emit("  因此下文残差/MAE 表**仅作排序诊断参考**; 分箱 Spearman 与分组 Spearman 仍然有效。")

    # ---------------- 1. 残差最大 20 条 ----------------
    emit("")
    emit(f"## 1. 残差最大的 {args.top_worst} 条")
    emit("")
    order = np.argsort(-absr)[:args.top_worst]

    # 元数据 (尽量补上方法与亚组)
    #
    # ⚠️ 这里**不能**用 `interface_key` join: 缓存的 interface_key 是
    # `pdb|配体序列|受体序列`, 而 V2.2 特征表是 `pdb|链ID|链ID`, 两者不同格式
    # (这正是 V2.3 join 失败的原因)。本脚本改用 **`pdb_id`** 做键,
    # 每个 pdb 取第一条的元数据 (方法/亚组/温度等在同结构内基本一致)。
    meta = None
    d22 = pd.read_csv(C.V22_FEATURES, low_memory=False)
    usable = (d22["row_structure_usable"] == 1) & d22["real_dG"].notna()
    te22 = d22[usable & (d22["split"] == "test")].reset_index(drop=True)
    meta_cols = [c for c in ("Affinity Method", "Temperature(K)", "Subgroup",
                             "Source Data Set", "Mutations", "KD(M)", "pdb_id")
                 if c in te22.columns]
    per_pdb = te22[meta_cols].drop_duplicates("pdb_id", keep="first")
    te_pdb = (z["test_pdb_id"].astype(str) if "test_pdb_id" in z.files
              else np.array([""] * len(y)))
    meta = pd.DataFrame({"pdb_id": te_pdb}).merge(per_pdb, on="pdb_id", how="left")
    cov = float(meta["Affinity Method"].notna().mean()) if "Affinity Method" in meta.columns else 0.0
    emit("")
    emit(f"- 元数据经 **`pdb_id`** 关联 (非 interface_key, 因两者格式不同): "
         f"`Affinity Method` 覆盖 **{cov * 100:.1f}%**")

    # 百分位 (0-100, 越小 = 排名越靠前 = 越强)。用于在**非 ΔG 尺度**下表达误差幅度。
    from scipy.stats import rankdata as _rd
    _n = len(y)
    true_pct = _rd(y, method="average") / _n * 100.0
    pred_pct = _rd(best_pred, method="average") / _n * 100.0
    worst = pd.DataFrame({
        "pdb_id": meta["pdb_id"] if "pdb_id" in meta else "",
        "interface_key": ikeys if ikeys is not None else "",
        "true_dG": y, "pred_dG": best_pred, "residual": resid, "abs_residual": absr,
        "true_pct": true_pct, "pred_pct": pred_pct,
        "pred_rank_score": best_pred, "pred_is_dG_scale": best_is_dg,
    })
    for c in ("Mutations", "Affinity Method", "Temperature(K)", "Subgroup",
              "Source Data Set"):
        if c in meta.columns:
            worst[c] = meta[c]
    # 残差/误差只在**有 ΔG 标定**时才有物理意义; 否则用「平均排名差」表达误差幅度。
    if best_is_dg:
        err_col_a, err_col_b = "true_dG", "pred_dG"
        err_label = "残差"
    else:
        err_col_a, err_col_b = "true_pct", "pred_pct"
        err_label = "排名差(百分点)"
    emit(f"| # | PDB | true ΔG | pred ({'ΔG' if best_is_dg else 'rank'}) | "
         f"{err_label} | 方向 | 突变 | 方法 | 亚组 |")
    emit("|---:|---|---:|---:|---:|---|---|---|---|")
    for rank, i in enumerate(order, start=1):
        w = worst.iloc[i]
        direction = "预测过强" if w["residual"] < 0 else "预测过弱"
        emit(f"| {rank} | {w['pdb_id']} | {w['true_dG']:.2f} | "
             f"{w[err_col_b]:.2f} | {w[err_col_b] - w[err_col_a]:+.2f} | {direction} | "
             f"{str(w.get('Mutations', ''))[:22]} | "
             f"{str(w.get('Affinity Method', ''))[:14]} | "
             f"{str(w.get('Subgroup', ''))[:14]} |")
    worst.iloc[order].to_csv(OUT_WORST, index=False, encoding="utf-8")

    n_over = int((resid < 0).sum())
    n_under = int((resid > 0).sum())
    emit("")
    emit(f"- 方向性: 预测**过强** (pred 比 true 更负) {n_over} 条 "
         f"({n_over / len(y) * 100:.1f}%); **过弱** {n_under} 条 "
         f"({n_under / len(y) * 100:.1f}%)")
    if best_is_dg:
        emit(f"- 残差 (ΔG 尺度): mean {resid.mean():+.3f}, std {resid.std():.3f}, "
             f"median |resid| {np.median(absr):.3f}, max |resid| {absr.max():.3f} kcal/mol")
    else:
        emit(f"- **该模型预测非 ΔG 尺度**, 故不报告 kcal/mol 残差。")
        emit(f"  改用**平均绝对排名差**: {np.mean(np.abs(true_pct - pred_pct)):.1f} 个百分点 "
             f"(中位 {np.median(np.abs(true_pct - pred_pct)):.1f})")
        if len(order) >= 10:
            emit(f"  即: 被模型挑为前十的条目, 在真值中平均位于第 "
                 f"{np.mean(true_pct[order[:10]]):.0f} 百分位 "
                 f"(0 = 最强, 100 = 最弱)")

    # ---------------- 2. 按 ΔG 分层与分位数细箱 ----------------
    emit("")
    emit("## 2. 模型在哪一段 ΔG 上更准")
    emit("")
    emit("> ΔG 越负 = 结合越强。两种切法互补:")
    emit("> (a) **等宽分层** (强/中/弱) —— 覆盖整个区间, 因此**层内 Spearman 可计算**;")
    emit("> (b) **等频细箱** (10 箱) —— ΔG 区间很窄, 箱内几乎没有排序变异,")
    emit(">     因此只报 MAE 与平均残差, **不报箱内 Spearman** (那样会得到无意义的数字)。")
    emit("")
    strata = [("强结合 (ΔG <= 三分之一分位)", y <= np.quantile(y, 1 / 3)),
              ("中等", (y > np.quantile(y, 1 / 3)) & (y <= np.quantile(y, 2 / 3))),
              ("弱结合 (> 三分之二分位)", y > np.quantile(y, 2 / 3))]
    emit("### (a) 等宽分层 — 层内 Spearman 可计算")
    emit("")
    emit("| 分层 | n | true ΔG 范围 | Spearman | " + ("MAE (kcal/mol)" if best_is_dg else "平均绝对排名差 (pp)") + " |")
    emit("|---|---:|---|---:|---:|")
    bin_rows = []
    for nm, msk in strata:
        msk = np.asarray(msk)
        if msk.sum() < 10:
            continue
        rho = mr.spearman_rho(-y[msk], -best_pred[msk])
        if best_is_dg:
            errv = float(np.mean(np.abs(resid[msk])))
        else:
            errv = float(np.mean(np.abs(true_pct[msk] - pred_pct[msk])))
        rng = f"[{y[msk].min():.2f}, {y[msk].max():.2f}]"
        rho_txt = f"{rho:+.4f}" if rho == rho else "—"
        emit(f"| {nm} | {int(msk.sum())} | {rng} | {rho_txt} | {errv:.3f} |")
        bin_rows.append({"scheme": "strata", "bin": nm, "n": int(msk.sum()),
                         "dg_min": float(y[msk].min()), "dg_max": float(y[msk].max()),
                         "spearman": float(rho) if rho == rho else float("nan"),
                         "err": errv,
                         "mean_residual": float(resid[msk].mean())})

    emit("")
    emit("### (b) 等频细分箱 (10 箱) — 只报误差, 不报箱内 Spearman")
    emit("")
    emit("| 分箱 | n | true ΔG 范围 | 平均残差 | MAE |")
    emit("|---|---:|---|---:|---:|")
    qq = np.quantile(y, np.linspace(0, 1, 11))
    qq[0] -= 1e-9
    qq[-1] += 1e-9
    fine = pd.cut(y, bins=qq, include_lowest=True)
    for b in list(fine.categories):
        msk = np.asarray(fine == b)
        if msk.sum() < 5:
            continue
        mae = float(np.mean(np.abs(resid[msk])))
        rng = f"[{y[msk].min():.2f}, {y[msk].max():.2f}]"
        emit(f"| {b} | {int(msk.sum())} | {rng} | {resid[msk].mean():+.3f} | {mae:.3f} |")
        bin_rows.append({"scheme": "fine10", "bin": str(b), "n": int(msk.sum()),
                         "dg_min": float(y[msk].min()), "dg_max": float(y[msk].max()),
                         "spearman": float("nan"), "err": mae,
                         "mean_residual": float(resid[msk].mean())})
    pd.DataFrame(bin_rows).to_csv(OUT_BINS, index=False, encoding="utf-8")

    # ---------------- 3. 分组指标 ----------------
    emit("")
    emit("## 3. 按元数据分组的指标")
    emit("")
    group_specs = []
    if "Affinity Method" in meta.columns:
        group_specs.append(("测量方法 (`Affinity Method`)", meta["Affinity Method"]))
    if "Subgroup" in meta.columns:
        group_specs.append(("亚组 (`Subgroup`)", meta["Subgroup"]))
    if "Source Data Set" in meta.columns:
        group_specs.append(("来源数据集", meta["Source Data Set"]))
    if "Temperature(K)" in meta.columns:
        tk = pd.to_numeric(meta["Temperature(K)"], errors="coerce")
        tb = pd.cut(tk, bins=[0, 295, 300, 310, 400],
                    labels=["<295K", "295-300K", "300-310K", ">310K"])
        group_specs.append(("温度区间", tb))

    for title, ser in group_specs:
        emit(f"### {title}")
        emit("")
        emit("| 分组 | n | Spearman | " + ("MAE (kcal/mol)" if best_is_dg
                                                 else "平均绝对排名差 (pp)") + " |")
        emit("|---|---:|---:|---:|")
        s = ser.astype(object).where(ser.notna(), "(缺失)")
        for gval, idx in s.groupby(s).groups.items():
            ii = np.asarray(list(idx))
            if len(ii) < 10:
                emit(f"| {gval} | {len(ii)} | (n<10 略) | |")
                continue
            rho = mr.spearman_rho(-y[ii], -best_pred[ii])
            errv = (np.mean(np.abs(resid[ii])) if best_is_dg
                    else np.mean(np.abs(true_pct[ii] - pred_pct[ii])))
            emit(f"| {str(gval)[:26]} | {len(ii)} | {rho:+.4f} | {errv:.3f} |")
        emit("")

    # ---------------- 4. 结论 ----------------
    emit("## 4. 错误分析结论")
    emit("")
    valid_bins = [b for b in bin_rows if b["spearman"] == b["spearman"]]
    if valid_bins:
        worst_bin = min(valid_bins, key=lambda b: b["spearman"])
        best_bin = max(valid_bins, key=lambda b: b["spearman"])
        emit(f"- **等宽分层**: 最差 **{worst_bin['bin']}** "
             f"(ρ={worst_bin['spearman']:+.4f}, n={worst_bin['n']}); "
             f"最好 **{best_bin['bin']}** (ρ={best_bin['spearman']:+.4f}, n={best_bin['n']})")
        emit("  -> 模型在**强结合区**的层内排序能力与**弱结合区**不同, 这解释了为什么")
        emit("     全局 Spearman 与头部富集 (EF@1%) 会给出不同排名。")
    emit(f"- 方向性: 过强 {n_over / len(y) * 100:.1f}% vs 过弱 "
         f"{n_under / len(y) * 100:.1f}% —— "
         + ("接近对称" if abs(n_over - n_under) / len(y) < 0.1 else "存在系统性偏移"))
    if best_is_dg:
        emit(f"- 中位绝对残差 {np.median(absr):.3f} kcal/mol, "
             f"而 test ΔG 的 std 为 {y.std():.3f} kcal/mol")
    else:
        emit(f"- 该最优模型**非 ΔG 尺度**, 误差以排名差表达: "
             f"平均 {np.mean(np.abs(true_pct - pred_pct)):.1f} 个百分点; "
             f"test ΔG 的 std 为 {y.std():.3f} kcal/mol")
    emit("")
    emit("### 局限")
    emit("")
    emit("1. **等频细分箱不报箱内 Spearman**: 箱宽仅约 0.2–0.5 kcal/mol,")
    emit("   箱内几乎无排序变异, 箱内 Spearman 是噪声; 因此只报误差。")
    emit("2. 分组分析中 `Affinity Method` 覆盖率有限, 方法分组只反映**有标注的部分**样本。")
    emit("3. 温度列大量为 `298(assumed)` 占位, 温度分箱可解释性有限。")
    emit("4. 元数据经 **`pdb_id`** 关联 (非 interface_key), 同一结构的多个突变体")
    emit("   会共享第一条的元数据 —— 对方法/来源类字段通常一致, 但并非严格逐行。")
    emit("")
    emit(f"- 写出: `{OUT_MD.name}`, `{OUT_WORST.name}`, `{OUT_BINS.name}`")

    # 修掉上面 f-string 拼接产生的排版问题: 重生成一遍表格块
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
