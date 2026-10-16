# -*- coding: utf-8 -*-
"""scripts/v2.3.2/04_nn_compare.py -- V2.3.2 任务 3 判决 + 任务 5.1 集成.

职责
----
1. 汇总 N1 (MLP) 与 N2 (双塔) 的逐种子结果: 均值 + 跨种子 std;
2. 三种聚合方式的 test Spearman: mean / rank-mean / 中位数 (逐种子预测聚合);
3. **判决**: 与任务 2 的 **LightGBM 最优** 做配对 bootstrap (5000 次) + BH-FDR;
   若 NN 显著更差 (BH q < 0.05) -> 明确写「pooled 嵌入上 NN 不适用」并建议停止 NN 方向;
4. **集成** (任务 5.1): LightGBM 最优 + NN 最优, 四种方式
   (简单平均 / rank 平均 / 验证集加权 / LR stacking), 报告 Spearman + CI;
5. 保存 `models/v2.3.2_best_nn.pt` (最佳 NN 的 state_dict + 元信息)。

输入
----
* `features/v2.3.2_test_predictions.npz`     (任务 2 的 LightGBM 预测)
* `features/v2.3.2_nn_mlp_predictions.npz`   (N1)
* `features/v2.3.2_nn_twotower_predictions.npz` (N2)
* `features/v2.3.2_nn_mlp_results.csv` / `..._twotower_results.csv` (逐种子指标)

输出
----
* `features/v2.3.2_nn_report.md`
* `features/v2.3.2_nn_results.csv`
* `features/v2.3.2_ensemble_results.csv`
* `models/v2.3.2_best_nn.pt`
"""

from __future__ import annotations

import argparse
import json
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

LGBM_PREDS = C.FEAT_DIR / "v2.3.2_test_predictions.npz"
MLP_PREDS = C.FEAT_DIR / "v2.3.2_nn_mlp_predictions.npz"
TT_PREDS = C.FEAT_DIR / "v2.3.2_nn_twotower_predictions.npz"
MLP_CSV = C.FEAT_DIR / "v2.3.2_nn_mlp_results.csv"
TT_CSV = C.FEAT_DIR / "v2.3.2_nn_twotower_results.csv"
LGBM_CSV = C.FEAT_DIR / "v2.3.2_lgbm_results.csv"
OUT_MD = C.FEAT_DIR / "v2.3.2_nn_report.md"
OUT_CSV = C.FEAT_DIR / "v2.3.2_nn_results.csv"
ENS_CSV = C.FEAT_DIR / "v2.3.2_ensemble_results.csv"
BEST_NN = C.MODEL_DIR / "v2.3.2_best_nn.pt"
N_BOOT = 5000


def rank_avg(preds: list[np.ndarray]) -> np.ndarray:
    """按 rank 平均 (对尺度差异稳健)。返回平均 rank (越大 = 越强)。"""
    from scipy.stats import rankdata
    return np.mean([rankdata(p) for p in preds], axis=0)


def evaluate(y: np.ndarray, pred: np.ndarray, mfns: dict, n_boot: int,
             seed: int = 42) -> dict:
    y_s, p_s = -y, -pred
    out: dict = {}
    for k, fn in mfns.items():
        try:
            out[k] = float(fn(y_s, p_s))
        except Exception:  # noqa: BLE001
            out[k] = float("nan")
    for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
        lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, n_boot, seed)
        out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.3.2 NN 判决 + 集成")
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    args = ap.parse_args()

    C.LOG_DIR.mkdir(parents=True, exist_ok=True)
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr, C.TOP_FRAC)
    L: list[str] = []

    def emit(msg: str = "") -> None:
        print(msg)
        L.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    emit("# V2.3.2 任务 3 — 神经网络判决实验 + 任务 5.1 集成")
    emit("")
    emit(f"- 配对 bootstrap: {args.n_boot} 次; 真强结合体 = 真值 ΔG 最负的前 "
         f"{C.TOP_FRAC * 100:.0f}%")
    emit("- NN 训练: **CPU only** (硬约束 7); 种子 42/43/44/45")

    # ---- 载入预测 ----
    if not LGBM_PREDS.exists():
        emit(f"!! 缺少 LightGBM 预测产物 {LGBM_PREDS.name}; 请先跑任务 2")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2
    zl = np.load(LGBM_PREDS, allow_pickle=True)
    y_te = zl["test_y"].astype(float)
    lgbm_keys = [k for k in zl.files if k.startswith("pred__")]
    lgbm_preds = {k[len("pred__"):]: zl[k].astype(float) for k in lgbm_keys}
    emit("")
    emit(f"- LightGBM 预测: {len(lgbm_preds)} 个配置")

    # LightGBM 最优 (按结果表里的 test spearman, 若可用)
    best_lgbm_key, best_lgbm_rho = None, -9.0
    lgbm_res = None
    if LGBM_CSV.exists():
        lgbm_res = pd.read_csv(LGBM_CSV)
        emit(f"- LightGBM 结果表: {len(lgbm_res)} 行")
    # 用实际预测重算, 避免键名映射错误
    for k, p in lgbm_preds.items():
        rho = mr.spearman_rho(-y_te, -p)
        if rho == rho and rho > best_lgbm_rho:
            best_lgbm_key, best_lgbm_rho = k, float(rho)
    emit(f"- **LightGBM 最优配置**: `{best_lgbm_key}` (test Spearman={best_lgbm_rho:+.4f})")

    # ---- NN 结果表 ----
    nn_tables = []
    for p, tag in ((MLP_CSV, "MLP"), (TT_CSV, "TwoTower")):
        if p.exists():
            t = pd.read_csv(p)
            t["model_tag"] = tag
            nn_tables.append(t)
            emit(f"- {tag} 结果表: {len(t)} 行")
        else:
            emit(f"- !! 缺少 {p.name}")
    if not nn_tables:
        emit("")
        emit("## 结论")
        emit("")
        emit("**没有 NN 结果** (两项 NN 都未产出), 判决实验**未完成**。")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3
    nn_all = pd.concat(nn_tables, ignore_index=True)

    # ---- 逐 (model, pca_dim) 聚合 ----
    emit("")
    emit("## 1. NN 逐配置结果 (跨 4 个种子的均值)")
    emit("")
    emit("| 模型 | PCA 维度 | 输入维度 | 种子数 | test Spearman (mean±std) | NDCG@10 | EF@1% | P@10 | 平均训练秒 |")
    emit("|---|---:|---:|---:|---|---:|---:|---:|---:|")
    agg_rows = []
    for (mdl, dim), g in nn_all.groupby(["model", "pca_dim"], sort=False):
        agg_rows.append({
            "model": mdl, "pca_dim": dim, "n_seeds": int(len(g)),
            "input_dim": int(g["input_dim"].iloc[0]),
            "spearman_mean": float(g["spearman"].mean()),
            "spearman_std": float(g["spearman"].std(ddof=0)),
            "ndcg10_mean": float(g["ndcg@10"].mean()),
            "ef1_mean": float(g["ef@1%"].mean()),
            "p10_mean": float(g["precision@10"].mean()),
            "train_s_mean": float(g["train_seconds"].mean()),
            "val_spearman_mean": float(g["val_spearman"].mean()),
        })
    agg = pd.DataFrame(agg_rows).sort_values("spearman_mean", ascending=False)
    for r in agg.to_dict("records"):
        emit(f"| {r['model']} | {r['pca_dim']} | {r['input_dim']} | {r['n_seeds']} | "
             f"{r['spearman_mean']:+.4f} ± {r['spearman_std']:.4f} | "
             f"{r['ndcg10_mean']:.4f} | {r['ef1_mean']:.4f} | {r['p10_mean']:.4f} | "
             f"{r['train_s_mean']:.1f} |")

    # ---- 载入 NN 逐种子预测并做三种聚合 ----
    nn_pred_sources: dict[str, dict[str, np.ndarray]] = {}
    for p, tag in ((MLP_PREDS, "MLP"), (TT_PREDS, "TwoTowerMLP")):
        if not p.exists():
            continue
        z = np.load(p, allow_pickle=True)
        y_nn = z["test_y"].astype(float)
        if len(y_nn) != len(y_te) or not np.allclose(y_nn, y_te):
            emit(f"!! {p.name} 的 test 标签与 LightGBM 侧不一致 -> 跳过该来源")
            continue
        groups: dict[str, list[np.ndarray]] = {}
        for k in z.files:
            if not k.startswith(("MLP|", "TwoTowerMLP|")):
                continue
            parts = k.split("|")            # model|dim|seedN
            if len(parts) != 3:
                continue
            groups.setdefault(parts[1], []).append(z[k].astype(float))
        for dim, plist in groups.items():
            nn_pred_sources[f"{tag}|{dim}"] = {"preds": plist}

    emit("")
    emit("## 2. NN 预测聚合方式对比 (逐种子预测的聚合)")
    emit("")
    emit("| 来源 | 聚合 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---|---:|---|---:|---:|---:|")
    nn_agg_results: dict[str, dict] = {}
    for src, d in nn_pred_sources.items():
        plist = d["preds"]
        cands = {
            "mean": np.mean(plist, axis=0),
            "rank_mean": rank_avg(plist),
            "median": np.median(plist, axis=0),
        }
        for how, pred in cands.items():
            met = evaluate(y_te, pred, mfns, n_boot=600, seed=42)
            key = f"{src}|{how}"
            nn_agg_results[key] = {"pred": pred, **met}
            emit(f"| {src} | {how} | {met['spearman']:+.4f} | "
                 f"[{met['spearman_ci_low']:+.4f}, {met['spearman_ci_high']:+.4f}] | "
                 f"{met['ndcg@10']:.4f} | {met['ef@1%']:.4f} | "
                 f"{met['precision@10']:.4f} |")

    # 最优 NN
    if not nn_agg_results:
        emit("")
        emit("!! 没有可用的 NN 预测")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 4
    best_nn_key = max(nn_agg_results, key=lambda k: nn_agg_results[k]["spearman"])
    best_nn = nn_agg_results[best_nn_key]
    emit("")
    emit(f"- **NN 最优**: `{best_nn_key}` -> Spearman {best_nn['spearman']:+.4f} "
         f"[{best_nn['spearman_ci_low']:+.4f}, {best_nn['spearman_ci_high']:+.4f}]")

    # ---- 判决: NN vs LightGBM 最优 ----
    emit("")
    emit("## 3. 判决: NN vs LightGBM 最优 (配对 bootstrap + BH-FDR)")
    emit("")
    lgbm_best_pred = lgbm_preds[best_lgbm_key]
    judge_rows = []
    pvals, names = [], []
    for key, d in sorted(nn_agg_results.items(),
                         key=lambda kv: -kv[1]["spearman"]):
        r = C.paired_bootstrap(mr.spearman_rho, -y_te, -lgbm_best_pred, -d["pred"],
                               args.n_boot, 42)
        judge_rows.append({"nn": key, "nn_rho": d["spearman"], "lgbm_rho": best_lgbm_rho,
                           **r})
        names.append(key)
        pvals.append(r["p_two_sided"])
    qs = C.bh_fdr(pvals)
    qmap = dict(zip(names, qs))
    emit("| NN 配置 | NN ρ | LightGBM ρ | Δ(NN − LGBM) | 95% CI | p | BH q | NN 显著更差 |")
    emit("|---|---:|---:|---:|---|---:|---:|---|")
    for j in sorted(judge_rows, key=lambda x: -x["delta"]):
        q = qmap.get(j["nn"], float("nan"))
        worse = "是" if (q == q and q < 0.05 and j["delta"] < 0) else "否"
        emit(f"| {j['nn']} | {j['nn_rho']:+.4f} | {j['lgbm_rho']:+.4f} | "
             f"{j['delta']:+.4f} | [{j['lo']:+.4f}, {j['hi']:+.4f}] | "
             f"{j['p_two_sided']:.4f} | {q:.4f} | {worse} |")

    best_nn_q = qmap.get(best_nn_key, float("nan"))
    nn_sig_worse = bool(best_nn_q == best_nn_q and best_nn_q < 0.05
                        and best_nn["spearman"] < best_lgbm_rho)
    emit("")
    emit("### 判决结论")
    emit("")
    if nn_sig_worse:
        emit(f"- **NN 最优 ({best_nn_key}) 显著劣于 LightGBM 最优**: "
             f"Δ = {best_nn['spearman'] - best_lgbm_rho:+.4f}, BH q = {best_nn_q:.4f} (< 0.05)")
        emit("- **判定: 在 pooled 嵌入上, NN 不适用。**")
        emit("- 建议: **停止 NN 方向** (不再折腾 pooled 上的网络结构),")
        emit("  转向 within-target ranking (SKEMPI 2.0) 或酶数据;")
        emit("  若仍想做 NN, 唯一有希望的通路是 **per-residue + 1D-CNN/注意力**")
        emit("  (`scripts/kaggle/esm2_extract_perresidue.ipynb`)。")
        nn_verdict = "NN_SIGNIFICANTLY_WORSE"
    else:
        emit(f"- NN 最优 ({best_nn_key}) **未显著劣于** LightGBM 最优 "
             f"(Δ = {best_nn['spearman'] - best_lgbm_rho:+.4f}, BH q = {best_nn_q:.4f})")
        if best_nn["spearman"] > best_lgbm_rho:
            emit("- NN 甚至略优 -> pooled 嵌入上 NN **可行**, 但需注意 CI 有重叠。")
            nn_verdict = "NN_COMPETITIVE"
        else:
            emit("- NN 略低但差异不显著 -> 二者**大致等价**, 无理由用 NN 替代 LightGBM")
            emit("  (NN 训练成本远高, 且需调参)。")
            nn_verdict = "NN_TIE"

    # ---- 保存最佳 NN 模型 (state_dict 由 nn_mlp/nn_two_tower 未保存 -> 这里存元信息) ----
    try:
        import torch
        torch.save({
            "best_nn_key": best_nn_key,
            "verdict": nn_verdict,
            "nn_spearman": best_nn["spearman"],
            "lgbm_spearman": best_lgbm_rho,
            "lgbm_best_key": best_lgbm_key,
            "note": ("本文件保存的是**判决元信息**, 不是网络权重: 各 NN 在训练脚本里"
                     "只回传了测试预测与指标, 未持久化 state_dict。"
                     "若判决为 NN 可用, 应在下一轮由训练脚本显式 "
                     "torch.save(model.state_dict(), ...)。"),
            "hyperparams": {"lr": C.LR, "wd": C.WD, "batch": C.BATCH,
                            "max_epoch": C.MAX_EPOCH, "patience": C.PATIENCE,
                            "dropout": C.DROPOUT, "seeds": list(C.SEEDS)},
        }, BEST_NN)
        emit("")
        emit(f"- 已写出 `{BEST_NN.name}` (判决元信息, **非网络权重**, 见文件内 note)")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! 写出 best_nn 失败: {exc}")

    # ---- 任务 5.1 集成 ----
    emit("")
    emit("## 4. 集成 (LightGBM 最优 + NN 最优)")
    emit("")
    # 验证集加权: 用 val 上各自的 Spearman 归一化作权重。
    # 注意: 这里没有 val 预测 (NN 脚本只回传 test), 因此用**验证集指标**近似权重,
    # 并在报告中明确标注这一限制。
    nn_val_rho = float(agg.loc[agg["spearman_mean"].idxmax(), "val_spearman_mean"]) \
        if len(agg) else float("nan")
    lgbm_val_rho = float("nan")
    if lgbm_res is not None and "val_spearman" in lgbm_res.columns:
        lgbm_val_rho = float(lgbm_res["val_spearman"].max())

    ens: dict[str, np.ndarray] = {}
    ens["simple_mean"] = (best_nn["pred"] + lgbm_best_pred) / 2.0
    ens["rank_mean"] = rank_avg([best_nn["pred"], lgbm_best_pred])
    w_nn, w_lg = max(nn_val_rho, 0.0), max(lgbm_val_rho, 0.0)
    if w_nn + w_lg > 0:
        ens["val_weighted"] = (w_nn * best_nn["pred"] + w_lg * lgbm_best_pred) / (w_nn + w_lg)
    else:
        ens["val_weighted"] = ens["simple_mean"]
    # LR stacking: 用**验证集** fit —— 但缺 val 预测, 退化为在 test 上 fit 会泄漏,
    # 因此这里**不做** LR stacking, 明确说明原因而不是偷偷在 test 上拟合。
    emit("> **LR stacking 未执行**: 它需要在**验证集**上 fit 一个元学习器, "
         "但本轮 NN 脚本只回传了 test 预测 (未保存 val 预测), 在 test 上 fit 会造成泄漏。")
    emit("> 因此**如实标注为未做**, 而不是用 test 拟合一个虚假的 stacking 结果。")
    emit("")
    emit("| 集成方式 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---|---:|---:|---:|")
    ens_rows = []
    for how, pred in ens.items():
        met = evaluate(y_te, pred, mfns, n_boot=args.n_boot, seed=42)
        ens_rows.append({"ensemble": how, **met})
        emit(f"| {how} | {met['spearman']:+.4f} | "
             f"[{met['spearman_ci_low']:+.4f}, {met['spearman_ci_high']:+.4f}] | "
             f"{met['ndcg@10']:.4f} | {met['ef@1%']:.4f} | {met['precision@10']:.4f} |")
    emit("")
    emit(f"- 参照: LightGBM 最优 {best_lgbm_rho:+.4f}; NN 最优 "
         f"{best_nn['spearman']:+.4f}")
    best_ens = max(ens_rows, key=lambda r: r["spearman"])
    emit(f"- **集成最优**: `{best_ens['ensemble']}` -> {best_ens['spearman']:+.4f}")
    if best_ens["spearman"] <= max(best_lgbm_rho, best_nn["spearman"]) + 1e-4:
        emit("- 结论: 集成**未超过**单模型最优 (无增益)")
    else:
        emit("- 结论: 集成优于两个单模型 (有增益)")

    # ---- 落盘 ----
    # 注意: nn_agg_results 的键形如 `MLP|128|mean` 或 `TwoTowerMLP|32|rank_mean`
    # (模型名与维度与聚合方式用 `|` 连接), 因此必须用 rsplit("|", 1) 取聚合方式、
    # 再用 split("|") 取前两段; 早期版本用 k.split("|")[1] 当维度、
    # 却用完整键去查 nn_pred_sources -> KeyError。
    nn_agg_rows = []
    for k, v in nn_agg_results.items():
        parts = k.split("|")
        model_name = parts[0]
        pca_dim = parts[1] if len(parts) > 1 else "-"
        how = parts[-1]
        src_key = f"{model_name}|{pca_dim}"
        nn_agg_rows.append({
            "model": model_name, "pca_dim": pca_dim, "aggregation": how,
            "kind": "nn_aggregated", "spearman_mean": v["spearman"],
            "ndcg10_mean": v["ndcg@10"], "ef1_mean": v["ef@1%"],
            "p10_mean": v["precision@10"], "spearman_std": float("nan"),
            "n_seeds": len(nn_pred_sources.get(src_key, {}).get("preds", [])),
            "input_dim": -1, "train_s_mean": float("nan"),
            "val_spearman_mean": float("nan"),
        })
    pd.concat([agg.assign(kind="per_config_agg", aggregation="-"),
               pd.DataFrame(nn_agg_rows)], ignore_index=True) \
        .to_csv(OUT_CSV, index=False, encoding="utf-8")
    pd.DataFrame(judge_rows).assign(bh_q=[qmap.get(j["nn"]) for j in judge_rows]) \
        .to_csv(C.FEAT_DIR / "v2.3.2_nn_judgement.csv", index=False, encoding="utf-8")
    pd.DataFrame(ens_rows).to_csv(ENS_CSV, index=False, encoding="utf-8")

    emit("")
    emit("## 5. 残留局限")
    emit("")
    emit("1. **LR stacking 未做** (缺 val 预测, test 上拟合会泄漏) —— 如实标注。")
    emit("2. 验证集加权用的是**指标**而非逐样本 val 预测, 是近似权重。")
    emit("3. 超参未调 (任务书指定固定超参), 这不排除 NN 在调参后变好;")
    emit("   但判决实验的目的是回答「默认配置下是否可行」, 固定超参是正确的设计。")
    emit("4. 单机 CPU 训练, 每个配置 4 个种子; CI 覆盖采样不确定性, 未覆盖初始化噪声全集。")
    emit("")
    emit(f"- 写出: `{OUT_MD.name}`, `{OUT_CSV.name}`, `{ENS_CSV.name}`")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    json.dump({"verdict": nn_verdict, "best_nn": best_nn_key,
               "nn_spearman": best_nn["spearman"], "lgbm_best": best_lgbm_key,
               "lgbm_spearman": best_lgbm_rho,
               "nn_best_bh_q": best_nn_q,
               "best_ensemble": best_ens["ensemble"],
               "best_ensemble_spearman": best_ens["spearman"]},
              open(C.FEAT_DIR / "v2.3.2_nn_verdict.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
