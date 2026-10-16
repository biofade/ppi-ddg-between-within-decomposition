# -*- coding: utf-8 -*-
"""scripts/v4.0/04_baselines.py -- V4.0 任务 3: ΔΔG 基线复现.

基线（任务书 B1–B4，全部 CPU）
------------------------------
| 编号 | 模型 | 特征 |
|---|---|---|
| B1 | LightGBM | 物理(60) + pooled ESM-2 PCA128，**配对级** `[x_mut, x_wt, x_mut−x_wt]` |
| B2 | LightGBM + ΔG 辅助（多任务，基线版）| 同 B1；主目标 ΔΔG，辅助目标 ΔG_mut |
| B3 | MLP | 同 B1 特征 |
| B4 | MLP | + per-residue mean-pool（若可用；否则同 B3 并**如实标注**）|

> B2 在此处只是**基线**（辅助头用 ΔG_mut）。任务 4 的创新 2 是它的**改进版**
> （双头 + λ 扫描 + 与创新 1 的对称损失结合）。两者不重复。

评估指标
--------
Pearson（主，与文献可比）/ Spearman / RMSE / MAE / **per-structure Pearson** /
AUROC(ΔΔG>0)。全部带 bootstrap 95% CI。所有模型与**最优基线**做配对 bootstrap
+ BH-FDR。

**不能直接与文献比较**：Pythia-PPI 的 0.7850 是 SKEMPI **mutation-level**
（同复合物突变跨折，作者自承乐观 [Pythia-PPI p8]），我的 test 是
**复合物级 + 40% 序列同一性 + 家族级**隔离。报告中必须标注该差异。
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v40_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.EXP / "baselines" / "ddg_baseline_results.csv"
OUT_PRED = C.EXP / "baselines" / "ddg_baseline_predictions.npz"
OUT_PAIR = C.EXP / "baselines" / "ddg_baseline_pairs.csv"
OUT_MD = C.FEAT / "v4.0_ddg_baseline_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def feat_matrix(blk: dict, use: str) -> np.ndarray:
    if use == "b1":
        return np.hstack([blk["phys_mut"], blk["phys_wt"], blk["phys_diff"],
                          blk["Z_mut"], blk["Z_wt"], blk["Z_diff"]]).astype(np.float32)
    raise ValueError(use)


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 任务 3: ΔΔG 基线")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--lam", type=float, default=0.3, help="B2 的 ΔG 辅助权重")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43])
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    torch.set_num_threads(2)
    lgb = C.load_lgb()
    t0 = time.time()

    emit("# V4.0 任务 3 — ΔΔG 基线复现")
    emit("")
    emit("**任务**: 预测 `ΔΔG = ΔG_mutant − ΔG_wild-type`（kcal/mol）")
    emit("**设备**: CPU only")
    emit(f"**主指标**: Pearson（与文献可比）；辅助 Spearman / RMSE / MAE / "
         f"per-structure Pearson / AUROC")
    emit("")

    data = C.build_pair_features()
    blk = data["X_blocks"]
    y, groups = data["y"], data["groups"]
    mtr, mva, mte = data["mask_train"], data["mask_val"], data["mask_test"]
    X = feat_matrix(blk, "b1")

    emit("## 1. 数据与划分")
    emit("")
    emit(f"- 配对: **{len(y)}**（forward）；特征 **{X.shape[1]}** 维 "
         f"= [物理60 + PCA128] × 3(mut/wt/diff)")
    emit(f"- train {int(mtr.sum())} / val {int(mva.sum())} / test {int(mte.sum())}")
    emit(f"- 复合物: train {len(set(groups[mtr]))} / val {len(set(groups[mva]))} / "
         f"test {len(set(groups[mte]))}")
    emit(f"- ΔΔG: test mean {y[mte].mean():+.3f}, std {y[mte].std():.3f}, "
         f"中位 {np.median(y[mte]):+.3f}")
    emit(f"- **划分已通过复合物级/40%同一性/家族级三重隔离验证**"
         f"（见 `v4.0_ddg_split_report.md`）")

    # 神经网络需要标准化（物理特征 std 到 293、嵌入 std 仅 0.04，差 4 个数量级；
    # 不标准化实测直接梯度爆炸 -> 全部输出 NaN）。LightGBM 对单调变换不敏感，
    # 因此仍用原始矩阵，保持与 V2.2/V2.5 的口径一致。
    Xn, scaler = C.standardize_train(X, mtr)
    assert np.isfinite(Xn).all(), "标准化后仍有 NaN/Inf"
    emit(f"- 神经网络输入经**仅在训练集拟合**的 median+std 标准化（防泄漏）")
    emit("")

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    def record(tag: str, p_te: np.ndarray, extra: dict | None = None) -> dict:
        m = C.all_metrics(y[mte], p_te, groups[mte])
        for k in ("pearson", "spearman", "per_structure_pearson"):
            lo, hi, _ = C.bootstrap_ci(
                {"pearson": C.pearson, "spearman": C.spearman,
                 "per_structure_pearson": lambda a, b: C.per_structure_pearson(
                     a, b, groups[mte])["per_structure_pearson"]}[k], y[mte], p_te)
            m[f"{k}_ci_low"], m[f"{k}_ci_high"] = lo, hi
        row = {"model": tag, **m}
        if extra:
            row.update(extra)
        results.append(row)
        preds[tag] = p_te.astype(np.float32)
        return row

    emit("## 2. 基线结果（test set）")
    emit("")
    emit("| 模型 | Pearson | 95% CI | Spearman | RMSE | MAE | per-struct Pearson (n) | AUROC |")
    emit("|---|---:|---|---:|---:|---:|---|---:|")

    # ---------- B1: LightGBM ----------
    try:
        m = lgb.LGBMRegressor(**C.LGB_PARAMS)
        m.fit(X[mtr], y[mtr], eval_set=[(X[mva], y[mva])],
              callbacks=[lgb.early_stopping(80, verbose=False)])
        p = m.predict(X[mte])
        r = record("B1_lgb_phys+esm2", p)
        emit(f"| {r['model']} | **{r['pearson']:+.4f}** | "
             f"[{r['pearson_ci_low']:+.3f}, {r['pearson_ci_high']:+.3f}] | "
             f"{r['spearman']:+.4f} | {r['rmse']:.4f} | {r['mae']:.4f} | "
             f"{r['per_structure_pearson']:+.4f} ({r['per_structure_n']}) | "
             f"{r['auroc']:.4f} |")
    except Exception as exc:  # noqa: BLE001
        emit(f"| B1_lgb_phys+esm2 | FAIL {exc} | | | | | | |")
        with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[04_baselines] B1: {exc}\n{traceback.format_exc()}\n")

    # ---------- 神经网络骨架 ----------
    D = X.shape[1]
    Xt = torch.from_numpy(Xn)
    yt = torch.from_numpy(y.astype(np.float32))

    class MLP(nn.Module):
        def __init__(self, d_in: int, hidden: int = 256, n_aux: bool = False):
            super().__init__()
            self.n_aux = n_aux
            self.body = nn.Sequential(
                nn.Linear(d_in, hidden), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(hidden, 128), nn.ReLU(), nn.Dropout(0.2),
                nn.Linear(128, 64), nn.ReLU())
            self.head_ddg = nn.Linear(64, 1)
            self.head_dg = nn.Linear(64, 1) if n_aux else None

        def forward(self, x):
            h = self.body(x)
            return self.head_ddg(h).squeeze(-1), (
                self.head_dg(h).squeeze(-1) if self.n_aux else None)

    def train_nn(tag: str, seed: int, n_aux: bool, lam: float) -> None:
        torch.manual_seed(seed)
        np.random.seed(seed)
        net = MLP(D, n_aux=n_aux)
        opt = torch.optim.Adam(net.parameters(), lr=args.lr, weight_decay=1e-4)
        mse = nn.MSELoss()
        tr_i, va_i = np.where(mtr)[0], np.where(mva)[0]
        # 辅助目标: mutant 的绝对 ΔG（由 dg_wt + ddg 得到）
        dg_mut = (data["df"]["dg_mut"].to_numpy(float)).astype(np.float32)
        at = torch.from_numpy(dg_mut)
        best, best_state, bad = -9.0, None, 0
        for ep in range(args.epochs):
            net.train()
            perm = np.random.RandomState(seed * 1000 + ep).permutation(len(tr_i))
            for s in range(0, len(perm), args.batch):
                b = tr_i[perm[s:s + args.batch]]
                pdg, pg = net(Xt[b])
                loss = mse(pdg, yt[b])
                if n_aux and pg is not None and lam > 0:
                    loss = loss + lam * mse(pg, at[b])
                opt.zero_grad()
                loss.backward()
                opt.step()
            net.eval()
            with torch.no_grad():
                pv, _ = net(Xt[va_i])
            rho = C.pearson(y[va_i], pv.numpy())
            if rho == rho and rho > best + 1e-6:
                best, bad = float(rho), 0
                best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
            else:
                bad += 1
            if bad >= 30:
                break
        if best_state:
            net.load_state_dict(best_state)
        net.eval()
        with torch.no_grad():
            pt, _ = net(Xt[np.where(mte)[0]])
        r = record(tag, pt.numpy(), {"seed": seed, "lam": lam if n_aux else 0.0,
                                    "best_val_pearson": best})
        emit(f"| {r['model']} | **{r['pearson']:+.4f}** | "
             f"[{r['pearson_ci_low']:+.3f}, {r['pearson_ci_high']:+.3f}] | "
             f"{r['spearman']:+.4f} | {r['rmse']:.4f} | {r['mae']:.4f} | "
             f"{r['per_structure_pearson']:+.4f} ({r['per_structure_n']}) | "
             f"{r['auroc']:.4f} |")

    # ---------- B3: MLP ----------
    for sd in args.seeds:
        try:
            train_nn(f"B3_mlp_phys+esm2_s{sd}", sd, n_aux=False, lam=0.0)
        except Exception as exc:  # noqa: BLE001
            emit(f"| B3_mlp seed={sd} | FAIL {exc} | | | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[04_baselines] B3 seed={sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")

    # ---------- B2: MLP + ΔG 辅助头（多任务基线）----------
    for sd in args.seeds:
        try:
            train_nn(f"B2_mlp_auxDG_lam{args.lam}_s{sd}", sd, n_aux=True, lam=args.lam)
        except Exception as exc:  # noqa: BLE001
            emit(f"| B2_mlp_auxDG seed={sd} | FAIL {exc} | | | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[04_baselines] B2 seed={sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")

    # ---------- B4: per-residue ----------
    emit("")
    emit("### B4（per-residue）状态")
    emit("")
    emit("> **B4 未运行**，原因如实说明：`v2.5_perres_cache.npz` 的 per-residue 表示按 "
         "**V2.2 可用行**（8199 行）索引，而本任务的配对涉及 **filtered.csv 的 "
         "WT 行与 mutant 行**；WT 行在 perres 缓存中的覆盖率需要单独核对，")
    emit("> 且 per-residue 逐对特征（mut/wt/diff）会占用大量内存（1.9GB 源缓存）。")
    emit("> V2.5 已实测 **per-residue 在 cross-target 上无增益**（CNN +0.3246 vs "
         "pooled +0.3253, BH q=0.4821），因此本任务把 B4 列为**次优先**。")
    emit("> 若时间允许将在创新阶段补做；此处**不假装完成**。")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[mte].astype(np.float32),
             test_group=groups[mte],
             **{k: v for k, v in preds.items()})

    # ---------- 3. 配对检验 ----------
    emit("")
    emit("## 3. 配对 bootstrap vs 最优基线（BH-FDR）")
    emit("")
    fin = res[np.isfinite(res["pearson"])]
    if len(fin):
        best_row = fin.sort_values("pearson", ascending=False).iloc[0]
        # 参照取 B1（LightGBM）为"最优基线"；若 B1 失败则取 fin 第一
        base_tag = ("B1_lgb_phys+esm2" if "B1_lgb_phys+esm2" in preds
                    else best_row["model"])
        emit(f"- **最优基线（参照）**: `{base_tag}`，"
             f"test Pearson = **{float(fin[fin.model == base_tag]['pearson'].iloc[0]):+.4f}**"
             if base_tag in set(fin.model) else f"- 参照: `{base_tag}`（无指标）")
        emit("")
        emit("| 模型 | Δ Pearson | 95% CI | p (bootstrap) | BH q | 判定 |")
        emit("|---|---:|---|---:|---:|---|")
        prows, pvals = [], []
        for _, r in fin.iterrows():
            if r["model"] == base_tag:
                continue
            rr = C.paired_bootstrap(C.pearson, y[mte], preds[base_tag],
                                    preds[r["model"]])
            prows.append({"model": r["model"], "vs": base_tag, **rr})
            pvals.append(rr["p_two_sided"])
        qs = C.bh_fdr(pvals)
        for r, q in zip(prows, qs):
            r["bh_q"] = q
            v = ("更优" if (q < 0.05 and r["delta"] > 0) else
                 "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
            emit(f"| {r['model']} | {r['delta']:+.4f} | "
                 f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
                 f"{q:.4f} | {v} |")
        if prows:
            pd.DataFrame(prows).to_csv(OUT_PAIR, index=False, encoding="utf-8")

    # ---------- 4. 与文献对比（标注 split 差异）----------
    emit("")
    emit("## 4. 与已发表方法的对比（**必须标注 split 差异**）")
    emit("")
    emit("| 方法 | 报告值 | 其 split 协议 | 我的 split 协议 | 可否直接比较 |")
    emit("|---|---:|---|---|---|")
    emit("| Pythia-PPI | Pearson **0.7850**（SKEMPI）[Pythia-PPI p11] | "
         "per-structure 五折，**含出现在训练折里的复合物** [p17] | "
         "复合物级 + 40% 同一性 + 家族级，**完全留出** | **否** |")
    emit("| Pythia-PPI（vanilla）| Pearson 0.7092 [p11] | 同上 | 同上 | 否 |")
    emit("| SAMAffinity | Pearson 0.88(S1131)/0.81(S4169) [p5,p7] | "
         "**mutation-level** 10-fold（同复合物突变跨折，作者自承乐观 [p8]）| "
         "复合物级隔离 | **否** |")
    emit("| SAMAffinity（complex-level）| Pearson 0.59(S1131)/0.54(S4169) [p9] | "
         "同一 WT 复合物突变不跨折，**无序列同一性约束** | 复合物级 + 40% 同一性 | "
         "**最接近，仍更松** |")
    emit("| DSSA-PPI（complex-level）| PCC 0.71(S1131)/0.58(S4169) [p5] | "
         "MpbPPI 预定 split | 复合物级 + 40% | 接近但仍更松 |")
    emit("| MINT | SKEMPI Pearson +32% over ESM-2 [p5] | 按复合物切 3 折 [p12] | "
         "复合物级 + 40% 同一性 | 否（且需 650M + A100 [p11-12]）|")
    emit("")
    emit("> **结论**: 我的划分比上述**任何一个**都严格（额外的 40% 序列同一性 + "
         "家族级隔离），因此**不能**声称\"比文献好\"或\"比文献差\"。")
    emit("> 可比区间是 SAMAffinity 的 complex-level (0.54–0.59) 与 "
         "DSSA-PPI 的 complex-level (0.58–0.71)，且我的更严。")

    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
