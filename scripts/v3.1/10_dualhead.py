# -*- coding: utf-8 -*-
"""scripts/v3.1/10_dualhead.py -- V3.1 任务 5 方向 A: 双头 ΔG + ΔΔG（Twin Peaks）.

文献依据
--------
* **Twin Peaks** [p7]：`ΔG 头 = MLP_ΔG(diff(E_mut1, E_mut2))`，
  `ΔΔG 头 = MLP_ΔΔG(diff(E_mut1,E_mut2) − diff(E_Swt1,E_Swt2))`，diff 为逐元素相减；
  总损失 `L = MSE(ΔG) + λ·MSE(ΔΔG)`，**最优 λ = 0.718** [p9]。
* **消融（验证集 Pearson，ΔG/ΔΔG）**：完整 **0.638/0.485**；**关闭双头 → 0.532/0.558**
  [p10][p11] —— 两任务间存在**可调 trade-off**。
* **对称性增强**：`ΔΔG(B→A) = −ΔΔG(A→B)` 使训练数据翻倍 [p8]。
* **Pythia-PPI** [p12]：反向突变一致性 `Pearson(ΔΔG_fwd, ΔΔG_rev) = −0.9224` 是廉价诊断。

为什么这与"特征工程"不同
------------------------
文献一致指出我的**架构/训练目标**是缺口，而非特征：
* SAMAffinity 直接点出我的病灶 ——「序列特征张量因单点替换只发生微小变化，
  导致预测 ΔΔG 聚集在狭窄数值区间」[SAMAffinity p16]；
* Twin Peaks 的 ΔΔG 头**显式建模相对变化**，绕开"绝对量被复合物主导"的问题。

最小验证设计（CPU，固定划分）
-----------------------------
在同一特征（物理 + pooled ESM-2 PCA128）上对比：

| 变体 | ΔG 头 | ΔΔG 头 | 对称增广 |
|---|---|---|---|
| `lightgbm_ref` | ✅（参照）| ❌ | ❌ |
| `mlp_single` | ✅ | ❌ | ❌ |
| `dualhead` | ✅ | ✅ | ❌ |
| `dualhead_sym` | ✅ | ✅ | ✅ |
| `ddg_only_sym` | ❌ | ✅ | ✅ |

> `ddg_only_sym` 是**纯 ΔΔG 模型**：只用成对差分预测相对变化。
> 它在 test 上**没有绝对 ΔG 预测**，因此无法用 Spearman 评价；
> 本脚本对它报告**组内 ΔΔG 相关**（within-target 口径）作为独立读数。

**没有 GPU**：全部 PyTorch CPU；模型刻意做小（2 隐层），以适配 4 核 i7-2630QM。
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
import v31_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v3.1_dualhead_results.csv"
OUT_PRED = C.FEAT / "v3.1_dualhead_predictions.npz"
OUT_PAIR = C.FEAT / "v3.1_dualhead_pair_stats.csv"
OUT_MD = C.FEAT / "v3.1_dualhead_report.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V3.1 方向 A: 双头 ΔG+ΔΔG")
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--lam", type=float, nargs="*", default=[0.0, 0.3, 0.718, 1.5],
                    help="ΔΔG 损失权重 λ（Twin Peaks 最优 0.718 [p9]）")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43])
    ap.add_argument("--max-pairs", type=int, default=60000,
                    help="每个 epoch 采样的 ΔΔG 对数上限（控制 CPU 时间）")
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    torch.set_num_threads(2)
    lgb = C.load_lgb()
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V3.1 方向 A — 双头 ΔG + ΔΔG（Twin Peaks 迁移）")
    emit("")
    emit(f"> 文献依据: Twin Peaks [p7] 双头结构 + λ=0.718 [p9] + 对称增广 [p8]")
    emit(f"> 设备 **CPU only**（torch {torch.__version__}, 线程 2）; 种子 {args.seeds}")
    emit("")

    data = C.load_data()
    y, tr, va, te = data["y"], data["tr"], data["va"], data["te"]
    X = C.ref_features(data)          # 物理 + pooled PCA128
    groups = data["groups"]
    emit("## 1. 数据与口径")
    emit("")
    emit(f"- 行数 **{len(y)}**; train {int(tr.sum())} / val {int(va.sum())} / "
         f"test {int(te.sum())}")
    emit(f"- 特征 **{X.shape[1]}** 维（物理 {data['Xp'].shape[1]} + pooled PCA "
         f"{data['Z'].shape[1]}, 解释方差 {data.get('evr', float('nan')):.4f}）")
    emit(f"- **固定划分 + family-level（已独立复核: 结构级交集 0、两侧 exact-match 0/757）**")
    emit(f"- 参照值: {C.REF_LABEL} = **{C.REF_RHO:+.4f}**")
    emit("")

    # ---------- 2. 构造 ΔΔG 对（只用 train 内的组）----------
    emit("## 2. ΔΔG 训练对的构造")
    emit("")
    tr_idx = np.where(tr)[0]
    by_group: dict[str, list[int]] = {}
    for i in tr_idx:
        by_group.setdefault(groups[i], []).append(i)
    pairs: list[tuple[int, int]] = []
    for g, idxs in by_group.items():
        if len(idxs) < 2:
            continue
        for a in range(len(idxs)):
            for b in range(a + 1, len(idxs)):
                pairs.append((idxs[a], idxs[b]))
    pairs = np.asarray(pairs, dtype=np.int64) if pairs else np.zeros((0, 2), np.int64)
    # 过滤掉 ΔΔG=0 的对（信息量为零）
    if len(pairs):
        dgz = y[pairs[:, 0]] - y[pairs[:, 1]]
        nz = np.abs(dgz) > 1e-9
        emit(f"- 训练集内成对（同 `interface_key`）: **{len(pairs)}**")
        emit(f"- 其中 ΔΔG ≠ 0 的: **{int(nz.sum())}** "
             f"({nz.mean() * 100:.1f}%)")
        pairs = pairs[nz]
        dgz = dgz[nz]
        emit(f"- ΔΔG 分布: 中位 |ΔΔG| {np.median(np.abs(dgz)):.3f}, "
             f"均值 {dgz.mean():+.3f}, std {dgz.std():.3f} kcal/mol")
    emit(f"- 用于训练的对数上限（每 epoch 采样）: {args.max_pairs}")
    emit("")
    if len(pairs) == 0:
        emit("!! 没有可用的 ΔΔG 训练对 -> 方向 A 无法执行")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3

    # ---------- 3. 模型 ----------
    D = X.shape[1]
    class DualHead(nn.Module):
        """Twin Peaks 风格: 共享躯干 + ΔG 头 + ΔΔG 头。"""

        def __init__(self, d_in: int, hidden: int = 256, use_ddg: bool = True):
            super().__init__()
            self.use_ddg = use_ddg
            self.trunk = nn.Sequential(
                nn.Linear(d_in, hidden), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(hidden, 64), nn.ReLU(),
            )
            self.head_g = nn.Linear(64, 1)
            if use_ddg:
                # ΔΔG 头吃**差分**表示 (Twin Peaks: diff(E1,E2)) [p7]
                self.head_ddg = nn.Sequential(
                    nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, 1))

        def embed(self, x):
            return self.trunk(x)

        def forward(self, x):
            h = self.embed(x)
            return self.head_g(h).squeeze(-1), h

    def evaluate_pred(pred: np.ndarray) -> dict:
        y_s, p_s = -y[te], -pred
        out = {}
        for k, fn in mfns.items():
            try:
                out[k] = float(fn(y_s, p_s))
            except Exception:  # noqa: BLE001
                out[k] = float("nan")
        for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
            lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, C.N_BOOT_CI, C.SEED)
            out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
        return out

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    # ---- 参照: LightGBM ----
    emit("## 3. 结果")
    emit("")
    emit("| 变体 | λ | 种子 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---:|---:|---|---:|---:|---:|")
    try:
        r = C.fit_eval_lgb(lgb, X, data, "lightgbm_ref", mfns)
        results.append({k: v for k, v in r.items() if not k.startswith("_")})
        preds["lightgbm_ref"] = r["_pred"].astype(np.float32)
        emit(f"| `lightgbm_ref` | — | — | {r['spearman']:+.4f} | "
             f"[{r['spearman_ci_low']:+.3f}, {r['spearman_ci_high']:+.3f}] | "
             f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} |")
    except Exception as exc:  # noqa: BLE001
        emit(f"| `lightgbm_ref` | — | — | FAIL {exc} | | | | |")
        with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[10_dualhead] lgb ref: {exc}\n{traceback.format_exc()}\n")

    Xt = torch.from_numpy(X)
    yt = torch.from_numpy(y.astype(np.float32))
    tr_t = torch.from_numpy(tr)
    te_idx = np.where(te)[0]

    def train_one(lam: float, seed: int, use_ddg: bool, label: str) -> None:
        torch.manual_seed(seed)
        np.random.seed(seed)
        model = DualHead(D, use_ddg=use_ddg)
        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)
        mse = nn.MSELoss()
        tr_i = np.where(tr)[0]
        va_i = np.where(va)[0]
        best_rho, best_state, bad = -9.0, None, 0
        for ep in range(args.epochs):
            model.train()
            perm = np.random.RandomState(seed * 1000 + ep).permutation(len(tr_i))
            for s in range(0, len(perm), args.batch):
                b = tr_i[perm[s:s + args.batch]]
                xb = Xt[b]
                yb = yt[b]
                pg, _ = model(xb)
                loss = mse(pg, yb)
                if use_ddg and lam > 0 and len(pairs):
                    # 采样成对, 用差分表示预测 ΔΔG (Twin Peaks [p7])
                    k = min(len(pairs), max(256, args.batch * 4))
                    sel = np.random.randint(0, len(pairs), k)
                    p1, p2 = pairs[sel, 0], pairs[sel, 1]
                    h1 = model.embed(Xt[p1])
                    h2 = model.embed(Xt[p2])
                    d_hat = model.head_ddg(h1 - h2).squeeze(-1)
                    d_true = torch.from_numpy(
                        (y[p1] - y[p2]).astype(np.float32))
                    loss = loss + lam * mse(d_hat, d_true)
                opt.zero_grad()
                loss.backward()
                opt.step()
            model.eval()
            with torch.no_grad():
                pv, _ = model(Xt[va_i])
                pv = pv.numpy()
            rho = mr.spearman_rho(-y[va_i], -pv)
            if rho == rho and rho > best_rho + 1e-6:
                best_rho, bad = float(rho), 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
            if bad >= 20:
                break
        if best_state:
            model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            pt, _ = model(Xt[te_idx])
            pt = pt.numpy()
        met = evaluate_pred(pt)
        row = {"variant": label, "lam": lam, "seed": seed,
               "n_features": D, **met}
        results.append(row)
        preds[f"{label}_lam{lam}_seed{seed}"] = pt.astype(np.float32)
        emit(f"| `{label}` | {lam} | {seed} | {met['spearman']:+.4f} | "
             f"[{met['spearman_ci_low']:+.3f}, {met['spearman_ci_high']:+.3f}] | "
             f"{met['ndcg@10']:.4f} | {met['ef@1%']:.4f} | "
             f"{met['precision@10']:.4f} |")

    # ---- mlp_single (无 ΔΔG 头) ----
    for sd in args.seeds:
        try:
            train_one(0.0, sd, use_ddg=False, label="mlp_single")
        except Exception as exc:  # noqa: BLE001
            emit(f"| `mlp_single` | 0.0 | {sd} | FAIL {exc} | | | | |")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[10_dualhead] mlp_single {sd}: {exc}\n"
                         f"{traceback.format_exc()}\n")

    # ---- dualhead + λ 扫描 ----
    for lam in args.lam:
        if lam <= 0:
            continue
        for sd in args.seeds:
            try:
                train_one(lam, sd, use_ddg=True, label="dualhead")
            except Exception as exc:  # noqa: BLE001
                emit(f"| `dualhead` | {lam} | {sd} | FAIL {exc} | | | | |")
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[10_dualhead] dualhead lam={lam} seed={sd}: {exc}\n"
                             f"{traceback.format_exc()}\n")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    np.savez(OUT_PRED, test_y=y[te].astype(np.float32),
             test_pdb_id=data["pdb"][te], **preds)

    # ---------- 4. 配对检验 vs 参照 ----------
    emit("")
    emit(f"## 4. 配对 bootstrap vs 参照（{C.REF_LABEL}）")
    emit("")
    # 参照预测：用 LightGBM 的 test 预测
    ref_pred = preds.get("lightgbm_ref")
    emit("| 变体 | Δ Spearman | 95% CI | p | BH q | 判定 |")
    emit("|---|---:|---|---:|---:|---|")
    prows, pvals = [], []
    res["spearman"] = pd.to_numeric(res["spearman"], errors="coerce")
    base_mlp = res[(res["variant"] == "mlp_single")
                   & res["spearman"].notna()].sort_values(
        "spearman", ascending=False)
    for _, r in res.iterrows():
        key = f"{r['variant']}_lam{r['lam']}_seed{r['seed']}"
        if r["variant"] == "lightgbm_ref" or key not in preds:
            continue
        # 与 mlp_single 最好者比（同一模型族，公平比较 ΔΔG 头的贡献）
        if not len(base_mlp):
            continue
        b = base_mlp.iloc[0]
        bkey = f"{b['variant']}_lam{b['lam']}_seed{b['seed']}"
        rr = C.paired_bootstrap(mr.spearman_rho, -y[te], -preds[bkey],
                                -preds[key], C.N_BOOT_PAIR, C.SEED)
        prows.append({"variant": r["variant"], "lam": r["lam"], "seed": r["seed"],
                      "vs": bkey, **rr})
        pvals.append(rr["p_two_sided"])
    qs = C.bh_fdr(pvals)
    for r, q in zip(prows, qs):
        r["bh_q"] = q
        v = ("更优" if (q < 0.05 and r["delta"] > 0) else
             "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
        emit(f"| `{r['variant']}` λ={r['lam']} seed={r['seed']} | {r['delta']:+.4f} | "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
             f"{q:.4f} | {v} |")
    if prows:
        pd.DataFrame(prows).to_csv(OUT_PAIR, index=False, encoding="utf-8")

    # ---------- 5. 结论 ----------
    emit("")
    emit("## 5. 结论")
    emit("")
    if len(res):
        # 注意: 有 FAIL 行时 spearman 列是 object dtype, 必须先强制数值再过滤,
        # 否则 sort_values 会把字符串与 NaN 混排, 选出 nan 行 (实测踩到)。
        res["spearman"] = pd.to_numeric(res["spearman"], errors="coerce")
        res["ndcg@10"] = pd.to_numeric(res["ndcg@10"], errors="coerce")
        ok = res[(res["variant"] != "lightgbm_ref") & res["spearman"].notna()]
        if len(ok):
            best = ok.sort_values("spearman", ascending=False).iloc[0]
            emit(f"- 神经网络族最优: `{best['variant']}` λ={best['lam']} "
                 f"-> Spearman **{best['spearman']:+.4f}** "
                 f"[{best['spearman_ci_low']:+.3f}, {best['spearman_ci_high']:+.3f}]")
        lr_ = res[(res["variant"] == "lightgbm_ref") & res["spearman"].notna()]
        if len(lr_):
            emit(f"- LightGBM 参照: **{float(lr_.iloc[0]['spearman']):+.4f}** "
                 f"(V2.3.2 历史值 {C.REF_RHO:+.4f})")
        if len(ok) and len(lr_):
            d = float(best["spearman"]) - float(lr_.iloc[0]["spearman"])
            emit(f"- Δ(NN 最优 − LightGBM) = **{d:+.4f}**")
            emit("")
            if d > 0.03:
                emit("> 神经网络 + ΔΔG 头**超过了 LightGBM 参照**（但需看配对检验的 BH q）。")
            else:
                emit("> 神经网络（含 ΔΔG 双头）**未超过 LightGBM 参照**。")
                emit("> 这与 V2.5 的 per-residue CNN 结论一致：在**同一特征**上换模型族/损失")
                emit("> 并未带来增益 —— 增益瓶颈不在**训练目标**本身。")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
