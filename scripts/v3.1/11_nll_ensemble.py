# -*- coding: utf-8 -*-
"""scripts/v3.1/11_nll_ensemble.py -- V3.1 任务 5 方向 C: NLL + 深度集成 + 不确定性分层.

文献依据
--------
* **CatPred** [p4][p5][p14]：用 **负对数似然（NLL）** 对 log10 标签拟合高斯分布，
  网络同时输出 μ 与 σ²；**集成 = 同一架构 10 个随机初始化副本** [p5]。
* **不确定性分解**：aleatoric = 集成各成员预测方差的算术平均再开方；
  epistemic = 集成均值的样本标准差；总不确定性 = 两者之和 [p10][p11]。
* **已验证的收益（这是关键）**：按不确定性分箱后指标**单调改善**，
  held-out **25 百分位**子集 R² 提升到 **0.78 / 0.76 / 0.61**（kcat/Km/Ki），
  MAE 相对 100 百分位下降约 **46% / 37% / 29%** [p11]。
* **为什么这对我的项目重要**：我的 EF@1%/P@10（头部富集）与全局 Spearman 存在张力。
  若不确定性分层有效，就能**只在高置信子集里做 top-K 挑选**，把"挑 10 个克隆"
  的成功率提上去 —— 这是**训练策略**层面的改进，不是特征。

最小验证设计（CPU，固定划分）
-----------------------------
对比（同一特征：物理 + pooled ESM-2 PCA128）：

| 变体 | 损失 | 集成 | 用途 |
|---|---|---|---|
| `lgb_ref` | L2（MSE）| ❌ | 参照 |
| `nn_mse_ens` | MSE | ✅ | 隔离"集成"本身的贡献 |
| `nn_nll_ens` | NLL（μ,σ²）| ✅ | 隔离"NLL 概率头"的贡献 |

然后对**最优变体**做不确定性分层：
* 按预测 σ 分位把 test 切成 4 档，报告每档的 Spearman / EF@1% / P@10；
* 检验是否**单调**（CatPred 报告的单调性是我要复现的关键结论）[p11]；
* 报告"最高置信 25% 子集"的 Spearman 相对全 test 的变化。

**注意**: 本任务的判据是"是否**显著优于 0.3890**（BH q<0.05）"。
不确定性分层**不会**改变全 test 的 Spearman，因此分层本身不能"胜出"；
它的价值在于**可用性**（能否挑出可信的子集）。这一点在报告中如实说明。
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

OUT_CSV = C.FEAT / "v3.1_nll_results.csv"
OUT_BIN = C.FEAT / "v3.1_nll_uncertainty_bins.csv"
OUT_PRED = C.FEAT / "v3.1_nll_predictions.npz"
OUT_MD = C.FEAT / "v3.1_nll_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V3.1 方向 C: NLL + 深度集成")
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--ens", type=int, default=8, help="集成成员数（CatPred 用 10 [p5]）")
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43, 44, 45])
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    torch.set_num_threads(2)
    lgb = C.load_lgb()
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V3.1 方向 C — NLL 概率回归 + 深度集成 + 不确定性分层（CatPred 迁移）")
    emit("")
    emit(f"> 文献依据: NLL 高斯头 [CatPred p4,p5,p14]；集成 = 多随机初始化 [p5]；")
    emit(f"> 不确定性分箱单调改善 [CatPred p11]")
    emit(f"> 设备 **CPU only**; 集成成员 {args.ens}; 种子 {args.seeds}")
    emit("")

    data = C.load_data()
    y, tr, va, te = data["y"], data["tr"], data["va"], data["te"]
    X = C.ref_features(data)
    emit("## 1. 数据与口径")
    emit("")
    emit(f"- 行数 {len(y)}; train {int(tr.sum())} / val {int(va.sum())} / "
         f"test {int(te.sum())}; 特征 {X.shape[1]} 维")
    emit(f"- 固定划分 + family-level（已独立复核）; 参照 {C.REF_LABEL} = "
         f"**{C.REF_RHO:+.4f}**")
    emit("")

    D = X.shape[1]
    Xt = torch.from_numpy(X)
    yt = torch.from_numpy(y.astype(np.float32))
    tr_i, va_i, te_i = np.where(tr)[0], np.where(va)[0], np.where(te)[0]

    class GaussNet(nn.Module):
        """输出 μ 与 log σ²（NLL 头）[CatPred p4]。"""

        def __init__(self, d_in: int, hidden: int = 256, nll: bool = True):
            super().__init__()
            self.nll = nll
            self.body = nn.Sequential(
                nn.Linear(d_in, hidden), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(hidden, 64), nn.ReLU(),
            )
            self.mu = nn.Linear(64, 1)
            self.logvar = nn.Linear(64, 1) if nll else None

        def forward(self, x):
            h = self.body(x)
            mu = self.mu(h).squeeze(-1)
            if self.nll:
                lv = self.logvar(h).squeeze(-1).clamp(-8.0, 8.0)
                return mu, lv
            return mu, None

    def nll_loss(mu, lv, target):
        # 0.5*(log σ² + (y-μ)²/σ²) 的高斯 NLL（去掉常数项）
        return (0.5 * (lv + (target - mu) ** 2 / torch.exp(lv))).mean()

    def train_net(seed: int, nll: bool):
        torch.manual_seed(seed)
        np.random.seed(seed)
        m = GaussNet(D, nll=nll)
        opt = torch.optim.Adam(m.parameters(), lr=args.lr, weight_decay=args.wd)
        best_rho, best_state, bad = -9.0, None, 0
        for ep in range(args.epochs):
            m.train()
            perm = np.random.RandomState(seed * 1000 + ep).permutation(len(tr_i))
            for s in range(0, len(perm), args.batch):
                b = tr_i[perm[s:s + args.batch]]
                mu, lv = m(Xt[b])
                if nll:
                    loss = nll_loss(mu, lv, yt[b])
                else:
                    loss = nn.functional.mse_loss(mu, yt[b])
                opt.zero_grad()
                loss.backward()
                opt.step()
            m.eval()
            with torch.no_grad():
                vmu, _ = m(Xt[va_i])
                vmu = vmu.numpy()
            rho = mr.spearman_rho(-y[va_i], -vmu)
            if rho == rho and rho > best_rho + 1e-6:
                best_rho, bad = float(rho), 0
                best_state = {k: v.detach().clone() for k, v in m.state_dict().items()}
            else:
                bad += 1
            if bad >= 20:
                break
        if best_state:
            m.load_state_dict(best_state)
        m.eval()
        with torch.no_grad():
            tmu, tlv = m(Xt[te_i])
            tmu = tmu.numpy()
            tsd = (torch.exp(0.5 * tlv).numpy() if tlv is not None
                   else np.zeros_like(tmu))
        return tmu, tsd, best_rho

    def run_variant(nll: bool, label: str) -> dict | None:
        mus, sds, rhos = [], [], []
        for k in range(args.ens):
            sd = args.seeds[k % len(args.seeds)] + 100 * (k // len(args.seeds))
            try:
                mu, sd_, rv = train_net(sd, nll)
                mus.append(mu)
                sds.append(sd_)
                rhos.append(rv)
            except Exception as exc:  # noqa: BLE001
                emit(f"  !! {label} 成员 seed={sd} 失败: {exc}")
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[11_nll_ensemble] {label} seed={sd}: {exc}\n"
                             f"{traceback.format_exc()}\n")
        if not mus:
            return None
        MU = np.stack(mus)                       # (M, n_te)
        SD = np.stack(sds)
        mu_ens = MU.mean(axis=0)
        # 不确定性分解 [CatPred p10-p11]
        alea = np.sqrt((SD ** 2).mean(axis=0))
        epis = MU.std(axis=0)
        tot = np.sqrt(alea ** 2 + epis ** 2)
        y_s, p_s = -y[te], -mu_ens
        met = {}
        for k, fn in mfns.items():
            try:
                met[k] = float(fn(y_s, p_s))
            except Exception:  # noqa: BLE001
                met[k] = float("nan")
        for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
            lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, C.N_BOOT_CI, C.SEED)
            met[f"{k}_ci_low"], met[f"{k}_ci_high"] = lo, hi
        emit(f"| `{label}` | {len(mus)} | {met['spearman']:+.4f} | "
             f"[{met['spearman_ci_low']:+.3f}, {met['spearman_ci_high']:+.3f}] | "
             f"{met['ndcg@10']:.4f} | {met['ef@1%']:.4f} | "
             f"{met['precision@10']:.4f} | "
             f"{np.nanmean(rhos):+.4f} | {np.nanmean(alea):.3f} |")
        return {"label": label, "metrics": met, "mu": mu_ens, "alea": alea,
                "epis": epis, "tot": tot, "member_rho_val": float(np.nanmean(rhos))}

    emit("## 2. 结果")
    emit("")
    emit("| 变体 | 成员数 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 | 成员 val ρ 均值 | 平均 aleatoric σ |")
    emit("|---|---:|---:|---|---:|---:|---:|---:|---:|")

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    # ---- LightGBM 参照 ----
    try:
        r = C.fit_eval_lgb(lgb, X, data, "lgb_ref", mfns)
        results.append({"label": "lgb_ref", **{k: v for k, v in r.items()
                                               if not k.startswith("_")}})
        preds["lgb_ref"] = r["_pred"].astype(np.float32)
        emit(f"| `lgb_ref` | — | {r['spearman']:+.4f} | "
             f"[{r['spearman_ci_low']:+.3f}, {r['spearman_ci_high']:+.3f}] | "
             f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} | — | — |")
    except Exception as exc:  # noqa: BLE001
        emit(f"| `lgb_ref` | — | FAIL {exc} | | | | | | |")

    out_mse = run_variant(False, "nn_mse_ens")
    out_nll = run_variant(True, "nn_nll_ens")
    for o in (out_mse, out_nll):
        if o is not None:
            results.append({"label": o["label"], "n_members": args.ens,
                            "member_rho_val": o["member_rho_val"],
                            **o["metrics"]})
            preds[o["label"]] = o["mu"].astype(np.float32)

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    if preds:
        save = {"test_y": y[te].astype(np.float32),
                "test_pdb_id": data["pdb"][te]}
        save.update({k: v.astype(np.float32) for k, v in preds.items()})
        if out_nll is not None:
            for k in ("alea", "epis", "tot"):
                save[f"nll_{k}"] = out_nll[k].astype(np.float32)
        np.savez(OUT_PRED, **save)

    # ---------- 3. 不确定性分层（核心：复现 CatPred 的单调性）----------
    emit("")
    emit("## 3. 不确定性分层（CatPred 结论复现 [p11]）")
    emit("")
    bins_rows: list[dict] = []
    src = out_nll if out_nll is not None else out_mse
    if src is not None:
        tot = src["tot"] if out_nll is not None else src["epis"]
        emit(f"### 按**总不确定性**分位切 4 档（用 `{src['label']}` 的预测）")
        emit("")
        emit("| 档位 | 不确定性子集 | n | Spearman | NDCG@10 | EF@1% | P@10 |")
        emit("|---|---|---:|---:|---:|---:|---:|")
        qs = np.quantile(tot, [0, 0.25, 0.5, 0.75, 1.0])
        qs[0] -= 1e-9; qs[-1] += 1e-9
        labels = ["Q1(最自信)", "Q2", "Q3", "Q4(最不自信)"]
        prev = None
        mono = True
        for i in range(4):
            msk = (tot > qs[i]) & (tot <= qs[i + 1])
            if msk.sum() < 10:
                emit(f"| {labels[i]} | [{qs[i]:.3f}, {qs[i+1]:.3f}] | {int(msk.sum())} "
                     f"| (n<10 略) | | | |")
                continue
            ys, ps = -y[te][msk], -src["mu"][msk]
            rho = mr.spearman_rho(ys, ps)
            nd = mr.ndcg_at_k(ys, ps, 10, top_frac=C.TOP_FRAC)
            ef = mr.enrichment_factor(ys, ps, 1.0, C.TOP_FRAC)
            p10 = mr.precision_at_k(ys, ps, 10, C.TOP_FRAC)
            emit(f"| {labels[i]} | [{qs[i]:.3f}, {qs[i+1]:.3f}] | {int(msk.sum())} | "
                 f"{rho:+.4f} | {nd:.4f} | {ef:.4f} | {p10:.4f} |")
            bins_rows.append({"scheme": "quartile", "bin": labels[i],
                              "n": int(msk.sum()), "unc_lo": float(qs[i]),
                              "unc_hi": float(qs[i + 1]), "spearman": float(rho),
                              "ndcg@10": float(nd), "ef@1%": float(ef),
                              "precision@10": float(p10)})
            if prev is not None and rho < prev - 1e-9:
                mono = False
            prev = rho
        # 最高置信 25%
        m25 = tot <= np.quantile(tot, 0.25)
        if m25.sum() >= 10:
            ys, ps = -y[te][m25], -src["mu"][m25]
            rho25 = mr.spearman_rho(ys, ps)
            full = src["metrics"]["spearman"]
            emit("")
            emit(f"- **最高置信 25% 子集**: n={int(m25.sum())}, Spearman = "
                 f"**{rho25:+.4f}** vs 全 test {full:+.4f} -> Δ = **{rho25 - full:+.4f}**")
            bins_rows.append({"scheme": "top25pct", "bin": "most_confident",
                              "n": int(m25.sum()), "spearman": float(rho25),
                              "delta_vs_full": float(rho25 - full)})
        emit("")
        emit(f"- **单调性**: {'✅ 成立（不确定性越低指标越好）' if mono else '❌ 不成立'}")
        emit("")
        if mono:
            emit("> ✅ **复现了 CatPred 的关键结论** [p11]：按不确定性分层后指标单调改善。")
            emit("> 实用含义：可以对候选**先按 σ 过滤**，只在最自信的子集里做 top-K 挑选。")
        else:
            emit("> ❌ 本数据上**未复现**单调性 —— 诚实报告，不宣称该收益。")
    else:
        emit("- 没有可用的神经网络结果，跳过分层分析")
    if bins_rows:
        pd.DataFrame(bins_rows).to_csv(OUT_BIN, index=False, encoding="utf-8")

    # ---------- 4. 结论 ----------
    emit("")
    emit("## 4. 结论")
    emit("")
    if len(res):
        emit("详见表。要点：")
        emit("")
        for lab in ("lgb_ref", "nn_mse_ens", "nn_nll_ens"):
            sub = res[res["label"] == lab]
            if len(sub) and np.isfinite(pd.to_numeric(sub["spearman"],
                                                      errors="coerce")).any():
                emit(f"- `{lab}`: Spearman {float(sub['spearman'].iloc[0]):+.4f}")
        emit("")
        emit("> **判据提醒**：本任务的胜出标准是「显著优于 0.3890（BH q<0.05）」。")
        emit("> 不确定性分层**不改变全 test 的 Spearman**，因此它**不能**'胜出'；")
        emit("> 它的价值是**可用性**（能否挑出可信子集），属训练策略而非精度提升。")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
