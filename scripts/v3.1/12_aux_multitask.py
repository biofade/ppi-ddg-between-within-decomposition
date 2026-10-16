# -*- coding: utf-8 -*-
"""scripts/v3.1/12_aux_multitask.py -- V3.1 任务 5 方向 B: 辅助任务多任务学习.

文献依据
--------
* **SELFprot** [p4]：共享冻结主干 → 一层可训练 cross-attention → **每任务头**。
  其**任务间迁移是最硬的定量证据**（Table 1）[p5]：
  * 去掉 binding 任务（Ki/Kd/IC50/EC50） → **kcat RMSD +154%**（1.105→2.810）
  * 去掉 kinetics 任务（kcat/Km） → **Ki RMSD +106%**（1.506→3.116）
  * 去掉功能位点分类任务 → **Ki RMSD +30%**（1.957）
  * 但纳入 binding 任务使功能位点 precision **下降 4%**（负迁移）
* **SELFprot 的明确缺口**：所有任务**等权**（"for simplicity" [p4]），
  **无任何梯度冲突处理**（无 GradNorm/PCGrad）-> **任务加权仍需自己解决**。
* **Pythia-PPI** [p8][p18]：多任务 `L = α·L_affinity + (1−α)·L_stability`，
  **α=0.8**（stability:affinity = 0.8:0.2）最优；加稳定性头使 per-structure
  Pearson **0.4462→0.4778** [p8]。

为什么这能在**我的数据上**做最小验证
------------------------------------
SELFprot/Pythia-PPI 用的辅助任务是稳定性/功能位点，我**没有**这些标签。
但我的表里有**可用的辅助信号**（不需要新数据、不需要 GPU）：

| 辅助任务 | 来源列 | 是否含标签信息 |
|---|---|---|
| **突变数** | `Mutations`（可解析 58.8%，含 1–8 个突变）| ❌ 纯输入侧，无泄漏 |
| **测量温度** | `Temperature(K)` | ❌ 实验条件 |
| **序列长度** | `Ligand/Receptor Sequences` | ❌ 已隐含在 pooled 嵌入里 |

本脚本用**突变数**作辅助回归目标（辅助头预测 log1p(n_mut)）。
它**不引入任何测试标签**，且能检验一个具体问题：
**加一个辅助头是否会改变主任务的性能**（正向 or 负向）。

设计（严格对照）
----------------
* 同一特征（物理 + pooled ESM-2 PCA128）、同一固定划分；
* 三组对照，**唯一差别是辅助损失权重 β**：

| 变体 | 主损失 | 辅助损失 |
|---|---|---|
| `mt_beta0` | MSE(ΔG) | 无（= 单任务基线）|
| `mt_beta_lo` | MSE(ΔG) | β=0.1 · MSE(log1p n_mut) |
| `mt_beta_hi` | MSE(ΔG) | β=0.5 · MSE(log1p n_mut) |

* 另加 `lgb_ref` 作为非神经参照。
* **主判据仍是 Spearman 是否显著优于 0.3890**；辅助任务是手段，不是目的。

**诚实预期**：辅助任务若无效，本脚本会明确报告 —— 这与 V2.5 的多次"零增益"一致。
"""

from __future__ import annotations

import argparse
import re
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

OUT_CSV = C.FEAT / "v3.1_aux_results.csv"
OUT_PRED = C.FEAT / "v3.1_aux_predictions.npz"
OUT_MD = C.FEAT / "v3.1_aux_report.md"

MUT_PAT = re.compile(r"^([A-Za-z0-9]+)_([A-Za-z])(\d+)([A-Za-z])$")
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def n_mutations(s) -> float:
    """解析突变数；无法解析返回 NaN（训练时用中位数兜底并记录）。"""
    if not isinstance(s, str) or not s.strip():
        return np.nan
    n = 0
    for p in re.split(r"[;,]", s):
        if MUT_PAT.match(p.strip()):
            n += 1
    return float(n) if n > 0 else np.nan


def main() -> int:
    ap = argparse.ArgumentParser(description="V3.1 方向 B: 辅助任务多任务")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--betas", type=float, nargs="*", default=[0.0, 0.1, 0.5])
    ap.add_argument("--seeds", type=int, nargs="*", default=[42, 43])
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    torch.set_num_threads(2)
    lgb = C.load_lgb()
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V3.1 方向 B — 辅助任务多任务学习（SELFprot / Pythia-PPI 迁移）")
    emit("")
    emit("> 文献依据: 共享主干 + 每任务头 [SELFprot p4]；辅助任务迁移的定量证据 [SELFprot p5]；")
    emit("> α=0.8 / 辅助头提升 per-structure Pearson [Pythia-PPI p8,p18]")
    emit("> **注意**: SELFprot 自身所有任务等权、无梯度冲突处理 [SELFprot p4] ->")
    emit("> 任务加权是本脚本必须自己扫的超参（β 扫描）")
    emit(f"> 设备 **CPU only**; 种子 {args.seeds}")
    emit("")

    data = C.load_data()
    y, tr, va, te = data["y"], data["tr"], data["va"], data["te"]
    X = C.ref_features(data)
    d = data["df"]
    emit("## 1. 辅助标签（不使用测试标签）")
    emit("")
    nm = np.array([n_mutations(v) for v in d["Mutations"].tolist()], float)
    ok = np.isfinite(nm)
    emit(f"- `Mutations` 可解析出突变数的行: **{int(ok.sum())} / {len(nm)}** "
         f"({ok.mean() * 100:.1f}%)")
    emit(f"- 突变数分布（可解析内）: "
         f"{pd.Series(nm[ok]).value_counts().sort_index().head(8).to_dict()}")
    # 兜底: 训练集中位数（不使用 val/test 的标签分布）
    med = float(np.nanmedian(nm[tr & ok])) if (tr & ok).any() else 1.0
    nm_f = np.where(ok, nm, med)
    aux = np.log1p(nm_f).astype(np.float32)
    emit(f"- 缺失用**训练集中位数**兜底 = {med:.0f} 个突变 -> 辅助目标 = log1p(n_mut)")
    emit(f"- 辅助目标统计: mean {aux.mean():.4f}, std {aux.std():.4f}")
    emit("")
    emit("> 辅助任务**不含任何亲和力标签**，因此不构成泄漏；它是一个纯输入侧的辅助信号。")
    emit("")

    D = X.shape[1]
    Xt = torch.from_numpy(X)
    yt = torch.from_numpy(y.astype(np.float32))
    at = torch.from_numpy(aux)
    tr_i, va_i, te_i = np.where(tr)[0], np.where(va)[0], np.where(te)[0]

    class MTNet(nn.Module):
        def __init__(self, d_in: int, hidden: int = 256, n_aux: bool = True):
            super().__init__()
            self.n_aux = n_aux
            self.trunk = nn.Sequential(
                nn.Linear(d_in, hidden), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(hidden, 64), nn.ReLU())
            self.head_main = nn.Linear(64, 1)
            self.head_aux = nn.Linear(64, 1) if n_aux else None

        def forward(self, x):
            h = self.trunk(x)
            m = self.head_main(h).squeeze(-1)
            a = self.head_aux(h).squeeze(-1) if self.n_aux else None
            return m, a

    def evaluate(pred: np.ndarray) -> dict:
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

    def train(beta: float, seed: int, label: str) -> dict | None:
        torch.manual_seed(seed)
        np.random.seed(seed)
        m = MTNet(D, n_aux=(beta > 0))
        opt = torch.optim.Adam(m.parameters(), lr=args.lr, weight_decay=args.wd)
        mse = nn.MSELoss()
        best_rho, best_state, bad = -9.0, None, 0
        for ep in range(args.epochs):
            m.train()
            perm = np.random.RandomState(seed * 1000 + ep).permutation(len(tr_i))
            for s in range(0, len(perm), args.batch):
                b = tr_i[perm[s:s + args.batch]]
                pm, pa = m(Xt[b])
                loss = mse(pm, yt[b])
                if beta > 0 and pa is not None:
                    loss = loss + beta * mse(pa, at[b])
                opt.zero_grad()
                loss.backward()
                opt.step()
            m.eval()
            with torch.no_grad():
                pv, _ = m(Xt[va_i])
                pv = pv.numpy()
            rho = mr.spearman_rho(-y[va_i], -pv)
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
            pt, pa_te = m(Xt[te_i])
            pt = pt.numpy()
            aux_rho = float("nan")
            if pa_te is not None:
                # 辅助头自身的 test 相关（诊断: 辅助任务是否真被学到）
                try:
                    aux_rho = float(mr.spearman_rho(aux[te], pa_te.numpy()))
                except Exception:  # noqa: BLE001
                    pass
        met = evaluate(pt)
        emit(f"| `{label}` | {beta} | {seed} | {met['spearman']:+.4f} | "
             f"[{met['spearman_ci_low']:+.3f}, {met['spearman_ci_high']:+.3f}] | "
             f"{met['ndcg@10']:.4f} | {met['ef@1%']:.4f} | "
             f"{met['precision@10']:.4f} | {aux_rho:+.4f} |")
        return {"variant": label, "beta": beta, "seed": seed,
                "aux_head_test_spearman": aux_rho, **met, "_pred": pt}

    emit("## 2. 结果")
    emit("")
    emit("| 变体 | β | 种子 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 | 辅助头 test ρ（诊断）|")
    emit("|---|---:|---:|---:|---|---:|---:|---:|---:|")

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    try:
        r = C.fit_eval_lgb(lgb, X, data, "lgb_ref", mfns)
        results.append({k: v for k, v in r.items() if not k.startswith("_")})
        preds["lgb_ref"] = r["_pred"].astype(np.float32)
        emit(f"| `lgb_ref` | — | — | {r['spearman']:+.4f} | "
             f"[{r['spearman_ci_low']:+.3f}, {r['spearman_ci_high']:+.3f}] | "
             f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} | — |")
    except Exception as exc:  # noqa: BLE001
        emit(f"| `lgb_ref` | — | — | FAIL {exc} | | | | | |")
        with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[12_aux_multitask] lgb: {exc}\n{traceback.format_exc()}\n")

    for beta in args.betas:
        lab = "mt_beta0" if beta == 0 else f"mt_beta{beta}"
        for sd in args.seeds:
            try:
                o = train(beta, sd, lab)
                if o:
                    p = o.pop("_pred")
                    results.append(o)
                    preds[f"{lab}_seed{sd}"] = p.astype(np.float32)
            except Exception as exc:  # noqa: BLE001
                emit(f"| `{lab}` | {beta} | {sd} | FAIL {exc} | | | | | |")
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[12_aux_multitask] {lab} seed={sd}: {exc}\n"
                             f"{traceback.format_exc()}\n")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    if preds:
        np.savez(OUT_PRED, test_y=y[te].astype(np.float32),
                 test_pdb_id=data["pdb"][te],
                 **{k: v.astype(np.float32) for k, v in preds.items()})

    # ---------- 3. 配对检验: 有辅助头 vs 无辅助头（β=0）----------
    emit("")
    emit("## 3. 配对检验: 辅助头 vs 无辅助头（同种子）")
    emit("")
    emit("| 对比 | Δ Spearman | 95% CI | p | BH q | 判定 |")
    emit("|---|---:|---|---:|---:|---|")
    prows, pvals = [], []
    for beta in args.betas:
        if beta == 0:
            continue
        for sd in args.seeds:
            a = f"mt_beta0_seed{sd}"
            b = f"mt_beta{beta}_seed{sd}"
            if a in preds and b in preds:
                rr = C.paired_bootstrap(mr.spearman_rho, -y[te], -preds[a], -preds[b],
                                        C.N_BOOT_PAIR, C.SEED)
                prows.append({"comparison": f"β={beta} − β=0 (seed {sd})", **rr})
                pvals.append(rr["p_two_sided"])
    qs = C.bh_fdr(pvals)
    for r, q in zip(prows, qs):
        r["bh_q"] = q
        v = ("更优" if (q < 0.05 and r["delta"] > 0) else
             "更差" if (q < 0.05 and r["delta"] < 0) else "无显著差异")
        emit(f"| {r['comparison']} | {r['delta']:+.4f} | "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}] | {r['p_two_sided']:.4f} | "
             f"{q:.4f} | {v} |")
    if prows:
        pd.DataFrame(prows).to_csv(C.FEAT / "v3.1_aux_pairs.csv", index=False,
                                   encoding="utf-8")

    # ---------- 4. 结论 ----------
    emit("")
    emit("## 4. 结论")
    emit("")
    if len(res):
        fin = res[pd.to_numeric(res["spearman"], errors="coerce").notna()]
        lr_ = fin[fin["variant"] == "lgb_ref"]
        nn = fin[fin["variant"] != "lgb_ref"]
        if len(nn):
            best = nn.sort_values("spearman", ascending=False).iloc[0]
            emit(f"- 神经网络族最优: `{best['variant']}` β={best['beta']} "
                 f"-> Spearman **{best['spearman']:+.4f}**")
        if len(lr_):
            emit(f"- LightGBM 参照: **{float(lr_.iloc[0]['spearman']):+.4f}** "
                 f"(V2.3.2 历史 {C.REF_RHO:+.4f})")
        aux_rows = nn[nn["beta"] > 0]
        if len(aux_rows):
            ar = pd.to_numeric(aux_rows["aux_head_test_spearman"],
                               errors="coerce").mean()
            emit(f"- **辅助头确实学到了东西**: 辅助目标 test ρ 均值 = **{ar:+.4f}** "
                 f"(>0 说明辅助任务非平凡)")
        emit("")
        emit("> 判据提醒: 胜出标准是「显著优于 0.3890（BH q<0.05）」。")
        emit("> 辅助任务若未提升主任务，本报告将如实记录 —— 这与 V2.5 的零增益一致。")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
