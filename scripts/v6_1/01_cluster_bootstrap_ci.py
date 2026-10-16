# -*- coding: utf-8 -*-
"""scripts/v6_1/01_cluster_bootstrap_ci.py -- V6.1 PHASE 1.1: cluster bootstrap CI.

对 CATH held-out（21 complex）上 5 个模型预测列，用 **cluster bootstrap**
（重采样单位 = complex_id）估计 `Macro per-PDB Spearman` 与 `Efficiency@10`
的 95% CI。

规格（V6.1 PHASE 1.1）
----------------------
* 数据源：`standardized_preds.csv` 中 `dataset_name == 'CATH'`
* 复合物范围：21 个（N≥5 过滤后）
* 重采样单位：complex_id
* 抽样次数：**5000**；seed = 42
* 若某次迭代 unique complex < 15 → 跳过该次
* 报告实际有效次数
* 模型：Flex-ddG / FoldX / DDGPred / ESM-IF / proteinMPNN
* 指标：Macro per-PDB Spearman, Efficiency@10

输出
----
* `experiments/v6_1/phase1_ci/bootstrap_ci_models.csv`
* `experiments/v6_1/phase1_ci/bootstrap_ci_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402  （口径一致、快 17 倍的实现）

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V61 = C.ROOT / "experiments" / "v6_1"
PH1 = V61 / "phase1_ci"
PH1.mkdir(parents=True, exist_ok=True)
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
OUT_CSV = PH1 / "bootstrap_ci_models.csv"
OUT_MD = PH1 / "bootstrap_ci_report.md"

N_BOOT = 5000
SEED = 42
MIN_CX = 15
MODELS = ["Flex-ddG", "FoldX", "DDGPred", "ESM-IF", "proteinMPNN"]
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def macro_spearman(y, p, g) -> float:
    return FM.macro_spearman(y, p, g, min_n=C.MIN_PS_N)


N_RAND_DRAWS = 200
N_RAND_BOOT = 20      # bootstrap 内部用较少随机抽样（已标注为近似）


def eff10(y, p, g, n_rand: int = N_RAND_DRAWS) -> float:
    """Efficiency@10（委托给 `fastmetrics`，与 `v5_common` 同公式）。"""
    return FM.eff10(y, p, g, n_rand=n_rand)


def boot(fn, y, p, g, n_boot=N_BOOT, seed=SEED, min_cx=MIN_CX):
    """cluster bootstrap：抽 complex 并携带其全部 rows。"""
    uniq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uniq]
    pt = float(fn(y, p, g))
    rs = np.random.RandomState(seed)
    vals, used, skipped = [], 0, 0
    for _ in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        u = len(np.unique(pick))
        if u < min_cx:
            skipped += 1
            continue
        rows = np.concatenate([cidx[j] for j in pick])
        # 重采样后 complex 标签重编号（同一 complex 被抽多次 = 多个副本）
        rep = np.concatenate([np.full(len(cidx[j]), i, dtype=np.int64)
                              for i, j in enumerate(pick)])
        if fn is eff10:
            v = fn(y[rows], p[rows], rep, n_rand=N_RAND_BOOT)
        else:
            v = fn(y[rows], p[rows], rep)
        if v == v:
            vals.append(v)
            used += 1
    if len(vals) < 20:
        return {"point": pt, "ci_low": float("nan"), "ci_high": float("nan"),
                "sd": float("nan"), "n_valid": len(vals), "n_skipped": skipped,
                "p_gt_zero": float("nan")}
    a = np.asarray(vals, float)
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": pt, "ci_low": float(lo), "ci_high": float(hi),
            "sd": float(a.std(ddof=1)), "n_valid": used, "n_skipped": skipped,
            "p_gt_zero": float((a > 0).mean())}


def main() -> int:
    t0 = time.time()
    emit("# V6.1 PHASE 1.1 — Cluster Bootstrap CI（CATH held-out）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 数据源: `experiments/v6_sota_audit/standardized_preds.csv` "
         f"(`dataset_name == 'CATH'`)")
    emit(f"- 重采样单位: **complex_id**；次数 **{N_BOOT}**；seed **{SEED}**")
    emit(f"- 若某次迭代 unique complex < **{MIN_CX}** → 跳过并计入 `n_skipped`")
    emit(f"- 指标实现: 复用 V5 `v5_common`（`per_structure_spearman` / `topk_metrics`）")
    emit("")

    if not STD.exists():
        emit("!! `standardized_preds.csv` 不存在 → 无法执行")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    d = pd.read_csv(STD, low_memory=False)
    ca = d[d["dataset_name"] == "CATH"].copy()
    n_cx = ca["complex_id"].nunique()
    emit(f"- CATH 数据: **{len(ca)} 行 / {n_cx} complex**；"
         f"模型 **{sorted(ca['model'].unique())}**")
    emit("")
    if n_cx < 30:
        emit(f"> ⚠️ **样本量警示**：仅 {n_cx} 个 complex（< 30）。"
             f"CI 会较宽，且 macro 统计量的抽样分布可能非正态 → "
             f"**结论需谨慎**（V6.1 全局约束 9）。")
        emit("")

    rows = []
    emit("## 结果")
    emit("")
    emit("| 模型 | n 行 | n complex | 指标 | point | 95% CI | SD | "
         "有效次数 | 跳过 | P(>0) |")
    emit("|---|---:|---:|---|---:|---|---:|---:|---:|---:|")
    for mo in MODELS:
        g0 = ca[ca["model"] == mo]
        if len(g0) < 30:
            emit(f"| {mo} | {len(g0)} | — | — | — | — | — | — | — | — |")
            continue
        y = g0["true_ddG"].to_numpy(float)
        p = g0["pred_ddG"].to_numpy(float)
        g = g0["complex_id"].astype(str).to_numpy()
        for mname, fn in (("macro_pdb_spearman", macro_spearman),
                          ("efficiency_at_10", eff10)):
            r = boot(fn, y, p, g)
            emit(f"| {mo} | {len(g0)} | {g0['complex_id'].nunique()} | {mname} | "
                 f"{r['point']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                 f"{r['sd']:.4f} | {r['n_valid']} | {r['n_skipped']} | "
                 f"{r['p_gt_zero']:.4f} |")
            rows.append({"dataset": "CATH", "model": mo, "metric": mname,
                         "n_rows": int(len(g0)),
                         "n_complex": int(g0["complex_id"].nunique()),
                         "point": round(r["point"], 4),
                         "ci_low": round(r["ci_low"], 4),
                         "ci_high": round(r["ci_high"], 4),
                         "sd": round(r["sd"], 4),
                         "n_valid": r["n_valid"], "n_skipped": r["n_skipped"],
                         "p_gt_zero": round(r["p_gt_zero"], 4),
                         "n_boot_requested": N_BOOT, "seed": SEED,
                         "min_complex_per_iter": MIN_CX})
    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")

    emit("")
    emit("## 判定（CI 是否显著 > 0）")
    emit("")
    emit("| 模型 | macro Spearman point | CI 下界 > 0? | Efficiency@10 point | CI 下界 > 0? |")
    emit("|---|---:|---|---:|---|")
    for mo in MODELS:
        s = res[(res["model"] == mo) & (res["metric"] == "macro_pdb_spearman")]
        e = res[(res["model"] == mo) & (res["metric"] == "efficiency_at_10")]
        if not len(s):
            continue
        s, e = s.iloc[0], e.iloc[0]
        emit(f"| {mo} | {s['point']:+.4f} | "
             f"{'✅ 是' if s['ci_low'] > 0 else '❌ 否（含 0）'} | "
             f"{e['point']:+.4f} | "
             f"{'✅ 是' if e['ci_low'] > 0 else '❌ 否（含 0）'} |")
    emit("")

    # 与 V5 baseline 并列（V5 的 CI 来自 V6 之前的 V5 主推断）
    emit("## 参考：V5 baseline 的对应 CI（来自 V5 主推断，抽样单位 PDB，n=61）")
    emit("")
    emit("| 指标 | point | 95% CI | n_clusters |")
    emit("|---|---:|---|---:|")
    emit("| macro per-PDB Spearman (seed42) | +0.1048 | [−0.0219, +0.2385] | 61 |")
    emit("| macro per-PDB Spearman (seed43) | +0.0905 | [−0.0379, +0.2180] | 61 |")
    emit("| macro per-PDB Spearman (seed44) | +0.0914 | [−0.0335, +0.2171] | 61 |")
    emit("")
    emit("> 该三行的原始产物：`experiments/v5_0/p0_benchmark/baseline_repro_report.md`。")
    emit("> **注意**：V5 与 CATH 数据集不同，此表**仅用于量级对照**，不是同集比较。")
    emit("")

    emit("## RBM 趋势性支持（V6.1 PHASE 1.4）")
    emit("")
    rb = d[d["dataset_name"] == "DSSA-PPI_RBM"]
    if len(rb):
        emit("| 模型 | n 行 | n 分箱 | macro Spearman（点估计，**无 CI**） |")
        emit("|---|---:|---:|---:|")
        for mo in ["FoldX", "DSSA-PPI", "SAAMBE-SEQ", "PPIFormer", "MpbPPI",
                   "ESM-2-FFN", "ProBASS", "MuTon"]:
            gg = rb[rb["model"] == mo]
            if len(gg) < 30:
                continue
            v = macro_spearman(gg["true_ddG"].to_numpy(float),
                               gg["pred_ddG"].to_numpy(float),
                               gg["complex_id"].astype(str).to_numpy())
            emit(f"| {mo} | {len(gg)} | {gg['complex_id'].nunique()} | {v:+.4f} |")
        emit("")
        emit("> ⚠️ **仅作趋势性支持**：RBM 只有 **1 个分箱**（人为按 `pos` 连续分块），")
        emit("> **无法计算 CI**，且其真值 `EXP_fitness` 是**适应度而非 ΔΔG**。")
        emit("> 该数据集的任何数值**不得**与 CATH / V5 的数字直接比较。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
