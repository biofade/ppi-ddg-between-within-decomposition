# -*- coding: utf-8 -*-
"""scripts/v6_1/05_pairwise_comparison.py -- V6.1 PHASE 3: CATH 内 pairwise 对比.

对 21 个 complex 做 **pairwise cluster bootstrap**：每对模型重采样 complex
（5000 次），计算 `Macro per-PDB Spearman` 之差的 95% CI 与 p 值。

输出
----
* `experiments/v6_1/phase3_pairwise/pairwise_results.csv`
* `experiments/v6_1/phase3_pairwise/pairwise_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH3 = C.ROOT / "experiments" / "v6_1" / "phase3_pairwise"
PH3.mkdir(parents=True, exist_ok=True)
OUT_CSV = PH3 / "pairwise_results.csv"
OUT_MD = PH3 / "pairwise_report.md"
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"

N_BOOT = 5000
SEED = 42
MIN_PS_N = 5
PAIRS = [("Flex-ddG", "ESM-IF"), ("Flex-ddG", "proteinMPNN"),
         ("Flex-ddG", "FoldX"), ("FoldX", "ESM-IF"), ("FoldX", "proteinMPNN"),
         ("ESM-IF", "proteinMPNN"), ("DDGPred", "ESM-IF"),
         ("DDGPred", "proteinMPNN"), ("DDGPred", "Flex-ddG"),
         ("DDGPred", "FoldX")]
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def per_cx_spearman(df: pd.DataFrame, model: str) -> dict[str, float]:
    out = {}
    for cx, gg in df[df["model"] == model].groupby("complex_id"):
        if len(gg) < MIN_PS_N:
            continue
        v = C.spearman(gg["true_ddG"].to_numpy(float),
                       gg["pred_ddG"].to_numpy(float))
        if v == v:
            out[str(cx)] = float(v)
    return out


def main() -> int:
    t0 = time.time()
    emit("# V6.1 PHASE 3 — CATH 内模型 pairwise 对比")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 单位: complex（per-complex Spearman，每 complex ≥{MIN_PS_N} 条突变）")
    emit(f"- 重采样: **cluster bootstrap，{N_BOOT} 次**，seed {SEED}")
    emit(f"- 指标: **Macro per-PDB Spearman** 的配对差")
    emit("")

    if not STD.exists():
        emit("!! `standardized_preds.csv` 不存在")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    d = pd.read_csv(STD, low_memory=False)
    ca = d[d["dataset_name"] == "CATH"].copy()
    models = sorted(ca["model"].unique())
    per = {m: per_cx_spearman(ca, m) for m in models}
    emit("| 模型 | 可用 complex | per-complex Spearman 均值 |")
    emit("|---|---:|---:|")
    for m in models:
        v = np.array(list(per[m].values()), float)
        emit(f"| {m} | {len(v)} | {v.mean():+.4f} |")
    emit("")

    rs = np.random.RandomState(SEED)
    rows = []
    emit("| A | B | n 配对 | Δ(A−B) point | 95% CI | p（双侧） | 显著? |")
    emit("|---|---|---:|---:|---|---:|---|")
    for a, b in PAIRS:
        if a not in per or b not in per:
            continue
        common = sorted(set(per[a]) & set(per[b]))
        if len(common) < 5:
            emit(f"| {a} | {b} | {len(common)} | — | — | — | 配对数不足 |")
            continue
        va = np.array([per[a][c] for c in common], float)
        vb = np.array([per[b][c] for c in common], float)
        dlt = va - vb
        obs = float(dlt.mean())
        n = len(dlt)
        boot = np.array([dlt[rs.randint(0, n, n)].mean() for _ in range(N_BOOT)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        dc = boot - boot.mean()
        pl = (np.sum(dc <= -abs(obs)) + 1) / (N_BOOT + 1)
        pr = (np.sum(dc >= abs(obs)) + 1) / (N_BOOT + 1)
        p = float(min(1.0, 2 * min(pl, pr)))
        sig = (lo > 0) or (hi < 0)
        emit(f"| {a} | {b} | {n} | {obs:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {p:.4f} | "
             f"{'✅ 显著' if sig else '❌ 不显著'} |")
        rows.append({"model_a": a, "model_b": b, "n_paired_complex": n,
                     "delta_mean": round(obs, 4),
                     "ci_low": round(float(lo), 4), "ci_high": round(float(hi), 4),
                     "p_two_sided": round(p, 4), "significant": bool(sig),
                     "bootstrap_unit": "complex", "n_boot": N_BOOT, "seed": SEED})
    res = pd.DataFrame(rows)
    if len(res):
        res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("")
    if len(res):
        nsig = int(res["significant"].sum())
        emit(f"## 汇总：**{nsig} / {len(res)} 对差异显著**（95% CI 排除 0）")
        emit("")
        emit("显著对（按 |Δ| 降序）：")
        emit("")
        emit("| A | B | Δ | 95% CI | p |")
        emit("|---|---|---:|---|---:|")
        for _, r in res[res["significant"]].sort_values(
                "delta_mean", key=lambda s: -s.abs()).iterrows():
            emit(f"| {r['model_a']} | {r['model_b']} | {r['delta_mean']:+.4f} | "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {r['p_two_sided']:.4f} |")
        emit("")
        # 类别内 / 类别间
        phys = {"Flex-ddG", "FoldX"}
        seq = {"ESM-IF", "proteinMPNN"}
        cross = res[((res["model_a"].isin(phys)) & (res["model_b"].isin(seq)))
                    | ((res["model_a"].isin(seq)) & (res["model_b"].isin(phys)))]
        within = res.drop(cross.index)
        emit("### 按对比类型拆分")
        emit("")
        emit("| 类型 | 对数 | 显著对数 | 所有 Δ 方向 |")
        emit("|---|---:|---:|---|")
        emit(f"| 跨类（物理 vs 序列） | {len(cross)} | "
             f"{int(cross['significant'].sum())} | "
             f"{'全部 A>B 为物理类占优' if len(cross) and (cross['delta_mean'] > 0).all() else '见 CSV'} |")
        emit(f"| 类内（同族互比） | {len(within)} | "
             f"{int(within['significant'].sum())} | — |")
        emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
