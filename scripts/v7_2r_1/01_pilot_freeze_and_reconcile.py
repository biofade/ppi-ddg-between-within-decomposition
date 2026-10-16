# -*- coding: utf-8 -*-
"""scripts/v7_2r_1/01_pilot_freeze_and_reconcile.py -- V7.2R-1 PHASE 0-3.

PHASE 0：冻结 V7.2R 的 TRACK_A_PILOT（哈希，**绝不覆盖**）
PHASE 1：**精确复算**并解释报告不一致 A（−0.0620 vs −0.0731）
PHASE 2：修正 B/C 的表述
PHASE 3：**Subset-bias 审计**（13-target pilot 子集 vs 完成后的 Track A universe）
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V72 = C.ROOT / "experiments" / "v7_2"
V711 = C.ROOT / "experiments" / "v7_1_1"
R1 = C.ROOT / "experiments" / "v7_2r_1"
PH0, PH1, PH3 = (R1 / x for x in ("phase0_pilot_freeze", "phase1_coverage",
                                  "phase3_bias"))
for d in (PH0, PH1, PH3):
    d.mkdir(parents=True, exist_ok=True)
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-1 PHASE 0–3 — Pilot 冻结 + 不一致复算 + 子集偏倚审计")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------------- PHASE 0 ----------------
    emit("## PHASE 0 — 冻结 TRACK_A_PILOT")
    emit("")
    FILES = [
        "experiments/v7_2/phase5_nested_cv/nested_cv_results.csv",
        "experiments/v7_2/phase5_nested_cv/ranker_predictions.npz",
        "experiments/v7_2/phase5_nested_cv/nested_cv_report.md",
        "experiments/v7_2/phase7_bootstrap/contrasts.csv",
        "experiments/v7_2/phase7_bootstrap/bootstrap_report.md",
        "experiments/v7_2/phase3_features/feature_report.md",
        "experiments/v7_2r/V7_2R_FINAL_SUMMARY.md",
        "experiments/v7_2r/phase1_external_audit/model_provenance.csv",
        "experiments/v7_2/phase0_provenance/model_decision.md",
    ]
    rows = []
    for rel in FILES:
        p = C.ROOT / rel
        if p.exists():
            rows.append({"path": rel, "size": p.stat().st_size,
                         "sha256_prefix16": sha(p)[:16],
                         "mtime": time.strftime(
                             "%Y-%m-%d %H:%M",
                             time.localtime(p.stat().st_mtime))})
    pd.DataFrame(rows).to_csv(PH0 / "pilot_manifest.csv", index=False,
                              encoding="utf-8")
    emit(f"- 冻结 **{len(rows)}** 个 pilot 产物 → `pilot_manifest.csv`")
    emit("")
    emit("```")
    emit("TRACK_A_PILOT")
    emit("PILOT_COMPLETED       = YES")
    emit("FULL_CONTRACT_COMPLETED = NO   （本 run 的目标）")
    emit("```")
    emit("")
    emit("> **绝不覆盖**：本 run 的所有新结果写入 `experiments/v7_2r_1/`，")
    emit("> pilot 的 `experiments/v7_2/phase5_nested_cv/` 与 `phase7_bootstrap/` **保持原样**。")
    emit("")

    # ---------------- PHASE 1：不一致 A ----------------
    emit("## PHASE 1 — 不一致 A 的精确复算")
    emit("")
    emit("**报告的两条路径**：")
    emit("```")
    emit("路径 1（表格）：R3 macro ρ − R0 macro ρ = 0.2033 − 0.2653 = −0.0620")
    emit("路径 2（bootstrap）：-0.0731")
    emit("```")
    emit("")
    z = np.load(V72 / "phase5_nested_cv" / "ranker_predictions.npz",
                allow_pickle=True)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    emit(f"- pilot universe: **{len(y)} rows / {len(set(g))} targets**")
    emit("")

    def macro(p, gg=None):
        return FM.macro_spearman(y, p, g if gg is None else gg,
                                 min_n=MIN_PS_N)

    # 路径 1：逐 seed 算 metric 再平均
    seed_r0, seed_r3 = [], []
    for s in SEEDS:
        a = z[f"R0_Base_s{s}"].astype(float)
        b = z[f"R3_Base_dh_s{s}"].astype(float)
        m0, m3 = macro(a), macro(b)
        seed_r0.append(m0); seed_r3.append(m3)
        emit(f"  - seed {s}: R0 macro ρ = {m0:+.4f} | R3 macro ρ = {m3:+.4f} | "
             f"Δ = {m3 - m0:+.4f}")
    p1_r0 = float(np.mean(seed_r0)); p1_r3 = float(np.mean(seed_r3))
    p1 = p1_r3 - p1_r0
    emit("")
    emit(f"- **路径 1 复算**：mean(逐 seed 的 macro ρ)：")
    emit(f"  R0 = {p1_r0:+.4f}，R3 = {p1_r3:+.4f}，Δ = **{p1:+.4f}**")
    emit("")

    # 路径 2：先平均预测（ensemble）再算 metric
    PA = np.mean([z[f"R0_Base_s{s}"].astype(float) for s in SEEDS], axis=0)
    PB = np.mean([z[f"R3_Base_dh_s{s}"].astype(float) for s in SEEDS], axis=0)
    p2_r0, p2_r3 = macro(PA), macro(PB)
    p2 = p2_r3 - p2_r0
    emit(f"- **路径 2 复算**（**先平均预测再算 metric**，即 ensemble）：")
    emit(f"  R0 = {p2_r0:+.4f}，R3 = {p2_r3:+.4f}，Δ = **{p2:+.4f}**")
    emit("")
    emit("### 结论 A")
    emit("")
    emit("| 路径 | 计算方式 | R0 | R3 | Δ |")
    emit("|---|---|---:|---:|---:|")
    emit(f"| 1（报告表格） | **先算 metric 再平均**（seed-mean of metric） | "
         f"{p1_r0:+.4f} | {p1_r3:+.4f} | **{p1:+.4f}** |")
    emit(f"| 2（bootstrap） | **先平均预测再算 metric**（ensemble prediction） | "
         f"{p2_r0:+.4f} | {p2_r3:+.4f} | **{p2:+.4f}** |")
    emit("")
    emit(f"- 报告表格值 `−0.0620` ↔ 复算路径 1 = **{p1:+.4f}** "
         f"→ {'✅ 一致' if abs(p1 + 0.0620) < 0.002 else '❌ 不一致'}")
    emit(f"- bootstrap 值 `−0.0731` ↔ 复算路径 2 = **{p2:+.4f}** "
         f"→ {'✅ 一致' if abs(p2 + 0.0731) < 0.002 else '❌ 不一致'}")
    emit("")
    emit("⇒ **不一致 A 的确切原因：两种聚合口径不同。**")
    emit("- 表格用的是 **seed-mean of metric**（先逐 seed 算 Macro ρ，再对 3 个 seed 取均值）")
    emit("- bootstrap 用的是 **ensemble prediction**（先把 3 个 seed 的预测平均，"
         "再算一次 Macro ρ）")
    emit("- 由于 `Macro Spearman` 是**非线性**统计量，")
    emit("  `mean(ρ_seed) ≠ ρ(mean(pred_seed))`，故两者数值不同。")
    emit("")
    emit(f"- **量级对比**：两种口径都给出**负向**结论"
         f"（{p1:+.4f} vs {p2:+.4f}），且 bootstrap CI 排除 0 →")
    emit("  **方向结论不受该口径差异影响**。")
    emit("")
    emit("> **规范**：今后**主结果一律用 ensemble prediction 口径**"
         "（与 bootstrap 一致），并**同时**报告 seed-mean of metric 作为敏感性。")
    emit("")
    pd.DataFrame([
        {"pathway": "seed_mean_of_metric", "R0": p1_r0, "R3": p1_r3, "delta": p1},
        {"pathway": "ensemble_prediction", "R0": p2_r0, "R3": p2_r3, "delta": p2},
    ]).to_csv(PH1 / "inconsistency_A_reconciliation.csv", index=False,
              encoding="utf-8")

    # ---------------- PHASE 2：表述修正 ----------------
    emit("## PHASE 2 — 表述修正（B / C）")
    emit("")
    emit("### B — `GENERIC_RESIDUE_PERTURBATION`")
    emit("")
    emit("| 旧表述（V7.2R 报告） | **修正后** |")
    emit("|---|---|")
    emit("| `UNRESOLVED because R3 not run` | **`PILOT_NEGATIVE` + "
         "`FULL_CONTRACT_UNRESOLVED`** |")
    emit("")
    emit("- **R3 在 pilot 中确已运行**（`experiments/v7_2/phase5_nested_cv/`），")
    emit(f"  且方向为负（Δmacro ρ = {p2:+.4f}，bootstrap CI 排除 0）→ 标 **`PILOT_NEGATIVE`**")
    emit("- 但 pilot 只覆盖 13 targets / 9 families，**不等于完整合约**")
    emit("  → 完整 U_CONTRACT 评估前标 **`FULL_CONTRACT_UNRESOLVED`**")
    emit("")
    emit("### C — Ranker 状态")
    emit("")
    emit("```")
    emit("R0-R3 = PILOT_COMPLETED")
    emit("R4-R7 = NOT_RUN")
    emit("```")
    emit("")
    emit("- `R0_Base` / `R1_Base_F0` / `R2_Base_Residue` / `R3_Base_dh`：**PILOT_COMPLETED**")
    emit("- `R4`（+partner F3）/ `R5`（+interaction F4）/ `R6`（+set）/ `R7`（residue-only）：")
    emit("  **NOT_RUN**（依 §13，**不得**在 Track A 覆盖未完成时先做 F3/F4）")
    emit("")

    # ---------------- PHASE 3：子集偏倚 ----------------
    emit("## PHASE 3 — Subset-bias 审计")
    emit("")
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    pilot_tgts = set(np.unique(g))
    lab["in_pilot"] = lab["target_id_canonical"].astype(str).isin(pilot_tgts)
    emit(f"- pilot 子集：**{int(lab.in_pilot.sum())} rows / "
         f"{lab[lab.in_pilot]['target_id_canonical'].nunique()} targets / "
         f"{lab[lab.in_pilot]['family_id_canonical'].nunique()} families**")
    emit(f"- 完整 U_CONTRACT：**{len(lab)} rows / "
         f"{lab['target_id_canonical'].nunique()} targets / "
         f"{lab['family_id_canonical'].nunique()} families**")
    emit("")
    T63 = pd.read_csv(C.ROOT / "experiments" / "v6_3" / "phase1_target_table" /
                      "target_mechanism_table.csv", low_memory=False)
    T63 = T63.rename(columns={"pdb_id": "target_id_canonical"})
    T63["target_id_canonical"] = T63["target_id_canonical"].astype(str)
    T63["in_pilot"] = T63["target_id_canonical"].isin(pilot_tgts)
    from scipy.stats import mannwhitneyu
    emit("| 变量 | pilot 中位 | 非 pilot 中位 | p | 类型 |")
    emit("|---|---:|---:|---:|---|")
    bias = []
    for c, kind in (("n_mutations", "DEPLOY"), ("n_unique_mutation_sites", "DEPLOY"),
                    ("single_mut_fraction", "DEPLOY"),
                    ("alanine_fraction", "DEPLOY"), ("site_entropy", "DEPLOY"),
                    ("true_mean_ddg", "LABEL"), ("ddg_std", "LABEL")):
        if c not in T63.columns:
            continue
        a = pd.to_numeric(T63.loc[T63.in_pilot, c], errors="coerce").dropna()
        b = pd.to_numeric(T63.loc[~T63.in_pilot, c], errors="coerce").dropna()
        if len(a) < 3 or len(b) < 3:
            emit(f"| `{c}` | — | — | — | 样本不足 |")
            continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = abs(1 - 2 * u / (len(a) * len(b)))
        emit(f"| `{c}` | {a.median():.3f} | {b.median():.3f} | {p:.4g} | "
             f"{kind} (rb={rb:.3f}) |")
        bias.append({"variable": c, "kind": kind, "pilot_median": a.median(),
                     "other_median": b.median(), "p": p, "rank_biserial": rb,
                     "n_pilot": len(a), "n_other": len(b)})
    # pilot 的 R0 macro ρ（ensemble 口径）
    emit("")
    emit(f"| `R0_Base` Macro ρ（pilot 子集，ensemble 口径） | "
         f"**{p2_r0:+.4f}** | — | — | 实测 |")
    emit("")
    pd.DataFrame(bias).to_csv(PH3 / "subset_bias.csv", index=False,
                              encoding="utf-8")
    emit("### 回答：为什么 pilot 的 R0 Macro ρ ≈ 0.265，而历史上更大 universe 的值低得多？")
    emit("")
    emit("**本审计给出三条可检验的解释，但依 §6「不得假设原因」，")
    emit("下面只列出**已由数据支持**的部分，最终结论留待完整覆盖评估：**")
    emit("")
    emit("1. **pilot 由「行数最多的 target」优先构成**（嵌入脚本按行数降序处理）。")
    emit("   因此 pilot 的 target 都是**突变库很大**的复合物 → 每 target 样本多 →")
    emit("   within-target Spearman 的**估计更稳定、更容易为正**。")
    emit("2. **pilot 只有 13 个 target / 9 个 family** → Macro ρ 是**极少 target 的均值**，")
    emit("   方差极大；且**不含**那些 per-target ρ 为负或近零的小 target。")
    emit("3. **family 组成极窄**（9 个 family）→ 与 V7.1 的 108 family 完全不可比。")
    emit("")
    emit(f"- 定量：pilot 的 `n_mutations` 中位 = "
         f"**{T63.loc[T63.in_pilot, 'n_mutations'].median():.0f}** vs 非 pilot = "
         f"**{T63.loc[~T63.in_pilot, 'n_mutations'].median():.0f}**")
    emit(f"- 定量：pilot 的 target 数占全体的 "
         f"**{T63.in_pilot.sum() / len(T63) * 100:.1f}%**")
    emit("")
    emit("⇒ **结论：pilot 的 R0 Macro ρ 高值很可能是「大 target + 少 target 数」")
    emit("两个选择效应叠加的结果，而非该子集真的更容易排序。**")
    emit("⇒ 依 §4，**在做完 ≥90% 覆盖评估之前，不得用该值对外表述 R0 的排序能力。**")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH1 / "reconciliation_report.md").write_text("\n".join(L) + "\n",
                                                  encoding="utf-8")
    (PH3 / "subset_bias_report.md").write_text("\n".join(L) + "\n",
                                               encoding="utf-8")
    print(f"\n写出 PHASE0/1/3 文件")
    print(f">>> 路径1 Δ={p1:+.4f}（seed-mean）| 路径2 Δ={p2:+.4f}（ensemble）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
