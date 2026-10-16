# -*- coding: utf-8 -*-
"""scripts/v8/14_remedy_case_study.py -- V8 PHASE 14 / §45/§46/§47：benchmark-guided remedy 案例。

**不重跑任何模型**（§46 硬规则）。只用**已冻结**的 V7.2 / V7.3 预测与 bootstrap 输出。

## 关键设计：先分层，再比较（§21 匹配宇宙）

不同 V7.x 阶段的预测是在**不同大小的宇宙**上生成的
（V7.2R1 有两种：174-target 与 166-target 子集；V7.3 为 166-target）。
**跨宇宙比较原始数值是不允许的**。因此本脚本：

  * 先按 `n_targets_ge5` 分层，只在**同层内**做 remedy 序列比较
  * 报告每层的 universe 大小与模型数
  * 明确标注哪些比较是**跨宇宙（不可直接比较）**

## 核心问题（§45）

Capability-resolved 评价能否：
  1. 识别出 V7.2 表示放大的**具体缺陷**（校准涨、排序不涨）？
  2. 区分**失败的 remedy**（8M / 650M / PPI context / static）与
     **成功的 remedy**（mutation × partner physics）？
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
REM = OUT / "remedy"
REM.mkdir(parents=True, exist_ok=True)

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
u = S[S.universe == "U_CANONICAL"].copy()
print(f"U_CANONICAL 模型 = {len(u)}")

# ---------- 1. 分层 ----------
print("\n=== 按 n_targets_ge5 分层（跨层不可直接比较）===")
u["stratum"] = u.n_targets_ge5.astype(int)
st = u.groupby("stratum").agg(models=("model_name", "size"),
                              families=("model_family", "nunique"))
print(st.to_string())
rows = []
for s_, g in u.groupby("stratum"):
    for fam, gg in g.groupby("model_family"):
        rows.append({"stratum_nge5": s_, "model_family": fam, "n": len(gg),
                     "overall_r": gg.overall_pearson.median(),
                     "macro": gg.macro_spearman.median(),
                     "C_W": gg.C_W.median(), "C_B": gg.C_B_row.median(),
                     "rmse": gg.rmse.median(), "eff10": gg["efficiency@10"].median()})
F = pd.DataFrame(rows).sort_values(["stratum_nge5", "macro"], ascending=[True, False])
F.to_csv(REM / "family_by_stratum.csv", index=False, encoding="utf-8")
print("\n=== 每层内各模型族三轴中位 ===")
for s_, g in F.groupby("stratum_nge5"):
    print(f"\n[stratum n_ge5={s_}]")
    print(g.drop(columns="stratum_nge5").round(4).to_string(index=False))

# ---------- 2. remedy 序列（同层内）----------
print("\n" + "=" * 70)
print("REMEDY SEQUENCE（同层内比较）")
print("=" * 70)

SEQ = [
    # (step, label, family, 需要匹配的 key 片段, 结局说明)
    ("R1", "8M PLM Δh（V7.2R1 R3）", "V7.2R1_representation", "R3_Base_dh"),
    ("R2", "650M 链隔离 Δr（V7.2R2 B2）", "V7.2R2_strong_repr", "B2_strong_dr"),
    ("R3", "PPI 条件化（V7.2R2 B3a/B4）", "V7.2R2_strong_repr", "B3a_static_PPI"),
    ("R4", "PPI 扰动（V7.2R2 B4）", "V7.2R2_strong_repr", "B4_PPI_perturb"),
    ("R5", "V7.3 F1 mutation-only", "V7.3_physics", "F1_Base_Mut"),
    ("R5", "V7.3 F2 static local", "V7.3_physics", "F2_Base_Static"),
    ("R5", "V7.3 F3 partner static", "V7.3_physics", "F3_Base_EnvPartner"),
    ("R5", "V7.3 F4 mutation×partner 扰动", "V7.3_physics", "F4_Base_PartnerPert"),
    ("R5", "V7.3 F5 intra 扰动", "V7.3_physics", "F5_Base_IntraPert"),
    ("R5", "V7.3 F8 full", "V7.3_physics", "F8_Full"),
]
rec = []
for step, label, fam, frag in SEQ:
    g = u[(u.model_family == fam) & (u.model_name.str.contains(frag, regex=False))]
    if not len(g):
        continue
    # 取 seed 42（若无则取全部中位）
    s42 = g[g.model_name.str.endswith("_s42")]
    gg = s42 if len(s42) else g
    rec.append({"step": step, "label": label, "family": fam,
                "stratum_nge5": int(gg.n_targets_ge5.median()),
                "n": len(gg),
                "overall_r": gg.overall_pearson.median(),
                "rmse": gg.rmse.median(),
                "macro_spearman": gg.macro_spearman.median(),
                "C_W": gg.C_W.median(),
                "C_B_row": gg.C_B_row.median(),
                "eff10": gg["efficiency@10"].median(),
                "ndcg10": gg["ndcg@10"].median()})
RS = pd.DataFrame(rec)
RS.to_csv(REM / "remedy_sequence.csv", index=False, encoding="utf-8")
print(f"{'step':<4}{'label':<34}{'n_ge5':>6}{'Overall r':>11}{'Macro ρ':>10}"
      f"{'C_W':>9}{'C_B':>9}{'Eff@10':>9}")
for r in RS.itertuples():
    print(f"{r.step:<4}{r.label:<34}{r.stratum_nge5:>6}{r.overall_r:>11.4f}"
          f"{r.macro_spearman:>10.4f}{r.C_W:>9.4f}{r.C_B_row:>9.4f}{r.eff10:>9.4f}")

# ---------- 3. 分解：全局增益 vs 排序增益 ----------
print("\n=== 关键分解：每个 remedy 步的增益落在哪一轴？===")
b = RS[RS.label.str.startswith("V7.3 F1")].iloc[0] if len(
    RS[RS.label.str.startswith("V7.3 F1")]) else None
if b is not None:
    print(f"基线 = V7.3 F1（mutation-only）: Overall r={b.overall_r:.4f}, "
          f"Macro ρ={b.macro_spearman:.4f}, C_W={b.C_W:.4f}")
    print(f"\n{'step':<34}{'ΔOverall r':>12}{'ΔMacro ρ':>11}{'ΔC_W':>10}{'ΔEff@10':>10}")
    dd = []
    for r in RS[RS.step == "R5"].itertuples():
        if r.label.startswith("V7.3 F1"):
            continue
        dd.append({"label": r.label,
                   "d_overall_r": r.overall_r - b.overall_r,
                   "d_macro": r.macro_spearman - b.macro_spearman,
                   "d_C_W": r.C_W - b.C_W,
                   "d_eff10": r.eff10 - b.eff10})
        print(f"{r.label:<34}{r.overall_r-b.overall_r:>+12.4f}"
              f"{r.macro_spearman-b.macro_spearman:>+11.4f}"
              f"{r.C_W-b.C_W:>+10.4f}{r.eff10-b.eff10:>+10.4f}")
    pd.DataFrame(dd).to_csv(REM / "v73_axis_gains.csv", index=False,
                            encoding="utf-8")

# ---------- 4. 判定 ----------
print("\n" + "=" * 70)
print("H5 判定依据")
print("=" * 70)
f_ok = RS[RS.label.str.contains("F4 mutation")]
if len(f_ok):
    r4 = f_ok.iloc[0]
    print(f"成功 remedy（F4 mutation×partner）: ΔMacro ρ = "
          f"{r4.macro_spearman-b.macro_spearman:+.4f}（相对 F1）")
# 失败的表示放大：V7.2R2 内部
b2 = u[(u.model_family == "V7.2R2_strong_repr") &
       (u.model_name.str.contains("B2_strong_dr", regex=False)) &
       (u.model_name.str.endswith("_s42"))]
b0 = u[(u.model_family == "V7.2R2_strong_repr") &
       (u.model_name.str.contains("B0_Base", regex=False)) &
       (u.model_name.str.endswith("_s42"))]
if len(b0) and len(b2):
    print(f"失败 remedy（650M 表示）: ΔOverall r = "
          f"{b2.overall_pearson.iloc[0]-b0.overall_pearson.iloc[0]:+.4f}, "
          f"ΔMacro ρ = {b2.macro_spearman.iloc[0]-b0.macro_spearman.iloc[0]:+.4f}")

lines = ["# V8 Benchmark-Guided Remedy Case Study（§45/§46/§47）\n",
         "> **未重跑任何模型**（§46）。全部数字来自已冻结的 V7.2/V7.3 预测制品。\n",
         "> ⚠️ **分层警告（§21）**：不同 V7.x 阶段的预测生成于**不同大小的宇宙**，",
         "> 跨层数值**不可直接比较**。下表已按 `n_targets_ge5` 分层。\n",
         "## 1. 分层\n", "```", st.to_string(), "```\n",
         "## 2. Remedy 序列（同层内）\n",
         "| step | label | n_ge5 | Overall r | Macro ρ | C_W | C_B | Eff@10 |",
         "|---|---|---:|---:|---:|---:|---:|---:|"]
for r in RS.itertuples():
    lines.append(f"| {r.step} | {r.label} | {r.stratum_nge5} | {r.overall_r:+.4f} | "
                 f"{r.macro_spearman:+.4f} | {r.C_W:+.4f} | {r.C_B_row:+.4f} | "
                 f"{r.eff10:+.4f} |")
lines += ["\n## 3. 增益落在哪一轴（相对 V7.3 F1 mutation-only）\n"]
if b is not None:
    lines += ["| step | ΔOverall r | ΔMacro ρ | ΔC_W | ΔEff@10 |",
              "|---|---:|---:|---:|---:|"]
    for r in RS[RS.step == "R5"].itertuples():
        if r.label.startswith("V7.3 F1"):
            continue
        lines.append(f"| {r.label} | {r.overall_r-b.overall_r:+.4f} | "
                     f"{r.macro_spearman-b.macro_spearman:+.4f} | "
                     f"{r.C_W-b.C_W:+.4f} | {r.eff10-b.eff10:+.4f} |")
(OUT / "BENCHMARK_GUIDED_REMEDY.md").write_text("\n".join(lines) + "\n",
                                                encoding="utf-8")
print("\n写出 remedy_sequence.csv / family_by_stratum.csv / v73_axis_gains.csv / "
      "BENCHMARK_GUIDED_REMEDY.md")
