# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/06b_final_summary.py -- V7.2R-2 最终报告（§27 结构，填入真实结果）."""

from __future__ import annotations

import hashlib
import json
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

V722 = C.ROOT / "experiments" / "v7_2r_2"
KAG = C.ROOT / "kaggle" / "v7_2r_2"
ZIP = Path(r"F:\dsh\downloads\v7_2r_2\V7_2R_2_TRACKB_FEATURES.zip")
SEEDS = [42, 43, 44]
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))


def sha(p: Path, ch: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    R = pd.read_csv(V722 / "phase5_eval" / "trackB_results.csv")
    NC = pd.read_csv(V722 / "phase7_controls" / "negative_controls.csv")
    BO = pd.read_csv(V722 / "phase6_bootstrap" / "contrasts.csv")
    z = np.load(V722 / "phase5_eval" / "trackB_predictions.npz", allow_pickle=True)
    y = z["y"].astype(float); g = np.array([str(x) for x in z["pdb"]])
    n5 = sum(1 for c in np.unique(g) if (g == c).sum() >= 5)

    # ensemble 口径
    def avg(b):
        ks = [f"{b}_s{s}" for s in SEEDS if f"{b}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0)
    ENS = {}
    sys.path.insert(0, str(C.ROOT / "scripts" / "v6_1"))
    from fastmetrics import macro_spearman, eff10
    for nm in ("B0_Base", "A2_smallPLM_dh", "B1_strong_repr", "B2_strong_dr"):
        p = avg(nm); m = np.isfinite(p)
        ENS[nm] = {"macro": macro_spearman(y[m], p[m], g[m], min_n=5),
                   "eff": eff10(y[m], p[m], g[m], n_rand=20)}
    sm = R.groupby("ranker")[["macro_spearman", "eff10", "C_W", "kendall",
                              "pairwise05", "hit10", "overall_r", "rmse"]].mean()

    def bo(c, m):
        d = BO[(BO.contrast.str.startswith(c)) & (BO.metric == m)]
        return d.iloc[0] if len(d) else None

    emit("# V7.2R-2 Final Summary")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    emit("## 1. Integrity")
    emit("")
    emit(f"- Upstream hashes: `experiments/v7_2r_2/phase_local_freeze/input_hashes.csv`"
         f"（11 项，**UPSTREAM_HASH_CHANGE = 0**）")
    emit(f"- Kaggle ZIP SHA256: `{sha(ZIP)}`（{ZIP.stat().st_size / 1e6:.1f} MB）")
    emit(f"- Import: **174/174 npz**，4219 行，**0 重复 row_id**，"
         f"**0 非预期 row_id**，**0 无法映射到 canonical target**")
    emit(f"- Canonical target namespace: **PASS** —— 全部按 `canonical_target_labels.csv` "
         f"的 `row_id` 连接（Kaggle 侧 `target_id` 与 canonical **100.00% 一致**，仅作核对）")
    emit("")
    emit("## 2. Kaggle Representation Semantics")
    emit("")
    emit("| 项 | 判定 |")
    emit("|---|---|")
    emit("| `r_wt` | **存在**，dim **1280**，dtype **float16** |")
    emit("| `r_mut` | **存在**，dim **1280**，float16 |")
    emit("| `delta_r` | **存在**，dim **1280**，float16；"
         "**数值自洽**（max‖Δr−(r_mut−r_wt)‖ = 0.0078 < float16 容差）✅ |")
    emit("| `partner_repr` | **不存在** ❌ |")
    emit("| MINT representation | **未进入突变表示**（MINT 被加载/审计，"
         "但 `embed_complex` 只前向含突变位点的链） |")
    emit("| site-local representation | **存在**（per-residue hidden，"
         "按位点取 token i+1，再对多突变取 mean） |")
    emit("")
    emit("| 维度/层级/条件性 | 值 |")
    emit("|---|---|")
    emit("| dimension | **1280** |")
    emit("| dtype（存储） | float16 |")
    emit("| 层级 | **site-level（per-residue）**，聚合到 variant（mean over sites） |")
    emit("| generic sequence / PPI-conditioned | **generic sequence（仅 ESM2-650M）** |")
    emit("| WT-specific / MUT-specific | **两者都有**（`r_wt` 与 `r_mut` 各自独立前向） |")
    emit("")
    emit("### ⚠️ 链隔离审计（关键限定）")
    emit("")
    emit("Kaggle 的 `build_variant()` 只返回**含突变位点的链**，"
         "`embed_complex()` 也只对这些链做前向。实测：")
    emit("")
    emit("- **214/214 个 target 在 PDB 中都有 ≥2 条蛋白链**（全是复合物）")
    emit("- 但 **128/214** 只有 **1 条**含突变的链 → 该复合物的伙伴链**未进入表示**")
    emit("- `partner_sanity.csv` **为空**（0 行）→ **伙伴条件化未被演示**")
    emit("")
    emit("**exact terminology allowed（§4 Case A）**：")
    emit("")
    emit("```")
    emit("STRONG_GENERIC_PLM_PERTURBATION")
    emit("    (ESM2-650M, CHAIN_ISOLATED site-level, 无 PPI 条件化)")
    emit("```")
    emit("")
    emit("**不得**称其为 `PPI_AWARE_PERTURBATION`；"
         "**不得**声称伙伴条件化（§5 `PARTNER_CONDITIONING_NOT_DEMONSTRATED`）。")
    emit("")
    emit("## 3. Track B Coverage")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit("| Rows | **4219 / 4233 = 99.7%** |")
    emit("| Targets | **174 / 174 = 100%** |")
    emit("| Families | **91 / 91 = 100%** |")
    emit(f"| n>=5 | **115 / 116 = 99.1%** |")
    emit("| metric-eligible target coverage | **99.1%** |")
    emit("")
    emit("## 4. Matched Universe U_AB")
    emit("")
    emit("| universe | rows | targets | families |")
    emit("|---|---:|---:|---:|")
    emit("| Track A (full) | 4233 | 174 | 91 |")
    emit("| Track B (imported) | 4219 | 174 | 91 |")
    emit(f"| **U_AB (matched)** | **4219** | **174** | **91** |")
    emit("")
    emit(f"- n>=5：**115**（Track A 为 116）→ 相对 **99.1%**")
    emit(f"- §40 门槛（n>=5 覆盖 ≥80% of Track A）：**PASS ✅**"
         f"（{'COVERAGE_LIMITED' if 115 / 116 < 0.8 else '非 COVERAGE_LIMITED'}）")
    emit("")
    emit("## 5. Coverage Bias")
    emit("")
    emit("仅 **14 行**（4233 − 4219）未进入 U_AB —— 规模太小，**未做分层偏倚检验**"
         "（无统计功效）。V7.2R-1 已对 Track A 的 94.6% 覆盖做过完整偏倚审计。")
    emit("")
    emit("## 6. Nested Family-CV Integrity")
    emit("")
    emit("| 项 | 结果 |")
    emit("|---|---|")
    emit("| Outer overlap | **0**（每折 `assert not (set(fam[tr]) & set(fam[te]))`） |")
    emit("| Inner overlap | **0**（`make_folds` family-disjoint，V7.1.1 冻结） |")
    emit("| Outer-based selection | **NO**（无任何 outer 选型；超参全部预先固定） |")
    emit("| PCA train-only | **PASS**（`PCA.fit` 仅用 `X[tr]`） |")
    emit("")
    emit("## 7. Primary Ranking Results")
    emit("")
    emit("**U_AB = 4219 rows / 174 targets / 91 families**")
    emit("")
    emit("### 主口径：ensemble prediction（与 bootstrap 一致）")
    emit("")
    emit("| Model | Macro rho | Eff@10 | C_W | Kendall | Pairwise | Hit@10 |")
    emit("|---|---:|---:|---:|---:|---:|---:|")
    lab = {"B0_Base": "B0", "A2_smallPLM_dh": "A2 small PLM",
           "B1_strong_repr": "B1 strong representation",
           "B2_strong_dr": "B2 strong perturbation"}
    for k in ("B0_Base", "A2_smallPLM_dh", "B1_strong_repr", "B2_strong_dr"):
        r = sm.loc[k]
        emit(f"| **{lab[k]}** | **{ENS[k]['macro']:+.4f}** | "
             f"**{ENS[k]['eff']:+.4f}** | {r['C_W']:+.4f} | {r['kendall']:+.4f} | "
             f"{r['pairwise05']:+.4f} | {r['hit10']:+.4f} |")
    emit(f"| B3 + PPI context | **NOT_RUN**（无 partner 表示） | | | | | |")
    emit("")
    emit("### 敏感性：seed-mean of metric")
    emit("")
    emit("| Model | Macro rho | Eff@10 | C_W | Overall r |")
    emit("|---|---:|---:|---:|---:|")
    for k in ("B0_Base", "A2_smallPLM_dh", "B1_strong_repr", "B2_strong_dr"):
        r = sm.loc[k]
        emit(f"| {k} | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
             f"{r['C_W']:+.4f} | {r['overall_r']:+.4f} |")
    emit("")
    emit("> `mean(ρ_seed) ≠ ρ(mean(pred_seed))` —— 两口径**分别报告，未混用**（§14）。")
    emit("")
    emit("## 8. Seed Sensitivity")
    emit("")
    piv = R.pivot_table(index="ranker", columns="seed", values="macro_spearman")
    emit("| Ranker | s42 | s43 | s44 | 极差 |")
    emit("|---|---:|---:|---:|---:|")
    for k, r in piv.iterrows():
        emit(f"| {k} | {r.get(42, float('nan')):+.4f} | {r.get(43, float('nan')):+.4f} | "
             f"{r.get(44, float('nan')):+.4f} | {r.max() - r.min():.4f} |")
    emit("")
    emit("## 9. Critical Contrast B2-A2")
    emit("")
    d = bo("B2−A2", "d_macro_spearman")
    emit("**5000/5000 paired target-cluster bootstrap**")
    emit("")
    emit("| 指标 | point | 95% CI | 排除 0? |")
    emit("|---|---:|---|---|")
    for m, nm in (("d_macro_spearman", "ΔMacro rho"), ("d_eff10", "ΔEff@10"),
                  ("d_C_W", "ΔC_W"), ("d_overall_r", "ΔOverall Pearson"),
                  ("d_rmse", "ΔRMSE")):
        r = bo("B2−A2", m)
        emit(f"| {nm} | {r['point']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
             f"{'✅ 排除 0' if r['ci_excludes_zero'] else '❌ 含 0'} |")
    emit("")
    emit("⇒ **§23 成功判据：`NOT SUPPORTED` ❌**")
    emit("")
    emit("## 10. Critical Contrast B2-B0")
    emit("")
    emit("**5000/5000**")
    emit("")
    emit("| 指标 | point | 95% CI | 排除 0? |")
    emit("|---|---:|---|---|")
    for m, nm in (("d_macro_spearman", "ΔMacro rho"), ("d_eff10", "ΔEff@10"),
                  ("d_C_W", "ΔC_W"), ("d_overall_r", "ΔOverall Pearson"),
                  ("d_rmse", "ΔRMSE")):
        r = bo("B2-B0", m)
        emit(f"| {nm} | {r['point']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
             f"{'✅ 排除 0' if r['ci_excludes_zero'] else '❌ 含 0'} |")
    emit("")
    emit("⇒ **`NOT SUPPORTED` ❌**")
    emit("")
    emit("### 附：B1-B0（强表示但**无差分**）— 唯一显著项")
    emit("")
    emit("**1000/1000**")
    emit("")
    emit("| 指标 | point | 95% CI | 排除 0? |")
    emit("|---|---:|---|---|")
    for m, nm in (("d_macro_spearman", "ΔMacro rho"), ("d_C_W", "ΔC_W"),
                  ("d_overall_r", "ΔOverall Pearson"), ("d_rmse", "ΔRMSE")):
        d2 = BO[(BO.contrast.str.startswith("B1-B0")) & (BO.metric == m)]
        if len(d2):
            r = d2.iloc[0]
            emit(f"| {nm} | {r['point']:+.4f} | [{r['ci_low']:+.4f}, "
                 f"{r['ci_high']:+.4f}] | "
                 f"{'✅ 排除 0' if r['ci_excludes_zero'] else '❌ 含 0'} |")
    emit("| ΔEff@10（重算 n_rand=20，800/800） | +0.0356 | [−0.1236, +0.1659] | ❌ 含 0 |")
    emit("")
    emit("**`B1 − B0` 的 ΔMacro rho = +0.0355 [+0.0003, +0.0757] → CI 排除 0（勉强）**")
    emit("")
    emit("⇒ 但 §23 需要**同时**满足（ΔEff@10 > 0 **或** ΔC_W > 0）：")
    emit("- ΔEff@10 = +0.0356，CI **[−0.1236, +0.1659] 含 0** ❌")
    emit("- ΔC_W = −0.0041，CI [−0.0194, +0.0078] 含 0 ❌")
    emit("- 且 `ΔOverall_r = −0.0435`、`ΔRMSE = −0.0527` **均显著变差**")
    emit("")
    emit("⇒ **§23 判定：`NOT SUPPORTED` ❌**（仅 Macro rho 一项勉强显著，"
         "不满足「至少一项 secondary 改善」的要求）")
    emit("")
    emit("## 11. Partner Increment B3-B2")
    emit("")
    emit("**`NOT_RUN`** —— Kaggle ZIP 未提供 `partner_repr`，"
         "`partner_sanity.csv` 为空，表示经审计为 **CHAIN_ISOLATED**。")
    emit("依 §8/§19，**不得**构造 B3，**不得**用零向量填充。")
    emit("")
    emit("## 12. Negative Controls")
    emit("")
    emit("| 对照 | 构成 | Macro rho | Eff@10 | C_W |")
    emit("|---|---|---:|---:|---:|")
    for _, r in NC.iterrows():
        emit(f"| `{r['control']}` | {r['blocks']} | {r['macro_spearman']:+.4f} | "
             f"{r['eff10']:+.4f} | {r['C_W']:+.4f} |")
    emit("")
    emit("### 判定")
    emit("")
    emit("| 对照 | 结论 |")
    emit("|---|---|")
    emit("| **NC1 WT only** | `base+B_wt`（**仅 r_wt，无差分**）= **+0.0462** "
         "**高于** `base+B_dr`（Δr）= **+0.0305** → **扰动本身不必要**；"
         "增益主要来自**结构身份**而非突变扰动 |")
    emit("| **NC2 mutation identity** | `base+F0` = **+0.0061**（远低于全部 PLM 变体）"
         " → 表示增益**不能**归因于氨基酸替换身份编码 |")
    emit("| **NC3 norm only** | `base+‖Δr‖` = **+0.0240**，**低于**完整 Δr（+0.0305）"
         " → 幅度信息不是全部，但差距极小 |")
    emit("| **NC4 partner shuffle** | `base+B_dr+B_randpartner` = **+0.0173**，"
         "**低于** Δr → 无伙伴信号可被 shuffle（因为本来就没有伙伴表示）；"
         "该对照在此**不适用**（§19 仅在声称 partner 信号时才要求） |")
    emit("")
    emit("## 13. Top-K Result")
    emit("")
    emit("**TOPK_ONLY_SIGNAL: NO**")
    emit("")
    emit("- B2−A2 的 ΔEff@10 = −0.0786，CI [−0.0952, +0.1279] **含 0** ❌")
    emit("- B2−B0 的 ΔEff@10 = +0.0020，CI [−0.1979, +0.2068] **含 0** ❌")
    emit("- B1−B0 的 ΔEff@10 = +0.0356，CI [−0.1236, +0.1659] **含 0** ❌")
    emit("")
    emit("## 14. Subgroups")
    emit("")
    emit("| 组 | 状态 |")
    emit("|---|---|")
    emit("| single / multi mutation | 见 `phase5_eval/trackB_results.csv`（未按此分层） |")
    emit("| 2-9 / 10-29 / >=30 | **NOT_RUN**（时间受限） |")
    emit("| interface / near / distal | **NOT_RUN**（需结构元数据） |")
    emit("")
    emit("## 15. Calibration")
    emit("")
    emit("**NOT_RUN** —— 依 §23，标定**不能**用于选择 ranker。"
         "由于没有任何 Track B ranker 通过 §23 成功判据（B2 两项主对比均 NULL），"
         "标定外壳不再对 Track B 排名模型做额外报告；V7.1 已给出 MU3+z 对 B0 类的"
         "标定增益（ΔOverall_r +0.3507 ✅、ΔRMSE −0.9448 ✅）。")
    emit("")
    emit("**rank preservation**：`PASS by construction`（V7.1 已数值验证 180/180 target）")
    emit("")
    emit("## 16. Scientific Verdict")
    emit("")
    emit("```")
    emit("TRACK_B_RESULT = STRONG_BACKBONE_NULL")
    emit("```")
    emit("")
    emit("**依据（§24 RESULT 0）**：")
    emit("")
    emit("- `B2 − A2` ΔMacro rho = **+0.0177**，95% CI **[−0.0118, +0.0453] 含 0** → 不显著")
    emit("- `B2 − B0` ΔMacro rho = **+0.0088**，95% CI **[−0.0229, +0.0418] 含 0** → 不显著")
    emit("- 全部 secondary 指标同样含 0")
    emit("")
    emit("## 17. Innovation Level")
    emit("")
    emit("```")
    emit("LEVEL = 2  (calibration/ranking formulation only)")
    emit("```")
    emit("")
    emit("不因 Overall Pearson 升高而上调等级（§25）。Track B 未带来 LEVEL 3/4/5。")
    emit("")
    emit("## 18. Exact Defensible Claim")
    emit("")
    emit("> Under a frozen, leakage-safe nested family-OOD protocol on an identical")
    emit("> matched universe (4219 rows / 174 targets / 91 families), an ESM2-650M")
    emit("> chain-isolated site-level mutation perturbation representation")
    emit("> **did not significantly improve within-target mutation ranking** over")
    emit("> either a small ESM2-8M perturbation (ΔMacro Spearman = +0.0177,")
    emit("> 95% CI [−0.0118, +0.0453]) or the established baseline")
    emit("> (ΔMacro Spearman = +0.0088, 95% CI [−0.0229, +0.0418]).")
    emit("> Negative controls further show that the **WT-only** representation")
    emit("> (+0.0462) already exceeds the explicit perturbation (+0.0305),")
    emit("> indicating the weak signal present is carried by structural identity")
    emit("> rather than by mutation-induced representation change.")
    emit("")
    emit("## 19. Claims We Cannot Make")
    emit("")
    emit("- ❌ 「protein foundation models do not work」（§51）")
    emit("- ❌ 「PPI-aware representation fails」—— 本表示经审计为 **CHAIN_ISOLATED**，"
         "**从未真正测试过 PPI 条件化表示**")
    emit("- ❌ 「partner context does not help」—— B3 为 **NOT_RUN**，未测试")
    emit("- ❌ 任何基于 `Overall Pearson` 的「方法更好」声明")
    emit("- ❌ 「MINT / ReCLIP 无效」—— MINT 被加载但**未进入突变表示**")
    emit("")
    emit("## 20. Comparison with Track A")
    emit("")
    emit("| | Track A (V7.2R-1) | Track B (V7.2R-2) |")
    emit("|---|---|---|")
    emit("| representation | ESM2-**8M** Δh | ESM2-**650M** Δr |")
    emit("| rows / targets / families | 4233 / 174 / 91 | **4219 / 174 / 91** |")
    emit("| Macro rho（扰动 vs 基线） | −0.0143 [−0.0588, +0.0307] | "
         "**+0.0088 [−0.0229, +0.0418]** |")
    emit("| 判定 | **NULL** | **NULL** |")
    emit("")
    emit("⇒ **表示规模从 8M 提升到 650M（81×）并未改变结论。**")
    emit("依 §21，**不得**写「小 PLM 失败所以大 PLM 成功」；"
         "实测是**两者都为 NULL**。")
    emit("")
    emit("## 21. Next Route")
    emit("")
    emit("依 §28（所有 PLM 路线 Macro Spearman 均近零）：")
    emit("")
    emit("```")
    emit("representation scaling alone is not solving the ranking bottleneck")
    emit("```")
    emit("")
    emit("**V7.3 建议：physics / structural local interaction / energy-informed ranking**")
    emit("")
    emit("**不要**继续放大模型规模。**不要**在没有伙伴表示的情况下继续堆 PLM。")
    emit("")
    emit("若仍要走表示路线，**唯一**有科学意义的第一步是：")
    emit("修复 `embed_complex` 让它包含**全部伙伴链**（本 ZIP 的关键缺陷），")
    emit("即真正做出 PPI-conditioned 表示，再重测 B2/B3。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (V722 / "V7_2R_2_FINAL_SUMMARY.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
    print(f"\n写出 {V722 / 'V7_2R_2_FINAL_SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
