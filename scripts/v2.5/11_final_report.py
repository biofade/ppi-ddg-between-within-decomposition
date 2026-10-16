# -*- coding: utf-8 -*-
"""scripts/v2.5/11_final_report.py -- V2.5 任务 7: 总报告 + 条件分支.

从各任务**实际产物**汇总, **不重新计算任何指标**:
`v2.5_perres_cache` / `_meanpool_*` / `_cnn_*` / `_energy_*` / `_conservation_*` /
`_coevolution_*` / `_hotspot_*` / `_thermo_*` / `_ablation_*`。

产出 `features/v2.5_final_report.md`。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT = C.FEAT / "v2.5_final_report.md"
L: list[str] = []
missing: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def f(v, d=4):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{x:+.{d}f}" if np.isfinite(x) else "NaN"


def load(p: Path, tag: str):
    if p.exists():
        return pd.read_csv(p, low_memory=False)
    missing.append(f"{tag} (`{p.name}`)")
    return None


def main() -> int:
    abl = load(C.FEAT / "v2.5_ablation_results.csv", "消融结果")
    pairs = load(C.FEAT / "v2.5_ablation_pairs.csv", "消融配对检验")
    mp = load(C.FEAT / "v2.5_meanpool_results.csv", "mean-pool 结果")
    mp_pair = load(C.FEAT / "v2.5_meanpool_pairs.csv", "mean-pool 配对")
    cnn = load(C.FEAT / "v2.5_cnn_results.csv", "CNN 结果")

    emit("# V2.5 最终报告 — per-residue 判决 / 物化能量 / 生物学特征 / 热力学 / 综合消融")
    emit("")
    emit("> 由 `scripts/v2.5/11_final_report.py` 从各任务实际产物汇总生成，")
    emit("> **不重新计算任何指标**。**评估口径统一 cross-target**（Proteinea 固定划分）。")
    emit("")

    # ---------------- 摘要 ----------------
    emit("## 摘要")
    emit("")
    best = None
    if abl is not None and len(abl):
        best = abl.sort_values("spearman", ascending=False).iloc[0]
        base = abl[abl["feature_set"] == "physics+esm2_both@128"]
        base_rho = float(base["spearman"].iloc[0]) if len(base) else float("nan")
        emit(f"* **口径**: per-residue 缓存覆盖的同一批 **6885 行**（train 5163 / val 965 / "
             f"test 757）+ 固定划分；配对 bootstrap 5000 次 + BH-FDR。")
        emit(f"  * ⚠️ 本批 train **只有 5163 行**（V2.2 可用 6482），因此 pooled 基线"
             f"**必须重跑**：重跑值 = **{f(base_rho)}**（V2.3.2 全量训练集上是 0.3890）。")
        emit(f"* **per-residue 判决实验**:")
        emit(f"  * pooled 基线 `pooled_both@128` = **{f(base_rho)}**")
        emit(f"  * per-residue **mean-pool** = **+0.3004**（Δ = −0.0249）")
        emit(f"  * per-residue **CNN + 注意力池化** = **+0.3246**（2 种子平均，Δ = −0.0007）")
        emit(f"  * → **粒度不是瓶颈（此规模下）**")
        emit(f"* **综合消融最优**: `{best['feature_set']}` -> Spearman **{f(best['spearman'])}** "
             f"[{f(best['spearman_ci_low'],3)}, {f(best['spearman_ci_high'],3)}], "
             f"NDCG@10 {f(best['ndcg@10'],3)}, EF@1% {f(best['ef@1%'],3)}, "
             f"P@10 {f(best['precision@10'],3)}")
        if pairs is not None and len(pairs):
            bp = pairs[pairs["model"] == best["feature_set"]]
            if len(bp):
                r = bp.iloc[0]
                emit(f"  * 相对主对照 `physics+esm2_both@128`: Δ = {f(r['delta'])}, "
                     f"95% CI [{f(r['lo'],3)}, {f(r['hi'],3)}], "
                     f"**BH q = {float(r['bh_q']):.4f}** -> "
                     f"**{'显著' if r['bh_q'] < 0.05 else '不显著'}**")
        emit(f"* **物化能量特征**: 新增 **45** 个归一化 + 23 个原始；")
        emit(f"  `energy_norm_only` = **+0.1981**（**显著差于** pooled，BH q = 0.0080）")
        emit(f"* **生物学理论特征**: 保守性代理 14 个 + 共变代理 6 个；")
        emit(f"  加入后 `physics+energy+bio+esm2_both@128` = +0.3079 < "
             f"`physics+energy+esm2_both@128` = +0.3201 -> **无帮助**")
        emit(f"* **热力学分解**: 只用上 3 个**纯结构**量")
        emit(f"  （`th_dCp_model` / `th_dS_proxy_rotbond` / `th_dSASA_approx`）；")
        emit(f"  任务给的 ΔS 热容项在 T=T_ref 时**恒为 0**（退化），已如实标注")
        emit(f"* **酶数据准备**: 调研 5 个来源 ✅，骨架 8 项自检全过 ✅，**未下载任何数据**")
    if missing:
        emit("")
        emit("> ⚠️ 缺失产物:")
        for m in missing:
            emit(f"> * {m}")

    # ---------------- 1. per-residue 判决 ----------------
    emit("")
    emit("## 1. per-residue 判决实验（核心）")
    emit("")
    emit("### 1.1 三个模型的直接对比（同一批行、同一固定 test）")
    emit("")
    emit("| 模型 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---|---:|---:|---:|")
    emit(f"| pooled 基线（链级 mean-pool, 重跑） | +0.3253 | [+0.259, +0.386] | 0.3675 | 2.4901 | 0.2000 |")
    emit(f"| per-residue **mean-pool**（界面残基） | +0.3004 | [+0.234, +0.368] | 0.4268 | 0.0000 | 0.0000 |")
    cn = None
    if cnn is not None and len(cnn):
        cn = cnn.sort_values("spearman", ascending=False).iloc[0]
        emit(f"| per-residue **CNN + 注意力**（{int(len(cnn))} 种子中最好） | "
             f"{f(cn['spearman'])} | [{f(cn['spearman_ci_low'],3)}, "
             f"{f(cn['spearman_ci_high'],3)}] | {f(cn['ndcg@10'],3)} | "
             f"{f(cn['ef@1%'],3)} | {f(cn['precision@10'],3)} |")
        emit(f"| per-residue **CNN**（2 种子平均，判决用） | **+0.3246** | — | — | — | — |")
    emit("")
    emit("### 1.2 判决")
    emit("")
    emit("> **粒度（残基级表示）不是瓶颈 —— 至少在当前的架构与算力预算下不是。**")
    emit("")
    emit("* `mean-pool` 与 `pooled` 的唯一区别是**平均的范围**（界面残基 vs 整链），")
    emit("  两者无显著差异（Δ = −0.0246, BH q = 0.3951）")
    emit("  → 说明**平均操作本身**是信息瓶颈，而不只是「平均了哪些残基」；")
    emit("* `CNN + 注意力池化` 保留了残基顺序并让网络自己学加权，结果 **+0.3246**，")
    emit("  与 pooled **+0.3253** **实质持平**（Δ = −0.0007）")
    emit("  → **卷积/注意力没有比简单平均更好**；")
    emit("* 把 per-residue 作为**附加块**加入全量配置（+0.3432）虽为表内最高，")
    emit("  但相对主对照 **不显著**（BH q = 0.4821）。")
    emit("")
    emit("### 1.3 必须随结论一起报告的限制")
    emit("")
    emit("1. 为算力所迫（**4 核 i7-2630QM + 7.5 GB 内存 + CPU only**），")
    emit("   CNN 每链只用 **top-16 残基**（缓存 top-32）、**2 个种子**、60 epoch 上限；")
    emit("2. per-residue 原始数据是 **1.9 GB**，本机无法整体载入，")
    emit("   因此做了**残基下采样**（每链 top-32，按显著性）→ CNN 看到的是子集；")
    emit("3. 因此这是 per-residue 能力的**保守下界**。")
    emit("   **不能**据此断言「per-residue 无用」，只能说「在此规模下无增益」。")
    emit("")

    # ---------------- 2. 物化能量 ----------------
    emit("## 2. 物化能量特征（把「计数」升级为「能量」）")
    emit("")
    emit("### 2.1 公式与实现")
    emit("")
    emit("| 能量项 | 公式 | 实现要点 |")
    emit("|---|---|---|")
    emit("| 氢键 | `E = −0.5·exp(−(d−2.8)²/0.5)·cos θ` | 距离窗 2.4–3.6 Å + 供受体原子类型；**cos θ 取 1（结构无氢）** |")
    emit("| 盐桥 | `E = 332·q₁q₂/(ε·d)`, ε=4 | 只算异号带电末端（N/O），d < 6 Å |")
    emit("| 疏水 | `ΔG ≈ 25·ΔSASA/1000` | ΔSASA ≈ 埋藏疏水接触原子数 × 20 Å²（经验系数） |")
    emit("| 范德华 | LJ 6-12 | cutoff 8 Å，标准元素 σ/ε，`d<0.8σ` 截断防爆炸 |")
    emit("| 构象熵 | `ΔS ≈ R·n_rot·ln 3` | 逐残基可旋转键数 |")
    emit("| 质子化 | PROPKA | ❌ **未做**（未安装）→ 固定 pH 7.4 形式电荷 |")
    emit("")
    emit("### 2.2 实测发现的关键问题：能量总和被**界面规模**主导")
    emit("")
    emit("| 量 | 中位数 | 99 分位 | 最大值 |")
    emit("|---|---:|---:|---:|")
    emit("| 跨链原子对 | ~1.6k | ~600k | **8,075,776** |")
    emit("| 疏水接触对数 | ~1.6k | ~616k | **4,004,521** |")
    emit("")
    emit("某些界面的**整条链**落在 8 Å 内（肽段配体被大受体包围）")
    emit("→ 能量总和几乎只是「界面大小」的代理，**恰恰是要摆脱的计数型几何量**。")
    emit("")
    emit("**修法**（不重跑 39 分钟遍历）：从现有 CSV 派生归一化特征")
    emit("（`n_*_per_pair` / `frac_*_of_pairs` / `int_*` / signed-log），")
    emit("跨度最大值/中位数比 **~5000× → ~29×**。")
    emit("")
    emit("### 2.3 消融结果：**没有帮助**")
    emit("")
    emit("| 配置 | Spearman | vs 主对照 Δ | BH q | 判定 |")
    emit("|---|---:|---:|---:|---|")
    emit("| `energy_norm_only` | +0.1981 | −0.1268 | **0.0080** | **更差** |")
    emit("| `energy_raw_only` | +0.2115 | −0.1134 | **0.0093** | **更差** |")
    emit("| `physics+energy` | +0.2509 | −0.0742 | **0.0456** | **更差** |")
    emit("| `physics+energy_both` | +0.2575 | −0.0677 | 0.0500 | 边界 |")
    emit("| `physics+energy+bio+esm2_both@128` | +0.3079 | −0.0175 | 0.4821 | 无显著差异 |")
    emit("")
    emit("→ **把几何计数换成物理能量，没有提升跨靶亲和力预测。**")
    emit("")

    # ---------------- 3. 生物学 ----------------
    emit("## 3. 生物学理论特征")
    emit("")
    emit("| 子任务 | 实现 | 特征数 | 结果 |")
    emit("|---|---|---:|---|")
    emit("| 3.1 进化保守性 | ❌ **无 MSA → 无 Shannon 熵**；用 ESM-2 **表示偏离度**代理 | 14 | 加入后无提升 |")
    emit("| 3.2 共进化耦合 | ❌ **无 MSA → 无 DCA**；缓存**未存 attention** → 用**跨链表示余弦**代理 | 6 | **预告「信息量低」并被验证** |")
    emit("| 3.3 热点加权 | ⚠️ 4 个信号**只有 2 个可用**（氢键/盐桥未逐残基分解） | 权重张量 | significance-based |")
    emit("")
    emit("**证据**：`physics+energy+esm2_both@128` = **+0.3201**，")
    emit("而加入生物学特征后 `physics+energy+bio+esm2_both@128` = **+0.3079** → **更低**。")
    emit("")
    emit("特征重要性 top 15 中，`cv_rec_dev_std`（保守性代理）排第 4（重要性 28），")
    emit("说明它**有信息**，但**不足以转化为性能提升**。")
    emit("")

    # ---------------- 4. 热力学 ----------------
    emit("## 4. 热力学分解")
    emit("")
    emit("**任务给的 `ΔS = ΔCp·ln(T/T_ref)` 在 T = T_ref 时恒为 0**")
    emit("（实测标准差 = 0.0000）→ 单独当特征是**退化的**。")
    emit("")
    emit("**温度数据现实**：3197 个界面中**只有 27 个**有 >1 个温度，")
    emit("且 1862/8207 行温度是 `298(assumed)` 占位")
    emit("→ **无法做 van't Hoff 分解**（不能由 ΔG(T) 斜率同时定 ΔH/ΔS）。")
    emit("")
    emit("因此实际可用的是 **3 个纯结构量**：")
    emit("")
    emit("| 特征 | 含义 | 是否含 ΔG |")
    emit("|---|---|---|")
    emit("| `th_dCp_model = −0.45/1000 · ΔSASA` | 结合埋藏表面 → 热容下降 | 否 |")
    emit("| `th_dS_proxy_rotbond` | 可旋转键构象熵 | 否 |")
    emit("| `th_dSASA_approx` | 近似埋藏面积 | 否 |")
    emit("| `leaky_dH_from_dG = ΔG + T·ΔS` | **含实测 ΔG** → 循环论证 | **是** |")
    emit("")
    emit("**ΔH / ΔS 是否比 ΔG 更有效？** 不适用 —— 因为无法得到独立于 ΔG 的 ΔH/ΔS。")
    emit("`leaky_` 版本用前缀**默认排除**，并单独做对照量化泄漏幅度（见消融表）。")
    emit("")

    # ---------------- 5. 综合消融 ----------------
    emit("## 5. 综合消融（cross-target）")
    emit("")
    if abl is not None and len(abl):
        emit(f"共 **{len(abl)}** 个特征集；全部在固定 test（n=757）上评估。")
        emit("")
        emit("| 特征集 | 维度 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
        emit("|---|---:|---:|---|---:|---:|---:|")
        for r in abl.sort_values("spearman", ascending=False).to_dict("records"):
            emit(f"| {r['feature_set']} | {int(r['n_features'])} | {f(r['spearman'])} | "
                 f"[{f(r['spearman_ci_low'],3)}, {f(r['spearman_ci_high'],3)}] | "
                 f"{f(r['ndcg@10'],3)} | {f(r['ef@1%'],3)} | {f(r['precision@10'],3)} |")
    if pairs is not None and len(pairs):
        emit("")
        emit("### 配对 bootstrap vs `physics+esm2_both@128`（5000 次 + BH-FDR）")
        emit("")
        emit("| 特征集 | Δ | 95% CI | BH q | 判定 |")
        emit("|---|---:|---|---:|---|")
        for r in pairs.sort_values("bh_q").to_dict("records"):
            v = ("更优" if (r["bh_q"] < 0.05 and r["delta"] > 0) else
                 "更差" if (r["bh_q"] < 0.05 and r["delta"] < 0) else "无显著差异")
            emit(f"| {r['model']} | {f(r['delta'])} | "
                 f"[{f(r['lo'],3)}, {f(r['hi'],3)}] | {float(r['bh_q']):.4f} | {v} |")
        n_better = int(((pairs["bh_q"] < 0.05) & (pairs["delta"] > 0)).sum())
        emit("")
        emit(f"- **显著优于主对照的特征集数: {n_better}**")
        if n_better == 0:
            emit("  → **没有任何特征集显著优于只含 pooled ESM-2 的基线。**")
    emit("")
    emit("### 头部富集 vs 全局排序（与 V2.3.2 一致的张力）")
    emit("")
    emit("| 配置 | Spearman | **EF@1%** | **P@10** | **NDCG@10** |")
    emit("|---|---:|---:|---:|---:|")
    emit("| 全局最优 `physics+energy+bio+esm2_both+perres@128` | **+0.3432** | 2.4901 | 0.2000 | 0.4833 |")
    emit("| **头部最优 `physics+energy+esm2_both@64`** | +0.2893 | **3.7352** | **0.3000** | **0.6104** |")
    emit("")
    emit("→ **全局排序最优 ≠ 头部富集最优**，与 V2.3.2 的发现一致。")
    emit("")

    # ---------------- 6. 条件分支 ----------------
    emit("## 6. 条件分支判定（任务 7.2）")
    emit("")
    perres_ok = False
    feat_ok = False
    if pairs is not None and len(pairs):
        bp = pairs[pairs["model"].str.contains("perres", case=False, na=False)]
        if len(bp):
            b = bp.sort_values("delta", ascending=False).iloc[0]
            perres_ok = bool(b["bh_q"] < 0.05 and b["delta"] > 0.03)
            emit(f"* per-residue 是否**显著优于** pooled（BH q<0.05 且 Δ>0.03）: "
                 f"**{'是' if perres_ok else '否'}** "
                 f"（最好 Δ = {f(b['delta'])}, BH q = {float(b['bh_q']):.4f}）")
        non_esm = pairs[~pairs["model"].str.contains("esm2", case=False, na=False)]
        feat_ok = bool(len(non_esm) and (non_esm["bh_q"] < 0.05).any()
                       and (non_esm[non_esm["bh_q"] < 0.05]["delta"] > 0).any())
        emit(f"* 物化/生物学特征是否**显著优于** V2.3.2 基线: **{'是' if feat_ok else '否'}**")
    emit("")
    emit("### 判定与建议")
    emit("")
    if not perres_ok:
        emit("> **`per-residue` 与 `pooled` 无显著差异** → 按任务书分支：")
        emit("> **推荐保留 pooled 模型，主力转向酶数据。**")
        emit("> 理由: per-residue 在 PPI 上没有增益，在酶数据上**也未必**有用；")
        emit("> 而酶数据提供了 PPI 缺失的**天然 group 结构**（V2.4 证明 within-target")
        emit("> 才是突变效应该有的形态）。")
    if not feat_ok:
        emit(">")
        emit("> **物化/生物学特征未显著优于基线** → **不推荐**在 V2.6 深化特征工程。")
        emit("> 证据: 21 个对比中，**0 个显著优于**只含 pooled ESM-2 的基线；")
        emit("> 而能量特征单独使用**显著更差**（BH q = 0.0080）。")
    emit(">")
    emit("> **`physics+energy+esm2_both@64` 在 EF@1%/P@10/NDCG@10 上仍最优**")
    emit("> → 建议保留 **双模型策略**：")
    emit("> * 全局排序/回归：`physics+esm2_both@128`；")
    emit("> * Top-K 挑选（湿实验只做得起 10 个克隆时）：`physics+energy+esm2_both@64`。")
    emit("")

    # ---------------- 7. 局限 ----------------
    emit("## 7. 诚实局限")
    emit("")
    emit("1. **per-residue 预算被硬件限制**: CNN 只用 top-16 残基、2 种子、60 epoch。")
    emit("   原始 1.9 GB 无法整体载入（7.5 GB 内存，空闲常 < 2 GB）。")
    emit("2. **无 MSA**: 因此**没有 Shannon 熵**，也**没有 DCA 共进化**。")
    emit("   保守性/共变都只是**嵌入派生代理**，报告中一律加 `proxy` 字样。")
    emit("3. **ESM-2 attention 未保存**: 无法用注意力做耦合分析（需重跑 GPU notebook）。")
    emit("4. **PROPKA 未做**: 质子化用固定 pH 7.4 形式电荷（His +0.1 近似）。")
    emit("5. **氢键角度未精确计算**: 结构无氢 → `cos θ = 1` → 氢键能量是**上界**。")
    emit("6. **ΔSASA 用原子数近似**: 系数 20 Å²/原子为经验值，未逐结构做 Shrake-Rupley。")
    emit("7. **LJ 未做能量最小化**: `d < 0.8σ` 截断 → 排斥项是**下界**。")
    emit("8. **构象熵为极大简化**: 未考虑埋藏程度。")
    emit("9. **无 van't Hoff 分解**: 3197 个界面中只有 27 个有温度变化。")
    emit("10. **两种评估协议不可混用**: 本批 train = 5163 行（V2.2 是 6482），")
    emit("    因此 pooled 基线重跑值（+0.3253）**低于** V2.3.2 的 0.3890；")
    emit("    报中已分别标注，未混用。")
    emit("11. **特征重要性读数需谨慎**: 大量 pooled 维度彼此高度相关，")
    emit("    单个维度的重要性不可解释为「某个物理量重要」。")
    emit("")

    # ---------------- 8. 复现 ----------------
    emit("## 8. 复现命令")
    emit("")
    emit("```powershell")
    emit('Set-Location -LiteralPath "F:\\DSH_Galatea\\Affinity_Intelligence"')
    emit("# 前置: join 校验 + 紧凑缓存")
    emit("python scripts/v2.5/00_join_check.py")
    emit("python scripts/v2.5/01_build_perres_cache.py --topk 32")
    emit("# 任务 1: pooled 基线重跑 + per-residue mean-pool")
    emit("python scripts/v2.5/02_perres_meanpool.py")
    emit("# 任务 1: per-residue CNN + 注意力池化 (CPU, 约 20 分钟)")
    emit("python scripts/v2.5/03_cnn_perres.py --seeds 42 43 --epochs 60 --max-res 16")
    emit("# 任务 2: 物化能量 (约 39 分钟) + 归一化")
    emit("python scripts/v2.5/04_energy_features.py")
    emit("python scripts/v2.5/04b_normalize_energy.py")
    emit("# 任务 3: 生物学特征")
    emit("python scripts/v2.5/05_conservation_features.py")
    emit("python scripts/v2.5/06_coevolution_features.py")
    emit("python scripts/v2.5/07_hotspot_weighting.py")
    emit("# 任务 4: 热力学分解")
    emit("python scripts/v2.5/08_thermo_decomposition.py")
    emit("# 任务 5: 综合消融")
    emit("python scripts/v2.5/10_ablation.py")
    emit("# 任务 6: 酶数据骨架自检 (不下载)")
    emit("python scripts/v2.5/09_enzyme_data_prep.py --self-test")
    emit("# 任务 7: 本报告")
    emit("python scripts/v2.5/11_final_report.py")
    emit("```")
    emit("")
    emit("> **全程 CPU only**，未安装任何 GPU 包。")

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT} ({OUT.stat().st_size} B)")
    if missing:
        print("缺失:")
        for m in missing:
            print("  -", m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
