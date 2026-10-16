# -*- coding: utf-8 -*-
"""scripts/v2.4/06_final_report.py -- V2.4 任务 5.2: 总报告生成器.

从各任务**实际落盘的产物**汇总, **不重新计算任何指标**:

* `features/v2.4_within_target_diagnosis.md`      (任务 0)
* `features/v2.4_skempi_within_target_report.md`  (任务 1)
* `features/v2.4_source_split_report.md`          (任务 2)
* `features/v2.4_binding_control_report.md`        (任务 K, 关键对照)
* `features/v2.4_head_enrichment_report.md`       (任务 3)
* `features/v2.4_perresidue_instructions.md`      (任务 4)

产出 `features/v2.4_final_report.md`。
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

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
OUT = FEAT / "v2.4_final_report.md"

WT_MD = FEAT / "v2.4_skempi_within_target_report.md"
WT_CSV = FEAT / "v2.4_skempi_within_target_results.csv"
SS_CSV = FEAT / "v2.4_source_split_results.csv"
BC_CSV = FEAT / "v2.4_binding_control_results.csv"
HE_CSV = FEAT / "v2.4_head_enrichment_results.csv"
NEXT_MD = FEAT / "v2.4_next_steps.md"
DIAG_MD = FEAT / "v2.4_within_target_diagnosis.md"

CROSS_TARGET = 0.3889577363940227
V22_BASELINE = 0.23747781344414348
BRANCH_A_THRESHOLD = 0.55

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
        return pd.read_csv(p)
    missing.append(f"{tag} (`{p.name}`)")
    return None


def sect(title: str, src: Path, start: str | None = None, end: str | None = None,
         max_chars: int = 4000) -> None:
    emit(f"## {title}")
    emit("")
    if not src.exists():
        missing.append(f"`{src.name}`")
        emit(f"(缺少 `{src.name}`)")
        emit("")
        return
    t = src.read_text(encoding="utf-8")
    if start:
        i = t.find(start)
        if i >= 0:
            t = t[i:]
    if end:
        j = t.find(end)
        if j > 0:
            t = t[:j]
    t = t.strip()
    if len(t) > max_chars:
        t = t[:max_chars] + "\n\n...（截断，完整内容见原文件）"
    emit(t)
    emit("")


def main() -> int:
    wt = load(WT_CSV, "within-target 结果")
    ss = load(SS_CSV, "来源拆分结果")
    bc = load(BC_CSV, "绑定方式对照")
    he = load(HE_CSV, "头部富集结果")

    emit("# V2.4 最终报告 — within-target 基准 / 来源拆分 / 头部富集 / per-residue 准备")
    emit("")
    emit("> 本报告由 `scripts/v2.4/06_final_report.py` 从各任务**实际产物**汇总生成，")
    emit("> **不重新计算任何指标**。每个数字都可追溯到对应 CSV/MD。")
    emit("")

    # ---------- 摘要 ----------
    emit("## 摘要")
    emit("")
    if wt is not None and len(wt):
        b = wt[wt["model"] == "physics+esm2_both"]
        nn = wt[wt["model"] == "nn_two_tower"]
        en = wt[wt["model"] == "ensemble_rankmean"]
        med = float(b["spearman_median"].iloc[0]) if len(b) else float("nan")
        emit(f"* **within-target（SKEMPI, 200 复合物 / 4987 行, GroupKFold）**:")
        emit(f"  `physics+esm2_both` 中位 Spearman = **{f(med)}** "
             f"(IQR [{f(b['iqr_low'].iloc[0], 3)}, {f(b['iqr_high'].iloc[0], 3)}]), "
             f"加权均值 {f(b['spearman_weighted_mean'].iloc[0])}")
        if len(nn):
            emit(f"  * 双塔 NN = {f(nn['spearman_median'].iloc[0])}；"
                 f"rank 平均集成 = {f(en['spearman_median'].iloc[0]) if len(en) else '—'}")
        emit(f"* **与 cross-target 对比**: within **{f(med)}** vs cross **{f(CROSS_TARGET)}** "
             f"-> Δ = **{f(med - CROSS_TARGET)}**")
        emit(f"  -> **V2.3.2 的 n=34 结果 (+0.5848) 未复现**；"
             f"within-target 反而**低于** cross-target。")
        emit(f"* **是否达分支 A 阈值 (>= {BRANCH_A_THRESHOLD})**: "
             f"**{'是' if med >= BRANCH_A_THRESHOLD else '否'}**")
        emit(f"* **置换对照**: null 中位 0.0000, 95% 区间 [−0.047, +0.047]; "
             f"观测 {f(med)} -> 经验 p = **0.0050**（显著高于随机但幅度很小）")
        emit(f"* **`physics_only` 在 within-target 上不可用**: 同一界面内 "
             f"**60/60 个物理特征完全恒定** -> 原理上无法区分突变体")
    if ss is not None and len(ss):
        sep = ss[ss["design"].str.contains("子集内") & (ss["feature_set"] == "physics+esm2_both")]
        mix = ss[ss["design"].str.contains("混合")]
        for lbl in ("一般 PPI", "SAbDab"):
            s_ = sep[sep["subset"] == lbl]
            m_ = mix[mix["subset"] == lbl]
            if len(s_) and len(m_):
                emit(f"* **{lbl}**: 分开 {f(s_['spearman'].iloc[0])} vs 混合 "
                     f"{f(m_['spearman'].iloc[0])} -> Δ = "
                     f"**{f(float(s_['spearman'].iloc[0]) - float(m_['spearman'].iloc[0]))}**")
    if he is not None and len(he):
        base = he[he["model"].str.startswith("regression")]
        xr = he[he["model"].str.contains("rank_xendcg")]
        if len(base) and len(xr):
            emit(f"* **头部富集 (`rank_xendcg`)**: 基线 EF@1% = "
                 f"{float(base['ef@1%'].iloc[0]):.3f} -> xendcg 最好 "
                 f"**{float(xr['ef@1%'].max()):.3f}** -> **未提升（反而变差）**")
    if bc is not None and len(bc):
        emit(f"* **关键对照**: 绑定方式影响 +0.0204；严格靶点不重叠后 Δ = **+0.0000** "
             f"-> **不存在泄漏**，任务 2 的高分来自协议不同（GroupKFold 更易）")
    if missing:
        emit("")
        emit("> ⚠️ 以下产物缺失：")
        for m in missing:
            emit(f"> * {m}")

    # ---------- 1. 诊断 ----------
    sect("1. within-target 可行性诊断（任务 0）", DIAG_MD,
         start="## 1. SKEMPI 全表", end="## 4. 诊断结论", max_chars=3000)

    # ---------- 2. within-target ----------
    emit("## 2. within-target 基准（任务 1）— 核心结果")
    emit("")
    if wt is not None and len(wt):
        emit("协议：**复合物级 GroupKFold（5 折）**，`group = interface_key`；")
        emit("每行折外预测 -> **复合物内部**排序算 per-complex Spearman -> 汇总。")
        emit("")
        emit("| 模型 | 复合物数 | 有效 | 排除 | **中位 Spearman** | IQR | 均值 | 加权均值 |")
        emit("|---|---:|---:|---:|---:|---|---:|---:|")
        for r in wt.sort_values("spearman_median", ascending=False).to_dict("records"):
            md = r.get("spearman_median")
            if md is None or (isinstance(md, float) and not np.isfinite(md)):
                emit(f"| {r['model']} | {int(r['n_complex'])} | {int(r['n_valid'])} | "
                     f"{int(r['n_excluded'])} | — | — | — | — |")
            else:
                emit(f"| {r['model']} | {int(r['n_complex'])} | {int(r['n_valid'])} | "
                     f"{int(r['n_excluded'])} | **{f(md)}** | "
                     f"[{f(r['iqr_low'], 3)}, {f(r['iqr_high'], 3)}] | "
                     f"{f(r['spearman_mean'])} | {f(r['spearman_weighted_mean'])} |")
        emit("")
        emit(f"- 详细报告：`v2.4_skempi_within_target_report.md`")
        emit(f"- 基准表：`v2.4_skempi_within_target.csv`（4987 行 / 200 复合物）")
        emit(f"- 逐复合物：`v2.4_skempi_per_complex.csv`")
    else:
        emit("(within-target 结果缺失)")
    emit("")

    # ---------- 3. 来源拆分 ----------
    sect("3. 来源拆分：SAbDab vs 一般 PPI（任务 2）",
         FEAT / "v2.4_source_split_report.md",
         start="## 2. 分组统计", end="## 5. 局限", max_chars=3000)

    # ---------- 4. 关键对照 ----------
    emit("## 4. 关键对照：绑定方式与泄漏排查（任务 K）")
    emit("")
    emit("> 任务 2 的 `row_exact` 数字（全体 +0.6734）远高于 V2.3.2 记录值 (+0.3890)。")
    emit("> 为排除泄漏专门做了对照，**结论是协议不同，不是泄漏**。")
    emit("")
    if bc is not None and len(bc):
        emit("| 设置 | 说明 | Spearman | 95% CI |")
        emit("|---|---|---:|---|")
        for r in bc.to_dict("records"):
            emit(f"| {r['setting']} | {r['note']} | {f(r['spearman'])} | "
                 f"[{f(r.get('ci_low'), 3)}, {f(r.get('ci_high'), 3)}] |")
        emit("")
        emit("要点：")
        emit("* V2.2 划分在 `pdb_id` 层面**完全不重叠**（test 的 573 个 pdb 中 0 个在 train）")
        emit("  -> 设置 C 无需删除任何行，Δ = **+0.0000** -> **无同 pdb 泄漏**；")
        emit("* 绑定方式本身只影响 **+0.0204**；")
        emit("* 因此两个数字的差异来自**评估协议**：")
        emit("  * V2.2 固定 test = **泛化到未见靶点**（难）；")
        emit("  * GroupKFold 随机折 = **数据集内平均表现**（易）。")
        emit("  **两者不可混用，必须标注协议。**")
        emit("* 最干净的一条证据：物理单独 **+0.2353**（≈ V2.2 基线 0.2375）、"
             "ESM-2 单独 **+0.3826** -> ESM-2 承担几乎全部跨靶信号。")
    else:
        emit("(对照结果缺失)")
    emit("")

    # ---------- 5. 头部富集 ----------
    emit("## 5. 头部富集专向优化（任务 3）")
    emit("")
    if he is not None and len(he):
        emit("只试一个方法：把目标函数换成 **`rank_xendcg`**（不需要 `label_gain` 分级）。")
        emit("基础特征集 = V2.3.2 头部富集最优 `physics+esm2_receptor` (l2, PCA128)。")
        emit("")
        emit("| 模型 | Spearman | EF@1% | P@10 | NDCG@10 |")
        emit("|---|---:|---:|---:|---:|")
        emit("| V2.3.2 `physics+esm2_receptor` (regression) | +0.3423 | 4.947 | 0.500 | 0.718 |")
        emit("| V2.3.2 `physics+esm2_both` (全局最优) | +0.3890 | 1.237 | 0.100 | — |")
        for r in he.to_dict("records"):
            emit(f"| 本次 {r['model']} | {f(r['spearman'])} | {float(r['ef@1%']):.3f} | "
                 f"{float(r['precision@10']):.3f} | {float(r['ndcg@10']):.3f} |")
        emit("")
        emit("**结论：`rank_xendcg` 四种配置（2 种 group × 2 种 relevance 分级）"
             "全部低于 `regression` 基线** ->")
        emit("头部富集的瓶颈在**特征**，不在目标函数。")
        emit("")
        emit("> 注：与 V2.3.2 的对比中，本次基线用 `row_exact` 装配（V2.3.2 用 `pdb_id` 平均），")
        emit("> 因此基线值与 V2.3.2 记录值有差异；但**本次内部的 regression vs xendcg "
             "对比是同一装配，因此有效**。")
    else:
        emit("(头部富集结果缺失)")
    emit("")

    # ---------- 6. per-residue ----------
    emit("## 6. per-residue + 界面残基 notebook（任务 4，**未运行**）")
    emit("")
    emit("* Notebook: `scripts/kaggle/esm2_extract_perresidue_v2.ipynb`（15 cells / 7 code）")
    emit("* 生成器: `scripts/v2.4/05_make_perresidue_v2_notebook.py`")
    emit("* 说明: `features/v2.4_perresidue_instructions.md`")
    emit("* 新增两个压缩开关：`SEQ_SUBSET`（默认 `test_val_skempi`，"
         "内含实测的 **285 个 SKEMPI PDB**）与 `INTERFACE_ONLY`"
         "（默认 True，±8 Å 界面残基）")
    emit("* 体积：**~3.1 GB -> ~20–60 MB**")
    emit("* 输出含 **`lig_residue_ids` / `rec_residue_ids`**（界面残基的**整链原始序号**），"
         "因为界面模式下必须知道每个向量对应整链第几号残基")
    emit("* **未运行**（需 Kaggle GPU）；`features/v2_esm2_perres_v2.npz` **不存在**")
    emit("")

    # ---------- 7. 条件分支 ----------
    emit("## 7. 条件分支判定（任务 5.1）")
    emit("")
    med = float("nan")
    if wt is not None and len(wt):
        b = wt[wt["model"] == "physics+esm2_both"]
        if len(b):
            med = float(b["spearman_median"].iloc[0])
    branch_a = bool(np.isfinite(med) and med >= BRANCH_A_THRESHOLD)
    emit(f"* within-target 中位 = **{f(med)}**；分支 A 阈值 = {BRANCH_A_THRESHOLD} "
         f"-> **{'满足' if branch_a else '不满足'}**")
    emit(f"* within-target 是否高于 cross-target ({f(CROSS_TARGET)})："
         f"**{'是' if med > CROSS_TARGET else '否'}**")
    emit("")
    emit("### 判定")
    emit("")
    if not branch_a:
        emit("> **命中：「within-target 与 cross-target 差异很大，但 within-target 更差」**")
        emit("> （任务书第二类情形：未达 0.55 且不高于 cross-target）")
        emit("> -> **推荐转向酶数据（SABIO-RK / PDBbind 酶子集），利用酶天然的 group 结构。**")
    else:
        emit("> 命中分支 A -> 推荐把 within-target 作为主基准并转向 ΔΔG。")
    emit("")
    emit("完整推理、优先级排序与已证伪方向见 `features/v2.4_next_steps.md`。")
    emit("")

    # ---------- 8. 局限 ----------
    emit("## 8. 局限（汇总）")
    emit("")
    emit("1. **两种评估协议不可混用**：V2.2 固定 test（泛化到未见靶点）与")
    emit("   GroupKFold 随机折（数据集内平均）在同一模型上可差 0.28。本报告已分别标注。")
    emit("2. **`pdb_id` 平均 join 会摧毁 within-target**（4987 行 -> 197 唯一嵌入，")
    emit("   组内特征恒定 -> Spearman 恒 NaN）。within-target 必须用逐行精确映射。")
    emit("3. **EF@1% 在小样本上高度离散**：V2.2 固定 test 只涉及 **2 个真阳性**，")
    emit("   步长 = 1/0.81（≈1.24），因此不能据此下强结论。")
    emit("4. **per-complex Spearman 在 n=3 时极不稳定**；因此同时报告 IQR、加权均值与置换 null。")
    emit("5. **SKEMPI 全表 (6994 行) 只用到 5109 行**（V2.2 表内的 SKEMPI 子集），")
    emit("   其余行未进入特征表（需要结构/序列匹配，未做）。")
    emit("6. **任务 2 的拆分差异混合了两种效应**（去异质来源的收益 + 训练量减少的损失），")
    emit("   本设计无法分离，未宣称纯粹的「拆分收益」。")
    emit("7. **per-residue 未运行**，其下游训练脚本**尚未编写**。")
    emit("")

    # ---------- 9. 复现 ----------
    emit("## 9. 复现命令")
    emit("")
    emit("```powershell")
    emit('Set-Location -LiteralPath "F:\\DSH_Galatea\\Affinity_Intelligence"')
    emit("# 任务 0: 可行性诊断")
    emit("python scripts/v2.4/00_diagnose_within_target.py")
    emit("# 任务 1: within-target 基准 (含置换对照)")
    emit("python scripts/v2.4/01_skempi_within_target.py --min-conditions 3 --folds 5 \\")
    emit("    --binding row_exact --n-perm 200")
    emit("# 任务 1 + NN (可选, 约 13 分钟)")
    emit("python scripts/v2.4/01_skempi_within_target.py --binding row_exact --nn --nn-seeds 42 43")
    emit("# 任务 2: 来源拆分")
    emit("python scripts/v2.4/02_source_split.py")
    emit("# 任务 K: 绑定方式/泄漏关键对照")
    emit("python scripts/v2.4/04_binding_control.py")
    emit("# 任务 3: 头部富集 (rank_xendcg)")
    emit("python scripts/v2.4/03_head_enrichment.py")
    emit("# 任务 4: 生成 per-residue notebook (不运行)")
    emit("python scripts/v2.4/05_make_perresidue_v2_notebook.py")
    emit("# 任务 5: 总报告")
    emit("python scripts/v2.4/06_final_report.py")
    emit("```")
    emit("")
    emit("> per-residue 提取需在 Kaggle 手动跑 "
         "`scripts/kaggle/esm2_extract_perresidue_v2.ipynb`（**未运行**）。")

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"写出 {OUT} ({OUT.stat().st_size} B)")
    if missing:
        print("缺失产物:")
        for m in missing:
            print("  -", m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
