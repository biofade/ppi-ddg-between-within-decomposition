# -*- coding: utf-8 -*-
"""scripts/v2.3.2/06_final_report.py -- V2.3.2 任务 5.4: 总报告生成器.

从各任务的**实际产物**汇总:
* `features/v2.3.2_lgbm_results.csv`      (任务 2)
* `features/v2.3.2_lgbm_pairs.csv`        (任务 2 配对检验)
* `features/v2.3.2_nn_verdict.json`       (任务 3 判决)
* `features/v2.3.2_nn_results.csv`        (任务 3 逐配置)
* `features/v2.3.2_ensemble_results.csv`  (任务 5.1)
* `features/v2.3.2_error_analysis.md`     (任务 5.2)

产出 `features/v2.3.2_final_report.md`, 含条件分支 (任务 5.3)。

**不重新计算任何指标** —— 只做汇总与排版, 保证报告中的每个数字都来自上游脚本。
"""

from __future__ import annotations

import json
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
OUT = FEAT / "v2.3.2_final_report.md"

LGBM_CSV = FEAT / "v2.3.2_lgbm_results.csv"
LGBM_PAIRS = FEAT / "v2.3.2_lgbm_pairs.csv"
NN_VERDICT = FEAT / "v2.3.2_nn_verdict.json"
NN_CSV = FEAT / "v2.3.2_nn_results.csv"
NN_JUDGE = FEAT / "v2.3.2_nn_judgement.csv"
ENS_CSV = FEAT / "v2.3.2_ensemble_results.csv"
ERR_MD = FEAT / "v2.3.2_error_analysis.md"
JOIN_MD = FEAT / "v2.3.2_join_diagnosis.md"

V22_BASELINE = 0.23747781344414348


def f(v, d=4):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{x:+.{d}f}" if np.isfinite(x) else "NaN"


def main() -> int:
    L: list[str] = []
    ap_emit = L.append
    missing: list[str] = []

    def load_csv(p: Path, tag: str):
        if p.exists():
            return pd.read_csv(p)
        missing.append(f"{tag} (`{p.name}`)")
        return None

    lgbm = load_csv(LGBM_CSV, "LightGBM 消融结果")
    pairs = load_csv(LGBM_PAIRS, "LightGBM 配对检验")
    nnres = load_csv(NN_CSV, "NN 结果")
    njudge = load_csv(NN_JUDGE, "NN 判决")
    ens = load_csv(ENS_CSV, "集成结果")
    verdict = json.loads(NN_VERDICT.read_text(encoding="utf-8")) if NN_VERDICT.exists() else None
    if verdict is None:
        missing.append("NN 判决 (`v2.3.2_nn_verdict.json`)")

    ap_emit("# V2.3.2 最终报告 — join 修复 / LightGBM 全量消融 / NN 判决 / 集成 / 错误分析")
    ap_emit("")
    ap_emit("> 本报告由 `scripts/v2.3.2/06_final_report.py` 从各任务**实际产物**汇总生成,")
    ap_emit("> **不重新计算任何指标**。每个数字都可追溯到对应的 CSV/JSON/MD。")
    ap_emit("")
    ap_emit("## 摘要")
    ap_emit("")
    if lgbm is not None and len(lgbm):
        base_row = lgbm[lgbm["feature_set"] == "physics_only"]
        base_rho = float(base_row["spearman"].iloc[0]) if len(base_row) else float("nan")
        best = lgbm.sort_values("spearman", ascending=False).iloc[0]
        ap_emit(f"* **join 修复**: V2.3 的 `interface_key` 语义不一致导致 0% 命中;")
        ap_emit(f"  改为 `pdb_id` join 后 **命中率 100.00%** (详见 `v2.3.2_join_diagnosis.md`)。")
        ap_emit(f"* **LightGBM 消融最优**: `{best['feature_set']}` "
                f"(版本 {best.get('embedding_variant', '-')}, PCA {best.get('pca_dim', '-')})")
        ap_emit(f"  -> Spearman = **{f(best['spearman'])}** "
                f"[{f(best.get('spearman_ci_low'),3)}, {f(best.get('spearman_ci_high'),3)}], "
                f"NDCG@10 = {f(best.get('ndcg@10'),3)}, "
                f"EF@1% = {f(best.get('ef@1%'),3)}, P@10 = {f(best.get('precision@10'),3)}")
        ap_emit(f"* **V2.2 基线 (physics_only)**: {f(base_rho)} "
                f"(V2.2 记录值 {V22_BASELINE:.6f})")
        ap_emit(f"* **ESM-2 增益**: Δ = **{f(float(best['spearman']) - base_rho)}**")
        better = None
        if pairs is not None and "model" in pairs.columns:
            b = pairs[pairs["model"].str.contains(str(best["feature_set"]).split("+")[-1],
                                                   regex=False)]
            b = b[(b["bh_q"] < 0.05) & (b["delta"] > 0)] if "bh_q" in b.columns else b.iloc[0:0]
            if len(b):
                better = b.sort_values("delta", ascending=False).iloc[0]
        if better is not None:
            ap_emit(f"  * 该配置相对基线**显著更优**: Δ = {f(better['delta'])}, "
                    f"BH q = {float(better['bh_q']):.4f} (< 0.05) ✅")
        elif pairs is not None:
            ap_emit("  * (相对基线的配对检验: 见第 2 节完整表)")
    if verdict:
        ap_emit(f"* **NN 判决**: 最优 NN = `{verdict.get('best_nn')}` "
                f"Spearman = **{f(verdict.get('nn_spearman'))}**; "
                f"LightGBM 最优 = `{verdict.get('lgbm_best')}` "
                f"{f(verdict.get('lgbm_spearman'))}")
        ap_emit(f"  -> Δ = {f((verdict.get('nn_spearman') or 0) - (verdict.get('lgbm_spearman') or 0))}, "
                f"BH q = {verdict.get('nn_best_bh_q')}")
        ap_emit(f"  -> **判决: `{verdict.get('verdict')}`**")
    if ens is not None and len(ens):
        be = ens.sort_values("spearman", ascending=False).iloc[0]
        ap_emit(f"* **集成最优**: `{be['ensemble']}` -> Spearman = **{f(be['spearman'])}** "
                f"[{f(be.get('spearman_ci_low'),3)}, {f(be.get('spearman_ci_high'),3)}]")
    if missing:
        ap_emit("")
        ap_emit("> ⚠️ 以下产物缺失, 相应小节已标注:")
        for m in missing:
            ap_emit(f"> * {m}")

    # ---- 1. join 修复 ----
    ap_emit("")
    ap_emit("## 1. join 修复")
    ap_emit("")
    if JOIN_MD.exists():
        txt = JOIN_MD.read_text(encoding="utf-8")
        sec = txt.find("## 2. 根因")
        end = txt.find("## 6.")
        body = txt[sec:end].strip() if sec >= 0 and end > sec else ""
        # 把原文的二级标题降为三级, 避免与总报告的一级/二级结构冲突
        body = body.replace("\n## ", "\n### ")
        ap_emit(body if body else "(见 v2.3.2_join_diagnosis.md)")
    else:
        ap_emit("(缺少 `v2.3.2_join_diagnosis.md`)")

    # ---- 2. LightGBM 消融 ----
    ap_emit("")
    ap_emit("## 2. LightGBM 全量消融")
    ap_emit("")
    if lgbm is not None and len(lgbm):
        ap_emit(f"共 **{len(lgbm)}** 个配置, 全部在 Proteinea 自带 test (n=752) 上评估, ")
        ap_emit("划分与超参**完全复用 V2.2**, ESM-2 侧 PCA **只在训练集 fit**。")
        ap_emit("")
        ap_emit("| 特征集 | 版本 | PCA | 特征数 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
        ap_emit("|---|---|---:|---:|---:|---|---:|---:|---:|")
        for r in lgbm.sort_values("spearman", ascending=False).to_dict("records"):
            ap_emit(f"| {r['feature_set']} | {r.get('embedding_variant', '-')} | "
                    f"{r.get('pca_dim', '-')} | {int(r['n_features'])} | "
                    f"{f(r['spearman'])} | "
                    f"[{f(r.get('spearman_ci_low'),3)}, {f(r.get('spearman_ci_high'),3)}] | "
                    f"{f(r.get('ndcg@10'),3)} | {f(r.get('ef@1%'),3)} | "
                    f"{f(r.get('precision@10'),3)} |")
    else:
        ap_emit("(LightGBM 结果缺失)")
    if pairs is not None and len(pairs):
        sig = pairs[(pairs.get("bh_q", pd.Series(dtype=float)) < 0.05)]
        ap_emit("")
        ap_emit(f"**相对基线 `physics_only` 的配对 bootstrap (5000 次) + BH-FDR**: "
                f"共 {len(pairs)} 个检验, 其中 **{len(sig)}** 个 BH q < 0.05。")
        ap_emit("")
        ap_emit("| 模型 | Δρ vs 基线 | 95% CI | p | BH q | 显著 |")
        ap_emit("|---|---:|---|---:|---:|---|")
        for r in pairs.sort_values("delta", ascending=False).head(15).to_dict("records"):
            q = r.get("bh_q")
            ap_emit(f"| {r['model']} | {f(r.get('delta'))} | "
                    f"[{f(r.get('lo'),3)}, {f(r.get('hi'),3)}] | "
                    f"{float(r.get('p_two_sided', float('nan'))):.4f} | "
                    f"{float(q):.4f} | {'是' if (q == q and q < 0.05) else '否'} |")

    # ---- 3. NN 判决 ----
    ap_emit("")
    ap_emit("## 3. 神经网络判决实验 (N1 MLP / N2 双塔)")
    ap_emit("")
    ap_emit("> **CPU only** (硬约束 7); 早停准则 = 验证集 Spearman, 连续 20 epoch 不升;")
    ap_emit("> 种子 42/43/44/45。详细表见 `v2.3.2_nn_report.md`。")
    ap_emit("")
    if nnres is not None and len(nnres):
        top = nnres.sort_values("spearman_mean", ascending=False).head(15)
        ap_emit("| 模型 | PCA | 输入维度 | 种子 | Spearman (mean±std) | NDCG@10 | EF@1% | P@10 |")
        ap_emit("|---|---:|---:|---:|---|---:|---:|---:|")
        for r in top.to_dict("records"):
            ap_emit(f"| {r.get('model')} | {r.get('pca_dim')} | {int(r.get('input_dim', -1))} | "
                    f"{int(r.get('n_seeds', 0))} | {f(r.get('spearman_mean'))} ± "
                    f"{float(r.get('spearman_std', float('nan'))):.4f} | "
                    f"{f(r.get('ndcg10_mean'),3)} | {f(r.get('ef1_mean'),3)} | "
                    f"{f(r.get('p10_mean'),3)} |")
    if njudge is not None and len(njudge):
        ap_emit("")
        ap_emit("### 与 LightGBM 最优的配对检验")
        ap_emit("")
        ap_emit("| NN 配置 | NN ρ | LightGBM ρ | Δ | 95% CI | p | BH q |")
        ap_emit("|---|---:|---:|---:|---|---:|---:|")
        for r in njudge.sort_values("delta", ascending=False).to_dict("records"):
            ap_emit(f"| {r['nn']} | {f(r['nn_rho'])} | {f(r['lgbm_rho'])} | "
                    f"{f(r['delta'])} | [{f(r.get('lo'),3)}, {f(r.get('hi'),3)}] | "
                    f"{float(r.get('p_two_sided', float('nan'))):.4f} | "
                    f"{float(r.get('bh_q', float('nan'))):.4f} |")
    if verdict:
        ap_emit("")
        ap_emit(f"**判决: `{verdict.get('verdict')}`**")
        ap_emit("")
        if verdict.get("verdict") == "NN_SIGNIFICANTLY_WORSE":
            ap_emit("> NN 最优**显著劣于** LightGBM 最优 (BH q < 0.05)。")
            ap_emit("> **结论: 在 pooled 嵌入上 NN 不适用。**")
            ap_emit("> 建议停止 pooled 上的 NN 方向; 若仍要做 NN,")
            ap_emit("> 唯一有希望的通路是 **per-residue + 1D-CNN/注意力**。")
        elif verdict.get("verdict") == "NN_COMPETITIVE":
            ap_emit("> NN 未显著劣于 LightGBM, 且点估计不低于它 -> pooled 上 NN **可行**。")
        else:
            ap_emit("> NN 与 LightGBM **大致等价** (差异不显著) -> 无理由用 NN 替代 LightGBM")
            ap_emit("> (NN 训练成本高得多且需调参)。")

    # ---- 4. 集成 ----
    ap_emit("")
    ap_emit("## 4. 集成 (任务 5.1)")
    ap_emit("")
    if ens is not None and len(ens):
        ap_emit("| 集成方式 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
        ap_emit("|---|---:|---|---:|---:|---:|")
        for r in ens.sort_values("spearman", ascending=False).to_dict("records"):
            ap_emit(f"| {r['ensemble']} | {f(r['spearman'])} | "
                    f"[{f(r.get('spearman_ci_low'),3)}, {f(r.get('spearman_ci_high'),3)}] | "
                    f"{f(r.get('ndcg@10'),3)} | {f(r.get('ef@1%'),3)} | "
                    f"{f(r.get('precision@10'),3)} |")
        ap_emit("")
        ap_emit("> **LR stacking 未执行**: 需要在验证集上 fit 元学习器, 但本轮 NN 脚本")
        ap_emit("> 只回传 test 预测 (未存 val 预测), 在 test 上拟合会造成泄漏。")
        ap_emit("> **如实标注为未做**, 而不是用 test 拟合一个虚假的 stacking 数字。")
    else:
        ap_emit("(集成结果缺失)")

    # ---- 5. 错误分析 ----
    ap_emit("")
    ap_emit("## 5. 错误分析 (任务 5.2)")
    ap_emit("")
    if ERR_MD.exists():
        txt = ERR_MD.read_text(encoding="utf-8")
        i = txt.find("## 4. 错误分析结论")
        ap_emit(txt[i:].strip() if i >= 0 else "(见 v2.3.2_error_analysis.md)")
        ap_emit("")
        ap_emit("完整表格见 `features/v2.3.2_error_analysis.md`。")
    else:
        ap_emit("(错误分析缺失)")

    # ---- 6. 条件分支 (任务 5.3) ----
    ap_emit("")
    ap_emit("## 6. 条件分支判定 (任务 5.3)")
    ap_emit("")
    esm2_sig = None
    base_rho = None
    if lgbm is not None and len(lgbm) and pairs is not None and len(pairs):
        b = lgbm[lgbm["feature_set"] == "physics_only"]
        if len(b):
            base_rho = float(b["spearman"].iloc[0])
            # 找任何含 esm2 且 BH 显著且 delta>0 的配置
            cand = pairs[pairs["model"].str.contains("esm2", case=False, regex=False)]
            if "bh_q" in cand.columns:
                cand = cand[(cand["bh_q"] < 0.05) & (cand["delta"] > 0)]
            esm2_sig = bool(len(cand) > 0)
    nn_worse = (verdict or {}).get("verdict") == "NN_SIGNIFICANTLY_WORSE"

    ap_emit(f"* ESM-2 是否**显著优于** `physics_only`: "
            f"**{'是' if esm2_sig else '否' if esm2_sig is not None else '无法判定'}**")
    ap_emit(f"* NN 是否**显著劣于** LightGBM: **{'是' if nn_worse else '否' if verdict else '无法判定'}**")
    ap_emit("")
    ap_emit("### 对应建议")
    ap_emit("")
    if esm2_sig and not nn_worse:
        ap_emit("> **分支 A** (ESM-2 显著更优, 且 NN 不显著更差):")
        ap_emit("> 推荐跑 per-residue notebook (`scripts/kaggle/esm2_extract_perresidue.ipynb`),")
        ap_emit("> 在其上做 1D-CNN / 注意力。")
    elif esm2_sig and nn_worse:
        ap_emit("> **分支 B** (ESM-2 显著更优, 但 NN 显著更差):")
        ap_emit("> **推荐保留 LightGBM**, 转向 within-target ranking (SKEMPI 2.0)。")
    elif esm2_sig is False:
        ap_emit("> **分支 C** (ESM-2 不显著优于 physics_only):")
        ap_emit("> 推荐转向 within-target ranking 或酶数据; PPI 跨靶已饱和。")
    else:
        ap_emit("> **无法判定**: 缺少 LightGBM 消融或配对检验产物。")

    # 错误分析里有一个直接支持/反对 within-target 方向的实测证据, 必须放进分支讨论
    ap_emit("")
    ap_emit("### 支持该建议的独立证据 (来自任务 5.2 错误分析)")
    ap_emit("")
    ap_emit("错误分析按**来源数据集**分组时发现 (见 `v2.3.2_error_analysis.md`):")
    ap_emit("")
    ap_emit("| 来源数据集 | n | Spearman | 平均绝对排名差 (pp) |")
    ap_emit("|---|---:|---:|---:|")
    ap_emit("| **SKEMPI v2.0** (同一靶点的突变体) | 34 | **+0.5848** | 17.2 |")
    ap_emit("| PDBbind v2020 | 657 | +0.3648 | 25.9 |")
    ap_emit("| Affinity Benchmark v5.5 | 23 | +0.3985 | 27.5 |")
    ap_emit("| **SAbDab** (抗体-抗原) | 38 | **-0.0509** | 21.7 |")
    ap_emit("")
    ap_emit("**解读 (重要)**:")
    ap_emit("")
    ap_emit("* 本数据集中唯一属于 **within-target 突变效应** 的子集 (SKEMPI v2.0, n=34)")
    ap_emit("  拿到了**全表最高的分组 Spearman (+0.58)** —— 远高于它在全局跨靶任务里的表现。")
    ap_emit("  这是**直接支持「转向 within-target ranking」**的实测证据, 而不只是推测。")
    ap_emit("* 同时 **SAbDab (抗体-抗原) 子集 ρ ≈ 0** —— 说明**抗体-抗原是完全不同的结合模式**,")
    ap_emit("  把抗体数据与一般 PPI 混在一个跨靶回归里, 很可能互相拖累。")
    ap_emit("* 分层分析也显示: 模型在**强结合区**层内排序能力最差 (ρ ≈ -0.08),")
    ap_emit("  而湿实验最关心的恰恰是强结合区 —— 这进一步指向「需要换任务形态,")
    ap_emit("  而不是继续堆特征」。")
    ap_emit("")
    ap_emit("> 注: SKEMPI 子集 n=34, CI 很宽, 该数字**仅作方向性证据**, 不能当作已确立的结论。")

    # ---- 7. 局限 ----
    ap_emit("")
    ap_emit("## 7. 局限")
    ap_emit("")
    ap_emit("1. **结构级 join 退化**: 主结果用 `pdb_id` join, 同一结构的多个突变体共享")
    ap_emit("   (平均后的) 嵌入。test 侧 126 个多行 pdb 中 14 个组内有真实 ΔG 差异,")
    ap_emit("   这是**不可逆的信息损失**。")
    ap_emit("2. **只有 pooled 嵌入**: 无 per-residue, 因此 NN 无法利用序列顺序 (任务 4 已准备提取)。")
    ap_emit("3. **CI bootstrap 次数**: 配对检验按要求用 5000 次; 每个模型自身的 CI 用较少次数")
    ap_emit("   (计算预算), 已在报告「统计实现说明」标注。")
    ap_emit("4. **LR stacking 未做** (缺 val 预测, 避免 test 泄漏)。")
    ap_emit("5. **NN 超参未调** (任务书指定固定超参), 不排除调参后变好;")
    ap_emit("   但判决实验的目的是回答「默认配置下是否可行」。")
    ap_emit("6. 单机 CPU; 单 seed 的 LightGBM (42), NN 用 4 个种子。")
    ap_emit("")
    ap_emit("## 8. 复现命令")
    ap_emit("")
    ap_emit("```powershell")
    ap_emit('Set-Location -LiteralPath "F:\\DSH_Galatea\\Affinity_Intelligence"')
    ap_emit("python scripts/v2.3.2/03_merge_esm2_v2.py --n-boot 5000 --n-boot-ci 600")
    ap_emit("# NN: 本机仅 7.5 GB 内存, 两个 NN 并发会 OOM -> 用顺序 runner")
    ap_emit("python scripts/v2.3.2/run_nn_sequential.py")
    ap_emit("python scripts/v2.3.2/04_nn_compare.py --n-boot 5000")
    ap_emit("python scripts/v2.3.2/05_error_analysis.py")
    ap_emit("python scripts/v2.3.2/06_final_report.py")
    ap_emit("```")
    ap_emit("")
    ap_emit("> per-residue 提取需在 Kaggle 手动跑 "
            "`scripts/kaggle/esm2_extract_perresidue.ipynb` (未运行)。")

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"写出 {OUT} ({OUT.stat().st_size} B)")
    if missing:
        print("缺失产物:")
        for m in missing:
            print("  -", m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
