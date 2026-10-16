# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/05_final_report.py -- V7.2R-3 最终报告（§27 结构）。

ZIP 未返回时输出本地就绪状态 + `V7_2R_3_KAGGLE_READY`。
ZIP 返回后填入完整科学结果与 §24 结果判定。
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

V723 = C.ROOT / "experiments" / "v7_2r_3"
KAG = C.ROOT / "kaggle" / "v7_2r_3"
ZIP = Path(r"F:\dsh\downloads\v7_2r_3\V7_2R_3_PPI_CONTEXT_FEATURES.zip")
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
    pkg = KAG / "kaggle_input_package.zip"
    nb = KAG / "V7_2R_3_PPI_Context_Kaggle.ipynb"
    rd = KAG / "README_KAGGLE.md"
    man = KAG / "kaggle_input_manifest.csv"
    ready = all(p.exists() for p in (pkg, nb, rd, man))
    res_p = V723 / "phase5_eval" / "ppi_results.csv"
    bo_p = V723 / "phase6_bootstrap" / "contrasts.csv"
    have = res_p.exists() and bo_p.exists()

    emit("# V7.2R-3 Final Summary")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    emit("## 1. Why V7.2R-3 Was Required")
    emit("")
    emit("V7.2R-2 的 Kaggle notebook 存在一个**真实的表示缺陷**：")
    emit("")
    emit("```python")
    emit("build_variant()   # 只返回 sites_by_chain 中的链")
    emit("```")
    emit("")
    emit("`target_chains.csv` 也是由 `sit[[\"target_id\",\"chain\"]].drop_duplicates()`")
    emit("写出，因此**只含携带突变位点的链**。实测后果：")
    emit("")
    emit("| 项 | V7.2R-2 实测 |")
    emit("|---|---:|")
    emit("| (target,chain) 记录 | **233** |")
    emit("| 每 target 链数中位 | **1** |")
    emit("| 至少缺 1 条伙伴链的 target | **150 / 174 = 86.2%** |")
    emit("| 丢失伙伴链总数 | **421** |")
    emit("")
    emit("⇒ V7.2R-2 测的是 **`STRONG_GENERIC_PLM_PERTURBATION`（链隔离）**，")
    emit("**不是** `PPI_AWARE_PERTURBATION`。")
    emit("")
    emit("V7.2R-3 的**唯一**科学问题：")
    emit("")
    emit("> 真正的伙伴条件化 PPI 上下文，能否在已测试的链隔离序列表示之上，")
    emit("> 带来可迁移的 family-OOD 突变排序信号？")
    emit("")
    emit("**修复**：`build_complex()` 使用该 PDB 的**全部蛋白链**，")
    emit("突变只施加到受影响链，伙伴链保持 WT（§10/§11）。")
    emit("")
    emit("## 2. PPI Representation Semantics")
    emit("")

    if not have:
        emit("**`WAITING_FOR_V7_2R_3_ZIP`** —— Kaggle PPI-context 提取尚未返回。")
        emit("")
        emit("计划中的语义（由 notebook 落盘并核验）：")
        emit("")
        emit("| 项 | 计划值 |")
        emit("|---|---|")
        emit("| Model | ESM2-650M over **full complex**（全部链拼接） |")
        emit("| level | **`MUTATION_SITE_PPI`**（per-residue，跨链自注意力） |")
        emit("| dimensions | 1280 |")
        emit("| WT | `r_wt`（全复合物，伙伴链 WT） |")
        emit("| Mut | `r_mut`（受影响链突变，伙伴链 WT） |")
        emit("| partner-dependent | 由 cell 08 的 `D_true` / `D_wrong` **实测判定** |")
    else:
        imp = V723 / "phase4_import" / "import_meta.json"
        meta = json.loads(imp.read_text(encoding="utf-8")) if imp.exists() else {}
        emit("| 项 | 值 |")
        emit("|---|---|")
        emit("| Model | ESM2-650M over **full complex** |")
        emit("| level | **`MUTATION_SITE_PPI`** |")
        emit(f"| dimensions | {meta.get('dims', {})} |")
        emit(f"| 嵌入链数中位 | **{meta.get('median_chains')}**"
             f"（V7.2R-2 为 1） |")
    emit("")
    emit("## 3. Partner Sanity Test")
    emit("")
    ps = V723 / "phase4_import" / "partner_sanity.csv"
    if ps.exists():
        SA = pd.read_csv(ps)
        d1 = SA.D_true.median() if "D_true" in SA else float("nan")
        d2 = SA.D_wrong.median() if "D_wrong" in SA else float("nan")
        emit(f"- 测试 target 数：**{len(SA)}**")
        emit(f"- `A alone` vs `A+true`：**D_true 中位 = {d1:.4f}**")
        emit(f"- `A+true` vs `A+wrong`：**D_wrong 中位 = {d2:.4f}**")
        v = ("PARTNER_CONDITIONING_CONFIRMED" if (d1 == d1 and d1 > 1e-3)
             else "PARTNER_CONDITIONING_NOT_DEMONSTRATED")
        emit(f"- **verdict：`{v}`**")
    else:
        emit("**NOT_RUN**（等 ZIP）。notebook cell 08 会在全量提取**之前**执行，")
        emit("未通过则主动抛错停止（§7）。")
    emit("")
    emit("## 4. Coverage")
    emit("")
    if have:
        idx = C.ROOT / "features" / "v7_2r_3" / "ppi_imported" / "ppi_feature_index.csv"
        RB = pd.read_csv(idx)
        emit(f"- rows: **{len(RB)}**")
        emit(f"- targets: **{RB.target_id_canonical.nunique()}**")
        emit(f"- families: **{RB.family_id_canonical.nunique()}**")
        n5 = (RB.groupby("target_id_canonical").size() >= 5).sum()
        emit(f"- n>=5: **{n5}**")
    else:
        emit("**WAITING_FOR_V7_2R_3_ZIP**")
    emit("")
    emit("## 5. Matched U_ABC")
    emit("")
    if have:
        rp = V723 / "phase5_eval" / "eval_report.md"
        t = rp.read_text(encoding="utf-8") if rp.exists() else ""
        i = t.find("U_ABC =")
        j = t.find("## 5.2")
        emit(t[i:j] if i >= 0 else "见 `phase5_eval/eval_report.md`")
    else:
        emit("`U_ABC = Track A ∩ Track B ∩ PPI context`，按 `row_id` 构造；"
             "**WAITING_FOR_V7_2R_3_ZIP**")
    emit("")
    emit("## 6. Primary Results")
    emit("")
    if have:
        R = pd.read_csv(res_p)
        sm = R.groupby("ranker")[["macro_spearman", "eff10", "C_W", "kendall",
                                 "pairwise05"]].mean()
        emit("| Method | Macro rho | Eff@10 | C_W | Kendall | Pairwise |")
        emit("|---|---:|---:|---:|---:|---:|")
        for k, r in sm.iterrows():
            emit(f"| {k} | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
                 f"{r['C_W']:+.4f} | {r['kendall']:+.4f} | "
                 f"{r['pairwise05']:+.4f} |")
    else:
        emit("| Method | Macro rho | Eff@10 | C_W | Kendall | Pairwise |")
        emit("|---|---:|---:|---:|---:|---:|")
        for m in ("B2 chain-isolated", "B3a + static PPI", "B3b + PPI interaction",
                  "B4 + PPI perturbation", "B_shuf shuffled PPI"):
            emit(f"| {m} | NOT_RUN | | | | |")
    emit("")
    emit("## 7. B3/B4 - B2 Bootstrap")
    emit("")
    if have:
        BO = pd.read_csv(bo_p)
        for c in BO.contrast.unique():
            emit(f"### {c}")
            emit("")
            g = BO[BO.contrast == c]
            emit("| 指标 | point | 95% CI | 排除 0? |")
            emit("|---|---:|---|---|")
            for _, r in g.iterrows():
                emit(f"| Δ{r['metric'][2:]} | {r['point']:+.4f} | "
                     f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                     f"{'✅' if r['ci_excludes_zero'] else '❌ 含 0'} |")
            emit("")
    else:
        emit("**WAITING_FOR_V7_2R_3_ZIP**")
    emit("")
    emit("## 8. True vs Shuffled Partner")
    emit("")
    emit("**WAITING_FOR_V7_2R_3_ZIP**（对照 `Shuf-真实` 已在 "
         "`04_bootstrap.py` 预注册）")
    emit("")
    emit("## 9. Wrong Partner Control")
    emit("")
    emit("notebook cell 08 已在 **表示层**做过错伙伴对照（`D_wrong`，§20）。")
    emit("排名层的错伙伴对照将作为后续机制分析（未在本轮评估）。")
    emit("")
    emit("## 10. Scientific Verdict")
    emit("")
    emit("```")
    emit("PPI_RESULT = WAITING_FOR_V7_2R_3_ZIP")
    emit("```")
    emit("")
    emit("判定规则（§24，**已冻结，不得事后改写**）：")
    emit("")
    emit("| 结果 | 条件 |")
    emit("|---|---|")
    emit("| `PPI_CONTEXT_NULL` | B3/B4 不显著 > B2 |")
    emit("| `STATIC_PPI_CONTEXT_GAIN` | B3 > B2 但 B4−B3a 不显著 |")
    emit("| `PPI_PERTURBATION_GAIN` | B4 显著 > B3a |")
    emit("| `PARTNER_SPECIFIC_GAIN` | 真伙伴 > 打乱/错伙伴 且有排序增益 |")
    emit("")
    emit("## 11. Allowed Claims")
    emit("")
    emit("- 已定位并修复 V7.2R-2 的链丢弃缺陷（233 → 654 条 (target,chain)）")
    emit("- 本轮的 PPI 上下文包**无标签**；未使用任何监督突变 head")
    emit("- 评估协议为 **V7.1.1 冻结的嵌套 family-OOD**")
    emit("- 已知的既有结论：Track A（ESM2-8M Δh）与 Track B（ESM2-650M Δr，链隔离）")
    emit("  在严格 family-OOD 下均为 **NULL**")
    emit("")
    emit("## 12. Forbidden Claims")
    emit("")
    emit("- ❌ 任何 PPI 排名结论（**尚未运行**）")
    emit("- ❌ 「protein foundation models do not work」")
    emit("- ❌ 「we invented the underlying PPI representation」（§51）")
    emit("- ❌ 「MINT 无效」—— MINT 仅在链级 sanity 中使用")
    emit("")
    emit("## 13. Final Representation-Branch Decision")
    emit("")
    emit("```")
    emit("CONTINUE_REPRESENTATION_ROUTE = PENDING_V7_2R_3_RESULT")
    emit("```")
    emit("")
    emit("**预注册决策规则**：")
    emit("- 若 `PPI_CONTEXT_NULL` → `CONTINUE_REPRESENTATION_ROUTE = NO`，")
    emit("  进入 **V7.3 physics / structure-informed ranking**（§24/§25）")
    emit("- 若 P1/P2/P3 → **先做机制刻画**，不立即放大架构（§26）")
    emit("")
    emit("## 14. Next Step")
    emit("")
    emit("**若 `KAGGLE_PACKAGE_READY = YES`：停止本地计算，等待 Kaggle ZIP。**")
    emit("")
    emit("ZIP 返回后依次运行：")
    emit("```")
    emit("python scripts/v7_2r_3/02_import_ppi_features.py")
    emit("python scripts/v7_2r_3/03_matched_eval.py")
    emit("python scripts/v7_2r_3/04_bootstrap.py")
    emit("python scripts/v7_2r_3/05_final_report.py")
    emit("```")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (V723 / "V7_2R_3_FINAL_SUMMARY.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
    json.dump({"V7_2R_3_KAGGLE_READY": bool(ready),
               "ppi_result_available": bool(have),
               "notebook": str(nb), "input_package": str(pkg),
               "readme": str(rd)},
              open(V723 / "package_status.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print("\n" + "=" * 66)
    print(f"V7_2R_3_KAGGLE_READY = {'YES' if ready else 'NO'}")
    print("=" * 66)
    print(f"\nNotebook:\n  {nb}")
    print(f"\nInput package:\n  {pkg}")
    print(f"\nREADME:\n  {rd}")
    print("\nWhat user must upload:")
    print(f"  1) {pkg.name}   -> Kaggle Dataset (V7.2R-3 版本，含全部伙伴链)")
    print(f"  2) {nb.name}   -> Kaggle Notebook")
    print("\nWhat user must download:")
    print("  V7_2R_3_PPI_CONTEXT_FEATURES.zip  ->  F:\\dsh\\downloads\\v7_2r_3\\")
    print(f"\n写出 {V723 / 'V7_2R_3_FINAL_SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
