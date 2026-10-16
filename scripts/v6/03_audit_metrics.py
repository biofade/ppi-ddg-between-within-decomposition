# -*- coding: utf-8 -*-
"""scripts/v6/03_audit_metrics.py -- V6.0 PHASE 2.2: 同协议重算 SOTA 指标.

对每个外部数据集的每个模型预测列，统一计算：
  1. Overall Pearson          （全局基准，SOTA 论文常报）
  2. Macro per-PDB Spearman   （靶点内排序能力；V6 核心）
  3. Hit@10 (ΔΔG<0)           （绝对命中率）
  4. Base Rate (ΔΔG<0)        （随机基准率）
  5. Efficiency@10            （(Random − Model)/(Random − Oracle)）
  6. n_complex / n_rows       （实际样本量）

过滤规则（V6 §2.2）：按 complex 分组，**丢弃组内突变数 < 5 的复合物**。
同时计算 **V5 baseline** 的对应指标作对照。

符号处理（V6 §0.2）：对每个预测列先按 `needs_flip`（V6 PHASE 2.1 判定）统一方向，
但**同时输出未反号版本**，不隐藏反号带来的人为提升。

输出
----
* `experiments/v6_sota_audit/sota_audit_metrics.csv`
* `experiments/v6_sota_audit/sota_audit_metrics.md`
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

OUT = C.ROOT / "experiments" / "v6_sota_audit"
OUT.mkdir(parents=True, exist_ok=True)
CATH_CSV = Path(r"F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv")
RBM_CSV = Path(r"F:\dsh\downloads\dssa_ppi\RBM_baselines_comparison_ddG.csv")
MIN_N = 5
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def eff10(y, p, g) -> dict:
    """Efficiency@10 = (Random − Model) / (Random − Oracle)，每 complex 内选 top-K。

    与 V5 的 `v5_common.topk_metrics` 口径一致（含 ≥5 成员的 complex）。
    """
    m = C.topk_metrics(np.asarray(y, float), np.asarray(p, float),
                       np.asarray(g))
    return {"efficiency_at_10": m.get("efficiency_at_10", float("nan")),
            "utility_model_mean": m.get("utility_model_mean", float("nan")),
            "utility_random_mean": m.get("utility_random_mean", float("nan")),
            "utility_oracle_mean": m.get("utility_oracle_mean", float("nan")),
            "n_complex_eff": m.get("n_complex_eff", 0),
            "p_at_10": m.get("p_at_10", float("nan")),
            "hit10_beneficial": m.get("hit10_beneficial", float("nan")),
            "hit10_ddg_m05": m.get("hit10_ddg_m05", float("nan")),
            "base_rate_best_true": m.get("base_rate_best_true", float("nan"))}


def audit_block(name: str, df: pd.DataFrame, group_col: str, truth_col: str,
                model_cols: list[str], flips: dict[str, bool],
                note: str) -> list[dict]:
    """对一组模型列算全套指标。"""
    rows = []
    d = df.copy()
    d[truth_col] = pd.to_numeric(d[truth_col], errors="coerce")
    d = d[d[truth_col].notna()]
    # 过滤：组内突变数 >= MIN_N
    cnt = d.groupby(group_col)[truth_col].transform("size")
    d = d[cnt >= MIN_N]
    n_cx = d[group_col].nunique()
    base_rate = float((d[truth_col] < 0).mean())
    emit(f"### {name}")
    emit("")
    emit(f"- 过滤后: **{len(d)} 行 / {n_cx} 个 {group_col}**"
         f"（已丢弃组内 < {MIN_N} 条的组）")
    emit(f"- 全局基准率 `ΔΔG<0` = **{base_rate:.4f}**")
    emit(f"- {note}")
    emit("")
    emit("| 模型 | Overall Pearson | Macro per-PDB Spearman | Hit@10(<0) | "
         "Base Rate | **Efficiency@10** | 反号? | n 行 | n 组 |")
    emit("|---|---:|---:|---:|---:|---:|---|---:|---:|")
    for col in model_cols:
        if col not in d.columns:
            continue
        pv = pd.to_numeric(d[col], errors="coerce")
        m = pv.notna()
        if m.sum() < 30:
            emit(f"| {col} | — | — | — | — | — | — | {int(m.sum())} | — |")
            continue
        dd = d[m]
        y = dd[truth_col].to_numpy(float)
        p_raw = pv[m].to_numpy(float)
        g = dd[group_col].astype(str).to_numpy()
        r_raw = C.pearson(y, p_raw)
        flip = bool(flips.get(col, r_raw < 0))
        p = (-p_raw) if flip else p_raw
        r = C.pearson(y, p)
        macro = C.per_structure_spearman(y, p, g)
        e = eff10(y, p, g)
        hit = float((y[np.argsort(p, kind="mergesort")[:10]] < 0).mean())
        rows.append({
            "dataset": name, "model": col, "needs_flip": flip,
            "overall_pearson": round(float(r), 4),
            "overall_pearson_noflip": round(float(r_raw), 4),
            "macro_pdb_spearman": round(float(macro), 4),
            "hit10_lt0": round(hit, 4), "base_rate_lt0": round(base_rate, 4),
            "efficiency_at_10": round(float(e["efficiency_at_10"]), 4),
            "utility_model": round(float(e["utility_model_mean"]), 4),
            "utility_random": round(float(e["utility_random_mean"]), 4),
            "utility_oracle": round(float(e["utility_oracle_mean"]), 4),
            "p_at_10": round(float(e["p_at_10"]), 4),
            "hit10_ddg_m05": round(float(e["hit10_ddg_m05"]), 4),
            "n_rows": int(m.sum()), "n_group": int(n_cx),
            "min_group_n": MIN_N, "note": note})
        emit(f"| {col} | {r:+.4f} | {macro:+.4f} | {hit:.4f} | {base_rate:.4f} | "
             f"**{e['efficiency_at_10']:+.4f}** | "
             f"{'✅' if flip else '—'} | {int(m.sum())} | {n_cx} |")
    emit("")
    return rows


def main() -> int:
    emit("# V6.0 PHASE 2.2 — 同协议 SOTA 指标重算")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 指标实现: **复用 V5 的 `v5_common`**（含已修正的 `efficiency_at_10`）")
    emit(f"- 过滤: 按 complex 分组，丢弃组内突变数 < **{MIN_N}** 的复合物")
    emit("- 约定: ΔΔG 越负越好；`Hit@10(ΔΔG<0)` = 预测前 10 中真实 ΔΔG<0 的比例")
    emit("")
    all_rows: list[dict] = []

    # ---------- A. CATH-ddG 官方 held-out CATH test ----------
    flip_map_cath = {"Flex-ddG": False, "DDGPred": False, "ddG_FoldX": False,
                     "ESM-IF": True, "proteinMPNN": False}
    if CATH_CSV.exists():
        ca = pd.read_csv(CATH_CSV, low_memory=False)
        te = ca[ca["cath_fold"].astype(str).str.strip().str.lower() == "val"].copy()
        emit("## A. CATH-ddG 官方 held-out CATH blind test（`cath_fold == val`）")
        emit("")
        emit(f"- 官方报告: **813 mutations / 53 PPIs**；本表该分区实测 "
             f"**{len(te)} 行 / {te['complex'].nunique()} complex**")
        emit("")
        all_rows += audit_block(
            "CATH-ddG/official_heldout_CATH", te, "complex", "ddG",
            ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"],
            flip_map_cath,
            "真值 = CATH 表 `ddG`（已由 V6 PHASE 2.1 用本地 SKEMPI 独立复算验证："
            "Pearson +0.9869、符号一致 94.99%）")

        # 用全表决定 flip，避免 test 小样本误判
        emit("> 说明：`needs_flip` 用 **CATH 全表**（train+test, n≈6500）判定"
             "（ESM-IF 反向、其余同向），此处套用到 test 分区。")
        emit("")
    else:
        emit("## A. CATH-ddG — `DATA_MISSING`")
        emit("")

    # ---------- B. DSSA-PPI RBM 案例集 ----------
    if RBM_CSV.exists():
        rb = pd.read_csv(RBM_CSV, low_memory=False)
        emit("## B. DSSA-PPI RBM 案例集（deep mutational scanning 风格）")
        emit("")
        emit(f"- 行数 **{len(rb)}**；`pos` 唯一值: "
             f"**{sorted(rb['pos'].dropna().unique().tolist())[:12]}**")
        emit(f"- ⚠️ 真值列是 `EXP_fitness`（实验**适应度**，非 ΔΔG）；"
             f"**符号语义未由仓库说明** →")
        emit(f"  本数据集**不进入 `Hit@10(ΔΔG<0)` 的绝对口径**，只用于 "
             f"Pearson / Spearman / within-position 排序对比。")
        emit("")
        # 过滤：该数据集只有一个 pos → 视为单一 target；为满足「组内≥5」按 pos 分组
        # 另做一个按 pos 分箱的版本（把 1702 条按位置连续分块）以产生可比的分组排序
        rb2 = rb.copy()
        rb2["pos"] = pd.to_numeric(rb2["pos"], errors="coerce")
        rb2["group_bin"] = (rb2["pos"].rank(method="dense") // 200).astype(int)
        cols = ["DSSA-PPI", "PPIFormer", "ESM-2-FFN", "MuTon", "MpbPPI",
                "ProBASS", "SAAMBE-SEQ", "FoldX"]
        flips = {"DSSA-PPI": False, "PPIFormer": True, "ESM-2-FFN": True,
                 "MuTon": False, "MpbPPI": False, "ProBASS": False,
                 "SAAMBE-SEQ": True, "FoldX": True}
        all_rows += audit_block(
            "DSSA-PPI/RBM_case_study", rb2, "group_bin", "EXP_fitness", cols,
            flips,
            "真值 = `EXP_fitness`（适应度，非 ΔΔG）；分组 = `pos` 连续分箱"
            "（每箱约 200 个位置序号）—— 这是**代理分组**，"
            "不等价于多个独立复合物，须谨慎解读")
    else:
        emit("## B. DSSA-PPI RBM — `DATA_MISSING`")
        emit("")

    # ---------- C. V5 baseline 对照 ----------
    emit("## C. V5 baseline 对照（同协议同指标实现）")
    emit("")
    p4 = C.V5 / "p4_conformal" / "p4_predictions.npz"
    if p4.exists():
        z = np.load(p4, allow_pickle=True)
        y = z["test_y"]
        g = np.array([str(x) for x in z["test_group"]])
        n_cx = len(set(g))
        base_rate = float((y < 0).mean())
        emit(f"- 数据: **{len(y)} 行 / {n_cx} 复合物**（V5 PPB-Affinity development benchmark）")
        emit(f"- 全局基准率 `ΔΔG<0` = **{base_rate:.4f}**")
        emit("")
        emit("| 模型 | Overall Pearson | Macro per-PDB Spearman | Hit@10(<0) | "
             "Base Rate | **Efficiency@10** | n 行 | n 组 |")
        emit("|---|---:|---:|---:|---:|---:|---:|---:|")
        for key in ("plain_s42", "plain_s43", "plain_s44"):
            if key not in z.files:
                continue
            p = z[key]
            r = C.pearson(y, p)
            macro = C.per_structure_spearman(y, p, g)
            e = eff10(y, p, g)
            hit = float((y[np.argsort(p, kind="mergesort")[:10]] < 0).mean())
            emit(f"| V5 baseline ({key}) | {r:+.4f} | {macro:+.4f} | {hit:.4f} | "
                 f"{base_rate:.4f} | **{e['efficiency_at_10']:+.4f}** | {len(y)} | "
                 f"{n_cx} |")
            all_rows.append({
                "dataset": "V5/PPB-Affinity_development", "model": f"V5_baseline({key})",
                "needs_flip": False, "overall_pearson": round(float(r), 4),
                "overall_pearson_noflip": round(float(r), 4),
                "macro_pdb_spearman": round(float(macro), 4),
                "hit10_lt0": round(hit, 4), "base_rate_lt0": round(base_rate, 4),
                "efficiency_at_10": round(float(e["efficiency_at_10"]), 4),
                "utility_model": round(float(e["utility_model_mean"]), 4),
                "utility_random": round(float(e["utility_random_mean"]), 4),
                "utility_oracle": round(float(e["utility_oracle_mean"]), 4),
                "p_at_10": round(float(e["p_at_10"]), 4),
                "hit10_ddg_m05": round(float(e["hit10_ddg_m05"]), 4),
                "n_rows": int(len(y)), "n_group": int(n_cx), "min_group_n": MIN_N,
                "note": "V5 baseline, PPB-Affinity triple-isolation split"})
        emit("")
    else:
        emit(f"- ⚠️ `{p4}` 不存在 → 无法对照")
        emit("")

    res = pd.DataFrame(all_rows)
    res.to_csv(OUT / "sota_audit_metrics.csv", index=False, encoding="utf-8")

    # ---------- 判定 ----------
    emit("## D. 假设判定（V6 §6 成功标准）")
    emit("")
    ext = res[~res["dataset"].str.startswith("V5/")]
    hall = ext[(ext["overall_pearson"] > 0.4)
               & (ext["macro_pdb_spearman"] < 0.15)]
    emit("| 判据 | 阈值 | 满足的外部模型数 |")
    emit("|---|---|---:|")
    emit(f"| Overall Pearson > 0.4 | >0.4 | "
         f"**{int((ext['overall_pearson'] > 0.4).sum())}** / {len(ext)} |")
    emit(f"| Macro per-PDB Spearman < 0.15 | <0.15 | "
         f"**{int((ext['macro_pdb_spearman'] < 0.15).sum())}** / {len(ext)} |")
    emit(f"| **两者同时成立（幻觉模式）** | — | **{len(hall)}** / {len(ext)} |")
    emit(f"| Macro per-PDB Spearman > 0.3（反驳） | >0.3 | "
         f"**{int((ext['macro_pdb_spearman'] > 0.3).sum())}** / {len(ext)} |")
    emit("")
    if len(hall) >= 1:
        emit(f"→ **假设 SUPPORTED**：至少 1 个外部模型同时满足"
             f"「全局高（>0.4）+ 局部低（<0.15）」。")
        emit("")
        emit("满足该模式的模型：")
        emit("")
        emit("| 模型 | Overall Pearson | Macro per-PDB Spearman | Hit@10 | Base Rate | Eff@10 |")
        emit("|---|---:|---:|---:|---:|---:|")
        for _, r in hall.iterrows():
            emit(f"| {r['model']} ({r['dataset']}) | {r['overall_pearson']:+.4f} | "
                 f"{r['macro_pdb_spearman']:+.4f} | {r['hit10_lt0']:.4f} | "
                 f"{r['base_rate_lt0']:.4f} | {r['efficiency_at_10']:+.4f} |")
    elif (ext["macro_pdb_spearman"] > 0.3).all() and len(ext):
        emit("→ **假设 REFUTED**：所有外部模型 Macro per-PDB Spearman > 0.3。")
    else:
        emit("→ **假设 INCONCLUSIVE**：外部模型数量或样本量不足。")
    emit("")
    emit("### 样本量警示（V6 §3 问题 3）")
    emit("")
    for ds, gg in res.groupby("dataset"):
        n = int(gg["n_group"].iloc[0])
        flag = "⚠️ **样本量不足，结论需谨慎**" if n < 30 else "✅ 样本量可用"
        emit(f"- `{ds}`: n_group = **{n}** → {flag}")
    emit("")

    (OUT / "sota_audit_metrics.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT / 'sota_audit_metrics.csv'}")
    print(f"写出 {OUT / 'sota_audit_metrics.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
