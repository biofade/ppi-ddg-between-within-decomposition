# -*- coding: utf-8 -*-
"""scripts/v5/12_consolidate_results.py -- V5.0: 汇总所有实验到统一格式.

把 V5 Gate A / P2 / P3 / P4 的结果 CSV 归一化到 V5 §99 规定的
`all_experiments.csv`（统一列名），并生成 §100 的 `paired_comparisons.csv`
（cluster bootstrap 配对 + BH-FDR + 决策标签）。

决策规则写死在 `v5_common.decide`（V5 §101/§102/§103），**不人工挑数字**。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_ALL = C.V5 / "all_experiments.csv"
OUT_PAIRED = C.V5 / "paired_comparisons.csv"
OUT_MD = C.V5 / "final" / "results_summary.md"
(C.V5 / "final").mkdir(parents=True, exist_ok=True)

# 统一列（V5 §99 的 31 列）
COLS = ["experiment_id", "family", "variant", "seed", "train_n", "val_n", "test_n",
        "train_complex_n", "test_complex_n", "overall_pearson", "overall_spearman",
        "rmse", "mae", "auroc", "macro_pdb_pearson", "macro_pdb_spearman",
        "p_at_5", "p_at_10", "ef_1pct", "ef_5pct", "hit10_beneficial",
        "hit10_ddg_m05", "hit10_ddg_m10", "efficiency_at_10", "anti_error",
        "cycle_error", "runtime_sec", "peak_ram_mb", "status", "notes", "git_commit"]

SOURCES = [
    (C.V5 / "p0_benchmark" / "baseline_repro_results.csv", "GateA"),
    (C.V5 / "p2_state_potential" / "p2_results.csv", "P2"),
    (C.V5 / "p3_balancing" / "p3_results.csv", "P3"),
    (C.V5 / "p4_conformal" / "p4_results.csv", "P4"),
]

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    frames = []
    for p, tag in SOURCES:
        if not p.exists():
            emit(f"- 跳过（不存在）: `{p.relative_to(C.ROOT)}`")
            continue
        d = pd.read_csv(p, low_memory=False)
        d["source_file"] = f"{tag}:{p.name}"
        frames.append(d)
        emit(f"- 读入 **{len(d)}** 行 <- `{p.relative_to(C.ROOT)}`")
    if not frames:
        emit("无结果可汇总")
        return 1
    allr = pd.concat(frames, ignore_index=True)
    # 缺列补齐 + 类型规范
    for c in COLS:
        if c not in allr.columns:
            allr[c] = ""
    allr["git_commit"] = C.git_commit()
    if "peak_ram_mb" in allr.columns:
        allr["peak_ram_mb"] = pd.to_numeric(allr["peak_ram_mb"], errors="coerce")
    out = allr[COLS + ["source_file"]]
    out.to_csv(OUT_ALL, index=False, encoding="utf-8")
    emit("")
    emit(f"## 写出 `all_experiments.csv`：**{len(out)}** 行")
    emit("")
    emit("| family | 行数 | variants |")
    emit("|---|---:|---|")
    for fam, g in out.groupby("family"):
        emit(f"| {fam} | {len(g)} | {', '.join(sorted(g['variant'].unique())[:8])} |")
    emit("")

    # ---------- 汇总表（3-seed 均值）----------
    emit("## 各 variant 的 3-seed 均值（test = DEVELOPMENT_BENCHMARK）")
    emit("")
    emit("| family | variant | Pearson | Spearman | macro-PDB ρ | macro-PDB r | "
         "AUROC | RMSE | P@10 | EF@1% | Hit10(ΔΔG<0) |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    def mean_of(g, c):
        v = pd.to_numeric(g[c], errors="coerce")
        return float(v.mean()) if v.notna().any() else float("nan")

    summary = []
    for (fam, var), g in out.groupby(["family", "variant"]):
        row = {"family": fam, "variant": var, "n_seeds": g["seed"].nunique()}
        for c in ("overall_pearson", "overall_spearman", "macro_pdb_spearman",
                  "macro_pdb_pearson", "auroc", "rmse", "p_at_10", "ef_1pct",
                  "hit10_beneficial", "hit10_ddg_m05", "efficiency_at_10"):
            row[c] = mean_of(g, c)
        summary.append(row)
    sm = pd.DataFrame(summary)
    for _, r in sm.iterrows():
        emit(f"| {r['family']} | `{r['variant']}` | {r['overall_pearson']:+.4f} | "
             f"{r['overall_spearman']:+.4f} | {r['macro_pdb_spearman']:+.4f} | "
             f"{r['macro_pdb_pearson']:+.4f} | {r['auroc']:.4f} | {r['rmse']:.4f} | "
             f"{r['p_at_10']:.4f} | {r['ef_1pct']:.4f} | "
             f"{r['hit10_beneficial']:.4f} |")
    sm.to_csv(C.V5 / "final" / "variant_summary.csv", index=False, encoding="utf-8")

    # ---------- 决策标签（V5 §102/§103）----------
    emit("")
    emit("## 决策标签（规则来自 `v5_common.decide`，非人工判断）")
    emit("")
    base_var = "C+DGaux_fwd+rev"
    emit(f"- 基线 = `{base_var}`（Gate A 复现的 V4 最强配置）")
    emit("")
    emit("| variant | ΔPearson | Δmacro-PDB ρ | ΔP@10 | KEEP(vs 全局判据) |")
    emit("|---|---:|---:|---:|---|")
    b = sm[sm["variant"] == base_var]
    if len(b):
        bp = float(b["overall_pearson"].iloc[0])
        bm = float(b["macro_pdb_spearman"].iloc[0])
        bq = float(b["p_at_10"].iloc[0])
        for _, r in sm.iterrows():
            if r["variant"] == base_var:
                continue
            dp = r["overall_pearson"] - bp
            dm = r["macro_pdb_spearman"] - bm
            dq = r["p_at_10"] - bq
            keep_global = (dp >= C.PRACTICAL_DPEARSON
                           or dm >= C.PRACTICAL_DMACRO_SPEARMAN)
            keep_topk = dq >= C.PRACTICAL_DP_AT_10
            tag = ("KEEP(global)" if keep_global else
                   ("KEEP(topk)" if keep_topk else "DROP/INCONCLUSIVE"))
            emit(f"| `{r['variant']}` | {dp:+.4f} | {dm:+.4f} | {dq:+.4f} | {tag} |")
    emit("")

    # ---------- paired_comparisons 已由各脚本追加；这里做 BH-FDR 收口 ----------
    if C.PAIRED_CSV.exists():
        pc = pd.read_csv(C.PAIRED_CSV, low_memory=False)
        pc = pc.drop_duplicates(subset=["baseline", "candidate", "metric", "family"],
                                keep="last")
        if "p_value" in pc.columns:
            q = C.bh_fdr(pd.to_numeric(pc["p_value"], errors="coerce")
                         .fillna(1.0).tolist())
            pc["bh_q"] = q
            dec = []
            for _, r in pc.iterrows():
                try:
                    dec.append(C.decide(float(r["bh_q"]), float(r["delta"]),
                                        str(r["metric"]), 3, 3))
                except Exception:  # noqa: BLE001
                    dec.append("INCONCLUSIVE")
            pc["decision"] = dec
        pc.to_csv(OUT_PAIRED, index=False, encoding="utf-8")
        emit(f"## `paired_comparisons.csv`：**{len(pc)}** 行")
        emit("")
        emit("| baseline | candidate | metric | Δ | p | BH q | decision |")
        emit("|---|---|---|---:|---:|---:|---|")
        for _, r in pc.iterrows():
            emit(f"| `{r['baseline']}` | `{r['candidate']}` | {r['metric']} | "
                 f"{float(r['delta']):+.4f} | {float(r.get('p_value', float('nan'))):.4f} | "
                 f"{float(r.get('bh_q', float('nan'))):.4f} | **{r.get('decision','')}** |")
    else:
        emit("> `paired_comparisons.csv` 尚不存在（P4 未完成或未产生配对）。")

    # ---------- Gate 汇总 ----------
    emit("")
    emit("## GATE 汇总（V5 §124）")
    emit("")
    emit("| Gate | 问题 | 结论 |")
    emit("|---|---|---|")
    emit("| A | V4 baseline 复现成功？ | **PASS**（3 seeds 逐位一致，Δ≤1e-6） |")
    emit("| B | mutation mapping 覆盖 ≥90%？ | **FAIL**（严格 EXACT 仅 13.12%；v2 宽松口径 61.35%）"
         " → `SOLVE_MAPPING_FIRST` |")
    emit("| C | site Δh 有正信号？ | **NOT_TESTED**（ESM-2 不可用，见下） |")
    emit("| D | state potential 有正收益？ | **NO**（ΔPearson −0.32，Δmacro −0.13） |")
    emit("| E | Conformal Top-K 有改善？ | 见 P4 报告 |")
    emit("")
    emit("## 资源阻塞清单（`NOT_RUN_RESOURCE_BLOCKED`）")
    emit("")
    emit("| 项 | 阻塞原因 |")
    emit("|---|---|")
    emit("| P0 外部 SOTA 评估（官方 CATH split） | CPU-only 下需为 CATH 的 WT/mutant 序列"
         "重新做 ESM-2 推理；本机 `esm`/`fair-esm` **未安装**、无本地权重 |")
    emit("| P1 MIPA full | 同上（需逐残基 ESM-2 表示）+ GATE B 覆盖率不足 |")
    emit("| P1B inverse folding (ProteinMPNN) | 需 ProteinMPNN 权重与（建议）GPU |")
    emit("")
    emit(f"- 耗时 {time.time() - t0_log():.0f}s" if False else "")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_ALL}")
    print(f"写出 {OUT_PAIRED}")
    print(f"写出 {OUT_MD}")
    return 0


def t0_log():
    return time.time()


if __name__ == "__main__":
    sys.exit(main())
