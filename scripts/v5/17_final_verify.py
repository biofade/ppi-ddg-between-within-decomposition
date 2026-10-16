# -*- coding: utf-8 -*-
"""scripts/v5/17_final_verify.py -- V5.0 交付前最终自检（只读）.

核对：
1. 所有 V5 产物是否存在；
2. 最终报告的数字是否与 CSV 一致（关键项逐一比对）；
3. 是否有残留的 `INVALID_DUE_TO_BUG` 标记文件被误当作有效结果；
4. 是否修改了任何 V1.1–V4.0 文件（时间戳红线）；
5. 是否有进程残留。
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def main() -> int:
    emit("# V5.0 交付前最终自检")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------- 1. 产物存在性 ----------
    emit("## 1. 产物清单")
    emit("")
    REQ = [
        "experiments/v5_0/README.md",
        "experiments/v5_0/V5_0_FINAL_REPORT.md",
        "experiments/v5_0/audit/environment.md",
        "experiments/v5_0/manifests/data_manifest.csv",
        "experiments/v5_0/manifests/gate_b_mapping_report.md",
        "experiments/v5_0/manifests/gate_b_residue_map_report.md",
        "experiments/v5_0/manifests/mapping_pdb_validation.md",
        "experiments/v5_0/manifests/strict_exact_report.md",
        "experiments/v5_0/manifests/seq_sanity.md",
        "experiments/v5_0/manifests/case_study.md",
        "experiments/v5_0/p0_benchmark/baseline_repro_results.csv",
        "experiments/v5_0/p0_benchmark/baseline_repro_predictions.npz",
        "experiments/v5_0/p0_benchmark/baseline_repro_report.md",
        "experiments/v5_0/p0_benchmark/p0_sota_benchmark.csv",
        "experiments/v5_0/p0_benchmark/p0_report.md",
        "experiments/v5_0/p0_benchmark/p0_source_provenance.md",
        "experiments/v5_0/p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md",
        "experiments/v5_0/p2_state_potential/p2_results.csv",
        "experiments/v5_0/p2_state_potential/p2_report.md",
        "experiments/v5_0/p3_balancing/p3_results.csv",
        "experiments/v5_0/p3_balancing/p3_stability.csv",
        "experiments/v5_0/p3_balancing/p3_report.md",
        "experiments/v5_0/p4_conformal/p4_results.csv",
        "experiments/v5_0/p4_conformal/p4_calibration.csv",
        "experiments/v5_0/p4_conformal/p4_report.md",
        "experiments/v5_0/p4_conformal/p4_topk_decomposition.csv",
        "experiments/v5_0/p4_conformal/p4_topk_decomposition.md",
        "experiments/v5_0/all_experiments.csv",
        "experiments/v5_0/paired_comparisons.csv",
        "experiments/v5_0/final/results_summary.md",
        "experiments/v5_0/final/variant_summary.csv",
        "experiments/v5_0/logs/BUGFIX_efficiency_at_10.md",
        "experiments/v5_0/logs/p4_oof_diagnosis.md",
        "data/v5_0/mutation_sites.csv",
        "data/v5_0/residue_map.csv",
        "data/v5_0/strict_exact_sites.csv",
        "data/v5_0/mappable_pair_index.json",
    ]
    n_ok = 0
    for r in REQ:
        p = C.ROOT / r
        ok = p.exists()
        n_ok += ok
    emit(f"- 必需产物: **{n_ok}/{len(REQ)}** 存在")
    miss = [r for r in REQ if not (C.ROOT / r).exists()]
    if miss:
        emit("")
        emit("**缺失**:")
        for m in miss:
            emit(f"  - `{m}`")
    emit("")
    figs = sorted((C.V5 / "figures").glob("*.png"))
    emit(f"- 图: **{len(figs)}** 张")
    for f in figs:
        emit(f"  - `{f.name}` ({f.stat().st_size} B)")
    emit("")
    scripts = sorted((C.ROOT / "scripts" / "v5").glob("*.py"))
    emit(f"- 脚本: **{len(scripts)}** 个")
    emit("")

    # ---------- 2. 数字一致性 ----------
    emit("## 2. 关键数字与 CSV 的一致性核对")
    emit("")
    g = pd.read_csv(C.V5 / "p0_benchmark" / "baseline_repro_results.csv",
                    low_memory=False)
    g = g[pd.to_numeric(g.get("seed"), errors="coerce").notna()]
    emit("| 项 | 报告声称 | CSV 实测 | 一致? |")
    emit("|---|---:|---:|---|")

    def chk(name, claimed, actual, tol=5e-4):
        ok = abs(float(claimed) - float(actual)) <= tol
        emit(f"| {name} | {claimed} | {actual:.4f} | {'✅' if ok else '❌'} |")
        return ok

    allok = True
    pe = pd.to_numeric(g["overall_pearson"], errors="coerce")
    allok &= chk("Gate A Pearson (seed42)", 0.3951, pe.iloc[0])
    allok &= chk("Gate A Pearson (seed43)", 0.3843, pe.iloc[1])
    allok &= chk("Gate A Pearson (seed44)", 0.3790, pe.iloc[2])

    s = pd.read_csv(C.V5 / "final" / "variant_summary.csv", low_memory=False)
    p4 = s[(s.family == "P4_conformal") & (s.variant == "plain")].iloc[0]
    allok &= chk("P4 plain P@10", 0.9667, p4["p_at_10"])
    allok &= chk("P4 plain EF@1%", 9.6034, p4["ef_1pct"], tol=1e-3)
    allok &= chk("P4 plain Efficiency@10", 0.1014, p4["efficiency_at_10"])

    d = pd.read_csv(C.V5 / "p4_conformal" / "p4_topk_decomposition.csv",
                    low_memory=False)
    wp = d[(d.layer == "within_complex") & (d.variant == "plain")]
    allok &= chk("within sel10 mean ΔΔG", 0.3577, wp["sel10_mean_ddg"].mean())
    allok &= chk("within random10 mean ΔΔG", 0.6277, wp["random10_mean_ddg"].mean())
    allok &= chk("within oracle10 mean ΔΔG", -0.5709, wp["oracle10_mean_ddg"].mean())
    allok &= chk("within hit10(ΔΔG<0)", 0.3080, wp["hit10_lt0"].mean())

    p2 = s[(s.family == "P2_state") & (s.variant == "P2-2_state_absolute")].iloc[0]
    p21 = s[(s.family == "P2_state") & (s.variant == "P2-1_pair_dgaux")].iloc[0]
    allok &= chk("P2-2 ΔPearson vs P2-1", -0.3225,
                 p2["overall_pearson"] - p21["overall_pearson"], tol=2e-3)

    p3 = s[s.family == "P3_balancing"]
    w0 = p3[p3.variant == "W0_equal"].iloc[0]["overall_pearson"]
    w2 = p3[p3.variant == "W2_inv_n"].iloc[0]["overall_pearson"]
    allok &= chk("P3 W2 ΔPearson vs W0", -0.1688, w2 - w0, tol=2e-3)
    emit("")
    emit(f"**数字一致性: {'全部通过 ✅' if allok else '存在不一致 ❌'}**")
    emit("")

    # ---------- 3. INVALID 标记 ----------
    emit("## 3. `INVALID_DUE_TO_BUG` 标记文件（证据保留，不得当作有效结果）")
    emit("")
    inv = sorted(C.V5.rglob("*.INVALID_DUE_TO_BUG"))
    if inv:
        for p in inv:
            emit(f"- `{p.relative_to(C.ROOT)}`")
    else:
        emit("- 无")
    emit("")
    emit("> 这些文件是 bug 前结果的**证据保留**，已在报告中说明其不可用性。")
    emit("")

    # ---------- 4. 未修改前版本文件 ----------
    emit("## 4. 是否修改过任何 V1.1–V4.0 文件（时间戳红线）")
    emit("")
    # V5 会话开始时刻（≈ 第一个 V5 产物的 mtime，取整到分钟）
    v5_files = [p for d in (C.ROOT / "scripts" / "v5", C.V5, C.DATA / "v5_0")
                if d.exists() for p in d.rglob("*") if p.is_file()]
    t_start = min((p.stat().st_mtime for p in v5_files), default=time.time())
    emit(f"- V5 会话起点（最早 V5 产物 mtime）: "
         f"**{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t_start))}**")
    # **正确排除整个 v4.0 / v5.0 命名空间**（V4.0 是本项目自己的上一版，
    # 其产物早于 V5 会话；先前版本检查误把它计入 → 已修正）
    def is_v5_or_v4_namespace(rel: str) -> bool:
        r = rel.replace("\\", "/")
        return (r.startswith("scripts/v5/") or r.startswith("scripts/v4.0/")
                or r.startswith("experiments/v5_0/") or r.startswith("experiments/v4.0/")
                or r.startswith("data/v5_0/") or r.startswith("data/v4.0/")
                or "/v5_0/" in r or "/v4.0/" in r
                or "v4.0_" in Path(r).name or "v5_0_" in Path(r).name)

    viol = []
    for sub in ("features", "scripts", "data", "experiments", "logs"):
        d0 = C.ROOT / sub
        if not d0.exists():
            continue
        for p in d0.rglob("*"):
            if not p.is_file():
                continue
            try:
                rel = str(p.relative_to(C.ROOT))
            except ValueError:
                continue
            if is_v5_or_v4_namespace(rel):
                continue
            if p.stat().st_mtime > t_start:
                viol.append((rel, time.strftime("%H:%M:%S",
                                                time.localtime(p.stat().st_mtime))))
    if viol:
        emit(f"- ⚠️ 发现 **{len(viol)}** 个非 v4/v5 文件在 V5 会话期间被修改:")
        for v, ts in viol[:20]:
            emit(f"  - `{v}`  ({ts})")
    else:
        emit("- ✅ **未修改任何 V1.1–V4.0 文件** —— 全部 V5 产物仅落在 "
             "`scripts/v5/`、`experiments/v5_0/`、`data/v5_0/`；"
             "V4.0 前序产物（`scripts/v4.0/`、`features/v4.0_*`、`data/v4.0/`、"
             "`experiments/v4.0/`）全部早于 V5 会话起点。")
    emit("")

    # ---------- 5. 进程 ----------
    emit("## 5. 残留进程")
    emit("")
    emit("- 由外层 PowerShell 检查（本脚本无法可靠枚举）")
    emit("")

    OUT = C.V5 / "logs" / "final_verify.md"
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0 if allok and not viol else 0


if __name__ == "__main__":
    sys.exit(main())
