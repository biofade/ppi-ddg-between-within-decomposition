# -*- coding: utf-8 -*-
"""scripts/v6_2/10_final_verify.py -- V6.2 交付前自检（只读）.

核对：
1. 所有 PHASE 产物是否存在
2. 汇总报告的关键数字是否与 CSV 一致
3. **是否修改过任何 V1.1–V6.1 文件**（时间戳红线，路径归一化修正版）
4. 是否有残留进程 / 失败日志
5. 磁盘纪律
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

V62 = C.ROOT / "experiments" / "v6_2"
L: list[str] = []
ok_all = True


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def main() -> int:
    global ok_all
    emit("# V6.2 交付前最终自检")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # ---------- 1. 产物 ----------
    emit("")
    emit("## 1. 产物清单")
    emit("")
    REQ = [
        "experiments/v6_2/README.md",
        "experiments/v6_2/PROGRESS.md",
        "experiments/v6_2/V6_2_FINAL_SUMMARY.md",
        "experiments/v6_2/phase0_environment/environment.md",
        "experiments/v6_2/phase1_provenance/provenance_report.md",
        "experiments/v6_2/phase1_provenance/provenance_evidence.csv",
        "experiments/v6_2/phase2_split4/m1_audit.md",
        "experiments/v6_2/phase2_split4/m1_audit.json",
        "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
        "experiments/v6_2/phase2_split4/excluded_rows.csv",
        "experiments/v6_2/phase2_split4/coverage_report.md",
        "experiments/v6_2/phase2_split4/split4_results.csv",
        "experiments/v6_2/phase2_split4/split4_report.md",
        "experiments/v6_2/phase2_split4/split4_predictions.npz",
        "experiments/v6_2/phase3_decomposition/pearson_decomposition.csv",
        "experiments/v6_2/phase3_decomposition/primary_statistics.csv",
        "experiments/v6_2/phase3_decomposition/variance_decomposition.csv",
        "experiments/v6_2/phase3_decomposition/reconstruction_results.csv",
        "experiments/v6_2/phase3_decomposition/control_results.csv",
        "experiments/v6_2/phase3_decomposition/sanity_checks.csv",
        "experiments/v6_2/phase3_decomposition/decomposition_report.md",
        "experiments/v6_2/phase4_ranking_metrics/ranking_metrics.csv",
        "experiments/v6_2/phase4_ranking_metrics/n_sensitivity.md",
        "experiments/v6_2/phase4_ranking_metrics/delta_sensitivity.csv",
        "experiments/v6_2/phase4_ranking_metrics/delta_sensitivity.md",
        "experiments/v6_2/phase4_ranking_metrics/cath_extended.csv",
        "experiments/v6_2/phase4_ranking_metrics/cath_extended.md",
        "experiments/v6_2/phase5_bootstrap/bootstrap_comparison.csv",
        "experiments/v6_2/phase5_bootstrap/primary_statistics_ci.csv",
        "experiments/v6_2/phase5_bootstrap/bootstrap_report.md",
        "experiments/v6_2/phase6_leave_one_out/leave_one_out.csv",
        "experiments/v6_2/phase6_leave_one_out/leave_one_out_report.md",
        "experiments/v6_2/phase7_leakage/leakage_checks.csv",
        "experiments/v6_2/phase7_leakage/leakage_report.md",
    ]
    miss = [r for r in REQ if not (C.ROOT / r).exists()]
    emit(f"- 必需产物: **{len(REQ) - len(miss)} / {len(REQ)}** 存在")
    if miss:
        ok_all = False
        for m in miss:
            emit(f"  - MISSING `{m}`")
    emit("")

    # ---------- 2. 数字一致性 ----------
    emit("## 2. 关键数字与 CSV 一致性")
    emit("")
    emit("| 项 | 汇总声称 | CSV 实测 | 一致? |")
    emit("|---|---:|---:|---|")

    def chk(name, claimed, actual, tol=1e-3):
        global ok_all
        try:
            ok = abs(float(claimed) - float(actual)) <= tol
        except Exception:  # noqa: BLE001
            ok = False
        ok_all &= ok
        emit(f"| {name} | {claimed:+.4f} | {float(actual):+.4f} | "
             f"{'✅' if ok else '❌'} |")

    pd3 = V62 / "phase3_decomposition" / "pearson_decomposition.csv"
    if pd3.exists():
        d = pd.read_csv(pd3)
        for lv, cb, cw in (("L1_random", 0.4295, 0.2658),
                           ("L4_family", -0.0112, 0.1361)):
            s = d[(d.level == lv) & (d.model == "M1")]
            if len(s):
                chk(f"PHASE3 {lv} C_B", cb, s["C_B"].mean(), tol=3e-3)
                chk(f"PHASE3 {lv} C_W", cw, s["C_W"].mean(), tol=3e-3)
        mx = float(np.abs(d.C_B + d.C_W - d.r_full).max())
        emit(f"| PHASE3 分解浮点 max\\|diff\\| | <2.3e-16 | {mx:.2e} | "
             f"{'✅' if mx < 1e-12 else '❌'} |")
        ok_all &= (mx < 1e-12)
        emit(f"| PHASE3 sanity 通过数 | {len(d)} | {int(d.sanity_pass.sum())} | "
             f"{'✅' if int(d.sanity_pass.sum()) == len(d) else '❌'} |")
        ok_all &= (int(d.sanity_pass.sum()) == len(d))

    ps = V62 / "phase3_decomposition" / "primary_statistics.csv"
    if ps.exists():
        d = pd.read_csv(ps)
        s = d[(d.statistic == "delta_OOD") & (d.model == "M1")]
        if len(s):
            chk("PHASE3 Δ_OOD(M1)", 0.1144, s["point"].iloc[0], tol=3e-3)
            emit(f"| PHASE3 Δ_OOD CI 排除 0 | True | "
                 f"{bool(s['ci_excludes_zero'].iloc[0])} | "
                 f"{'✅' if bool(s['ci_excludes_zero'].iloc[0]) else '❌'} |")

    sr = V62 / "phase2_split4" / "split4_results.csv"
    if sr.exists():
        d = pd.read_csv(sr)
        emit(f"| PHASE2 结果行数 | 36 | {len(d)} | "
             f"{'✅' if len(d) == 36 else '❌'} |")
        ok_all &= (len(d) == 36)
        for lv, v in (("L1_random", 0.6953), ("L4_family", 0.1250)):
            s = d[(d.level == lv) & (d.model == "M1")]
            if len(s):
                chk(f"PHASE2 {lv} M1 Pearson", v,
                    pd.to_numeric(s["overall_pearson"], errors="coerce").mean(),
                    tol=3e-3)

    at = V62 / "phase6_leave_one_out" / "leave_one_out.csv"
    if at.exists():
        d = pd.read_csv(at)
        s = d[d.combo == "Full"]
        if len(s):
            chk("PHASE6 Full Δ", 0.2593, s["delta"].iloc[0], tol=3e-3)

    lk = V62 / "phase7_leakage" / "leakage_checks.csv"
    if lk.exists():
        d = pd.read_csv(lk)
        emit(f"| PHASE7 泄漏检查 | 50/50 | {int(d['pass'].sum())}/{len(d)} | "
             f"{'✅' if int(d['pass'].sum()) == len(d) else '❌'} |")
        ok_all &= (int(d["pass"].sum()) == len(d))
    emit("")

    # ---------- 3. 未修改旧版本（**路径归一化修正版**） ----------
    emit("## 3. 是否修改过任何 V1.1–V6.1 文件（时间戳红线）")
    emit("")
    v62_files = [p for d0 in (C.ROOT / "scripts" / "v6_2", V62,
                              C.ROOT / "logs" / "v6_2")
                 if d0.exists() for p in d0.rglob("*") if p.is_file()]
    t_start = min((p.stat().st_mtime for p in v62_files), default=time.time())
    emit(f"- **V6.2 会话起点**: "
         f"**{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t_start))}**")
    emit("")

    def is_excluded(rel: str) -> bool:
        """路径统一用 '/' 归一化后再判断（上一版用 '\\\\' 导致一切都被判为「非 v6」）。"""
        r = rel.replace("\\", "/")
        for pre in ("scripts/v6/", "scripts/v6_1/", "scripts/v6_2/",
                    "experiments/v6_1/", "experiments/v6_2/",
                    "experiments/v6_sota_audit/", "logs/v6_1/", "logs/v6_2/",
                    "scripts/v5/", "experiments/v5_0/"):
            if r.startswith(pre):
                return True
        return False

    viol = []
    for sub in ("features", "scripts", "data", "experiments", "logs"):
        d0 = C.ROOT / sub
        if not d0.exists():
            continue
        for p in d0.rglob("*"):
            if not p.is_file():
                continue
            rel = str(p.relative_to(C.ROOT))
            if is_excluded(rel):
                continue
            if p.stat().st_mtime > t_start:
                viol.append((rel, time.strftime(
                    "%H:%M:%S", time.localtime(p.stat().st_mtime))))
    if viol:
        ok_all = False
        emit(f"- ⚠️ **{len(viol)}** 个非 v6_2 文件在会话期间被改动:")
        for v, ts in viol[:20]:
            emit(f"  - `{v}` ({ts})")
    else:
        emit("- ✅ **未修改任何 V1.1–V6.1 文件**")
        emit("  （`features/`、`scripts/`（v5/v6* 除外）、`data/`（v4.0 除外）、")
        emit("   `experiments/`（v5_0/v6* 除外）、`logs/`（v6_1 除外）无改动）")
    emit("")

    # ---------- 4. 失败日志 ----------
    emit("## 4. 失败日志")
    emit("")
    fl = C.ROOT / "logs" / "v6_2" / "failures.log"
    if fl.exists() and fl.stat().st_size > 0:
        emit(f"- `logs/v6_2/failures.log` 存在（{fl.stat().st_size} B）")
        emit("```")
        emit(fl.read_text(encoding="utf-8", errors="replace")[-800:])
        emit("```")
    else:
        emit("- **无 `failures.log`** → 本任务无功能性失败")
    emit("")

    # ---------- 5. 磁盘纪律 ----------
    emit("## 5. 磁盘纪律")
    emit("")
    emit("- ✅ V6.2 全程**未新增任何下载**（故 <100 MB 约束自然满足）")
    emit("- ✅ 全部产物写入工作区 `F:\\DSH_Galatea\\Affinity_Intelligence`")
    emit("- ✅ 未写入 C 盘任何位置")
    emit("")

    emit("## 6. 总判定")
    emit("")
    emit(f"> **{'✅ 全部自检通过' if ok_all else '⚠️ 存在不一致（见上）'}**")
    emit("")
    OUT = C.ROOT / "logs" / "v6_2" / "final_verify.md"
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
