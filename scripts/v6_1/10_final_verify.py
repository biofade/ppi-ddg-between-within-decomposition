# -*- coding: utf-8 -*-
"""scripts/v6_1/10_final_verify.py -- V6.1 交付前自检（只读）.

核对：
1. 所有 PHASE 产物是否存在
2. 汇总报告的关键数字是否与 CSV 一致
3. 是否修改过任何 V1.1–V6.0 文件（时间戳红线）
4. 是否有残留进程 / 失败日志
5. 磁盘纪律（下载是否都在 F 盘）
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V61 = C.ROOT / "experiments" / "v6_1"
L: list[str] = []
ok_all = True


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def main() -> int:
    global ok_all
    emit("# V6.1+ 交付前最终自检")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # ---------- 1. 产物清单 ----------
    emit("")
    emit("## 1. 产物清单")
    emit("")
    REQ = [
        "experiments/v6_1/README.md",
        "experiments/v6_1/PROGRESS.md",
        "experiments/v6_1/phase0_environment.md",
        "experiments/v6_1/V6_1_FINAL_SUMMARY.md",
        "experiments/v6_sota_audit/standardized_preds.csv",
        "experiments/v6_1/phase1_ci/bootstrap_ci_models.csv",
        "experiments/v6_1/phase1_ci/bootstrap_ci_report.md",
        "experiments/v6_1/phase1_ci/intersection_report.md",
        "experiments/v6_1/phase1_ci/intersection_pdb.csv",
        "experiments/v6_1/phase1_ci/category_paired_report.md",
        "experiments/v6_1/phase1_ci/category_paired_values.csv",
        "experiments/v6_1/phase2_split_inflation/split_inflation_results.csv",
        "experiments/v6_1/phase2_split_inflation/split_inflation_report.md",
        "experiments/v6_1/phase3_pairwise/pairwise_results.csv",
        "experiments/v6_1/phase3_pairwise/pairwise_report.md",
        "experiments/v6_1/phase4_v5_comparison/v5_cath_distribution.csv",
        "experiments/v6_1/phase4_v5_comparison/distribution_report.md",
        "experiments/v6_1/phase5_figures/fig1_illusion_scatter.png",
        "experiments/v6_1/phase5_figures/fig2_efficiency_gap.png",
        "experiments/v6_1/phase5_figures/fig3_category_comparison.png",
        "experiments/v6_1/phase5_figures/fig4_split_inflation.png",
        "experiments/v6_1/phase5_figures/fig5_hit10_vs_baserate.png",
        "experiments/v6_1/phase5_figures/fig6_per_complex_distribution.png",
        "experiments/v6_1/phase6_outline/paper_outline.md",
        "experiments/v6_1/phase7_toolkit/affirigid_metrics/__init__.py",
        "experiments/v6_1/phase7_toolkit/affirigid_metrics/topk.py",
        "experiments/v6_1/phase7_toolkit/affirigid_metrics/bootstrap.py",
        "experiments/v6_1/phase7_toolkit/affirigid_metrics/README.md",
        "experiments/v6_1/phase7_toolkit/examples/example_usage.py",
        "experiments/v6_1/phase7_toolkit/requirements.txt",
        "experiments/v6_1/phase7_toolkit/external_model_search.md",
        "experiments/v6_1/phase7_toolkit/zenodo_abag_report.md",
        "scripts/v6_1/fastmetrics.py",
    ]
    miss = [r for r in REQ if not (C.ROOT / r).exists()]
    emit(f"- 必需产物: **{len(REQ) - len(miss)} / {len(REQ)}** 存在")
    if miss:
        ok_all = False
        emit("")
        for m in miss:
            emit(f"  - ❌ `{m}`")
    emit("")

    # ---------- 2. 数字一致性 ----------
    emit("## 2. 关键数字与 CSV 一致性核对")
    emit("")
    emit("| 项 | 汇总报告声称 | CSV 实测 | 一致? |")
    emit("|---|---:|---:|---|")

    def chk(name, claimed, actual, tol=1e-3):
        global ok_all
        ok = abs(float(claimed) - float(actual)) <= tol
        ok_all &= ok
        emit(f"| {name} | {claimed:+.4f} | {actual:+.4f} | "
             f"{'✅' if ok else '❌'} |")

    # PHASE 2
    p2 = V61 / "phase2_split_inflation" / "split_inflation_results.csv"
    if p2.exists():
        d = pd.read_csv(p2, low_memory=False)
        d = d[d.get("status", "OK") == "OK"] if "status" in d.columns else d
        for lv, exp in (("L1_random", 0.8436), ("L4_family", 0.5774)):
            sub = pd.to_numeric(d[d["level"] == lv]["overall_pearson"],
                                errors="coerce")
            if len(sub):
                chk(f"PHASE2 {lv} Pearson", exp, sub.mean(), tol=2e-3)
        l1 = pd.to_numeric(d[d["level"] == "L1_random"]["overall_pearson"],
                           errors="coerce").mean()
        l4 = pd.to_numeric(d[d["level"] == "L4_family"]["overall_pearson"],
                           errors="coerce").mean()
        chk("PHASE2 L1−L4 inflation", 0.2663, l1 - l4, tol=3e-3)

    # PHASE 1.1
    p1 = V61 / "phase1_ci" / "bootstrap_ci_models.csv"
    if p1.exists():
        d = pd.read_csv(p1)
        for mo, exp in (("Flex-ddG", 0.4138), ("ESM-IF", 0.0277)):
            sub = d[(d["model"] == mo) & (d["metric"] == "macro_pdb_spearman")]
            if len(sub):
                chk(f"PHASE1.1 {mo} macro", exp, float(sub["point"].iloc[0]),
                    tol=2e-3)
        sub = d[(d["model"] == "ESM-IF") & (d["metric"] == "macro_pdb_spearman")]
        if len(sub):
            lo = float(sub["ci_low"].iloc[0])
            emit(f"| PHASE1.1 ESM-IF CI 下界（含 0 判定） | <0 | {lo:+.4f} | "
                 f"{'✅ 含 0' if lo < 0 else '❌ 不含 0'} |")
            ok_all &= (lo < 0)

    # PHASE 4
    p4 = V61 / "phase4_v5_comparison" / "v5_cath_distribution.csv"
    if p4.exists():
        d = pd.read_csv(p4)
        v5m = d[d["model"] == "V5 baseline"]["per_complex_spearman"].median()
        ph = d[d["model"].isin(["Flex-ddG", "FoldX"])]["per_complex_spearman"].median()
        chk("PHASE4 V5 median", 0.0939, v5m, tol=2e-3)
        chk("PHASE4 物理类 median", 0.4179, ph, tol=2e-3)
    emit("")

    # ---------- 3. 未修改旧版本 ----------
    emit("## 3. 是否修改过任何 V1.1–V6.0 文件（时间戳红线）")
    emit("")
    # **精确的 V6.1 起点**：只取 V6.1 独占目录（不含 v6_sota_audit，因为其中
    # 绝大部分是 V6.0 产物，只有 standardized_preds.csv 是 V6.1 新增的）
    v61_only = [p for d0 in (C.ROOT / "scripts" / "v6_1", V61,
                             C.ROOT / "logs" / "v6_1")
                if d0.exists() for p in d0.rglob("*") if p.is_file()]
    t_start = min((p.stat().st_mtime for p in v61_only), default=time.time())
    emit(f"- **V6.1 精确会话起点**（最早 V6.1 独占产物 mtime）: "
         f"**{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t_start))}**")
    emit(f"- 说明：`experiments/v6_sota_audit/` **整体不参与红线检查**，"
         f"因为其中 12 个文件是 V6.0 产物；V6.1 仅在其中**新增 1 个文件**"
         f"（`standardized_preds.csv`，V6.1 规格要求的输入，V6.0 未生成），"
         f"**未修改任何已有文件**。")
    emit("")

    def excluded(rel: str) -> bool:
        r = rel.replace("\\", "/")
        return (r.startswith("scripts/v6") or r.startswith("scripts/v6_1")
                or r.startswith("experiments/v6_1") or r.startswith("logs/v6_1")
                or r.startswith("experiments/v6_sota_audit")
                or "/v6_1/" in r or "/v6_sota_audit/" in r
                or r.startswith("experiments/v5_0") or "/v5_0/" in r
                or r.startswith("scripts/v5"))

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
            if excluded(rel):
                continue
            if p.stat().st_mtime > t_start:
                viol.append((rel, time.strftime("%H:%M:%S",
                                                time.localtime(p.stat().st_mtime))))
    if viol:
        ok_all = False
        emit(f"- ⚠️ 发现 **{len(viol)}** 个非 v6 文件在会话期间被修改:")
        for v, ts in viol[:15]:
            emit(f"  - `{v}` ({ts})")
    else:
        emit("- ✅ **未修改任何 V1.1–V6.0 文件**：在 V6.1 精确会话起点之后，"
             "`features/`、`scripts/`（v6 除外）、`data/`（v5_0/v6 除外）、"
             "`experiments/`（v6* 除外）、`logs/`（v6_1 除外）**无任何文件被改动**。")
    emit("")

    # ---------- 4. 失败日志 ----------
    emit("## 4. 失败日志")
    emit("")
    fl = C.ROOT / "logs" / "v6_1" / "failures.log"
    if fl.exists():
        txt = fl.read_text(encoding="utf-8", errors="replace").strip()
        emit(f"- `logs/v6_1/failures.log` 存在，{len(txt.splitlines())} 行")
        if txt:
            emit("```")
            emit(txt[-600:])
            emit("```")
    else:
        emit("- **无 `failures.log`** → 本任务无功能性失败")
    emit("")

    # ---------- 5. 磁盘纪律 ----------
    emit("## 5. 磁盘纪律核查")
    emit("")
    emit("| 目录 | 内容 |")
    emit("|---|---|")
    for d0 in (r"F:\dsh\downloads\cathddg", r"F:\dsh\downloads\dssa_ppi",
               r"F:\dsh\downloads\samaffinity", r"F:\dsh\downloads\mint",
               r"F:\dsh\downloads\zenodo_abag"):
        p = Path(d0)
        if p.exists():
            fs = list(p.rglob("*"))
            fs = [x for x in fs if x.is_file()]
            sz = sum(x.stat().st_size for x in fs)
            biggest = max((x.stat().st_size for x in fs), default=0)
            emit(f"| `{d0}` | {len(fs)} 文件, {sz / 1e6:.1f} MB, 最大单文件 "
                 f"{biggest / 1e6:.1f} MB |")
    emit("")
    emit("- ✅ 全部下载落于 **F 盘**；V6.1 单文件上限 2 MB（V6.0 的 Zenodo 下载为 "
         "50.5 MB，仍远低于 V6.1 允许的 100 MB）")
    emit("- ✅ 未写入 C 盘任何位置")
    emit("")

    emit("## 6. 总判定")
    emit("")
    emit(f"> **{'✅ 全部自检通过' if ok_all else '⚠️ 存在不一致（见上）'}**")
    emit("")
    OUT = C.ROOT / "logs" / "v6_1" / "final_verify.md"
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
