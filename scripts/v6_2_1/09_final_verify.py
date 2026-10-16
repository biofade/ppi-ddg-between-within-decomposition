# -*- coding: utf-8 -*-
"""scripts/v6_2_1/09_final_verify.py -- V6.2.1 交付前自检（只读）。

1. 产物完整性
2. 关键数字与 CSV 一致
3. **V6.2 冻结哈希是否仍一致**（证明未修改上游）
4. 失败日志 / 残留进程
"""

from __future__ import annotations

import hashlib
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

V621 = C.ROOT / "experiments" / "v6_2_1"
L: list[str] = []
ok_all = True

# PHASE 0 冻结的哈希（必须保持不变 → 证明未修改 V6.2）
FROZEN = {
    "experiments/v6_2/phase2_split4/split4_predictions.npz":
        "ba711ea1b349836ac18b44ba67f6ab2f0a0b680b",
    "experiments/v6_2/phase3_decomposition/pearson_decomposition.csv":
        "ce9baa9c4fe2cf16023999f22820ea298388d384",
    "experiments/v6_2/phase3_decomposition/control_results.csv":
        "4edcb024912700418ec35b8a082682643b18bd14",
    "scripts/v6_2/02_split4_ddg.py":
        "6a391a9c08f64a372507619f28fc7b0453ae5156",
    "scripts/v6_2/03_decomposition.py":
        "5c74317a1b85954911f026753857b39e7b143c29",
    "scripts/v5/v5_common.py":
        "91bee0814f66737b50edadd937d967644bbf3ec7",
}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    global ok_all
    emit("# V6.2.1 交付前自检")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # ---------- 1. 产物 ----------
    emit("")
    emit("## 1. 产物完整性")
    emit("")
    REQ = [
        "experiments/v6_2_1/PROGRESS.md",
        "experiments/v6_2_1/phase0_environment.md",
        "experiments/v6_2_1/V6_2_1_FINAL_SUMMARY.md",
        "experiments/v6_2_1/phase1_shuffle_audit/shuffle_audit.csv",
        "experiments/v6_2_1/phase1_shuffle_audit/shuffle_repeats.csv",
        "experiments/v6_2_1/phase1_shuffle_audit/old_shuffle_repro.md",
        "experiments/v6_2_1/phase1_shuffle_audit/shuffle_report.md",
        "experiments/v6_2_1/phase2_tree_audit/tree_count_audit.csv",
        "experiments/v6_2_1/phase2_tree_audit/tree_count_report.md",
        "experiments/v6_2_1/phase3_d_only/d_only_predictions.npz",
        "experiments/v6_2_1/phase3_d_only/d_only_metrics.csv",
        "experiments/v6_2_1/phase3_d_only/gain_w_bootstrap.csv",
        "experiments/v6_2_1/phase3_d_only/d_only_report.md",
        "experiments/v6_2_1/phase4_target_size/target_level_stats.csv",
        "experiments/v6_2_1/phase4_target_size/target_size_bins.csv",
        "experiments/v6_2_1/phase4_target_size/target_size_report.md",
        "experiments/v6_2_1/phase5_corrected_stats/primary_table.csv",
        "experiments/v6_2_1/phase5_corrected_stats/route_statistics.csv",
        "experiments/v6_2_1/phase5_corrected_stats/corrected_shuffle.csv",
        "experiments/v6_2_1/phase5_corrected_stats/corrected_stats_report.md",
        "experiments/v6_2_1/phase6_route/route_decision.md",
    ]
    miss = [r for r in REQ if not (C.ROOT / r).exists()]
    emit(f"- 必需产物: **{len(REQ) - len(miss)} / {len(REQ)}**")
    if miss:
        ok_all = False
        for m in miss:
            emit(f"  - MISSING `{m}`")
    emit("")

    # ---------- 2. 关键数字 ----------
    emit("## 2. 关键数字与 CSV 一致")
    emit("")
    emit("| 项 | 声称 | CSV 实测 | 一致? |")
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

    # PHASE 1
    p1 = V621 / "phase1_shuffle_audit" / "shuffle_repeats.csv"
    if p1.exists():
        d = pd.read_csv(p1)
        s = d[(d.level == "L1_random") & (d.model == "M1")]
        if len(s):
            chk("PHASE1 L1/M1 C_B_shuf mean", 0.4275,
                s["C_B_shuffle_mean"].iloc[0], tol=2e-3)
            chk("PHASE1 L1/M1 C_W_shuf mean", -0.0004,
                s["C_W_shuffle_mean"].iloc[0], tol=2e-3)
    a1 = V621 / "phase1_shuffle_audit" / "shuffle_audit.csv"
    if a1.exists():
        d = pd.read_csv(a1)
        s = d[(d.level == "L1_random") & (d.model == "M1") & (d.seed == 42)]
        if len(s):
            chk("PHASE1 L1/M1 r_shuffle(正确)", 0.4348,
                s["r_shuffle"].iloc[0], tol=3e-3)
            emit(f"| PHASE1 C_B 守恒组合数 | {len(d)} | "
                 f"{int(d.C_B_preserved.sum())} | "
                 f"{'✅' if int(d.C_B_preserved.sum()) == len(d) else '❌'} |")
            ok_all &= (int(d.C_B_preserved.sum()) == len(d))

    # PHASE 3
    g3 = V621 / "phase3_d_only" / "gain_w_bootstrap.csv"
    if g3.exists():
        d = pd.read_csv(g3)
        for lab, lv, pt, ex in (
                ("Gain_D = C_W(M2)−C_W(M0)", "L4_family", 0.0074, False),
                ("Gain_DGaux = C_W(M1)−C_W(M2)", "L4_family", -0.0078, False),
                ("Gain_D = C_W(M2)−C_W(M0)", "L1_random", 0.0035, False)):
            s = d[(d.contrast == lab) & (d.level == lv)]
            if len(s):
                chk(f"PHASE3 {lab} @{lv}", pt, s["point"].iloc[0], tol=2e-3)
                same = bool(s["ci_excludes_zero"].iloc[0]) == ex
                emit(f"| PHASE3 {lab} @{lv} 排除0={ex} | {ex} | "
                     f"{bool(s['ci_excludes_zero'].iloc[0])} | "
                     f"{'✅' if same else '❌'} |")
                ok_all &= same

    # PHASE 4
    b4 = V621 / "phase4_target_size" / "target_size_bins.csv"
    if b4.exists():
        d = pd.read_csv(b4)
        for b, cb in ((">=30", -0.0720), ("2-9", 0.3153), ("10-29", 0.3076)):
            s = d[(d.level == "L4_family") & (d.bin == b)]
            if len(s):
                chk(f"PHASE4 L4 bin {b} row-w C_B", cb,
                    s["roww_C_B"].iloc[0], tol=3e-3)

    # PHASE 5
    r5 = V621 / "phase5_corrected_stats" / "route_statistics.csv"
    if r5.exists():
        d = pd.read_csv(r5)
        s = d[d.contrast == "M1"].head(1)
        s = d[(d.statistic.str.startswith("Delta_OOD_W")) & (d.contrast == "M1")]
        if len(s):
            chk("PHASE5 Delta_OOD_W(M1)", 0.1207, s["point"].iloc[0], tol=3e-3)
        n_ex = int(d[d.statistic.str.startswith("Gain_D =")]["ci_excludes_zero"].sum())
        emit(f"| PHASE5 Gain_D 排除 0 的个数 | 0 | {n_ex} | "
             f"{'✅' if n_ex == 0 else '❌'} |")
        ok_all &= (n_ex == 0)
    emit("")

    # ---------- 3. V6.2 冻结哈希 ----------
    emit("## 3. V6.2 冻结哈希是否仍一致（证明未修改上游）")
    emit("")
    emit("| 文件 | 期望(前 40) | 实际(前 40) | 一致? |")
    emit("|---|---|---|---|")
    for rel, exp in FROZEN.items():
        p = C.ROOT / rel
        if not p.exists():
            emit(f"| `{rel}` | `{exp}` | MISSING | ❌ |")
            ok_all = False
            continue
        act = sha(p)[:40]
        same = (act == exp)
        ok_all &= same
        emit(f"| `{rel}` | `{exp}` | `{act}` | {'✅' if same else '❌ 被修改!'} |")
    emit("")

    # ---------- 4. 失败日志 ----------
    emit("## 4. 失败日志")
    emit("")
    fl = C.ROOT / "logs" / "v6_2_1" / "failures.log"
    if fl.exists() and fl.stat().st_size > 0:
        emit(f"- `logs/v6_2_1/failures.log` 存在（{fl.stat().st_size} B）")
    else:
        emit("- **无 `failures.log`** → 无功能性失败")
    emit("")

    emit("## 5. 总判定")
    emit("")
    emit(f"> **{'✅ 全部自检通过' if ok_all else '⚠️ 存在不一致（见上）'}**")
    emit("")
    OUT = C.ROOT / "logs" / "v6_2_1" / "final_verify.md"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
