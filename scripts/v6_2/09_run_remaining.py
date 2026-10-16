# -*- coding: utf-8 -*-
"""scripts/v6_2/09_run_remaining.py -- V6.2 顺序驱动剩余 PHASE（3→5→6→7→4）.

用途：PHASE 2（4 档 × 3 seeds）耗时较长。本脚本在 PHASE 2 完成后
**顺序**执行其余 PHASE，每步记录耗时与退出码，失败不中断。

顺序说明
--------
* PHASE 3（Pearson 分解 + Primary Statistics）**依赖 PHASE 2 的 npz**，必须最先
* PHASE 5（bootstrap 修正）依赖 PHASE 2 npz + PHASE 3 的 primary_statistics.csv
* PHASE 6（leave-one-out）依赖 `standardized_preds.csv`（独立）
* PHASE 7（泄漏检查）依赖 PHASE 2 的 matched set + PHASE 3 的 sanity_checks.csv
* PHASE 4（排序 + δ）依赖 PHASE 2 npz（kendall/pairwise 较慢，放最后）
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
LOG = ROOT / "logs" / "v6_2"
LOG.mkdir(parents=True, exist_ok=True)

STEPS = [
    ("PHASE 3 — Pearson 分解 + Primary Statistics",
     ["scripts/v6_2/03_decomposition.py"]),
    ("PHASE 5 — Bootstrap 修正",
     ["scripts/v6_2/05_bootstrap_fixed.py"]),
    ("PHASE 6 — Leave-one-model-out + BH-FDR",
     ["scripts/v6_2/06_leave_one_out.py"]),
    ("PHASE 7 — 分级泄漏检查 + 浮点 sanity",
     ["scripts/v6_2/07_leakage_assertions.py"]),
    ("PHASE 4 — 排序指标 + δ 敏感性",
     ["scripts/v6_2/04_ranking_metrics.py"]),
]


def main() -> int:
    print("=" * 74)
    print("V6.2 剩余 PHASE 顺序驱动")
    print("=" * 74)
    results = []
    for name, cmd in STEPS:
        print()
        print("-" * 74)
        print(f">>> {name}")
        print(f"    python -u {' '.join(cmd)}")
        t = time.time()
        r = subprocess.run([sys.executable, "-u"] + cmd, cwd=str(ROOT),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace")
        dt = time.time() - t
        out = (r.stdout or "") + "\n--- STDERR ---\n" + (r.stderr or "")
        lp = LOG / (Path(cmd[0]).stem + ".driver.log")
        lp.write_text(out, encoding="utf-8")
        ok = r.returncode == 0
        print(f"    -> {'OK' if ok else f'FAIL rc={r.returncode}'}  ({dt:.0f}s)  "
              f"日志: {lp.name}")
        if not ok:
            tail = (r.stderr or "")[-700:]
            print("    stderr 尾部:")
            for x in tail.splitlines()[-10:]:
                print("      " + x)
        results.append((name, ok, dt))
    print()
    print("=" * 74)
    print("汇总")
    print("=" * 74)
    for n, ok, dt in results:
        print(f"  [{'OK  ' if ok else 'FAIL'}] {n}  ({dt:.0f}s)")
    allok = all(o for _, o, _ in results)
    print()
    print("全部成功 ✅" if allok else "存在失败 ❌（见上，已记录日志）")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
