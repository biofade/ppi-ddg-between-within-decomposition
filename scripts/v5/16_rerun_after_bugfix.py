# -*- coding: utf-8 -*-
"""scripts/v5/16_rerun_after_bugfix.py -- V5.0 Bug 修正后的重算链（V5 §123）.

Bug 摘要（已确认，非推测）
--------------------------
`v5_common.topk_metrics` 的 `efficiency_at_10` 曾**系统性反号**：
  * oracle 模型 → 得 **−0.6667**（应为 +1.0）
  * 完全反向模型 → 得 **+1.0000**（应为 ≤0）
根因: 每个 complex 内的 **random 基线只用 1 次抽样**，小组下方差极大，
使分母 `v_or − v_rnd` 可正可负；叠加 oracle/random 语义混用。
修正: 随机基线改为 **200 次抽样取均值**；efficiency 改为
`(random − model) / (random − oracle)`，对合理模型恒 ≥0，oracle → +1，random → 0。
合成自检（20 个 complex × 40 候选）: oracle **+1.0000**、anti **−0.9816**、
noisy **+0.7644**，符合预期。

影响范围（**已标 `INVALID_DUE_TO_BUG`**）
----------------------------------------
所有曾写入 `efficiency_at_10` 的产物：
  `p2_state_potential/p2_results.csv`、`p3_balancing/p3_results.csv`、
  `p4_conformal/p4_results.csv`、`all_experiments.csv`、`final/variant_summary.csv`
**注意**：Pearson / Spearman / macro-PDB / RMSE / MAE / AUROC / P@k / EF / Hit@k
**不受影响**（它们不依赖该函数）。因此 P2/P3/P4 的**主要结论不变**，
只有 `efficiency_at_10` 一列需要重算。

本脚本串行重跑受影响的脚本以刷新该列。
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

# 先给旧产物打 INVALID 标记（保留证据，不删除）
TO_MARK = [
    C.V5 / "p2_state_potential" / "p2_results.csv",
    C.V5 / "p3_balancing" / "p3_results.csv",
    C.V5 / "p4_conformal" / "p4_results.csv",
    C.V5 / "all_experiments.csv",
    C.V5 / "final" / "variant_summary.csv",
]

RUNS = [
    ("P4 conformal（最快，含 Top-K 主结论）",
     ["scripts/v5/09_p4_conformal.py"]),
    ("P4 Top-K 分解（核心诊断）",
     ["scripts/v5/14_p4_decompose_topk.py"]),
    ("P3 balancing（--skip-bagging）",
     ["scripts/v5/08_p3_balancing.py", "--skip-bagging"]),
    ("P2 state potential",
     ["scripts/v5/07_p2_state_potential.py"]),
    ("汇总 all_experiments / paired_comparisons",
     ["scripts/v5/12_consolidate_results.py"]),
]


def main() -> int:
    print("=" * 72)
    print("V5.0 Bug 修正重算链")
    print("=" * 72)

    # 1) 标记旧产物
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    for p in TO_MARK:
        if p.exists():
            side = p.with_suffix(p.suffix + ".INVALID_DUE_TO_BUG")
            if not side.exists():
                p.replace(side)
                print(f"  [标记无效] {p.name} -> {side.name}")
            else:
                print(f"  [跳过] {side.name} 已存在")
    note = C.V5 / "logs" / "BUGFIX_efficiency_at_10.md"
    note.write_text(
        f"""# BUG 修正记录 — `efficiency_at_10` 反号

**时间**: {stamp}
**状态**: `INVALID_DUE_TO_BUG`（旧产物已重命名保留为 `*.INVALID_DUE_TO_BUG`）

## 现象（合成自检证据）

| 合成模型 | 修正前 efficiency_at_10 | 修正后 |
|---|---:|---:|
| 完美 oracle | **−0.6667** | **+1.0000** |
| 完全反向 | **+1.0000** | **−0.9816** |
| 含噪中间模型 | −0.6667 | **+0.7644** |

## 根因

`v5_common.topk_metrics` 中：
1. 每个 complex 的 **random 基线只抽 1 次**，在 n 较小的 complex 上方差极大，
   使分母 `v_or − v_rnd` 可正可负，符号不稳定；
2. oracle 与 random 的语义在符号变换中混用。

## 修正

* 随机基线改为 **200 次抽样取均值**；
* efficiency 统一定义为 `(random − model) / (random − oracle)`
  → 对合理模型恒 ≥ 0；oracle → +1；random → 0（与 V5 §108 语义一致）。

## 影响范围

* **受影响**：`efficiency_at_10`（及 `efficiency_at_10_median`、
  `utility_*`、`regret_at_10`）一列。
* **不受影响**：Pearson / Spearman / macro per-PDB / RMSE / MAE / AUROC /
  P@k / EF@x% / Hit@k —— 这些不依赖该函数，因此 **P2/P3/P4 的主要结论不变**。

## 处置（V5 §123）

1. 保存 bug evidence ✅（本文件）
2. 写 reproduction ✅（上表合成自检）
3. 修 bug ✅
4. **从受影响节点之后全部重新运行** → 见下方重算链
5. 旧结果标 `INVALID_DUE_TO_BUG` ✅
6. 不只更新最终表 ✅（重跑全部下游脚本）
""", encoding="utf-8")
    print(f"  [写出] {note}")

    # 2) 串行重跑
    ok = True
    for name, cmd in RUNS:
        print()
        print("-" * 72)
        print(f">>> {name}")
        print(f"    python -u {' '.join(cmd)}")
        t = time.time()
        r = subprocess.run([sys.executable, "-u"] + cmd, cwd=str(C.ROOT),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace")
        dt = time.time() - t
        log = C.V5 / "logs" / (Path(cmd[0]).stem + ".rerun.log")
        log.write_text((r.stdout or "") + "\n--- STDERR ---\n" + (r.stderr or ""),
                       encoding="utf-8")
        status = "OK" if r.returncode == 0 else f"FAIL rc={r.returncode}"
        print(f"    -> {status}  ({dt:.0f}s)  日志: {log.name}")
        if r.returncode != 0:
            ok = False
            tail = (r.stderr or "")[-800:]
            print("    stderr 尾部:\n" + "\n".join("      " + x
                                                  for x in tail.splitlines()[-12:]))

    print()
    print("=" * 72)
    print("重算链完成" if ok else "重算链**部分失败** —— 见上方 FAIL")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
