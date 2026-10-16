# -*- coding: utf-8 -*-
"""scripts/v2.3.2/run_nn_sequential.py -- 顺序跑 N1 -> N2 (避免并发 OOM).

背景
----
本机只有 7.5 GB 物理内存。并发跑两个 NN 脚本会在加载 ESM-2 缓存 +
PCA (640x640 特征值分解) 时 OOM (实测: `_ArrayMemoryError`, 连 20 MiB 都分配不出)。
因此改为**顺序执行**, 并在进入前显式限制 BLAS 线程数以降低内存峰值。

用法
----
    python scripts/v2.3.2/run_nn_sequential.py [--raw640] [--seeds 42 43 44 45]

行为
----
1. 先跑 `nn_mlp.py`  (N1)
2. 再跑 `nn_two_tower.py` (N2)
3. 任一步失败**不中断**后一步 (硬约束 4: 单步失败不中断全流程),
   失败信息写入 `logs/v2.3.2/failures.log`。
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
LOG_DIR = ROOT / "logs" / "v2.3.2"
FAIL_LOG = LOG_DIR / "failures.log"
PY = sys.executable


def run(script: str, args: list[str], log_name: str) -> int:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [PY, "-u", str(ROOT / "scripts" / "v2.3.2" / script)] + args
    print(f"\n{'='*78}\n[{time.strftime('%H:%M:%S')}] RUN {script} {' '.join(args)}\n{'='*78}")
    env = dict(os.environ)
    # 限制 BLAS/OpenMP 线程: 降低内存峰值, 也避免与其他任务争抢
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        env[k] = "2"
    env["PYTHONIOENCODING"] = "utf-8"
    log_path = LOG_DIR / log_name
    with open(log_path, "w", encoding="utf-8", errors="replace") as fh:
        proc = subprocess.run(cmd, cwd=str(ROOT), env=env,
                              stdout=fh, stderr=subprocess.STDOUT)
    print(f"[{time.strftime('%H:%M:%S')}] {script} exit={proc.returncode} "
          f"(log: {log_path.name})")
    # 回显最后几行, 便于在主控窗口看到进度
    try:
        tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-8:]
        for ln in tail:
            print("   ", ln)
    except Exception:  # noqa: BLE001
        pass
    if proc.returncode != 0:
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[run_nn_sequential] {script} exit={proc.returncode}\n")
    return proc.returncode


def main() -> int:
    ap = argparse.ArgumentParser(description="顺序跑 N1 -> N2")
    ap.add_argument("--raw640", action="store_true")
    ap.add_argument("--seeds", type=int, nargs="*", default=None)
    args = ap.parse_args()

    common: list[str] = []
    if args.raw640:
        common.append("--raw640")
    if args.seeds:
        common += ["--seeds"] + [str(s) for s in args.seeds]

    rc1 = run("nn_mlp.py",
              common + ["--out-csv", "features/v2.3.2_nn_mlp_results.csv",
                        "--preds-npz", "features/v2.3.2_nn_mlp_predictions.npz"],
              "nn_mlp_console.log")
    rc2 = run("nn_two_tower.py",
              common + ["--out-csv", "features/v2.3.2_nn_twotower_results.csv",
                        "--preds-npz", "features/v2.3.2_nn_twotower_predictions.npz"],
              "nn_twotower_console.log")
    print(f"\n=== 完成: N1 exit={rc1}, N2 exit={rc2} ===")
    return 0 if (rc1 == 0 or rc2 == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
