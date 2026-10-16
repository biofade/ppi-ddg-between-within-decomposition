# -*- coding: utf-8 -*-
"""scripts/v5/00_env_snapshot.py -- V5.0 PHASE 0: 环境快照（只读，不修改任何已有文件）.

输出: experiments/v5_0/audit/environment.md
若文件已存在则追加时间戳段落（**不覆盖历史快照**）。
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
OUTDIR = ROOT / "experiments" / "v5_0" / "audit"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = OUTDIR / "environment.md"


def sh(cmd: list[str], cwd: Path | None = None) -> str:
    try:
        r = subprocess.run(cmd, cwd=str(cwd or ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace",
                           timeout=60)
        return (r.stdout or "").strip() + (("\n" + r.stderr.strip()) if r.stderr else "")
    except Exception as exc:  # noqa: BLE001
        return f"<failed: {exc}>"


def modver(name: str) -> str:
    try:
        m = __import__(name)
        return str(getattr(m, "__version__", "?"))
    except Exception:  # noqa: BLE001
        return "<not installed>"


def main() -> int:
    L: list[str] = []
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    L.append(f"## 环境快照 @ {ts}")
    L.append("")

    L.append("### OS / 硬件")
    L.append("")
    L.append("| 项 | 值 |")
    L.append("|---|---|")
    L.append(f"| OS | {platform.platform()} |")
    L.append(f"| 机器 | {platform.machine()} |")
    try:
        import psutil
        phys = psutil.cpu_count(logical=False)
        logi = psutil.cpu_count(logical=True)
        vm = psutil.virtual_memory()
        L.append(f"| CPU 物理核 / 逻辑核 | {phys} / {logi} |")
        L.append(f"| RAM 总量 | {vm.total / 1e9:.2f} GB |")
        L.append(f"| RAM 可用 | {vm.available / 1e9:.2f} GB |")
    except Exception:  # noqa: BLE001
        L.append("| CPU / RAM | psutil 未安装；见下方 PowerShell 记录 |")
    L.append(f"| CPU 型号 | {os.environ.get('PROCESSOR_IDENTIFIER', 'n/a')} |")
    L.append(f"| 磁盘可用 (项目盘) | {shutil.disk_usage(ROOT).free / 1e9:.1f} GB |")
    L.append("")

    L.append("### GPU / CUDA")
    L.append("")
    L.append(f"- `nvidia-smi`: {sh(['nvidia-smi', '--query-gpu=name,memory.total', '--format=csv'])[:200]}")
    try:
        import torch
        L.append(f"- torch: {torch.__version__} | CUDA available: {torch.cuda.is_available()} | "
                 f"device_count: {torch.cuda.device_count()}")
    except Exception as exc:  # noqa: BLE001
        L.append(f"- torch: <not importable: {exc}>")
    L.append("> 本项目全程 **CPU only**（见 V1.1–V4.0 记录）。")
    L.append("")

    L.append("### Python 环境")
    L.append("")
    L.append(f"- Python: {sys.version.split()[0]}")
    L.append(f"- 解释器: `{sys.executable}`")
    L.append(f"- 虚拟环境: {os.environ.get('VIRTUAL_ENV', '<none>')}")
    L.append("")
    L.append("| 包 | 版本 |")
    L.append("|---|---|")
    for m in ("numpy", "pandas", "scipy", "sklearn", "lightgbm", "torch",
              "Bio", "matplotlib", "pypdf", "psutil"):
        L.append(f"| {m} | {modver(m)} |")
    L.append("")

    L.append("### 版本控制")
    L.append("")
    if (ROOT / ".git").exists():
        L.append(f"- git HEAD: `{sh(['git', 'rev-parse', 'HEAD'])}`")
        L.append(f"- git branch: `{sh(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])}`")
        porc = sh(["git", "status", "--porcelain"])
        L.append(f"- git status --porcelain 行数: {len(porc.splitlines())}")
    else:
        L.append("- **不是 git 仓库**（无 `.git`）→ commit hash 不可用，标记 `GIT_NOT_AVAILABLE`")
    L.append("")

    L.append("### 时间与运行成本")
    L.append("")
    L.append(f"- 快照时刻: {ts}")
    L.append(f"- 项目根: `{ROOT}`")
    L.append("")

    header = ""
    if not OUT.exists():
        header = ("# V5.0 PHASE 0 — 环境快照\n\n"
                  "> 由 `scripts/v5/00_env_snapshot.py` 生成。"
                  "**每次运行追加一段，不覆盖历史。**\n\n")
    with open(OUT, "a", encoding="utf-8") as fh:
        fh.write(header + "\n".join(L) + "\n\n---\n\n")
    print("\n".join(L))
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
