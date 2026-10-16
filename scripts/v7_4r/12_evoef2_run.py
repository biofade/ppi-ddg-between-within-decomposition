# -*- coding: utf-8 -*-
"""scripts/v7_4r/12_evoef2_run.py -- V7.4R Track B：EvoEF2 执行（restart-safe）.

在 **WSL2 / Linux / Colab CPU** 下调用 `EvoEF2` 二进制。**Windows 原生不推荐**
（上游面向 Linux/g++）。GPU 不需要（§34）。

流程（§35）
----------
1. 每个 target 只 `RepairStructure` **一次**（缓存 `<T>_protein_Repair.pdb`）
2. 逐变体 `BuildMutant`（`individual_list.txt` 每行一个 token；输出 `Model_1.pdb`）
3. WT 与 MUT 各 `ComputeBinding` 一次 → `dE_bind = E_MUT - E_WT`
4. 逐变体落盘 + manifest（可续跑）

用法（在 WSL2/Linux 下）
------------------------
    python scripts/v7_4r/12_evoef2_run.py --evoef2 /path/to/EvoEF2 --smoke
    python scripts/v7_4r/12_evoef2_run.py --evoef2 /path/to/EvoEF2
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "evoef2" / "v7_4r"
STRUCT = OUT / "structures"
RUNS = OUT / "runs"
RUNS.mkdir(parents=True, exist_ok=True)
PH3 = R / "experiments" / "v7_4r" / "phase3_evoef2"
PH3.mkdir(parents=True, exist_ok=True)

# EvoEF2 输出里的能量行（宽松匹配，避免对其格式做死假设）
NUM = re.compile(r"[-+]?\d+\.\d+")


def run(cmd, cwd=None, timeout=1800):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=timeout)
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def parse_total_energy(text: str):
    """从 ComputeBinding/ComputeStability 输出中提取 'Total' 能量。

    不假定列顺序：找含 'total' 的行，取最后一个浮点数。
    """
    best = None
    for ln in text.splitlines():
        low = ln.lower()
        if "total" not in low:
            continue
        nums = NUM.findall(ln)
        if nums:
            try:
                best = float(nums[-1])
            except ValueError:
                pass
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evoef2", required=True, help="EvoEF2 可执行文件路径")
    ap.add_argument("--smoke", action="store_true", help="只做少量变体")
    ap.add_argument("--limit", type=int, default=0, help="额外限制变体数")
    ap.add_argument("--num-runs-repair", type=int, default=3)
    ap.add_argument("--num-runs-mutant", type=int, default=10)
    a = ap.parse_args()

    exe = Path(a.evoef2)
    if not exe.exists():
        print(f"EvoEF2 不存在: {exe}"); return 2
    idx = pd.read_csv(OUT / "evoef2_variant_index.csv")
    if a.limit:
        idx = idx.head(a.limit)
    if a.smoke:
        idx = idx.head(5)
    print(f"EvoEF2: {exe}")
    print(f"变体数: {len(idx)} | targets: {idx.target_id.nunique()}")

    man_path = RUNS / "evoef2_manifest.csv"
    man = pd.read_csv(man_path) if man_path.exists() else pd.DataFrame(
        columns=["row_id", "target_id", "WT_status", "MUT_status", "E_bind_WT",
                 "E_bind_MUT", "dE_bind", "runtime", "error_code"])
    done = set(man.loc[man.WT_status == "COMPLETE", "row_id"].astype(int)) \
        if len(man) else set()
    print("已完成:", len(done))

    repaired: dict[str, Path] = {}
    t0 = time.time()
    for t, g in idx.groupby("target_id"):
        pdb = STRUCT / f"{t}_protein.pdb"
        if not pdb.exists():
            print(f"  缺结构 {pdb}"); continue
        # 1) repair 一次
        rep = STRUCT / f"{t}_protein_Repair.pdb"
        if not rep.exists():
            rc, out = run([str(exe), "--command=RepairStructure",
                           f"--pdb={pdb.name}", f"--num_of_runs={a.num_runs_repair}"],
                          cwd=str(STRUCT))
            if rc != 0 or not rep.exists():
                print(f"  repair FAIL {t}: rc={rc} {out[-200:]}")
                continue
        repaired[t] = rep

        # 2) 逐变体 build + 3) binding
        for r in g.itertuples():
            rid = int(r.row_id)
            if rid in done:
                continue
            rec = {"row_id": rid, "target_id": t, "WT_status": "RUNNING",
                   "MUT_status": "PENDING", "E_bind_WT": None, "E_bind_MUT": None,
                   "dE_bind": None, "runtime": 0.0, "error_code": ""}
            try:
                tb = time.time()
                ml = RUNS / f"{t}_{rid}_individual_list.txt"
                ml.write_text(str(r.evoef2_mutation) + "\n", encoding="utf-8")
                rc, out = run([str(exe), "--command=BuildMutant",
                               f"--pdb={rep.name}", f"--mutant_file={ml.name}",
                               f"--num_of_runs={a.num_runs_mutant}"],
                              cwd=str(RUNS))
                model = RUNS / f"{rep.stem}_Model_1.pdb"
                if rc != 0 or not model.exists():
                    raise RuntimeError(f"BuildMutant rc={rc} {out[-160:]}")
                # WT binding
                rc1, o1 = run([str(exe), "--command=ComputeBinding",
                               f"--pdb={rep.name}"], cwd=str(STRUCT))
                e_wt = parse_total_energy(o1)
                # MUT binding
                rc2, o2 = run([str(exe), "--command=ComputeBinding",
                               f"--pdb={model.name}"], cwd=str(RUNS))
                e_mu = parse_total_energy(o2)
                if e_wt is None or e_mu is None:
                    raise RuntimeError("energy parse failed")
                rec.update(WT_status="COMPLETE", MUT_status="COMPLETE",
                           E_bind_WT=e_wt, E_bind_MUT=e_mu, dE_bind=e_mu - e_wt,
                           runtime=round(time.time() - tb, 1))
                done.add(rid)
            except Exception as e:  # noqa: BLE001
                rec.update(WT_status="FAIL_BUILD",
                           error_code=f"{type(e).__name__}:{str(e)[:80]}")
            man = man[man.row_id.astype(int) != rid] if len(man) else man
            man = pd.concat([man, pd.DataFrame([rec])], ignore_index=True)
            if len(man) % 10 == 0:
                man.to_csv(man_path, index=False)
        if (time.time() - t0) / 60 > 1:
            print(f"  [{t}] rows={len(man)} "
                  f"elapsed={(time.time() - t0) / 60:.1f} min")
    man.to_csv(man_path, index=False)
    ok = man[man.WT_status == "COMPLETE"]
    print(f"\n完成 {len(ok)}/{len(idx)}")
    print(f"manifest -> {man_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
