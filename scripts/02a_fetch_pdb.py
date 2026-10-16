# -*- coding: utf-8 -*-
"""
02a_fetch_pdb.py -- Affinity Intelligence / 阶段二-结构抓取

从 RCSB 按需下载 mvp_structures.txt 中的 PDB 结构到 data/pdb/。
纪律: 并发数强制 2, 每个文件独立 try/except, 403/404 直接放弃并记日志,
      不重试超过 1 次, 已存在且非空则跳过。
"""

import argparse
import os
import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
DATA_DIR = os.path.join(ROOT, "data")
PDB_DIR = os.path.join(DATA_DIR, "pdb")
MVP_LIST = os.path.join(DATA_DIR, "mvp_structures.txt")
FULL_LIST = os.path.join(DATA_DIR, "full_structures.txt")
BENCH = os.path.join(DATA_DIR, "benchmark.csv")
LOG_PATH = os.path.join(DATA_DIR, "fetch_log.txt")

os.makedirs(PDB_DIR, exist_ok=True)

logger = logging.getLogger("fetch")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)

HEADERS = {"User-Agent": "AffinityIntelligence/0.1 (academic pipeline; contact: local)"}
SOURCES = [
    "https://files.rcsb.org/download/{pid}.pdb",
    "https://files.rcsb.org/download/{pid}.cif",
]
TIMEOUT = 45


def fetch_one(pid):
    """返回 (pid, path, bytes, source) 或 (pid, None, 0, err)。"""
    for ext_try in SOURCES:
        url = ext_try.format(pid=pid)
        dest = os.path.join(PDB_DIR, url.rsplit("/", 1)[-1])
        if os.path.exists(dest) and os.path.getsize(dest) > 500:
            return pid, dest, os.path.getsize(dest), "cached"
        for attempt in (1, 2):
            try:
                r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
                if r.status_code == 200 and len(r.content) > 500:
                    tmp = dest + ".part"
                    with open(tmp, "wb") as fh:
                        fh.write(r.content)
                    os.replace(tmp, dest)
                    return pid, dest, len(r.content), url.split("/")[-1]
                if r.status_code in (404, 403):
                    break            # 该格式不存在, 换下一种格式
            except Exception as exc:                          # noqa: BLE001
                if attempt == 2:
                    logger.warning("下载失败 %s (%s): %s", pid, url.split("/")[-1], exc)
                time.sleep(1.0)
    return pid, None, 0, "not_found"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", type=str, default="",
                    help="结构 ID 清单文件 (默认 data/mvp_structures.txt)")
    ap.add_argument("--all", action="store_true",
                    help="从 benchmark.csv 取全部非 MVP PDB (写入 data/full_structures.txt)")
    args = ap.parse_args()

    t0 = time.time()
    list_path = args.list or MVP_LIST
    if args.all:
        import pandas as pd
        bench = pd.read_csv(BENCH, dtype={"pdb_id": str})
        rest = bench[bench["is_mvp"] != True]                    # noqa: E712
        mvp = set(bench[bench["is_mvp"] == True]["pdb_id"])       # noqa: E712
        pids = sorted(set(rest["pdb_id"]) - mvp)
        with open(FULL_LIST, "w", encoding="utf-8") as fh:
            fh.write("\n".join(pids) + "\n")
        list_path = FULL_LIST
        logger.info("从 benchmark.csv 生成全量清单: %d 个 PDB -> %s",
                    len(pids), FULL_LIST)
    else:
        with open(list_path, "r", encoding="utf-8") as fh:
            pids = [ln.strip() for ln in fh if ln.strip()]

    need = [p for p in pids if not _present(p)]
    logger.info("清单 %d 个 PDB, 其中 %d 个待抓取 -> %s",
                len(pids), len(need), PDB_DIR)
    if not need:
        print("全部结构已存在, 无需抓取")
        return

    ok, fail = 0, []
    done = 0
    with ThreadPoolExecutor(max_workers=2) as pool:     # 并发强制 = 2
        futs = {pool.submit(fetch_one, p): p for p in need}
        for fut in as_completed(futs):
            pid = futs[fut]
            done += 1
            try:
                pid, path, size, src = fut.result()
            except Exception as exc:                      # noqa: BLE001
                fail.append((pid, str(exc)[:80]))
                continue
            if path:
                ok += 1
            else:
                fail.append((pid, src))
            if done % 100 == 0 or done == len(need):
                logger.info("进度 %d/%d  成功=%d 失败=%d  已用时 %.0fs",
                            done, len(need), ok, len(fail), time.time() - t0)

    if fail:
        logger.warning("失败清单 (%d): %s", len(fail), fail[:40])
    print("\n抓取完成: 成功 %d / %d, 失败 %d, 用时 %.0fs"
          % (ok, len(need), len(fail), time.time() - t0))
    print("目录: %s" % PDB_DIR)


def _present(pid):
    for ext in (".pdb", ".cif", ".ent"):
        p = os.path.join(PDB_DIR, pid + ext)
        if os.path.exists(p) and os.path.getsize(p) > 500:
            return True
    return False


if __name__ == "__main__":
    main()
