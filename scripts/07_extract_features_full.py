# -*- coding: utf-8 -*-
"""
07_extract_features_full.py -- Affinity Intelligence V1.0 / 全量特征提取 (断点续算)

输入 : data/benchmark.csv                 (全量 4762 行; 默认排除 is_mvp 的 400 行)
       data/pdb/*.pdb|*.cif               (结构文件, 由 02a_fetch_pdb.py 预先抓取)
输出 : features/features_full.csv         (每 --exit-every 条结构落盘一次; 可续算)
       features/extract_checkpoint.json   (进度 + 失败清单; 重启即续算)
       features/features_full_log.txt

断点续算机制:
  * 进度以 "结构 (PDB ID)" 为单位记录: done_pdb / failed_pdb / skipped_pdb。
  * 每处理 --exit-every(=100) 条结构, 立即 (a) 追加写出新结果,
    (b) 原子替换 (os.replace) 写检查点。进程被杀最多丢最近不到一批的结果。
  * 启动时读检查点跳过已完成/已失败结构; 同时读回既有 CSV, 以 interface_key
    去重 (保留最后一条), 保证输出文件不重复、不丢行。
  * CSV 与检查点不一致时 (例如检查点落后), 以 CSV 实际内容为准补记进度。
  * 失败结构写入 failed_pdb 并带原因, 重启不再重试 (不硬撑、不反复重试)。

速度策略 (纯 CPU):
  * --pair-sasa auto (默认): 单链路对时逐链对量在数学上等于整体量, 精确直填,
    省掉 2~3 次 FreeSASA 调用 (零信息损失); 多链路对才真正计算, 上限 8 个链对。
  * --pair-sasa N: 强制上限 N; --pair-sasa 0: 完全跳过链路对计算 (最快)。
  * 每批结束 gc.collect(); 每批打印进度、内存与 ETA。
"""

import argparse
import gc
import json
import logging
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# ---------------------------------------------------------------- 路径 / 常量
ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
DATA_DIR = os.path.join(ROOT, "data")
FEAT_DIR = os.path.join(ROOT, "features")
PDB_DIR = os.path.join(DATA_DIR, "pdb")
BENCH = os.path.join(DATA_DIR, "benchmark.csv")

OUT_CSV = os.path.join(FEAT_DIR, "features_full.csv")
CHECKPOINT = os.path.join(FEAT_DIR, "extract_checkpoint.json")
LOG_PATH = os.path.join(FEAT_DIR, "features_full_log.txt")

MAX_PAIR_SASA = 8          # 用户指令: 逐链对 ΔSASA 最多算前 8 个链对
EXIT_EVERY = 100           # 默认每 100 条结构落盘一次
MEM_WARN_MB = 4000

os.makedirs(FEAT_DIR, exist_ok=True)

logger = logging.getLogger("full")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)


def get_rss_mb():
    """当前进程常驻内存 (MB); 失败返回 NaN。"""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return pmc.WorkingSetSize / 1024 / 1024
    except Exception:                                           # noqa: BLE001
        return float("nan")


def resolve_struct_path(pdb_id):
    """定位结构文件; 不存在返回 None。"""
    for ext in (".pdb", ".cif", ".ent"):
        p = os.path.join(PDB_DIR, pdb_id + ext)
        if os.path.exists(p) and os.path.getsize(p) > 500:
            return p
    return None


def _load_feature_module():
    """复用 02_extract_features.py 的特征函数 (单一实现来源, 避免两份代码漂移)。"""
    import importlib.util
    p = os.path.join(ROOT, "scripts", "02_extract_features.py")
    spec = importlib.util.spec_from_file_location("feat02", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_FEAT02 = _load_feature_module()


def fill_single_pair_features(feat):
    """单链路对: 逐链对量在数学上等于整体量, 精确直填 (零信息损失)。"""
    tot = feat.get("dsasa_total", np.nan)
    feat["weakest_chainpair_dsasa"] = tot
    feat["strongest_chainpair_dsasa"] = tot
    feat["mean_chainpair_dsasa"] = tot
    feat["std_chainpair_dsasa"] = 0.0
    feat["weakest_over_total_dsasa"] = 1.0
    feat["pair_sasa_complete"] = 1
    return feat


def nan_pair_features(feat):
    """多链路对但未计算链路对量: 置 NaN 并标注不完整 (不伪装成有效值)。"""
    for k in ("weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
              "mean_chainpair_dsasa", "std_chainpair_dsasa",
              "weakest_over_total_dsasa"):
        feat[k] = np.nan
    feat["pair_sasa_complete"] = 0
    return feat


# ---------------------------------------------------------------- 检查点
def load_checkpoint(path):
    if not os.path.exists(path):
        return {"done_pdb": [], "failed_pdb": {}, "skipped_pdb": {},
                "started_at": None, "last_update": None, "n_rows_written": 0}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            cp = json.load(fh)
        for k, v in (("done_pdb", []), ("failed_pdb", {}), ("skipped_pdb", {})):
            cp.setdefault(k, v)
        return cp
    except Exception as exc:                                    # noqa: BLE001
        logger.error("检查点损坏 (%s); 将从 CSV 重建进度", exc)
        return {"done_pdb": [], "failed_pdb": {}, "skipped_pdb": {},
                "started_at": None, "last_update": None, "n_rows_written": 0,
                "corrupted": True}


def save_checkpoint(path, cp):
    cp["last_update"] = time.strftime("%Y-%m-%d %H:%M:%S")
    tmp = path + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(cp, fh, ensure_ascii=False, indent=1)
        os.replace(tmp, path)
    except Exception as exc:                                    # noqa: BLE001
        logger.error("检查点写出失败: %s", exc)


def load_existing_csv(path):
    """读回已写结果并按 interface_key 去重 (保留最后一条)。"""
    if not os.path.exists(path) or os.path.getsize(path) < 10:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path)
        if "interface_key" in df.columns:
            n0 = len(df)
            df = df.drop_duplicates(subset=["interface_key"], keep="last")
            if len(df) != n0:
                logger.warning("既有 CSV 去重: %d -> %d 行", n0, len(df))
        return df
    except Exception as exc:                                    # noqa: BLE001
        logger.error("读取既有 CSV 失败 (%s); 将重建", exc)
        return pd.DataFrame()


def append_rows(path, rows):
    """追加写出; 返回写入行数。列集合变化时按并集重写, 保证可续算。

    注意: 先对每批做 V0.2 的派生特征计算 (derive_normalized_features), 保证
    全量表与 MVP 特征表的列集合完全一致, 否则训练阶段特征对不齐。
    """
    if not rows:
        return 0
    try:
        df = _FEAT02.derive_normalized_features(pd.DataFrame(rows))
    except Exception as exc:                                    # noqa: BLE001
        logger.error("派生特征计算失败 (%s); 该批按原始列写出", exc)
        df = pd.DataFrame(rows)
    fresh = not os.path.exists(path) or os.path.getsize(path) < 10
    try:
        if fresh:
            df.to_csv(path, index=False, encoding="utf-8")
            return len(df)
        cols = pd.read_csv(path, nrows=0).columns.tolist()
        missing = [c for c in cols if c not in df.columns]
        extra = [c for c in df.columns if c not in cols]
        if missing or extra:
            logger.warning("列集合变化 (缺 %s / 多 %s); 合并重写该批", missing, extra)
            old = pd.read_csv(path)
            merged = pd.concat([old, df], ignore_index=True)
            if "interface_key" in merged.columns:
                merged = merged.drop_duplicates(subset=["interface_key"], keep="last")
            merged.to_csv(path, index=False, encoding="utf-8")
            return len(df)
        df[cols].to_csv(path, mode="a", header=False, index=False, encoding="utf-8")
        return len(df)
    except Exception as exc:                                    # noqa: BLE001
        logger.error("追加落盘失败: %s", exc)
        return 0


def summarize(out_csv, ckpt_path, cp, pdb_ids, todo_rows, t0):
    df = load_existing_csv(out_csv)
    logger.info("=" * 78)
    print("\n" + "=" * 80)
    print("全量特征提取汇总")
    print("=" * 80)
    print("输出文件      : %s" % out_csv)
    print("检查点        : %s" % ckpt_path)
    print("目标结构数    : %d" % len(pdb_ids))
    print("已完成结构    : %d" % len(cp.get("done_pdb", [])))
    print("失败结构      : %d" % len(cp.get("failed_pdb", {})))
    print("跳过结构      : %d  %s" % (len(cp.get("skipped_pdb", {})),
                                     list(cp.get("skipped_pdb", {}).items())[:5]))
    print("特征表行数    : %d" % len(df))
    if not df.empty:
        rate = 100.0 * len(df) / max(todo_rows, 1)
        print("行级成功率    : %.1f%%  (%d / %d 行)" % (rate, len(df), todo_rows))
        if "pdb_id" in df.columns:
            print("覆盖唯一结构  : %d" % df["pdb_id"].nunique())
        if "pair_sasa_complete" in df.columns:
            print("链路对特征完整: %d / %d" % (int((df["pair_sasa_complete"] == 1).sum()),
                                              len(df)))
        if "self_interface" in df.columns:
            print("自界面标记行  : %d" % int(df["self_interface"].sum()))
    print("耗时          : %.1f 分钟" % ((time.time() - t0) / 60))
    print("=" * 80)
    logger.info("汇总: 行数=%d, 完成结构=%d, 失败=%d, 跳过=%d, 耗时=%.1fmin",
                len(df), len(cp.get("done_pdb", [])), len(cp.get("failed_pdb", {})),
                len(cp.get("skipped_pdb", {})), (time.time() - t0) / 60)


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0,
                    help="只处理前 N 条结构 (验证用; 产物自动加 _test 后缀)")
    ap.add_argument("--include-mvp", action="store_true",
                    help="包含 is_mvp 的 400 行 (默认排除)")
    ap.add_argument("--tag", type=str, default="",
                    help="产物后缀 (默认: 有 --limit 时自动加 _test)")
    ap.add_argument("--pair-sasa", type=str, default="auto",
                    help="逐链对 ΔSASA 策略: auto(单链路对精确直填, 多链对才算) / "
                         "整数上限 (如 8) / 0(全部跳过)")
    ap.add_argument("--exit-every", type=int, default=EXIT_EVERY,
                    help="每处理多少条结构落盘一次 (默认 100)")
    args = ap.parse_args()

    max_pair = -1 if args.pair_sasa == "auto" else int(args.pair_sasa)
    exit_every = max(1, args.exit_every)

    tag = args.tag or ("_test" if args.limit else "")
    out_csv = OUT_CSV.replace(".csv", "%s.csv" % tag)
    ckpt_path = CHECKPOINT.replace(".json", "%s.json" % tag)

    t0 = time.time()
    logger.info("=" * 78)
    logger.info("全量特征提取启动 | 输出=%s | 检查点=%s | pair-sasa=%s | 落盘间隔=%d",
                out_csv, ckpt_path, args.pair_sasa, exit_every)
    feat02 = _FEAT02

    from Bio.PDB import MMCIFParser, PDBParser

    bench = pd.read_csv(BENCH, dtype={"pdb_id": str})
    todo = bench if args.include_mvp else bench[bench["is_mvp"] != True]   # noqa: E712
    logger.info("benchmark 总行 %d -> 待处理 %d 行", len(bench), len(todo))

    by_pdb = {}
    for rec in todo.to_dict("records"):
        by_pdb.setdefault(rec["pdb_id"], []).append(rec)
    pdb_ids = sorted(by_pdb)
    if args.limit:
        pdb_ids = pdb_ids[:args.limit]
    todo_rows = sum(len(by_pdb[p]) for p in pdb_ids)
    logger.info("唯一结构数: %d (对应 %d 行)", len(pdb_ids), todo_rows)

    cp = load_checkpoint(ckpt_path)
    existing = load_existing_csv(out_csv)
    done = set(cp["done_pdb"])
    failed = dict(cp["failed_pdb"])
    skipped = dict(cp["skipped_pdb"])

    if not existing.empty and "pdb_id" in existing.columns:
        csv_pdbs = set(existing["pdb_id"].astype(str).unique())
        newly = csv_pdbs - done
        if newly:
            logger.warning("CSV 中有 %d 个结构不在检查点, 补记为已完成", len(newly))
            done |= newly
    if cp.get("started_at") is None:
        cp["started_at"] = time.strftime("%Y-%m-%d %H:%M:%S")

    remaining = [p for p in pdb_ids
                 if p not in done and p not in failed and p not in skipped]
    logger.info("已完成=%d 已失败=%d 已跳过=%d 剩余=%d | 既有 CSV 行数=%d",
                len(done), len(failed), len(skipped), len(remaining), len(existing))
    if not remaining:
        logger.info("无剩余结构, 直接汇总")
        cp["done_pdb"] = sorted(done)
        save_checkpoint(ckpt_path, cp)
        summarize(out_csv, ckpt_path, cp, pdb_ids, todo_rows, t0)
        print("全量特征提取完成，等待下一步指令")
        return

    parser_pdb, parser_cif = PDBParser(QUIET=True), MMCIFParser(QUIET=True)
    buf, n_ok, n_fail, n_skip = [], 0, 0, 0

    for i, pid in enumerate(remaining, start=1):
        path = resolve_struct_path(pid)
        if path is None:
            skipped[pid] = "no_structure_file"
            n_skip += 1
            continue

        recs = by_pdb[pid]
        try:
            parser = parser_cif if path.lower().endswith(".cif") else parser_pdb
            structure = parser.get_structure(pid, path)
            model = structure[0]
            chain_ids = [c.id for c in model]
            ss_map = feat02.build_ss_map(path)
        except Exception as exc:                                # noqa: BLE001
            failed[pid] = "parse_error: %s" % str(exc)[:80]
            n_fail += 1
            logger.warning("解析失败(跳过, 不再重试) %s: %s", pid, str(exc)[:70])
            gc.collect()
            continue

        for rec in recs:
            lig = [c.strip() for c in str(rec["ligand_chains"]).split(",") if c.strip()]
            rcp = [c.strip() for c in str(rec["receptor_chains"]).split(",") if c.strip()]
            try:
                lig_p = [c for c in lig if c in chain_ids]
                rcp_p = [c for c in rcp if c in chain_ids]
                if not lig_p or not rcp_p:
                    raise ValueError("链组完全缺失 (lig=%s rec=%s)" % (lig, rcp))
                n_missing = len(lig) + len(rcp) - len(lig_p) - len(rcp_p)

                xa, ea, na, ra, rna, ssa = feat02.extract_chain_data(model, lig_p, ss_map)
                xb, eb, nb, rb, rnb, ssb = feat02.extract_chain_data(model, rcp_p, ss_map)
                if len(xa) == 0 or len(xb) == 0:
                    raise ValueError("链组无有效蛋白原子")

                feat = {"pdb_id": pid,
                        "interface_key": rec["interface_key"],
                        "ligand_chains": rec["ligand_chains"],
                        "receptor_chains": rec["receptor_chains"],
                        "real_dG": float(rec["real_dG"]),
                        "split": rec.get("split", ""),
                        "split_v0.1": rec.get("split_v0.1", rec.get("split", "")),
                        "family_id": rec.get("family_id", ""),
                        "subgroup": rec.get("subgroup", ""),
                        "source_dataset": rec.get("source_dataset", ""),
                        "n_missing_chain": int(n_missing),
                        "chain_map_ok": int(n_missing == 0),
                        "self_interface": int(set(lig_p) == set(rcp_p))}

                gf = feat02.group_features(xa, ea, na, ra, rna, ssa,
                                           xb, eb, nb, rb, rnb, ssb)
                pair_contacts = gf.pop("_pair_contacts", {})
                pair_hbonds = gf.pop("_pair_hbonds", {})
                feat.update(gf)

                with tempfile.TemporaryDirectory() as td:
                    ncp = int(feat.get("n_chainpair_contacts", 1) or 1)
                    if ncp <= 1:
                        feat.update(feat02.sasa_features(
                            model, lig_p, rcp_p, td,
                            pair_contacts={}, pair_hbonds={}, max_pair_sasa=0))
                        fill_single_pair_features(feat)
                    elif max_pair == 0:
                        feat.update(feat02.sasa_features(
                            model, lig_p, rcp_p, td,
                            pair_contacts={}, pair_hbonds={}, max_pair_sasa=0))
                        nan_pair_features(feat)
                    else:
                        cap = MAX_PAIR_SASA if max_pair == -1 else max_pair
                        feat.update(feat02.sasa_features(
                            model, lig_p, rcp_p, td,
                            pair_contacts=pair_contacts, pair_hbonds=pair_hbonds,
                            max_pair_sasa=cap))

                feat["predicted_dG"] = np.nan
                feat["prodigy_ic"] = np.nan
                feat["feature_version"] = "v1.0-full-freesasa-Biopython"
                buf.append(feat)
                n_ok += 1
            except Exception as exc:                            # noqa: BLE001
                logger.warning("特征失败(跳过) %s %s|%s: %s",
                               pid, rec["ligand_chains"], rec["receptor_chains"],
                               str(exc)[:80])
                failed["%s|%s|%s" % (pid, rec["ligand_chains"],
                                     rec["receptor_chains"])] = \
                    "feature_error: %s" % str(exc)[:60]
                n_fail += 1

        del structure, model
        done.add(pid)

        if i % exit_every == 0 or i == len(remaining):
            written = append_rows(out_csv, buf)
            cp["done_pdb"] = sorted(done)
            cp["failed_pdb"] = failed
            cp["skipped_pdb"] = skipped
            cp["n_rows_written"] = int(cp.get("n_rows_written", 0)) + written
            save_checkpoint(ckpt_path, cp)
            buf = []
            gc.collect()
            elapsed = time.time() - t0
            rate = elapsed / max(i, 1)
            eta = rate * (len(remaining) - i) / 60
            logger.info("进度 %d/%d | 成功=%d 失败=%d 跳过=%d | 本批写入=%d | "
                        "内存=%.0fMB | 已用=%.1fmin | ETA=%.1fmin",
                        i, len(remaining), n_ok, n_fail, n_skip, written,
                        get_rss_mb(), elapsed / 60, eta)

    cp["done_pdb"] = sorted(done)
    cp["failed_pdb"] = failed
    cp["skipped_pdb"] = skipped
    save_checkpoint(ckpt_path, cp)
    gc.collect()
    summarize(out_csv, ckpt_path, cp, pdb_ids, todo_rows, t0)
    print("全量特征提取完成，等待下一步指令")


if __name__ == "__main__":
    main()
