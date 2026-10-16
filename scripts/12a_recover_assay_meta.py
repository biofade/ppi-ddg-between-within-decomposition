# -*- coding: utf-8 -*-
"""12a_recover_assay_meta.py -- 从源表恢复测量条件元数据 (方法 / 温度 / 亚组).

背景 (为什么需要这个脚本)
------------------------
`01_clean_data.py` 在把 12062 条原始读数按 `interface_key` 合并成 4762 条时,
对 `affinity_method` 与 `temperature_k` 只取了 `"first"` (见 01_clean_data.py:241-242)。
结果是 benchmark.csv 里这两列的覆盖率只有 22.6% / 14.0% —— 信息在聚合时被丢弃了,
而源表本身覆盖率高得多。本脚本把源表按同一 `interface_key` 重新聚合, 给出:

* `method_majority`   同一界面下出现最多的测量方法 (票数最多者, 平票取字典序最小)
* `method_set`        该界面出现过的全部方法 (``|`` 连接)
* `method_n_votes`    参与投票的读数条数
* `method_conflict`   1 = 同一界面存在 >1 种方法 (测量异质性标签)
* `temperature_median_k`  该界面全部温度读数的中位数 (K)
* `temperature_min_k` / `temperature_max_k`
* `temperature_n`     有温度读数的条数 (0 = 无温度信息)
* `temperature_spread_k`  max - min (0 = 单温度, >0 = 温度混杂)
* `temp_imputed`      1 = 无实验温度, 主流程按 298.15 K 兜底
* `temp_assumed_flag` 源表中形如 ``298(assumed)`` 的占位温度比例
* `subgroup_any`      源表 Subgroup 非空则取之

本脚本**只读** PPB-Affinity.xlsx 与 benchmark.csv, 只写
`data/v2_assay_meta.csv` 与 `logs/v2/12a_recover_assay_meta.log`,
不改动任何 V1.1 产物。

纪律: 纯 CPU; 无网络; 固定随机种子 (本脚本无随机过程); 逐行容错不中断。
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
DATA_DIR = ROOT / "data"
FEAT_DIR = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2"

XLSX = DATA_DIR / "PPB-Affinity.xlsx"
BENCHMARK = DATA_DIR / "benchmark.csv"
OUT_CSV = DATA_DIR / "v2_assay_meta.csv"
LOG_PATH = LOG_DIR / "12a_recover_assay_meta.log"

SEED = 42
IMPUTE_TEMP_K = 298.15

_TEMP_NUM = re.compile(r"(-?\d+(?:\.\d+)?)")

logger = logging.getLogger("recover_assay_meta")


# --------------------------------------------------------------------- helpers
def _setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)


def norm_chain(value: Any) -> str:
    """把 ``'A, B'`` / ``'a;b'`` / ``'A:B'`` 规范成 ``'A,B'`` (大写/去重/排序)。

    与 ``01_clean_data.py::norm_chain`` 保持完全一致, 否则 interface_key 对不上。
    """
    if value is None:
        return ""
    txt = str(value).replace(";", ",").replace("|", ",")
    parts = [p.strip().upper() for p in txt.split(",") if p.strip()]
    parts = [p.split(":")[0].strip() for p in parts]
    parts = [p for p in parts if p]
    return ",".join(sorted(dict.fromkeys(parts)))


def parse_temperature(raw: Any) -> tuple[float, bool]:
    """``'298(assumed)'`` -> ``(298.0, True)``; ``'ITC'``/``''`` -> ``(nan, False)``."""
    if raw is None:
        return float("nan"), False
    txt = str(raw).strip()
    if not txt or txt.lower() in ("nan", "none"):
        return float("nan"), False
    m = _TEMP_NUM.search(txt)
    if not m:
        return float("nan"), False
    try:
        val = float(m.group(1))
    except ValueError:
        return float("nan"), False
    if not np.isfinite(val) or val <= 0 or val > 400:
        return float("nan"), False
    return val, "assumed" in txt.lower()


def norm_method(raw: Any) -> str:
    if raw is None:
        return ""
    txt = str(raw).strip()
    if not txt or txt.lower() in ("nan", "none"):
        return ""
    return txt


# --------------------------------------------------------------------- main
def build_interface_key(frame: pd.DataFrame, colmap: dict[str, str]) -> pd.Series:
    pdb = frame[colmap["pdb"]].fillna("").astype(str).str.strip().str.upper()
    lig = frame[colmap["lig"]].map(norm_chain)
    rec = frame[colmap["rec"]].map(norm_chain)
    return pdb + "|" + lig + "|" + rec


def recover(source: pd.DataFrame, benchmark: pd.DataFrame) -> pd.DataFrame:
    """按 interface_key 聚合源表的测量条件。"""
    colmap = {
        "pdb": "PDB",
        "lig": "Ligand Chains",
        "rec": "Receptor Chains",
        "method": "Affinity Method",
        "temp": "Temperature(K)",
        "subgroup": "Subgroup",
    }
    missing = [v for v in colmap.values() if v not in source.columns]
    if missing:
        raise RuntimeError(f"源表缺少必需列: {missing}")

    src = pd.DataFrame(
        {
            "interface_key": build_interface_key(source, colmap),
            "method": source[colmap["method"]].map(norm_method),
            "subgroup": source[colmap["subgroup"]].map(norm_method),
        }
    )
    temps = source[colmap["temp"]].map(parse_temperature)
    src["temperature_k"] = [t[0] for t in temps]
    src["temp_assumed"] = [t[1] for t in temps]

    n_src = len(src)
    src = src[src["interface_key"].str.count(r"\|") == 2]
    src = src[src["interface_key"].str.split("|").str[0].ne("")]
    logger.info("源表读数 %d 条, 可用于聚合 %d 条 (唯一界面 %d)",
                n_src, len(src), src["interface_key"].nunique())

    records: list[dict[str, Any]] = []
    for ik, grp in src.groupby("interface_key", sort=True):
        votes = Counter(m for m in grp["method"] if m)
        if votes:
            top = max(votes.items(), key=lambda kv: (kv[1], -ord(kv[0][0]) if kv[0] else 0))
            # 平票时按方法名字典序取最小, 保证确定性
            best_n = top[1]
            method_majority = sorted(k for k, v in votes.items() if v == best_n)[0]
        else:
            method_majority = ""
        tvals = grp["temperature_k"].dropna().to_numpy(dtype=float)
        n_temp = int(tvals.size)
        if n_temp:
            t_med = float(np.median(tvals))
            t_min, t_max = float(tvals.min()), float(tvals.max())
            spread = t_max - t_min
        else:
            t_med = t_min = t_max = float("nan")
            spread = float("nan")
        sub = [s for s in grp["subgroup"] if s]
        records.append(
            {
                "interface_key": ik,
                "method_majority": method_majority,
                "method_set": "|".join(sorted(votes)) if votes else "",
                "method_n_votes": int(sum(votes.values())),
                "method_conflict": int(len(votes) > 1),
                "affinity_method_recovered": method_majority,
                "temperature_median_k": t_med,
                "temperature_min_k": t_min,
                "temperature_max_k": t_max,
                "temperature_n": n_temp,
                "temperature_spread_k": spread,
                "temp_imputed": int(n_temp == 0),
                "temp_assumed_frac": (
                    float(np.mean(grp.loc[grp["temperature_k"].notna(), "temp_assumed"]))
                    if n_temp
                    else float("nan")
                ),
                "subgroup_any": sub[0] if sub else "",
            }
        )
    meta = pd.DataFrame(records)

    merged = benchmark[["interface_key"]].merge(meta, on="interface_key", how="left")
    n_missing = int(merged["method_majority"].isna().sum())
    if n_missing:
        logger.warning("%d 条 benchmark 界面在源表中无聚合结果", n_missing)

    # 温度兜底: 无实验温度 -> 298.15 K (并在主流程标注 temp_imputed)
    merged["temperature_k_used"] = merged["temperature_median_k"].fillna(IMPUTE_TEMP_K)
    return merged


def summarise(meta: pd.DataFrame) -> list[str]:
    lines = ["", "=" * 78, "测量条件恢复汇总", "=" * 78]
    n = len(meta)
    lines.append(f"benchmark 界面数            : {n}")
    cov_m = float((meta['method_majority'].fillna('') != '').mean())
    cov_t = float((meta["temperature_n"].fillna(0) > 0).mean())
    lines.append(f"方法覆盖 (源表恢复)          : {cov_m * 100:.1f}%  ({int(cov_m * n)} 条)")
    lines.append(f"温度覆盖 (源表恢复)          : {cov_t * 100:.1f}%  ({int(cov_t * n)} 条)")
    lines.append(f"温度需兜底 298.15 K          : {int(meta['temp_imputed'].sum())} 条")
    lines.append(f"同界面方法冲突 (测量异质性)   : {int(meta['method_conflict'].sum())} 条")
    spread = meta["temperature_spread_k"].dropna()
    lines.append(f"同界面温度跨度 >5 K          : {int((spread > 5).sum())} 条")
    lines.append("")
    lines.append("方法分布 (众数投票):")
    vc = meta["method_majority"].replace("", "(无方法信息)").value_counts()
    for k, v in vc.items():
        lines.append(f"  {k:<20s} {v:6d}  ({v / n * 100:5.1f}%)")
    lines.append("")
    lines.append("温度区间分布 (中位数温度):")
    tb = pd.cut(
        meta["temperature_median_k"],
        bins=[0, 290, 295, 299, 303, 310, 400],
        labels=["<290", "290-295", "295-299", "299-303", "303-310", ">310"],
    )
    for k, v in tb.value_counts().sort_index().items():
        lines.append(f"  {str(k):<12s} {v:6d}  ({v / n * 100:5.1f}%)")
    lines.append("=" * 78)
    return lines


def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="从 PPB-Affinity 源表恢复测量条件元数据")
    ap.add_argument("--xlsx", type=Path, default=XLSX, help="源工作簿路径")
    ap.add_argument("--benchmark", type=Path, default=BENCHMARK, help="benchmark.csv 路径")
    ap.add_argument("--out", type=Path, default=OUT_CSV, help="输出 CSV 路径")
    ap.add_argument("--seed", type=int, default=SEED, help="随机种子 (本脚本无随机过程, 仅记录)")
    ap.add_argument("--dry-run", action="store_true", help="只打印汇总, 不写文件")
    args = ap.parse_args()

    _setup_logging(LOG_PATH)
    logger.info("seed=%d (本脚本为确定性聚合, 无随机过程)", args.seed)

    if not args.xlsx.exists():
        logger.error("源表不存在: %s", args.xlsx)
        return 2
    if not args.benchmark.exists():
        logger.error("benchmark 不存在: %s", args.benchmark)
        return 2

    logger.info("读取源表: %s", args.xlsx)
    source = pd.read_excel(args.xlsx, sheet_name="Sheet1", dtype=str)
    logger.info("源表 %d 行 x %d 列", *source.shape)
    benchmark = pd.read_csv(args.benchmark, dtype={"pdb_id": str})
    logger.info("benchmark %d 行", len(benchmark))

    meta = recover(source, benchmark)
    for line in summarise(meta):
        logger.info("%s", line)

    if args.dry_run:
        logger.info("--dry-run: 不写文件")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + ".tmp")
    meta.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out)
    logger.info("写出 %s (%d 行 x %d 列)", args.out, *meta.shape)
    logger.info("耗时 %.1fs", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
