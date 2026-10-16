# -*- coding: utf-8 -*-
"""scripts/v5/01_data_manifest.py -- V5.0 PHASE 0: 数据哈希清单（只读）.

对 V4 及更早的关键产物记录: path / size / sha256 / row count / col count / mtime。
不修改任何被扫描文件（只读打开）。

输出: experiments/v5_0/manifests/data_manifest.csv
      experiments/v5_0/manifests/data_manifest.md（可读版）
"""

from __future__ import annotations

import csv
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
MAN = ROOT / "experiments" / "v5_0" / "manifests"
MAN.mkdir(parents=True, exist_ok=True)
OUT_CSV = MAN / "data_manifest.csv"
OUT_MD = MAN / "data_manifest.md"

# 需要登记的关键产物（缺失则记 MISSING，不报错）
KEY_FILES: list[tuple[str, str]] = [
    # --- 原始数据 ---
    ("data/ppb_affinity_proteinea/filtered.csv", "原始 PPB-Affinity filtered（划分来源）"),
    ("data/ppb_affinity_proteinea/raw.csv", "原始 raw"),
    ("data/v2.4/skempi_v2_cleaned.csv", "SKEMPI v2.0 清洗表"),
    # --- 特征表 ---
    ("features/v2.2_proteinea_features.csv", "V2.2 物理+元数据特征表（V4/V5 主表）"),
    ("features/features_v1.1.csv", "V1.1 特征表（历史口径，标签口径不同）"),
    # --- 缓存 ---
    ("features/v2_esm2_cache.npz", "pooled ESM-2 缓存（8207 行 × 640 × 2 链）"),
    ("features/v2_esm2_perres_v2.npz", "per-residue 原始缓存（1.9 GB）"),
    ("features/v2.5_perres_cache.npz", "per-residue 紧凑缓存（582 MB）"),
    # --- V2.5 衍生特征 ---
    ("features/v2.5_energy_features_norm.csv", "V2.5 归一化能量特征"),
    ("features/v2.5_conservation_features.csv", "V2.5 保守性代理"),
    ("features/v2.5_coevolution_features.csv", "V2.5 共变代理"),
    ("features/v2.5_thermo_features.csv", "V2.5 热力学分解"),
    ("features/v2.5_hotspot_weights.npz", "V2.5 热点权重"),
    # --- V4 配对 / 划分 ---
    ("data/v4.0/ddg_pairs.csv", "V4 ΔΔG 配对表（含 forward+reverse）"),
    ("data/v4.0/ddg_groups.json", "V4 每复合物 WT/mutant 结构"),
    ("data/v4.0/ddg_split.json", "V4 三重隔离 split（V5 起点）"),
    ("data/v4.0/ddg_split_stratified.json", "V4 诊断用均值对齐变体（非主划分）"),
    ("data/v4.0/ddg_benchmark/ddg_benchmark_pairs.csv", "V4 benchmark 配对（7526 行）"),
    ("data/v4.0/ddg_benchmark/split.json", "V4 benchmark split"),
    # --- V4 结果 ---
    ("experiments/v4.0/innovation_results.csv", "V4 创新 1 结果"),
    ("experiments/v4.0/innovation_predictions.npz", "V4 创新 1 逐样本预测"),
    ("experiments/v4.0/innovation2_results.csv", "V4 创新 2 结果"),
    ("experiments/v4.0/innovation2_predictions.npz", "V4 创新 2 逐样本预测 + OOF ΔG 辅助"),
    ("experiments/v4.0/baselines/ddg_baseline_results.csv", "V4 基线结果（含早停 bug 版）"),
    ("experiments/v4.0/baselines/ddg_baseline_predictions.npz", "V4 基线预测"),
    ("experiments/v4.0/baselines/ddg_lgb_curve.csv", "V4 LightGBM 学习曲线"),
    ("experiments/v4.0/baselines/ddg_split_variants.csv", "V4 划分变体诊断"),
    # --- 历史报告（交叉核对用） ---
    ("features/v4.0_final_report.md", "V4 最终报告"),
    ("features/v4.0_all_results.csv", "V4 全结果汇总"),
    ("features/v3.1_paper_summaries.md", "V3.1 文献结构化（带页码）"),
]


def sha256(p: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def describe(p: Path) -> tuple[object, object, str]:
    """返回 (n_rows, n_cols, note)。只读。"""
    suf = p.suffix.lower()
    if suf == ".csv":
        try:
            d = pd.read_csv(p, low_memory=False, nrows=5)
            ncol = d.shape[1]
            # 行数：快速数换行（避免大内存）
            with open(p, "rb") as fh:
                n = sum(1 for _ in fh)
            return n - 1, ncol, ""
        except Exception as exc:  # noqa: BLE001
            return "", "", f"read_failed:{type(exc).__name__}"
    if suf == ".npz":
        try:
            z = np.load(p, allow_pickle=True)
            keys = list(z.files)
            shapes = []
            for k in keys[:8]:
                try:
                    shapes.append(f"{k}{tuple(z[k].shape)}")
                except Exception:  # noqa: BLE001
                    shapes.append(f"{k}<?>")
            nmax = ""
            try:
                nmax = int(z[keys[0]].shape[0])
            except Exception:  # noqa: BLE001
                pass
            return nmax, len(keys), ";".join(shapes)[:300]
        except Exception as exc:  # noqa: BLE001
            return "", "", f"npz_read_failed:{type(exc).__name__}"
    if suf == ".json":
        try:
            import json
            j = json.loads(p.read_text(encoding="utf-8"))
            if isinstance(j, dict):
                return len(j), len(j), "dict_keys=" + ",".join(list(j)[:6])
            if isinstance(j, list):
                return len(j), "", "list"
        except Exception as exc:  # noqa: BLE001
            return "", "", f"json_read_failed:{type(exc).__name__}"
    if suf == ".md":
        return "", "", f"{p.stat().st_size} bytes"
    return "", "", ""


def main() -> int:
    t0 = time.time()
    rows = []
    print(f"{'exists':6s} {'size':>12s}  {'rows':>8s} {'cols':>6s}  path")
    for rel, desc in KEY_FILES:
        p = ROOT / rel
        rec = {"path": rel, "description": desc, "exists": int(p.exists()),
               "size_bytes": "", "sha256": "", "n_rows": "", "n_cols": "",
               "mtime": "", "note": ""}
        if p.exists():
            st = p.stat()
            rec["size_bytes"] = st.st_size
            rec["mtime"] = time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(st.st_mtime))
            rec["sha256"] = sha256(p)
            nr, nc, note = describe(p)
            rec["n_rows"], rec["n_cols"], rec["note"] = nr, nc, note
        else:
            rec["note"] = "MISSING"
        rows.append(rec)
        print(f"{rec['exists']:^6d} {str(rec['size_bytes']):>12s}  "
              f"{str(rec['n_rows']):>8s} {str(rec['n_cols']):>6s}  {rel}"
              + (f"   [{rec['note'][:60]}]" if rec["note"] else ""))

    with open(OUT_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_ex = sum(r["exists"] for r in rows)
    L = ["# V5.0 PHASE 0 — 数据哈希清单", "",
         f"- 生成: {time.strftime('%Y-%m-%d %H:%M:%S')}",
         f"- 登记文件: **{len(rows)}**，存在 **{n_ex}**，缺失 **{len(rows) - n_ex}**",
         f"- 说明: 只读扫描，未修改任何被登记文件。", "",
         "| exists | size (B) | rows | cols | mtime | sha256(前16) | path | note |",
         "|---:|---:|---:|---:|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['exists']} | {r['size_bytes']} | {r['n_rows']} | {r['n_cols']} | "
                 f"{r['mtime']} | `{str(r['sha256'])[:16]}` | `{r['path']}` | "
                 f"{r['note'][:70]} |")
    L += ["", "### 缺失文件（如实登记，不猜测内容）", ""]
    miss = [r for r in rows if not r["exists"]]
    if miss:
        for r in miss:
            L.append(f"- `{r['path']}` — {r['description']}")
    else:
        L.append("- 无")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")

    print(f"\n存在 {n_ex}/{len(rows)}；耗时 {time.time() - t0:.0f}s")
    print(f"写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
