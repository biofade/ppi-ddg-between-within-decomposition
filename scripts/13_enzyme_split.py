# -*- coding: utf-8 -*-
"""13_enzyme_split.py -- 阶段 4: 酶数据接入与三层隔离 (EC / CATH / 序列).

背景与动机
----------
学长建议「往酶方向靠」的技术理由是: **酶数据天然有 group 结构**
(同一酶 × 多底物, 或同一酶 × 多突变), 组内 ΔG 各不相同, 正是 Learning-to-Rank
需要的形态。阶段 3 已实测证明当前 PPI 数据的 group 里 61.7% 的组内 ΔG 完全相同
(同一测量的链组合冗余), 排序梯度恒为 0。

本脚本做两件事
--------------
**A. 三层隔离键构建器** (目标 schema, 可直接复用)

.. code-block:: text

    source, pdb_id, chain, sequence, ec_number, cath_id,
    ligand_id, label_type, value, unit, temperature_K, ph,
    assay_method, real_dG, group_id, family_id

三层隔离优先级: **EC > CATH > 序列(40% 同一性)**, 最后附加 ``PDB_<id>`` 指纹兜底。
硬校验: ``isolation_key`` 与 ``pdb_id`` 跨侧数量必须为 0, 否则抛错停止。

**B. 本地可得数据的诚实审计**

提示词列出的数据源 (SABIO-RK / PDBbind 酶子集 / SKEMPI 2.0 / BRENDA / BindingDB)
都需要**联网下载**。本脚本**不伪造**这些数据, 而是:

1. 检查本地是否存在 EC / CATH 元数据文件 (``data/ec_annotations.csv`` /
   ``data/cath_annotations.csv``), 存在则启用对应层;
2. 不存在则**明确标记该层为「未接入」**, 只启用序列层, 并在报告中给出
   可执行的下游接入指引 (要跑哪个脚本、要拿到什么文件);
3. 从本地 PDB **HEADER 记录**提取大分子类别 (HYDROLASE / TRANSFERASE /
   OXIDOREDUCTASE / LYASE / ISOMERASE / LIGASE), 作为"酶样"的**弱证据**,
   并显式声明这只是 PDB 沉积注释, **不是 EC 号, 不能替代 EC 层**。

序列层是本脚本的核心可用部分: 用 40% 序列同一性 (k-mer Jaccard 代理) 贪心聚类,
既是隔离键, 也顺便回答"本地数据里有多少真正多成员的酶样 group"。

纪律: 纯 CPU; 无网络; seed=42; 逐结构容错; 不改动 V1.1 产物。
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT_DIR = ROOT / "features"
DATA_DIR = ROOT / "data"
PDB_DIR = DATA_DIR / "pdb"
LOG_DIR = ROOT / "logs" / "v2"

FEATURES_CSV = FEAT_DIR / "features_v1.1.csv"
BENCHMARK_CSV = DATA_DIR / "benchmark.csv"
ASSAY_META_CSV = DATA_DIR / "v2_assay_meta.csv"
SEQ_CACHE = DATA_DIR / "pdb_seqs.json"
EC_FILE = DATA_DIR / "ec_annotations.csv"
CATH_FILE = DATA_DIR / "cath_annotations.csv"

OUT_DATASET = FEAT_DIR / "v2_enzyme_dataset.csv"
OUT_SPLIT = FEAT_DIR / "v2_enzyme_split.csv"
REPORT_MD = FEAT_DIR / "v2_enzyme_dataset_report.md"
JSON_OUT = FEAT_DIR / "v2_enzyme_dataset_summary.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
IDENTITY_THRESHOLD = 0.40        # 提示词要求: 40% 序列同一性
KMER_K = 3                       # k-mer Jaccard 代理 (与 01/09 同口径)
PDB_FINGERPRINT_PREFIX = "PDB_"

# PDB HEADER 的大分子类别 -> 对应 EC 一级类 (仅作"酶样"弱证据, 不是 EC 号)
EC_CLASS_FROM_HEADER = {
    "OXIDOREDUCTASE": "EC_1",
    "TRANSFERASE": "EC_2",
    "HYDROLASE": "EC_3",
    "LYASE": "EC_4",
    "ISOMERASE": "EC_5",
    "LIGASE": "EC_6",
    "TRANSLOCASE": "EC_7",
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}

logger = logging.getLogger("enzyme_split")


# ===================================================================== 隔离键工具
def ec_layer(ec: str | None) -> str:
    """EC 层: 取前 2 段 (``1.1.1.1`` -> ``EC_1.1``); 缺失返回空串。"""
    if not ec:
        return ""
    txt = str(ec).strip().replace("EC:", "").replace("EC_", "")
    if not txt or txt.lower() in ("nan", "none", "-"):
        return ""
    parts = [p for p in txt.split(".") if p.strip()]
    if len(parts) < 2:
        return f"EC_{parts[0]}" if parts else ""
    return f"EC_{parts[0]}.{parts[1]}"


def cath_layer(cath: str | None) -> str:
    """CATH 层: 取前 3 段 (``1.10.8.10`` -> ``CATH_1.10.8``); 缺失返回空串。"""
    if not cath:
        return ""
    txt = str(cath).strip().replace("CATH_", "")
    if not txt or txt.lower() in ("nan", "none", "-"):
        return ""
    parts = [p for p in txt.split(".") if p.strip()]
    if len(parts) < 3:
        return f"CATH_{'.'.join(parts)}" if parts else ""
    return f"CATH_{'.'.join(parts[:3])}"


def sequence_cluster_key(seq: str, reps: list[tuple[set, str]],
                         threshold: float = IDENTITY_THRESHOLD,
                         k: int = KMER_K) -> tuple[str, int]:
    """40% 同一性贪心代表聚类 -> (簇 ID, 是否新建)。

    用 **k-mer 集合 Jaccard** 作为序列同一性的可复算代理 (与
    ``01_clean_data.py``/``09_family_split_full.py`` 同口径, 见 v1.1 报告残留局限第 3 条:
    "同源性为代理指标, 簇比单链接更细更保守")。

    为什么不在这里跑 CD-HIT/MMseqs2: 二进制度依赖会破坏"纯 CPU 可复现"纪律,
    且 3000+ 序列的 O(n²) 全比对在本机不可行。该近似的局限已在报告中显式声明。
    """
    if not seq:
        return "", -1
    ks = {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else {seq}
    for rep_ks, uid in reps:
        union = len(ks | rep_ks)
        if union and len(ks & rep_ks) / union >= threshold:
            return uid, 0
    uid = f"SEQ_{len(reps):05d}"
    reps.append((ks, uid))
    return uid, 1


def build_isolation_key(ec_key: str, cath_key: str, seq_key: str,
                        pdb_id: str) -> tuple[str, str]:
    """合并优先级 EC > CATH > 序列 > PDB 指纹。返回 (isolation_key, 层级来源)。"""
    if ec_key:
        return f"{ec_key}|{seq_key}" if seq_key else ec_key, "EC"
    if cath_key:
        return f"{cath_key}|{seq_key}" if seq_key else cath_key, "CATH"
    if seq_key:
        return seq_key, "sequence"
    return f"{PDB_FINGERPRINT_PREFIX}{pdb_id}", "pdb_fingerprint"


def hard_check_no_leak(df: pd.DataFrame, key_cols: list[str]) -> dict[str, int]:
    """硬校验: 每个隔离键 / pdb_id 绝不能跨 train-test。返回违规计数。"""
    out = {}
    for col in key_cols:
        if col not in df.columns or "split" not in df.columns:
            out[col] = -1
            continue
        sub = df[df[col].astype(str).str.len() > 0]
        if len(sub) == 0:
            out[col] = 0
            continue
        out[col] = int((sub.groupby(col)["split"].nunique() > 1).sum())
    return out


# ===================================================================== 结构解析
def header_class(pdb_path: Path) -> str:
    """从 PDB HEADER 记录提取大分子类别 (仅作"酶样"弱证据, 不是 EC 号)。"""
    try:
        with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("HEADER"):
                    txt = line[10:50].strip().upper()
                    for cls in EC_CLASS_FROM_HEADER:
                        if cls in txt:
                            return cls
                    return ""
                if line.startswith("ATOM"):
                    break
    except Exception:  # noqa: BLE001
        return ""
    return ""


def seq_from_structure(pdb_path: Path, chains: list[str]) -> dict[str, str]:
    """从 ATOM 记录提取链序列 (缓存链 ID 不匹配时的兜底)。"""
    want = set(chains)
    per: dict[str, list[tuple[int, str]]] = {c: [] for c in chains}
    seen: dict[str, set] = {c: set() for c in chains}
    try:
        with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line[:6] != "ATOM  ":
                    continue
                ch = line[21]
                if ch not in want:
                    continue
                aa = AA3_TO_1.get(line[17:20].strip().upper())
                if aa is None:
                    continue
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                key = (resseq, line[26])
                if key in seen[ch]:
                    continue
                seen[ch].add(key)
                per[ch].append((resseq, aa))
    except Exception:  # noqa: BLE001
        return {c: "" for c in chains}
    return {c: "".join(a for _r, a in sorted(v, key=lambda t: t[0])) for c, v in per.items()}


def get_chain_seq(pdb_id: str, chains: list[str], seqs: dict) -> str:
    cs = seqs.get(pdb_id, {}) or {}
    parts = [cs.get(c, {}).get("seq", "") for c in chains]
    if all(parts):
        return "".join(parts)
    p = PDB_DIR / f"{pdb_id}.pdb"
    if not p.exists():
        p = PDB_DIR / f"{pdb_id}.ent"
    if p.exists():
        fb = seq_from_structure(p, chains)
        parts2 = [fb.get(c, "") for c in chains]
        merged = [a if len(a) >= len(b) else b for a, b in zip(parts, parts2)]
        if any(merged):
            return "".join(merged)
    return "".join(p for p in parts if p)


# ===================================================================== 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 4: 酶数据三层隔离 (EC/CATH/序列)")
    ap.add_argument("--features", type=Path, default=FEATURES_CSV)
    ap.add_argument("--benchmark", type=Path, default=BENCHMARK_CSV)
    ap.add_argument("--assay-meta", type=Path, default=ASSAY_META_CSV)
    ap.add_argument("--ec-file", type=Path, default=EC_FILE)
    ap.add_argument("--cath-file", type=Path, default=CATH_FILE)
    ap.add_argument("--out-dataset", type=Path, default=OUT_DATASET)
    ap.add_argument("--out-split", type=Path, default=OUT_SPLIT)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--identity", type=float, default=IDENTITY_THRESHOLD)
    ap.add_argument("--test-frac", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "13_enzyme_split.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))

    emit("# Affinity Intelligence V2.0 — 阶段 4: 酶数据接入与三层隔离报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}")
    emit(f"- 序列同一性阈值: {args.identity:.2f} (k={KMER_K} k-mer Jaccard 代理)")
    emit(f"- 隔离优先级: **EC > CATH > 序列 > PDB 指纹**")

    # ---- 0. 数据源可用性审计 (诚实: 不联网就明说) -----------------------
    emit("")
    emit("## 0. 数据源可用性审计")
    emit("")
    emit("| 数据源 | 本地是否可得 | 用途 | 状态 |")
    emit("|---|---|---|---|")
    audit: dict[str, Any] = {}
    src_rows = [
        ("PPB-Affinity (PDBbind/SKEMPI/SAbDab/ATLAS 汇总)", BENCHMARK_CSV.exists(),
         "Kd -> ΔG 主标签", "已接入"),
        ("测量条件元数据 (方法/温度)", ASSAY_META_CSV.exists(),
         "标签净化 + 特征", "已接入" if ASSAY_META_CSV.exists() else "未接入"),
        ("EC 号注释", args.ec_file.exists(), "三层隔离第 1 层",
         "已接入" if args.ec_file.exists() else "**未接入 (需联网获取)**"),
        ("CATH 结构域分类", args.cath_file.exists(), "三层隔离第 2 层",
         "已接入" if args.cath_file.exists() else "**未接入 (需联网获取)**"),
        ("SABIO-RK (Km/kcat/Ki/Kd)", False, "酶动力学标签 + EC", "**未接入 (需联网)**"),
        ("BRENDA", False, "酶标签 (许可严格)", "**未接入 (需联网 + 许可)**"),
        ("BindingDB", False, "Ki/IC50 扩展", "**未接入 (需联网)**"),
        ("SKEMPI 2.0 突变 ΔΔG", False, "突变效应 (已在 PPB 内部分覆盖)", "**未接入 (需联网)**"),
        ("PDB 结构文件 (本地)", PDB_DIR.exists(), "序列 / HEADER 类别", "已接入"),
    ]
    for name, have, use, status in src_rows:
        emit(f"| {name} | {'是' if have else '否'} | {use} | {status} |")
    audit["sources"] = [{"name": n, "available": h, "use": u, "status": s}
                        for n, h, u, s in src_rows]
    emit("")
    emit("> **重要**: 本机离线。EC / CATH / SABIO-RK / BRENDA / BindingDB 层**均未接入**,")
    emit("> 报告中的任何 EC / CATH 统计都只会在对应文件存在时出现。**不伪造**这些数据。")
    emit("> 接入方式见第 6 节「下游接入指引」。")

    # ---- 1. 载入并用 PDB HEADER 做"酶样"弱分类 -------------------------
    if not args.features.exists():
        emit(f"!! 特征表不存在: {args.features}")
        return 2
    df = pd.read_csv(args.features, low_memory=False)
    if args.limit:
        df = df.head(args.limit).copy()
    emit("")
    emit("## 1. 本地结构与酶样分类 (弱证据)")
    emit("")
    emit(f"- 输入: {args.features.name} — {len(df)} 行 / {df['pdb_id'].nunique()} 结构")

    seqs: dict = {}
    if SEQ_CACHE.exists():
        seqs = json.loads(SEQ_CACHE.read_text(encoding="utf-8"))
    emit(f"- 序列缓存: {len(seqs)} 个 PDB")

    ec_map: dict[str, str] = {}
    if args.ec_file.exists():
        ecdf = pd.read_csv(args.ec_file)
        key = "pdb_id" if "pdb_id" in ecdf.columns else ecdf.columns[0]
        val = "ec_number" if "ec_number" in ecdf.columns else ecdf.columns[-1]
        ec_map = dict(zip(ecdf[key].astype(str).str.upper(), ecdf[val].astype(str)))
        emit(f"- EC 注释: {len(ec_map)} 条 (来自 {args.ec_file.name})")
    cath_map: dict[str, str] = {}
    if args.cath_file.exists():
        cdf = pd.read_csv(args.cath_file)
        key = "pdb_id" if "pdb_id" in cdf.columns else cdf.columns[0]
        val = "cath_id" if "cath_id" in cdf.columns else cdf.columns[-1]
        cath_map = dict(zip(cdf[key].astype(str).str.upper(), cdf[val].astype(str)))
        emit(f"- CATH 注释: {len(cath_map)} 条 (来自 {args.cath_file.name})")

    # 每个唯一 PDB 只读一次 HEADER
    uniq_pdbs = sorted(df["pdb_id"].astype(str).str.upper().unique())
    hdr: dict[str, str] = {}
    t_hdr = time.time()
    for i, pid in enumerate(uniq_pdbs):
        p = PDB_DIR / f"{pid}.pdb"
        if not p.exists():
            p = PDB_DIR / f"{pid}.ent"
        hdr[pid] = header_class(p) if p.exists() else ""
        if (i + 1) % 500 == 0:
            emit(f"  ... HEADER 扫描 {i + 1}/{len(uniq_pdbs)} ({time.time() - t_hdr:.0f}s)")
    hdr_counts = Counter(v for v in hdr.values() if v)
    emit(f"- HEADER 扫描完成: {len(uniq_pdbs)} 结构, {sum(hdr_counts.values())} 个带酶类注释")
    emit("")
    emit("| HEADER 类别 | 结构数 | 占比 | 对应 EC 一级类 |")
    emit("|---|---:|---:|---|")
    for k, v in hdr_counts.most_common():
        emit(f"| {k} | {v} | {v / len(uniq_pdbs) * 100:.1f}% | {EC_CLASS_FROM_HEADER[k]} |")
    emit("")
    emit("> **声明**: HEADER 类别是 PDB 沉积时的大分子注释, 用来筛"
         "「这个结构是不是酶」**只是弱证据**。")
    emit("> 它**不是 EC 号**, 一个 HEADER=HYDROLASE 的结构可能是 EC 3.1.1.1 也可能是")
    emit("> EC 3.4.21.4, 精度不足以支撑 EC 层隔离。**因此 HEADER 只用于统计画像,")
    emit("> 不参与 isolation_key 构建。**")
    audit["header_class_counts"] = dict(hdr_counts)

    if args.dry_run:
        emit("")
        emit("--dry-run: 不构建数据集")
        return 0

    # ---- 2. 构造目标 schema -------------------------------------------
    emit("")
    emit("## 2. 目标 schema 构建 (三层隔离)")
    emit("")
    df["_pdb"] = df["pdb_id"].astype(str).str.upper()
    if ASSAY_META_CSV.exists():
        meta = pd.read_csv(ASSAY_META_CSV)
        df = df.merge(meta, on="interface_key", how="left", validate="one_to_one")
        emit(f"- 已并入测量条件元数据 ({args.assay_meta.name})")
    else:
        for c in ("method_majority", "temperature_k_used", "temperature_median_k",
                  "temp_imputed"):
            df[c] = np.nan
        emit("- !! 测量条件元数据缺失, 相关列填 NaN")

    # 序列 (逐结构缓存, 避免重复读文件)
    seq_cache_local: dict[str, str] = {}
    n_seq_ok = 0
    for i, r in enumerate(df.to_dict("records")):
        pid = r["_pdb"]
        if pid in seq_cache_local:
            continue
        lig = [c.strip() for c in str(r["ligand_chains"]).split(",") if c.strip()]
        rec = [c.strip() for c in str(r["receptor_chains"]).split(",") if c.strip()]
        s = get_chain_seq(pid, lig, seqs) + get_chain_seq(pid, rec, seqs)
        seq_cache_local[pid] = s
        n_seq_ok += int(bool(s))
    df["_seq"] = df["_pdb"].map(seq_cache_local)
    emit(f"- 序列可得结构: {n_seq_ok}/{df['_pdb'].nunique()}")

    # EC / CATH 层
    df["ec_number"] = df["_pdb"].map(lambda p: ec_map.get(p, ""))
    df["cath_id"] = df["_pdb"].map(lambda p: cath_map.get(p, ""))
    df["_ec_key"] = df["ec_number"].map(ec_layer)
    df["_cath_key"] = df["cath_id"].map(cath_layer)
    ec_cov = float((df["_ec_key"] != "").mean())
    cath_cov = float((df["_cath_key"] != "").mean())
    emit(f"- EC 层覆盖率: {ec_cov * 100:.1f}%  "
         f"({'已启用' if ec_cov > 0 else '**未接入 -> 该层停用**'})")
    emit(f"- CATH 层覆盖率: {cath_cov * 100:.1f}%  "
         f"({'已启用' if cath_cov > 0 else '**未接入 -> 该层停用**'})")

    # 序列层: 40% 同一性贪心聚类 (只对唯一序列算一次)
    uniq_seqs = sorted({s for s in df["_seq"].dropna().unique() if s})
    emit(f"- 唯一序列数: {len(uniq_seqs)}")
    t_seq = time.time()
    reps: list[tuple[set, str]] = []
    seq_to_cluster: dict[str, str] = {}
    n_new = 0
    for i, s in enumerate(uniq_seqs):
        cid, is_new = sequence_cluster_key(s, reps, args.identity, KMER_K)
        seq_to_cluster[s] = cid
        n_new += is_new
        if (i + 1) % 500 == 0:
            emit(f"  ... 序列聚类 {i + 1}/{len(uniq_seqs)} "
                 f"({time.time() - t_seq:.0f}s, 簇 {n_new})")
    emit(f"- 序列聚类完成: {n_new} 个簇 (耗时 {time.time() - t_seq:.0f}s)")

    df["_seq_key"] = df["_seq"].map(lambda s: seq_to_cluster.get(s, ""))

    # 合并隔离键
    triples = [
        build_isolation_key(e, c, s, p)
        for e, c, s, p in zip(df["_ec_key"], df["_cath_key"], df["_seq_key"], df["_pdb"])
    ]
    df["isolation_key"] = [t[0] for t in triples]
    df["isolation_source"] = [t[1] for t in triples]
    emit("")
    emit("隔离键来源分布:")
    for k, v in df["isolation_source"].value_counts().items():
        emit(f"  - {k}: {v} ({v / len(df) * 100:.1f}%)")

    # ---- 3. 划分 (整族分侧) + 硬校验 ------------------------------------
    emit("")
    emit("## 3. 划分与红线硬校验")
    emit("")
    df["_h"] = df["isolation_key"].map(
        lambda k: int(__import__("hashlib").md5(
            f"ai-v2-enzyme::{k}".encode()).hexdigest()[:8], 16) / 0xFFFFFFFF)
    fam_stat = (df.groupby("isolation_key")
                  .agg(n=("real_dG", "size"), h=("_h", "first"))
                  .reset_index().sort_values(["h", "isolation_key"]))
    fam_stat["cum"] = fam_stat["n"].cumsum()
    cut = (1 - args.test_frac) * len(df)
    best = int((fam_stat["cum"] - cut).abs().idxmin())
    best = max(0, min(best, len(fam_stat) - 2)) if len(fam_stat) > 1 else 0
    train_keys = set(fam_stat.loc[:best, "isolation_key"])
    df["split"] = np.where(df["isolation_key"].isin(train_keys), "train", "test")

    checks = hard_check_no_leak(df, ["isolation_key", "pdb_id", "_seq_key", "_ec_key",
                                     "_cath_key"])
    emit("| 隔离键 | 跨侧违规数 | 要求 | 结果 |")
    emit("|---|---:|---|---|")
    for col, n_viol in checks.items():
        if col.startswith("_") and col not in ("_seq_key", "_ec_key", "_cath_key"):
            continue
        emit(f"| `{col}` | {n_viol if n_viol >= 0 else 'n/a'} | 必须为 0 | "
             f"{'PASS' if n_viol == 0 else ('FAIL' if n_viol > 0 else 'n/a')} |")
    fatal = {k: v for k, v in checks.items()
             if k in ("isolation_key", "pdb_id") and v != 0}
    if fatal:
        emit("")
        emit(f"!! **红线失败**: {fatal} —— 拒绝输出数据集")
        return 3
    emit("")
    emit("- **红线硬校验通过**: `isolation_key` 与 `pdb_id` 跨侧数量均为 0。")

    tr, te = df[df["split"] == "train"], df[df["split"] == "test"]
    emit(f"- 划分: train {len(tr)} 行 / {tr['isolation_key'].nunique()} 隔离组; "
         f"test {len(te)} 行 / {te['isolation_key'].nunique()} 隔离组")

    # ---- 4. train/test 可比性 (提示词明确要求) --------------------------
    emit("")
    emit("## 4. train / test 可比性检查")
    emit("")
    emit("### 4.1 ΔG 分布")
    emit("")
    emit("| 子集 | n | mean | std | min | q25 | median | q75 | max |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for tag, s in (("train", tr["real_dG"]), ("test", te["real_dG"])):
        emit(f"| {tag} | {len(s)} | {s.mean():.3f} | {s.std():.3f} | {s.min():.3f} | "
             f"{s.quantile(.25):.3f} | {s.median():.3f} | {s.quantile(.75):.3f} | "
             f"{s.max():.3f} |")
    from scipy.stats import ks_2samp, mannwhitneyu
    ks = ks_2samp(tr["real_dG"], te["real_dG"])
    mw = mannwhitneyu(tr["real_dG"], te["real_dG"])
    emit("")
    emit(f"- KS 检验: D={ks.statistic:.4f}, p={ks.pvalue:.4f} -> "
         f"{'分布可比' if ks.pvalue > 0.05 else '**分布不可比**'}")
    emit(f"- Mann-Whitney U: p={mw.pvalue:.4f} -> "
         f"{'中位数可比' if mw.pvalue > 0.05 else '**中位数不可比**'}")

    emit("")
    emit("### 4.2 隔离组规模与组内变异 (LTR 可用性的核心指标)")
    emit("")
    emit("| 子集 | 隔离组数 | 多条目组 | 单条目组 | 多条目组覆盖行 | 组内 ΔG 有差异的组 |")
    emit("|---|---:|---:|---:|---:|---:|")
    grp_rows = []
    for tag, sub in (("train", tr), ("test", te), ("全体", df)):
        g = sub.groupby("isolation_key")["real_dG"]
        sz, nu = g.size(), g.nunique()
        multi = sz > 1
        grp_rows.append({
            "split": tag, "n_groups": int(len(sz)), "n_multi": int(multi.sum()),
            "n_single": int((~multi).sum()), "rows_multi": int(sz[multi].sum()),
            "n_varying": int((multi & (nu > 1)).sum()),
        })
        emit(f"| {tag} | {len(sz)} | {int(multi.sum())} | {int((~multi).sum())} | "
             f"{int(sz[multi].sum())} | {int((multi & (nu > 1)).sum())} |")
    g_tr = grp_rows[0]
    eff = float(g_tr["rows_multi"] / max(len(tr), 1))
    emit("")
    emit(f"- **训练侧可用于 group 排序的行占比 = {eff * 100:.1f}%** (多条目组覆盖行)")
    emit(f"- 训练侧组内 ΔG 真有差异的组: {g_tr['n_varying']} / {g_tr['n_multi']} 多条目组")

    emit("")
    emit("### 4.3 隔离组 / 家族数可比性")
    emit("")
    emit("| 子集 | isolation_key 数 | 序列簇数 | EC 一级类数 | 结构数 | 来源数据集数 |")
    emit("|---|---:|---:|---:|---:|---:|")
    for tag, sub in (("train", tr), ("test", te)):
        emit(f"| {tag} | {sub['isolation_key'].nunique()} | "
             f"{sub['_seq_key'].replace('', np.nan).nunique()} | "
             f"{sub['_ec_key'].replace('', np.nan).nunique()} | "
             f"{sub['_pdb'].nunique()} | {sub['source_dataset'].nunique()} |")
    emit("")
    emit("> **外推声明**: 若是家族/隔离组级划分, 测试集里的隔离组在训练集中**完全不可见**。")
    emit("> 因此本表上的任何指标都是**对未见酶/未见家族的外推能力**, 不能当作普遍性能。")
    emit(f"> 本数据 test 侧共 {te['isolation_key'].nunique()} 个隔离组, 全部为训练集未见。")

    # ---- 5. 落盘 -------------------------------------------------------
    schema_cols = [
        "source", "pdb_id", "chain", "sequence", "ec_number", "cath_id",
        "ligand_id", "label_type", "value", "unit", "temperature_K", "ph",
        "assay_method", "real_dG", "group_id", "family_id",
    ]
    out = pd.DataFrame({
        "source": df.get("source_dataset", ""),
        "pdb_id": df["_pdb"],
        "chain": df["ligand_chains"].astype(str) + "|" + df["receptor_chains"].astype(str),
        "sequence": df["_seq"],
        "ec_number": df["ec_number"],
        "cath_id": df["cath_id"],
        "ligand_id": df.get("ligand_name", ""),
        "label_type": "Kd",                       # 源表只有 KD(M), 见 AUTO_ASSUMPTIONS A3
        "value": df.get("kd_m", np.nan),
        "unit": "M",
        "temperature_K": df.get("temperature_k_used", np.nan),
        "ph": np.nan,                             # 源表无 pH 列
        "assay_method": df.get("method_majority", ""),
        "real_dG": df["real_dG"],
        "group_id": df["_seq_key"],
        "family_id": df.get("family_id", ""),
        "isolation_key": df["isolation_key"],
        "isolation_source": df["isolation_source"],
        "split": df["split"],
        "interface_key": df["interface_key"],
        "header_enzyme_class": df["_pdb"].map(hdr),
        "temperature_source": np.where(
            pd.to_numeric(df.get("temp_imputed", np.nan), errors="coerce") == 1,
            "imputed_298.15K", "experimental"),
    })
    args.out_dataset.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out_dataset.with_suffix(args.out_dataset.suffix + ".tmp")
    out.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out_dataset)

    split_df = out[["pdb_id", "interface_key", "isolation_key", "isolation_source",
                    "group_id", "family_id", "real_dG", "split"]].copy()
    tmp2 = args.out_split.with_suffix(args.out_split.suffix + ".tmp")
    split_df.to_csv(tmp2, index=False, encoding="utf-8")
    tmp2.replace(args.out_split)

    emit("")
    emit("## 5. 输出")
    emit("")
    emit(f"- 数据集: `{args.out_dataset.name}` — {len(out)} 行 x {out.shape[1]} 列")
    emit(f"- 划分表: `{args.out_split.name}` — {len(split_df)} 行")
    emit(f"- 目标 schema 字段: {', '.join(schema_cols)}")
    emit("")
    emit("### 未能填充的 schema 字段 (诚实披露)")
    emit("")
    emit("| 字段 | 状态 | 原因 |")
    emit("|---|---|---|")
    emit("| `ph` | **全空** | PPB-Affinity 源表无 pH 列 (19 列已全量枚举) |")
    emit("| `ec_number` | "
         f"{'已填充' if ec_cov > 0 else '**全空**'} | "
         f"{'来自本地 EC 注释文件' if ec_cov > 0 else '需联网获取 ExPASy/SIFTS 映射'} |")
    emit("| `cath_id` | "
         f"{'已填充' if cath_cov > 0 else '**全空**'} | "
         f"{'来自本地 CATH 注释文件' if cath_cov > 0 else '需联网获取 CATH Gene3D 映射'} |")
    emit("| `label_type` | 恒为 `Kd` | 源表只有 KD(M) 一列, 无 Ki/IC50/Km |")
    emit("| `ligand_id` | 部分 | 来自 ligand_name, 非标准化 CHEBI/CCD 编码 |")

    # ---- 6. 下游接入指引 -----------------------------------------------
    emit("")
    emit("## 6. 下游接入指引 (把 EC / CATH / 酶动力学数据接进来)")
    emit("")
    emit("本脚本已实现完整的三层隔离逻辑, **只要提供下列文件即自动启用对应层** ——")
    emit("不需要改代码:")
    emit("")
    emit("1. **EC 层**: 生成 `data/ec_annotations.csv`, 至少含列")
    emit("   `pdb_id,ec_number` (如 `1A22,3.4.21.4`)。")
    emit("   来源: PDB `DBREF` + UniProt 的 EC 号, 或 SABIO-RK 导出, 或 ExPASy ENZYME。")
    emit("2. **CATH 层**: 生成 `data/cath_annotations.csv`, 至少含列")
    emit("   `pdb_id,cath_id` (如 `1A22,3.40.50.300`)。")
    emit("   来源: CATH 官方 `cath-b-s35` / Gene3D 的 PDB 映射文件。")
    emit("3. **酶动力学标签**: SABIO-RK 可导出 Km/kcat/Ki/Kd + 实验条件;")
    emit("   接入后 `label_type` 不再恒为 `Kd`, 阶段 1 的标签净化变体才真正可分")
    emit("   (当前因源表只有 Kd, `Kd_Ki_only` / `Km_aux` 退化为 no-op)。")
    emit("")
    emit("重跑本脚本即自动启用新层:")
    emit("")
    emit("```bash")
    emit("python scripts/13_enzyme_split.py")
    emit("```")
    emit("")
    emit("### 技术理由回顾 (为什么值得接入酶数据)")
    emit("")
    emit(f"- 当前 PPI 数据: 训练侧多条目隔离组 {g_tr['n_multi']} 个, 其中组内 ΔG 真有差异的"
         f"只有 {g_tr['n_varying']} 个")
    emit("  (阶段 3 已独立验证: 61.7% 的'多条目组'是同一测量的链组合冗余, 排序梯度恒为 0)。")
    emit("- 酶数据的天然结构是「同一酶 × 多底物」或「同一酶 × 多突变」,")
    emit("  每个 group 内 ΔG **天然各不相同**, 正是 Learning-to-Rank 需要的形态。")
    emit("")
    emit("## 7. 残留局限")
    emit("")
    emit("1. **序列同一性是 k-mer Jaccard 代理**, 不是真正的比对同一性 (无 CD-HIT/MMseqs2,")
    emit("   以保持纯 CPU 可复现)。该代理**更保守** (簇更细), 见 v1.1 报告残留局限第 3 条。")
    emit("2. **EC / CATH 层未启用**, 当前 isolation_key 实际只由序列层构成;")
    emit("   报告中的 `isolation_source` 列如实标注了每一行的真实来源。")
    emit("3. HEADER 酶类注释只用于画像, **不参与隔离** (精度不足, 一个类别对多个 EC 号)。")
    emit("4. `ph` 字段无法填充; 温度缺失值按 298.15 K 兜底并在")
    emit("   `temperature_source` 列标记为 `imputed_298.15K`。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "identity_threshold": args.identity, "kmer_k": KMER_K,
        "n_rows": int(len(df)), "n_structures": int(df["_pdb"].nunique()),
        "ec_coverage": ec_cov, "cath_coverage": cath_cov,
        "isolation_source_counts": df["isolation_source"].value_counts().to_dict(),
        "leak_checks": checks, "split": {"train": int(len(tr)), "test": int(len(te))},
        "group_stats": grp_rows,
        "header_class_counts": dict(hdr_counts),
        "sources": audit["sources"],
        "unfilled_fields": {"ph": "源表无 pH 列",
                            "ec_number": "未接入" if ec_cov == 0 else "已填充",
                            "cath_id": "未接入" if cath_cov == 0 else "已填充"},
    }
    JSON_OUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str),
                        encoding="utf-8")
    args.report.write_text("\n".join(run_log) + "\n", encoding="utf-8")
    emit("")
    emit(f"报告: {args.report}")
    emit(f"摘要 JSON: {JSON_OUT}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[13_enzyme_split] {traceback.format_exc()}\n")
        raise
