# -*- coding: utf-8 -*-
"""
01_clean_data.py -- Affinity Intelligence / 阶段一：数据基建与清洗

输入 : data/PPB-Affinity.xlsx     (Zenodo record 13054646, 公开)
输出 : data/benchmark.csv         (清洗后主表, 含 train/test 划分与 MVP 标记)
       data/mvp_structures.txt    (MVP 子集需要抓取的唯一 PDB ID 清单)
       data/clean_log.txt         (清洗日志)

关键纪律:
  * 纯 CPU 低内存: 不下载 3.1GB 全量结构包 (实测 350KB/s), 结构按需从 RCSB 取。
  * 标签换算: dG = 0.592 * ln(Kd[M]), 单位 kcal/mol (用户指定公式)。
  * 划分红线: 按靶点家族整体分侧。家族 = entity 序列单链接聚类 (有 RCSB 序列缓存时),
    否则退回「结构身份 (PDB 级)」作为最保守代理, 保证同一结构绝不跨界。
  * 类型纪律: xlsx 空单元格读作 float('nan'), 所有文本字段统一走 sstr() 兜底,
    单行异常只跳过该行并记日志, 不重试不修复; 循环内周期性 gc.collect()。
"""

import gc
import hashlib
import logging
import os
import random
import sys
import time

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- 常量 / 路径
ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
DATA_DIR = os.path.join(ROOT, "data")
XLSX = os.path.join(DATA_DIR, "PPB-Affinity.xlsx")
SEQ_CACHE = os.path.join(DATA_DIR, "pdb_seqs.json")
OUT_CSV = os.path.join(DATA_DIR, "benchmark.csv")
OUT_MVP = os.path.join(DATA_DIR, "mvp_structures.txt")
LOG_PATH = os.path.join(DATA_DIR, "clean_log.txt")

RT_LITERAL = 0.592            # kcal/mol, 用户指定
FAMILY_IDENTITY = 0.90        # 3-mer Jaccard 阈值, >= 判为同一靶点家族
TRAIN_FRACTION = 0.80
HASH_SALT = "affinity-intelligence-v0.1"
SEED = 20240101
MVP_ROWS = 400                # 降级快速验证集规模 (用户要求 200-500)
MVP_MAX_PER_STRUCT = 1        # 同一结构(界面)在 MVP 中最多占几行

# ---------------------------------------------------------------- 日志
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

logger = logging.getLogger("clean")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)


# ---------------------------------------------------------------- 工具函数
def sstr(v):
    """类型安全的字符串转换: NaN/None/float 一律安全。"""
    if v is None:
        return ""
    if isinstance(v, float):
        if np.isnan(v):
            return ""
        return repr(v)
    return str(v)


def norm_chain(v):
    """'A, B' / 'a;b' -> 'A,B' (大写、去重、排序)。"""
    txt = sstr(v).replace(";", ",").replace("|", ",")
    parts = [p.strip().upper() for p in txt.split(",") if p.strip()]
    parts = [p.split(":")[0].strip() for p in parts]
    parts = [p for p in parts if p]
    return ",".join(sorted(dict.fromkeys(parts)))


def dg_from_kd(kd_m):
    """dG = 0.592 * ln(Kd[M]) -> kcal/mol (用户指定公式)。"""
    return RT_LITERAL * float(np.log(kd_m))


def hash_unit(key):
    """确定性映射到 [0,1)。"""
    h = hashlib.md5((HASH_SALT + "::" + str(key)).encode("utf-8")).hexdigest()
    return int(h[:8], 16) / 0xFFFFFFFF


def cluster_families(sequences, threshold=FAMILY_IDENTITY):
    """单链接贪心聚类: 3-mer 集合 Jaccard >= threshold 视为同家族。

    纯 CPU / 低内存代理, 避免 3000+ 序列的 O(n^2) 全比对。
    """
    def kmers(seq, k=3):
        if len(seq) < k:
            return {seq} if seq else set()
        return {seq[i:i + k] for i in range(len(seq) - k + 1)}

    reps, assign = [], {}
    for seq in sequences:
        ks = kmers(seq)
        hit = None
        for rep_ks, uid in reps:
            union = len(ks | rep_ks)
            if union and len(ks & rep_ks) / union >= threshold:
                hit = uid
                break
        if hit is None:
            hit = "fam%04d" % len(reps)
            reps.append((ks, hit))
        assign[seq] = hit
    return assign


# ---------------------------------------------------------------- 主流程
def main():
    t0 = time.time()
    summary = {}

    # ---- 1. 读取原始标签表 ------------------------------------------------
    logger.info("读取标签表: %s", XLSX)
    raw = pd.read_excel(XLSX, sheet_name="Sheet1", dtype=str)
    logger.info("原始行数: %d, 列: %d", len(raw), raw.shape[1])
    summary["raw_rows"] = len(raw)

    colmap = {
        "Source Data Set": "source_dataset",
        "Complex ID": "complex_id",
        "PDB": "pdb_id",
        "Mutations": "mutations",
        "Ligand Chains": "ligand_chains",
        "Receptor Chains": "receptor_chains",
        "Ligand Name": "ligand_name",
        "Receptor Name": "receptor_name",
        "KD(M)": "kd_m",
        "Affinity Method": "affinity_method",
        "Temperature(K)": "temperature_k",
        "Resolution(Å)": "resolution_a",
        "Subgroup": "subgroup",
    }
    missing = [c for c in colmap if c not in raw.columns]
    if missing:
        raise RuntimeError("标签表缺少必需列: %s" % missing)
    df = raw[list(colmap)].rename(columns=colmap).copy()
    del raw
    gc.collect()

    # ---- 2. 逐行清洗 (类型安全; 坏行跳过并记日志) ------------------------
    rows, bad = [], {"no_pdb": 0, "no_chains": 0, "bad_kd": 0, "err": 0}
    for i, rec in enumerate(df.to_dict("records")):
        try:
            pdb_id = sstr(rec["pdb_id"]).strip().upper()
            if not pdb_id or pdb_id in ("NAN", "NONE"):
                bad["no_pdb"] += 1
                continue

            lig = norm_chain(rec["ligand_chains"])
            rec_c = norm_chain(rec["receptor_chains"])
            if not lig or not rec_c:
                bad["no_chains"] += 1
                continue

            try:
                kd = float(sstr(rec["kd_m"]).strip())
            except (TypeError, ValueError):
                bad["bad_kd"] += 1
                continue
            if not np.isfinite(kd) or kd <= 0:
                bad["bad_kd"] += 1
                continue

            rows.append({
                "pdb_id": pdb_id,
                "receptor_chains": rec_c,
                "ligand_chains": lig,
                "real_dG": dg_from_kd(kd),
                "kd_m": kd,
                "source_dataset": sstr(rec["source_dataset"]).strip(),
                "complex_id": sstr(rec["complex_id"]).strip(),
                "mutations": sstr(rec["mutations"]).strip(),
                "ligand_name": sstr(rec["ligand_name"]).strip(),
                "receptor_name": sstr(rec["receptor_name"]).strip(),
                "affinity_method": sstr(rec["affinity_method"]).strip(),
                "temperature_k": sstr(rec["temperature_k"]).strip(),
                "resolution_a": sstr(rec["resolution_a"]).strip(),
                "subgroup": sstr(rec["subgroup"]).strip(),
            })
        except Exception as exc:                       # noqa: BLE001 - 坏行绝不阻断
            bad["err"] += 1
            if bad["err"] <= 3:
                import traceback
                logger.error("坏行 #%d (row=%d) 真实异常:\n%s",
                             bad["err"], i, traceback.format_exc())
        if (i + 1) % 2500 == 0:
            gc.collect()

    df = pd.DataFrame(rows)
    del rows
    gc.collect()
    logger.info("清洗后行数: %d (丢弃: 无PDB=%d, 缺链=%d, 坏Kd=%d, 其他异常=%d)",
                len(df), bad["no_pdb"], bad["no_chains"], bad["bad_kd"], bad["err"])
    summary.update(bad)
    summary["clean_rows"] = len(df)
    if df.empty:
        raise RuntimeError("清洗后无数据, 终止")

    # ---- 3. 派生键 --------------------------------------------------------
    # interface_key = 具体界面的测量条目 (同 PDB 同链组合); 突变体会共享它
    df["interface_key"] = (df["pdb_id"] + "|" + df["ligand_chains"] + "|"
                           + df["receptor_chains"])
    # structure_key = 结构身份, 用于 MVP 结构级封顶
    df["structure_key"] = df["pdb_id"]

    df["pdb_path"] = df["pdb_id"].map(
        lambda p: os.path.join("data", "pdb", "%s.pdb" % p))

    # ---- 4. 合并同一测量条目的重复读数 (取均值, 降噪) --------------------
    before = len(df)
    agg = (df.sort_values("complex_id")
             .groupby("interface_key", as_index=False)
             .agg(pdb_id=("pdb_id", "first"),
                  structure_key=("structure_key", "first"),
                  receptor_chains=("receptor_chains", "first"),
                  ligand_chains=("ligand_chains", "first"),
                  real_dG=("real_dG", "mean"),
                  kd_m=("kd_m", "mean"),
                  n_measure=("real_dG", "size"),
                  source_dataset=("source_dataset", "first"),
                  complex_id=("complex_id", "first"),
                  mutations=("mutations", "first"),
                  has_mutation=("mutations", lambda s: int((s != "").any())),
                  ligand_name=("ligand_name", "first"),
                  receptor_name=("receptor_name", "first"),
                  affinity_method=("affinity_method", "first"),
                  temperature_k=("temperature_k", "first"),
                  resolution_a=("resolution_a", "first"),
                  subgroup=("subgroup", "first")))
    agg["pdb_path"] = agg["pdb_id"].map(
        lambda p: os.path.join("data", "pdb", "%s.pdb" % p))
    logger.info("按界面条目合并重复读数: %d -> %d 条 (唯一界面 %d, 唯一结构 %d)",
                before, len(agg), agg["interface_key"].nunique(),
                agg["pdb_id"].nunique())
    summary["entries"] = len(agg)
    summary["unique_structures"] = int(agg["pdb_id"].nunique())
    del df
    gc.collect()

    # ---- 5. 靶点家族 (序列聚类优先, 结构身份兜底) -----------------------
    seqs = {}
    if os.path.exists(SEQ_CACHE):
        try:
            import json
            with open(SEQ_CACHE, "r", encoding="utf-8") as fh:
                seqs = json.load(fh)
            logger.info("载入 RCSB 序列缓存: %d 个 PDB", len(seqs))
        except Exception as exc:                        # noqa: BLE001
            logger.warning("序列缓存读取失败, 退回结构身份分组: %s", exc)
            seqs = {}

    if seqs:
        uniq_receptor_seqs = sorted({s for s in seqs.values() if s})
        fam_map = cluster_families(uniq_receptor_seqs, FAMILY_IDENTITY)
        agg["family_id"] = agg["pdb_id"].map(
            lambda p: fam_map.get(seqs.get(p, "")) or ("fam_struct_" + p))
        fam_src = "RCSB entity 序列单链接聚类 (3-mer Jaccard >= %.2f)" % FAMILY_IDENTITY
    else:
        agg["family_id"] = "fam_struct_" + agg["pdb_id"]
        fam_src = "结构身份 (每个 PDB 一个家族) + 家族级哈希兜底"
    logger.info("划分依据: %s; 家族数: %d", fam_src, agg["family_id"].nunique())
    summary["families"] = int(agg["family_id"].nunique())
    summary["family_source"] = fam_src

    # ---- 6. 按家族整体分侧 (家族级排序取前 80% 行 -> 保证不跨界) ---------
    fam_stat = (agg.groupby("family_id")
                   .agg(n=("real_dG", "size"))
                   .reset_index())
    fam_stat["h"] = fam_stat["family_id"].map(hash_unit)
    fam_stat = fam_stat.sort_values(["h", "family_id"]).reset_index(drop=True)
    fam_stat["cum"] = fam_stat["n"].cumsum()
    cut = TRAIN_FRACTION * len(agg)
    # 找到累积行数最接近 80% 的切点, 且保证两侧都非空
    best_idx = int((fam_stat["cum"] - cut).abs().idxmin())
    n_fam = len(fam_stat)
    best_idx = max(0, min(best_idx, n_fam - 2)) if n_fam > 1 else 0
    train_fams = set(fam_stat.loc[:best_idx, "family_id"])
    agg["split"] = np.where(agg["family_id"].isin(train_fams), "train", "test")

    # ---- 7. 红线校验: 同一家族 / 同一 PDB 绝不跨界 ----------------------
    cross_fam = int((agg.groupby("family_id")["split"].nunique() > 1).sum())
    cross_pdb = int((agg.groupby("pdb_id")["split"].nunique() > 1).sum())
    if cross_fam or cross_pdb:
        raise RuntimeError("数据泄漏! 跨界家族=%d, 跨界PDB=%d" % (cross_fam, cross_pdb))
    logger.info("红线校验通过: 0 个家族 / 0 个 PDB 同时出现在 train 与 test")

    # ---- 8. MVP 快速验证集 (降级; 结构级封顶保证多样性) -----------------
    rng = random.Random(SEED)
    agg["is_mvp"] = False
    work = agg.copy()
    # 同一界面优先保留"无突变(野生型)"读数, 再按 complex_id 稳定排序
    work["_wt"] = (work["has_mutation"] == 0).astype(int)
    work = work.sort_values(["structure_key", "_wt", "complex_id"],
                            ascending=[True, False, True])
    picked, per_struct = [], {}
    for rec in work.to_dict("records"):
        sk = rec["structure_key"]
        if per_struct.get(sk, 0) >= MVP_MAX_PER_STRUCT:
            continue
        # 第一遍: 每个结构取满 1 条; 第二遍再补 (本处 1 条/结构)
        picked.append(rec["interface_key"])
        per_struct[sk] = per_struct.get(sk, 0) + 1
        if len(picked) >= MVP_ROWS:
            break
    if len(picked) < MVP_ROWS:
        # 额度没用完 (结构数不足), 放开封顶, 按突变体补齐
        rest = [r for r in work.to_dict("records")
                if r["interface_key"] not in set(picked)]
        rng.shuffle(rest)
        for rec in rest:
            if len(picked) >= MVP_ROWS:
                break
            picked.append(rec["interface_key"])
    mvp_keys = set(picked)
    agg.loc[agg["interface_key"].isin(mvp_keys), "is_mvp"] = True
    mvp = agg[agg["is_mvp"]]

    # ---- 9. 汇总 ----------------------------------------------------------
    tr, te = agg[agg["split"] == "train"], agg[agg["split"] == "test"]
    summary.update({
        "train_rows": int(len(tr)), "test_rows": int(len(te)),
        "train_ratio": round(len(tr) / max(len(agg), 1), 4),
        "train_structures": int(tr["pdb_id"].nunique()),
        "test_structures": int(te["pdb_id"].nunique()),
        "mvp_rows": int(len(mvp)),
        "mvp_structures": int(mvp["pdb_id"].nunique()),
        "dg_min": round(float(agg["real_dG"].min()), 3),
        "dg_max": round(float(agg["real_dG"].max()), 3),
        "dg_mean": round(float(agg["real_dG"].mean()), 3),
        "dg_std": round(float(agg["real_dG"].std()), 3),
    })

    # ---- 10. 落盘 ---------------------------------------------------------
    out_cols = ["pdb_id", "receptor_chains", "ligand_chains", "real_dG", "pdb_path",
                "kd_m", "interface_key", "family_id", "split", "is_mvp",
                "source_dataset", "complex_id", "mutations", "has_mutation",
                "subgroup", "affinity_method", "temperature_k", "resolution_a",
                "n_measure", "ligand_name", "receptor_name"]
    required = ["pdb_id", "receptor_chains", "ligand_chains", "real_dG", "pdb_path"]
    assert all(c in out_cols for c in required)
    agg[out_cols].to_csv(OUT_CSV, index=False, encoding="utf-8")
    logger.info("写出 %s (%d 行)", OUT_CSV, len(agg))

    mvp_pdbs = sorted(mvp["pdb_id"].unique())
    with open(OUT_MVP, "w", encoding="utf-8") as fh:
        fh.write("\n".join(mvp_pdbs) + "\n")
    logger.info("写出 %s (%d 个唯一 PDB)", OUT_MVP, len(mvp_pdbs))

    # ---- 11. 报告 ---------------------------------------------------------
    print("\n" + "=" * 80)
    print("阶段一完成")
    print("=" * 80)
    print("benchmark.csv 行数 : %d    (原始 %d 行 -> 清洗/合并重复读数后)"
          % (len(agg), summary["raw_rows"]))
    print("唯一界面 / 唯一结构 : %d / %d"
          % (agg["interface_key"].nunique(), summary["unique_structures"]))
    print("划分比例           : train %d 行 (%.1f%%, %d 结构) / test %d 行 (%.1f%%, %d 结构)"
          % (len(tr), 100 * summary["train_ratio"], summary["train_structures"],
             len(te), 100 * (1 - summary["train_ratio"]), summary["test_structures"]))
    print("划分依据           : %s" % fam_src)
    print("家族数             : %d  (整族分侧, 0 跨界)" % summary["families"])
    print("红线校验           : PASS — 0 个 PDB、0 个家族同时出现在训练/测试集")
    print("MVP 快速验证集      : %d 行 / %d 个唯一结构 (降级集)"
          % (summary["mvp_rows"], summary["mvp_structures"]))
    print("real_dG 分布       : %.3f ~ %.3f kcal/mol (mean %.3f, std %.3f)"
          % (summary["dg_min"], summary["dg_max"], summary["dg_mean"], summary["dg_std"]))
    print("必需列             : %s" % ", ".join(required))
    print("\n前 5 行 (必需列):")
    print(agg[required].head(5).to_string(index=False))
    print("\n前 5 行 (划分与来源):")
    print(agg[["pdb_id", "ligand_chains", "receptor_chains", "real_dG",
               "family_id", "split", "is_mvp", "subgroup", "source_dataset"]]
          .head(5).to_string(index=False))
    print("\n标签换算: real_dG = %.3f * ln(Kd[M])   [kcal/mol]" % RT_LITERAL)
    print("日志: %s" % LOG_PATH)
    print("耗时: %.1fs" % (time.time() - t0))
    print("=" * 80)
    print("阶段一完成，等待下一步指令")


if __name__ == "__main__":
    main()
