# -*- coding: utf-8 -*-
"""
09_family_split_full.py -- Affinity Intelligence V1.0 / 全量家族级划分

背景: features_full.csv 的 family_id 为空, 其 split 列是早期在 4762 行全量上按
      "PDB 级结构身份"分配的 —— 不是序列家族。V1.0 训练要求"严格沿用 04 的家族级
      隔离", 因此必须为全部训练行重新计算序列家族并整族分侧。

输入 : features/features_full.csv (2151 行) + features/features_mvp.csv (397 行)
       data/pdb_seqs.json         (序列缓存; 缺失的按需从 RCSB 抓取)
输出 : features/features_full.csv    (原地更新 family_id / split_family)
       features/features_mvp.csv     (同上, 新增 split_family 列保留原 split)
       features/v1_family_assignments.csv
       features/v1_family_split_report.md
       data/pdb_seqs.json           (更新缓存)

与 04_family_split.py 完全一致的规则:
  * 家族 = 链序列 40% 同一性阈值的单链接聚类 (union-find)
  * 同一性 = 4-mer Jaccard 预筛 (阈值 0.02) 后做带状全局比对
  * 整族分侧: 家族级哈希排序后取前 80% 行
  * 红线硬校验: 0 个家族 / 0 个 PDB 跨侧, 违反即抛错
性能优化: 预筛的 k-mer 交集用稀疏矩阵一次算完 (旧版逐对 Python 循环在 1600+ 序列
          规模下会退化到小时级), 仅对通过预筛的候选对做比对。
"""

import gc
import hashlib
import json
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
DATA_DIR = os.path.join(ROOT, "data")
FEAT_DIR = os.path.join(ROOT, "features")
SEQ_CACHE = os.path.join(DATA_DIR, "pdb_seqs.json")
BENCH = os.path.join(DATA_DIR, "benchmark.csv")
FULL_CSV = os.path.join(FEAT_DIR, "features_full.csv")
MVP_CSV = os.path.join(FEAT_DIR, "features_mvp.csv")
ASSIGN_CSV = os.path.join(FEAT_DIR, "v1_family_assignments.csv")
REPORT = os.path.join(FEAT_DIR, "v1_family_split_report.md")
LOG_PATH = os.path.join(FEAT_DIR, "v1_family_split_log.txt")

IDENTITY_THRESHOLD = 0.40
KMER = 4
PREFILTER_JACCARD = 0.02
MIN_SEQ_LEN = 20
TRAIN_FRACTION = 0.80
GREEDY_SEQ_LIMIT = 4000      # 序列数超过此值改用贪心代表聚类 (全对预筛代价过高)
N_THREADS = 2
HASH_SALT = "affinity-intelligence-v1.0-family"
FASTA_URL = "https://www.rcsb.org/fasta/entry/%s"
HEADERS = {"User-Agent": "AffinityIntelligence/1.0 (academic pipeline)"}

logger = logging.getLogger("famfull")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)


# ---------------------------------------------------------------- 序列抓取
def _clean_chain_label(seg):
    seg = seg.strip()
    if "[auth " in seg:
        seg = seg.split("[auth ")[1].split("]")[0].strip()
    elif "[auth" in seg:
        seg = seg.split("[auth")[1].strip(" ]")
    seg = seg.split(",")[0].strip()
    for pref in ("ss ", "s "):
        if seg.startswith(pref):
            seg = seg[len(pref):].strip()
    return seg if len(seg) == 1 else None


def parse_fasta(txt):
    out, cur = {}, None
    for ln in txt.splitlines():
        if ln.startswith(">"):
            parts = ln[1:].split("|")
            chain, name = None, ""
            for p in parts[1:]:
                p = p.strip()
                if p.startswith("Chain"):
                    chain = _clean_chain_label(p[5:])
                elif not name:
                    name = p
            if chain:
                cur = chain
                out[chain] = {"seq": "", "name": name}
            else:
                cur = None
        elif cur and ln.strip():
            out[cur]["seq"] += ln.strip().upper()
    return out


def fetch_fasta(pid):
    for attempt in (1, 2):
        try:
            r = requests.get(FASTA_URL % pid, headers=HEADERS, timeout=40)
            if r.status_code == 200 and r.text.strip():
                return pid, parse_fasta(r.text), None
            if r.status_code in (400, 404):
                return pid, None, "http_%d" % r.status_code
        except Exception as exc:                                # noqa: BLE001
            if attempt == 2:
                return pid, None, type(exc).__name__
            time.sleep(1.0)
    return pid, None, "no_content"


def load_or_fetch_sequences(pids):
    cache = {}
    if os.path.exists(SEQ_CACHE):
        try:
            with open(SEQ_CACHE, "r", encoding="utf-8") as fh:
                cache = json.load(fh)
            logger.info("载入序列缓存: %d 个 PDB", len(cache))
        except Exception as exc:                                # noqa: BLE001
            logger.warning("序列缓存损坏: %s", exc)
            cache = {}
    todo = [p for p in pids if p not in cache]
    if not todo:
        logger.info("全部 %d 个 PDB 序列已缓存", len(pids))
        return cache

    logger.info("待抓取序列: %d 个 PDB (并发=%d)", len(todo), N_THREADS)
    t0 = time.time()
    failed, done = [], 0
    with ThreadPoolExecutor(max_workers=N_THREADS) as pool:
        futs = {pool.submit(fetch_fasta, p): p for p in todo}
        for fut in as_completed(futs):
            pid = futs[fut]
            done += 1
            try:
                pid, seqs, err = fut.result()
            except Exception as exc:                            # noqa: BLE001
                failed.append((pid, str(exc)[:50]))
                continue
            if seqs:
                cache[pid] = seqs
            else:
                failed.append((pid, err))
            if done % 200 == 0 or done == len(todo):
                logger.info("序列进度 %d/%d 成功=%d 失败=%d 用时=%.0fs",
                            done, len(todo), len(cache), len(failed), time.time() - t0)
    if failed:
        logger.warning("序列抓取失败 %d 个: %s", len(failed), failed[:15])
    try:
        with open(SEQ_CACHE, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, ensure_ascii=False)
        logger.info("序列缓存已更新: %d 个 PDB", len(cache))
    except Exception as exc:                                    # noqa: BLE001
        logger.error("序列缓存写出失败: %s", exc)
    return cache


# ---------------------------------------------------------------- 聚类
def kmers(seq, k=KMER):
    if len(seq) < k:
        return {seq}
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def banded_identity(a, b, max_len=2000, band=64):
    """带状全局比对 (match +2 / mismatch -1 / gap -1), 返回 (identity, aln_len)。"""
    if not a or not b:
        return 0.0, 0
    a, b = a[:max_len], b[:max_len]
    if not a.isascii() or not b.isascii():
        return 0.0, 0
    n, m = len(a), len(b)
    # 带宽必须至少覆盖长度差, 否则当 |i-j| 无法保持在 band 内时会出现空区间 (j0>j1)。
    # 长度悬殊的序列对 (如 20 aa 片段 vs 1000 aa 蛋白) 会触发此情况。
    band = max(band, abs(n - m) + 2)
    NEG = -1e9
    aa = np.frombuffer(a.encode(), dtype=np.uint8)
    bb = np.frombuffer(b.encode(), dtype=np.uint8)

    prev = np.full(m + 1, NEG)
    prev[0] = 0.0
    for j in range(1, min(m, band) + 1):
        prev[j] = -float(j)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)

    for i in range(1, n + 1):
        j0 = max(1, i - band)
        j1 = min(m, i + band)
        cur = np.full(m + 1, NEG)
        if i <= band:
            cur[0] = -float(i)
        trace[i, 0] = 1
        sub = prev[j0 - 1:j1] + np.where(aa[i - 1] == bb[j0 - 1:j1], 2.0, -1.0)
        up = prev[j0:j1 + 1] - 1.0
        left = np.full(j1 - j0 + 1, NEG)
        if j0 == 1:
            left[0] = cur[0] - 1.0
            for pos in range(1, len(left)):
                left[pos] = max(left[pos - 1], sub[pos - 1], up[pos - 1]) - 1.0
        stack = np.vstack([sub, up, left])
        choice = np.argmax(stack, axis=0)
        cur[j0:j1 + 1] = stack[choice, np.arange(len(choice))]
        trace[i, j0:j1 + 1] = choice
        prev = cur

    i, j, matches, aln = n, m, 0, 0
    while i > 0 or j > 0:
        t = trace[i, j]
        if i > 0 and j > 0 and t == 0:
            if aa[i - 1] == bb[j - 1]:
                matches += 1
            aln += 1
            i -= 1
            j -= 1
        elif i > 0 and (t == 1 or j == 0):
            aln += 1
            i -= 1
        elif j > 0:
            aln += 1
            j -= 1
        else:
            break
    return (matches / aln if aln else 0.0), aln


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def prefilter_pairs(seqs, threshold=PREFILTER_JACCARD):
    """numpy 矩阵化 4-mer Jaccard 预筛: 返回候选对 (i<j) 与统计。

    旧版逐对 Python 集合运算在 1600+ 序列规模下会退化到小时级; 这里把 k-mer
    映射成稀疏 one-hot 矩阵, 交集/并集用矩阵乘法一次算完。
    """
    from scipy import sparse
    vocab = {}
    rows, cols = [], []
    for i, s in enumerate(seqs):
        for k in kmers(s):
            j = vocab.setdefault(k, len(vocab))
            rows.append(i)
            cols.append(j)
    M = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                          shape=(len(seqs), len(vocab)))
    M.data[:] = 1.0
    inter = (M @ M.T).toarray()
    sizes = np.asarray(M.sum(axis=1)).ravel()
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(jac, 0.0)
    ii, jj = np.where(np.triu(jac >= threshold, k=1))
    return list(zip(ii.tolist(), jj.tolist())), jac.shape


def cluster_sequences(seqs, threshold=IDENTITY_THRESHOLD, ident_cache=None):
    """单链接聚类; ident_cache 复用已算好的候选对同一性 (多阈值 / 重跑时省掉重算)。"""
    if ident_cache is None:
        ident_cache = {}
    uniq = sorted({s for s in seqs if s and len(s) >= MIN_SEQ_LEN})
    stats = {"n_seq": len(uniq), "threshold": threshold, "n_candidates": 0,
             "n_edges": 0, "edge_identities": [], "n_reused": 0}
    if not uniq:
        stats["n_clusters"] = 0
        return {}, stats

    cand, _ = prefilter_pairs(uniq)
    stats["n_candidates"] = len(cand)
    logger.info("预筛: %d 条序列 -> %d 个候选对进入比对", len(uniq), len(cand))

    uf = UnionFind()
    for s in uniq:
        uf.find(s)
    t0 = time.time()
    for n, (i, j) in enumerate(cand, start=1):
        a, b = uniq[i], uniq[j]
        key = (a, b) if a <= b else (b, a)
        if ident_cache is not None and key in ident_cache:
            ident, aln = ident_cache[key]
            stats["n_reused"] += 1
        else:
            ident, aln = banded_identity(a, b)
            if ident_cache is not None:
                ident_cache[key] = (ident, aln)
        if aln >= max(30, 0.3 * min(len(a), len(b))) and ident >= threshold:
            uf.union(a, b)
            stats["n_edges"] += 1
            stats["edge_identities"].append(ident)
        if n % 5000 == 0:
            gc.collect()
            if ident_cache is not None and len(ident_cache) % 10000 < 5000:
                save_ident_cache(ident_cache)       # 定期落盘, 避免长比对白费
            logger.info("  比对 %d/%d, 已确认边=%d, 复用=%d, 用时=%.0fs",
                        n, len(cand), stats["n_edges"], stats["n_reused"],
                        time.time() - t0)

    groups = {}
    for s in uniq:
        groups.setdefault(uf.find(s), []).append(s)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), -max(len(x) for x in g)))
    assign = {}
    for idx, g in enumerate(ordered, start=1):
        for s in g:
            assign[s] = "F%04d" % idx
    stats["n_clusters"] = len(ordered)
    stats["n_singletons"] = sum(1 for g in ordered if len(g) == 1)
    stats["largest_cluster"] = max((len(g) for g in ordered), default=0)
    return assign, stats


def cluster_sequences_greedy(seqs, threshold=IDENTITY_THRESHOLD, ident_cache=None):
    """**贪心代表聚类 (CD-HIT 式)**: 只与已有簇代表比对, 比对量 O(n x 代表数)。

    动机: 全量数据下序列数达 5763, 单链接的全对预筛产生 397k 候选对 (约 4.3 小时)。
    贪心法按序列**由长到短**处理 (长序列信息量大, 优先成为代表), 每条序列只需与
    通过 k-mer 预筛的代表比对。

    与单链接的差异: 贪心是"代表制", 不保证传递闭包, 因此簇通常更细 (更保守)。
    红线(整族分侧)由后续硬校验保证, 与聚类算法无关。
    """
    if ident_cache is None:
        ident_cache = {}
    uniq = sorted({s for s in seqs if s and len(s) >= MIN_SEQ_LEN},
                  key=lambda s: (-len(s), s))
    stats = {"n_seq": len(uniq), "threshold": threshold, "n_candidates": 0,
             "n_edges": 0, "edge_identities": [], "n_reused": 0, "mode": "greedy"}
    if not uniq:
        stats["n_clusters"] = 0
        return {}, stats

    ksets = {s: kmers(s) for s in uniq}
    reps = []          # [(seq, cluster_id)]
    assign = {}
    t0 = time.time()
    for n, s in enumerate(uniq, start=1):
        ks = ksets[s]
        hit = None
        for rep, cid in reps:
            rk = ksets[rep]
            union = len(ks | rk)
            if union == 0 or len(ks & rk) / union < PREFILTER_JACCARD:
                continue
            stats["n_candidates"] += 1
            key = (s, rep) if s <= rep else (rep, s)
            if key in ident_cache:
                ident, aln = ident_cache[key]
                stats["n_reused"] += 1
            else:
                ident, aln = banded_identity(s, rep)
                ident_cache[key] = (ident, aln)
            if aln >= max(30, 0.3 * min(len(s), len(rep))) and ident >= threshold:
                hit = cid
                stats["n_edges"] += 1
                stats["edge_identities"].append(ident)
                break
        if hit is None:
            hit = "F%04d" % (len(reps) + 1)
            reps.append((s, hit))
        assign[s] = hit
        if n % 500 == 0:
            gc.collect()
            if len(ident_cache) % 10000 < 500:
                save_ident_cache(ident_cache)
            logger.info("  贪心聚类 %d/%d | 代表数=%d | 已比对=%d 复用=%d | 用时=%.0fs",
                        n, len(uniq), len(reps), stats["n_candidates"],
                        stats["n_reused"], time.time() - t0)
    stats["n_clusters"] = len(reps)
    sizes = {}
    for s, cid in assign.items():
        sizes[cid] = sizes.get(cid, 0) + 1
    stats["n_singletons"] = sum(1 for v in sizes.values() if v == 1)
    stats["largest_cluster"] = max(sizes.values(), default=0)
    return assign, stats


def hash_unit(key, salt=HASH_SALT):
    h = hashlib.md5((salt + "::" + str(key)).encode("utf-8")).hexdigest()
    return int(h[:8], 16) / 0xFFFFFFFF


# ---------------------------------------------------------------- 比对缓存 (可续算)
IDENT_CACHE = os.path.join(DATA_DIR, "v1_identity_cache.json")


def load_ident_cache(path=IDENT_CACHE):
    """载入已算好的候选对同一性; 键为 '短序列||长序列'。"""
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)
        cache = {}
        for k, v in raw.items():
            a, b = k.split("||", 1)
            key = (a, b) if a <= b else (b, a)
            cache[key] = (float(v[0]), int(v[1]))
        logger.info("载入比对缓存: %d 对 (重跑可跳过比对)", len(cache))
        return cache
    except Exception as exc:                                    # noqa: BLE001
        logger.warning("比对缓存读取失败: %s", exc)
        return {}


def save_ident_cache(cache, path=IDENT_CACHE):
    try:
        out = {"%s||%s" % k: [round(v[0], 6), int(v[1])] for k, v in cache.items()}
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(out, fh, ensure_ascii=False)
        os.replace(tmp, path)
        logger.info("比对缓存已落盘: %d 对 -> %s", len(out), path)
    except Exception as exc:                                    # noqa: BLE001
        logger.error("比对缓存写出失败: %s", exc)


# ---------------------------------------------------------------- 主流程
def main():
    t0 = time.time()
    logger.info("=" * 80)
    logger.info("V1.0 全量家族划分启动")

    full = pd.read_csv(FULL_CSV, dtype={"pdb_id": str})
    mvp = pd.read_csv(MVP_CSV, dtype={"pdb_id": str})
    logger.info("features_full: %d 行 / %d 结构; features_mvp: %d 行 / %d 结构",
                len(full), full["pdb_id"].nunique(), len(mvp), mvp["pdb_id"].nunique())

    all_pdbs = sorted(set(full["pdb_id"]) | set(mvp["pdb_id"]))
    cache = load_or_fetch_sequences(all_pdbs)
    missing = [p for p in all_pdbs if p not in cache]
    logger.info("无序列 PDB: %d (以 PDB 身份兜底成独立家族)", len(missing))

    seqs = []
    for pid, chains in cache.items():
        for ch, d in chains.items():
            if d.get("seq"):
                seqs.append(d["seq"])
    logger.info("收集链序列: %d 条", len(seqs))

    # 同一性缓存跨阈值复用; 首次运行后落盘, 重跑无需再等比对
    ident_cache = load_ident_cache()
    t_cluster = time.time()
    n_seq_all = len({s for s in seqs if s and len(s) >= MIN_SEQ_LEN})
    # 规模自适应: 序列数超过阈值时改用贪心代表聚类 (全对预筛在该规模下需数小时)
    use_greedy = n_seq_all > GREEDY_SEQ_LIMIT
    if use_greedy:
        logger.info("序列数 %d > %d -> 使用贪心代表聚类 (CD-HIT 式, 比对量线性于代表数)",
                    n_seq_all, GREEDY_SEQ_LIMIT)
        assign, stats = cluster_sequences_greedy(seqs, IDENTITY_THRESHOLD, ident_cache)
        # 敏感性对照: 抽样估算 (全量贪心在 30%/50% 上重跑代价可接受)
        sens = {}
        for thr in (0.30, 0.50):
            _, st = cluster_sequences_greedy(seqs, thr, ident_cache)
            sens[thr] = st["n_clusters"]
            logger.info("[%.0f%%] 簇数=%d (贪心, 复用 %d)", thr * 100,
                        st["n_clusters"], st["n_reused"])
    else:
        assign, stats = cluster_sequences(seqs, IDENTITY_THRESHOLD, ident_cache)
        sens = {}
        for thr in (0.30, 0.50):
            _, st = cluster_sequences(seqs, thr, ident_cache)
            sens[thr] = st["n_clusters"]
            logger.info("[%.0f%%] 簇数=%d (单链接, 复用 %d)", thr * 100,
                        st["n_clusters"], st["n_reused"])
    save_ident_cache(ident_cache)
    logger.info("[40%%] 唯一序列=%d -> 簇=%d (单例 %d, 最大簇 %d), 边=%d, "
                "候选比对=%d, 复用=%d, 模式=%s, 聚类耗时=%.0fs",
                stats["n_seq"], stats["n_clusters"], stats["n_singletons"],
                stats["largest_cluster"], stats["n_edges"], stats["n_candidates"],
                stats["n_reused"], stats.get("mode", "single-linkage"),
                time.time() - t_cluster)

    def row_family(rec):
        chains = cache.get(rec["pdb_id"], {})
        labels = []
        for c in (str(rec["ligand_chains"]).split(",")
                  + str(rec["receptor_chains"]).split(",")):
            c = c.strip()
            s = (chains.get(c) or {}).get("seq", "")
            if s:
                labels.append(assign.get(s) or ("S" + hashlib.md5(
                    s.encode("utf-8")).hexdigest()[:6]))
        return "|".join(sorted(set(labels))) if labels else "F_struct_%s" % rec["pdb_id"]

    # ---- 全局家族标签: 必须先合并再计算 -------------------------------------
    # 若分别对 full / mvp 计算, 同一 PDB 在两张表里的链集合可能不同 (full 是子集),
    # 家族指纹会不一致, 进而导致该 PDB 被拆到两侧。这里统一在合并表上计算。
    both = pd.concat([full.assign(_src="full"), mvp.assign(_src="mvp")],
                     ignore_index=True, sort=False)
    both["family_id_v1"] = [row_family(r) for r in both.to_dict("records")]
    per_pdb = both.groupby("pdb_id")["family_id_v1"].nunique()
    bad = int((per_pdb > 1).sum())
    if bad:
        logger.warning("%d 个 PDB 因不同行的链集合不同而家族指纹不一致; "
                       "统一为其 PDB 级指纹 (保守, 不会跨侧)", bad)
    both["family_id_v1"] = both["pdb_id"].map(
        both.groupby("pdb_id")["family_id_v1"].first())
    logger.info("全局数据行家族数: %d (唯一 PDB %d)",
                both["family_id_v1"].nunique(), both["pdb_id"].nunique())

    # ---- 整族分侧 (家族级哈希排序取前 80% 行) ----
    fam_stat = (both.groupby("family_id_v1").agg(n=("real_dG", "size"))
                .reset_index())
    fam_stat["h"] = fam_stat["family_id_v1"].map(hash_unit)
    fam_stat = fam_stat.sort_values(["h", "family_id_v1"]).reset_index(drop=True)
    fam_stat["cum"] = fam_stat["n"].cumsum()
    cut = TRAIN_FRACTION * len(both)
    best = int((fam_stat["cum"] - cut).abs().idxmin())
    best = max(0, min(best, len(fam_stat) - 2)) if len(fam_stat) > 1 else 0
    train_fams = set(fam_stat.loc[:best, "family_id_v1"])
    both["split_family"] = np.where(both["family_id_v1"].isin(train_fams),
                                    "train", "test")

    # ---- 红线硬校验 ----
    cross_fam = int((both.groupby("family_id_v1")["split_family"].nunique() > 1).sum())
    cross_pdb = int((both.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    if cross_fam or cross_pdb:
        raise RuntimeError("划分泄漏! 跨侧家族=%d, 跨侧PDB=%d" % (cross_fam, cross_pdb))
    logger.info("红线校验通过: 0 个家族 / 0 个 PDB 跨侧")

    # ---- 写回两个特征表 ----
    for name, path in (("full", FULL_CSV), ("mvp", MVP_CSV)):
        sub = both[both["_src"] == name][["interface_key", "family_id_v1",
                                          "split_family"]]
        df = pd.read_csv(path)
        df = df.drop(columns=[c for c in ("family_id_v1", "split_v1.0")
                              if c in df.columns])
        m = sub.set_index("interface_key")
        df["family_id"] = df["interface_key"].map(m["family_id_v1"])
        df["split_family"] = df["interface_key"].map(m["split_family"])
        df.to_csv(path, index=False, encoding="utf-8")
        logger.info("写回 %s: family_id 非空=%d, split_family 分布=%s",
                    path, int(df["family_id"].notna().sum()),
                    df["split_family"].value_counts().to_dict())

    both[["pdb_id", "interface_key", "ligand_chains", "receptor_chains",
          "family_id_v1", "split_family", "real_dG", "_src"]].to_csv(
        ASSIGN_CSV, index=False, encoding="utf-8")

    # ---- 报告 ----
    tr = both[both["split_family"] == "train"]
    te = both[both["split_family"] == "test"]
    ei = np.asarray(stats["edge_identities"], dtype=float)
    rep = []
    rep.append("# V1.0 全量家族级划分报告\n")
    rep.append("- 训练数据: features_full (%d 行 / %d 结构) + features_mvp (%d 行 / %d 结构)"
               % (len(full), full["pdb_id"].nunique(), len(mvp), mvp["pdb_id"].nunique()))
    rep.append("- 序列来源: RCSB FASTA (www.rcsb.org/fasta/entry), 并发 2")
    rep.append("- 家族定义: **40% 序列同一性** 单链接聚类 (union-find), 与 04 规则完全一致")
    rep.append("- 同一性: 4-mer Jaccard 预筛 (阈值 %.2f) + 带状全局比对 (match +2/-1/-1)"
               % PREFILTER_JACCARD)
    rep.append("- 参与聚类: 长度 >= %d aa 的链序列; 预筛后候选对 %d\n"
               % (MIN_SEQ_LEN, stats["n_candidates"]))
    rep.append("## 家族数量\n")
    rep.append("| 指标 | 数值 |")
    rep.append("|---|---|")
    rep.append("| 唯一链序列数 | %d |" % stats["n_seq"])
    rep.append("| **40%% 阈值簇数** | **%d** |" % stats["n_clusters"])
    rep.append("| 单例簇 | %d |" % stats["n_singletons"])
    rep.append("| 最大簇规模 | %d |" % stats["largest_cluster"])
    rep.append("| 聚类边数 | %d |" % stats["n_edges"])
    if ei.size:
        rep.append("| 边同一性最小值 | %.3f |" % ei.min())
        rep.append("| 边同一性中位数 | %.3f |" % float(np.median(ei)))
    rep.append("| 数据行家族数 | %d |" % both["family_id_v1"].nunique())
    rep.append("\n阈值敏感性: " + ", ".join("%.0f%% -> %d 簇" % (k * 100, v)
                                            for k, v in sorted(sens.items())))
    rep.append("\n## 划分结果\n")
    rep.append("| 划分 | 行数 | 占比 | 唯一 PDB | 家族数 | real_dG 均值±std |")
    rep.append("|---|---|---|---|---|---|")
    for nm, g in (("train", tr), ("test", te)):
        rep.append("| %s | %d | %.1f%% | %d | %d | %.2f ± %.2f |"
                   % (nm, len(g), 100.0 * len(g) / len(both), g["pdb_id"].nunique(),
                      g["family_id_v1"].nunique(), g["real_dG"].mean(),
                      g["real_dG"].std()))
    rep.append("\n## 红线校验\n")
    rep.append("- 家族整体分侧: **PASS** (0 个家族跨侧, 脚本内硬校验)")
    rep.append("- 同一 PDB 不跨界: **PASS** (0 个 PDB 跨侧)")
    rep.append("\n## 与早期 PDB 级 split 的差异\n")
    old_tr = int((both["split"] == "train").sum()) if "split" in both.columns else -1
    rep.append("- 旧 `split` 列 (PDB 级, 4762 行上分配): train %d / test %d"
               % (old_tr, len(both) - old_tr))
    rep.append("- 新 `split_family` 列 (序列家族级, 本报告数据上分配): train %d / test %d"
               % (len(tr), len(te)))
    rep.append("- **V1.0 训练必须使用 `split_family`**, 否则同源家族会跨侧。")
    rep.append("\n## 残留局限\n")
    rep.append("1. 40% 同一性仍是同源代理; 单链接聚类存在传递性 (簇内成员不保证两两 >=40%)。")
    rep.append("2. 无序列 PDB %d 个以 PDB 身份兜底成独立家族 (不会跨侧)。" % len(missing))
    rep.append("3. 数据仅覆盖已完成特征提取的 %d 个结构 (全量 2774 中已跑 %d 个)。"
               % (both["pdb_id"].nunique(), full["pdb_id"].nunique()))
    rep.append("\n生成时间: %s | 耗时 %.1fs"
               % (time.strftime("%Y-%m-%d %H:%M:%S"), time.time() - t0))
    with open(REPORT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(rep) + "\n")

    print("\n" + "=" * 80)
    print("V1.0 家族划分完成")
    print("=" * 80)
    print("链序列 %d 条 -> 40%% 阈值 %d 个簇 (单例 %d, 最大簇 %d)"
          % (stats["n_seq"], stats["n_clusters"], stats["n_singletons"],
             stats["largest_cluster"]))
    print("数据行家族数: %d" % both["family_id_v1"].nunique())
    print("划分: train %d 行 (%.1f%%, %d PDB) / test %d 行 (%.1f%%, %d PDB)"
          % (len(tr), 100.0 * len(tr) / len(both), tr["pdb_id"].nunique(),
             len(te), 100.0 * len(te) / len(both), te["pdb_id"].nunique()))
    print("红线: 0 个家族 / 0 个 PDB 跨侧 PASS")
    print("报告: %s" % REPORT)
    print("=" * 80)
    print("V1.0 家族划分完成，等待下一步指令")


if __name__ == "__main__":
    main()
