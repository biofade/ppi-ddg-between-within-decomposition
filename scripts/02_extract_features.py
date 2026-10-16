# -*- coding: utf-8 -*-
"""
02_extract_features.py -- Affinity Intelligence / 阶段二：物理特征提取层

输入 : data/benchmark.csv        (阶段一产物, 取 is_mvp==True 的行)
       data/pdb/*.pdb|*.cif      (阶段二 02a 抓取的 MVP 结构)
输出 : features/features_mvp.csv (特征主表)
       features/features_mvp_log.txt

特征口径 (全部记录在 CSV 的 feature_version 里):
  几何 (Biopython + scipy KD-tree):
    n_contact_pairs       跨链重原子距离 < 5.0 Å 的原子对数量
    n_contact_res_pairs   跨链接触的残基对数量 (去重)
    n_iface_res_lig/rec   界面残基个数 (任一重原子参与接触)
    iface_res_frac_*      界面残基占该组总残基的比例
    min_interatomic_dist  最近跨链重原子距离
    n_hbond               N/O 间距离 <= 3.5 Å 的跨链氢键 (含供体/受体元素约束)
    n_salt_bridge         带电侧链间距离 <= 4.0 Å 的盐桥 (去重残基对)
    n_hydrophobic_contact 双方均为疏水残基的接触残基对
    n_chainpair_contacts  逐链对的接触数 (取最大/最小, 表征界面分布)
  面积 (FreeSASA, Lee-Richards, 去氢去 HETATM 的统一原子集):
    dsasa_total / dsasa_lig / dsasa_rec   界面埋藏面积 (Å^2)
    dsasa_frac_rec / dsasa_frac_lig       埋藏面积占单体 SASA 比例
    sasa_complex / sasa_lig_mono / sasa_rec_mono
  序列 (Kyte-Doolittle):
    kd_hydro_*            界面残基平均疏水指数
    frac_charged_iface / frac_polar_iface / frac_hydrophobic_iface
  二级结构 (PDB HELIX/SHEET 记录, 不依赖 DSSP):
    frac_helix_iface / frac_sheet_iface
  规模:
    n_res_lig / n_res_rec / n_atoms_lig / n_atoms_rec
  能量:
    predicted_dG          PRODIGY 预测 (环境不可用时整列缺失, 并写日志)

纪律: 纯 CPU; 每个 PDB 独立 try/except, 失败即跳过并记日志, 绝不重试/修复;
      循环内周期性 gc.collect(); 超过 DEADLINE 分钟仍未跑完则自动降级取前 200 条。
"""

import gc
import json
import logging
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

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
OUT_CSV = os.path.join(FEAT_DIR, "features_mvp.csv")
OUT_JSON = os.path.join(FEAT_DIR, "features_mvp.json")
LOG_PATH = os.path.join(FEAT_DIR, "features_mvp_log.txt")

FEATURE_VERSION = "v0.1-freesasa-Biopython"
CONTACT_CUTOFF = 5.0        # Å, 跨链接触
HBOND_CUTOFF = 3.5          # Å, N/O 氢键
SALT_CUTOFF = 4.0           # Å, 盐桥
DEADLINE_MIN = 20           # 超时自动降级
FALLBACK_N = 200            # 降级后前 N 条
MEM_WARN_MB = 3000          # 内存预警阈值

os.makedirs(FEAT_DIR, exist_ok=True)

logger = logging.getLogger("features")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)

# Kyte-Doolittle 疏水指数
KD = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
    "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}
CHARGED = {"LYS", "ARG", "HIS", "ASP", "GLU"}
POSITIVE_ATOMS = {"LYS": {"NZ"}, "ARG": {"NH1", "NH2", "NE"},
                  "HIS": {"ND1", "NE2"}}
NEGATIVE_ATOMS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "TRP", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}
AA3 = set(KD)


# ---------------------------------------------------------------- 工具
def get_rss_mb():
    """无 psutil 依赖的常驻内存读数 (Windows)。"""
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
    except Exception:                                       # noqa: BLE001
        return float("nan")


def resolve_struct_path(pdb_id):
    for ext in (".pdb", ".cif", ".ent"):
        p = os.path.join(PDB_DIR, pdb_id + ext)
        if os.path.exists(p) and os.path.getsize(p) > 500:
            return p
    return None


def build_ss_map(path):
    """直接解析 PDB 的 HELIX/SHEET 固定列记录, 返回 {(chain,resseq): 1|2}。

    注意: Biopython 1.88 的 Structure 对象不暴露 helix/sheet, Residue 也没有
    is_helix()/is_sheet(), 因此只能读原始记录 (mmCIF 暂不支持, 记 0)。
    """
    ss = {}
    if not path.lower().endswith((".pdb", ".ent")):
        return ss
    try:
        with open(path, "r", errors="replace") as fh:
            for ln in fh:
                rec = ln[:6]
                if rec == "HELIX ":
                    try:
                        ch = ln[19]
                        s, e = int(ln[21:25]), int(ln[33:37])
                        if s > e:
                            s, e = e, s
                        for rs in range(s, e + 1):
                            ss[(ch, rs)] = 1
                    except (ValueError, IndexError):
                        continue
                elif rec == "SHEET ":
                    try:
                        ch = ln[21]
                        s, e = int(ln[22:26]), int(ln[33:37])
                        if s > e:
                            s, e = e, s
                        for rs in range(s, e + 1):
                            ss[(ch, rs)] = 2
                    except (ValueError, IndexError):
                        continue
    except OSError:
        return ss
    return ss


def extract_chain_data(structure, chain_ids, ss_map):
    """返回该链组的 (coords, elements, atom_names, res_keys, res_names, ss_flags)。"""
    coords, elements, anames, rkeys, rnames, ss = [], [], [], [], [], []
    present = {c.id for c in structure}
    for cid in chain_ids:
        if cid not in present:
            continue
        chain = structure[cid]
        for res in chain:
            if res.id[0] != " ":            # 跳过水/配体等 HETATM
                continue
            rname = res.get_resname().strip().upper()
            if rname not in AA3:
                continue
            ss_flag = ss_map.get((cid, res.id[1]), 0)
            for atom in res.get_atoms():
                el = (atom.element or "").strip().upper()
                if el in ("H", "D", ""):
                    continue
                coords.append(atom.get_coord())
                elements.append(el)
                anames.append(atom.get_name().strip().upper())
                rkeys.append("%s:%s%s" % (cid, res.id[1], res.id[2].strip()))
                rnames.append(rname)
                ss.append(ss_flag)
    return (np.asarray(coords, dtype=np.float64),
            np.asarray(elements, dtype=object),
            np.asarray(anames, dtype=object),
            np.asarray(rkeys, dtype=object),
            np.asarray(rnames, dtype=object),
            np.asarray(ss, dtype=np.int8))


def contact_pairs(xyz_a, xyz_b, cutoff):
    """返回跨链 (i,j) 索引对, 重原子距离 < cutoff。"""
    if len(xyz_a) == 0 or len(xyz_b) == 0:
        return np.empty((0, 2), dtype=int), float("nan")
    tree = cKDTree(xyz_a)
    pairs = tree.query_ball_point(xyz_b, cutoff)
    ii, jj = [], []
    for j, lst in enumerate(pairs):
        if lst:
            for i in lst:
                ii.append(i)
                jj.append(j)
    if not ii:
        d, _ = tree.query(xyz_b)
        return np.empty((0, 2), dtype=int), float(np.min(d)) if len(d) else float("nan")
    return np.column_stack([ii, jj]), float("nan")


def group_features(xa, ea, na, ra, rna, ssa,
                   xb, eb, nb, rb, rnb, ssb):
    """计算两个链组之间的界面特征。"""
    f = {}
    pairs, _ = contact_pairs(xa, xb, CONTACT_CUTOFF)
    n_pairs = len(pairs)
    f["n_contact_pairs"] = n_pairs
    f["n_atoms_lig"] = int(len(xa))
    f["n_atoms_rec"] = int(len(xb))
    f["n_res_lig"] = int(len(np.unique(ra)))
    f["n_res_rec"] = int(len(np.unique(rb)))

    if n_pairs == 0:
        # 无任何跨链接触: 全零 + 最近距离
        if len(xa) and len(xb):
            tree = cKDTree(xa)
            d, _ = tree.query(xb)
            f["min_interatomic_dist"] = float(np.min(d))
        else:
            f["min_interatomic_dist"] = float("nan")
        for k in ("n_contact_res_pairs", "n_iface_res_lig", "n_iface_res_rec",
                  "n_hbond", "n_salt_bridge", "n_hydrophobic_contact",
                  "iface_res_frac_lig", "iface_res_frac_rec", "kd_hydro_lig",
                  "kd_hydro_rec", "kd_hydro_all", "frac_charged_iface",
                  "frac_polar_iface", "frac_hydrophobic_iface", "frac_helix_iface",
                  "frac_sheet_iface", "contact_density", "n_chainpair_max",
                  "n_chainpair_min"):
            f[k] = 0.0
        return f

    ia, ib = pairs[:, 0], pairs[:, 1]

    # 残基对 (numpy unique(axis=0) 不支持 object 数组, 用集合去重)
    rp_u = sorted({(ra[i], rb[j]) for i, j in zip(ia, ib)})
    f["n_contact_res_pairs"] = int(len(rp_u))
    iface_lig = {ra[i] for i in ia}
    iface_rec = {rb[j] for j in ib}
    f["n_iface_res_lig"] = int(len(iface_lig))
    f["n_iface_res_rec"] = int(len(iface_rec))
    f["iface_res_frac_lig"] = float(len(iface_lig) / max(f["n_res_lig"], 1))
    f["iface_res_frac_rec"] = float(len(iface_rec) / max(f["n_res_rec"], 1))
    f["min_interatomic_dist"] = float(np.min(np.linalg.norm(
        xa[ia] - xb[ib], axis=1)))

    # 氢键: N/O 之间 <= 3.5 Å
    is_no_a = np.isin(ea, ["N", "O"])
    is_no_b = np.isin(eb, ["N", "O"])
    ia_no = np.where(is_no_a)[0]
    ib_no = np.where(is_no_b)[0]
    hb_pairs, _ = contact_pairs(xa[is_no_a], xb[is_no_b], HBOND_CUTOFF)
    if len(hb_pairs):
        # 去重到残基对, 避免同一残基对多原子重复计数
        f["n_hbond"] = int(len({(ra[ia_no[i]], rb[ib_no[j]])
                                for i, j in hb_pairs}))
    else:
        f["n_hbond"] = 0

    # 盐桥: 带电侧链原子对 <= 4.0 Å
    pos_a = np.array([(rna[i] in POSITIVE_ATOMS and na[i] in POSITIVE_ATOMS[rna[i]])
                      for i in range(len(xa))])
    neg_a = np.array([(rna[i] in NEGATIVE_ATOMS and na[i] in NEGATIVE_ATOMS[rna[i]])
                      for i in range(len(xa))])
    pos_b = np.array([(rnb[i] in POSITIVE_ATOMS and nb[i] in POSITIVE_ATOMS[rnb[i]])
                      for i in range(len(xb))])
    neg_b = np.array([(rnb[i] in NEGATIVE_ATOMS and nb[i] in NEGATIVE_ATOMS[rnb[i]])
                      for i in range(len(xb))])
    salt_res = []
    for sel_a, sel_b in ((pos_a, neg_b), (neg_a, pos_b)):
        if sel_a.any() and sel_b.any():
            sp, _ = contact_pairs(xa[sel_a], xb[sel_b], SALT_CUTOFF)
            for i, j in sp:
                salt_res.append((ra[np.where(sel_a)[0][i]], rb[np.where(sel_b)[0][j]]))
    f["n_salt_bridge"] = int(len(set(salt_res)))

    # 疏水接触残基对
    fa = {x: i for i, x in enumerate(ra)}
    fb = {x: i for i, x in enumerate(rb)}
    hyd = set()
    for i, j in rp_u:
        if rna[fa[i]] in HYDROPHOBIC and rnb[fb[j]] in HYDROPHOBIC:
            hyd.add((i, j))
    f["n_hydrophobic_contact"] = int(len(hyd))

    # 界面残基序列/结构特征 (按残基键聚合, 不能用按原子的掩码)
    def res_level(rkeys, rnames, ssflags):
        d = {}
        for k, n, s in zip(rkeys, rnames, ssflags):
            if k not in d:
                d[k] = (n, int(s))
        return d

    def iface_stats(per_res, iface_keys):
        names, ss = [], []
        for k in iface_keys:
            if k in per_res:
                nm, sf = per_res[k]
                names.append(nm)
                ss.append(sf)
        if not names:
            return {"kd": float("nan"), "charged": 0.0, "polar": 0.0,
                    "hydro": 0.0, "helix": 0.0, "sheet": 0.0}
        kd_vals = np.array([KD.get(x, np.nan) for x in names], dtype=float)
        kd_mean = float(np.nan) if np.all(np.isnan(kd_vals)) else float(np.nanmean(kd_vals))
        ss_arr = np.asarray(ss)
        return {
            "kd": kd_mean,
            "charged": float(np.mean([x in CHARGED for x in names])),
            "polar": float(np.mean([x in POLAR for x in names])),
            "hydro": float(np.mean([x in HYDROPHOBIC for x in names])),
            "helix": float(np.mean(ss_arr == 1)),
            "sheet": float(np.mean(ss_arr == 2)),
        }

    sa = iface_stats(res_level(ra, rna, ssa), iface_lig)
    sb = iface_stats(res_level(rb, rnb, ssb), iface_rec)
    f["kd_hydro_lig"] = sa["kd"]
    f["kd_hydro_rec"] = sb["kd"]
    both = [v for v in (sa["kd"], sb["kd"]) if v == v]
    f["kd_hydro_all"] = float(np.mean(both)) if both else float("nan")
    f["frac_charged_iface"] = 0.5 * (sa["charged"] + sb["charged"])
    f["frac_polar_iface"] = 0.5 * (sa["polar"] + sb["polar"])
    f["frac_hydrophobic_iface"] = 0.5 * (sa["hydro"] + sb["hydro"])
    f["frac_helix_iface"] = 0.5 * (sa["helix"] + sb["helix"])
    f["frac_sheet_iface"] = 0.5 * (sa["sheet"] + sb["sheet"])

    # 界面分布: 逐链对的接触数 / 氢键数 (链对键为排序后的 (受体链, 配体链))
    pair_counts = {}
    for i, j in zip(ra[ia], rb[ib]):
        key = tuple(sorted((i.split(":")[0], j.split(":")[0])))
        pair_counts[key] = pair_counts.get(key, 0) + 1
    hb_pair_counts = {}
    if len(hb_pairs):
        for i, j in hb_pairs:
            key = tuple(sorted((ra[ia_no[i]].split(":")[0],
                                rb[ib_no[j]].split(":")[0])))
            hb_pair_counts[key] = hb_pair_counts.get(key, 0) + 1
    vals = list(pair_counts.values()) or [0]
    f["n_chainpair_max"] = int(max(vals))
    f["n_chainpair_min"] = int(min(vals))
    f["n_chainpair_contacts"] = len(vals)

    # 接触密度: 每界面残基平均接触原子对
    f["contact_density"] = float(n_pairs / max(f["n_iface_res_lig"] + f["n_iface_res_rec"], 1))
    f["_pair_contacts"] = {k: int(v) for k, v in pair_counts.items()}
    f["_pair_hbonds"] = {k: int(v) for k, v in hb_pair_counts.items()}
    return f


def _write_subset(model, out_path, chains=None):
    """写出 protein-only (可选链过滤) 的临时 PDB。"""
    from Bio.PDB import PDBIO, Select

    class Pick(Select):
        def accept_residue(self, residue):
            ok = residue.id[0] == " "
            if chains is not None:
                ok = ok and residue.get_parent().id in chains
            return ok

    io = PDBIO()
    io.set_structure(model)
    io.save(out_path, Pick())
    return out_path


def sasa_features(model, lig_chains, rec_chains, tmpdir, pair_contacts=None,
                  pair_hbonds=None, max_pair_sasa=40):
    """FreeSASA: 统一原子集下算单体/复合物 SASA + **逐链对界面埋藏面积**。

    V0.3 新增: 对每个真实存在接触的 (配体链, 受体链) 对单独算
        dsasa_pair = SASA(L单体) + SASA(R单体) - SASA(L+R 双链复合物)
    给出"最弱界面"物理量, 让模型关注界面强度而非蛋白体型。
    计时纪律: 单体 SASA 按链缓存复用; 链对数超过 max_pair_sasa 时按接触数
    降序只算前 max_pair_sasa 个, 并以 pair_sasa_complete 标注是否完整。
    """
    import freesasa

    freesasa.setVerbosity(freesasa.silent)
    model_chains = {c.id for c in model}

    out = {}
    all_path = _write_subset(model, os.path.join(tmpdir, "all.pdb"))
    st_all = freesasa.Structure(all_path)
    res_all = freesasa.calc(st_all)
    out["sasa_complex"] = float(res_all.totalArea())
    areas = res_all.residueAreas()

    def chain_area(ch):
        """该链在复合物中的 SASA (从复合物单次计算的 residueAreas 取得)。"""
        if ch in areas:
            return float(sum(r.total for r in areas[ch].values()))
        p = _write_subset(model, os.path.join(tmpdir, "ch_%s.pdb" % ch), chains={ch})
        return float(freesasa.calc(freesasa.Structure(p)).totalArea())

    # 单体 SASA 按链缓存 (同一链在多个链对里复用, 避免重复计算)
    mono_cache = {}

    def mono_of(ch):
        if ch not in mono_cache:
            p = _write_subset(model, os.path.join(tmpdir, "mono_%s.pdb" % ch),
                              chains={ch})
            mono_cache[ch] = float(freesasa.calc(freesasa.Structure(p)).totalArea())
        return mono_cache[ch]

    def pair_sasa(la, rb_):
        """两链双元复合物 SASA; 失败返回 None (不强行修复)。"""
        p = _write_subset(model, os.path.join(tmpdir, "pair_%s_%s.pdb" % (la, rb_)),
                          chains={la, rb_})
        try:
            return float(freesasa.calc(freesasa.Structure(p)).totalArea())
        except Exception as exc:                                # noqa: BLE001
            logger.debug("链对 SASA 失败 %s-%s: %s", la, rb_, exc)
            return None

    def group_mono(chains, tag):
        p = _write_subset(model, os.path.join(tmpdir, "%s.pdb" % tag), chains=set(chains))
        return float(freesasa.calc(freesasa.Structure(p)).totalArea())

    try:
        out["sasa_lig_mono"] = group_mono(lig_chains, "lig")
    except Exception as exc:                                    # noqa: BLE001
        logger.debug("单体 SASA 计算失败 (lig): %s", exc)
        out["sasa_lig_mono"] = float("nan")
    try:
        out["sasa_rec_mono"] = group_mono(rec_chains, "rec")
    except Exception as exc:                                    # noqa: BLE001
        logger.debug("单体 SASA 计算失败 (rec): %s", exc)
        out["sasa_rec_mono"] = float("nan")

    lig_present = [c for c in lig_chains if c in model_chains]
    rec_present = [c for c in rec_chains if c in model_chains]
    lig_in = float(sum(chain_area(c) for c in lig_present))
    rec_in = float(sum(chain_area(c) for c in rec_present))
    out["sasa_lig_in_complex"] = lig_in
    out["sasa_rec_in_complex"] = rec_in
    out["dsasa_lig"] = out["sasa_lig_mono"] - lig_in
    out["dsasa_rec"] = out["sasa_rec_mono"] - rec_in
    out["dsasa_total"] = out["dsasa_lig"] + out["dsasa_rec"]
    out["dsasa_frac_lig"] = out["dsasa_lig"] / max(out["sasa_lig_mono"], 1e-6)
    out["dsasa_frac_rec"] = out["dsasa_rec"] / max(out["sasa_rec_mono"], 1e-6)

    # ---- V0.3: 逐链对 ΔSASA ----
    pair_contacts = pair_contacts or {}
    # 只对"有接触的链对"算埋藏面积 (无接触的链对 ΔSASA≈0, 不占算力)
    cand = []
    for (ca, cb), ncnt in pair_contacts.items():
        # 键是排序后的链 ID, 需要判断哪条属于配体组
        if ca in lig_present and cb in rec_present:
            cand.append((ca, cb, ncnt))
        elif cb in lig_present and ca in rec_present:
            cand.append((cb, ca, ncnt))
    cand.sort(key=lambda x: (-x[2], x[0], x[1]))     # 接触多的优先, 确定性排序
    truncated = len(cand) > max_pair_sasa
    cand = cand[:max_pair_sasa]

    pair_dsasa = {}
    for la, rb_, _ncnt in cand:
        try:
            s_pair = pair_sasa(la, rb_)
            if s_pair is None:
                continue
            pair_dsasa["%s-%s" % (la, rb_)] = float(
                mono_of(la) + mono_of(rb_) - s_pair)
        except Exception as exc:                                # noqa: BLE001
            logger.debug("链对 ΔSASA 失败 %s-%s: %s", la, rb_, exc)

    dv = [v for v in pair_dsasa.values() if v > 0]
    out["n_chainpairs_with_contact"] = len(cand)
    out["n_chainpairs_sasa_ok"] = len(pair_dsasa)
    out["pair_sasa_complete"] = int(not truncated)
    if dv:
        out["weakest_chainpair_dsasa"] = float(min(dv))
        out["strongest_chainpair_dsasa"] = float(max(dv))
        out["mean_chainpair_dsasa"] = float(np.mean(dv))
        out["std_chainpair_dsasa"] = float(np.std(dv))
        out["weakest_over_total_dsasa"] = float(min(dv) / max(out["dsasa_total"], 1e-6))
    else:
        for k in ("weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
                  "mean_chainpair_dsasa", "std_chainpair_dsasa",
                  "weakest_over_total_dsasa"):
            out[k] = float("nan")
    return out


def try_prodigy(pdb_path, lig_chains, rec_chains):
    """PRODIGY 预测 ΔG; 环境不可用时返回 (None, 原因)。"""
    try:
        from prodigy_prot.modules import predict_IC as _pic      # noqa: F401
    except Exception:                                            # noqa: BLE001
        pass
    try:
        from prodigy_prot.predict_IC import Predictor
        pred = Predictor()
        temp = {"A": "A", "B": "B"}
        pred.predict(pdb_path, temp, chain_ids=[",".join(rec_chains),
                                                ",".join(lig_chains)])
        return {"predicted_dG": getattr(pred, "dG", None),
                "prodigy_ic": getattr(pred, "ic", None)}, None
    except Exception as exc:                                     # noqa: BLE001
        return None, "%s: %s" % (type(exc).__name__, str(exc)[:120])


def derive_normalized_features(df):
    """V0.2 新增: 体积/规模归一化特征 + 标签相关派生量。

    重要区分 (避免目标泄漏):
      * 无标签归一化特征 (可安全入模): interface_residue_ratio_* 等 —— 分子分母都来自结构。
      * **标签派生量** (分子含 real_dG, 入模即泄漏): dG_per_contact / dG_per_dSASA。
        仅保留作逐界面"能量效率"诊断, 列名以 _leaky 结尾, 模型脚本必须排除。
    """
    out = df.copy()
    num = out.select_dtypes(include=[np.number])

    def safe_div(a, b):
        a = pd.to_numeric(a, errors="coerce")
        b = pd.to_numeric(b, errors="coerce")
        return a / b.replace(0, np.nan)

    # ---- 无标签: 界面残基比例 (用户的 interface_residue_ratio) ----
    out["interface_residue_ratio_lig"] = safe_div(out["n_iface_res_lig"], out["n_res_lig"])
    out["interface_residue_ratio_rec"] = safe_div(out["n_iface_res_rec"], out["n_res_rec"])
    out["interface_residue_ratio"] = safe_div(
        num.get("n_iface_res_lig", 0) + num.get("n_iface_res_rec", 0),
        num.get("n_res_lig", 0) + num.get("n_res_rec", 0))

    # ---- 无标签: 体积/规模归一化的界面强度 (替代原始规模量的候选) ----
    out["contact_per_iface_res"] = safe_div(out["n_contact_pairs"],
                                            num.get("n_iface_res_lig", 0)
                                            + num.get("n_iface_res_rec", 0))
    out["contacts_per_res_total"] = safe_div(out["n_contact_pairs"],
                                             num.get("n_res_lig", 0) + num.get("n_res_rec", 0))
    out["dsasa_per_iface_res"] = safe_div(out["dsasa_total"],
                                          num.get("n_iface_res_lig", 0)
                                          + num.get("n_iface_res_rec", 0))
    out["dsasa_per_contact"] = safe_div(out["dsasa_total"], out["n_contact_pairs"])
    out["hbond_per_contact"] = safe_div(out["n_hbond"], out["n_contact_pairs"])
    out["hbond_per_iface_res"] = safe_div(out["n_hbond"],
                                          num.get("n_iface_res_lig", 0)
                                          + num.get("n_iface_res_rec", 0))
    out["saltbridge_per_iface_res"] = safe_div(out["n_salt_bridge"],
                                               num.get("n_iface_res_lig", 0)
                                               + num.get("n_iface_res_rec", 0))
    out["iface_res_per_contact"] = safe_div(
        num.get("n_iface_res_lig", 0) + num.get("n_iface_res_rec", 0),
        out["n_contact_pairs"])
    out["dsasa_frac_total"] = safe_div(
        out["dsasa_total"], num.get("sasa_complex", 0) + out["dsasa_total"])

    # ---- V0.3: 链对级物理量 (用"每链对/每接触"的强度替代绝对规模) ----
    ncp = num.get("n_chainpair_contacts", 1)
    out["chainpair_contact_density"] = safe_div(out["n_contact_pairs"], ncp)
    out["chainpair_dsasa_mean_from_total"] = safe_div(out["dsasa_total"], ncp)
    # 最弱界面 vs 平均界面: 衡量结合是否由单一界面主导
    out["weakest_to_mean_chainpair_dsasa"] = safe_div(
        out["weakest_chainpair_dsasa"], out["mean_chainpair_dsasa"])
    out["mean_chainpair_dsasa_per_contact"] = safe_div(
        out["mean_chainpair_dsasa"], out["n_contact_pairs"])
    # 单个接触的平均埋藏面积 (界面紧密程度的形状代理)
    out["dsasa_per_contactpair"] = safe_div(out["dsasa_total"], out["n_contact_pairs"])
    out["hbond_density_periface_area"] = safe_div(
        out["n_hbond"], out["dsasa_total"])
    out["saltbridge_density_periface_area"] = safe_div(
        out["n_salt_bridge"], out["dsasa_total"])

    # ---- 标签派生量 (仅诊断, 禁止入模) ----
    out["dG_per_contact__leaky"] = safe_div(out["real_dG"], out["n_contact_pairs"])
    out["dG_per_dSASA__leaky"] = safe_div(out["real_dG"], out["dsasa_total"])
    out["dG_per_iface_res__leaky"] = safe_div(
        out["real_dG"], num.get("n_iface_res_lig", 0) + num.get("n_iface_res_rec", 0))

    return out.replace([np.inf, -np.inf], np.nan)


# ---------------------------------------------------------------- 主流程
def main():
    t0 = time.time()
    from Bio.PDB import PDBParser, MMCIFParser

    logger.info("特征版本: %s", FEATURE_VERSION)
    bench = pd.read_csv(BENCH, dtype={"pdb_id": str})
    mvp = bench[bench["is_mvp"] == True].copy()                  # noqa: E712
    logger.info("benchmark 总行 %d, MVP 行 %d, 唯一结构 %d",
                len(bench), len(mvp), mvp["pdb_id"].nunique())

    # 结构可用性过滤
    mvp["struct_path"] = mvp["pdb_id"].map(resolve_struct_path)
    n_missing = int(mvp["struct_path"].isna().sum())
    if n_missing:
        logger.warning("缺少结构文件的行: %d (将计入解析失败)", n_missing)
    todo = mvp[mvp["struct_path"].notna()].copy()
    logger.info("待提取: %d 行", len(todo))

    # 以结构为单位组织 (同一 PDB 只解析一次)
    by_pdb = {}
    for rec in todo.to_dict("records"):
        by_pdb.setdefault(rec["pdb_id"], []).append(rec)
    pdb_ids = sorted(by_pdb)
    logger.info("唯一结构数: %d", len(pdb_ids))

    # ---- PRODIGY 可用性一次性排查 (失败即整列跳过, 绝不死磕) ----
    prodigy_ok = False
    prodigy_reason = ""
    probe = todo.iloc[0]
    out_p, reason = try_prodigy(probe["struct_path"],
                                str(probe["ligand_chains"]).split(","),
                                str(probe["receptor_chains"]).split(","))
    if out_p and out_p.get("predicted_dG") is not None:
        prodigy_ok = True
        logger.info("PRODIGY 可用, predicted_dG 已纳入特征")
    else:
        prodigy_reason = reason
        logger.warning("PRODIGY不可用，使用Biopython基础特征 (原因: %s)", reason)

    # ---- 逐结构提取 (超时/超内存则降级为前 200 条, 已算结果不丢弃) ----
    parser_pdb, parser_cif = PDBParser(QUIET=True), MMCIFParser(QUIET=True)
    rows, failed, downgraded, stopped_at = [], [], False, None

    def extract_range(pdb_subset, deadline):
        """返回 (是否被中断)。已算结果追加到 rows。"""
        for k, pid in enumerate(pdb_subset):
            if time.time() > deadline:
                logger.warning("已超过 %d 分钟: 在结构 %d/%d (%s) 处中断",
                               DEADLINE_MIN, k, len(pdb_subset), pid)
                return True
            rss = get_rss_mb()
            if rss == rss and rss > MEM_WARN_MB:
                logger.warning("内存 %.0f MB 超阈值 %d MB", rss, MEM_WARN_MB)
                return True

            rec0 = by_pdb[pid][0]
            path = rec0["struct_path"]
            try:
                parser = parser_cif if path.lower().endswith(".cif") else parser_pdb
                structure = parser.get_structure(pid, path)
                model = structure[0]
                chain_ids = [c.id for c in model]
                ss_map = build_ss_map(path)
            except Exception as exc:                            # noqa: BLE001
                logger.warning("解析失败(跳过) %s: %s: %s",
                               pid, type(exc).__name__, str(exc)[:80])
                failed.extend([(r["interface_key"], "parse_error")
                               for r in by_pdb[pid]])
                gc.collect()
                continue

            for rec in by_pdb[pid]:
                lig = [c.strip() for c in str(rec["ligand_chains"]).split(",") if c.strip()]
                rcp = [c.strip() for c in str(rec["receptor_chains"]).split(",") if c.strip()]
                try:
                    # 链级容错: 标注链在结构文件中缺失时, 用可得的链计算并标注 (不强行修复)
                    lig_p = [c for c in lig if c in chain_ids]
                    rcp_p = [c for c in rcp if c in chain_ids]
                    n_missing_chain = len(lig) + len(rcp) - len(lig_p) - len(rcp_p)
                    if not lig_p or not rcp_p:
                        raise ValueError("链组完全缺失 (lig=%s rec=%s, 结构含 %s)"
                                         % (lig, rcp, chain_ids))

                    xa, ea, na, ra, rna, ssa = extract_chain_data(model, lig_p, ss_map)
                    xb, eb, nb, rb, rnb, ssb = extract_chain_data(model, rcp_p, ss_map)
                    if len(xa) == 0 or len(xb) == 0:
                        raise ValueError("链组无有效蛋白原子")

                    feat = {"pdb_id": pid,
                            "interface_key": rec["interface_key"],
                            "ligand_chains": rec["ligand_chains"],
                            "receptor_chains": rec["receptor_chains"],
                            "real_dG": float(rec["real_dG"]),
                            "split": rec["split"],
                            "split_v0.1": rec.get("split_v0.1", rec["split"]),
                            "family_id": rec.get("family_id", ""),
                            "subgroup": rec.get("subgroup", ""),
                            "source_dataset": rec.get("source_dataset", ""),
                            "n_missing_chain": int(n_missing_chain),
                            "chain_map_ok": int(n_missing_chain == 0)}
                    gf = group_features(xa, ea, na, ra, rna, ssa,
                                        xb, eb, nb, rb, rnb, ssb)
                    pair_contacts = gf.pop("_pair_contacts", {})
                    pair_hbonds = gf.pop("_pair_hbonds", {})
                    feat.update(gf)
                    # 数据质量标记: 配体链组 == 受体链组 时, 接触/面积是"同一条链与
                    # 自己"的比较 (自界面重复计数), 如 1TZN:A|A 的 102002 个接触对。
                    # 不删除该行, 而是给模型一个可隔离它的显式标记。
                    feat["self_interface"] = int(set(lig_p) == set(rcp_p))
                    with tempfile.TemporaryDirectory() as td:
                        feat.update(sasa_features(model, lig_p, rcp_p, td,
                                                  pair_contacts=pair_contacts,
                                                  pair_hbonds=pair_hbonds))

                    if prodigy_ok:
                        p, _rsn = try_prodigy(path, lig, rcp)
                        feat["predicted_dG"] = (p or {}).get("predicted_dG")
                        feat["prodigy_ic"] = (p or {}).get("prodigy_ic")
                    else:
                        feat["predicted_dG"] = np.nan
                        feat["prodigy_ic"] = np.nan

                    feat["feature_version"] = FEATURE_VERSION
                    rows.append(feat)
                except Exception as exc:                        # noqa: BLE001
                    logger.warning("特征失败(跳过) %s %s|%s: %s: %s",
                                   pid, rec["ligand_chains"], rec["receptor_chains"],
                                   type(exc).__name__, str(exc)[:90])
                    failed.append((rec["interface_key"], "%s" % type(exc).__name__))

            del structure, model
            gc.collect()
            if (k + 1) % 25 == 0 or (k + 1) == len(pdb_subset):
                logger.info("进度 %d/%d 成功=%d 失败=%d 内存=%.0fMB 用时=%.0fs",
                            k + 1, len(pdb_subset), len(rows), len(failed),
                            get_rss_mb(), time.time() - t0)
        return False

    deadline = t0 + DEADLINE_MIN * 60
    interrupted = extract_range(pdb_ids, deadline)

    if interrupted:
        done = {r["pdb_id"] for r in rows}
        remaining = [p for p in pdb_ids if p not in done]
        if len(remaining) > FALLBACK_N:
            downgraded = True
            keep = remaining[:FALLBACK_N]
            logger.warning("降级: 仅继续前 %d 条剩余结构 (原剩余 %d 条)",
                           FALLBACK_N, len(remaining))
            stopped_at = len(pdb_ids) - len(remaining) + FALLBACK_N
            extract_range(keep, t0 + (DEADLINE_MIN + 15) * 60)
        else:
            logger.info("剩余 %d 条 <= 降级阈值 %d, 继续跑完", len(remaining), FALLBACK_N)
            extract_range(remaining, t0 + (DEADLINE_MIN + 15) * 60)

    # ---- 落盘 ----
    if not rows:
        raise RuntimeError("没有任何特征被成功提取")
    df = pd.DataFrame(rows)
    df = derive_normalized_features(df)          # V0.2: 归一化特征 + 标签派生量
    if not prodigy_ok:
        logger.info("predicted_dG 整列缺失 (PRODIGY不可用，使用Biopython基础特征)")
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, ensure_ascii=False, default=str)

    total_attempt = len(todo) + n_missing
    rate = len(df) / max(total_attempt, 1) * 100
    summary = {
        "attempted_rows": total_attempt,
        "success_rows": len(df),
        "failed_rows": len(failed) + n_missing,
        "missing_structures": n_missing,
        "success_rate": round(rate, 2),
        "unique_pdb_ok": int(df["pdb_id"].nunique()),
        "prodigy_available": prodigy_ok,
        "prodigy_reason": prodigy_reason if not prodigy_ok else "",
        "downgraded": downgraded,
        "feature_version": FEATURE_VERSION,
        "elapsed_s": round(time.time() - t0, 1),
        "rss_mb": round(get_rss_mb(), 1),
    }
    with open(os.path.join(FEAT_DIR, "features_mvp_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    if failed:
        fail_stats = pd.Series([f[1] for f in failed]).value_counts().to_dict()
    else:
        fail_stats = {}

    # ---- 报告 ----
    print("\n" + "=" * 80)
    print("阶段二完成")
    print("=" * 80)
    print("特征表            : %s" % OUT_CSV)
    print("尝试提取          : %d 行 -> 成功 %d 行, 跳过 %d 行, 成功率 %.1f%%"
          % (total_attempt, len(df), total_attempt - len(df), rate))
    print("跳过原因分布      : %s" % (fail_stats or "无"))
    print("缺少结构文件      : %d" % n_missing)
    print("成功覆盖唯一结构  : %d / %d" % (summary["unique_pdb_ok"], len(pdb_ids)))
    print("PRODIGY           : %s"
          % ("可用, predicted_dG 已纳入" if prodigy_ok
             else "不可用 -> 已跳过, 仅用 Biopython 几何/序列特征 (原因: %s)" % prodigy_reason))
    print("降级触发          : %s" % ("是 (超时/内存)" if downgraded else "否"))
    print("峰值内存          : %.0f MB" % summary["rss_mb"])
    norm_cols = [c for c in df.columns if c.startswith(("interface_residue_ratio",
                                                        "contact_per", "contacts_per",
                                                        "dsasa_per", "hbond_per",
                                                        "saltbridge_per", "iface_res_per",
                                                        "dsasa_frac_total", "chainpair_",
                                                        "weakest_", "strongest_",
                                                        "mean_chainpair",
                                                        "std_chainpair"))]
    leaky_cols = [c for c in df.columns if c.endswith("__leaky")]
    print("归一化特征 (V0.2)  : %d 个 -> %s" % (len(norm_cols), ", ".join(norm_cols)))
    print("标签派生量 (禁止入模): %d 个 -> %s" % (len(leaky_cols), ", ".join(leaky_cols)))
    meta_cols = ("pdb_id", "interface_key", "ligand_chains", "receptor_chains",
                 "real_dG", "split", "split_v0.1", "family_id", "subgroup",
                 "source_dataset", "feature_version")
    feat_cols = [c for c in df.columns if c not in meta_cols]
    print("特征数            : %d 列 (含 %d 个归一化, 不含 %d 个标签派生量)"
          % (len(feat_cols), len(norm_cols), len(leaky_cols)))
    print("耗时              : %.1fs" % summary["elapsed_s"])

    feat_cols = [c for c in df.columns if c not in
                 ("pdb_id", "interface_key", "ligand_chains", "receptor_chains",
                  "real_dG", "split", "subgroup", "source_dataset", "feature_version")]
    print("\n特征列: %s" % ", ".join(feat_cols))
    print("\n前 5 行特征数据 (核心列):")
    show = ["pdb_id", "ligand_chains", "receptor_chains", "real_dG",
            "n_contact_pairs", "n_contact_res_pairs", "n_hbond", "n_salt_bridge",
            "dsasa_total", "kd_hydro_all"]
    print(df[[c for c in show if c in df.columns]].head(5).to_string(index=False))
    print("\n特征描述统计 (节选):")
    print(df[["n_contact_pairs", "n_hbond", "n_salt_bridge", "dsasa_total",
              "kd_hydro_all", "iface_res_frac_lig", "frac_hydrophobic_iface"]]
          .describe().round(3).to_string())
    print("\n日志: %s" % LOG_PATH)
    print("=" * 80)
    print("阶段二完成，等待下一步指令")


if __name__ == "__main__":
    main()
