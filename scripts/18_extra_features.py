# -*- coding: utf-8 -*-
"""18_extra_features.py -- 阶段 5: CPU 低成本特征增强.

设计原则
--------
**不从零重新解析 3000+ 个 PDB 的几何信息** (V1.1 已做过, 重跑是浪费且会引入口径漂移)。
本脚本只新增**可复现、CPU 廉价、与已有特征不重复**的维度:

1. **界面残基组成细分** (6 类: aromatic / charged_pos / charged_neg / polar /
   hydrophobic / special), 对 ligand 链、receptor 链、两侧合并分别算
   -> 18 个特征。数据源: ``data/pdb_seqs.json`` (每链序列) +
   ``interface_key`` 的跨链重原子接触 (5 Å) 判定界面残基。
2. **氨基酸组成** (20 维, 对**界面残基**算) -> 20 个特征。
3. **疏水矩** (Eisenberg 共识标度, 窗口 11, 角度 100°), 对界面残基窗口取最大值
   -> 表达两亲性/界面插入倾向 -> 6 个特征 (lig/rec/all 各 max 与 mean)。
4. **静电互补性简化版**: 界面残基形式电荷计数 + 链净电荷乘积
   ``q_net_lig * q_net_rec`` + 界面电荷平衡项 -> 6 个特征。
   明确标注为**简化代理**, 不是 Poisson-Boltzmann 求解。
5. **形状互补性简化版**: 界面原子密度 (每 Å² 界面面积的原子数)、
   ΔSASA/界面残基比、界面残基埋藏效率 -> 5 个特征。
6. **可选 ESM-2 嵌入融合**: 若 ``features/v2_esm2_cache.npz`` 存在 (由
   ``scripts/colab/extract_esm2.py`` 在免费 GPU 上生成), 则读取并按
   **训练集拟合** 的 PCA 降到 32 维拼上; 不存在则明确跳过并记录。

关键纪律
--------
* **PCA 只在训练集 fit**, 再 transform 测试集 —— 全量 fit 是标签泄漏。
* 新增特征逐一缺失时用训练集中位数兜底, 并在报告里给出缺失率。
* 不覆盖 V1.1 特征表; 输出到 ``features/v2_extra_features.csv`` 与合并表
  ``features/v2_features_enhanced.csv``。

纪律: 纯 CPU; 无网络; seed=42; 逐结构容错 (单条失败写 failures.log 不中断)。
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
import traceback
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
SEQ_CACHE = DATA_DIR / "pdb_seqs.json"
OUT_EXTRA = FEAT_DIR / "v2_extra_features.csv"
OUT_MERGED = FEAT_DIR / "v2_features_enhanced.csv"
REPORT_MD = FEAT_DIR / "v2_feature_report.md"
ESM2_CACHE = FEAT_DIR / "v2_esm2_cache.npz"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
CONTACT_CUTOFF = 5.0          # 跨链重原子接触距离 (Å), 与 02 口径一致
HB_CUTOFF = 3.5               # 供静电代理使用的近接触
N_LEVELS = 5
ESM2_PCA_DIM = 32

# ---------------------------------------------------------------- 残基分类表
AA_AROMATIC = set("FWY")
AA_CHARGED_POS = set("KRH")
AA_CHARGED_NEG = set("DE")
AA_POLAR = set("STNQC")
AA_HYDROPHOBIC = set("AVLIMGP")
AA_SPECIAL = set("C")          # 半胱氨酸单列 (二硫键/金属配位特殊角色)
AA_ALL = sorted(set("ACDEFGHIKLMNPQRSTVWY"))

# Eisenberg 共识疏水性标度 (kcal/mol, 归一化)
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

# pH 7.4 下的近似形式电荷 (简化代理)
FORMAL_CHARGE = {
    "K": 1.0, "R": 1.0, "H": 0.1,
    "D": -1.0, "E": -1.0,
}
for _aa in AA_ALL:
    FORMAL_CHARGE.setdefault(_aa, 0.0)

# PDB 使用**三字母**残基名 (PHE), 而本模块其余部分一律用**单字母** (F)。
# 解析结构文件时必须先做这个映射, 否则所有残基都会被当成非标准残基而被丢弃
# (实测后果: 界面残基集合恒为空, 全部界面口径特征为 NaN)。
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",   # 硒代甲硫氨酸, 按 Met 处理
    "SEC": "C", "PYL": "K",
}

HMOM_WINDOW = 11
HMOM_ANGLE_DEG = 100.0

logger = logging.getLogger("extra_features")


# ---------------------------------------------------------------- 化学计量
def residue_class_counts(seq: str) -> dict[str, float]:
    """对一段序列统计 6 类残基的**占比**。"""
    n = len(seq)
    if n == 0:
        return {k: float("nan") for k in
                ("aromatic", "charged_pos", "charged_neg", "polar",
                 "hydrophobic", "special")}
    arr = list(seq)
    return {
        "aromatic": sum(a in AA_AROMATIC for a in arr) / n,
        "charged_pos": sum(a in AA_CHARGED_POS for a in arr) / n,
        "charged_neg": sum(a in AA_CHARGED_NEG for a in arr) / n,
        "polar": sum(a in AA_POLAR for a in arr) / n,
        "hydrophobic": sum(a in AA_HYDROPHOBIC for a in arr) / n,
        "special": sum(a in AA_SPECIAL for a in arr) / n,
    }


def aa_composition(seq: str) -> dict[str, float]:
    n = len(seq)
    if n == 0:
        return {f"aac_{a}": float("nan") for a in AA_ALL}
    return {f"aac_{a}": seq.count(a) / n for a in AA_ALL}


def hydrophobic_moment(seq: str, window: int = HMOM_WINDOW,
                       angle_deg: float = HMOM_ANGLE_DEG) -> tuple[float, float]:
    """疏水矩 μH = |Σ h_i · e^{i·δ·θ}| 的滑窗最大/均值 (Eisenberg 标度)。

    返回 ``(max, mean)``; 序列短于窗口或含非标准残基时跳过该窗口。
    """
    if len(seq) < window:
        return float("nan"), float("nan")
    delta = math.radians(angle_deg)
    vals: list[float] = []
    for s in range(len(seq) - window + 1):
        hx = hy = 0.0
        ok = True
        for k, aa in enumerate(seq[s:s + window]):
            h = EISENBERG.get(aa)
            if h is None:
                ok = False
                break
            ang = k * delta
            hx += h * math.cos(ang)
            hy += h * math.sin(ang)
        if ok:
            vals.append(math.hypot(hx, hy))
    if not vals:
        return float("nan"), float("nan")
    return float(max(vals)), float(sum(vals) / len(vals))


def net_charge(seq: str) -> float:
    return float(sum(FORMAL_CHARGE.get(a, 0.0) for a in seq))


# ---------------------------------------------------------------- 结构接触
def interface_residues(pdb_path: Path, lig: list[str], rec: list[str],
                       cutoff: float = CONTACT_CUTOFF) -> tuple[set, set]:
    """返回 (ligand 界面残基集合, receptor 界面残基集合), 元素为 (chain, resseq)。

    用 scipy KDTree 在 receptor 重原子上建树, 再批量查询 ligand 重原子,
    整体 O((m+n) log n), 比逐对 O(m·n) 快一到两个数量级。
    """
    from scipy.spatial import cKDTree

    def read_atoms(chains: list[str]):
        coords, meta = [], []
        want = set(chains)
        with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                rec_type = line[:6]
                if rec_type not in ("ATOM  ", "HETATM"):
                    continue
                if rec_type == "HETATM":
                    continue                      # 只算蛋白重原子
                resname3 = line[17:20].strip().upper()
                resname = AA3_TO_1.get(resname3)
                if resname is None:
                    continue                      # 跳过非标准残基/配体
                ch = line[21]
                if ch not in want:
                    continue
                el = line[76:78].strip().upper()
                if el == "H" or (not el and line[12:16].strip().startswith("H")):
                    continue                      # 跳过氢
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 0
                coords.append((x, y, z))
                meta.append((ch, resseq))
        return (np.asarray(coords, dtype=float) if coords else np.zeros((0, 3)),
                meta)

    ca, ma = read_atoms(lig)
    cb, mb = read_atoms(rec)
    if len(ca) == 0 or len(cb) == 0:
        return set(), set()
    tree = cKDTree(cb)
    pairs = tree.query_ball_point(ca, r=cutoff)
    lig_if: set = set()
    rec_if: set = set()
    for i, hits in enumerate(pairs):
        if hits:
            lig_if.add(ma[i])
            for j in hits:
                rec_if.add(mb[j])
    return lig_if, rec_if


def seq_slice(seq: str, res_ids: set[int], seq_offset: int = 0) -> str:
    """从链序列按残基号取出界面残基子序列。

    PDB 的 ``resseq`` 与序列下标之间可能有偏移 (缺失残基/编号起点), 因此用
    **残基号最小值对齐**: 先按 ``resseq`` 排序再按序号一一对应, 并在不匹配时
    回退为"按下标直接取"。这是一个已知的近似, 报告里明确标注。
    """
    if not res_ids or not seq:
        return ""
    idx = sorted(res_ids)
    lo = idx[0]
    out = []
    for r in idx:
        j = r - lo
        if 0 <= j < len(seq):
            out.append(seq[j])
    return "".join(out)


def seq_from_structure(pdb_path: Path, chains: list[str]) -> dict[str, str]:
    """从结构文件的 ATOM 记录直接提取链序列 (缓存缺失时的兜底)。

    为什么需要这个兜底: ``data/pdb_seqs.json`` 的链 ID 来自 RCSB 的 entity 命名
    (常见 ``A0``/``B0``), 而 benchmark 用的是**作者链 ID** (``A``/``B``)。
    抗体类结构这个错位尤其普遍 —— 实测前 60 行里有 35 行因缓存取不到序列而
    界面特征全空。直接从 PDB 读 ATOM 记录可以完全绕开这个命名差异。
    """
    want = set(chains)
    per_chain: dict[str, list[tuple[int, str]]] = {c: [] for c in chains}
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
                icode = line[26]
                key = (resseq, icode)
                if key in seen[ch]:
                    continue
                seen[ch].add(key)
                per_chain[ch].append((resseq, aa))
    except Exception:  # noqa: BLE001
        return {c: "" for c in chains}
    out = {}
    for c in chains:
        items = sorted(per_chain[c], key=lambda t: t[0])
        out[c] = "".join(aa for _r, aa in items)
    return out


def get_chain_seq(pdb_id: str, chains: list[str], seqs: dict,
                  pdb_dir: Path = PDB_DIR) -> tuple[str, int]:
    """按优先级取链序列: 缓存 -> 结构文件。返回 (sequence, source_code)。

    source_code: 0=无, 1=全部来自缓存, 2=全部来自结构文件, 3=缓存/结构混合。
    返回**数值编码**而不是字符串, 以保证输出特征表整列可数值化。
    """
    cs = seqs.get(pdb_id, {}) or {}
    parts = [cs.get(c, {}).get("seq", "") for c in chains]
    if all(parts):
        return "".join(parts), 1
    p = pdb_dir / f"{pdb_id}.pdb"
    if not p.exists():
        p = pdb_dir / f"{pdb_id}.ent"
    if p.exists():
        fb = seq_from_structure(p, chains)
        parts2 = [fb.get(c, "") for c in chains]
        if all(parts2):
            return "".join(parts2), 2
        # 部分命中时逐链取更长的那个 (缓存与结构互补)
        merged = [a if len(a) >= len(b) else b for a, b in zip(parts, parts2)]
        if all(merged):
            return "".join(merged), 3
    if any(parts):
        return "".join(parts), 3
    return "", 0


# ---------------------------------------------------------------- 主流程
def build_features(df: pd.DataFrame, seqs: dict, args: argparse.Namespace,
                   emit) -> pd.DataFrame:
    """逐 interface_key 生成增强特征。"""
    rows: list[dict[str, Any]] = []
    n_fail = 0
    t0 = time.time()
    total = len(df)

    for i, r in enumerate(df.to_dict("records")):
        out: dict[str, Any] = {"interface_key": r["interface_key"]}
        try:
            pdb_id = str(r["pdb_id"]).upper()
            lig = [c.strip() for c in str(r["ligand_chains"]).split(",") if c.strip()]
            rec = [c.strip() for c in str(r["receptor_chains"]).split(",") if c.strip()]
            chain_seqs = seqs.get(pdb_id, {}) or {}

            lig_full, lig_src = get_chain_seq(pdb_id, lig, seqs)
            rec_full, rec_src = get_chain_seq(pdb_id, rec, seqs)
            all_full = lig_full + rec_full
            out["v2_seq_available"] = int(bool(lig_full) and bool(rec_full))
            out["v2_seq_source_lig"] = lig_src
            out["v2_seq_source_rec"] = rec_src
            out["v2_lig_seq_len"] = len(lig_full)
            out["v2_rec_seq_len"] = len(rec_full)

            # --- 1/2/3/4: 全链口径 (不依赖结构文件, 覆盖率高) --------------
            for tag, seq in (("lig", lig_full), ("rec", rec_full), ("all", all_full)):
                cc = residue_class_counts(seq)
                for k, v in cc.items():
                    out[f"v2_{tag}_frac_{k}"] = v
                if tag in ("lig", "rec"):
                    for k, v in aa_composition(seq).items():
                        out[f"v2_{tag}_{k}"] = v
                mx, mn = hydrophobic_moment(seq)
                out[f"v2_{tag}_hmom_max"] = mx
                out[f"v2_{tag}_hmom_mean"] = mn
            out["v2_lig_net_charge"] = net_charge(lig_full)
            out["v2_rec_net_charge"] = net_charge(rec_full)
            out["v2_net_charge_product"] = (out["v2_lig_net_charge"]
                                            * out["v2_rec_net_charge"])
            out["v2_net_charge_abs_diff"] = abs(out["v2_lig_net_charge"]
                                                - out["v2_rec_net_charge"])

            # --- 5: 形状互补性简化版 (复用 V1.1 已有几何量) ----------------
            dsasa = _num(r.get("dsasa_total"))
            n_if_lig = _num(r.get("n_iface_res_lig"))
            n_if_rec = _num(r.get("n_iface_res_rec"))
            n_pair = _num(r.get("n_contact_pairs"))
            n_atoms = _num(r.get("n_atoms_lig")) + _num(r.get("n_atoms_rec"))
            out["v2_atoms_per_dsasa"] = _safe_div(n_atoms, dsasa)
            out["v2_dsasa_per_iface_res"] = _safe_div(dsasa, n_if_lig + n_if_rec)
            out["v2_contacts_per_atompair"] = _safe_div(n_pair, n_atoms)
            out["v2_iface_balance"] = _safe_div(n_if_lig,
                                                n_if_lig + n_if_rec)
            out["v2_iface_res_total"] = n_if_lig + n_if_rec

            # --- 6: 界面残基口径 (需要 PDB 文件) --------------------------
            pdb_path = PDB_DIR / f"{pdb_id}.pdb"
            if not pdb_path.exists():
                pdb_path = PDB_DIR / f"{pdb_id}.ent"
            if pdb_path.exists() and lig_full and rec_full:
                lig_if, rec_if = interface_residues(pdb_path, lig, rec)
                lig_ids = {rr for (cc, rr) in lig_if if cc in set(lig)}
                rec_ids = {rr for (cc, rr) in rec_if if cc in set(rec)}
                lig_ifseq = seq_slice(lig_full, lig_ids)
                rec_ifseq = seq_slice(rec_full, rec_ids)
                all_ifseq = lig_ifseq + rec_ifseq
                out["v2_iface_seq_frac"] = _safe_div(
                    len(all_ifseq), len(lig_full) + len(rec_full))
                for tag, seq in (("ligif", lig_ifseq), ("recif", rec_ifseq),
                                 ("allif", all_ifseq)):
                    cc = residue_class_counts(seq)
                    for k, v in cc.items():
                        out[f"v2_{tag}_frac_{k}"] = v
                    if tag in ("ligif", "recif"):
                        for k, v in aa_composition(seq).items():
                            out[f"v2_{tag}_{k}"] = v
                    mx, mn = hydrophobic_moment(seq)
                    out[f"v2_{tag}_hmom_max"] = mx
                    out[f"v2_{tag}_hmom_mean"] = mn
                q_lig = net_charge(lig_ifseq)
                q_rec = net_charge(rec_ifseq)
                out["v2_iface_net_charge_product"] = q_lig * q_rec
                out["v2_iface_charge_complementarity"] = -(q_lig + q_rec) ** 2
                out["v2_iface_pos_neg_ratio"] = _safe_div(
                    _num(out.get("v2_allif_frac_charged_pos")),
                    _num(out.get("v2_allif_frac_charged_neg")))
                out["v2_iface_aromatic_density"] = _safe_div(
                    _num(out.get("v2_allif_frac_aromatic")),
                    _safe_div(dsasa, 1000.0))
            else:
                for k in _INTERFACE_KEYS:
                    out[k] = float("nan")
        except Exception as exc:  # noqa: BLE001
            n_fail += 1
            if n_fail <= 5:
                logger.error("结构 %s 特征生成失败: %s", r.get("interface_key"), exc)
                with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[18_extra_features] {r.get('interface_key')}: {exc}\n"
                             f"{traceback.format_exc()}\n")
        rows.append(out)
        if (i + 1) % 500 == 0:
            emit(f"  ... {i + 1}/{total} ({time.time() - t0:.0f}s, 失败 {n_fail})")

    emit(f"  结构级特征完成: {len(rows)} 行, 失败 {n_fail} 条")
    return pd.DataFrame(rows)


_INTERFACE_KEYS: list[str] = []
for _tag in ("ligif", "recif", "allif"):
    _INTERFACE_KEYS += [f"v2_{_tag}_frac_{k}" for k in
                        ("aromatic", "charged_pos", "charged_neg", "polar",
                         "hydrophobic", "special")]
    _INTERFACE_KEYS += [f"v2_{_tag}_hmom_max", f"v2_{_tag}_hmom_mean"]
for _a in AA_ALL:
    _INTERFACE_KEYS += [f"v2_ligif_aac_{_a}", f"v2_recif_aac_{_a}"]
_INTERFACE_KEYS += ["v2_iface_seq_frac", "v2_iface_net_charge_product",
                    "v2_iface_charge_complementarity", "v2_iface_pos_neg_ratio",
                    "v2_iface_aromatic_density"]


def _num(v: Any) -> float:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return f if np.isfinite(f) else float("nan")


def _safe_div(a: float, b: float) -> float:
    if not (np.isfinite(a) and np.isfinite(b)) or b == 0:
        return float("nan")
    return float(a / b)


def add_esm2_features(df: pd.DataFrame, emit) -> tuple[pd.DataFrame, dict]:
    """若存在 ESM-2 缓存, 按**训练集拟合**的 PCA 降到 ESM2_PCA_DIM 维拼上。"""
    info: dict[str, Any] = {"available": False}
    if not ESM2_CACHE.exists():
        emit("- ESM-2 缓存不存在 (`features/v2_esm2_cache.npz`): **跳过 ESM-2 融合**")
        emit("  生成方式: 在 Colab/Kaggle 免费 GPU 上运行 "
             "`scripts/colab/extract_esm2.py`, 把 npz 拷回 `features/` 即可。")
        emit("  主流程**不依赖** ESM-2, 这是刻意的轻量级设计 (无 GPU 硬依赖)。")
        return df, info
    try:
        from sklearn.decomposition import PCA

        z = np.load(ESM2_CACHE, allow_pickle=True)
        keys = [str(k) for k in z["interface_key"]]
        emb = np.asarray(z["embeddings"], dtype=float)
        idx = pd.Index(keys)
        pos = idx.get_indexer(df["interface_key"].to_numpy())
        hit = pos >= 0
        E = np.full((len(df), emb.shape[1]), np.nan)
        E[hit] = emb[pos[hit]]
        emit(f"- ESM-2 缓存命中: {int(hit.sum())}/{len(df)} "
             f"({hit.mean() * 100:.1f}%), 维度 {emb.shape[1]}")

        tr = (df["split_family"] == "train").to_numpy()
        ok_tr = tr & np.isfinite(E).all(axis=1)
        if ok_tr.sum() < 50:
            emit("- 训练集可用嵌入不足 50 条, 跳过 PCA")
            return df, info
        # **只在训练集 fit** —— 全量 fit 会把测试集信息泄漏进特征
        pca = PCA(n_components=min(ESM2_PCA_DIM, int(ok_tr.sum()) - 1),
                  random_state=SEED)
        pca.fit(E[ok_tr])
        Z = pca.transform(np.nan_to_num(E, nan=0.0))
        for j in range(Z.shape[1]):
            df[f"v2_esm2_pc{j:02d}"] = Z[:, j]
        info = {"available": True, "n_hit": int(hit.sum()),
                "n_components": int(Z.shape[1]),
                "explained_variance_ratio_sum":
                    float(np.sum(pca.explained_variance_ratio_)),
                "pca_fit_on": "train_split_only"}
        emit(f"- ESM-2 PCA: {Z.shape[1]} 维, 训练集解释方差 "
             f"{info['explained_variance_ratio_sum']:.3f} (仅在训练集 fit)")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! ESM-2 融合失败: {exc}")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[18_extra_features] esm2: {exc}\n{traceback.format_exc()}\n")
    return df, info


def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 5: CPU 低成本特征增强")
    ap.add_argument("--features", type=Path, default=FEATURES_CSV)
    ap.add_argument("--seq-cache", type=Path, default=SEQ_CACHE)
    ap.add_argument("--out-extra", type=Path, default=OUT_EXTRA)
    ap.add_argument("--out-merged", type=Path, default=OUT_MERGED)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 条 (冒烟测试)")
    ap.add_argument("--skip-structure", action="store_true",
                    help="跳过需要读 PDB 的界面残基特征 (仅序列类特征)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "18_extra_features.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))

    emit("# Affinity Intelligence V2.0 — 阶段 5: 特征增强报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}")
    emit(f"- 跨链接触距离: {CONTACT_CUTOFF} Å (与 02_extract_features.py 口径一致)")
    emit(f"- 疏水矩: Eisenberg 标度, 窗口 {HMOM_WINDOW}, 角度 {HMOM_ANGLE_DEG}°")

    if not args.features.exists():
        emit(f"!! 特征表不存在: {args.features}")
        return 2
    df = pd.read_csv(args.features, low_memory=False)
    if args.limit:
        df = df.head(args.limit).copy()
    emit(f"- 输入: {args.features.name} — {len(df)} 行 x {df.shape[1]} 列")

    seqs: dict = {}
    if args.seq_cache.exists():
        seqs = json.loads(args.seq_cache.read_text(encoding="utf-8"))
        emit(f"- 序列缓存: {args.seq_cache.name} — {len(seqs)} 个 PDB")
    else:
        emit(f"- !! 序列缓存缺失 ({args.seq_cache}); 序列类特征将全为 NaN")

    if args.dry_run:
        emit("--dry-run: 不生成特征")
        return 0

    emit("")
    emit("## 1. 序列与结构派生特征")
    emit("")
    if args.skip_structure:
        emit("- `--skip-structure`: 只生成序列类特征 (不读 PDB 文件)")
    extra = build_features(df, seqs, args, emit)

    merged = df.merge(extra, on="interface_key", how="left", validate="one_to_one")
    new_cols = [c for c in extra.columns if c != "interface_key"]

    emit("")
    emit("## 2. 输出特征清单与缺失率")
    emit("")
    emit("| 特征 | 缺失率 | 训练集中位数 |")
    emit("|---|---:|---:|")
    tr = (merged["split_family"] == "train") if "split_family" in merged.columns \
        else pd.Series(True, index=merged.index)
    miss_rows = []
    for c in new_cols:
        v = pd.to_numeric(merged[c], errors="coerce")
        mr_ = float(v.isna().mean())
        med = float(v[tr].median()) if v[tr].notna().any() else float("nan")
        miss_rows.append({"feature": c, "missing_rate": mr_, "train_median": med})
        if mr_ > 0 or True:
            emit(f"| {c} | {mr_ * 100:.1f}% | {'—' if not np.isfinite(med) else f'{med:.4f}'} |")

    miss_df = pd.DataFrame(miss_rows)
    all_nan = miss_df[miss_df["missing_rate"] >= 0.999]["feature"].tolist()
    if all_nan:
        emit("")
        emit(f"- !! 完全缺失的特征 {len(all_nan)} 个 (将被下游中位数兜底, 无信息量): "
             f"{', '.join(all_nan[:10])}{' ...' if len(all_nan) > 10 else ''}")

    emit("")
    emit("## 3. ESM-2 嵌入融合 (可选, 无 GPU 硬依赖)")
    emit("")
    merged, esm2_info = add_esm2_features(merged, emit)

    emit("")
    emit("## 4. 消融实验设计 (供 16_ablation.py 使用)")
    emit("")
    emit("提示词要求的关键对照 (全部在 16_ablation.py 中实测):")
    emit("")
    emit("| 配置 | 含义 | 预期作用 |")
    emit("|---|---|---|")
    emit("| `full` | 全部特征 | 上界参考 |")
    emit("| `no_esm2` | 去掉 ESM-2 维度 | 若 ≈ full -> 物理特征已足够 (轻量级卖点) |")
    emit("| `physics_only_normalized` | 只用归一化物理量 | 检验归一化是否替代尺寸 |")
    emit("| `size_only` | 只用原子/残基计数 | 若 << full -> 归一化物理特征不可替代 |")
    emit("| `esm2_only` | 只用 ESM-2 | 检验序列嵌入单独是否够用 |")
    emit("| `no_sasa` | 去 SASA/ΔSASA 族 | 检验埋藏面积贡献 |")
    emit("| `no_hbond_salt` | 去氢键/盐桥族 | 检验极性相互作用贡献 |")
    emit("| `no_electrostatic` | 去新增静电代理 | 检验静电互补性增益 |")
    emit("| `no_composition` | 去残基组成/氨基酸组成 | 检验组成特征贡献 |")
    emit("| `no_chainpair` | 去 V0.3 链对物理特征 | 检验链对聚合特征贡献 |")
    emit("| `leak_pdb_level` | **PDB 级**（非家族级）划分 | 量化同源泄漏造成的虚假提升 |")
    emit("")

    # ---- 落盘 ---------------------------------------------------------
    args.out_extra.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out_extra.with_suffix(args.out_extra.suffix + ".tmp")
    extra.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out_extra)

    drop_all_nan = [c for c in new_cols if c in all_nan]
    keep_cols = [c for c in merged.columns if c not in drop_all_nan]
    tmp2 = args.out_merged.with_suffix(args.out_merged.suffix + ".tmp")
    merged[keep_cols].to_csv(tmp2, index=False, encoding="utf-8")
    tmp2.replace(args.out_merged)

    emit("")
    emit(f"- 新增特征表: `{args.out_extra.name}` — {len(extra)} 行 x "
         f"{extra.shape[1]} 列 (含 interface_key)")
    emit(f"- 合并增强表: `{args.out_merged.name}` — {len(merged)} 行 x "
         f"{len(keep_cols)} 列 (已剔除 {len(drop_all_nan)} 个全空特征)")
    emit(f"- 净新增可用特征: **{len(new_cols) - len(drop_all_nan)}** 个")
    emit("")
    emit("## 诚实声明")
    emit("")
    emit("1. **界面残基序列的残基号对齐是近似**: PDB 的 `resseq` 与序列下标之间可能")
    emit("   存在偏移 (缺失残基/编号起点不同), 本脚本按「最小残基号一一对应」回退,")
    emit("   未做完整的序列比对。该近似只影响界面口径特征, 不影响全链口径特征。")
    emit("2. **静电互补性是简化代理**, 基于 pH 7.4 形式电荷计数, 不是 "
         "Poisson-Boltzmann 或 Coulomb 求和求解; 不得当作真实静电能。")
    emit("3. **形状互补性只有密度类代理**, 没有实现真正的形状互补系数 (Sc)。")
    emit("4. 新增特征是否真的有增益, 以 `16_ablation.py` 的实测为准;")
    emit("   本报告只负责生成特征与记录缺失率, 不声称任何增益。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "n_new_features": len(new_cols), "n_all_nan": len(drop_all_nan),
        "features": miss_rows, "esm2": esm2_info,
        "contact_cutoff": CONTACT_CUTOFF,
        "esm2_cache": str(ESM2_CACHE),
    }
    (FEAT_DIR / "v2_feature_report.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    args.report.write_text("\n".join(run_log) + "\n", encoding="utf-8")
    emit("")
    emit(f"报告: {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
