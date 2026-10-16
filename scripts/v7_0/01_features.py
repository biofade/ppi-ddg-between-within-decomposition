# -*- coding: utf-8 -*-
"""scripts/v7_0/01_features.py -- V7.0 PHASE 0 + PHASE 1: 冻结 + 特征块构造.

PHASE 0：冻结输入哈希 + 定义 baseline B0（复用历史 M0，不重写）
PHASE 1：构造三个特征块（**必须分 block，禁止一开始混合**）
  * **BLOCK S** —— Set context（V6.3 的 deployable mutation-set descriptors）
  * **BLOCK L** —— Mutation-local（**从本地 PDB 坐标真实构造**）
  * **BLOCK G** —— Global target context（deployable 结构描述符）

⚠️ 可用性审计结论（必须如实记录）
--------------------------------
1. **`features/v2_esm2_cache.npz` 只有 pooled 链级嵌入**（640 维/链），
   **没有 residue-level 嵌入** → V7.0 §7 BLOCK L 的
   `Δh_m = h_mut_m − h_wt_m` **【不可构造】**
   ⇒ 按 §7 允许的回退方案：使用 **WT 残基性质 + mutation identity 向量 +
   局部 partner 环境**，并在报告中显式说明（**不假装有 residue embedding**）。
2. **本地有 3031 个 PDB，覆盖 230/230 targets（100%）** → 真实局部结构环境
   （半径内接触数、partner 链接触、到界面距离）**可以构造**。
3. `data/v5_0/mutation_sites.csv` 提供 mutation→PDB 残基映射
   （`mapping_status == EXACT`，canonical 覆盖 89.6%）。

**禁止**把 pooled `Z_mut − Z_wt` 当作 Local block（§7 明令）。
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V63 = C.ROOT / "experiments" / "v6_3"
PH0 = C.ROOT / "experiments" / "v7_0" / "phase0_baseline"
PH1 = C.ROOT / "experiments" / "v7_0" / "phase1_features"
for d in (PH0, PH1):
    d.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
L: list[str] = []

# ---------------- 氨基酸理化性质（用于 BLOCK L 的 mutation identity） ----------------
AA1 = list("ACDEFGHIKLMNPQRSTVWY")
# Kyte-Doolittle 疏水性
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# 体积（Å³）
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
       "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
       "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
       "Y": 193.6, "V": 140.0}
# 电荷（pH 7.4 近似）
CHG = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}
# 极性
POL = {"R": 1, "N": 1, "D": 1, "Q": 1, "E": 1, "H": 1, "K": 1, "S": 1, "T": 1,
       "Y": 1, "C": 0.5}
# 芳香性
ARO = {"F": 1, "W": 1, "Y": 1, "H": 0.5}
# 柔性（B-factor 代理：Gly/Ser 高，Pro 低）
FLEX = {"G": 1.0, "S": 0.9, "P": 0.2, "A": 0.6, "T": 0.7, "N": 0.7, "D": 0.7,
        "C": 0.5, "V": 0.5, "I": 0.4, "L": 0.4, "M": 0.4, "F": 0.4, "Y": 0.4,
        "W": 0.4, "H": 0.5, "K": 0.7, "R": 0.7, "Q": 0.7, "E": 0.7}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


# ============================================================ PDB 解析
def parse_pdb(path: Path):
    """返回 (residues, ca_map)。

    residues: list of dict(chain, resseq, icode, resname, ca=(x,y,z), atoms=[coords])
    只取 ATOM 记录（标准氨基酸）。
    """
    res = {}
    order = []
    with open(path, "r", errors="replace") as f:
        for ln in f:
            if not ln.startswith("ATOM"):
                continue
            resname = ln[17:20].strip().upper()
            if resname not in ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                               "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                               "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"):
                continue
            ch = ln[21]
            try:
                seq = int(ln[22:26])
            except ValueError:
                continue
            ic = ln[26].strip()
            key = (ch, seq, ic)
            if key not in res:
                res[key] = {"chain": ch, "resseq": seq, "icode": ic,
                            "resname": "MET" if resname == "MSE" else resname,
                            "ca": None, "atoms": []}
                order.append(key)
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError:
                continue
            res[key]["atoms"].append((x, y, z))
            if ln[12:16].strip() == "CA":
                res[key]["ca"] = (x, y, z)
    out = [res[k] for k in order if res[k]["ca"] is not None]
    ca_map = {(r["chain"], r["resseq"], r["icode"]): r["ca"] for r in out}
    return out, ca_map


def build_local_features(pdb_id: str, sites: list[dict], cache: dict):
    """对给定 target 的全部 mutation site 计算 BLOCK L 特征。

    半径：primary 8 Å（同时给 6/10 Å 以便敏感性分析）
    """
    key = pdb_id.upper()
    if key not in cache:
        p = PDBDIR / f"{key}.pdb"
        if not p.exists():
            cache[key] = None
        else:
            cache[key] = parse_pdb(p)
    got = cache[key]
    if got is None:
        return None
    res, ca_map = got
    if len(res) < 10:
        return None
    coords = np.array([r["ca"] for r in res], float)
    chains = np.array([r["chain"] for r in res])
    # ---- 界面：跨链 CA 距离 < 10 Å 的残基集合 ----
    uchain = sorted(set(chains.tolist()))
    iface = np.zeros(len(res), dtype=bool)
    if len(uchain) >= 2:
        from scipy.spatial import cKDTree
        trees = {c: cKDTree(coords[chains == c]) for c in uchain}
        for i, c in enumerate(uchain):
            for c2 in uchain[i + 1:]:
                d, _ = trees[c].query(coords[chains == c2], k=1,
                                      distance_upper_bound=10.0)
                m = np.isfinite(d)
                if m.any():
                    iface[chains == c2] |= m
                d2, _ = trees[c2].query(coords[chains == c], k=1,
                                        distance_upper_bound=10.0)
                m2 = np.isfinite(d2)
                if m2.any():
                    iface[chains == c] |= m2
    # ---- 每残基的局部接触（同链 + 异链），三个半径 ----
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    res_index = {(r["chain"], r["resseq"], r["icode"]): i
                 for i, r in enumerate(res)}
    feats = []
    for st in sites:
        k = (st["pdb_chain"], st["pdb_resseq"], st["icode"])
        idx = res_index.get(k)
        if idx is None:
            # 允许 icode 不匹配时按 (chain,resseq) 兜底
            cand = [i for i, r in enumerate(res)
                    if r["chain"] == st["pdb_chain"]
                    and r["resseq"] == st["pdb_resseq"]]
            if not cand:
                feats.append(None)
                continue
            idx = cand[0]
        d_all = np.linalg.norm(coords - coords[idx], axis=1)
        f = {"mut_dist_to_iface": float(np.min(d_all[iface])) if iface.any()
             else float(d_all.max()),
             "mut_is_iface": int(iface[idx])}
        for rad in (6.0, 8.0, 10.0):
            sel = (d_all <= rad) & (np.arange(len(res)) != idx)
            f[f"mut_n_contacts_r{int(rad)}"] = int(sel.sum())
            partner = sel & (chains != chains[idx])
            f[f"mut_n_partner_contacts_r{int(rad)}"] = int(partner.sum())
            f[f"mut_mean_dist_r{int(rad)}"] = (float(d_all[sel].mean())
                                               if sel.any() else float(rad))
            # partner 的理化性质分布
            if partner.any():
                pn = [res[i]["resname"] for i in np.where(partner)[0]]
                aa = [THREE2ONE.get(x, "X") for x in pn]
                f[f"partner_mean_kd_r{int(rad)}"] = float(
                    np.mean([KD.get(a, 0.0) for a in aa]))
                f[f"partner_mean_vol_r{int(rad)}"] = float(
                    np.mean([VOL.get(a, 0.0) for a in aa]))
                f[f"partner_mean_chg_r{int(rad)}"] = float(
                    np.mean([CHG.get(a, 0.0) for a in aa]))
                f[f"partner_frac_charged_r{int(rad)}"] = float(
                    np.mean([1.0 if a in "DEKR" else 0.0 for a in aa]))
                f[f"partner_frac_hydrophobic_r{int(rad)}"] = float(
                    np.mean([1.0 if a in "AVLIMFW" else 0.0 for a in aa]))
                f[f"partner_frac_aromatic_r{int(rad)}"] = float(
                    np.mean([1.0 if a in "FWY" else 0.0 for a in aa]))
            else:
                for kk in ("partner_mean_kd", "partner_mean_vol",
                           "partner_mean_chg", "partner_frac_charged",
                           "partner_frac_hydrophobic", "partner_frac_aromatic"):
                    f[f"{kk}_r{int(rad)}"] = 0.0
            # 同链接触
            same = sel & (chains == chains[idx])
            f[f"mut_n_samechain_contacts_r{int(rad)}"] = int(same.sum())
        feats.append(f)
    return feats


THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
             "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
             "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
             "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def mut_identity(wt: str, mu: str) -> dict:
    """BLOCK L 的 mutation identity（§7.1）。"""
    d = {}
    d["id_d_charge"] = CHG.get(mu, 0.0) - CHG.get(wt, 0.0)
    d["id_d_hydro"] = KD.get(mu, 0.0) - KD.get(wt, 0.0)
    d["id_d_volume"] = VOL.get(mu, 0.0) - VOL.get(wt, 0.0)
    d["id_d_polarity"] = POL.get(mu, 0.0) - POL.get(wt, 0.0)
    d["id_d_aromatic"] = ARO.get(mu, 0.0) - ARO.get(wt, 0.0)
    d["id_d_flex"] = FLEX.get(mu, 0.0) - FLEX.get(wt, 0.0)
    d["id_wt_hydro"] = KD.get(wt, 0.0)
    d["id_wt_volume"] = VOL.get(wt, 0.0)
    d["id_mut_hydro"] = KD.get(mu, 0.0)
    d["id_mut_volume"] = VOL.get(mu, 0.0)
    d["id_is_ala"] = int(mu == "A")
    d["id_is_gly"] = int(mu == "G")
    d["id_is_pro"] = int(mu == "P")
    # one-hot WT / MUT（20 维各）
    for a in AA1:
        d[f"id_wt_{a}"] = int(wt == a)
        d[f"id_mut_{a}"] = int(mu == a)
    return d


def main() -> int:
    t0 = time.time()
    emit("# V7.0 PHASE 0 + 1 — 冻结 + 特征块构造")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ================= PHASE 0 =================
    emit("## PHASE 0 — 输入冻结")
    emit("")
    FILES = [
        ("matched_analysis_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
        ("split4_predictions", "experiments/v6_2/phase2_split4/split4_predictions.npz"),
        ("split4_results", "experiments/v6_2/phase2_split4/split4_results.csv"),
        ("d_only_predictions", "experiments/v6_2_1/phase3_d_only/d_only_predictions.npz"),
        ("v63_target_table", "experiments/v6_3/phase1_target_table/target_mechanism_table.csv"),
        ("v63_feature_dict", "experiments/v6_3/phase1_target_table/feature_dictionary.csv"),
        ("ddg_split", "data/v4.0/ddg_split.json"),
        ("mutation_sites", "data/v5_0/mutation_sites.csv"),
        ("v22_features", "features/v2.2_proteinea_features.csv"),
        ("esm2_cache", "features/v2_esm2_cache.npz"),
        ("v5_common", "scripts/v5/v5_common.py"),
        ("standardized_preds", "experiments/v6_sota_audit/standardized_preds.csv"),
    ]
    rows = []
    for name, rel in FILES:
        p = C.ROOT / rel
        if p.exists():
            rows.append({"name": name, "path": rel, "size": p.stat().st_size,
                         "sha256_prefix40": sha(p)[:40]})
            emit(f"- `{name}`: {p.stat().st_size:,} B  `{sha(p)[:16]}…`")
        else:
            rows.append({"name": name, "path": rel, "size": 0,
                         "sha256_prefix40": "MISSING"})
            emit(f"- `{name}`: **MISSING**")
    pd.DataFrame(rows).to_csv(PH0 / "input_manifest.csv", index=False,
                              encoding="utf-8")
    emit("")

    # B0 定义
    emit("### Baseline B0（= 历史 M0，**复用不重写**）")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit("| 定义 | pooled baseline：`[phys_diff(60), Z_diff(128)]` + fwd/rev 增广 |")
    emit("| 实现 | `scripts/v6_2/02_split4_ddg.py`（V6.2 §2.0 已审计） |")
    emit("| 预测来源 | `experiments/v6_2/phase2_split4/split4_predictions.npz` 的 `L4_family_M0_s{42,43,44}` |")
    emit("| LightGBM | `LGB_PARAMS` + `n_estimators=300` + `random_state=seed` |")
    emit("| 特征维度 | 188（60 物理 + 128 pooled ESM-2 PCA） |")
    emit("")
    (PH0 / "baseline_definition.md").write_text("\n".join(L) + "\n",
                                                encoding="utf-8")
    emit("")

    # ================= PHASE 1 =================
    emit("## PHASE 1 — 特征块")
    emit("")
    T63 = pd.read_csv(V63 / "phase1_target_table" / "target_mechanism_table.csv",
                      low_memory=False)
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    sites = site[site.pair_row.isin(can.row.values)
                 & (site.mapping_status == "EXACT")].copy()
    emit(f"### BLOCK L 的 mutation-site 映射")
    emit(f"- `mutation_sites.csv`（`mapping_status == EXACT`）覆盖 canonical: "
         f"**{sites.pair_row.nunique()} / {len(can)}** "
         f"({sites.pair_row.nunique() / len(can) * 100:.1f}%)")
    emit("")

    cache: dict = {}
    recs = []
    n_nopdb = 0
    pdb_ids = sorted(can.pdb_id.astype(str).unique())
    emit(f"### 构造中（{len(pdb_ids)} targets，逐 target 解析 PDB）")
    for i, pdb in enumerate(pdb_ids):
        sub = sites[sites.pdb_id.astype(str).str.upper() == pdb.upper()]
        if not len(sub):
            n_nopdb += 1
            continue
        gp = sub.groupby("pair_row")
        per_site = []
        for pr, g in gp:
            lst = [{"pdb_chain": r.pdb_chain if hasattr(r, "pdb_chain")
                    else r.ann_chain,
                    "pdb_resseq": int(r.pdb_resseq),
                    "icode": getattr(r, "icode", "") or "",
                    "wt_aa": str(r.wt_aa), "mut_aa": str(r.mut_aa),
                    "side": str(r.side)} for r in g.itertuples()]
            per_site.append((pr, lst))
        allsites = [s for _, lst in per_site for s in lst]
        lf = build_local_features(pdb, allsites, cache)
        if lf is None:
            n_nopdb += 1
            continue
        # 回填到 pair_row
        cur = 0
        for pr, lst in per_site:
            here = lf[cur:cur + len(lst)]
            cur += len(lst)
            ok = [h for h in here if h is not None]
            if not ok:
                continue
            agg = {}
            keys = set(ok[0].keys())
            for k in keys:
                vals = [h[k] for h in ok if k in h]
                agg[k] = float(np.mean(vals))
                agg[f"{k}__max"] = float(np.max(vals))
            # mutation identity（每 variant 的**平均** wt→mut 身份向量）
            idv = [mut_identity(s["wt_aa"], s["mut_aa"]) for s in lst
                   if s["wt_aa"] in AA1 and s["mut_aa"] in AA1]
            if idv:
                for k in idv[0]:
                    agg[k] = float(np.mean([d[k] for d in idv]))
            agg["pair_row"] = int(pr)
            agg["n_sites_mapped"] = len(lst)
            recs.append(agg)
        if (i + 1) % 40 == 0:
            emit(f"  - 已处理 {i + 1}/{len(pdb_ids)} targets，"
                 f"累计 rows {len(recs)}")
    emit("")
    LF = pd.DataFrame(recs)
    emit(f"### BLOCK L 覆盖")
    emit(f"- 构造出 BLOCK L 的 canonical rows: **{len(LF)} / {len(can)}** "
         f"({len(LF) / len(can) * 100:.1f}%)")
    emit(f"- 无映射/无 PDB 的 target 数: **{n_nopdb}**")
    emit(f"- BLOCK L 特征列数: **{LF.shape[1] - 2}**")
    emit("")
    LF.to_csv(PH1 / "block_L_local.csv", index=False, encoding="utf-8")

    # ---------- BLOCK S ----------
    SET_COLS = ["log_n_mutations", "n_mutations", "n_unique_mutation_sites",
                "mutations_per_site", "max_mutations_per_site", "site_entropy",
                "mutation_pair_entropy", "WT_residue_entropy",
                "Mut_residue_entropy", "alanine_fraction", "glycine_fraction",
                "charged_to_neutral_fraction", "hydrophobic_to_polar_fraction",
                "conservative_fraction", "nonconservative_fraction",
                "single_mut_fraction", "multi_mut_fraction",
                "mean_mutations_per_variant", "max_mutations_per_variant",
                "site_count_gini", "mutation_count_gini", "site_coverage_all",
                "site_coverage_interface", "protein_length",
                "interface_residue_count"]
    haveS = [c for c in SET_COLS if c in T63.columns]
    S = T63[["pdb_id"] + haveS].copy()
    if "multi_mut_fraction" not in S.columns:
        S["multi_mut_fraction"] = 1.0 - T63["single_mut_fraction"]
        haveS.append("multi_mut_fraction")
    emit(f"### BLOCK S 覆盖")
    emit(f"- 特征列: **{len(haveS)}**；target 覆盖 **{len(S)}** / "
         f"{can.pdb_id.nunique()}")
    S.to_csv(PH1 / "block_S_set.csv", index=False, encoding="utf-8")
    emit("")

    # ---------- BLOCK G ----------
    G_COLS = [c for c in T63.columns if c.startswith("struct_")]
    G_COLS = [c for c in G_COLS if pd.api.types.is_numeric_dtype(
        pd.to_numeric(T63[c], errors="coerce")) and T63[c].notna().sum() >= 150]
    G = T63[["pdb_id"] + G_COLS].copy()
    emit(f"### BLOCK G 覆盖")
    emit(f"- 特征列: **{len(G_COLS)}**；target 覆盖 **{len(G)}**")
    G.to_csv(PH1 / "block_G_global.csv", index=False, encoding="utf-8")
    emit("")

    # ---------- 覆盖与排除报告 ----------
    emit("## 覆盖 / 排除报告（§5 要求）")
    emit("")
    lf_rows = set(LF.pair_row.astype(int))
    allrows = set(can.row.astype(int))
    excl = sorted(allrows - lf_rows)
    emit(f"- canonical rows: **{len(allrows)}**")
    emit(f"- BLOCK L 可用 rows: **{len(lf_rows)}**")
    emit(f"- **excluded rows（无 BLOCK L）: {len(excl)}**"
         f"（{len(excl) / len(allrows) * 100:.1f}%）")
    excl_tg = sorted(can[can.row.isin(excl)].pdb_id.astype(str).unique())
    emit(f"- **excluded targets: {len(excl_tg)}** → {excl_tg[:20]}")
    emit("")
    emit("> ⚠️ 按 §5：**禁止拿不同 row universe 的模型直接比较**。")
    emit("> 因此 PHASE 2 的所有模型将在 **BLOCK L 可用的 matched-coverage 子集**")
    emit("> 上训练与评估；同时另报 A0 在全 7523 行上的读数作为参考。")
    emit("")
    json.dump({"canonical_rows": len(allrows),
               "blockL_rows": len(lf_rows),
               "excluded_rows": len(excl),
               "excluded_targets": excl_tg},
              open(PH1 / "coverage.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "feature_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH0 / 'input_manifest.csv'} / baseline_definition.md")
    print(f"写出 {PH1}/ block_L_local.csv, block_S_set.csv, block_G_global.csv, "
          f"feature_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
