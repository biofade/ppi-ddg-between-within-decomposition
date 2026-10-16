# -*- coding: utf-8 -*-
"""scripts/v7_3/02_build_physics_features.py -- V7.3 PHASE 3.

构建 **mutation-conditioned local physical interaction perturbation proxies**
（§11–§24）。**术语纪律**：这些是 `fixed-backbone physics-derived perturbation
proxies`，**不是**精确自由能；禁止称 `kcal/mol energy` / `binding energy`。

特征族（§26，全部在**看 outer 结果之前**冻结）
--------------------------------------------
| 族 | 构成 |
|---|---|
| `F1`  | mutation property deltas（`F_MUT`，§13） |
| `F2`  | **STATIC_LOCAL_CONTROL** = V7.0 BLOCK L（§25，只读复用） |
| `F3`  | partner static environment（`F_ENV_PARTNER`，§17） |
| `F4`  | **mutation × partner 扰动**（核心，§19） |
| `F5`  | mutation × intra 扰动（`F_ENV_INTRA`，§18） |
| `F6`  | partner + intra 扰动 |
| `F7`  | F6 + burial/SASA 条件化扰动（§20） |
| `F8`  | F7 + 多突变位点耦合几何（`F_EPI_GEOM`，§24） |
| `F9`  | **PHYSICS-ONLY**（无历史 Base） |

关键实现约束
------------
* **环境分离**（§14）：`ENV_PARTNER`（跨链）与 `ENV_INTRA`（同链）
  **绝不合并后分析**
* **重原子最小距离**（§15）：不用 CA 代替
* 半径 **8 Å**、权重 `exp(-d/4)`（§16，冻结）
* 多突变聚合：**SUM / MEAN / MAX_ABS**（§23，冻结）

输出：`features/v7_3/physics_features.parquet` + `feature_dict.csv`
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V73 = C.ROOT / "experiments" / "v7_3"
FEAT = C.ROOT / "features" / "v7_3"
CACHE = FEAT / "pdb_environment_cache"
PH3 = V73 / "phase3_features"
PH3.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
R_PRIMARY = 8.0
L: list[str] = []

# ---------------- 20-AA 物理属性表（§12，冻结，记录来源） ----------------
# 来源全部为**已发表的标准标度**，在 `aa_physics_table.csv` 中逐列标注。
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
# Kyte-Doolittle hydropathy (J. Mol. Biol. 157:105, 1982)
KD = dict(zip(AA_LIST, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                        1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
# Side-chain volume (A^3) — Zamyatnin 1972 / standard tabulation
VOL = dict(zip(AA_LIST, [88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7,
                         168.6, 166.7, 162.9, 114.1, 143.8, 112.6, 85.0, 89.0,
                         116.1, 193.6, 113.9, 140.0]))
# Formal charge at pH 7 (side chain)
CHG = {a: 0.0 for a in AA_LIST}
CHG["D"] = -1.0; CHG["E"] = -1.0
CHG["K"] = 1.0; CHG["R"] = 1.0; CHG["H"] = 0.1      # His 部分质子化
# Polarity (Grantham 1974 polarity scale)
GRANTHAM_POL = dict(zip(AA_LIST, [8.1, 5.5, 13.0, 12.3, 5.2, 3.0, 0.0, 10.4,
                                  5.7, 4.9, 5.2, 11.6, 10.5, 10.5, 9.2, 8.6,
                                  8.0, 5.9, 11.6, 5.4]))
# H-bond donor / acceptor capacity（侧链可参与氢键的数目，粗粒度）
DONOR = dict(zip(AA_LIST, [0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
DONOR["S"] = 1; DONOR["T"] = 1; DONOR["Y"] = 1; DONOR["W"] = 1; DONOR["H"] = 1
DONOR["N"] = 1; DONOR["Q"] = 1; DONOR["K"] = 1; DONOR["R"] = 1
ACCEPTOR = dict(zip(AA_LIST, [0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                             1, 0, 1, 1]))
# flags
AROMATIC = {a: float(a in "FWYH") for a in AA_LIST}
CATION = {a: float(a in "KRH") for a in AA_LIST}
ANION = {a: float(a in "DE") for a in AA_LIST}
ALIPHATIC = {a: float(a in "AVLIM") for a in AA_LIST}
GLY = {a: float(a == "G") for a in AA_LIST}
PRO = {a: float(a == "P") for a in AA_LIST}

PROPS = {"charge": CHG, "hydropathy": KD, "volume": VOL, "polarity": GRANTHAM_POL,
         "donor": DONOR, "acceptor": ACCEPTOR, "aromatic": AROMATIC,
         "cation": CATION, "anion": ANION, "aliphatic": ALIPHATIC,
         "gly": GLY, "pro": PRO}
PROP_SRC = {
    "charge": "formal charge at pH 7.4 (His 取 0.1 部分质子化)",
    "hydropathy": "Kyte & Doolittle, J. Mol. Biol. 157:105 (1982)",
    "volume": "side-chain volume (Zamyatnin 1972 / standard tabulation)",
    "polarity": "Grantham, Science 185:862 (1974) polarity scale",
    "donor": "侧链 H-bond donor 能力（粗粒度计数）",
    "acceptor": "侧链 H-bond acceptor 能力（粗粒度计数）",
    "aromatic": "flag: F/W/Y/H",
    "cation": "flag: K/R/H",
    "anion": "flag: D/E",
    "aliphatic": "flag: A/V/L/I/M",
    "gly": "flag: G", "pro": "flag: P",
}
CHEMS = ["positive", "negative", "polar", "hydrophobic", "aromatic", "glypro"]


def chem_of(a: str) -> str:
    if a in "KR": return "positive"
    if a in "DE": return "negative"
    if a in "FWY": return "aromatic"
    if a in "GP": return "glypro"
    if a in "STNQCY": return "polar"
    if a in "AVLIM": return "hydrophobic"
    return "other"


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-targets", type=int, default=10 ** 9)
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.3 PHASE 3 — 物理扰动代理特征")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    emit("> **术语纪律（§11）**：以下全部是 `fixed-backbone physics-derived "
         "perturbation proxies`，**不是**精确自由能。")
    emit("")

    # ---------- AA 属性表 ----------
    rows = []
    for aa in AA_LIST:
        r = {"aa": aa}
        for k, d in PROPS.items():
            r[k] = d[aa]
        r["chem_class"] = chem_of(aa)
        rows.append(r)
    AAT = pd.DataFrame(rows)
    AAT.attrs = {}
    AAT.to_csv(FEAT / "aa_physics_table.csv", index=False, encoding="utf-8")
    emit(f"## 3.1 20-AA 物理属性表（§12）")
    emit("")
    emit("| 属性 | 来源/标度 |")
    emit("|---|---|")
    for k, v in PROP_SRC.items():
        emit(f"| `{k}` | {v} |")
    emit("")
    emit(f"→ `features/v7_3/aa_physics_table.csv`（{len(AAT)} 行 × "
         f"{AAT.shape[1]} 列，**在看 outer 结果之前冻结**）")
    emit("")

    # ---------- 载入 sites / cache ----------
    with open(CACHE / "_sites.pkl", "rb") as f:
        keep = pickle.load(f)
    labs = pd.DataFrame([{k: v for k, v in r.items() if k != "sites"}
                         for r in keep])
    emit(f"- U_PHYS_STRICT：**{len(keep)}** 行 / "
         f"{labs.target_id_canonical.nunique()} targets")
    emit("")

    # ---------- 历史 Base（复用 V7.0 物理+Z 块，与之前所有实验一致）----------
    d = C.build_pair_features(verbose=False)
    blk = d["X_blocks"]
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    all_rows = can.row.to_numpy(int)
    BASE_ALL = np.hstack([blk["phys_diff"], blk["Z_diff"]]).astype(np.float32)
    pos = pd.Series(np.arange(len(all_rows)), index=all_rows)

    # ---------- BLOCK L（STATIC_LOCAL_CONTROL，§25）----------
    BL = pd.read_csv(C.ROOT / "experiments" / "v7_0" / "phase1_features" /
                     "block_L_local.csv", low_memory=False)
    bl_cols = [c for c in BL.columns if c not in ("pair_row", "pdb_id")]
    BL_i = BL.set_index("pair_row")[bl_cols]

    # ---------- 逐 target 构建 ----------
    # 逐 target 构建（每 target 只解析一次 PDB + 建一次 KD 树）
    F = {}          # row_id -> dict
    n_done = 0

    # 提前载入 v7f（parse_heavy 需要 THREE）
    import importlib.util as _ilu
    _sp = _ilu.spec_from_file_location(
        "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
    v7f = _ilu.module_from_spec(_sp)
    sys.modules["v7f"] = v7f
    _sp.loader.exec_module(v7f)
    THREE = v7f.THREE2ONE

    def parse_heavy(p: Path):
        ch_res, ch_atoms = {}, {}
        from collections import defaultdict
        ch_res = defaultdict(list); ch_atoms = defaultdict(list)
        for ln in open(p, "r", errors="replace"):
            if ln[:6] not in ("ATOM  ", "HETATM"):
                continue
            elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if elem == "H":
                continue
            rn = ln[17:20].strip().upper()
            if rn not in THREE:
                continue
            try:
                x, y, zz = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            except ValueError:
                continue
            c = ln[21]
            try:
                rs = int(ln[22:26])
            except ValueError:
                continue
            ic = ln[26].strip()
            key = (rs, ic, rn)
            if not ch_res[c] or ch_res[c][-1] != key:
                ch_res[c].append(key); ch_atoms[c].append([])
            ch_atoms[c][-1].append((x, y, zz))
        return ch_res, ch_atoms

    # 逐 target
    for ti, t in enumerate(sorted(labs.target_id_canonical.astype(str).unique())):
        if ti >= a.max_targets:
            break
        if (time.time() - t0) / 60 > 70:
            emit(f"- ⏱ 软时限，已处理 {ti} targets"); break
        cp = CACHE / f"{t}.npz"
        if not cp.exists():
            continue
        z = np.load(cp, allow_pickle=True)
        chains = z["chain"].astype(str)
        aas = z["aa"].astype(str)
        xyz = z["res_xyz"].astype(np.float64)
        sasa_c, sasa_i, bur = z["sasa_complex"], z["sasa_isolated"], z["burial"]
        env = z["env"].astype(np.float32)
        env_names = [str(x) for x in z["env_names"]]
        ENVI = {nm: i for i, nm in enumerate(env_names)}
        nres = len(aas)
        pp = PDBDIR / f"{t.upper()}.pdb"
        ch_res, ch_atoms = parse_heavy(pp)
        # 重原子数组 + 全局残基下标
        atom_xyz, atom_res = [], []
        chain_res_pos = {}
        gi = 0
        for c in sorted(ch_res):
            for k in range(len(ch_res[c])):
                chain_res_pos[(c, k)] = gi
                for (x, y, zz) in ch_atoms[c][k]:
                    atom_xyz.append((x, y, zz)); atom_res.append(gi)
                gi += 1
        atom_xyz = np.asarray(atom_xyz, np.float64)
        atom_res = np.asarray(atom_res, np.int32)
        tree = cKDTree(atom_xyz) if len(atom_xyz) else None

        def res_env_fields(site_gidx, exclude_gidx, cross_chain):
            """返回该位点周围（跨链或同链）的加权场 dict（§17/§18/§19）。"""
            out = {}
            if tree is None:
                return out
            cand = tree.query_ball_point(xyz[site_gidx], r=R_PRIMARY + 6.0)
            best = {}
            for ai in cand:
                j = int(atom_res[ai])
                if j == site_gidx or j in exclude_gidx:
                    continue
                iscross = chains[j] != chains[site_gidx]
                if iscross != cross_chain:
                    continue
                d = float(np.linalg.norm(atom_xyz[ai] - xyz[site_gidx]))
                if j not in best or d < best[j]:
                    best[j] = d
            if not best:
                return {"present": 0.0}
            js = list(best.keys())
            ds = np.array([best[j] for j in js], np.float64)
            keep = ds <= R_PRIMARY
            if not keep.any():
                return {"present": 0.0}
            js = [js[i] for i in np.where(keep)[0]]
            ds = ds[keep]
            w = np.exp(-ds / 4.0)
            q = np.array([CHG[aas[j]] for j in js])
            h = np.array([KD[aas[j]] for j in js])
            v = np.array([VOL[aas[j]] for j in js])
            po = np.array([GRANTHAM_POL[aas[j]] for j in js])
            dn = np.array([DONOR[aas[j]] for j in js])
            ac = np.array([ACCEPTOR[aas[j]] for j in js])
            ar = np.array([AROMATIC[aas[j]] for j in js])
            ca = np.array([CATION[aas[j]] for j in js])
            an = np.array([ANION[aas[j]] for j in js])
            out = {
                "present": 1.0,
                "n_contact_4": float((ds <= 4.0).sum()),
                "n_contact_6": float((ds <= 6.0).sum()),
                "n_contact_8": float((ds <= 8.0).sum()),
                "min_dist": float(ds.min()),
                "mean_dist": float(ds.mean()),
                "w_charge": float((w * q).sum()),
                "w_charge_over_d": float((w * q / (ds + 1.0)).sum()),
                "w_hydropathy": float((w * h).sum()),
                "w_volume": float((w * v).sum()),
                "w_density": float(w.sum()),
                "w_polarity": float((w * po).sum()),
                "w_donor": float((w * dn).sum()),
                "w_acceptor": float((w * ac).sum()),
                "w_aromatic": float((w * ar).sum()),
                "w_cation": float((w * ca).sum()),
                "w_anion": float((w * an).sum()),
            }
            for cn in CHEMS:
                m = np.array([chem_of(aas[j]) == cn for j in js])
                out[f"c_{cn}"] = float(m.sum())
                out[f"wc_{cn}"] = float(w[m].sum())
            return out

        # 该 target 的行
        trew = [r for r in keep if str(r["target_id_canonical"]) == t]
        for rec in trew:
            rid = int(rec["row_id"])
            site_feats = []
            for s in rec["sites"]:
                gidx = int(s["gidx"])
                wt, mu = s["wt"], s["mut"]
                if wt not in CHG or mu not in CHG:
                    continue
                dp = {k: PROPS[k][mu] - PROPS[k][wt] for k in PROPS}
                # 双环境（§14）
                ep = res_env_fields(gidx, set(), cross_chain=True)
                ei = res_env_fields(gidx, set(), cross_chain=False)
                f = {"present": 1.0}
                # --- F_MUT（§13）---
                for k, v in dp.items():
                    f[f"d_{k}"] = v
                for k in ("charge", "volume", "hydropathy", "polarity"):
                    f[f"abs_d_{k}"] = abs(dp[k])
                # --- F_ENV_PARTNER 静态（§17）---
                for k, v in ep.items():
                    f[f"pt_{k}"] = v
                # --- F_ENV_INTRA 静态（§18）---
                for k, v in ei.items():
                    f[f"in_{k}"] = v
                # --- 核心：mutation × partner 扰动（§19）---
                if ep.get("present", 0) > 0:
                    f["F_elec_partner"] = dp["charge"] * ep["w_charge_over_d"]
                    f["F_hyd_partner"] = dp["hydropathy"] * ep["w_hydropathy"]
                    f["F_pack_partner"] = dp["volume"] * ep["w_density"]
                    f["F_packabs_partner"] = abs(dp["volume"]) * ep["w_density"]
                    f["F_hbond_partner"] = (dp["donor"] * ep["w_acceptor"]
                                            + dp["acceptor"] * ep["w_donor"])
                    f["F_salt_partner"] = (dp["cation"] * ep["w_anion"]
                                           + dp["anion"] * ep["w_cation"])
                    f["F_arom_partner"] = (dp["aromatic"] * ep["w_aromatic"]
                                           + dp["cation"] * ep["w_aromatic"]
                                           + dp["aromatic"] * ep["w_cation"])
                    f["F_pol_partner"] = dp["polarity"] * ep["w_polarity"]
                    f["F_chg_x_pos"] = dp["charge"] * ep.get("wc_positive", 0.0)
                    f["F_chg_x_neg"] = dp["charge"] * ep.get("wc_negative", 0.0)
                    f["F_hyd_x_pho"] = dp["hydropathy"] * ep.get(
                        "wc_hydrophobic", 0.0)
                    f["F_arom_x_arom"] = dp["aromatic"] * ep.get(
                        "wc_aromatic", 0.0)
                # --- intra 扰动 ---
                if ei.get("present", 0) > 0:
                    f["F_elec_intra"] = dp["charge"] * ei["w_charge_over_d"]
                    f["F_hyd_intra"] = dp["hydropathy"] * ei["w_hydropathy"]
                    f["F_pack_intra"] = dp["volume"] * ei["w_density"]
                    f["F_packabs_intra"] = abs(dp["volume"]) * ei["w_density"]
                    f["F_hbond_intra"] = (dp["donor"] * ei["w_acceptor"]
                                          + dp["acceptor"] * ei["w_donor"])
                    f["F_salt_intra"] = (dp["cation"] * ei["w_anion"]
                                         + dp["anion"] * ei["w_cation"])
                # --- 埋藏/SASA 条件化（§20）---
                ib = float(bur[gidx]) if np.isfinite(bur[gidx]) else 0.0
                sc = float(sasa_c[gidx]) if np.isfinite(sasa_c[gidx]) else np.nan
                f["site_sasa_complex"] = sc
                f["site_burial"] = ib
                f["F_burial_hyd"] = dp["hydropathy"] * ib
                f["F_burial_pol"] = dp["polarity"] * ib
                f["F_burial_chg"] = dp["charge"] * ib
                f["F_burial_vol"] = abs(dp["volume"]) * ib
                # --- 界面状态（§21，冻结定义）---
                dmin = ep.get("min_dist", np.nan)
                if ep.get("present", 0) > 0 and dmin <= 5.0:
                    st = 0
                elif ep.get("present", 0) > 0 and dmin <= 10.0:
                    st = 1
                else:
                    st = 2
                f["iface_state"] = float(st)
                site_feats.append(f)
            if not site_feats:
                continue
            # --- 多突变聚合（§23，冻结）---
            keys = sorted({k for s in site_feats for k in s})
            agg = {}
            for k in keys:
                vals = np.array([s.get(k, np.nan) for s in site_feats],
                                np.float64)
                v = vals[np.isfinite(vals)]
                if len(v) == 0:
                    continue
                if k in ("iface_state",):
                    agg[k] = float(np.mean(v))
                    continue
                agg[f"{k}__sum"] = float(v.sum())
                agg[f"{k}__mean"] = float(v.mean())
                agg[f"{k}__maxabs"] = float(v[np.argmax(np.abs(v))])
            agg["n_mut_sites"] = float(len(site_feats))
            agg["n_iface_sites"] = float(sum(
                1 for s in site_feats if s.get("iface_state", 2) == 0))
            # --- F_EPI_GEOM（§24）---
            gs = [int(s["gidx"]) for s in rec["sites"]]
            if len(gs) > 1:
                dd, same = [], 0
                for i in range(len(gs)):
                    for j in range(i + 1, len(gs)):
                        dd.append(float(np.linalg.norm(xyz[gs[i]] - xyz[gs[j]])))
                        same += int(chains[gs[i]] == chains[gs[j]])
                dd = np.array(dd)
                agg["epi_min_dist"] = float(dd.min())
                agg["epi_mean_dist"] = float(dd.mean())
                agg["epi_n_pairs_8"] = float((dd <= 8.0).sum())
                agg["epi_n_pairs_12"] = float((dd <= 12.0).sum())
                agg["epi_samechain_frac"] = float(same / len(dd))
            else:
                agg["epi_min_dist"] = np.nan
                agg["epi_mean_dist"] = np.nan
                agg["epi_n_pairs_8"] = 0.0
                agg["epi_n_pairs_12"] = 0.0
                agg["epi_samechain_frac"] = np.nan
            F[rid] = agg
        n_done += 1
        if (ti + 1) % 20 == 0:
            emit(f"  - [{ti + 1}] {t}: n_rows={len(F)} 累计, "
                 f"{(time.time() - t0) / 60:.1f} min")

    emit("")
    emit(f"- 完成 **{n_done}** targets / **{len(F)}** rows")
    # ---------- 组装矩阵 ----------
    rids = sorted(F.keys())
    all_keys = sorted({k for r in rids for k in F[r]})
    FM = {k: np.full(len(rids), np.nan, np.float32) for k in all_keys}
    for i, r in enumerate(rids):
        for k, v in F[r].items():
            FM[k][i] = v
    PH = pd.DataFrame({"row_id": rids})
    for k in all_keys:
        PH[k] = FM[k]
    # 拼接 Base 与 BLOCK L（按 row_id）
    bi = pos.reindex(rids).to_numpy()
    ok = np.isfinite(bi)
    base_m = np.full((len(rids), BASE_ALL.shape[1]), np.nan, np.float32)
    base_m[ok] = BASE_ALL[bi[ok].astype(int)]
    for j in range(BASE_ALL.shape[1]):
        PH[f"base_{j}"] = base_m[:, j]
    bl = BL_i.reindex(rids)
    for c in bl_cols:
        PH[f"bl_{c}"] = bl[c].to_numpy(np.float32)
    PH.to_parquet(FEAT / "physics_features.parquet", index=False)
    # feature dict
    fd = []
    for k in all_keys:
        kind = ("F_MUT" if k.startswith(("d_", "abs_d_")) else
                "F_ENV_PARTNER" if k.startswith("pt_") else
                "F_ENV_INTRA" if k.startswith("in_") else
                "F_ELEC" if "elec" in k else
                "F_HYD" if "hyd" in k else
                "F_PACK" if "pack" in k else
                "F_HBOND" if "hbond" in k else
                "F_SALT" if "salt" in k else
                "F_AROM" if "arom" in k else
                "F_BURIAL" if "burial" in k or "sasa" in k else
                "F_EPI_GEOM" if k.startswith("epi_") else
                "F_IFACE" if "iface" in k else "OTHER")
        fd.append({"feature": k, "family": kind,
                   "agg": (k.split("__")[-1] if "__" in k else "")})
    FDD = pd.DataFrame(fd)
    FDD.to_csv(PH3 / "feature_dict.csv", index=False, encoding="utf-8")
    emit("")
    emit("## 3.2 物理特征族统计")
    emit("")
    emit("| family | n_features |")
    emit("|---|---:|")
    for fam, g in FDD.groupby("family"):
        emit(f"| `{fam}` | {len(g)} |")
    emit(f"| **合计物理特征** | **{len(FDD)}** |")
    emit("")
    emit(f"- 输出：`features/v7_3/physics_features.parquet` "
         f"（{PH.shape[0]} 行 × {PH.shape[1]} 列）")
    emit(f"- 含 `base_*`（历史 Base）与 `bl_*`（STATIC_LOCAL_CONTROL）")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH3 / "phase3_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 physics_features.parquet + aa_physics_table.csv + feature_dict.csv")
    print(f">>> rows={PH.shape[0]} cols={PH.shape[1]} 物理特征={len(FDD)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
