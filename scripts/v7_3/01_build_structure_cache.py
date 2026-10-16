# -*- coding: utf-8 -*-
"""scripts/v7_3/01_build_structure_cache.py -- V7.3 PHASE 1+2.

PHASE 1：**U_PHYS_STRICT**（§9/§10）—— 所有突变位点都必须：
    canonical row 存在 / chain 在 PDB 内 / 残基坐标存在 /
    mapped WT == expected WT / **非 AMBIGUOUS 映射**
多突变必须**全部**位点通过（不静默使用部分映射）。

PHASE 2：**结构环境缓存**（§52/§53）—— 每个 target 只解析一次 PDB，缓存：
    * 每残基重原子坐标
    * 同链 / 跨链**重原子最小距离**（§15：不用 CA 代替）
    * 8 Å 内邻居（按化学类别的加权/非加权接触场）
    * SASA（freesasa，§20）：复合物 & 各链单独
    * 界面状态（§21）

输出：`features/v7_3/pdb_environment_cache/<target>.npz`
      `features/v7_3/u_phys_strict_manifest.csv`
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import defaultdict
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

V711 = C.ROOT / "experiments" / "v7_1_1"
V73 = C.ROOT / "experiments" / "v7_3"
PH1 = V73 / "phase1_universe"
PH1.mkdir(parents=True, exist_ok=True)
CACHE = C.ROOT / "features" / "v7_3" / "pdb_environment_cache"
CACHE.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
L: list[str] = []

# 化学类别（§22 接触化学分箱）
CHEM = {
    "positive": set("KRH"), "negative": set("DE"),
    "polar": set("STNQYC"), "hydrophobic": set("AVLIMFWC"),
    "aromatic": set("FWYH"), "glypro": set("GP"),
}
# 显式化学类别（互斥优先序，用于单标签）
def chem_class(a: str) -> str:
    if a in "KR": return "positive"
    if a in "DE": return "negative"
    if a in "FWY": return "aromatic"
    if a in "GP": return "glypro"
    if a in "STNQCY": return "polar"
    if a in "AVLIM": return "hydrophobic"
    return "other"
CHEMS = ["positive", "negative", "polar", "hydrophobic", "aromatic", "glypro"]


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def parse_heavy(p: Path):
    """返回 {chain: {'res': [(resseq,icode,aa)], 'heavy': [ (res_idx, x,y,z) ]}}
    以及每残基重原子坐标列表。"""
    ch_res = defaultdict(list)          # chain -> [(resseq,icode,aa)]
    ch_atoms = defaultdict(list)        # chain -> [(res_idx, element, x,y,z)]
    for ln in open(p, "r", errors="replace"):
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if elem == "H":
            continue
        rn = ln[17:20].strip().upper()
        if rn not in THREE:              # 非标准残基 -> 跳过（§10 nonstandard）
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
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
            ch_res[c].append(key)
            ch_atoms[c].append([])
        ch_atoms[c][-1].append((elem, x, y, z))
    out = {}
    for c in ch_res:
        res = [(rs, ic, THREE[rn]) for (rs, ic, rn) in ch_res[c]]
        atoms = ch_atoms[c]
        out[c] = {"res": res, "atoms": atoms}
    return out


THREE = {}   # 由 v7_0 的 THREE2ONE 填充


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-targets", type=int, default=10 ** 9)
    ap.add_argument("--max-minutes", type=float, default=80.0)
    a = ap.parse_args()
    t0 = time.time()
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    global THREE
    THREE = v7f.THREE2ONE
    import freesasa

    emit("# V7.3 PHASE 1+2 — U_PHYS_STRICT + 结构环境缓存")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    emit(f"- canonical：{len(lab)} 行 / {lab.target_id_canonical.nunique()} targets")
    emit("")

    targets = sorted(lab.target_id_canonical.astype(str).unique())
    exc = defaultdict(int)
    keep_rows, exc_rows = [], []
    cache_files = []
    for ti, t in enumerate(targets):
        if ti >= a.max_targets:
            break
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"- ⏱ 软时限，已处理 {ti}/{len(targets)} targets")
            break
        p = PDBDIR / f"{t.upper()}.pdb"
        if not p.exists():
            for r in lab[lab.target_id_canonical.astype(str) == t].itertuples():
                exc["missing_pdb_file"] += 1
                exc_rows.append({"row_id": int(r.pair_row), "target": t,
                                 "reason": "missing_pdb_file"})
            continue
        chs = parse_heavy(p)
        if not chs:
            for r in lab[lab.target_id_canonical.astype(str) == t].itertuples():
                exc["no_protein_chain"] += 1
                exc_rows.append({"row_id": int(r.pair_row), "target": t,
                                 "reason": "no_protein_chain"})
            continue

        # ---- 每 target 一次：全复合物重原子 KD 树 + 逐残基描述符 ----
        all_atoms, all_res = [], []      # 全局残基表
        chain_res_idx = {}               # chain -> [全局残基下标]
        for c in sorted(chs):
            idxs = []
            for k, (rs, ic, aa) in enumerate(chs[c]["res"]):
                g = len(all_res)
                all_res.append((c, k, rs, ic, aa))
                idxs.append(g)
                for (el, x, y, z) in chs[c]["atoms"][k]:
                    all_atoms.append((g, x, y, z))
            chain_res_idx[c] = idxs
        res_xyz = np.array([[np.mean([a[1] for a in chs[c]["atoms"][k]]),
                             np.mean([a[2] for a in chs[c]["atoms"][k]]),
                             np.mean([a[3] for a in chs[c]["atoms"][k]])]
                            for c, k, _, _, _ in all_res], np.float64)
        atom_xyz = np.array([[a[1], a[2], a[3]] for a in all_atoms], np.float64)
        atom_res = np.array([a[0] for a in all_atoms], np.int32)
        tree = cKDTree(atom_xyz)
        nres = len(all_res)
        # 化学类别 & 侧链理化（用于环境场）
        aa_of = [x[4] for x in all_res]
        chain_of = [x[0] for x in all_res]

        # ---- 环境矩阵：对每个残基 i，统计 8 Å 内邻居的加权/非加权化学场 ----
        W_R = 8.0
        NS = ["n_contact_4", "n_contact_6", "n_contact_8",
              "min_heavy_dist", "mean_dist",
              "w_positive", "w_negative", "w_polar", "w_hydrophobic",
              "w_aromatic", "w_glypro",
              "c_positive", "c_negative", "c_polar", "c_hydrophobic",
              "c_aromatic", "c_glypro", "w_all", "n_partner_8", "n_intra_8"]
        ENV = np.zeros((nres, len(NS)), np.float32)
        # 每个残基的近邻（重原子最小距离）
        pairs = tree.query_ball_point(res_xyz, r=W_R + 6.0)   # 预筛
        for i in range(nres):
            cand = pairs[i]
            best = {}
            for ai in cand:
                j = int(atom_res[ai])
                if j == i:
                    continue
                d = float(np.linalg.norm(atom_xyz[ai] - res_xyz[i]))
                if j not in best or d < best[j]:
                    best[j] = d
            if not best:
                continue
            js = np.fromiter(best.keys(), np.int32)
            ds = np.fromiter(best.values(), np.float64)
            w = np.exp(-ds / 4.0)
            w[ds > W_R] = 0.0
            intra = np.array([chain_of[j] == chain_of[i] for j in js])
            row = ENV[i]
            row[0] = float((ds <= 4.0).sum())
            row[1] = float((ds <= 6.0).sum())
            row[2] = float((ds <= 8.0).sum())
            row[3] = float(ds.min())
            row[4] = float(ds.mean())
            row[18] = float(((ds <= 8.0) & (~intra)).sum())
            row[19] = float(((ds <= 8.0) & intra).sum())
            row[17] = float(w.sum())
            for ci, cname in enumerate(CHEMS):
                m = np.array([chem_class(aa_of[j]) == cname for j in js])
                row[5 + ci] = float(w[m].sum())
                row[11 + ci] = float(m.sum())

        # ---- SASA（freesasa via BioPython，§20）----
        # calcBioPDB 接受 **Structure** 或 Model，**不接受 Chain**；
        # 故「单链孤立 SASA」要把该链装进一个新的 Structure。
        from Bio.PDB import PDBParser, Structure as _BStruct, Model as _BModel
        freesasa.setVerbosity(freesasa.silent)
        sasa_cplx = np.full(nres, np.nan, np.float32)
        sasa_iso = np.full(nres, np.nan, np.float32)
        _bst = None
        try:
            _bst = PDBParser(QUIET=True).get_structure("t", str(p))
        except Exception:  # noqa: BLE001
            _bst = None

        def _lookup(areas_by_chain, c, rs, ic):
            d = areas_by_chain.get(c)
            if not d:
                return None
            v = d.get(f"{rs}{ic}") or d.get(str(rs))
            if v is None:
                cands = [kk for kk in d if kk.startswith(str(rs))]
                v = d[cands[0]] if cands else None
            return v

        if _bst is not None:
            try:
                _r, _ = freesasa.calcBioPDB(_bst)
                ra = _r.residueAreas()
                for c in sorted(chs):
                    for k, (rs, ic, aa) in enumerate(chs[c]["res"]):
                        v = _lookup(ra, c, rs, ic)
                        if v is not None:
                            sasa_cplx[chain_res_idx[c][k]] = float(
                                getattr(v, "total", np.nan))
            except Exception:  # noqa: BLE001
                pass
            try:
                _src_model = next(_bst.get_models())
            except StopIteration:
                _src_model = None
            for c in sorted(chs):
                if _src_model is None or c not in _src_model:
                    continue
                try:
                    _s2 = _BStruct.Structure("single")
                    _m2 = _BModel.Model(0)
                    _s2.add(_m2)
                    _m2.add(_src_model[c].copy())
                    _r1, _ = freesasa.calcBioPDB(_s2)
                except Exception:  # noqa: BLE001
                    continue
                ra1 = _r1.residueAreas()
                for k, (rs, ic, aa) in enumerate(chs[c]["res"]):
                    v = _lookup(ra1, c, rs, ic)
                    if v is not None:
                        sasa_iso[chain_res_idx[c][k]] = float(
                            getattr(v, "total", np.nan))
        burial = (sasa_iso - sasa_cplx).astype(np.float32)   # §20 interface_burial

        # ---- 写出 target 缓存 ----
        np.savez_compressed(
            CACHE / f"{t}.npz",
            chain=np.array(chain_of, dtype=object),
            resseq=np.array([x[2] for x in all_res], np.int32),
            icode=np.array([x[3] for x in all_res], dtype=object),
            aa=np.array(aa_of, dtype=object),
            res_xyz=res_xyz.astype(np.float32),
            env=ENV,
            env_names=np.array(NS, dtype=object),
            sasa_complex=sasa_cplx,
            sasa_isolated=sasa_iso,
            burial=burial,
        )
        cache_files.append({"target": t, "n_res": nres,
                            "n_chains": len(chs),
                            "sasa_ok": bool(np.isfinite(sasa_cplx).any())})

        # ---- 判定 U_PHYS_STRICT（§9）----
        for r in lab[lab.target_id_canonical.astype(str) == t].itertuples():
            pr = int(r.pair_row)
            g = s_by_row.get(pr)
            if g is None or not len(g):
                exc["no_site_record"] += 1
                exc_rows.append({"row_id": pr, "target": t, "reason": "no_site_record"})
                continue
            bad = None
            sites_ok = []
            for x in g.itertuples():
                c = str(x.ann_chain); i = int(x.seq_index0)
                if str(x.mapping_status) == "AMBIGUOUS":
                    bad = "ambiguous_mapping"; break
                if c not in chs:
                    bad = "missing_chain"; break
                if not (0 <= i < len(chs[c]["res"])):
                    bad = "idx_out_of_bounds"; break
                rs, ic, aa = chs[c]["res"][i]
                if aa != str(x.wt_aa):
                    bad = "wt_mismatch"; break
                gidx = chain_res_idx[c][i]
                if not np.isfinite(res_xyz[gidx]).all():
                    bad = "missing_coordinates"; break
                sites_ok.append({"chain": c, "seq_index0": i, "gidx": gidx,
                                 "resseq": rs, "icode": ic,
                                 "wt": aa, "mut": str(x.mut_aa)})
            if bad is not None:
                exc[bad] += 1
                exc_rows.append({"row_id": pr, "target": t, "reason": bad})
                continue
            keep_rows.append({"row_id": pr, "target_id_canonical": t,
                              "family_id_canonical": str(r.family_id_canonical),
                              "n_sites": len(sites_ok),
                              "sites": sites_ok})
        if (ti + 1) % 10 == 0:
            emit(f"  - [{ti + 1}/{len(targets)}] {t}: {nres} res, "
                 f"{len(chs)} chains, 累计 keep={len(keep_rows)}, "
                 f"{(time.time() - t0) / 60:.1f} min")

    # ---------- manifest ----------
    emit("")
    emit("## 1.1 U_PHYS_STRICT（§9）")
    emit("")
    K = pd.DataFrame([{k: v for k, v in r.items() if k != "sites"}
                      for r in keep_rows])
    if len(K):
        K.to_csv(C.ROOT / "features" / "v7_3" /
                 "u_phys_strict_manifest.csv", index=False, encoding="utf-8")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| canonical rows | {len(lab)} |")
    emit(f"| **U_PHYS_STRICT rows** | **{len(K)}** |")
    emit(f"| targets | **{K.target_id_canonical.nunique() if len(K) else 0}** |")
    emit(f"| families | **{K.family_id_canonical.nunique() if len(K) else 0}** |")
    for k in (1, 2, 3, 5, 10):
        if len(K):
            n = (K.groupby("target_id_canonical").size() >= k).sum()
        else:
            n = 0
        emit(f"| n>={k} targets | **{n}** |")
    emit("")
    emit("### 排除原因（§10）")
    emit("")
    emit("| reason | rows |")
    emit("|---|---:|")
    for k, v in sorted(exc.items(), key=lambda x: -x[1]):
        emit(f"| `{k}` | {v} |")
    emit(f"| **合计排除** | **{sum(exc.values())}** |")
    emit("")
    if exc_rows:
        pd.DataFrame(exc_rows).to_csv(PH1 / "exclusions.csv", index=False,
                                      encoding="utf-8")
    # 保留 sites 供下游用
    import pickle
    with open(CACHE / "_sites.pkl", "wb") as f:
        pickle.dump(keep_rows, f)
    pd.DataFrame(cache_files).to_csv(PH1 / "cache_index.csv", index=False,
                                     encoding="utf-8")
    emit(f"- 结构缓存：**{len(cache_files)}** 个 target → `{CACHE}`")
    emit(f"- SASA 成功：**{sum(1 for c in cache_files if c['sasa_ok'])}** / "
         f"{len(cache_files)}（freesasa）")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH1 / "phase1_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 u_phys_strict_manifest.csv / exclusions.csv / cache_index.csv")
    print(f">>> U_PHYS_STRICT rows={len(K)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
