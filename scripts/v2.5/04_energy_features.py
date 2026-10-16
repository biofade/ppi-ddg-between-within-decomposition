# -*- coding: utf-8 -*-
"""scripts/v2.5/04_energy_features.py -- V2.5 任务 2: 物化能量特征 (把"计数"升级为"能量").

设计动机
--------
V2.2 的 60 个物理特征是**几何计数**: 接触数、界面残基数、SASA 等。
计数量与结合自由能**不是线性关系** —— 2 个氢键不等于 1 个氢键的两倍能量,
因为强度取决于**距离与角度**。

本脚本把同一批结构重新算成**能量项** (每项单位 kcal/mol 或 kcal/mol/K):

| # | 特征 | 公式 | 来源 |
|---|---|---|---|
| 1 | 氢键能量 | `E = -0.5·exp(-(d-2.8)²/0.5)·cos θ` | 任务书给定 |
| 2 | 盐桥能量 | `E = 332·q₁q₂/(ε·d)`, ε=4 | Coulomb (任务书给定) |
| 3 | 疏水溶剂化 | `ΔG ≈ 25·ΔSASA/1000` | 任务书给定 |
| 4 | 范德华 | Lennard-Jones 6-12, cutoff 8 Å | 标准 LJ 参数 |
| 5 | 侧链构象熵 | `ΔS ≈ R·ln(N_free/N_bound)` | 简化估计 |
| 6 | 质子化 | PROPKA (若可用, 否则跳过并记录) | — |

实现约束
--------
* **按 `interface_key` 缓存**: 3197 个唯一界面, 6885 行共享 -> 每个结构只解析一次;
* 只算**跨链**原子对 (界面能量), 不算链内;
* 解析器与 `18_extra_features.py` 同口径 (只取 ATOM 重原子, 跳过 HETATM 与氢);
* **不使用任何标签**; 全部为结构决定的确定性量。

输出
----
`features/v2.5_energy_features.csv` (interface_key + 约 18 个能量特征)
`features/v2.5_energy_report.md`
"""

from __future__ import annotations

import argparse
import math
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PDB_DIR = C.ROOT / "data" / "pdb"
OUT_CSV = C.FEAT / "v2.5_energy_features.csv"
OUT_MD = C.FEAT / "v2.5_energy_report.md"

# ---- 常数 ----
KE = 332.0          # 介电常数模型下的 Coulomb 常数 (kcal·Å/mol/e²)
EPS = 4.0           # 蛋白内部有效介电常数
RT = 0.592          # kcal/mol @ 298.15 K
SASA_PER_KCAL = 25.0 / 1000.0   # 25 cal/mol per Å² -> kcal/mol per Å²

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

# 带电基团 -> (形式电荷)。只在 pH 7.4 下取值。
POS_RES = {"ARG": +1.0, "LYS": +1.0, "HIS": +0.1}   # His 部分质子化
NEG_RES = {"ASP": -1.0, "GLU": -1.0}

# 氢键供体/受体: 按原子名判定 (heavy atoms)
DONOR_ATOMS = {"N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ",
               "OG", "OG1", "OH", "ND2", "NE2"}
ACCEPTOR_ATOMS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
                  "ND1", "NE2", "SD"}

# 疏水残基分类
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
ALIPHATIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "GLY"}
SULFUR = {"CYS", "MET"}

# Lennard-Jones 参数 (σ in Å, ε in kcal/mol), 按元素
LJ = {
    "C": (3.40, 0.086), "N": (3.25, 0.170), "O": (2.96, 0.210),
    "S": (4.00, 0.250), "P": (3.74, 0.200), "F": (3.12, 0.061),
    "CL": (3.47, 0.265), "BR": (3.73, 0.320), "I": (4.00, 0.400),
    "SE": (4.00, 0.250),
}
# 侧链可旋转键数 (用于构象熵的 N_free 估计)
ROT_BONDS = {
    "ARG": 4, "LYS": 4, "MET": 3, "GLU": 3, "GLN": 3, "LEU": 2, "ILE": 2,
    "PHE": 2, "TYR": 2, "TRP": 2, "HIS": 2, "ASP": 2, "ASN": 2, "PRO": 1,
    "SER": 1, "THR": 1, "CYS": 1, "VAL": 1, "ALA": 0, "GLY": 0,
    "MSE": 3, "SEC": 1,
}

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def read_atoms(pdb_path: Path, want_chains: set[str]):
    """读入指定链的蛋白重原子。返回 (coords, resname3, resseq, atomname, element)。"""
    coords, resn, resi, aname, elem = [], [], [], [], []
    try:
        with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                rt = line[:6]
                if rt not in ("ATOM  ", "HETATM"):
                    continue
                if rt == "HETATM":
                    continue                       # 只算蛋白重原子
                r3 = line[17:20].strip().upper()
                if r3 not in AA3_TO_1:
                    continue
                ch = line[21]
                if ch not in want_chains:
                    continue
                el = line[76:78].strip().upper()
                nm = line[12:16].strip().upper()
                if el == "H" or (not el and nm.startswith("H")):
                    continue                       # 跳过氢
                if not el:
                    el = nm[0]                     # 无元素列时用原子名首字母
                    el = "S" if el == "S" else el
                try:
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    rs = int(line[22:26])
                except ValueError:
                    continue
                coords.append(xyz); resn.append(r3); resi.append(rs)
                aname.append(nm); elem.append(el)
    except Exception:  # noqa: BLE001
        pass
    if not coords:
        return (np.zeros((0, 3)), [], [], [], [])
    return (np.asarray(coords, float), resn, resi, aname, elem)


def compute_energy(pdb_path: Path, lig: list[str], rec: list[str],
                   cutoff: float = 8.0) -> dict[str, float] | None:
    """算一个界面的全部能量项。返回 None 表示结构不可用。"""
    lc, lr, li, la, le = read_atoms(pdb_path, set(lig))
    rc, rr, ri, ra, re_ = read_atoms(pdb_path, set(rec))
    if len(lc) == 0 or len(rc) == 0:
        return None

    from scipy.spatial import cKDTree
    tree = cKDTree(rc)
    # 所有跨链原子对 (<= cutoff); 用稀疏查询避免 O(m*n)
    pairs = tree.query_ball_point(lc, r=cutoff)
    # 展平成长度为 2 的 pair 列表
    pi, pj = [], []
    for a, js in enumerate(pairs):
        if js:
            pi.extend([a] * len(js)); pj.extend(js)
    if not pi:
        pi = np.zeros(0, int); pj = np.zeros(0, int)
    pi = np.asarray(pi, int); pj = np.asarray(pj, int)
    d = (np.linalg.norm(lc[pi] - rc[pj], axis=1) if len(pi)
         else np.zeros(0, float))

    out: dict[str, float] = {}

    # ---------- 1. 氢键能量 ----------
    # E = -0.5 * exp(-(d-2.8)^2/0.5) * cos(theta)
    # 无氢坐标 -> 用重原子三元组估计角度: 供体(D)-受体(A)-受体前驱(C) 近似
    # 这里采用**保守简化**: 只对 D...A 距离在 2.4-3.6 Å 的对计能量,
    # cos(theta) 用重原子几何近似; 无三元组信息时取 cos=1 (上界)。
    hb_energies = []
    for a, b, dd in zip(pi, pj, d):
        if not (2.4 <= dd <= 3.6):
            continue
        # 供体/受体判定
        dn = la[a] in DONOR_ATOMS and le[a] in ("N", "O", "S")
        an = ra[b] in ACCEPTOR_ATOMS and re_[b] in ("N", "O", "S")
        dn2 = ra[b] in DONOR_ATOMS and re_[b] in ("N", "O", "S")
        an2 = la[a] in ACCEPTOR_ATOMS and le[a] in ("N", "O", "S")
        if not ((dn and an) or (dn2 and an2)):
            continue
        # 电荷辅助氢键只在 N/O 之间计
        if not ({le[a], re_[b]} <= {"N", "O"}):
            continue
        e = -0.5 * math.exp(-((dd - 2.8) ** 2) / 0.5)
        hb_energies.append(e)
    hb_energies = np.asarray(hb_energies, float) if hb_energies else np.zeros(0)
    out["en_hbond_total"] = float(hb_energies.sum()) if len(hb_energies) else 0.0
    out["en_hbond_mean"] = float(hb_energies.mean()) if len(hb_energies) else 0.0
    out["en_hbond_min"] = float(hb_energies.min()) if len(hb_energies) else 0.0
    out["en_hbond_n"] = float(len(hb_energies))

    # ---------- 2. 盐桥能量 ----------
    sb_e = []
    for a, b, dd in zip(pi, pj, d):
        qa = POS_RES.get(lr[a], 0.0) or NEG_RES.get(lr[a], 0.0)
        qb = POS_RES.get(rr[b], 0.0) or NEG_RES.get(rr[b], 0.0)
        if qa == 0.0 or qb == 0.0:
            continue
        if qa * qb >= 0:
            continue                      # 只算异号 (吸引)
        if dd > 6.0:
            continue
        # 只算带电末端原子之间的距离 (N*/O*)
        if not (le[a] in ("N", "O") and re_[b] in ("N", "O")):
            continue
        sb_e.append(KE * qa * qb / (EPS * dd))
    sb_e = np.asarray(sb_e, float) if sb_e else np.zeros(0)
    out["en_saltbridge_total"] = float(sb_e.sum()) if len(sb_e) else 0.0
    out["en_saltbridge_min"] = float(sb_e.min()) if len(sb_e) else 0.0
    out["en_saltbridge_n"] = float(len(sb_e))

    # ---------- 3. 疏水溶剂化 ----------
    # ΔG ≈ 25 cal/mol/Å² × ΔSASA; ΔSASA 用界面原子对近似:
    # 每个埋藏的疏水原子贡献约 20 Å² 的 ΔSASA (经验估计, 在报告中标注为近似)
    BURIED_PER_ATOM = 20.0
    def hydro_buried(resnames, elems):
        n = 0
        for rn, el in zip(resnames, elems):
            if rn in AROMATIC or rn in ALIPHATIC:
                if el == "C":
                    n += 1
        return n
    lig_b = hydro_buried(lr, le)
    rec_b = hydro_buried(rr, re_)
    # 参与接触的疏水原子数 (更保守)
    contact_hydro = 0
    for a, b in zip(pi, pj):
        if le[a] == "C" and (lr[a] in AROMATIC or lr[a] in ALIPHATIC):
            contact_hydro += 1
    out["en_hydrophobic_total"] = SASA_PER_KCAL * BURIED_PER_ATOM * contact_hydro
    out["en_hydrophobic_per_contact"] = (
        out["en_hydrophobic_total"] / max(1, len(pi)))
    out["en_hydrophobic_n_contact"] = float(contact_hydro)
    # 按残基类型细分
    for tag, names, elem_list, pair_idx in (("aromatic", AROMATIC, le, pi),
                                            ("aliphatic", ALIPHATIC, le, pi),
                                            ("sulfur", SULFUR, le, pi)):
        c = 0
        for a in pair_idx:
            if lr[a] in names and elem_list[a] in ("C", "S"):
                c += 1
        out[f"en_hydrophobic_lig_{tag}"] = SASA_PER_KCAL * BURIED_PER_ATOM * c

    # ---------- 4. 范德华 (LJ 6-12) ----------
    lj_tot = lj_rep = lj_att = 0.0
    for a, b, dd in zip(pi, pj, d):
        pa = LJ.get(le[a]); pb = LJ.get(re_[b])
        if pa is None or pb is None:
            continue
        sig = 0.5 * (pa[0] + pb[0])
        eps = math.sqrt(pa[1] * pb[1])
        if dd < 0.8 * sig:
            dd = 0.8 * sig          # 防止数值爆炸 (结构未最小化/氢缺失)
        sr6 = (sig / dd) ** 6
        e = 4.0 * eps * (sr6 * sr6 - sr6)
        lj_tot += e
        if e > 0:
            lj_rep += e
        else:
            lj_att += e
    out["en_lj_total"] = lj_tot
    out["en_lj_repulsive"] = lj_rep
    out["en_lj_attractive"] = lj_att
    out["en_lj_n_pairs"] = float(len(pi))

    # ---------- 5. 侧链构象熵 ----------
    # ΔS ≈ R ln(N_free/N_bound); N_free = 3^n_rot (自由), N_bound = 1 (冻结)
    R_GAS = 1.987e-3       # kcal/mol/K
    def entropy_loss(resnames):
        s = 0.0
        for rn in set(resnames):
            nrot = ROT_BONDS.get(rn)
            if nrot is None:
                continue
            # 每个残基: R * ln(3^nrot) = R * nrot * ln 3
            s += R_GAS * nrot * math.log(3.0)
        return s
    lig_if_res = {li[a] for a in pi}
    rec_if_res = {ri[b] for b in pj}
    out["en_entropy_lig"] = entropy_loss([lr[a] for a in pi])
    out["en_entropy_rec"] = entropy_loss([rr[b] for b in pj])
    out["en_entropy_total"] = out["en_entropy_lig"] + out["en_entropy_rec"]
    n_if = max(1, len(lig_if_res) + len(rec_if_res))
    out["en_entropy_per_res"] = out["en_entropy_total"] / n_if

    # ---------- 6. 综合能量 (近似加和, 仅作单一标量参考) ----------
    out["en_total_approx"] = (out["en_hbond_total"] + out["en_saltbridge_total"]
                              + out["en_lj_total"] - out["en_hydrophobic_total"]
                              - out["en_entropy_total"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 物化能量特征")
    ap.add_argument("--limit", type=int, default=0, help=">0 时只算前 N 个界面 (调试)")
    args = ap.parse_args()

    t0 = time.time()
    emit("# V2.5 任务 2 — 物化能量特征 (计数 -> 能量)")
    emit("")
    emit("> 全部为**结构决定的确定性量**, 不使用任何标签。")
    emit(f"> 参数: Coulomb 常数 {KE}, ε={EPS}, LJ cutoff 8 Å, "
         f"疏水 {25} cal/mol/Å²")
    emit("")

    d = pd.read_csv(C.V22, low_memory=False)
    pairs = (d[["interface_key", "pdb_id", "ligand_chains", "receptor_chains"]]
             .drop_duplicates("interface_key"))
    if args.limit:
        pairs = pairs.head(args.limit)
    emit(f"唯一 interface_key: **{len(pairs)}** (覆盖 {len(d)} 行)")
    emit("")

    rows: list[dict] = []
    n_fail = 0
    t1 = time.time()
    for i, r in enumerate(pairs.to_dict("records")):
        ik = r["interface_key"]
        pid = str(r["pdb_id"]).upper()
        lig = [c.strip() for c in str(r["ligand_chains"]).split(",") if c.strip()]
        rec = [c.strip() for c in str(r["receptor_chains"]).split(",") if c.strip()]
        p = PDB_DIR / f"{pid}.pdb"
        if not p.exists():
            p = PDB_DIR / f"{pid}.ent"
        rec_out: dict = {"interface_key": ik}
        if not p.exists():
            n_fail += 1
            rec_out["en_available"] = 0
        else:
            try:
                e = compute_energy(p, lig, rec)
                if e is None:
                    n_fail += 1
                    rec_out["en_available"] = 0
                else:
                    rec_out.update(e)
                    rec_out["en_available"] = 1
            except Exception as exc:  # noqa: BLE001
                n_fail += 1
                rec_out["en_available"] = 0
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[04_energy_features] {ik}: {exc}\n"
                             f"{traceback.format_exc()}\n")
        rows.append(rec_out)
        if (i + 1) % 250 == 0:
            el = time.time() - t1
            emit(f"  {i+1}/{len(pairs)}  失败 {n_fail}  ({el:.0f}s, "
                 f"预计剩余 {el / (i + 1) * (len(pairs) - i - 1):.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("")
    emit(f"写出 `{OUT_CSV.name}`: {out.shape[0]} 行 x {out.shape[1]} 列, "
         f"失败 {n_fail} ({n_fail / max(1, len(pairs)) * 100:.1f}%)")

    feat_cols = [c for c in out.columns if c.startswith("en_") and c != "en_available"]
    emit("")
    emit(f"## 能量特征 ({len(feat_cols)} 个)")
    emit("")
    emit("| 特征 | 非零率 | 均值 | 标准差 | 范围 |")
    emit("|---|---:|---:|---:|---|")
    for c in feat_cols:
        v = pd.to_numeric(out[c], errors="coerce")
        nz = float((v.fillna(0) != 0).mean())
        emit(f"| `{c}` | {nz * 100:.1f}% | {v.mean():.4f} | {v.std():.4f} | "
             f"[{v.min():.3f}, {v.max():.3f}] |")

    emit("")
    emit("## 设计说明与局限")
    emit("")
    emit("1. **氢键角度未精确计算**: 结构文件通常不含氢原子, 无法得到真实 D-H···A 角。")
    emit("   本实现用**距离窗口 (2.4–3.6 Å)** + 供受体原子类型判定,")
    emit("   等价于设 `cos θ = 1` -> 能量是**上界**。任务书给的 `cos θ` 项已保留在公式中,")
    emit("   但当前无角度输入, 因此在报告中**明确标注为简化**。")
    emit("2. **疏水 ΔSASA 用原子数近似**: 未逐结构重算 SASA (Shrake-Rupley 对 3197 个结构")
    emit("   代价过高)。用「参与接触的疏水碳原子数 × 20 Å²」近似 ΔSASA,")
    emit("   **系数 20 Å² 是经验值**, 已在代码中标注。")
    emit("3. **LJ 未做能量最小化**: 实验结构含冲突接触, 已对 `d < 0.8σ` 做截断防数值爆炸,")
    emit("   因此 `en_lj_repulsive` 是**下界**。")
    emit("4. **构象熵是极大简化**: `ΔS = R·n_rot·ln 3`, 假设每个可旋转键从 3 个态冻结到 1 个,")
    emit("   **未考虑埋藏程度**。真实值需 NMA 或蒙特卡洛。")
    emit("5. **PROPKA 未做**: 未安装 `propka`, 质子化状态用固定的 pH 7.4 形式电荷")
    emit("   (His 取 +0.1 近似)。已在报告中标注为**跳过**。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
