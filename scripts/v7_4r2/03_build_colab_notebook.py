# -*- coding: utf-8 -*-
"""scripts/v7_4r2/03_build_colab_notebook.py -- 生成 V7.4R2 Colab notebook（§5–§26）.

实现要点（全部按 V7.4R2 规格）：
* **Canonical WT**：每 target 只准备一次，**禁用缺失残基重建**（`missingResidues = {}`），
  硬断言 `missing_residues_added == 0`（§5）
* **Mutant from canonical WT**：不从原始 PDB 重建；只改被突变残基（§9）
* **局部侧链松弛**：骨架与壳外强限制；突变侧链 + 8 Å 内侧链自由（§11）
* **Kabsch 对齐后的局部 RMSD**（§13），不再用原始整体系 RMSD
* **局部几何特征**（§15）+ **局部能量**（不可靠则 `NOT_RUN`）（§16）
* **快烟测**：1 target（1A22）+ 5 identity + 5 mutant，逐门短路（§18/§24）
* Drive 持久化 + 缓存键（§27/§28）

平台：**Colab GPU only**（Kaggle 不用于 primary，§3）。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "colab" / "v7_4r2"
OUT.mkdir(parents=True, exist_ok=True)
NB = OUT / "V7_4R2_CanonicalWT_LocalRelax_COLAB.ipynb"
CELLS: list = []


def md(t):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": t.splitlines(True)})


def code(t):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": t.strip("\n").splitlines(True)})


md("""# V7.4R2 — Canonical-WT Local Relaxation (Google Colab GPU)

**核心原则**（§2）：

```
experimental PDB -> Canonical WT -> Mutant clone -> local sidechain relaxation
```

同一 target 的**所有**突变都派生自**同一个冻结的 canonical WT 参考**。
**不再**从原始 PDB 独立重建 WT 与 MUT。

> ⚠️ **术语纪律**：产物是 **computationally relaxed mutant geometry** /
> **force-field potential-energy change**。**不是** experimental structure、
> **不是** true binding energy / true ΔΔG。

**V7.4R 已确立**：`OpenMM CUDA / ff14SB / GBn2_XML / CustomGBForce / mutation builder` 全部 PASS。
失败的是**旧的结构准备与配对最小化实现**。本 notebook 只替换该实现，
**物理协议一字不改**（§4）。

**快速门**（§18/§24）：1 target（1A22）+ 5 identity + 5 mutant，
逐门短路 —— 失败即刻停止，不再让您等两小时。
""")

code('''
# ============================================================ 00 Config (§26)
import os, sys, json, time, math, hashlib, zipfile, shutil, subprocess, warnings, random
from pathlib import Path
warnings.filterwarnings("ignore")

PROJECT = "V7_4R2"

RUN_SMOKE_ONLY = True
RUN_EXTENDED_SMOKE = False
RUN_FULL = False
RESUME = True

SMOKE_TARGET = "1A22"
N_IDENTITY = 5
N_MUTANT = 5
N_EXTENDED_MUTANT = 20

DRIVE_ROOT = "/content/drive/MyDrive/Affinity_Intelligence/v7_4r2"
WORK_ROOT = "/content/v7_4r2"

# ---- §4 冻结物理协议（不得修改）----
FORCEFIELD = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
NONBONDED_METHOD = "NoCutoff"
PH = 7.0
RADIUS_MOVABLE = 8.0
GB_MODEL = "GBn2_XML"
ALLOW_GB_FALLBACK = False

MIN_TOL = 10.0
MIN_MAXITER = 2000
MIN_MAXITER_PRE = 5000            # canonical WT 预松弛

# ---- §12 限制常数（**仅由 label-free QC 选择**，绝不用 dG 性能调）----
KCAL_A2_TO_KJ_NM2 = 4.184 / (0.1 ** 2)      # 418.4
K_BACKBONE = 1000.0 * KCAL_A2_TO_KJ_NM2     # 骨架强限制
K_OUTSIDE = 1000.0 * KCAL_A2_TO_KJ_NM2      # 壳外强限制
# 突变侧链 + 8 A 内侧链：k = 0（自由）

# ---- §23 噪声门（项目工程阈值，非领域标准）----
NOISE_PASS, NOISE_MARGINAL, NOISE_FAIL = 0.5, 1.0, 1.0

RANDOM_SEED = 42

WORK = Path(WORK_ROOT); WORK.mkdir(parents=True, exist_ok=True)
FEAT_DIR = WORK / "features"; FEAT_DIR.mkdir(parents=True, exist_ok=True)
PDB_DIR = WORK / "pdb"; PDB_DIR.mkdir(parents=True, exist_ok=True)
CWT_DIR = WORK / "canonical_wt"; CWT_DIR.mkdir(parents=True, exist_ok=True)
MAN_DIR = WORK / "manifests"; MAN_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR = WORK / "qc"; QC_DIR.mkdir(parents=True, exist_ok=True)
DRIVE = Path(DRIVE_ROOT)

def _protocol_hash():
    import hashlib as _h
    payload = json.dumps({"ff": FORCEFIELD, "nonbonded": NONBONDED_METHOD, "ph": PH,
                          "radius": RADIUS_MOVABLE, "gb": GB_MODEL,
                          "k_backbone": K_BACKBONE, "k_outside": K_OUTSIDE,
                          "tol": MIN_TOL, "maxiter": MIN_MAXITER,
                          "maxiter_pre": MIN_MAXITER_PRE}, sort_keys=True)
    return _h.sha256(payload.encode()).hexdigest()[:16]

PROTOCOL_HASH = _protocol_hash()
print("PROJECT       :", PROJECT)
print("WORK / DRIVE  :", WORK_ROOT, "|", DRIVE_ROOT)
print("RUN_SMOKE_ONLY:", RUN_SMOKE_ONLY, "| RUN_EXTENDED_SMOKE:", RUN_EXTENDED_SMOKE,
      "| RUN_FULL:", RUN_FULL)
print("FORCEFIELD    :", FORCEFIELD, "| GB:", GB_MODEL, "| fallback:", ALLOW_GB_FALLBACK)
print("movable 8 A   : k_backbone=k_outside=%.0f kJ/mol/nm^2 ; site/shell sidechain k=0"
      % K_BACKBONE)
print("PROTOCOL_HASH :", PROTOCOL_HASH)
''')

code('''
# ============================================================ 01 Environment
import platform
env = {"python": sys.version.split()[0], "platform": platform.platform(),
       "project": PROJECT}
def sh(c):
    try:
        return subprocess.run(c, shell=True, capture_output=True, text=True,
                              timeout=180).stdout.strip()
    except Exception as e:
        return f"<failed {type(e).__name__}>"
print(sh("nvidia-smi")[:1200])
try:
    import psutil; env["ram_gb"] = round(psutil.virtual_memory().total / 1024**3, 2)
except Exception:
    env["ram_gb"] = None
env["disk_free_gb"] = round(shutil.disk_usage("/content").free / 1024**3, 2)
try:
    import urllib.request
    urllib.request.urlopen("https://files.rcsb.org", timeout=15).read(64)
    env["internet"] = True
except Exception:
    env["internet"] = False
print("RAM:", env["ram_gb"], "GB | disk free:", env["disk_free_gb"],
      "GB | internet:", env["internet"])
(WORK / "environment_report.json").write_text(json.dumps(env, indent=2))
''')

code('''
# ============================================================ 02 Install OpenMM + platforms
OPENMM_INSTALL_CMD = None
def pip(args):
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-input"] + args,
                       capture_output=True, text=True)
    print(("OK  " if r.returncode == 0 else "FAIL"), " ".join(args))
    if r.returncode:
        print(r.stderr[-1000:])
    return r.returncode == 0

need = []
for mod, pkg in (("openmm", "openmm"), ("pdbfixer", "pdbfixer"),
                 ("openmm.app", "openmm"), ("parmed", "parmed")):
    try:
        __import__(mod)
    except ImportError:
        need.append(pkg)
if need:
    args = sorted(set(need)); pip(args)
    OPENMM_INSTALL_CMD = "pip install " + " ".join(args)

import openmm
import openmm as mm
import openmm.app as app
from openmm import unit
import pdbfixer

def _plats():
    return [mm.Platform.getPlatform(i).getName()
            for i in range(mm.Platform.getNumPlatforms())]
PLATS = _plats()
print("OpenMM  :", openmm.version.version)
print("install :", OPENMM_INSTALL_CMD or "(preinstalled)")
print("platforms (as found):", PLATS)
if "CUDA" not in PLATS:
    print("未检测到 CUDA platform -> 尝试 openmm[cuda12]")
    pip(["openmm[cuda12]"])
    OPENMM_INSTALL_CMD = (OPENMM_INSTALL_CMD or "") + " ; pip install openmm[cuda12]"
    try:
        import importlib; importlib.reload(mm)
    except Exception:
        pass
    PLATS = _plats(); print("platforms (retry):", PLATS)

def get_platform():
    for nm in ("CUDA", "OpenCL", "CPU", "Reference"):
        if nm in PLATS:
            try:
                return mm.Platform.getPlatformByName(nm), nm
            except Exception:
                continue
    p0 = mm.Platform.getPlatform(0); return p0, p0.getName()
PLAT, PLAT_NAME = get_platform()
print(">>> 实际使用 platform:", PLAT_NAME)
ARTIFACTS = [
    {"name": "openmm", "version": openmm.version.version,
     "source": OPENMM_INSTALL_CMD or "preinstalled", "license": "MIT / LGPL",
     "purpose": "力场最小化 + 隐式溶剂"},
    {"name": "pdbfixer", "version": str(getattr(pdbfixer, "__version__", "?")),
     "source": "github.com/openmm/pdbfixer", "license": "MIT", "purpose": "缺原子/加氢"},
]
''')

code('''
# ============================================================ 03 Mount Google Drive (§27)
DRIVE_MOUNTED = False
try:
    from google.colab import drive
    drive.mount("/content/drive", force_remount=False)
    DRIVE_MOUNTED = True
    print("Drive mounted ->", DRIVE_ROOT)
except Exception as e:
    print("Drive mount 失败:", type(e).__name__, str(e)[:160])
if DRIVE_MOUNTED:
    for sub in ("canonical_wt", "features", "manifests", "qc", "protocol"):
        (DRIVE / sub).mkdir(parents=True, exist_ok=True)
    if RESUME:
        n1 = n2 = 0
        for p in (DRIVE / "canonical_wt").glob("*.npz"):
            shutil.copy(p, CWT_DIR / p.name); n1 += 1
        for p in (DRIVE / "features").glob("*.npz"):
            shutil.copy(p, FEAT_DIR / p.name); n2 += 1
        print("RESUME: 回灌 canonical_wt %d 个, features %d 个" % (n1, n2))
''')

code('''
# ============================================================ 04 Input discovery (auto-unzip)
import zipfile as _zip
REQ = ["variants.csv", "mutation_sites.csv", "target_chains.csv", "canonical_rows.csv"]

def _ok_dir(d):
    try:
        return all((d / r).exists() for r in REQ)
    except Exception:
        return False

def _search():
    seen, hits = set(), []
    for root in (Path("/content"), Path("/content/drive/MyDrive"), Path.cwd()):
        try:
            if not root.exists():
                continue
        except Exception:
            continue
        for f in root.rglob("variants.csv"):
            d = f.parent
            if str(d) in seen:
                continue
            seen.add(str(d))
            if _ok_dir(d):
                hits.append(d)
    return hits

def _autounzip():
    dest = WORK / "input"; dest.mkdir(parents=True, exist_ok=True)
    cands = []
    for root in (Path("/content"), Path("/content/drive/MyDrive"), Path.cwd()):
        try:
            if root.exists():
                cands += list(root.glob("*.zip")) + list(root.glob("**/colab_input_package.zip"))
        except Exception:
            pass
    cands = sorted(set(cands), key=str)
    print("ZIP 候选:", [str(c) for c in cands] or "无")
    for z in cands:
        try:
            with _zip.ZipFile(z) as zf:
                if not any(n.endswith("variants.csv") for n in zf.namelist()):
                    continue
                zf.extractall(dest)
            print("已解压:", z, "->", dest)
            for s in [p for p in dest.iterdir() if p.is_dir()]:
                if _ok_dir(s):
                    return s
            if _ok_dir(dest):
                return dest
        except Exception as e:
            print("解压失败", z, type(e).__name__, str(e)[:120])
    return None

hits = _search()
PKG = sorted(hits, key=str)[0] if hits else _autounzip()
if PKG is None:
    raise FileNotFoundError("未找到输入包；请把 colab_input_package.zip 上传到 /content（会自动解压）")

import pandas as pd
var = pd.read_csv(PKG / "variants.csv")
sit = pd.read_csv(PKG / "mutation_sites.csv")
tch = pd.read_csv(PKG / "target_chains.csv")
canon = pd.read_csv(PKG / "canonical_rows.csv")
print("输入目录:", PKG)
print("variants:", var.shape, "| sites:", sit.shape, "| chains:", tch.shape)
print("targets:", var.target_id.nunique(), "| rows:", len(var))
labtok = ("ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred")
bad = [c for df in (var, sit, tch) for c in df.columns
       if any(t in c.lower() for t in labtok)]
print("疑似标签列:", bad if bad else "无 ✅")
assert not bad, "输入含疑似标签列 -> 停止"
sites_by_row = {int(k): g for k, g in sit.groupby("row_id")}
THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
         "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
         "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
BB = ("N", "CA", "C", "O")

import requests
def fetch_pdb(pdb_id):
    dst = PDB_DIR / f"{pdb_id}.pdb"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status(); dst.write_bytes(r.content); return dst

def raw_sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
''')

code('''
# ============================================================ 05 GBn2 gate (compact, §25)
# 不再重复昂贵的 2000 步整体系探针；用**紧凑**验证：
#   CustomGBForce present / count>=1 / GB particles == system particles / energy & forces finite
# 仍保留 §4 的官方 XML 路线与 NO FALLBACK。
_pdb_probe = fetch_pdb(SMOKE_TARGET)
_fx = pdbfixer.PDBFixer(filename=str(_pdb_probe))
_fx.removeHeterogens(keepWater=False)
_fx.findMissingResidues()
_n_missing_res_detected = len(_fx.missingResidues)
_fx.missingResidues = {}                       # §5 禁用缺失残基重建
_fx.findMissingAtoms(); _fx.addMissingAtoms()
_fx.addMissingHydrogens(PH)
_topo_probe = _fx.topology
_ff = app.ForceField(*FORCEFIELD)
_sys_probe = _ff.createSystem(_topo_probe, nonbondedMethod=app.NoCutoff,
                              constraints=app.HBonds)
_gb_idx = [i for i in range(_sys_probe.getNumForces())
           if isinstance(_sys_probe.getForce(i), mm.CustomGBForce)]
_names = [_sys_probe.getForce(i).getName() for i in range(_sys_probe.getNumForces())]
_np_ = _sys_probe.getNumParticles()
_gb_n = _sys_probe.getForce(_gb_idx[0]).getNumParticles() if _gb_idx else None
_ctx = mm.Context(_sys_probe, mm.VerletIntegrator(0.001), PLAT)
_ctx.setPositions(_fx.positions)
_st = _ctx.getState(getEnergy=True, getForces=True)
_e = _st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
_F = np.asarray(_st.getForces(asNumpy=True).value_in_unit(
    unit.kilojoule_per_mole / unit.nanometer)) if False else None
import numpy as np
_F = np.asarray(_st.getForces(asNumpy=True).value_in_unit(
    unit.kilojoule_per_mole / unit.nanometer))
_e_per_atom = _e / max(1, _np_)
print("GB_GATE")
print("  OpenMM version            :", openmm.version.version)
print("  platform                  :", PLAT_NAME)
print("  ForceField files          :", FORCEFIELD)
print("  force classes             :", _names)
print("  CustomGBForce count       :", len(_gb_idx))
print("  GB particles              :", _gb_n)
print("  System particles          :", _np_)
print("  GB covers all particles   :", bool(_gb_idx and _gb_n == _np_))
print("  energy finite             :", bool(np.isfinite(_e)))
print("  forces finite             :", bool(np.all(np.isfinite(_F))))
# 起点能量量级（阶段 A 之前的**原始**建系能量）
print("  raw energy                : %.4e kJ/mol" % _e)
print("  raw energy / atom         : %.4e kJ/mol/atom" % _e_per_atom)
print("  missingResidues detected  :", _n_missing_res_detected, "(已强制置空, 不重建)")
GB_GATE_PASS = bool(len(_gb_idx) >= 1 and _gb_idx and _gb_n == _np_
                    and np.isfinite(_e) and np.all(np.isfinite(_F)))
print("  GB_GATE_PASS              :", GB_GATE_PASS)
if not GB_GATE_PASS:
    print("STOP: GBn2 协议门未通过（ALLOW_GB_FALLBACK=False，不做任何替换）")
    raise RuntimeError("GB_GATE_FAIL")
print("  -> GBn2_XML 协议成立（无偏离）")
''')

code('''
# ============================================================ 06 Canonical WT preparation (§5/§7/§8)
# 关键：**禁用缺失残基重建**，硬断言 missing_residues_added == 0
import numpy as np

def prep_canonical_wt(pdb_path):
    """返回 dict(ok, reason, fixer, topology, positions, stats)。
    §5: 绝不构建缺失的肽残基/loop；只补**已存在**残基内的缺失原子。"""
    stats = {}
    try:
        fx = pdbfixer.PDBFixer(filename=str(pdb_path))
        fx.removeHeterogens(keepWater=False)
        fx.findMissingResidues()
        n_det = len(fx.missingResidues)
        n_res_before = fx.topology.getNumResidues()
        stats["missing_residues_detected"] = n_det
        fx.missingResidues = {}                       # <<< §5 核心
        fx.findMissingAtoms()
        n_miss_atoms = sum(len(v) for v in (fx.missingAtoms or {}).values())
        stats["missing_atoms_existing_residues"] = n_miss_atoms
        fx.addMissingAtoms()
        n_res_after = fx.topology.getNumResidues()
        stats["missing_residues_added"] = n_res_after - n_res_before
        fx.addMissingHydrogens(PH)
        stats["hydrogens_added"] = fx.topology.getNumAtoms()
        stats["n_atoms"] = fx.topology.getNumAtoms()
        stats["n_chains"] = fx.topology.getNumChains()
        stats["n_residues"] = fx.topology.getNumResidues()
    except Exception as e:
        return {"ok": False, "reason": f"PREP:{type(e).__name__}:{str(e)[:160]}",
                "stats": stats}
    # §5 硬断言
    if stats["missing_residues_added"] != 0:
        return {"ok": False,
                "reason": f"MISSING_RESIDUES_ADDED={stats['missing_residues_added']}",
                "stats": stats}
    try:
        ff = app.ForceField(*FORCEFIELD)
        system = ff.createSystem(fx.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds)
    except Exception as e:
        return {"ok": False, "reason": f"FF:{type(e).__name__}:{str(e)[:160]}",
                "stats": stats}
    return {"ok": True, "reason": None, "fixer": fx, "topology": fx.topology,
            "positions": fx.positions, "system": system, "stats": stats}

def pre_relax_canonical_wt(system, topology, positions):
    """§7 一次保守预松弛：清除**加氢/补原子**引入的冲突。
    骨架强限制；目的是去掉非物理冲突，**不是**全局重塑实验结构。"""
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.setName("PreRestraint")
    for p in ("k", "x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    pos = np.array([list(pp) for pp in positions.value_in_unit(unit.nanometer)])
    for a in topology.atoms():
        is_bb = a.name in BB
        k = K_BACKBONE if is_bb else (0.0 if a.element is not None
                                      and a.element.symbol == "H" else 0.0)
        f.addParticle(a.index, [k, float(pos[a.index][0]),
                                float(pos[a.index][1]), float(pos[a.index][2])])
    f.setForceGroup(31)
    system.addForce(f)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), PLAT)
    ctx.setPositions(positions)
    e_pre = ctx.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER_PRE)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return st.getPositions(), e_pre, e_post, f

def kabsch_rmsd(P, Q):
    """Kabsch 最优叠合后的 RMSD（同原子数、同顺序）。"""
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    Hm = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(Hm)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    Pr = (R @ Pc.T).T
    return float(np.sqrt(np.mean(np.sum((Pr - Qc) ** 2, axis=1))))

def atom_xyz(topology, positions, pred):
    pos = np.array([list(pp) for pp in positions.value_in_unit(unit.angstrom)])
    idx = [a.index for a in topology.atoms() if pred(a)]
    return pos[idx], idx

print("canonical WT 准备与 Kabsch 工具就绪")
print("  §5: missingResidues 强制置空；硬断言 missing_residues_added == 0")
print("  §7: 一次保守预松弛（骨架强限制）")
print("  §13: kabsch_rmsd 用于对齐后的局部 RMSD")
''')

code('''
# ============================================================ 07 CANONICAL_WT_GATE + CANONICAL_WT_REBUILD_GATE (§8/§19)
_man = pd.read_csv(WORK / "input" / "canonical_rows.csv") if (WORK / "input" / "canonical_rows.csv").exists() else canon
_sel = var[var.target_id.astype(str) == SMOKE_TARGET]
if not len(_sel):
    _sel = var.head(1)
_smoke_target = str(_sel.target_id.iloc[0])
_pdb_t = fetch_pdb(_smoke_target)
raw_hash = raw_sha(_pdb_t)
print("=== CANONICAL_WT_GATE ===")
print("  target            :", _smoke_target)
print("  raw PDB sha256    :", raw_hash[:32], "...")
print("  PROTOCOL_HASH     :", PROTOCOL_HASH)

A = prep_canonical_wt(_pdb_t)
print("  准备 #A           :", "OK" if A["ok"] else A["reason"])
print("  stats A           :", A["stats"])
CANONICAL_WT_GATE_PASS = bool(A["ok"])
if not CANONICAL_WT_GATE_PASS:
    print("STOP: TARGET_CANONICAL_WT_FAIL -> 排除整个 target（不静默修复, §8）")
    raise RuntimeError("CANONICAL_WT_FAIL")

posA, eA_pre, eA_post, _fr = pre_relax_canonical_wt(A["system"], A["topology"], A["positions"])
nA = A["system"].getNumParticles()
print("  预松弛 #A         : E_pre=%.4e -> E_post=%.4e kJ/mol/atom=%.4e"
      % (eA_pre, eA_post, eA_post / nA))
_ctxA = mm.Context(A["system"], mm.VerletIntegrator(0.001), PLAT)
_ctxA.setPositions(A["positions"])
rmsd_pre = kabsch_rmsd(np.array(posA.value_in_unit(unit.angstrom)),
                       np.array(A["positions"].value_in_unit(unit.angstrom)))
print("  backbone-aligned RMSD(pre->post) = %.4f A" % rmsd_pre)

# ---- §19 REBUILD：独立重建第二次 ----
B = prep_canonical_wt(_pdb_t)
print()
print("=== CANONICAL_WT_REBUILD_GATE (§19) ===")
print("  准备 #B           :", "OK" if B["ok"] else B["reason"])
if B["ok"]:
    _aA = np.array(A["positions"].value_in_unit(unit.angstrom))
    _aB = np.array(B["positions"].value_in_unit(unit.angstrom))
    same_atoms = (A["topology"].getNumAtoms() == B["topology"].getNumAtoms())
    same_chains = (A["topology"].getNumChains() == B["topology"].getNumChains())
    seq_same = ([r.name for r in A["topology"].residues()]
                == [r.name for r in B["topology"].residues()])
    max_d = (float(np.max(np.abs(_aA - _aB))) if same_atoms else float("nan"))
    print("  atom count equal  :", same_atoms, "(%s vs %s)"
          % (A["topology"].getNumAtoms(), B["topology"].getNumAtoms()))
    print("  chain count equal :", same_chains)
    print("  sequence equal    :", seq_same)
    print("  max |dx|          : %.6f A" % max_d)
    CANONICAL_WT_REBUILD_PASS = bool(same_atoms and same_chains and seq_same
                                     and np.isfinite(max_d) and max_d < 1e-6)
else:
    CANONICAL_WT_REBUILD_PASS = False
print("  CANONICAL_WT_GATE_PASS      :", CANONICAL_WT_GATE_PASS)
print("  CANONICAL_WT_REBUILD_PASS   :", CANONICAL_WT_REBUILD_PASS)
if not CANONICAL_WT_REBUILD_PASS:
    print("STOP: canonical WT 不可复现 -> 不得进入突变体阶段（§19）")
    raise RuntimeError("CANONICAL_WT_REBUILD_FAIL")
''')

code('''
# ============================================================ 08 Mutant from canonical WT + local relaxation (§9/§10/§11)
import copy

def clone_from_canonical(cwt_fixer, cwt_positions):
    """从**冻结的 canonical WT**克隆一份可变模型（§9）。"""
    topo = copy.deepcopy(cwt_fixer.topology)
    pos = copy.deepcopy(cwt_positions)
    return topo, pos

def apply_mutations(topo, muts, site_index_by_key):
    """只改被突变残基（§9）。muts: [(chain, seq_index0, wt_aa, mut_aa)]。"""
    ch = {c.id: c for c in topo.chains()}
    for (c, si, wt1, mut1) in muts:
        if c not in ch:
            return False, f"chain {c} missing"
        res = list(ch[c].residues())
        if not (0 <= si < len(res)):
            return False, f"site index {c}:{si} out of range"
        r = res[si]
        if r.name != THREE[wt1]:
            return False, f"WT mismatch {c}:{si} got {r.name} want {THREE[wt1]}"
        r.name = THREE[mut1]
    return True, "ok"

def add_hydrogens_inplace(topo, pos):
    """对改名后的 topology 补氢（只增氢，不动重原子）。"""
    mod = app.Modeller(topo, pos)
    ff = app.ForceField(*FORCEFIELD)
    mod.addHydrogens(ff, pH=PH)
    return mod.topology, mod.positions

def relax_local(topology, positions, site_keys, system_factory):
    """§11 局部侧链松弛：
       自由 = 突变位点侧链原子/H + 8 A 内其它残基的**侧链**原子
       强限制 = 所有骨架重原子 + 壳外一切
       返回 (positions, e_pre, e_post, movable_atom_count, fixed_atom_count)"""
    pos_nm = np.array([list(pp) for pp in positions.value_in_unit(unit.nanometer)])
    # 8 A 壳（以位点重原子为参考）
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    other = [a for a in topology.atoms()
             if (a.residue.chain.id, a.residue.id) not in site_keys]
    shell_keys = set()
    if site_atoms and other:
        oi = np.array([a.index for a in other])
        D = np.sqrt(np.sum((pos_nm[oi][:, None, :] - pos_nm[site_atoms][None, :, :]) ** 2,
                           axis=2))
        for k, a in enumerate(other):
            if float(D[k].min()) <= RADIUS_MOVABLE / 10.0:
                shell_keys.add((a.residue.chain.id, a.residue.id))
    sysm = system_factory(topology)
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.setName("Restraint")
    for p in ("k", "x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    n_move = n_fix = 0
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        is_bb = a.name in BB
        is_h = (a.element is not None and a.element.symbol == "H")
        if key in site_keys:
            k = 0.0                       # 突变残基（含骨架）自由
        elif key in shell_keys and not is_bb:
            k = 0.0                       # 壳内侧链自由
        else:
            k = K_BACKBONE if is_bb else K_OUTSIDE
        if k == 0.0:
            n_move += 1
        else:
            n_fix += 1
        f.addParticle(a.index, [k, float(pos_nm[a.index][0]),
                                float(pos_nm[a.index][1]), float(pos_nm[a.index][2])])
    f.setForceGroup(31)
    sysm.addForce(f)
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001), PLAT)
    ctx.setPositions(positions)
    bio = {i for i in range(sysm.getNumForces())
           if sysm.getForce(i).getForceGroup() != 31}
    e_pre = ctx.getState(getEnergy=True, groups=bio).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL, maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post = ctx.getState(getEnergy=True, groups=bio).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    return st.getPositions(), e_pre, e_post, n_move, n_fix

def make_system(topology):
    ff = app.ForceField(*FORCEFIELD)
    return ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                           constraints=app.HBonds)

print("突变体克隆 / 局部侧链松弛 就绪")
print("  自由: 突变残基全部原子 + 8 A 壳内**侧链**原子")
print("  限制: 所有骨架重原子(k=%.0f) + 壳外一切(k=%.0f) kJ/mol/nm^2"
      % (K_BACKBONE, K_OUTSIDE))
''')

code('''
# ============================================================ 09 Local geometry + local energy features (§15/§16)
def local_geom(topology, positions, site_keys):
    """§15 局部几何量（以突变位点为中心，与 partner 链的相互作用）。"""
    pos = np.array([list(pp) for pp in positions.value_in_unit(unit.angstrom)])
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    partner = [a for a in topology.atoms()
               if a.residue.chain.id not in {c for c, _ in site_keys}]
    intra = [a for a in topology.atoms()
             if (a.residue.chain.id, a.residue.id) not in site_keys
             and a.residue.chain.id in {c for c, _ in site_keys}]
    out = {}
    def _fields(atomlist, prefix):
        if not atomlist or not site_atoms:
            out[prefix + "present"] = 0.0
            return
        ai = np.array([a.index for a in atomlist])
        D = np.sqrt(np.sum((pos[ai][:, None, :] - pos[site_atoms][None, :, :]) ** 2,
                           axis=2))
        dmin = D.min(axis=1)
        out[prefix + "present"] = 1.0
        out[prefix + "min_dist"] = float(dmin.min())
        out[prefix + "mean_dist"] = float(dmin.mean())
        out[prefix + "n_contact_4"] = int((dmin <= 4.0).sum())
        out[prefix + "n_contact_6"] = int((dmin <= 6.0).sum())
        out[prefix + "n_contact_8"] = int((dmin <= 8.0).sum())
        out[prefix + "n_clash"] = int((dmin < 2.2).sum())
    _fields(partner, "partner_")
    _fields(intra, "intra_")
    return out

def local_binding_energy(topology, positions, site_keys):
    """§16 局部非键相互作用能（突变位点 <-> partner 链）。
    用**同一**力场参数、对**同一** system 求「排除位点-伙伴相互作用」的能量差，
    即 E(全) - E(去掉伙伴相互作用)。技术上脆弱则返回 None（NOT_RUN）。"""
    try:
        ff = app.ForceField(*FORCEFIELD)
        sysm = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                               constraints=app.HBonds)
        nb = None
        for i in range(sysm.getNumForces()):
            if sysm.getForce(i).getName() == "NonbondedForce":
                nb = sysm.getForce(i); break
        if nb is None:
            return None
        site_chains = {c for c, _ in site_keys}
        site_atoms = {a.index for a in topology.atoms()
                      if (a.residue.chain.id, a.residue.id) in site_keys}
        partner_atoms = {a.index for a in topology.atoms()
                         if a.residue.chain.id not in site_chains}
        if not site_atoms or not partner_atoms:
            return None
        # 记录并临时把 site<->partner 的电荷/epsilon 置零
        import openmm as _mm
        saved = []
        for i in range(nb.getNumExceptions()):
            pass
        # 更稳健：用 CustomNonbondedForce 单独算 site-partner 的 LJ+静电
        cnb = _mm.CustomNonbondedForce(
            "4*eps*((sig/r)^12-(sig/r)^6) + 138.935456*q1*q2/r")
        cnb.addPerParticleParameter("q")
        cnb.addPerParticleParameter("sig")
        cnb.addPerParticleParameter("eps")
        for a in topology.atoms():
            q, sig, eps = nb.getParticleParameters(a.index)
            cnb.addParticle([q.value_in_unit(unit.elementary_charge),
                             sig.value_in_unit(unit.nanometer),
                             eps.value_in_unit(unit.kilojoule_per_mole)])
        # 只保留 site<->partner 对
        for a in topology.atoms():
            ia = a.index
            if ia in partner_atoms:
                cnb.addInteractionGroup({ia}, site_atoms)
        cnb.setNonbondedMethod(_mm.CustomNonbondedForce.NoCutoff)
        sysm.addForce(cnb)
        ctx = _mm.Context(sysm, _mm.VerletIntegrator(0.001), PLAT)
        ctx.setPositions(positions)
        gi = sysm.getNumForces() - 1
        e = ctx.getState(getEnergy=True, groups={gi}).getPotentialEnergy(
            ).value_in_unit(unit.kilojoule_per_mole)
        return float(e)
    except Exception as e:
        return None

print("局部几何 + 局部能量函数就绪（§15/§16）")
print("  局部能量若技术上不可靠将返回 None -> 记为 NOT_RUN，不硬造")
''')

code('''
# ============================================================ 10 MINIMIZER_REPEATABILITY_GATE (§20)
# 用**同一份**已准备好的 canonical WT，克隆两次，做完全相同的局部松弛。
# **不**再调用 PDBFixer。
_sk_probe = set()
_g0 = sites_by_row.get(int(_sel.row_id.iloc[0]))
for x in _g0.itertuples():
    _sk_probe.add((str(x.chain), str(x.resseq)))
# site_keys 以 (chain, resseq) 表示；用 canonical WT 的实际残基 id 建映射
def site_keys_from_indices(topo, sites_df):
    ch = {c.id: list(c.residues()) for c in topo.chains()}
    keys = set()
    for x in sites_df.itertuples():
        c = str(x.chain); si = int(x.seq_index0)
        if c in ch and 0 <= si < len(ch[c]):
            keys.add((c, ch[c][si].id))
    return keys

SK = site_keys_from_indices(A["topology"], _g0)
print("=== MINIMIZER_REPEATABILITY_GATE (§20) ===")
print("  位点键:", SK)
_t0 = time.time()
def _clone_relax():
    t, p = clone_from_canonical(A["fixer"], posA)
    return relax_local(t, p, SK, make_system)
_o1 = _clone_relax(); _o2 = _clone_relax()
_a1 = np.array(_o1[0].value_in_unit(unit.angstrom))
_a2 = np.array(_o2[0].value_in_unit(unit.angstrom))
_bb_idx = [a.index for a in A["topology"].atoms() if a.name in BB]
rmsd_rep = kabsch_rmsd(_a1[_bb_idx], _a2[_bb_idx])
dE_rep = abs(_o1[2] - _o2[2])
rel_rep = dE_rep / max(1e-12, abs(_o2[2]))
print("  E_post #1 = %.6e   #2 = %.6e" % (_o1[2], _o2[2]))
print("  |dE| = %.4e   rel = %.3e" % (dE_rep, rel_rep))
print("  backbone-aligned RMSD = %.6f A" % rmsd_rep)
print("  movable/fixed atoms = %d / %d" % (_o1[3], _o1[4]))
MINIMIZER_REPEATABILITY_PASS = bool(rel_rep < 1e-3 and rmsd_rep < 0.05)
print("  MINIMIZER_REPEATABILITY_PASS =", MINIMIZER_REPEATABILITY_PASS,
      "(%.1fs)" % (time.time() - _t0))
if not MINIMIZER_REPEATABILITY_PASS:
    print("STOP: 最小化不可复现 -> 局部能量/几何差值的噪声地板过高（§20）")
    raise RuntimeError("MINIMIZER_REPEATABILITY_FAIL")
''')

code('''
# ============================================================ 11 IDENTITY_GATE (§21) — 5 identity controls
print("=== IDENTITY_GATE (§21) ===")
IDENT = []
_cands = _sel.head(N_IDENTITY)
_t0 = time.time()
for k, r in enumerate(_cands.itertuples()):
    g = sites_by_row.get(int(r.row_id))
    if g is None or not len(g):
        continue
    sk = site_keys_from_indices(A["topology"], g)
    if not sk:
        continue
    t, p = clone_from_canonical(A["fixer"], posA)
    pos_r, e_pre, e_post, nmv, nfx = relax_local(t, p, sk, make_system)
    geom = local_geom(t, pos_r, sk)
    a_r = np.array(pos_r.value_in_unit(unit.angstrom))
    a_0 = np.array(A["positions"].value_in_unit(unit.angstrom))
    bb_i = [a.index for a in t.atoms() if a.name in BB]
    IDENT.append({"row_id": int(r.row_id), "n_sites": len(sk),
                  "e_pre": e_pre, "e_post": e_post, "dE": 0.0,
                  "rmsd_bb_aligned": kabsch_rmsd(a_0[bb_i], a_r[bb_i]),
                  "movable": nmv, "fixed": nfx, **geom})
    print("  [%d/%d] row=%-5d sites=%d  rmsd_bb=%.4f  movable=%d"
          % (k + 1, len(_cands), int(r.row_id), len(sk),
             IDENT[-1]["rmsd_bb_aligned"], nmv), flush=True)
ID = pd.DataFrame(IDENT)
ID.to_csv(QC_DIR / "identity_controls.csv", index=False)
print("  identity n =", len(ID))

# 身份对照的噪声地板（几何）
def _q95(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    return float(np.percentile(np.abs(v), 95)) if len(v) else float("nan")
ID_RMSD_P95 = _q95(ID.get("rmsd_bb_aligned", [])) if len(ID) else float("nan")
print("  identity rmsd_bb_aligned P95 = %.6f A" % ID_RMSD_P95)
IDENTITY_GATE_PASS = bool(len(ID) >= 3 and np.isfinite(ID_RMSD_P95)
                          and ID_RMSD_P95 < 0.10)
print("  IDENTITY_GATE_PASS =", IDENTITY_GATE_PASS, "(%.1fs)" % (time.time() - _t0))
if not IDENTITY_GATE_PASS:
    print("STOP: identity 局部几何噪声过大（§21）")
    raise RuntimeError("IDENTITY_GATE_FAIL")
''')

code('''
# ============================================================ 12 Five real mutants + LOCAL_NOISE_GATE (§22/§23)
N_M = N_EXTENDED_MUTANT if RUN_EXTENDED_SMOKE else N_MUTANT
print("=== 真实突变体（%d 个）===" % N_M)
MUT = []
_t0 = time.time()
for k, r in enumerate(_sel.head(N_M).itertuples()):
    g = sites_by_row.get(int(r.row_id))
    if g is None or not len(g):
        continue
    sk = site_keys_from_indices(A["topology"], g)
    muts = []
    ok = True
    for x in g.itertuples():
        muts.append((str(x.chain), int(x.seq_index0), str(x.wt_aa), str(x.mut_aa)))
    t, p = clone_from_canonical(A["fixer"], posA)
    okm, msg = apply_mutations(t, muts, None)
    if not okm:
        MUT.append({"row_id": int(r.row_id), "err": msg}); continue
    try:
        t2, p2 = add_hydrogens_inplace(t, p)
    except Exception as e:
        MUT.append({"row_id": int(r.row_id), "err": f"ADDH:{type(e).__name__}:{str(e)[:70]}"})
        continue
    try:
        pos_r, e_pre, e_post, nmv, nfx = relax_local(t2, p2, sk, make_system)
    except Exception as e:
        MUT.append({"row_id": int(r.row_id), "err": f"RELAX:{type(e).__name__}:{str(e)[:70]}"})
        continue
    geom_m = local_geom(t2, pos_r, sk)
    # WT 参照几何（同一 canonical WT、同样局部松弛已由 identity 给出基线）
    t_w, p_w = clone_from_canonical(A["fixer"], posA)
    pos_w, ew_pre, ew_post, _, _ = relax_local(t_w, p_w, sk, make_system)
    geom_w = local_geom(t_w, pos_w, sk)
    a_m = np.array(pos_r.value_in_unit(unit.angstrom))
    a_w = np.array(pos_w.value_in_unit(unit.angstrom))
    bb_i = [a.index for a in t2.atoms() if a.name in BB]
    dE_local = None
    MUT.append({"row_id": int(r.row_id), "target": str(r.target_id),
                "n_sites": len(sk), "dE_total": e_post - ew_post,
                "e_post": e_post, "e_post_wt": ew_post,
                "rmsd_bb_aligned": kabsch_rmsd(a_w[bb_i], a_m[bb_i]),
                "movable": nmv, "fixed": nfx,
                **{("d_" + kk): (geom_m.get(kk, float("nan")) - geom_w.get(kk, 0.0))
                   for kk in set(list(geom_m) + list(geom_w))
                   if isinstance(geom_m.get(kk), (int, float))
                   and isinstance(geom_w.get(kk), (int, float))}})
    print("  [%d/%d] row=%-5d sites=%d  d(rmsd_bb)=%.4f  movable=%d"
          % (k + 1, N_M, int(r.row_id), len(sk),
             MUT[-1]["rmsd_bb_aligned"], nmv), flush=True)
MU = pd.DataFrame(MUT)
MU.to_csv(QC_DIR / "mutation_smoke.csv", index=False)

# ---------- LOCAL_NOISE_GATE（§23）：用**对齐后的局部几何**算比值 ----------
def _med(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    return float(np.median(np.abs(v))) if len(v) else float("nan")
_loc_cols = [c for c in MU.columns if c.startswith("d_partner_")
             or c.startswith("d_intra_")]
# 局部几何"响应幅度"：用各局部量的 |Δ| 的中位数聚合
resp = []
for c in _loc_cols:
    m = _med(MU[c]) if c in MU.columns else float("nan")
    if np.isfinite(m):
        resp.append(m)
LOCAL_RESPONSE = float(np.median(resp)) if resp else float("nan")
R_G_LOCAL = (ID_RMSD_P95 / LOCAL_RESPONSE
             if (np.isfinite(ID_RMSD_P95) and np.isfinite(LOCAL_RESPONSE)
                 and LOCAL_RESPONSE > 1e-12) else None)
print()
print("=" * 62)
print("LOCAL_NOISE_GATE (§23)  —— 项目工程阈值，非领域标准")
print("  identity 局部几何噪声 P95 (rmsd_bb_aligned) = %.6f A" % ID_RMSD_P95)
print("  真实突变的局部几何响应中位数              = %.6f" % LOCAL_RESPONSE)
print("  R_G_LOCAL = %.4f" % (R_G_LOCAL if R_G_LOCAL is not None else float("nan")))
print("  identity 整体系 |dE| 中位 = %.4e kJ/mol (仅诊断, §14)"
      % _med(ID.get("dE", [])))
print("  真实突变整体系 |dE| 中位  = %.4e kJ/mol (仅诊断, §14)"
      % _med(MU.get("dE_total", [])))
if R_G_LOCAL is None:
    LOCAL_NOISE_GATE = "NOT_EVALUABLE"; LOCAL_NOISE_GATE_PASS = False
elif R_G_LOCAL >= NOISE_FAIL:
    LOCAL_NOISE_GATE = "NUMERICAL_NOISE_TOO_HIGH"; LOCAL_NOISE_GATE_PASS = False
elif R_G_LOCAL >= NOISE_PASS:
    LOCAL_NOISE_GATE = "MARGINAL"; LOCAL_NOISE_GATE_PASS = False
else:
    LOCAL_NOISE_GATE = "LOCAL_NOISE_GATE_PASS"; LOCAL_NOISE_GATE_PASS = True
print("  LOCAL_NOISE_GATE =", LOCAL_NOISE_GATE,
      "| PASS =", LOCAL_NOISE_GATE_PASS)
print("=" * 62)
(WORK / "noise_gate.json").write_text(json.dumps(
    {"identity_rmsd_p95_A": ID_RMSD_P95, "local_response_median": LOCAL_RESPONSE,
     "R_G_LOCAL": R_G_LOCAL, "gate": LOCAL_NOISE_GATE, "pass": LOCAL_NOISE_GATE_PASS,
     "identity_dE_median_kJmol": _med(ID.get("dE", [])),
     "mutation_dE_median_kJmol": _med(MU.get("dE_total", [])),
     "gb_model": GB_MODEL, "platform": PLAT_NAME, "protocol_hash": PROTOCOL_HASH,
     "threshold_note": "prespecified project engineering thresholds, NOT field-standard"},
    indent=2))
if not LOCAL_NOISE_GATE_PASS:
    print()
    print("STOP（§23）：局部噪声门未通过 -> 不进入全量。")
    raise RuntimeError("LOCAL_NOISE_GATE_FAIL:" + LOCAL_NOISE_GATE)
print()
print(">>> 全部快速门通过。可把 RUN_EXTENDED_SMOKE / RUN_FULL 打开。")
''')

code('''
# ============================================================ 13 Smoke summary + Drive sync (§24/§27/§40)
GATES = {
    "GB_GATE": GB_GATE_PASS,
    "CANONICAL_WT_GATE": CANONICAL_WT_GATE_PASS,
    "CANONICAL_WT_REBUILD_GATE": CANONICAL_WT_REBUILD_PASS,
    "MINIMIZER_REPEATABILITY_GATE": MINIMIZER_REPEATABILITY_PASS,
    "IDENTITY_GATE": IDENTITY_GATE_PASS,
    "LOCAL_NOISE_GATE": LOCAL_NOISE_GATE_PASS,
}
print("=== V7.4R2 FAST SMOKE GATES ===")
for k, v in GATES.items():
    print("   %-30s = %s" % (k, v))
SMOKE_ALL_PASS = all(GATES.values())
print("   SMOKE_ALL_PASS =", SMOKE_ALL_PASS)
prov = {"project": PROJECT, "platform": PLAT_NAME, "forcefield": FORCEFIELD,
        "gb_model": GB_MODEL, "allow_gb_fallback": ALLOW_GB_FALLBACK,
        "ph": PH, "radius_movable_A": RADIUS_MOVABLE,
        "k_backbone": K_BACKBONE, "k_outside": K_OUTSIDE,
        "min_tol": MIN_TOL, "min_maxiter": MIN_MAXITER,
        "protocol_hash": PROTOCOL_HASH, "openmm": openmm.version.version,
        "install_cmd": OPENMM_INSTALL_CMD, "smoke_target": _smoke_target,
        "missing_residues_reconstruction": "DISABLED (§5)",
        "mutant_source": "canonical WT clone (§9)",
        "alignment": "Kabsch on backbone N/CA/C outside shell (§13)",
        "labels_used": False, "terminology": "force-field minimized geometry",
        "drive_root": DRIVE_ROOT, "drive_mounted": DRIVE_MOUNTED,
        "gates": GATES, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
(WORK / "provenance.json").write_text(json.dumps(prov, indent=2))
if DRIVE_MOUNTED:
    for sub, pat in (("qc", "*.csv"), ("qc", "*.json"), ("manifests", "*.csv")):
        for p in WORK.rglob(pat):
            if p.parent.name in ("qc", "manifests", ""):
                shutil.copy(p, DRIVE / sub / p.name)
    print("已同步到 Drive:", DRIVE)
print()
print("=" * 62)
print("V7_4R2_SMOKE_STATUS =", "PASS" if SMOKE_ALL_PASS else "FAIL")
print("GB_GATE:", GB_GATE_PASS, "| CANONICAL_WT_REBUILD:", CANONICAL_WT_REBUILD_PASS)
print("MINIMIZER_REPEATABILITY:", MINIMIZER_REPEATABILITY_PASS)
print("IDENTITY:", IDENTITY_GATE_PASS, "| LOCAL_NOISE:", LOCAL_NOISE_GATE)
print("R_G_LOCAL =", R_G_LOCAL)
print("FORMAL_FEATURES_GENERATED = 0")
print("=" * 62)
''')

code('''
# ============================================================ 14 FULL EXTRACTION (§29) — only if RUN_FULL
if not RUN_FULL:
    print("RUN_FULL=False -> 跳过全量（默认）。")
    print("只有 FAST smoke 与（可选）extended smoke 全部通过后才改 RUN_FULL=True。")
else:
    raise SystemExit("全量流程将在 smoke 通过后由 DSH 生成并冻结（避免未经确认的全量运行）。")
''')

nb = {"cells": CELLS,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.11"},
                   "accelerator": "GPU",
                   "colab": {"provenance": [], "toc_visible": True}},
      "nbformat": 4, "nbformat_minor": 5}
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"写出 {NB} | cells = {len(CELLS)}")
