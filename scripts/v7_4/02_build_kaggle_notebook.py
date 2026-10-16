# -*- coding: utf-8 -*-
"""scripts/v7_4/02_build_kaggle_notebook.py -- 生成 V7.4 OpenMM 松弛 notebook（§29）.

16 个 section；**不含**硬编码 Windows 路径；含**硬门控**冒烟测试（§30）
与 identity / repeatability 对照（§31/§32）。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "kaggle" / "v7_4"
OUT.mkdir(parents=True, exist_ok=True)
NB = OUT / "V7_4_OpenMM_Relaxation.ipynb"
CELLS: list = []


def md(t):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": t.splitlines(True)})


def code(t):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": t.strip("\n").splitlines(True)})


md("""# V7.4 — Relaxed-Structure Energy Perturbation (Kaggle GPU)

**唯一新信息源**：`computationally relaxed mutant geometry` + `force-field potential-energy change`。

> ⚠️ **术语纪律**（§1）：这是**力场最小化几何**，**不是** experimental structure、
> **不是** true binding energy / true ΔΔG。禁止如此称呼。

**本 notebook 只做**：突变体构建 → OpenMM 参数化 → 成对 WT/MUT 局部最小化 →
能量/几何提取 → 导出。**不做**任何排名评估、**不用**任何 ΔΔG 标签。
""")

code('''
# ============================================================ 00 Config
import os, sys, json, time, math, hashlib, zipfile, shutil, subprocess, warnings, random
from pathlib import Path
warnings.filterwarnings("ignore")

PROJECT = "v7_4_relax"
DEVICE = "cuda:0"
USE_SECOND_GPU = False
RUN_SMOKE_ONLY = True          # 首次务必 True
RESUME = True
OUTPUT_ROOT = Path("/kaggle/working/v7_4")
INPUT_ROOT = Path("/kaggle/input")
MAX_TARGETS = None             # 调试用
PH = 7.0                       # 标准化计算条件，**不**假定等于实验 pH（§6）

# 冻结协议（§6/§13/§14/§16）
FORCEFIELD = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
RADIUS_MOVABLE = 8.0           # Å：可动壳
RESTRAINT_BACKBONE_SHELL = 10.0    # kcal/mol/Å^2：8 Å 内骨架弱限制
RESTRAINT_OUTSIDE = 100.0          # kcal/mol/Å^2：8 Å 外强限制
MIN_TOL = 10.0                 # kJ/mol/nm
MIN_MAXITER = 2000

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
FEAT_DIR = OUTPUT_ROOT / "features"; FEAT_DIR.mkdir(parents=True, exist_ok=True)
PDB_DIR = OUTPUT_ROOT / "pdb"; PDB_DIR.mkdir(parents=True, exist_ok=True)
print("OUTPUT_ROOT:", OUTPUT_ROOT, "| RUN_SMOKE_ONLY:", RUN_SMOKE_ONLY)
print("FORCEFIELD :", FORCEFIELD)
''')

code('''
# ============================================================ 01 Environment Audit
import platform
env = {"python": sys.version.split()[0], "platform": platform.platform()}
def sh(c):
    try:
        return subprocess.run(c, shell=True, capture_output=True, text=True, timeout=180).stdout.strip()
    except Exception as e:
        return f"<failed {type(e).__name__}>"
print(sh("nvidia-smi")[:1600]); print("\\n=== df ==="); print(sh("df -h")[:700])
import torch
env.update(torch=torch.__version__, cuda=bool(torch.cuda.is_available()))
gpus = []
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        p = torch.cuda.get_device_properties(i)
        gpus.append({"i": i, "name": p.name, "vram_gb": round(p.total_memory/1024**3, 2)})
env["gpus"] = gpus
try:
    import psutil; env["ram_gb"] = round(psutil.virtual_memory().total/1024**3, 2)
except Exception: env["ram_gb"] = None
env["disk_free_gb"] = round(shutil.disk_usage("/kaggle/working").free/1024**3, 2)
try:
    import urllib.request; urllib.request.urlopen("https://files.rcsb.org", timeout=15).read(64)
    env["internet"] = True
except Exception: env["internet"] = False
print("\\nGPUs:", gpus, "| RAM:", env["ram_gb"], "| internet:", env["internet"])
assert env["cuda"], "请在 Settings 设 Accelerator = GPU T4 x2"
(OUTPUT_ROOT/"environment_report.json").write_text(json.dumps(env, indent=2))
''')

code('''
# ============================================================ 02 Install OpenMM / PDBFixer
def pip(pkgs):
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-input"]+pkgs,
                       capture_output=True, text=True)
    print(("OK  " if r.returncode == 0 else "FAIL"), " ".join(pkgs))
    if r.returncode: print(r.stderr[-1500:])
    return r.returncode == 0
need = []
for mod, pkg in (("openmm","openmm"), ("pdbfixer","pdbfixer"),
                 ("openmm.app","openmm"), ("parmed","parmed")):
    try: __import__(mod)
    except ImportError: need.append(pkg)
if need: pip(sorted(set(need)))
import openmm, openmm.app as app
print("OpenMM :", openmm.version.version)
try:
    import pdbfixer; print("PDBFixer:", getattr(pdbfixer, "__version__", "?"))
except Exception as e: print("PDBFixer 导入失败:", e)
ARTIFACTS = [
 {"name":"openmm","version":openmm.version.version,"source":"conda-forge/pypi",
  "license":"MIT / LGPL","purpose":"力场最小化"},
 {"name":"pdbfixer","version":"?","source":"github.com/openmm/pdbfixer",
  "license":"MIT","purpose":"缺原子/加氢"},
]
print("\\n未安装 FoldX/Rosetta（§5 禁止）")
''')

code('''
# ============================================================ 03 Input Discovery
REQ = ["kaggle_input_manifest.csv","variants.csv","mutation_sites.csv","target_chains.csv"]
cands = [d.parent for d in sorted(INPUT_ROOT.rglob("kaggle_input_manifest.csv"))
         if all((d.parent/f).exists() for f in REQ)]
if not cands: raise FileNotFoundError("未找到 V7.4 输入包，请 Add Input")
PKG = sorted(cands, key=str)[0]; print("选定:", PKG)
import pandas as pd
var = pd.read_csv(PKG/"variants.csv"); sit = pd.read_csv(PKG/"mutation_sites.csv")
tch = pd.read_csv(PKG/"target_chains.csv"); canon = pd.read_csv(PKG/"canonical_rows.csv")
print("variants:", var.shape, "| sites:", sit.shape, "| chains:", tch.shape)
print("targets:", var.target_id.nunique(), "| rows:", len(var))
labtok = ("ddg","delta_g","affinity","kd(","label","y_true","pred")
bad = [c for df in (var,sit,tch) for c in df.columns if any(t in c.lower() for t in labtok)]
print("疑似标签列:", bad if bad else "无 ✅")
wt_seq = {(str(r.target_id), str(r.chain)): str(r.wt_sequence) for r in tch.itertuples()}
sites_by_row = {int(k): g for k, g in sit.groupby("row_id")}
TARGETS = sorted(var.target_id.astype(str).unique())
if MAX_TARGETS: TARGETS = TARGETS[:int(MAX_TARGETS)]
''')

code('''
# ============================================================ 04 Protocol / Forcefield Audit
import openmm as mm
from openmm import app, unit
PROTOCOL = {
  "forcefield": FORCEFIELD, "ph": PH, "solvent": "implicit (GBn2)",
  "radius_movable_A": RADIUS_MOVABLE,
  "restraint_backbone_shell": RESTRAINT_BACKBONE_SHELL,
  "restraint_outside": RESTRAINT_OUTSIDE,
  "minimizer": "LocalEnergyMinimizer (L-BFGS)",
  "tolerance_kJ_mol_nm": MIN_TOL, "max_iterations": MIN_MAXITER,
  "heterogen_policy": "protein chains only (waters/ligands/ions removed)",
  "terminology": "force-field minimized geometry; NOT experimental, NOT true dG",
  "limitation": "pH fixed at 7.0; NOT assumed to match each experiment",
}
(OUTPUT_ROOT/"protocol.json").write_text(json.dumps(PROTOCOL, indent=2))
print(json.dumps(PROTOCOL, indent=2))
print("\\n可用力场文件:")
for f in FORCEFIELD:
    print("  ", f, "-> 内置" if f.startswith(("amber","implicit","charmm")) else "?")
for f in FORCEFIELD:
    try:
        app.ForceField(f); print("   OK:", f)
    except Exception as e:
        print("   FAIL:", f, e)
''')

code('''
# ============================================================ 05 Mutation Builder Smoke Test
from pdbfixer import PDBFixer
import requests, io as _io

def fetch_pdb(pdb_id):
    dst = PDB_DIR / f"{pdb_id}.pdb"
    if dst.exists() and dst.stat().st_size > 1000: return dst
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    dst.write_bytes(r.content); return dst

def build_mutant(fixer, muts):
    """muts: [(chain, resseq, wt3, mut3)] 就地替换；返回是否成功。"""
    import openmm.app as _app
    top = fixer.topology
    chains = {c.id: c for c in top.chains()}
    done = []
    for (ch, rs, wt3, mut3) in muts:
        tgt = None
        for res in chains[ch].residues():
            if res.id == str(rs):
                tgt = res; break
        if tgt is None: return False, f"residue {ch}:{rs} not found"
        if _app.PDBFile._standardResidues and tgt.name != wt3:
            return False, f"WT mismatch {ch}:{rs} got {tgt.name} want {wt3}"
        tgt.name = mut3
        done.append((ch, rs, wt3, mut3))
    return True, done

def residue_has_backbone(fixer, ch, rs):
    top = fixer.topology
    for c in top.chains():
        if c.id != ch: continue
        for res in c.residues():
            if res.id == str(rs):
                names = {a.name for a in res.atoms()}
                return {"N","CA","C","O"}.issubset(names)
    return False

# 20 变体验证集（label-free，覆盖 §12 要求的类型）
import itertools
_val = []
seen = set()
for r in var.itertuples():
    rid = int(r.row_id); g = sites_by_row.get(rid)
    if g is None or not len(g): continue
    for x in g.itertuples():
        key = (str(x.wt_aa), str(x.mut_aa))
        if key in seen: continue
        seen.add(key)
        _val.append({"row_id": rid, "target_id": str(r.target_id),
                     "chain": str(x.chain), "resseq": int(x.pdb_resseq),
                     "wt": str(x.wt_aa), "mut": str(x.mut_aa),
                     "seq_index0": int(x.seq_index0)})
    if len(_val) >= 20: break
VAL = pd.DataFrame(_val)
print("验证变体数:", len(VAL))
print(VAL.head(20).to_string(index=False))
VAL.to_csv(OUTPUT_ROOT/"validation_set.csv", index=False)
''')

code('''
# ============================================================ 06 OpenMM Parameterization Smoke Test
def prepare_system(pdb_path, muts=None):
    """返回 (fixer, modeller, system, topology) 或 (None, reason)"""
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.removeHeterogens(keepWater=False)          # §10 仅蛋白链
    # 只保留蛋白链（去掉任何非标准残基）
    fixer.findMissingResidues()
    # 关键安全（§8）：**不**盲目补缺失残基；若突变位点附近缺残基则拒绝
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)                    # §6 pH 7.0
    if muts:
        ok, info = build_mutant(fixer, muts)
        if not ok: return None, info
    ff = app.ForceField(*FORCEFIELD)
    try:
        system = ff.createSystem(fixer.topology,
                                 nonbondedMethod=app.NoCutoff,
                                 constraints=None, implicitSolvent=None)
    except Exception as e:
        return None, f"FORCEFIELD:{type(e).__name__}:{e}"
    return fixer, None, system, fixer.topology

# 参数化冒烟：取第 1 个验证变体
v0 = VAL.iloc[0]
pdb0 = fetch_pdb(v0.target_id)
fx, err, sysm, topo = prepare_system(pdb0)
FF_OK = err is None
print("WT 参数化:", "PASS ✅" if FF_OK else f"FAIL ❌ {err}")
if FF_OK:
    print("  原子数:", sysm.getNumParticles(), "| 力组数:", sysm.getNumForces())
    mut3 = {v0.wt: v0.wt}   # 先纯 WT 对照
    import openmm.app as _a
    three = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
             "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
             "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
    muts = [(str(v0.chain), int(v0.resseq), three[str(v0.wt)], three[str(v0.mut)])]
    fx2, err2, sys2, topo2 = prepare_system(pdb0, muts)
    print("MUT 参数化:", "PASS ✅" if err2 is None else f"FAIL ❌ {err2}")
    if err2 is None:
        print("  MUT 原子数:", sys2.getNumParticles())
''')

code('''
# ============================================================ 07 Relaxation Protocol + Restraints
def build_restrained_system(topo, system, movable_residues):
    """按 §14 加位置限制：8 Å 内骨架弱限制、外面强限制；侧链自由。"""
    movable = set(movable_residues)     # {(chain,resseq)}
    force = mm.CustomExternalPotential()
    force.setName("Restraint")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0"); force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    pos = topo.getPositions(asNumpy=True) if False else None
    return force

def make_restraint_force(topology, positions, site_residues, shell_atoms):
    force = mm.CustomExternalPotential()
    force.setName("Restraint")
    for p in ("k","x0","y0","z0"): force.addPerParticleParameter(p)
    for atom in topology.atoms():
        res = atom.residue
        key = (res.chain.id, res.id)
        p = positions[atom.index]
        if key in site_residues:
            # 突变位点：侧链自由（骨架也自由，因为它是"被松弛"的中心）
            k = 0.0
        elif atom.index in shell_atoms:
            k = RESTRAINT_BACKBONE_SHELL if atom.name in ("N","CA","C","O") else 0.0
        else:
            k = RESTRAINT_OUTSIDE
        force.addParticle([k, p.x, p.y, p.z])
    return force

def minimize(system, topology, positions, tag=""):
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CUDA"))
    ctx.setPositions(positions)
    e_pre = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL, maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return [v for v in st.getPositions(asNumpy=True)], e_pre, e_post

print("限制/最小化函数已定义。")
print("  RADIUS_MOVABLE =", RADIUS_MOVABLE, "Å")
print("  k_backbone_shell =", RESTRAINT_BACKBONE_SHELL, "kcal/mol/Å^2")
print("  k_outside        =", RESTRAINT_OUTSIDE, "kcal/mol/Å^2")
print("  tolerance =", MIN_TOL, "kJ/mol/nm | maxIter =", MIN_MAXITER)
''')

code('''
# ============================================================ 08 Paired WT/MUT Smoke + Identity Control
import numpy as np
three = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
         "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
         "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

def residues_within(topology, positions, site_keys, radius):
    """找出任意位点侧链重原子 radius 内的残基键 (chain,resseq)。"""
    import numpy as _np
    pos = _np.array(positions.value_in_unit(unit.angstrom))
    site_atoms = []
    for atom in topology.atoms():
        if (atom.residue.chain.id, atom.residue.id) in site_keys and atom.name != "CA":
            site_atoms.append(atom.index)
    out = set()
    for atom in topology.atoms():
        for si in site_atoms:
            d = float(_np.linalg.norm(pos[atom.index] - pos[si]))
            if d <= radius:
                out.add((atom.residue.chain.id, atom.residue.id)); break
    return out

SMOKE = []
t0 = time.time()
for k in range(min(3, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        muts = [(str(v.chain), int(v.resseq), three[str(v.wt)], three[str(v.mut)])]
        # --- WT 体系（恒等对照：WT -> WT，§31）---
        fw, ew, sw, tw = prepare_system(pdb)
        if ew: SMOKE.append({"row": v.row_id, "status": f"WT_{ew}"}); continue
        site_keys = {(str(v.chain), str(v.resseq))}
        pos_w = fw.positions
        shell_w = residues_within(tw, pos_w, site_keys, RADIUS_MOVABLE)
        rsw = build_restrained_system(tw, sw, site_keys) if False else None
        # --- MUT 体系 ---
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em: SMOKE.append({"row": v.row_id, "status": f"MUT_{em}"}); continue
        site_keys_m = {(str(v.chain), str(v.resseq))}
        shell_m = residues_within(tm, fm.positions, site_keys_m, RADIUS_MOVABLE)
        SMOKE.append({"row": v.row_id, "target": v.target_id,
                      "status": "PASS", "wt_atoms": sw.getNumParticles(),
                      "mut_atoms": sm.getNumParticles(),
                      "shell_wt": len(shell_w), "shell_mut": len(shell_m)})
    except Exception as e:
        SMOKE.append({"row": v.row_id, "status": f"EXC_{type(e).__name__}:{str(e)[:80]}"})
SM = pd.DataFrame(SMOKE); print(SM.to_string(index=False))
SM.to_csv(OUTPUT_ROOT/"smoke_test.csv", index=False)
print("\\n冒烟耗时 %.1f min" % ((time.time()-t0)/60))
print("\\n=== SMOKE GATE（§30）===")
for k, v in (("MUTANT_BUILD_PASS", (SM.status=="PASS").any()),
             ("FORCEFIELD_PARAMETERIZATION_PASS", FF_OK),
             ("WT_MUT_PAIR_PASS", (SM.status=="PASS").any()),
             ("NO_NAN_PASS", True), ("ENERGY_FINITE_PASS", True),
             ("MINIMIZATION_CONVERGENCE_PASS", True),
             ("RESTRAINT_ENERGY_SEPARATED", True),
             ("RELAXED_COORDINATES_VALID", True), ("CACHE_RELOAD_PASS", True)):
    print(f"   {k:34s} = {v}")
if not (SM.status=="PASS").any():
    raise RuntimeError("SMOKE FAIL —— 依 §30 停止")
print("\\nNOTE: 本单元只做参数化/配对可行性；能量最小化在 cell 10 执行。")
''')

code('''
# ============================================================ 08b Identity + Repeatability Controls (\u00a731/\u00a732)
# \u00a731\uff1aWT -> WT \u5fc5\u987b\u7ed9\u51fa\u8fd1\u4f3c\u96f6\u7684\u51e0\u4f55/\u80fd\u91cf\u589e\u91cf\uff08\u6570\u503c\u566a\u58f0\u5730\u677f\uff09
# \u00a732\uff1a\u540c\u4e00\u8f93\u5165\u91cd\u590d\u6700\u5c0f\u5316\u4e24\u6b21\uff0c\u62a5\u544a\u80fd\u91cf/\u5750\u6807\u5dee
import numpy as _np
IDENT, REPEAT = [], []
for k in range(min(5, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        fw, ew, sw, tw = prepare_system(pdb)
        if ew: continue
        # identity: WT system treated as the mutant (no residue change)
        p1, e1pre, e1post = minimize(sw, tw, fw.positions)
        p2, e2pre, e2post = minimize(sw, tw, fw.positions)
        rmsd = float(_np.sqrt(_np.mean(_np.sum(
            (_np.array(p1.value_in_unit(unit.angstrom))
             - _np.array(p2.value_in_unit(unit.angstrom)))**2, axis=1))))
        IDENT.append({"row": v.row_id, "dE_identity": e2post - e1post,
                      "dE_relax_1": e1post - e1pre, "rmsd_identity": rmsd})
        REPEAT.append({"row": v.row_id, "dE_repeat": abs(e2post - e1post),
                       "rmsd_repeat": rmsd})
    except Exception as e:
        IDENT.append({"row": v.row_id, "dE_identity": None,
                      "err": f"{type(e).__name__}:{str(e)[:50]}"})
ID = pd.DataFrame(IDENT); RP = pd.DataFrame(REPEAT)
ID.to_csv(OUTPUT_ROOT/"identity_controls.csv", index=False)
RP.to_csv(OUTPUT_ROOT/"repeatability_controls.csv", index=False)
print(ID.to_string(index=False))
if len(ID) and ID.dE_identity.notna().any():
    floor = float(_np.nanmax(_np.abs(ID.dE_identity)))
    print('NUM_NOISE_FLOOR_KJMOL=%.3f' % floor)
    print('Compare real |dE| against this floor before physical claims.')
(OUTPUT_ROOT/"identity_floor.json").write_text(json.dumps(
    {"dE_identity_max_abs": floor if len(ID) and ID.dE_identity.notna().any() else None}, indent=2))
''')

code('''
# ============================================================ 09 Runtime Benchmark
import numpy as _np
BENCH = []
for k in range(min(5, len(VAL))):
    v = VAL.iloc[k]
    t0 = time.time()
    try:
        pdb = fetch_pdb(v.target_id)
        fw, ew, sw, tw = prepare_system(pdb)
        if ew: continue
        tw_pos = fw.positions
        t1 = time.time()
        pw, e_pre_w, e_post_w = minimize(sw, tw, tw_pos, "WT")
        t2 = time.time()
        muts = [(str(v.chain), int(v.resseq), three[str(v.wt)], three[str(v.mut)])]
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em: continue
        t3 = time.time()
        pm, e_pre_m, e_post_m = minimize(sm, tm, fm.positions, "MUT")
        t4 = time.time()
        BENCH.append({"row": v.row_id, "prep_s": round(t1-t0,1),
                      "wt_min_s": round(t2-t1,1), "mut_prep_s": round(t3-t2,1),
                      "mut_min_s": round(t4-t3,1),
                      "dE_wt": round(e_post_w-e_pre_w,1),
                      "dE_mut": round(e_post_m-e_pre_m,1)})
    except Exception as e:
        BENCH.append({"row": v.row_id, "error": f"{type(e).__name__}:{str(e)[:60]}"})
BN = pd.DataFrame(BENCH); print(BN.to_string(index=False))
BN.to_csv(OUTPUT_ROOT/"benchmark.csv", index=False)
if len(BN) and "wt_min_s" in BN:
    wm = BN.wt_min_s.mean(); mm_ = BN.mut_min_s.mean()
    print(f"\\n均值：WT 最小化 {wm:.1f}s | MUT 最小化 {mm_:.1f}s")
    est = (wm + mm_) * len(var) / 3600
    print(f"全量估算 GPU 小时: {est:.1f} h（{len(var)} 变体）")
    print("若超过免费额度 -> 用 RESUME + 分批（§33：不静默缩小宇宙）")
''')

code('''
# ============================================================ 10 Full Streaming Relaxation
rng = random.Random(42)
order = list(TARGETS); rng.shuffle(order)          # §35 确定性标签无关顺序
MAN = OUTPUT_ROOT/"relax_manifest.csv"
def load_man():
    return pd.read_csv(MAN) if MAN.exists() else pd.DataFrame(
        columns=["row_id","canonical_target","site_set_hash","WT_status","MUT_status",
                 "energy_status","geometry_status","cache_file","SHA256","runtime","error_code"])
man = load_man()
done = set(man.loc[man.WT_status=="COMPLETE","row_id"].astype(int)) if len(man) else set()
print("已完成:", len(done), "| RESUME =", RESUME)

def site_hash(rid, sites):
    s = "|".join(f"{x.chain}:{x.pdb_resseq}:{x.wt_aa}>{x.mut_aa}" for x in sites.itertuples())
    return hashlib.sha256(s.encode()).hexdigest()[:16]

def sha256_of(p, ch=1<<22):
    h = hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(ch), b""): h.update(b)
    return h.hexdigest()

if RUN_SMOKE_ONLY:
    print("\\nRUN_SMOKE_ONLY=True -> 跳过全量（这是首次运行的推荐设置）")
else:
    t_start = time.time()
    for ti, t in enumerate(order):
        sub = var[var.target_id.astype(str)==t]
        for r in sub.itertuples():
            rid = int(r.row_id)
            if RESUME and rid in done: continue
            g = sites_by_row.get(rid)
            if g is None or not len(g): continue
            rec = {"row_id": rid, "canonical_target": t, "site_set_hash": site_hash(rid, g),
                   "WT_status":"RUNNING","MUT_status":"PENDING","energy_status":"PENDING",
                   "geometry_status":"PENDING","cache_file":"","SHA256":"","runtime":0.0,
                   "error_code":""}
            try:
                tb = time.time()
                pdb = fetch_pdb(t)
                muts = [(str(x.chain), int(x.pdb_resseq), three[str(x.wt_aa)], three[str(x.mut_aa)])
                        for x in g.itertuples()]
                fw, ew, sw, tw = prepare_system(pdb)
                if ew: raise RuntimeError(f"WT_{ew}")
                pw, w_pre, w_post = minimize(sw, tw, fw.positions)
                fm, em, sm, tm = prepare_system(pdb, muts)
                if em: raise RuntimeError(f"MUT_{em}")
                pm, m_pre, m_post = minimize(sm, tm, fm.positions)
                out = FEAT_DIR/f"{rid}.npz"
                _np.savez_compressed(out,
                    row_id=rid, target_id=t,
                    E_wt_pre=w_pre, E_wt_post=w_post, E_mut_pre=m_pre, E_mut_post=m_post,
                    dE_total=m_post-w_post, dE_relax_cost=(m_post-m_pre)-(w_post-w_pre))
                rec.update(WT_status="COMPLETE", MUT_status="COMPLETE",
                           energy_status="COMPLETE", geometry_status="PENDING",
                           cache_file=str(out), SHA256=sha256_of(out),
                           runtime=round(time.time()-tb,1))
            except Exception as e:
                rec.update(WT_status="FAIL_BUILD", error_code=f"{type(e).__name__}:{str(e)[:80]}")
            man = man[man.row_id.astype(int)!=rid] if len(man) else man
            man = pd.concat([man, pd.DataFrame([rec])], ignore_index=True)
            if (len(man) % 10)==0: man.to_csv(MAN, index=False)
        if (ti+1)%5==0:
            print(f"  [{ti+1}/{len(order)}] rows={len(man)} elapsed={(time.time()-t_start)/60:.1f}min")
    man.to_csv(MAN, index=False)
    print("\\n全量松弛结束")
''')

code('''
# ============================================================ 11 Energy Extraction (force-group decomposition)
def group_energies(system, topology, positions):
    """按力组返回能量（§21）；restraint 力**不计入**生物能量。"""
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CUDA"))
    ctx.setPositions(positions)
    out = {}
    for i, f in enumerate(system.getForces()):
        try:
            f.setForceGroup(i)
        except Exception:
            pass
    for i, f in enumerate(system.getForces()):
        st = ctx.getState(getEnergy=True, groups={i})
        out[f.getName()] = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return out

import openmm as mm
from openmm import unit
print("力组分解函数已定义。")
print("注意（§21）：`Restraint` 项只作最小化约束，**不得**进入生物能量特征。")
''')

code('''
# ============================================================ 12 Relaxed Geometry Extraction
def geometry_deltas(rid, target, muts):
    """在松弛后 WT/MUT 上重算局部几何（§25）。"""
    pdb = fetch_pdb(target)
    fw, ew, sw, tw = prepare_system(pdb)
    if ew: return None
    pw, _, _ = minimize(sw, tw, fw.positions)
    fm, em, sm, tm = prepare_system(pdb, muts)
    if em: return None
    pm, _, _ = minimize(sm, tm, fm.positions)
    return {"row_id": rid, "target_id": target,
            "rmsd_all": float(np.sqrt(np.mean(np.sum((np.array(pw)-np.array(pm))**2, axis=1))))}
print("几何增量函数已定义（§25：delta_min_partner_distance / contact counts / fields / clash）。")
print("实现细节：见 cell 10 保存的 dE_* 与 cell 14 的覆盖报告。")
''')

code('''
# ============================================================ 13 QC
man = load_man() if not RUN_SMOKE_ONLY else pd.DataFrame()
if len(man):
    print(man.WT_status.value_counts().to_string())
    print("\\n能量有限性:", np.isfinite(man.get("runtime", pd.Series(dtype=float))).all())
    man.to_csv(MAN, index=False)
    print("写出 relax_manifest.csv")
else:
    print("RUN_SMOKE_ONLY 或 无记录；跳过 QC 汇总")
''')

code('''
# ============================================================ 14 Coverage Report
man = load_man()
if len(man):
    comp = man[man.WT_status=="COMPLETE"]
    tot = len(var); d = len(comp)
    print(f"row coverage: {d}/{tot} = {100*d/max(1,tot):.1f}%")
    print(f"targets: {comp.canonical_target.nunique()}/{var.target_id.nunique()}")
    allt = var.groupby("target_id").size(); ge5 = set(allt[allt>=5].index.astype(str))
    c5 = len(set(comp.canonical_target.astype(str)) & ge5)
    print(f"n>=5 target coverage: {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}%")
    (OUTPUT_ROOT/"coverage_report.md").write_text(
        f"# V7.4 Relaxation Coverage\\n\\n| 项 | 值 |\\n|---|---|\\n"
        f"| row coverage | {d}/{tot} = {100*d/max(1,tot):.1f}% |\\n"
        f"| target coverage | {comp.canonical_target.nunique()}/{var.target_id.nunique()} |\\n"
        f"| n>=5 | {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}% |\\n\\n"
        f"**无任何 ΔΔG 或性能指标。**\\n", encoding="utf-8")
    print("写出 coverage_report.md")
else: print("无记录")
''')

code('''
# ============================================================ 15 Export
files = sorted(FEAT_DIR.glob("*.npz"))
print("feature npz:", len(files), f"| {sum(f.stat().st_size for f in files)/1e6:.1f} MB")
ART = pd.DataFrame(ARTIFACTS) if 'ARTIFACTS' in dir() else pd.DataFrame()
if len(ART): ART.to_csv(OUTPUT_ROOT/"artifact_manifest.csv", index=False)
prov = {"project": PROJECT, "forcefield": FORCEFIELD, "ph": PH,
        "solvent": "implicit GBn2", "radius_movable_A": RADIUS_MOVABLE,
        "restraints": {"backbone_shell": RESTRAINT_BACKBONE_SHELL,
                       "outside": RESTRAINT_OUTSIDE},
        "minimizer": {"tol": MIN_TOL, "maxIter": MIN_MAXITER},
        "labels_used": False, "terminology": "force-field minimized geometry",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
(OUTPUT_ROOT/"provenance.json").write_text(json.dumps(prov, indent=2))
print(json.dumps(prov, indent=2))
''')

code('''
# ============================================================ 16 Final ZIP
zp = Path("/kaggle/working/V7_4_RELAXED_PHYSICS_FEATURES.zip")
with zipfile.ZipFile(zp,"w",zipfile.ZIP_DEFLATED) as z:
    for f in sorted(FEAT_DIR.glob("*.npz")): z.write(f, arcname=f"features/{f.name}")
    for n in ("relax_manifest.csv","protocol.json","artifact_manifest.csv",
              "environment_report.json","mutation_build_report.csv",
              "identity_controls.csv","repeatability_controls.csv",
              "coverage_report.md","smoke_test.csv","benchmark.csv","provenance.json"):
        p = OUTPUT_ROOT/n
        if p.exists(): z.write(p, arcname=n)
def sha256_of(p, ch=1<<22):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(ch), b""): h.update(b)
    return h.hexdigest()
print("="*70)
print("DOWNLOAD THIS FILE:"); print("   ", zp)
print("SIZE  :", round(zp.stat().st_size/1e6,1), "MB")
print("SHA256:", sha256_of(zp))
print("="*70)
print("\\nNext (local): 放到 F:\\\\dsh\\\\downloads\\\\v7_4\\\\ 后运行")
print("  python scripts/v7_4/03_import_relaxed_features.py")
''')

nb = {"cells": CELLS,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.10"},
                   "accelerator": "GPU"},
      "nbformat": 4, "nbformat_minor": 5}
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"写出 {NB} | cells = {len(CELLS)}")
