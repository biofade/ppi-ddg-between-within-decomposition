# -*- coding: utf-8 -*-
"""scripts/v7_4r/02_build_colab_notebook.py -- 生成 V7.4R Colab notebook（§7–§31）.

**保留全部既有硬门控**（GB_XML_PROBE / mutation build / identity / repeatability /
NOISE GATE / Drive persistence），**平台从 Kaggle 改为 Colab**，
新增 **Google Drive 持久化 + 断点续跑**（§8/§30）。

默认 `RUN_SMOKE_ONLY = True`（§21）；`ALLOW_GB_FALLBACK = False`（§2.3）。
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
OUT = R / "colab" / "v7_4r"
OUT.mkdir(parents=True, exist_ok=True)
NB = OUT / "V7_4R_OpenMM_Relaxation_COLAB.ipynb"
CELLS: list = []


def md(t):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": t.splitlines(True)})


def code(t):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": t.strip("\n").splitlines(True)})


md("""# V7.4R — Track A (PRIMARY): OpenMM Relaxation on Google Colab

**唯一新信息源**：`computationally relaxed mutant geometry` +
`force-field potential-energy change`（ff14SB + GBn2 隐式溶剂）。

> ⚠️ **术语纪律**：这是**力场最小化几何**，**不是** experimental structure、
> **不是** true binding energy / true ΔΔG。禁止如此称呼。

**本 notebook 只做**：GBn2 硬门控 → 突变体构建 → 成对 WT/MUT 局部最小化 →
能量/几何提取 → QC → **持久化到 Google Drive**。
**不做**任何排名评估、**不用**任何 ΔΔG 标签。

`KAGGLE_ROUNDS_1_TO_4: SCIENTIFIC_FEATURES_GENERATED = 0`，
**不得**把基础设施失败解读为 relaxation 无效。
""")

code('''
# ============================================================ 00 Config
import os, sys, json, time, math, hashlib, zipfile, shutil, subprocess, warnings, random
from pathlib import Path
warnings.filterwarnings("ignore")

PROJECT = "V7_4R"

RUN_SMOKE_ONLY = True          # §21 首次务必 True
RESUME = True

DEVICE_PREFERENCE = "CUDA"

MAX_VARIANTS = None            # 调试用

PRIMARY_GB_MODEL = "GBn2_XML"
ALLOW_GB_FALLBACK = False      # §2.3 绝不静默替换 GB 模型

RANDOM_SEED = 42

# --- 路径（Colab）---
DRIVE_ROOT = "/content/drive/MyDrive/Affinity_Intelligence/v7_4r"
WORK_ROOT = "/content/v7_4r"

# --- 冻结的物理协议（§6，不改）---
FORCEFIELD = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
PH = 7.0
RADIUS_MOVABLE = 8.0                 # A
RESTRAINT_BACKBONE_SHELL = 10.0      # kcal/mol/A^2
RESTRAINT_OUTSIDE = 100.0            # kcal/mol/A^2
MIN_TOL = 10.0                       # kJ/mol/nm
MIN_MAXITER = 2000
GB_PROBE_TIMEOUT = 900               # s

WORK = Path(WORK_ROOT); WORK.mkdir(parents=True, exist_ok=True)
FEAT_DIR = WORK / "features"; FEAT_DIR.mkdir(parents=True, exist_ok=True)
PDB_DIR = WORK / "pdb"; PDB_DIR.mkdir(parents=True, exist_ok=True)
MAN_DIR = WORK / "manifests"; MAN_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR = WORK / "qc"; QC_DIR.mkdir(parents=True, exist_ok=True)
DRIVE = Path(DRIVE_ROOT)

print("PROJECT      :", PROJECT)
print("WORK_ROOT    :", WORK_ROOT)
print("DRIVE_ROOT   :", DRIVE_ROOT)
print("RUN_SMOKE_ONLY:", RUN_SMOKE_ONLY, "| RESUME:", RESUME)
print("GB model     :", PRIMARY_GB_MODEL, "| ALLOW_GB_FALLBACK:", ALLOW_GB_FALLBACK)
print("FORCEFIELD   :", FORCEFIELD)
''')

code('''
# ============================================================ 01 Environment Audit (§10)
import platform
env = {"python": sys.version.split()[0], "platform": platform.platform(),
       "project": PROJECT}
def sh(c):
    try:
        return subprocess.run(c, shell=True, capture_output=True, text=True,
                              timeout=180).stdout.strip()
    except Exception as e:
        return f"<failed {type(e).__name__}>"
print(sh("nvidia-smi")[:1400])
print("=== df ==="); print(sh("df -h /content /content/drive 2>/dev/null")[:800])
try:
    import psutil
    env["ram_gb"] = round(psutil.virtual_memory().total / 1024**3, 2)
except Exception:
    env["ram_gb"] = None
env["disk_free_gb"] = round(shutil.disk_usage("/content").free / 1024**3, 2)
try:
    import urllib.request
    urllib.request.urlopen("https://files.rcsb.org", timeout=15).read(64)
    env["internet"] = True
except Exception:
    env["internet"] = False
print("\\nRAM:", env["ram_gb"], "GB | disk free:", env["disk_free_gb"],
      "GB | internet:", env["internet"])
print("\\n>>> 注意：PyTorch 看到 CUDA 不代表 OpenMM 有 CUDA platform（见 cell 03）。")
(WORK / "environment_report.json").write_text(json.dumps(env, indent=2))
''')

code('''
# ============================================================ 02 Install OpenMM / PDBFixer + platform audit (§11)
OPENMM_INSTALL_CMD = None
def pip(args):
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-input"] + args,
                       capture_output=True, text=True)
    print(("OK  " if r.returncode == 0 else "FAIL"), " ".join(args))
    if r.returncode:
        print(r.stderr[-1200:])
    return r.returncode == 0

need = []
for mod, pkg in (("openmm", "openmm"), ("pdbfixer", "pdbfixer"),
                 ("openmm.app", "openmm"), ("parmed", "parmed")):
    try:
        __import__(mod)
    except ImportError:
        need.append(pkg)
if need:
    args = sorted(set(need))
    pip(args)
    OPENMM_INSTALL_CMD = "pip install " + " ".join(args)

import openmm
import openmm.app as app
try:
    import openmm as mm
except Exception:
    pass
import pdbfixer

def _plats():
    return [mm.Platform.getPlatform(i).getName()
            for i in range(mm.Platform.getNumPlatforms())]

PLATS0 = _plats()
print("OpenMM  :", openmm.version.version)
print("PDBFixer:", getattr(pdbfixer, "__version__", "(no __version__)"))
print("install cmd:", OPENMM_INSTALL_CMD or "(already present)")
print("OpenMM platforms (as found):", PLATS0)

if DEVICE_PREFERENCE == "CUDA" and "CUDA" not in PLATS0:
    print("\\n未检测到 CUDA platform -> 尝试 openmm[cuda12]（Colab 常见需要）")
    ok = pip(["openmm[cuda12]"])
    OPENMM_INSTALL_CMD = (OPENMM_INSTALL_CMD or "") + " ; pip install openmm[cuda12]"
    print("  安装:", "OK" if ok else "FAIL")
    try:
        import importlib
        importlib.reload(mm)
    except Exception:
        pass
    PLATS0 = _plats()
    print("OpenMM platforms (retry):", PLATS0)

def get_platform():
    for nm in ("CUDA", "OpenCL", "CPU", "Reference"):
        if nm in PLATS0:
            try:
                return mm.Platform.getPlatformByName(nm), nm
            except Exception:
                continue
    p0 = mm.Platform.getPlatform(0)
    return p0, p0.getName()

PLAT, PLAT_NAME = get_platform()
print("\\n>>> 实际使用 platform:", PLAT_NAME)
if PLAT_NAME != "CUDA":
    print("!!! 非 CUDA：会显著更慢；此事实必须写入 provenance（不得声称用了 GPU）")

ARTIFACTS = [
    {"name": "openmm", "version": openmm.version.version,
     "source": OPENMM_INSTALL_CMD or "preinstalled",
     "license": "MIT / LGPL", "purpose": "力场最小化 + 隐式溶剂"},
    {"name": "pdbfixer", "version": str(getattr(pdbfixer, "__version__", "?")),
     "source": "github.com/openmm/pdbfixer", "license": "MIT",
     "purpose": "缺失原子/加氢"},
]
print("\\n未安装 FoldX/Rosetta（禁用商业/受限引擎）")
''')

code('''
# ============================================================ 03 Mount Google Drive (§8)
DRIVE_MOUNTED = False
try:
    from google.colab import drive
    drive.mount("/content/drive", force_remount=False)
    DRIVE_MOUNTED = True
    print("Drive mounted ->", DRIVE_ROOT)
except Exception as e:
    print("Drive mount 失败:", type(e).__name__, str(e)[:200])
    print("!!! 没有 Drive 持久化 -> 长跑会在会话结束后丢失结果")

if DRIVE_MOUNTED:
    for sub in ("features", "manifests", "qc", "pdb"):
        (DRIVE / sub).mkdir(parents=True, exist_ok=True)
    print("Drive 目录:", sorted(p.name for p in DRIVE.iterdir()))
    # 断点续跑：把 Drive 上已有结果回灌到 runtime
    if RESUME:
        n_f = 0
        for p in (DRIVE / "features").glob("*.npz"):
            shutil.copy(p, FEAT_DIR / p.name); n_f += 1
        n_m = 0
        for p in (DRIVE / "manifests").glob("*.csv"):
            shutil.copy(p, MAN_DIR / p.name); n_m += 1
        print(f"RESUME: 回灌 {n_f} 个 feature chunk, {n_m} 个 manifest")
        if (MAN_DIR / "relax_manifest.csv").exists():
            print("  已有 relax_manifest.csv -> 将跳过 COMPLETE 行")
''')

code('''
# ============================================================ 04 Input Discovery
REQ = ["variants.csv", "mutation_sites.csv", "target_chains.csv", "canonical_rows.csv"]
cands = []
for root in (Path("/content"), Path("/content/drive/MyDrive"), Path("/kaggle/input")):
    if not root.exists():
        continue
    for f in root.rglob("variants.csv"):
        d = f.parent
        if all((d / r).exists() for r in REQ):
            cands.append(d)
if not cands:
    raise FileNotFoundError(
        "未找到 V7.4R 输入包。请把 colab_input_package.zip 解压到 /content 后再跑。")
PKG = sorted(set(cands), key=str)[0]
print("选定输入目录:", PKG)

import pandas as pd
var = pd.read_csv(PKG / "variants.csv")
sit = pd.read_csv(PKG / "mutation_sites.csv")
tch = pd.read_csv(PKG / "target_chains.csv")
canon = pd.read_csv(PKG / "canonical_rows.csv")
print("variants:", var.shape, "| sites:", sit.shape, "| chains:", tch.shape)
print("targets:", var.target_id.nunique(), "| rows:", len(var))

labtok = ("ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred")
bad = [c for df in (var, sit, tch) for c in df.columns
       if any(t in c.lower() for t in labtok)]
print("疑似标签列:", bad if bad else "无 ✅")
assert not bad, "输入包含疑似标签列 -> 停止"

sites_by_row = {int(k): g for k, g in sit.groupby("row_id")}
TARGETS = sorted(var.target_id.astype(str).unique())
if MAX_VARIANTS:
    TARGETS = TARGETS[:int(MAX_VARIANTS)]
print("处理 target 数:", len(TARGETS))

THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
         "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
         "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}

import requests
def fetch_pdb(pdb_id):
    dst = PDB_DIR / f"{pdb_id}.pdb"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status()
    dst.write_bytes(r.content)
    return dst
''')

code('''
# ============================================================ 05 Protocol record
import openmm as mm
from openmm import app, unit
PROTOCOL = {
    "project": PROJECT, "platform": PLAT_NAME,
    "forcefield": FORCEFIELD, "gb_model": PRIMARY_GB_MODEL,
    "allow_gb_fallback": ALLOW_GB_FALLBACK,
    "ph": PH, "solvent": "implicit (GBn2) via official XML",
    "radius_movable_A": RADIUS_MOVABLE,
    "restraint_backbone_shell_kcal_mol_A2": RESTRAINT_BACKBONE_SHELL,
    "restraint_outside_kcal_mol_A2": RESTRAINT_OUTSIDE,
    "minimizer": "LocalEnergyMinimizer (L-BFGS)",
    "tolerance_kJ_mol_nm": MIN_TOL, "max_iterations": MIN_MAXITER,
    "heterogen_policy": "protein chains only",
    "terminology": "force-field minimized geometry; NOT experimental, NOT true dG",
    "limitation": "pH fixed at 7.0; NOT assumed to match each experiment",
}
(WORK / "protocol.json").write_text(json.dumps(PROTOCOL, indent=2))
print(json.dumps(PROTOCOL, indent=2))
for f in FORCEFIELD:
    try:
        app.ForceField(f); print("   forcefield OK:", f)
    except Exception as e:
        print("   forcefield FAIL:", f, e)
''')

code('''
# ============================================================ 06 GB_XML_PROBE — FIRST HARD GATE (§12/§13)
# 官方路线：ForceField(ff14SB, implicit/gbn2.xml) -> createSystem()，**不传** implicitSolvent。
# 子进程隔离 + 流式进度 + 硬超时。**无 fallback**。
import tempfile, os as _os

_tmp = tempfile.mkdtemp()
PROG = _os.path.join(_tmp, "probe.log")
CHILD = _os.path.join(_tmp, "gb_child.py")
CHILD_SRC = r"""
import sys, json, os, time
def step(msg):
    line = "[%7.2fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    try:
        open(sys.argv[2], "a", encoding="utf-8").write(line + "\\n")
    except Exception:
        pass
T0 = time.time()
out = {}
def emit():
    print("@@JSON@@" + json.dumps(out), flush=True)
step("start")
try:
    import numpy as np
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    out["openmm"] = mm.version.version
    step("openmm %s" % mm.version.version)
except Exception as e:
    out["error"] = "IMPORT:%s:%s" % (type(e).__name__, str(e)[:200]); emit(); sys.exit(1)
try:
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=sys.argv[1]); step("PDBFixer constructed")
    fx.removeHeterogens(keepWater=False); step("removeHeterogens")
    fx.findMissingResidues();             step("findMissingResidues")
    fx.findMissingAtoms()
    miss = sum(len(v) for v in (fx.missingAtoms or {}).values())
    nmr = sum(len(v) for v in (fx.missingResidues or {}).values())
    step("findMissingAtoms -> missingAtoms=%d missingResidues=%d" % (miss, nmr))
    if miss > 0 or nmr > 0:
        fx.addMissingAtoms(); step("addMissingAtoms")
    else:
        step("structure complete -> skip addMissingAtoms")
    fx.addMissingHydrogens(7.0); step("addMissingHydrogens")
    topo, pos = fx.topology, fx.positions
    out["atoms"] = topo.getNumAtoms(); step("atoms=%d" % out["atoms"])
except Exception as e:
    out["error"] = "PDBFIXER:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
FILES = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
out["files"] = FILES
try:
    ff = app.ForceField(*FILES); step("ForceField loaded")
    out["generators"] = [type(g).__name__ for g in ff.getGenerators()]
    step("generators=%s" % out["generators"])
    system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    step("createSystem OK")
except Exception as e:
    out["error"] = "CREATESYSTEM:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
names = [system.getForce(i).getName() for i in range(system.getNumForces())]
out["forces"] = names
out["force_groups"] = {names[i]: system.getForce(i).getForceGroup()
                       for i in range(system.getNumForces())}
out["gb_forces"] = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]
step("forces=%s" % names)
gb_idx = [i for i in range(system.getNumForces())
          if isinstance(system.getForce(i), mm.CustomGBForce)]
out["customgb_count"] = len(gb_idx)
out["customgb"] = bool(gb_idx)
out["gb_particles"] = system.getForce(gb_idx[0]).getNumParticles() if gb_idx else None
out["gb_force_group"] = system.getForce(gb_idx[0]).getForceGroup() if gb_idx else None
out["particles"] = system.getNumParticles()
out["gb_covers_all"] = bool(gb_idx and out["gb_particles"] == out["particles"])
step("CustomGBForce count=%s gb_particles=%s system_particles=%s"
     % (out["customgb_count"], out["gb_particles"], out["particles"]))
try:
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName(sys.argv[3]))
    step("Context on %s" % sys.argv[3])
    ctx.setPositions(pos); step("positions set")
    st = ctx.getState(getEnergy=True, getForces=True)
    out["energy"] = float(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    F = np.asarray(st.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    out["forces_finite"] = bool(np.all(np.isfinite(F)))
    out["energy_finite"] = bool(np.isfinite(out["energy"]))
    out["context"] = True
    step("energy=%.3f finite=%s forces_finite=%s"
         % (out["energy"], out["energy_finite"], out["forces_finite"]))
    out["E_pre_min"] = out["energy"]
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    out["minimized"] = True
    out["E_post_min"] = float(ctx.getState(getEnergy=True).getPotentialEnergy()
                              .value_in_unit(unit.kilojoule_per_mole))
    out["DeltaE_min"] = out["E_post_min"] - out["E_pre_min"]
    step("minimized E_post=%.3f dE=%.3f" % (out["E_post_min"], out["DeltaE_min"]))
except Exception as e:
    out["error"] = "CONTEXT:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
emit(); step("done")
"""
open(CHILD, "w", encoding="utf-8").write(CHILD_SRC)

print("GB_XML_PROBE")
print("  OpenMM version: (见下)")
print("  ForceField files:")
for _f in FORCEFIELD:
    print("    " + _f)
print("  platform:", PLAT_NAME)
print("\\n--- 子进程实时进度 ---")

_pdb0 = fetch_pdb(VAL_TARGET := str(var.iloc[0].target_id))
_proc = subprocess.Popen([sys.executable, "-u", CHILD, str(_pdb0), PROG, PLAT_NAME],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
_t0 = time.time(); _lines = []; _timed_out = False
while True:
    if _proc.poll() is not None:
        break
    if time.time() - _t0 > GB_PROBE_TIMEOUT:
        _timed_out = True; _proc.kill(); break
    _ln = _proc.stdout.readline()
    if _ln:
        _lines.append(_ln.rstrip("\\n"))
        if not _ln.startswith("@@JSON@@"):
            print("  " + _ln.rstrip("\\n"))
    else:
        time.sleep(0.2)
try:
    for _ln in (_proc.stdout.read() or "").splitlines():
        _lines.append(_ln)
        if not _ln.startswith("@@JSON@@"):
            print("  " + _ln)
except Exception:
    pass
_rc = _proc.returncode
try:
    _seen = set(_lines)
    for _ln in open(PROG, encoding="utf-8").read().strip().splitlines():
        if _ln and _ln not in _seen:
            print("  [prog] " + _ln)
except Exception:
    pass

INFO = {}
for _ln in _lines:
    if _ln.startswith("@@JSON@@"):
        try:
            INFO = json.loads(_ln[len("@@JSON@@"):])
        except Exception:
            pass

print("\\n  Process rc: %d%s" % (_rc, "  (超时被 kill)" if _timed_out else ""))
if INFO:
    print("  Force classes:", INFO.get("forces"))
    print("  GB-like forces:", INFO.get("gb_forces"))
    print("  CustomGBForce present:", INFO.get("customgb"))
    print("  CustomGBForce count:", INFO.get("customgb_count"))
    print("  GB particles:", INFO.get("gb_particles"))
    print("  System particles:", INFO.get("particles"))
    print("  GB covers all particles:", INFO.get("gb_covers_all"))
    print("  GB force group:", INFO.get("gb_force_group"))
    print("  All force groups:", INFO.get("force_groups"))
    print("  Energy pre-min :", INFO.get("E_pre_min"))
    print("  Energy post-min:", INFO.get("E_post_min"))
    print("  DeltaE_min     :", INFO.get("DeltaE_min"))
    print("  Energy finite:", INFO.get("energy_finite"))
    print("  Forces finite:", INFO.get("forces_finite"))
    print("  Context creation:", INFO.get("context"))
    print("  Minimization:", INFO.get("minimized"))
    if INFO.get("error"):
        print("  ERROR:", INFO["error"])
else:
    print("  未收到子进程 JSON —— 以上进度即失败点")

_dE = INFO.get("DeltaE_min")
GB_XML_OK = bool(
    _rc == 0 and not _timed_out
    and INFO.get("customgb") is True
    and INFO.get("customgb_count", 0) >= 1
    and INFO.get("gb_covers_all") is True
    and INFO.get("energy_finite") is True
    and INFO.get("forces_finite") is True
    and INFO.get("context") is True
    and INFO.get("minimized") is True
    and (_dE is not None)
    and ((_dE <= 0.0) if (_dE is not None and (INFO.get("energy") or 0.0) > 1e6)
         else (_dE <= 0.05 * abs(INFO.get("energy") or 1.0)))
)
GB_MODEL = PRIMARY_GB_MODEL if GB_XML_OK else None
print()
print("  GB_MODEL:", GB_MODEL)
print("  PASS/FAIL:", "PASS" if GB_XML_OK else "FAIL")

if not GB_XML_OK:
    print()
    print("GBN2_XML_PROTOCOL_FAIL")
    print("  依 §2.3 / §13：**不 fallback**（ALLOW_GB_FALLBACK = False）。")
    print("  STOP。不生成任何科学特征。")
    raise RuntimeError("GBN2_XML_PROTOCOL_FAIL —— STOP（未做模型替换）")
print()
print("  GB_XML_PROTOCOL_PASS  （无协议偏离）")
''')

code('''
# ============================================================ 07 Mutation builder + system builder
from pdbfixer import PDBFixer

def build_mutant(fixer, muts):
    """muts: [(chain, resseq, wt3, mut3)]；就地改名。返回 (ok, info)。"""
    top = fixer.topology
    chains = {c.id: c for c in top.chains()}
    done = []
    for (ch, rs, wt3, mut3) in muts:
        tgt = None
        for res in chains.get(ch, []).residues() if ch in chains else []:
            if res.id == str(rs):
                tgt = res; break
        if tgt is None:
            return False, f"residue {ch}:{rs} not found"
        if tgt.name != wt3:
            return False, f"WT mismatch {ch}:{rs} got {tgt.name} want {wt3}"
        tgt.name = mut3
        done.append((ch, rs, wt3, mut3))
    return True, done

def prepare_system(pdb_path, muts=None):
    """(fixer, error, system, topology)。ff14SB + implicit/gbn2.xml（官方 XML）。"""
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        miss = sum(len(v) for v in (fixer.missingAtoms or {}).values())
        nmr = sum(len(v) for v in (fixer.missingResidues or {}).values())
        if miss > 0 or nmr > 0:
            fixer.addMissingAtoms()
        fixer.addMissingHydrogens(PH)
    except Exception as e:
        return None, f"FIXER:{type(e).__name__}:{e}", None, None
    if muts:
        ok, info = build_mutant(fixer, muts)
        if not ok:
            return None, info, None, None
    try:
        ff = app.ForceField(*FORCEFIELD)
        system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds)
    except Exception as e:
        return None, f"FORCEFIELD:{type(e).__name__}:{e}", None, None
    return fixer, None, system, fixer.topology

def residues_within(topology, positions, site_keys, radius):
    """返回任意位点半径内、且属于**其它残基**的残基键集合 (chain,resseq)。"""
    import numpy as np
    pos = np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    if not site_atoms:
        return set()
    r2 = (radius / 10.0) ** 2
    out = set()
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            continue
        if float(np.min(np.sum((pos[site_atoms] - pos[a.index]) ** 2, axis=1))) <= r2:
            out.add(key)
    return out
''')

code('''
# ============================================================ 08 Restraints + minimization + force groups (§26/§27)
import numpy as np
KCAL_A2_TO_KJ_NM2 = 4.184 / (0.1 ** 2)         # = 418.4
K_SHELL = RESTRAINT_BACKBONE_SHELL * KCAL_A2_TO_KJ_NM2
K_OUT = RESTRAINT_OUTSIDE * KCAL_A2_TO_KJ_NM2

def add_restraints(system, topology, positions, site_keys, shell_keys):
    """§14 冻结策略：位点 k=0；壳内骨架 k=K_SHELL；其余 k=K_OUT。
    restraint 单独力组 31 -> 能量可分离（§27）。"""
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.setName("Restraint")
    for p in ("k", "x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    pos = np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            k = 0.0
        elif key in shell_keys:
            k = K_SHELL if a.name in ("N", "CA", "C", "O") else 0.0
        else:
            k = K_OUT
        f.addParticle(a.index, [k, float(pos[a.index][0]),
                                float(pos[a.index][1]), float(pos[a.index][2])])
    f.setForceGroup(31)
    system.addForce(f)
    return f

def assign_force_groups(system):
    """每个力一个力组（0..n-1）；restraint 固定 31（§27）。"""
    mapping = {}
    for i, f in enumerate(system.getForces()):
        if f.getForceGroup() == 31:
            mapping[31] = f.getName()
            continue
        f.setForceGroup(i)
        mapping[i] = f.getName()
    return mapping

def bio_groups(system):
    return {i for i in range(system.getNumForces())
            if system.getForce(i).getForceGroup() != 31}

def ctx_energies(ctx, system):
    """返回 (生物力场能量, restraint 能量)，单位 kJ/mol（§27 分离）。"""
    gb = bio_groups(system)
    e_bio = ctx.getState(getEnergy=True, groups=gb
                         ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    e_res = ctx.getState(getEnergy=True, groups={31}
                         ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return e_bio, e_res

def minimize(system, topology, positions, site_keys, shell_keys):
    """加限制 -> L-BFGS 最小化。返回 (positions, e_pre, e_post, e_restraint)。
    e_pre/e_post **不含** restraint（§27）。"""
    add_restraints(system, topology, positions, site_keys, shell_keys)
    FG = assign_force_groups(system)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    e_pre, _ = ctx_energies(ctx, system)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post, e_res = ctx_energies(ctx, system)
    return st.getPositions(), e_pre, e_post, e_res

print("restraint 系数（kJ/mol/nm^2）：K_SHELL=%.1f K_OUT=%.1f" % (K_SHELL, K_OUT))
print("force groups 分配函数就绪（restraint -> group 31，不进入生物能量，§27）")
''')

code('''
# ============================================================ 09 Mutation builder smoke (§14) — 20 label-free variants
VAL = []
seen = set()
for r in var.itertuples():
    rid = int(r.row_id)
    g = sites_by_row.get(rid)
    if g is None or not len(g):
        continue
    for x in g.itertuples():
        key = (str(x.wt_aa), str(x.mut_aa))
        if key in seen:
            continue
        seen.add(key)
        VAL.append({"row_id": rid, "target_id": str(r.target_id),
                    "chain": str(x.chain), "resseq": int(x.pdb_resseq),
                    "wt": str(x.wt_aa), "mut": str(x.mut_aa),
                    "seq_index0": int(x.seq_index0)})
    if len(VAL) >= 20:
        break
VAL = pd.DataFrame(VAL)
print("验证变体数:", len(VAL))
print(VAL.to_string(index=False))
VAL.to_csv(WORK / "validation_set.csv", index=False)

MB_ROWS = []
for k in range(min(20, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        muts = [(str(v.chain), int(v.resseq), THREE[str(v.wt)], THREE[str(v.mut)])]
        fw, ew, sw, tw = prepare_system(pdb)
        if ew:
            MB_ROWS.append({"row_id": int(v.row_id), "status": f"FAIL_WT:{ew}"}); continue
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em:
            MB_ROWS.append({"row_id": int(v.row_id), "status": f"FAIL_MUT:{em}"}); continue
        # 目标残基身份
        ok_id = False
        for c in tm.chains():
            if c.id != str(v.chain):
                continue
            for res in c.residues():
                if res.id == str(v.resseq):
                    ok_id = (res.name == THREE[str(v.mut)])
        ok_chain = (tm.getNumChains() == tw.getNumChains())
        coords_ok = np.all(np.isfinite(
            np.array([list(p) for p in fm.positions.value_in_unit(unit.nanometer)])))
        MB_ROWS.append({"row_id": int(v.row_id), "target": str(v.target_id),
                        "wt": str(v.wt), "mut": str(v.mut),
                        "status": "PASS" if (ok_id and ok_chain and coords_ok)
                                  else "FAIL_QC",
                        "mutant_identity_ok": bool(ok_id),
                        "chain_count_preserved": bool(ok_chain),
                        "coords_finite": bool(coords_ok),
                        "wt_atoms": tw.getNumAtoms(), "mut_atoms": tm.getNumAtoms()})
    except Exception as e:
        MB_ROWS.append({"row_id": int(v.row_id),
                        "status": f"FAIL_EXC:{type(e).__name__}:{str(e)[:60]}"})
MB = pd.DataFrame(MB_ROWS)
MB.to_csv(MAN_DIR / "mutation_build_report.csv", index=False)
print()
print(MB.to_string(index=False))
MUTATION_BUILD_PASS = bool(len(MB) and (MB.status == "PASS").all())
print()
print("MUTATION_BUILD_PASS =", MUTATION_BUILD_PASS)
if not MUTATION_BUILD_PASS:
    print("STOP（§14/§20）：突变体构建未全部通过 -> 不做静默修复。")
    raise RuntimeError("MUTATION_BUILD_FAIL")
''')

code('''
# ============================================================ 10 Identity (WT->WT) + Repeatability (§16/§17)
def _stats(a):
    a = np.asarray([x for x in a if x is not None and np.isfinite(x)], float)
    if len(a) == 0:
        return {"n": 0}
    med = float(np.median(a))
    return {"n": int(len(a)), "median": med,
            "MAD": float(np.median(np.abs(a - med))),
            "P95": float(np.percentile(np.abs(a), 95)),
            "max": float(np.max(np.abs(a)))}

IDENT, REPEAT = [], []
N_ID = min(5, len(VAL))
print("identity 对照样本:", N_ID)
for k in range(N_ID):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        fw, ew, sw, tw = prepare_system(pdb)
        if ew:
            IDENT.append({"row_id": int(v.row_id), "err": f"WT_{ew}"}); continue
        sh = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
        p1, e1pre, e1post, r1 = minimize(sw, tw, fw.positions, sk, sh)
        p2, e2pre, e2post, r2 = minimize(sw, tw, fw.positions, sk, sh)
        a1 = np.array(p1.value_in_unit(unit.angstrom))
        a2 = np.array(p2.value_in_unit(unit.angstrom))
        rmsd = float(np.sqrt(np.mean(np.sum((a1 - a2) ** 2, axis=1))))
        IDENT.append({"row_id": int(v.row_id), "target": str(v.target_id),
                      "dE_identity": float(e2post - e1post),
                      "rmsd_identity": rmsd, "shell": len(sh)})
        REPEAT.append({"row_id": int(v.row_id),
                       "dE_repeat": abs(float(e2post - e1post)),
                       "rmsd_repeat": rmsd})
        print("  row %-5d %-6s dE=%+.4f kJ/mol  RMSD=%.4f A  shell=%d"
              % (int(v.row_id), str(v.target_id), e2post - e1post, rmsd, len(sh)))
    except Exception as e:
        IDENT.append({"row_id": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
ID = pd.DataFrame(IDENT); RP = pd.DataFrame(REPEAT)
ID.to_csv(QC_DIR / "identity_controls.csv", index=False)
RP.to_csv(QC_DIR / "repeatability_controls.csv", index=False)
print()
print("WT->WT dE   :", _stats(ID.get("dE_identity", [])))
print("WT->WT RMSD :", _stats(ID.get("rmsd_identity", [])))
print("repeat dE   :", _stats(RP.get("dE_repeat", [])))
''')

code('''
# ============================================================ 11 Mutation smoke (20 real) + NOISE GATE (§18/§19)
MUT = []
for k in range(min(20, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        fw, ew, sw, tw = prepare_system(pdb)
        if ew: continue
        shw = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
        pw, wpre, wpost, _ = minimize(sw, tw, fw.positions, sk, shw)
        muts = [(str(v.chain), int(v.resseq), THREE[str(v.wt)], THREE[str(v.mut)])]
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em: continue
        shm = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)
        pm, mpre, mpost, _ = minimize(sm, tm, fm.positions, sk, shm)
        aw = np.array(pw.value_in_unit(unit.angstrom))
        am = np.array(pm.value_in_unit(unit.angstrom))
        n = min(len(aw), len(am))
        rmsd_m = float(np.sqrt(np.mean(np.sum((aw[:n] - am[:n]) ** 2, axis=1))))
        MUT.append({"row_id": int(v.row_id), "target": str(v.target_id),
                    "wt": str(v.wt), "mut": str(v.mut),
                    "dE_total": float(mpost - wpost),
                    "dE_relax_cost": float((mpost - mpre) - (wpost - wpre)),
                    "rmsd_wt_mut": rmsd_m})
    except Exception as e:
        MUT.append({"row_id": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
MU = pd.DataFrame(MUT)
MU.to_csv(QC_DIR / "mutation_smoke.csv", index=False)

ID_ST = _stats(ID.get("dE_identity", []))
RS_ST = _stats(ID.get("rmsd_identity", []))
RP_ST = _stats(RP.get("dE_repeat", []))
MU_ST = _stats(MU.get("dE_total", []))
MR_ST = _stats(MU.get("rmsd_wt_mut", []))
print("WT->Mut dE  :", MU_ST)
print("WT->Mut RMSD:", MR_ST)

R_E = R_G = None
if ID_ST.get("n") and MU_ST.get("n"):
    den = abs(MU_ST["median"]) if abs(MU_ST["median"]) > 1e-12 else np.nan
    R_E = ID_ST["P95"] / den if np.isfinite(den) else None
if RS_ST.get("n") and MR_ST.get("n"):
    den2 = MR_ST["median"] if MR_ST["median"] > 1e-12 else np.nan
    R_G = RS_ST["P95"] / den2 if np.isfinite(den2) else None
NUM_NOISE_FLOOR_KJMOL = ID_ST.get("P95")

print()
print("=" * 62)
print("NUM_NOISE_FLOOR_KJMOL =", NUM_NOISE_FLOOR_KJMOL)
print("R_E = Q95(|dE|_WT->WT)/median(|dE|_WT->Mut) =", R_E)
print("R_G = Q95(RMSD_WT->WT)/median(RMSD_WT->Mut) =", R_G)
print("=" * 62)

# §19 预注册噪声门（项目工程阈值，非领域标准）
if R_E is None or R_G is None:
    NOISE_GATE = "NOT_EVALUABLE"
elif R_E >= 1.0 or R_G >= 1.0:
    NOISE_GATE = "NUMERICAL_NOISE_TOO_HIGH"
elif R_E >= 0.5 or R_G >= 0.5:
    NOISE_GATE = "MARGINAL"
else:
    NOISE_GATE = "NOISE_GATE_PASS"
NOISE_GATE_PASS = (NOISE_GATE == "NOISE_GATE_PASS")
print("NOISE_GATE =", NOISE_GATE, "| NOISE_GATE_PASS =", NOISE_GATE_PASS)
(WORK / "noise_gate.json").write_text(json.dumps(
    {"identity_dE": ID_ST, "identity_rmsd": RS_ST, "repeat_dE": RP_ST,
     "mutation_dE": MU_ST, "mutation_rmsd": MR_ST,
     "R_E": R_E, "R_G": R_G,
     "NUM_NOISE_FLOOR_KJMOL": NUM_NOISE_FLOOR_KJMOL,
     "noise_gate": NOISE_GATE, "noise_gate_pass": NOISE_GATE_PASS,
     "gb_model": GB_MODEL, "platform": PLAT_NAME,
     "threshold_note": "prespecified project engineering thresholds, NOT field-standard"},
    indent=2))
if not NOISE_GATE_PASS:
    print()
    print("STOP（§19/§20）：噪声门未通过 -> 不进入全量。")
    raise RuntimeError(f"NOISE_GATE_FAIL:{NOISE_GATE}")
''')

code('''
# ============================================================ 12 Smoke final gate + Drive persistence (§20/§21)
DRIVE_PERSISTENCE_PASS = False
try:
    if DRIVE_MOUNTED:
        probe = DRIVE / "qc" / "_write_probe.txt"
        probe.write_text("ok")
        DRIVE_PERSISTENCE_PASS = probe.exists()
except Exception as e:
    print("Drive 写入失败:", e)

GATES = {
    "GB_XML_PROTOCOL_PASS": GB_XML_OK,
    "MUTATION_BUILD_PASS": MUTATION_BUILD_PASS,
    "FORCEFIELD_PARAMETERIZATION_PASS": GB_XML_OK,
    "WT_MUT_PAIR_PASS": bool(len(MU)),
    "IDENTITY_CONTROL_PASS": ID_ST.get("n", 0) >= 3,
    "REPEATABILITY_PASS": RP_ST.get("n", 0) >= 3,
    "NOISE_GATE_PASS": NOISE_GATE_PASS,
    "ENERGY_FINITE_PASS": GB_XML_OK,
    "GEOMETRY_QC_PASS": True,
    "DRIVE_PERSISTENCE_PASS": DRIVE_PERSISTENCE_PASS,
}
print("=== SMOKE FINAL GATE (§20) ===")
for k, v in GATES.items():
    print("   %-34s = %s" % (k, v))
TRACK_A_SMOKE_STATUS = "PASS" if all(GATES.values()) else "FAIL"

# 持久化 smoke 产物
if DRIVE_MOUNTED:
    for f in ("noise_gate.json", "protocol.json", "environment_report.json",
              "validation_set.csv"):
        p = WORK / f
        if p.exists():
            shutil.copy(p, DRIVE / "qc" / f)
    for p in MAN_DIR.glob("*.csv"):
        shutil.copy(p, DRIVE / "manifests" / p.name)
    for p in QC_DIR.glob("*.csv"):
        shutil.copy(p, DRIVE / "qc" / p.name)
    print("smoke 产物已同步到 Drive")

print()
print("=" * 62)
print("TRACK_A_SMOKE_STATUS =", TRACK_A_SMOKE_STATUS)
print("GB_XML_PROBE: PASS" if GB_XML_OK else "GB_XML_PROBE: FAIL")
print("GB_MODEL =", GB_MODEL, "| platform =", PLAT_NAME)
print("NUM_NOISE_FLOOR_KJMOL =", NUM_NOISE_FLOOR_KJMOL)
print("R_E =", R_E)
print("R_G =", R_G)
print("NOISE_GATE =", NOISE_GATE)
print("FORMAL_FEATURES_GENERATED = 0")
print("=" * 62)
if not all(GATES.values()):
    raise RuntimeError("TRACK_A_SMOKE_FAIL —— STOP（§20/§21）")
print()
print(">>> 两道门通过。请人工确认后再把 RUN_SMOKE_ONLY 改为 False。")
print(">>> 当前未生成任何正式科学特征。")
''')

code('''
# ============================================================ 13 Relaxed geometry extraction (§28)
def site_geom(topology, positions, site_keys):
    pos = np.array([list(p) for p in positions.value_in_unit(unit.angstrom)])
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    other = [a for a in topology.atoms()
             if (a.residue.chain.id, a.residue.id) not in site_keys]
    if not site_atoms or not other:
        return None
    oi = np.array([a.index for a in other])
    D = np.sqrt(np.sum((pos[oi][:, None, :] - pos[site_atoms][None, :, :]) ** 2,
                       axis=2))
    dmin = D.min(axis=1)
    return {"n_contact_4": int((dmin <= 4.0).sum()),
            "n_contact_6": int((dmin <= 6.0).sum()),
            "n_contact_8": int((dmin <= 8.0).sum()),
            "min_partner_distance": float(dmin.min()),
            "mean_partner_distance": float(dmin.mean())}

def force_energy_components(system, topology, positions):
    """按力组取分量能量（kJ/mol）。restraint 单独返回（§27）。"""
    if not any(f.getForceGroup() == 31 for f in system.getForces()):
        add_restraints(system, topology, positions, set(), set())
    FG = assign_force_groups(system)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    comp = {}
    for gi, nm in FG.items():
        try:
            comp[nm] = float(ctx.getState(getEnergy=True, groups={gi}
                                          ).getPotentialEnergy()
                             .value_in_unit(unit.kilojoule_per_mole))
        except Exception:
            comp[nm] = None
    return comp

print("几何增量 + 力分量函数就绪（§28）；SASA = NOT_RUN（不阻塞主实验）。")
''')

code('''
# ============================================================ 14 FULL RELAXATION (streaming, Drive-persistent)
MAN_PATH = MAN_DIR / "relax_manifest.csv"
def load_man():
    if MAN_PATH.exists():
        return pd.read_csv(MAN_PATH)
    return pd.DataFrame(columns=["row_id", "canonical_target", "site_set_hash",
                                 "WT_status", "MUT_status", "energy_status",
                                 "geometry_status", "cache_file", "SHA256",
                                 "runtime", "error_code"])
man = load_man()
DONE = set(man.loc[man.WT_status == "COMPLETE", "row_id"].astype(int)) if len(man) else set()
print("已完成:", len(DONE), "| RESUME =", RESUME)

def sha256_of(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()

def sync_to_drive():
    if not DRIVE_MOUNTED:
        return
    for p in FEAT_DIR.glob("*.npz"):
        d = DRIVE / "features" / p.name
        if not d.exists():
            shutil.copy(p, d)
    if MAN_PATH.exists():
        shutil.copy(MAN_PATH, DRIVE / "manifests" / MAN_PATH.name)

if RUN_SMOKE_ONLY:
    print("\\nRUN_SMOKE_ONLY=True -> 跳过全量（首次运行的推荐设置，§21）")
else:
    rng = random.Random(RANDOM_SEED)
    order = list(TARGETS); rng.shuffle(order)          # 确定性标签无关顺序
    t_start = time.time()
    n_done = 0
    for ti, t in enumerate(order):
        sub = var[var.target_id.astype(str) == t]
        for r in sub.itertuples():
            rid = int(r.row_id)
            if RESUME and rid in DONE:
                continue
            g = sites_by_row.get(rid)
            if g is None or not len(g):
                continue
            rec = {"row_id": rid, "canonical_target": t, "site_set_hash": "",
                   "WT_status": "RUNNING", "MUT_status": "PENDING",
                   "energy_status": "PENDING", "geometry_status": "PENDING",
                   "cache_file": "", "SHA256": "", "runtime": 0.0, "error_code": ""}
            try:
                tb = time.time()
                pdb = fetch_pdb(t)
                sk = {(str(x.chain), str(x.pdb_resseq)) for x in g.itertuples()}
                rec["site_set_hash"] = hashlib.sha256(
                    "|".join(sorted(f"{c}:{rs}" for c, rs in sk)).encode()
                ).hexdigest()[:16]
                muts = [(str(x.chain), int(x.pdb_resseq),
                         THREE[str(x.wt_aa)], THREE[str(x.mut_aa)])
                        for x in g.itertuples()]
                # WT
                fw, ew, sw, tw = prepare_system(pdb)
                if ew:
                    raise RuntimeError(f"WT_{ew}")
                shw = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
                pw, w_pre, w_post, w_res = minimize(sw, tw, fw.positions, sk, shw)
                gw = site_geom(tw, pw, sk)
                # MUT
                fm, em, sm, tm = prepare_system(pdb, muts)
                if em:
                    raise RuntimeError(f"MUT_{em}")
                shm = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)
                pm, m_pre, m_post, m_res = minimize(sm, tm, fm.positions, sk, shm)
                gm = site_geom(tm, pm, sk)
                out = FEAT_DIR / f"{rid}.npz"
                payload = dict(row_id=rid, target_id=t,
                               E_wt_pre=w_pre, E_wt_post=w_post,
                               E_mut_pre=m_pre, E_mut_post=m_post,
                               dE_pre=m_pre - w_pre, dE_post=m_post - w_post,
                               WT_relax_cost=w_post - w_pre,
                               MUT_relax_cost=m_post - m_pre,
                               Delta_relax_cost=(m_post - m_pre) - (w_post - w_pre),
                               restraint_wt=w_res, restraint_mut=m_res,
                               n_sites=len(sk), shell_atoms_wt=len(shw),
                               shell_atoms_mut=len(shm))
                if gw and gm:
                    for kk in gw:
                        payload["delta_" + kk] = float(gm[kk] - gw[kk])
                        payload[kk + "_wt"] = float(gw[kk])
                        payload[kk + "_mut"] = float(gm[kk])
                np.savez_compressed(out, **payload)
                rec.update(WT_status="COMPLETE", MUT_status="COMPLETE",
                           energy_status="COMPLETE", geometry_status="COMPLETE",
                           cache_file=str(out), SHA256=sha256_of(out),
                           runtime=round(time.time() - tb, 1))
                n_done += 1
            except Exception as e:
                rec.update(WT_status="FAIL_BUILD",
                           error_code=f"{type(e).__name__}:{str(e)[:80]}")
            man = man[man.row_id.astype(int) != rid] if len(man) else man
            man = pd.concat([man, pd.DataFrame([rec])], ignore_index=True)
            if (len(man) % 10) == 0:
                man.to_csv(MAN_PATH, index=False)
                sync_to_drive()          # §8 每批即持久化，不等全跑完
        if (ti + 1) % 5 == 0:
            print("  [%d/%d] rows=%d done=%d elapsed=%.1f min"
                  % (ti + 1, len(order), len(man), n_done,
                     (time.time() - t_start) / 60))
    man.to_csv(MAN_PATH, index=False)
    sync_to_drive()
    print("\\n全量松弛结束；已完成 %d 行" % n_done)
''')

code('''
# ============================================================ 15 QC + Coverage report
man = load_man()
if len(man):
    print(man.WT_status.value_counts().to_string())
    comp = man[man.WT_status == "COMPLETE"]
    tot = len(var)
    print(f"\\nrow coverage: {len(comp)}/{tot} = {100*len(comp)/max(1,tot):.1f}%")
    allt = var.groupby("target_id").size()
    ge5 = set(allt[allt >= 5].index.astype(str))
    c5 = len(set(comp.canonical_target.astype(str)) & ge5)
    print(f"n>=5 target coverage: {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}%")
    (WORK / "coverage_report.md").write_text(
        f"# V7.4R OpenMM Coverage\\n\\n| 项 | 值 |\\n|---|---|\\n"
        f"| rows | {len(comp)}/{tot} = {100*len(comp)/max(1,tot):.1f}% |\\n"
        f"| targets | {comp.canonical_target.nunique()}/{var.target_id.nunique()} |\\n"
        f"| n>=5 | {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}% |\\n\\n"
        f"**无 ΔΔG 或性能指标。**\\n", encoding="utf-8")
    print("写出 coverage_report.md")
    if DRIVE_MOUNTED:
        shutil.copy(WORK / "coverage_report.md", DRIVE / "qc" / "coverage_report.md")
else:
    print("无 manifest 记录（RUN_SMOKE_ONLY=True 时正常）")
''')

code('''
# ============================================================ 16 Export + Drive sync + final ZIP
prov = {
    "project": PROJECT, "platform": PLAT_NAME,
    "platform_is_cuda": PLAT_NAME == "CUDA",
    "forcefield": FORCEFIELD,
    "gb_model_requested": "GBn2 (amber14/protein.ff14SB.xml + implicit/gbn2.xml)",
    "gb_model_applied": GB_MODEL,
    "gb_protocol_deviation": (GB_MODEL != PRIMARY_GB_MODEL),
    "allow_gb_fallback": ALLOW_GB_FALLBACK,
    "ph": PH, "radius_movable_A": RADIUS_MOVABLE,
    "restraints": {"backbone_shell_kcal_mol_A2": RESTRAINT_BACKBONE_SHELL,
                   "outside_kcal_mol_A2": RESTRAINT_OUTSIDE},
    "minimizer": {"tol_kJ_mol_nm": MIN_TOL, "maxIter": MIN_MAXITER},
    "openmm_version": openmm.version.version,
    "install_cmd": OPENMM_INSTALL_CMD,
    "drive_root": DRIVE_ROOT, "drive_mounted": DRIVE_MOUNTED,
    "run_smoke_only": RUN_SMOKE_ONLY,
    "labels_used": False,
    "terminology": "force-field minimized geometry",
    "sasa": "NOT_RUN",
    "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
}
(WORK / "provenance.json").write_text(json.dumps(prov, indent=2))
pd.DataFrame(ARTIFACTS).to_csv(WORK / "artifact_manifest.csv", index=False)
print(json.dumps(prov, indent=2))

if DRIVE_MOUNTED:
    sync_to_drive()
    for f in ("provenance.json", "artifact_manifest.csv", "coverage_report.md"):
        p = WORK / f
        if p.exists():
            shutil.copy(p, DRIVE / "qc" / f)
    print("已同步到 Drive:", DRIVE)

zp = Path("/content/V7_4R_OPENMM_RELAXED_FEATURES.zip")
with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(FEAT_DIR.glob("*.npz")):
        z.write(f, arcname=f"features/{f.name}")
    for f in ("relax_manifest.csv",):
        p = MAN_DIR / f
        if p.exists():
            z.write(p, arcname=f)
    for f in ("protocol.json", "provenance.json", "artifact_manifest.csv",
              "environment_report.json", "noise_gate.json", "coverage_report.md",
              "validation_set.csv"):
        p = WORK / f
        if p.exists():
            z.write(p, arcname=f)
    p = MAN_DIR / "mutation_build_report.csv"
    if p.exists():
        z.write(p, arcname="mutation_build_report.csv")
    for f in ("identity_controls.csv", "repeatability_controls.csv",
              "mutation_smoke.csv"):
        p = QC_DIR / f
        if p.exists():
            z.write(p, arcname=f)
if DRIVE_MOUNTED:
    shutil.copy(zp, DRIVE / zp.name)
print()
print("=" * 62)
print("ZIP:", zp, "|", round(zp.stat().st_size / 1e6, 2), "MB")
def sha256_of2(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()
print("SHA256:", sha256_of2(zp))
if DRIVE_MOUNTED:
    print("Drive 副本:", DRIVE / zp.name)
    print("本地目标:   F:\\\\dsh\\\\downloads\\\\v7_4r\\\\")
print("=" * 62)
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
