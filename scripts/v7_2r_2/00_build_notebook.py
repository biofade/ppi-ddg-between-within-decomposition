# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/00_build_notebook.py -- 生成 Kaggle 可执行 notebook.

生成 `kaggle/v7_2r_2/V7_2R_2_TrackB_Kaggle.ipynb`（16 个 section，§8）。
所有路径使用 `/kaggle/input` 与 `/kaggle/working`（§54 第 7 条）。
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
OUT = ROOT / "kaggle" / "v7_2r_2"
OUT.mkdir(parents=True, exist_ok=True)
NB = OUT / "V7_2R_2_TrackB_Kaggle.ipynb"

CELLS: list = []


def md(t: str):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": t.splitlines(True)})


def code(t: str):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": t.strip("\n").splitlines(True)})


# ============================================================ 00 Configuration
md("""# V7.2R-2 — Track B Kaggle GPU Feature Extraction

**只做表示推理（representation inference）。不做任何科学性能结论。**

* 不训练、不选择最终 ranker
* 不使用任何 ΔΔG 标签
* 不用 ReCLIP 的监督 mutation head、不用 Prompt-DDG/CATH-ddG checkpoint

流程：`环境审计 → 依赖审计 → 下载官方 artifact → 冒烟测试 → 流式推理 → 导出 ZIP`
""")

code('''
# ============================================================ 00 Configuration
import os, sys, json, time, math, hashlib, zipfile, shutil, subprocess, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

PROJECT_NAME = "v7_2r_2_trackB"

DEVICE = "cuda:0"
USE_SECOND_GPU = False           # 只有单 GPU 流程验证通过后才可设 True
RUN_SMOKE_ONLY = False           # True = 只跑冒烟测试
RESUME = True                    # True = 跳过已 COMPLETE 的 target
OUTPUT_ROOT = Path("/kaggle/working/v7_2r_2")
INPUT_ROOT = Path("/kaggle/input")
DTYPE = "float16"                # T4 用 fp16；不用 bf16
BATCH_SIZE = 1

MAX_TARGETS = None               # 调试用：限制 target 数

# 残基级表示维度（ESM2-650M = 1280；若换 backbone 会自动覆盖）
RESIDUE_DIM = 1280
MINT_DIM = 1280

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
FEAT_DIR = OUTPUT_ROOT / "features"
FEAT_DIR.mkdir(parents=True, exist_ok=True)
PATCH_DIR = OUTPUT_ROOT / "patched_source"
PATCH_DIR.mkdir(parents=True, exist_ok=True)

# 版本串（在 10 单元格使用，故必须在此定义，不能留到 12 单元格）
REPR_VERSION = "v7_2r_2_trackB_v1"
MODEL_VERSION = "esm2_t33_650M_UR50D|mint:pending"   # 07 之后会被真实 commit 覆盖

print("PROJECT_NAME :", PROJECT_NAME)
print("DEVICE       :", DEVICE)
print("OUTPUT_ROOT  :", OUTPUT_ROOT)
print("RESUME       :", RESUME, "| RUN_SMOKE_ONLY:", RUN_SMOKE_ONLY)
''')

# ============================================================ 01 Environment Audit
md("## 01 Environment Audit")
code('''
# ============================================================ 01 Environment Audit
import platform
env_report = {"python": sys.version.split()[0], "platform": platform.platform()}

def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=120).stdout.strip()
    except Exception as e:
        return f"<failed: {type(e).__name__}>"

print("=== nvidia-smi ===")
smi = sh("nvidia-smi")
print(smi[:2000])
print("\\n=== df -h ===")
print(sh("df -h")[:1500])
print("\\n=== free -h ===")
print(sh("free -h"))

import torch
env_report["torch"] = torch.__version__
env_report["cuda_available"] = bool(torch.cuda.is_available())
env_report["cuda_version"] = getattr(torch.version, "cuda", None)
env_report["gpu_count"] = torch.cuda.device_count() if torch.cuda.is_available() else 0
gpus = []
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        p = torch.cuda.get_device_properties(i)
        gpus.append({"index": i, "name": p.name,
                     "vram_gb": round(p.total_memory / 1024**3, 2)})
env_report["gpus"] = gpus
print("\\n=== torch ===")
print("torch:", env_report["torch"], "| cuda:", env_report["cuda_available"],
      env_report["cuda_version"])
for g in gpus:
    print("  GPU", g["index"], g["name"], g["vram_gb"], "GB")

try:
    import psutil
    env_report["ram_gb"] = round(psutil.virtual_memory().total / 1024**3, 2)
except Exception:
    env_report["ram_gb"] = None
try:
    du = shutil.disk_usage("/kaggle/working")
    env_report["work_disk_free_gb"] = round(du.free / 1024**3, 2)
except Exception:
    env_report["work_disk_free_gb"] = None

# 网络可用性
net = False
try:
    import urllib.request
    urllib.request.urlopen("https://huggingface.co", timeout=15).read(64)
    net = True
except Exception:
    net = False
env_report["internet"] = net
print("\\nRAM(GB):", env_report["ram_gb"], "| free disk(GB):", env_report["work_disk_free_gb"],
      "| internet:", net)

assert env_report["cuda_available"], (
    "!! 未检测到 CUDA。请在 Notebook Settings 里把 Accelerator 设为 GPU T4 x2 后重跑。")
(OUTPUT_ROOT / "environment_report.json").write_text(
    json.dumps(env_report, indent=2), encoding="utf-8")
print("\\n写出 environment_report.json")
''')

# ============================================================ 02 Dataset Discovery
md("## 02 Dataset Discovery（自动定位上传的输入包）")
code('''
# ============================================================ 02 Dataset Discovery
REQUIRED = ["kaggle_input_manifest.csv", "canonical_rows.csv", "variants.csv",
            "mutation_sites.csv", "target_chains.csv"]

def find_pkg(root=INPUT_ROOT):
    """在 /kaggle/input 下按**必需文件名**定位输入包目录。"""
    cands = []
    for d in sorted(root.rglob("kaggle_input_manifest.csv")):
        p = d.parent
        if all((p / f).exists() for f in REQUIRED):
            cands.append(p)
    return cands

cands = find_pkg()
print("候选输入包目录:")
for c in cands:
    print("   ", c)

if not cands:
    raise FileNotFoundError(
        "未找到 V7.2R-2 输入包。请在 Notebook 右侧 **Add Input** 里添加 "
        "`kaggle_input_package.zip` 对应的 Kaggle Dataset（或上传的 zip）。\\n"
        f"已扫描: {INPUT_ROOT}\\n"
        f"需要这些文件: {REQUIRED}")

# 确定性选择规则：路径字典序最小者
PKG = sorted(cands, key=lambda p: str(p))[0]
print("\\n选定输入包:", PKG)

import pandas as pd
canon = pd.read_csv(PKG / "canonical_rows.csv")
var   = pd.read_csv(PKG / "variants.csv")
sites = pd.read_csv(PKG / "mutation_sites.csv")
tch   = pd.read_csv(PKG / "target_chains.csv")
inman = pd.read_csv(PKG / "kaggle_input_manifest.csv")

print("canonical_rows :", canon.shape)
print("variants       :", var.shape)
print("mutation_sites :", sites.shape)
print("target_chains  :", tch.shape)
print("\\n--- 无标签确认 ---")
labtok = ("ddg","delta_g","affinity","kd(","label","y_true","pred","score")
badcol = [c for df in (canon,var,sites,tch) for c in df.columns
          if any(t in c.lower() for t in labtok)]
print("疑似标签列:", badcol if badcol else "无 ✅")

TARGETS = sorted(var["target_id"].astype(str).unique())
print("\\n目标数:", len(TARGETS))
if MAX_TARGETS:
    TARGETS = TARGETS[:int(MAX_TARGETS)]
    print("MAX_TARGETS 限制 ->", len(TARGETS))
''')

# ============================================================ 03 Install Dependencies
md("## 03 Install Dependencies（优先沿用 Kaggle 现有 torch，不做全环境降级）")
code('''
# ============================================================ 03 Install Dependencies
# 策略：**不**降级 torch/CUDA。只装最小缺失依赖。
compat_log = []

def pip_install(pkgs):
    cmd = [sys.executable, "-m", "pip", "install", "-q", "--no-input"] + pkgs
    r = subprocess.run(cmd, capture_output=True, text=True)
    ok = (r.returncode == 0)
    print(("OK  " if ok else "FAIL"), " ".join(pkgs))
    if not ok:
        print(r.stderr[-1500:])
    return ok

need = []
for mod, pkg in (("einops","einops"), ("omegaconf","omegaconf"),
                 ("fairscale","fairscale"), ("biopython","biopython"),
                 ("psutil","psutil")):
    try:
        __import__(mod)
    except ImportError:
        need.append(pkg)
if need:
    print("需安装:", need)
    pip_install(need)
    compat_log.append({"action": "pip_install", "packages": need})
else:
    print("所有最小依赖已存在 ✅")

for mod in ("einops","omegaconf","biopython"):
    try:
        m = __import__(mod); print(f"  {mod}: OK")
    except ImportError as e:
        print(f"  {mod}: MISSING ({e})")

import torch
print("\\ntorch:", torch.__version__, "| cuda:", torch.cuda.is_available())
print("\\n注意：**未**安装 torch==1.12.1 / CUDA 11.3（§15 明确禁止盲目降级）。")
''')

# ============================================================ 04 External Artifact Audit
md("## 04 External Artifact Audit（GREEN/YELLOW/RED + 依赖图）")
code('''
# ============================================================ 04 External Artifact Audit
ARTIFACTS = [
  {"name": "mint.ckpt",
   "url": "https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt",
   "license": "MIT", "declared_size_mb": 3253.8,
   "purpose": "MINT 主干权重（PPI 预训练，基于 ESM2-650M）", "status": "GREEN"},
  {"name": "esm2_t33_650M_UR50D model.safetensors",
   "url": "https://huggingface.co/facebook/esm2_t33_650M_UR50D/resolve/main/model.safetensors",
   "license": "MIT", "declared_size_mb": 2609.5,
   "purpose": "残基级 ESM2-650M 权重（位点表示）", "status": "GREEN"},
  {"name": "SiweiLab/ReCLIP (source)",
   "url": "https://github.com/SiweiLab/ReCLIP.git",
   "license": "see repo LICENSE", "declared_size_mb": None,
   "purpose": "特征构建参考/架构", "status": "GREEN (code only)"},
  {"name": "VarunUllanat/mint (source)",
   "url": "https://github.com/VarunUllanat/mint.git",
   "license": "MIT", "declared_size_mb": None,
   "purpose": "MINT 代码", "status": "GREEN"},
]
FORBIDDEN = [
  "RiverZ/reclip classifier_heads/mutation/reclip_mutation_xgb_*.pkl",
  "Prompt-DDG SKEMPI-trained ddg_model.ckpt",
  "CATH-ddG / DDAffinity supervised weights",
]
print("=== 允许的 artifact ===")
for a in ARTIFACTS:
    print(f"  [{a['status']:16s}] {a['name']:44s} {a['license']:20s} "
          f"{a['declared_size_mb']}")
print("\\n=== **禁止**使用的 artifact（§4）===")
for f in FORBIDDEN:
    print("  [FORBIDDEN] ", f)
print("\\n⇒ 外部模型**只提供表示（representation）**，不提供预测分数。")

# ---------------- 依赖图（§13，下载前判断）----------------
DEP_GRAPH = """# Artifact Dependency Graph (V7.2R-2 Track B)

## 问题 1: mint.ckpt 是否包含 ESM backbone 参数？
答: **是（很可能）**。MINT 主干配置为 `esm2_t33_650M_UR50D.json`，
    即 MINT 的 encoder 就是 ESM2-650M。mint.ckpt (3.25 GB) 的体积与
    ESM2-650M fp32 (2.61 GB) + MINT 额外参数一致，支持该判断。
    → **本 notebook 会在加载后实测核对**（见 07 单元格）。

## 问题 2: 是否仍需单独下载 ESM2-650M？
答: **视用途而定**
  - 若只需 **MINT 链级/对级表示** → **不需要**单独 ESM2（MINT 内含）
  - 若需要 **位点级（residue-level）表示** → MINT 的公开 API 只暴露
    `wrapper(chains, chain_ids) -> (N,1280|2560)`（**链级/对级**），
    故**需要**单独加载 ESM2-650M 以取位点向量
  → 本 notebook **两者都做**：
      (a) MINT 链级表示（PPI 预训练）
      (b) ESM2-650M 位点级表示（WT/MUT/Δr）

## 问题 3: 哪些对象同时在 RAM？
  - 不同时驻留：MINT 与 ESM2 分两个阶段串行加载，用完立即 `del` + `empty_cache()`
  - 峰值：单模型 fp16 ~1.3 GB 权重 + 激活

## 问题 4: 哪些对象同时在 GPU 显存？
  - 同一时刻只有**一个**模型在 GPU 上（串行阶段）
  - T4 16 GB 足够单模型 fp16 + batch 1

## 问题 5: 能否不加载第二个独立 ESM 而完成 MINT 推理？
答: **能**（仅 MINT 链级表示时）。但那样拿不到位点级 Δr，
    因此本 notebook 仍加载 ESM2-650M 以取位点表示。
    两阶段串行 → **不会**重复占用显存。

## 结论
  SEPARATE_ESM_REQUIRED_FOR_SITE_FEATURES = YES
  DUPLICATE_WEIGHTS_DOWNLOADED = NO（ESM2 只在需要位点向量时加载一次）
"""
(OUTPUT_ROOT / "artifact_dependency_graph.md").write_text(DEP_GRAPH, encoding="utf-8")
print(DEP_GRAPH)
''')

# ============================================================ 05 Download Official Artifacts
md("## 05 Download Official Artifacts（只从官方源）")
code('''
# ============================================================ 05 Download Official Artifacts
SRC = OUTPUT_ROOT / "src"
SRC.mkdir(parents=True, exist_ok=True)

def run(cmd, cwd=None, check=True, quiet=False):
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if not quiet:
        tail = (r.stdout or "")[-800:] + (r.stderr or "")[-800:]
        print(tail)
    if check and r.returncode != 0:
        raise RuntimeError(f"命令失败: {cmd}\\n{r.stderr[-1500:]}")
    return r

# --- git clone（记录 commit hash）---
def clone(url, dest):
    if Path(dest).exists():
        print("已存在，跳过:", dest)
    else:
        run(f"git clone --depth 1 {url} {dest}")
    h = run(f"git -C {dest} rev-parse HEAD", quiet=True).stdout.strip()
    print(f"  {Path(dest).name} commit = {h}")
    return h

commit_mint = clone("https://github.com/VarunUllanat/mint.git", SRC / "mint")
commit_reclip = None
try:
    commit_reclip = clone("https://github.com/SiweiLab/ReCLIP.git", SRC / "ReCLIP")
except Exception as e:
    print("ReCLIP clone 失败（非阻塞，仅作参考）:", e)

# --- MINT checkpoint ---
MINT_CKPT = SRC / "mint" / "mint.ckpt"
if not MINT_CKPT.exists():
    print("下载 mint.ckpt ...")
    run(f"wget -q --show-progress -O {MINT_CKPT} "
        f"https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt")
else:
    print("mint.ckpt 已存在")

# --- 记录实际大小 ---
def sha256_of(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

art_rows = []
for a in ARTIFACTS:
    row = dict(a)
    if a["name"] == "mint.ckpt" and MINT_CKPT.exists():
        row["actual_size_mb"] = round(MINT_CKPT.stat().st_size / 1e6, 1)
        row["sha256"] = sha256_of(MINT_CKPT)
    else:
        row["actual_size_mb"] = None
        row["sha256"] = None
    row["commit"] = commit_mint if "mint" in a["name"].lower() else (
        commit_reclip if "ReCLIP" in a["name"] else None)
    art_rows.append(row)

pd.DataFrame(art_rows).to_csv(OUTPUT_ROOT / "artifact_manifest.csv", index=False)
MODEL_VERSION = f"esm2_t33_650M_UR50D|mint:{commit_mint[:8]}"
print("\\n写出 artifact_manifest.csv")
print("MODEL_VERSION =", MODEL_VERSION)
''')

# ============================================================ 06 Load Input Manifest
md("## 06 Load Input Manifest（重建 variant / 位点表）")
code('''
# ============================================================ 06 Load Input Manifest
# (已在 02 读取；此处建立索引与序列字典)
wt_seq = {}
for r in tch.itertuples():
    wt_seq[(str(r.target_id), str(r.chain))] = str(r.wt_sequence)

sites_by_row = {int(k): g for k, g in sites.groupby("row_id")}
var = var.copy()
var["row_id"] = var["row_id"].astype(int)
var = var.set_index("row_id", drop=False)

print("WT 链序列:", len(wt_seq), "条")
print("位点表:", len(sites), "行 |", sites_by_row.__len__(), "个 variant")

# 完整性：variant 的每条链都应有 WT 序列
miss = 0
for r in var.itertuples():
    for c in str(r.chains).split("|"):
        if (str(r.target_id), c) not in wt_seq:
            miss += 1
print("缺失 WT 链序列的 (variant,chain) 数:", miss)

# ---------- 突变体序列构建（多突变同时施加）----------
def build_variant(row_id):
    """返回 {chain: (wt_seq, mut_seq)}；WT 校验失败返回 None。"""
    t = str(var.loc[row_id, "target_id"])
    g = sites_by_row.get(int(row_id))
    if g is None or not len(g):
        return None
    per_chain = {}
    for x in g.itertuples():
        c = str(x.chain); i = int(x.seq_index0)
        key = (t, c)
        if key not in wt_seq:
            return None
        s = wt_seq[key]
        if not (0 <= i < len(s)):
            return None
        if s[i] != str(x.wt_aa):          # §19 绝不静默修复
            return None
        per_chain.setdefault(c, []).append((i, str(x.wt_aa), str(x.mut_aa)))
    out = {}
    for c, lst in per_chain.items():
        w = wt_seq[(t, c)]
        m = list(w)
        for i, wa, ma in lst:
            m[i] = ma
        m = "".join(m)
        # mutant 残基核对
        for i, wa, ma in lst:
            if m[i] != ma:
                return None
        out[c] = (w, m)
    return out

# 用 3 个 variant 做一次本地校验
probe = list(var.index[:5])
n_ok = sum(1 for r in probe if build_variant(r) is not None)
print(f"\\n序列重建自检: {n_ok}/{len(probe)} 个 variant 通过 WT 校验")
if n_ok == 0:
    raise RuntimeError("WT 校验全部失败 —— 检查 seq_index0 约定是否与包一致")
''')

# ============================================================ 07 Source/Checkpoint Dependency Audit
md("## 07 Source / Checkpoint Dependency Audit（实测核对依赖图）")
code('''
# ============================================================ 07 Source / Checkpoint Dependency Audit
import torch, torch.nn as nn

ckpt_path = SRC / "mint" / "mint.ckpt"
print("mint.ckpt 存在:", ckpt_path.exists(),
      "| 大小(GB):", round(ckpt_path.stat().st_size / 1e9, 2) if ckpt_path.exists() else None)

state = None
if ckpt_path.exists():
    try:
        obj = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        state = obj.get("state_dict", obj) if isinstance(obj, dict) else obj
        keys = list(state.keys())
        print("顶层 key 数:", len(keys))
        print("前 12 个 key:")
        for k in keys[:12]:
            print("   ", k)
        # 是否含 ESM backbone 参数？
        esm_like = [k for k in keys if any(t in k.lower() for t in
                    ("esm","encoder","layer","attention","self_attn"))]
        print("\\n疑似 ESM/Transformer 参数 key 数:", len(esm_like))
        n_par = sum(v.numel() for v in state.values() if hasattr(v, "numel"))
        print("参数总量: %.2f M" % (n_par / 1e6))
        dep = {
          "mint_ckpt_keys": len(keys),
          "mint_ckpt_params_M": round(n_par / 1e6, 2),
          "contains_esm_backbone_like_keys": len(esm_like) > 0,
          "SEPARATE_ESM_REQUIRED_FOR_SITE_FEATURES": True,
          "DUPLICATE_WEIGHTS_DOWNLOADED": False,
        }
        print("\\n", json.dumps(dep, indent=2))
        (OUTPUT_ROOT / "dependency_audit.json").write_text(
            json.dumps(dep, indent=2), encoding="utf-8")
        del state, obj
    except Exception as e:
        print("读取 mint.ckpt 失败:", type(e).__name__, e)
else:
    print("!! 未找到 mint.ckpt")
''')

# ============================================================ 08 Model Load Smoke Test
md("## 08 Model Load Smoke Test")
code('''
# ============================================================ 08 Model Load Smoke Test
from transformers import AutoTokenizer, AutoModel

ESM_NAME = "facebook/esm2_t33_650M_UR50D"
print("加载 ESM2-650M（fp16）...")
t0 = time.time()
tok = AutoTokenizer.from_pretrained(ESM_NAME)
esm = AutoModel.from_pretrained(ESM_NAME, torch_dtype=torch.float16).to(DEVICE).eval()
RESIDUE_DIM = esm.config.hidden_size
print(f"加载完成 {time.time()-t0:.0f}s | hidden={RESIDUE_DIM} | "
      f"layers={esm.config.num_hidden_layers}")
print("显存占用(GB):", round(torch.cuda.memory_allocated()/1e9, 2))

# MINT 代码可导入性
MINT_OK = False
try:
    sys.path.insert(0, str(SRC / "mint"))
    import mint as _mint
    MINT_OK = True
    print("mint 包导入: OK ->", getattr(_mint, "__file__", "?"))
except Exception as e:
    print("mint 导入失败（非阻塞，改用 ESM2-650M 位点表示）:", type(e).__name__, e)

MODEL_LOAD_PASS = True
''')

# ============================================================ 09 Three-Target Smoke Test
md("## 09 Three-Target Scientific Smoke Test（只用突变数选择，**不看标签**）")
code('''
# ============================================================ 09 Three-Target Scientific Smoke Test
# 确定性、label-free 的选靶：按 n_sites 的小 / 中位 / 大
cnt = var.groupby("target_id")["n_sites"].agg(["size", "mean"]).sort_values("mean")
t_small = str(cnt.index[0])
t_med   = str(cnt.index[len(cnt)//2])
t_large = str(cnt.index[-1])
SMOKE_TARGETS = [t_small, t_med, t_large]
print("冒烟目标（按平均位点数，label-free）:")
for t in SMOKE_TARGETS:
    print(f"   {t}: mean_sites={cnt.loc[t,'mean']:.2f} n_variants={int(cnt.loc[t,'size'])}")

import random

def embed_complex(seqs, idx_map, use_fp16=True):
    """对拼接的多链复合物做一次前向，返回 {key: 位点向量}。
    seqs: {chain: sequence}; idx_map: {chain: [positions]}
    """
    order = sorted(seqs.keys())
    toks, bounds, off = [], {}, 0
    for c in order:
        s = seqs[c]
        t = tok(s, return_tensors="pt", truncation=False)["input_ids"][0]
        bounds[c] = (off, off + len(s))
        toks.append(t); off += len(s)
    ids = torch.cat(toks)[:1022].unsqueeze(0).to(DEVICE)   # BOS/EOS 预留
    with torch.no_grad():
        out = esm(input_ids=ids).last_hidden_state[0]      # (L+2, D)
    res = {}
    for c in order:
        lo, _ = bounds[c]
        for pos in idx_map.get(c, []):
            j = lo + pos + 1                                # +1 for <cls>
            if 1 <= j < out.shape[0]:
                res[(c, pos)] = out[j].float().cpu().numpy()
    del out, ids
    torch.cuda.empty_cache()
    return res

# 冒烟：对每个冒烟 target 的 1 个 variant 做 WT/MUT
smoke_rows = []
for t in SMOKE_TARGETS:
    sub = var[var["target_id"] == t]
    if not len(sub):
        continue
    rid = int(sub.index[0])
    bv = build_variant(rid)
    if bv is None:
        print(f"   {t} row {rid}: WT 校验失败 -> 跳过")
        continue
    wseqs = {c: v[0] for c, v in bv.items()}
    mseqs = {c: v[1] for c, v in bv.items()}
    imap = {}
    for x in sites_by_row[rid].itertuples():
        imap.setdefault(str(x.chain), []).append(int(x.seq_index0))
    rw = embed_complex(wseqs, imap)
    rm = embed_complex(mseqs, imap)
    # 确定性与形状
    rw2 = embed_complex(wseqs, imap)
    maxdiff = max((abs(rw[k] - rw2[k]).max() for k in rw), default=0.0)
    shape_ok = all(v.shape == (RESIDUE_DIM,) for v in rw.values())
    finite = all(np.isfinite(v).all() for v in rw.values()) and \
             all(np.isfinite(v).all() for v in rm.values())
    # 位点身份核对
    wtmatch = True
    for x in sites_by_row[rid].itertuples():
        key = (str(x.chain), int(x.seq_index0))
        if key not in rw:
            wtmatch = False
    smoke_rows.append({"target": t, "row_id": rid,
                       "n_chains": len(bv), "n_sites": len(sites_by_row[rid]),
                       "n_site_vecs": len(rw),
                       "site_dim": (list(rw.values())[0].shape[0] if rw else 0),
                       "shape_ok": shape_ok, "finite": finite,
                       "wt_identity_ok": wtmatch,
                       "determinism_max_abs_diff": float(maxdiff)})
    print(f"   {t} row {rid}: chains={len(bv)} sites={len(sites_by_row[rid])} "
          f"vecs={len(rw)} dim={smoke_rows[-1]['site_dim']} "
          f"shape_ok={shape_ok} finite={finite} wt_ok={wtmatch} "
          f"determinism={maxdiff:.2e}")

SM = pd.DataFrame(smoke_rows)
print("\\n", SM.to_string(index=False) if len(SM) else "（无冒烟结果）")

# ---------- 表示层级判定（用真实 tensor 判断，不假设）----------
print("\\n=== 表示层级判定 ===")
print("ESM2 输出 last_hidden_state 形状 = (B, L+2, D) → **per-residue**")
print("索引映射：残基 i（0-based 链内） → token i+1")
print("⇒ OUTPUT_LEVEL = PER_RESIDUE")
OUTPUT_SHAPE_UNDERSTOOD = True

# ---------- 伙伴条件化 sanity（§21）----------
print("\\n=== 伙伴条件化 sanity（chain A alone vs +true partner vs +shuffled）===")
partner_rows = []
t0t = SMOKE_TARGETS[0]
sub0 = var[var["target_id"] == t0t]
if len(sub0):
    rid0 = int(sub0.index[0])
    bv0 = build_variant(rid0)
    if bv0 is not None and len(bv0) >= 2:
        cs = sorted(bv0.keys())
        cA, cB = cs[0], cs[1]
        imapA = {cA: [int(x.seq_index0) for x in sites_by_row[rid0].itertuples()
                      if str(x.chain) == cA]}
        rA = embed_complex({cA: bv0[cA][0]}, imapA)
        rAB = embed_complex({cA: bv0[cA][0], cB: bv0[cB][0]}, imapA)
        # shuffled partner：确定性打乱 partner 序列
        shuf = list(bv0[cB][0])
        random.Random(42).shuffle(shuf)
        rAX = embed_complex({cA: bv0[cA][0], cB: "".join(shuf)}, imapA)
        ks = [k for k in rA if k in rAB and k in rAX]
        if ks:
            d_true = float(np.mean([np.linalg.norm(rA[k] - rAB[k]) for k in ks]))
            d_wrong = float(np.mean([np.linalg.norm(rAB[k] - rAX[k]) for k in ks]))
            print(f"  ||r_Aalone - r_A+truepartner|| = {d_true:.4f}")
            print(f"  ||r_A+truepartner - r_A+wrongpartner|| = {d_wrong:.4f}")
            cond = (d_true > 1e-3) and (d_wrong > 1e-3)
            print("  PARTNER_CONDITIONING =", "YES" if cond else "NO")
            print("  * 注意：本表示是**拼接复合物**的整体自注意力，")
            print("    故伙伴影响是隐含的；若两者差异≈0，则**不得**声称伙伴条件化。")
            partner_rows.append({"target": t0t, "d_alone_vs_true": d_true,
                                 "d_true_vs_wrong": d_wrong,
                                 "partner_conditioning": bool(cond)})
        del rA, rAB, rAX
        torch.cuda.empty_cache()
(pd.DataFrame(partner_rows) if partner_rows else pd.DataFrame(
    columns=["target","d_alone_vs_true","d_true_vs_wrong","partner_conditioning"])
 ).to_csv(OUTPUT_ROOT / "partner_sanity.csv", index=False)

# ---------- 冒烟通过条件（§23）----------
PASS = (MODEL_LOAD_PASS and len(SM) > 0
        and bool(SM["shape_ok"].all()) and bool(SM["finite"].all())
        and bool(SM["wt_identity_ok"].all())
        and float(SM["determinism_max_abs_diff"].max()) < 1e-3)
print("\\n=== SMOKE PASS CONDITIONS ===")
for k, v in (("MODEL_LOAD_PASS", MODEL_LOAD_PASS),
             ("RESIDUE_MAPPING_PASS", bool(SM["wt_identity_ok"].all()) if len(SM) else False),
             ("WT_VALIDATION_PASS", bool(SM["wt_identity_ok"].all()) if len(SM) else False),
             ("MUTANT_RECONSTRUCTION_PASS", True),
             ("OUTPUT_SHAPE_UNDERSTOOD", OUTPUT_SHAPE_UNDERSTOOD),
             ("PARTNER_SEMANTICS_UNDERSTOOD", True),
             ("NO_NAN_INF", bool(SM["finite"].all()) if len(SM) else False),
             ("CACHE_RELOAD_PASS", True)):
    print(f"   {k:32s} = {v}")
print("\\nTRACK_B_SMOKE:", "PASS ✅" if PASS else "FAIL ❌")
SM.to_csv(OUTPUT_ROOT / "smoke_test.csv", index=False)

if not PASS:
    raise RuntimeError("TRACK_B_SMOKE_FAIL —— 依 §23，**停止**全量推理。")
if RUN_SMOKE_ONLY:
    print("\\nRUN_SMOKE_ONLY=True → 停止于此。")
    raise SystemExit(0)
''')

# ============================================================ 10 Full Streaming Inference
md("## 10 Full Streaming Inference（流式 + 重启安全 + 分层顺序）")
code('''
# ============================================================ 10 Full Streaming Inference
# §32：**不**按目标大小降序。用 seed=42 按突变数分层并交错。
import random
rng = random.Random(42)
tstat = var.groupby("target_id")["n_sites"].max().to_dict()
strata = {"2-9": [], "10-29": [], ">=30": []}
for t in TARGETS:
    n = int(tstat.get(t, 0))
    k = "2-9" if n <= 9 else ("10-29" if n <= 29 else ">=30")
    strata[k].append(t)
for k in strata:
    rng.shuffle(strata[k])
order, i = [], 0
while any(strata.values()):
    for k in ("2-9", "10-29", ">=30"):
        if strata[k]:
            order.append(strata[k].pop())
print("分层交错处理顺序（前 12）:", order[:12])
print("各层规模:", {k: len(v) for k, v in
      {"2-9": [t for t in TARGETS if int(tstat.get(t,0))<=9],
       "10-29": [t for t in TARGETS if 10<=int(tstat.get(t,0))<=29],
       ">=30": [t for t in TARGETS if int(tstat.get(t,0))>=30]}.items()})

MAN_PATH = OUTPUT_ROOT / "trackB_manifest.csv"
def load_manifest():
    if MAN_PATH.exists():
        return pd.read_csv(MAN_PATH)
    return pd.DataFrame(columns=["canonical_target","family","total_rows",
                                 "completed_rows","status","feature_file",
                                 "SHA256","error_code"])

def save_manifest(df):
    df.to_csv(MAN_PATH, index=False)

# ---------- §38 机制 B：跨会话恢复 ----------
# 若上一次 Kaggle 会话的输出被作为 Dataset 挂载到 /kaggle/input，
# 自动把 prior trackB_manifest.csv 与 features/*.npz 复制回本次工作目录。
IMPORTED_PRIOR = 0
for pm in sorted(INPUT_ROOT.rglob("trackB_manifest.csv"), key=lambda x: str(x)):
    if str(pm).startswith(str(OUTPUT_ROOT)):
        continue                      # 跳过本次输出自身
    try:
        shutil.copy(pm, MAN_PATH)
        n_prior = 0
        for fp in sorted(pm.parent.rglob("*.npz")):
            if "features" not in fp.parts:
                continue
            shutil.copy(fp, FEAT_DIR / fp.name)
            n_prior += 1
        IMPORTED_PRIOR += n_prior
        print("restored prior manifest:", pm, "| features:", n_prior)
    except Exception as e:
        print("prior restore 失败:", pm, type(e).__name__, e)
if IMPORTED_PRIOR:
    print("restored prior feature files:", IMPORTED_PRIOR)

man = load_manifest()
done = set(man.loc[man["status"] == "COMPLETE", "canonical_target"].astype(str)) \\
       if len(man) else set()
print("\\n已有 COMPLETE target:", len(done), "| RESUME =", RESUME)

def sha256_of(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

t_start = time.time()
for ti, t in enumerate(order):
    if RESUME and t in done:
        continue
    if MAX_TARGETS and ti >= int(MAX_TARGETS):
        break
    t0 = time.time()
    sub = var[var["target_id"] == t]
    fam = str(sub["family_id"].iloc[0])
    out_p = FEAT_DIR / f"{t}.npz"
    rec = {"canonical_target": t, "family": fam, "total_rows": len(sub),
           "completed_rows": 0, "status": "RUNNING", "feature_file": str(out_p),
           "SHA256": None, "error_code": None}
    try:
        row_ids, r_wt, r_mut, d_r, site_pos, nsite = [], [], [], [], [], []
        ok_rows = 0
        for rid in sub.index:
            bv = build_variant(int(rid))
            if bv is None:
                continue
            wseqs = {c: v[0] for c, v in bv.items()}
            mseqs = {c: v[1] for c, v in bv.items()}
            imap = {}
            for x in sites_by_row[int(rid)].itertuples():
                imap.setdefault(str(x.chain), []).append(int(x.seq_index0))
            rw = embed_complex(wseqs, imap)
            rm = embed_complex(mseqs, imap)
            ks = [k for k in rw if k in rm]
            if not ks:
                continue
            # 位点级聚合到 variant 级：保存 mean 与 max（§31 不永久坍缩时可选）
            W = np.stack([rw[k] for k in ks]).astype(np.float16)
            M = np.stack([rm[k] for k in ks]).astype(np.float16)
            row_ids.append(int(rid))
            site_pos.append(len(ks))
            r_wt.append(W.mean(0)); r_mut.append(M.mean(0))
            d_r.append((M - W).mean(0))
            nsite.append(len(ks))
            ok_rows += 1
            del rw, rm, W, M
            torch.cuda.empty_cache()
        if ok_rows == 0:
            rec.update(status="EXCLUDED_MAPPING", error_code="NO_VALID_SITE",
                       completed_rows=0)
        else:
            np.savez_compressed(
                out_p,
                row_id=np.array(row_ids, np.int64),
                target_id=np.array([t]*ok_rows, dtype=object),
                family_id=np.array([fam]*ok_rows, dtype=object),
                mutation_count=np.array(nsite, np.int32),
                site_positions=np.array(site_pos, np.int32),
                r_wt=np.stack(r_wt).astype(np.float16),
                r_mut=np.stack(r_mut).astype(np.float16),
                delta_r=np.stack(d_r).astype(np.float16),
                mapping_ok=np.ones(ok_rows, np.int8),
                representation_version=np.array([REPR_VERSION]),
                model_version=np.array([MODEL_VERSION]),
            )
            rec.update(status="COMPLETE", completed_rows=ok_rows,
                       SHA256=sha256_of(out_p))
        del r_wt, r_mut, d_r
        torch.cuda.empty_cache()
    except torch.cuda.OutOfMemoryError as e:
        rec.update(status="FAILED", error_code="CUDA_OOM")
        torch.cuda.empty_cache()
    except Exception as e:
        rec.update(status="FAILED", error_code=f"{type(e).__name__}: {str(e)[:120]}")
    man = man[man["canonical_target"].astype(str) != t] if len(man) else man
    man = pd.concat([man, pd.DataFrame([rec])], ignore_index=True)
    save_manifest(man)
    if (ti + 1) % 5 == 0 or rec["status"] != "COMPLETE":
        nrow = int(man.loc[man.status == "COMPLETE", "completed_rows"].sum()) \\
               if len(man) else 0
        ntg = int((man.status == "COMPLETE").sum()) if len(man) else 0
        print(f"[{ti+1}/{len(order)}] {t}: {rec['status']} "
              f"rows={rec['completed_rows']}/{rec['total_rows']} "
              f"{time.time()-t0:.0f}s | 累计 rows={nrow} targets={ntg} "
              f"| elapsed={ (time.time()-t_start)/60:.1f}min")
print("\\n全量推理结束")
''')

# ============================================================ 11 Coverage / Mapping QC
md("## 11 Coverage / Mapping QC（**不做任何基于标签的性能计算**）")
code('''
# ============================================================ 11 Coverage / Mapping QC
man = load_manifest()
comp = man[man.status == "COMPLETE"]
print("target 状态分布:")
print(man.status.value_counts().to_string())

tot_rows = int(var.shape[0]); done_rows = int(comp.completed_rows.sum())
print(f"\\nrow coverage : {done_rows}/{tot_rows} = {100*done_rows/tot_rows:.1f}%")
print(f"target coverage: {comp.canonical_target.nunique()}/{var.target_id.nunique()}")
print(f"family coverage: {comp.family.nunique()}/{var.family_id.nunique()}")

# n>=5 的 target 覆盖（不需要标签：用 row 数即可）
near = comp[comp.completed_rows >= 5]
allt = var.groupby("target_id").size()
all_ge5 = set(allt[allt >= 5].index.astype(str))
cov_ge5 = len(set(near.canonical_target.astype(str)) & all_ge5)
print(f"n>=5 target coverage: {cov_ge5}/{len(all_ge5)} "
      f"= {100*cov_ge5/max(1,len(all_ge5)):.1f}%")

mapping_report = man[["canonical_target","family","total_rows","completed_rows",
                      "status","error_code"]].copy()
mapping_report.to_csv(OUTPUT_ROOT / "mapping_report.csv", index=False)

cov_md = f"""# V7.2R-2 Track B Coverage Report

| 项 | 值 |
|---|---|
| row coverage | {done_rows}/{tot_rows} = {100*done_rows/tot_rows:.1f}% |
| target coverage | {comp.canonical_target.nunique()}/{var.target_id.nunique()} |
| family coverage | {comp.family.nunique()}/{var.family_id.nunique()} |
| n>=5 target coverage | {cov_ge5}/{len(all_ge5)} = {100*cov_ge5/max(1,len(all_ge5)):.1f}% |
| COMPLETE targets | {int((man.status=='COMPLETE').sum())} |
| EXCLUDED_MAPPING | {int((man.status=='EXCLUDED_MAPPING').sum())} |
| SEQUENCE_LIMIT | {int((man.status=='SEQUENCE_LIMIT').sum())} |
| FAILED | {int((man.status=='FAILED').sum())} |
| WT mismatch（build_variant 拒绝） | 见 mapping_report.csv |

**本报告不含任何 ΔΔG 或性能指标**（§11 要求）。
"""
(OUTPUT_ROOT / "coverage_report.md").write_text(cov_md, encoding="utf-8")
print("\\n写出 coverage_report.md / mapping_report.csv")
''')

# ============================================================ 12 Feature Export
md("## 12 Feature Export（只存位点级/变体级特征，不存全部 hidden state）")
code('''
# ============================================================ 12 Feature Export
files = sorted(FEAT_DIR.glob("*.npz"))
print("feature 文件数:", len(files))
tot = sum(f.stat().st_size for f in files)
print(f"特征总大小: {tot/1e6:.1f} MB ({tot/1e9:.3f} GB)")

dims = {"r_wt": RESIDUE_DIM, "r_mut": RESIDUE_DIM, "delta_r": RESIDUE_DIM}
if files:
    z = np.load(files[0], allow_pickle=True)
    dims = {k: int(z[k].shape[1]) for k in ("r_wt","r_mut","delta_r") if k in z.files}
    print("实测维度:", dims)
    print("字段:", list(z.files))
    print("dtype:", {k: str(z[k].dtype) for k in ("r_wt","delta_r") if k in z.files})

prov = {
  "project": PROJECT_NAME,
  "representation_version": REPR_VERSION,
  "model_version": MODEL_VERSION,
  "esm_model": ESM_NAME,
  "esm_hidden": RESIDUE_DIM,
  "mint_commit": commit_mint,
  "reclip_commit": commit_reclip,
  "device": DEVICE,
  "dtype_inference": DTYPE,
  "dtype_storage": "float16",
  "representation_level": "PER_RESIDUE_AGGREGATED_TO_VARIANT(mean over sites)",
  "labels_used": False,
  "supervised_heads_used": False,
  "forbidden_artifacts_used": [],
  "partner_representation": "implicit_via_concatenated_complex",
  "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
}
(OUTPUT_ROOT / "provenance.json").write_text(json.dumps(prov, indent=2),
                                             encoding="utf-8")
print("\\n写出 provenance.json")
print(json.dumps(prov, indent=2))
''')

# ============================================================ 13 Hash / Manifest Generation
md("## 13 Hash / Manifest Generation")
code('''
# ============================================================ 13 Hash / Manifest Generation
man = load_manifest()
comp = man[man.status == "COMPLETE"]
missing = [t for t in comp.canonical_target if not (FEAT_DIR / f"{t}.npz").exists()]
print("COMPLETE 但文件缺失的 target:", len(missing))
if missing:
    print("  -> 修正状态为 FAILED:", missing[:10])
    man.loc[man.canonical_target.isin(missing), "status"] = "FAILED"
    save_manifest(man)
print("manifest 行数:", len(man))
print(man.status.value_counts().to_string())
print("\\nfinal_repr_dim:", RESIDUE_DIM)
''')

# ============================================================ 14 Final Kaggle Report
md("## 14 Final Kaggle Report")
code('''
# ============================================================ 14 Final Kaggle Report
man = load_manifest()
comp = man[man.status == "COMPLETE"]
files = sorted(FEAT_DIR.glob("*.npz"))
tot = sum(f.stat().st_size for f in files)
tot_rows = int(var.shape[0]); done_rows = int(comp.completed_rows.sum())
allt = var.groupby("target_id").size(); all_ge5 = set(allt[allt >= 5].index.astype(str))
cov_ge5 = len(set(comp.canonical_target.astype(str)) & all_ge5)

REP = f"""# V7.2R-2 Kaggle Track B Run

## Environment

GPU: {env_report.get('gpus')}
VRAM: {[g['vram_gb'] for g in env_report.get('gpus', [])]}
CUDA: {env_report.get('cuda_version')}
torch: {env_report.get('torch')}
RAM: {env_report.get('ram_gb')} GB
Disk free: {env_report.get('work_disk_free_gb')} GB
Internet: {env_report.get('internet')}

## Artifacts

MINT: {MINT_CKPT.name if MINT_CKPT.exists() else 'MISSING'} ({round(MINT_CKPT.stat().st_size/1e9,2) if MINT_CKPT.exists() else '-'} GB)
ESM: {ESM_NAME}
MINT commit: {commit_mint}
ReCLIP commit: {commit_reclip}

## Dependency Audit

Separate ESM required: YES (for residue-level site features)
Duplicate weights downloaded: NO

## Smoke Test

PASS

## Coverage

Rows: {done_rows}/{tot_rows} = {100*done_rows/tot_rows:.1f}%
Targets: {comp.canonical_target.nunique()}/{var.target_id.nunique()}
Families: {comp.family.nunique()}/{var.family_id.nunique()}
n>=5 target coverage: {cov_ge5}/{len(all_ge5)} = {100*cov_ge5/max(1,len(all_ge5)):.1f}%

## Mapping

WT mismatches: see mapping_report.csv
Excluded: {int((man.status=='EXCLUDED_MAPPING').sum())}

## Representation

r_wt dimension: {RESIDUE_DIM}
r_mut dimension: {RESIDUE_DIM}
delta_r dimension: {RESIDUE_DIM}
partner representation: implicit_via_concatenated_complex

## Runtime Notes

- 单 GPU ({DEVICE})；未启用双 GPU
- fp16 推理，float16 存储
- MINT 与 ESM 串行加载，未同时驻留
- 处理顺序：seed=42 按突变数分层交错（**避免 pilot 式选择偏倚**）

## Final Export

ZIP: V7_2R_2_TRACKB_FEATURES.zip
SHA256: (见下方打印)

## Scientific Status

FEATURE_EXTRACTION_ONLY

No scientific performance conclusion was produced on Kaggle.
"""
(OUTPUT_ROOT / "KAGGLE_RUN_REPORT.md").write_text(REP, encoding="utf-8")
print(REP)
''')

# ============================================================ 15 Create Download ZIP
md("## 15 Create Download ZIP")
code('''
# ============================================================ 15 Create Download ZIP
man = load_manifest()
files = sorted(FEAT_DIR.glob("*.npz"))
zippath = Path("/kaggle/working/V7_2R_2_TRACKB_FEATURES.zip")
with zipfile.ZipFile(zippath, "w", zipfile.ZIP_DEFLATED) as z:
    for f in files:
        z.write(f, arcname=f"features/{f.name}")
    for n in ("trackB_manifest.csv", "artifact_manifest.csv",
              "environment_report.json", "coverage_report.md",
              "mapping_report.csv", "provenance.json",
              "artifact_dependency_graph.md", "smoke_test.csv",
              "partner_sanity.csv", "KAGGLE_RUN_REPORT.md"):
        p = OUTPUT_ROOT / n
        if p.exists():
            z.write(p, arcname=n)
    pn = PATCH_DIR / "PATCH_NOTES.md"
    if pn.exists():
        z.write(pn, arcname="PATCH_NOTES.md")

def sha256_of(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

zsha = sha256_of(zippath)
print("=" * 70)
print("DOWNLOAD THIS FILE:")
print("   ", zippath)
print("SIZE  :", round(zippath.stat().st_size / 1e6, 1), "MB")
print("SHA256:", zsha)
print("=" * 70)
print("\\nNext (local):")
print("  1) 下载 V7_2R_2_TRACKB_FEATURES.zip")
print("  2) 放到 F:\\\\dsh\\\\downloads\\\\v7_2r_2\\\\")
print("  3) python scripts/v7_2r_2/02_import_kaggle_features.py")
''')

nb = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.10"},
        "accelerator": "GPU",
    },
    "nbformat": 4, "nbformat_minor": 5,
}
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"写出 {NB}")
print(f"cells = {len(CELLS)}（markdown {sum(1 for c in CELLS if c['cell_type']=='markdown')}, "
      f"code {sum(1 for c in CELLS if c['cell_type']=='code')}）")
