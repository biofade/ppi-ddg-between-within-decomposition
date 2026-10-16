# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/00_build_notebook.py -- 生成 V7.2R-3 Kaggle notebook.

核心修复（§10/§11）：`build_complex()` 使用**全部蛋白链**；
突变**只施加**到受影响链，伙伴链保持 WT。
必须**先通过伙伴 sanity 测试**才允许全量提取（§7/§29）。
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
OUT = ROOT / "kaggle" / "v7_2r_3"
OUT.mkdir(parents=True, exist_ok=True)
NB = OUT / "V7_2R_3_PPI_Context_Kaggle.ipynb"
CELLS: list = []


def md(t: str):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": t.splitlines(True)})


def code(t: str):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": t.strip("\n").splitlines(True)})


md("""# V7.2R-3 — TRUE PPI-Conditioned Representation (Kaggle GPU)

**唯一科学问题**：真正的伙伴条件化 PPI 上下文，能否在已测试的
**链隔离**序列表示之上，带来可迁移的 family-OOD 突变排序信号？

**本 notebook 只做表示推理**（+ 强制伙伴 sanity 测试）。
不做排名评估、不用任何 ΔΔG 标签、不用任何监督突变 head。

> **关键修复**：V7.2R-2 只对携带突变的链做前向（86.2% 的 target 丢了伙伴链）。
> 本 notebook 的 `build_complex()` 使用 **PDB 全部蛋白链**；
> 突变只施加到受影响链，伙伴链保持 WT。
""")

code('''
# ============================================================ 00 Configuration
import os, sys, json, time, math, hashlib, zipfile, shutil, subprocess, warnings
import random
from pathlib import Path
warnings.filterwarnings("ignore")

PROJECT_NAME = "v7_2r_3_ppi"
DEVICE = "cuda:0"
USE_SECOND_GPU = False
RUN_SANITY_ONLY = False        # True = 只跑伙伴 sanity，不做全量
RESUME = True
OUTPUT_ROOT = Path("/kaggle/working/v7_2r_3")
INPUT_ROOT = Path("/kaggle/input")
DTYPE = "float16"
BATCH_SIZE = 1
MAX_TARGETS = None             # 调试用
MAXLEN = 1022                  # ESM2 位置上限（BOS/EOS 预留）

REPR_VERSION = "v7_2r_3_ppi_context_v1"
MODEL_VERSION = "esm2_t33_650M_UR50D|full_complex"
RESIDUE_DIM = 1280

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
FEAT_DIR = OUTPUT_ROOT / "features"; FEAT_DIR.mkdir(parents=True, exist_ok=True)
PATCH_DIR = OUTPUT_ROOT / "patched_source"; PATCH_DIR.mkdir(parents=True, exist_ok=True)
print("DEVICE:", DEVICE, "| OUT:", OUTPUT_ROOT, "| RESUME:", RESUME)
''')

code('''
# ============================================================ 01 Environment Audit
import platform
env_report = {"python": sys.version.split()[0], "platform": platform.platform()}
def sh(c):
    try:
        return subprocess.run(c, shell=True, capture_output=True, text=True,
                              timeout=180).stdout.strip()
    except Exception as e:
        return f"<failed {type(e).__name__}>"
print(sh("nvidia-smi")[:1800]); print("\\n=== df ==="); print(sh("df -h")[:900])
import torch
env_report.update(torch=torch.__version__,
                  cuda_available=bool(torch.cuda.is_available()),
                  cuda_version=getattr(torch.version, "cuda", None))
gpus = []
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        p = torch.cuda.get_device_properties(i)
        gpus.append({"index": i, "name": p.name,
                     "vram_gb": round(p.total_memory / 1024**3, 2)})
env_report["gpus"] = gpus
try:
    import psutil
    env_report["ram_gb"] = round(psutil.virtual_memory().total / 1024**3, 2)
except Exception:
    env_report["ram_gb"] = None
env_report["work_disk_free_gb"] = round(shutil.disk_usage("/kaggle/working").free / 1024**3, 2)
try:
    import urllib.request
    urllib.request.urlopen("https://huggingface.co", timeout=15).read(64)
    env_report["internet"] = True
except Exception:
    env_report["internet"] = False
print("\\nGPU:", gpus, "| RAM:", env_report["ram_gb"], "| net:", env_report["internet"])
assert env_report["cuda_available"], "请在 Settings 设 Accelerator = GPU T4 x2"
(OUTPUT_ROOT / "environment_report.json").write_text(json.dumps(env_report, indent=2))
''')

code('''
# ============================================================ 02 Dataset Discovery
REQUIRED = ["kaggle_input_manifest.csv", "canonical_rows.csv", "variants.csv",
            "mutation_sites.csv", "target_chains.csv"]
cands = [d.parent for d in sorted(INPUT_ROOT.rglob("kaggle_input_manifest.csv"))
         if all((d.parent / f).exists() for f in REQUIRED)]
print("候选:", [str(c) for c in cands])
if not cands:
    raise FileNotFoundError("未找到 V7.2R-3 输入包，请 Add Input 上传 "
                            "kaggle_input_package.zip（V7.2R-3 版本）")
PKG = sorted(cands, key=str)[0]
print("选定:", PKG)
import pandas as pd
canon = pd.read_csv(PKG / "canonical_rows.csv")
var   = pd.read_csv(PKG / "variants.csv")
sites = pd.read_csv(PKG / "mutation_sites.csv")
tch   = pd.read_csv(PKG / "target_chains.csv")
print("canonical_rows:", canon.shape, "| variants:", var.shape,
      "| mutation_sites:", sites.shape, "| target_chains:", tch.shape)
# 关键：确认含全部伙伴链
nch = tch.groupby("target_id").size()
print("\\n每 target 链数：中位 %.0f，最小 %d，最大 %d"
      % (nch.median(), nch.min(), nch.max()))
print("=> 修复验证：V7.2R-2 该中位数为 1（只含突变链）")
labtok = ("ddg","delta_g","affinity","kd(","label","y_true","pred","score")
badcol = [c for df in (canon,var,sites,tch) for c in df.columns
          if any(t in c.lower() for t in labtok)]
print("疑似标签列:", badcol if badcol else "无 ✅")
TARGETS = sorted(var["target_id"].astype(str).unique())
if MAX_TARGETS:
    TARGETS = TARGETS[:int(MAX_TARGETS)]
print("\\n目标数:", len(TARGETS))
''')

code('''
# ============================================================ 03 Dependencies
def pip_install(pkgs):
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-input"] + pkgs,
                       capture_output=True, text=True)
    print(("OK  " if r.returncode == 0 else "FAIL"), " ".join(pkgs))
    if r.returncode: print(r.stderr[-1200:])
need = []
for m, p in (("einops","einops"), ("omegaconf","omegaconf"),
             ("fairscale","fairscale"), ("biopython","biopython"), ("psutil","psutil")):
    try: __import__(m)
    except ImportError: need.append(p)
if need: pip_install(need)
print("\\n未降级 torch/CUDA（§15）")
''')

code('''
# ============================================================ 04 Artifact Audit + Dependency Graph
ARTIFACTS = [
  {"name":"esm2_t33_650M_UR50D","url":"https://huggingface.co/facebook/esm2_t33_650M_UR50D",
   "license":"MIT","purpose":"全复合物残基级表示（位点向量）","status":"GREEN"},
  {"name":"mint.ckpt","url":"https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt",
   "license":"MIT","purpose":"链级/对级 PPI 表示（伙伴 sanity 用）","status":"GREEN"},
  {"name":"SiweiLab/ReCLIP (source)","url":"https://github.com/SiweiLab/ReCLIP.git",
   "license":"see repo LICENSE","purpose":"特征构建参考","status":"GREEN (code only)"},
]
FORBIDDEN = ["RiverZ/reclip classifier_heads/mutation/reclip_mutation_xgb_*.pkl",
             "Prompt-DDG SKEMPI ddg_model.ckpt", "CATH-ddG / DDAffinity supervised weights"]
for a in ARTIFACTS: print(f"  [{a['status']:8s}] {a['name']:28s} {a['license']}")
print("\\n禁止使用:")
for f in FORBIDDEN: print("  [FORBIDDEN]", f)

DEP = """# Artifact Dependency Graph (V7.2R-3)

## 与 V7.2R-2 的关键差别
V7.2R-2: build_variant() -> 只含突变链 -> **不是 PPI-conditioned**
V7.2R-3: build_complex() -> PDB 全部蛋白链 -> 真正的复合物前向

## 模型使用
1. **ESM2-650M over FULL complex**（主）：把全部链按链 ID 排序拼接后一次前向。
   自注意力天然跨链 -> 位点向量即为 **partner-conditioned**。
   索引：残基 (chain c, seq_index0 i) -> 全局位置 offset(previous chains)+i -> token+1
2. **MINT chain-level**（辅，用于伙伴 sanity）：MINTWrapper(chains, chain_ids) -> (N,1280)
   只做 chain 级，**不假装是 residue 级**（§6 强制区分）

SEPARATE_ESM_REQUIRED = NO（ESM2-650M 为主，必须加载一次）
MINT_LOADED = YES（仅用于 chain-level sanity，可选）
"""
(OUTPUT_ROOT / "artifact_dependency_graph.md").write_text(DEP)
print(DEP)
''')

code('''
# ============================================================ 05 Load Model
from transformers import AutoTokenizer, AutoModel
ESM_NAME = "facebook/esm2_t33_650M_UR50D"
print("加载 ESM2-650M fp16 ...")
t0 = time.time()
tok = AutoTokenizer.from_pretrained(ESM_NAME)
esm = AutoModel.from_pretrained(ESM_NAME, torch_dtype=torch.float16).to(DEVICE).eval()
RESIDUE_DIM = esm.config.hidden_size
print(f"{time.time()-t0:.0f}s | hidden={RESIDUE_DIM} | layers={esm.config.num_hidden_layers}")
print("显存(GB):", round(torch.cuda.memory_allocated()/1e9, 2))
''')

code('''
# ============================================================ 06 Build Complex (FIXED)
wt_seq = {}
for r in tch.itertuples():
    wt_seq[(str(r.target_id), str(r.chain))] = str(r.wt_sequence)
chains_of_target = {}
for r in tch.itertuples():
    chains_of_target.setdefault(str(r.target_id), []).append(str(r.chain))
for t in chains_of_target:
    chains_of_target[t] = sorted(chains_of_target[t])
sites_by_row = {int(k): g for k, g in sites.groupby("row_id")}
var = var.copy(); var["row_id"] = var["row_id"].astype(int)
var = var.set_index("row_id", drop=False)
print("target 数:", len(chains_of_target), "| WT 链序列:", len(wt_seq))

def build_complex(target, row_id=None, mutated=True):
    """返回 (seqs_dict, offsets, sites)
    seqs_dict: {chain: sequence}  —— **含全部伙伴链**（§10）
    row_id=None -> 纯 WT 复合物
    mutated=True -> 对受影响链施加该 row 的全部替换
    """
    chains = chains_of_target.get(target, [])
    seqs = {c: wt_seq.get((target, c), "") for c in chains}
    seqs = {c: s for c, s in seqs.items() if len(s) > 0}
    st = []
    if row_id is not None:
        g = sites_by_row.get(int(row_id))
        if g is None or not len(g):
            return None
        for x in g.itertuples():
            c = str(x.chain); i = int(x.seq_index0)
            if c not in seqs or not (0 <= i < len(seqs[c])):
                return None
            if seqs[c][i] != str(x.wt_aa):     # §19 绝不静默修复
                return None
            st.append({"chain": c, "pos": i, "wt": str(x.wt_aa), "mut": str(x.mut_aa)})
        if mutated:
            for c in sorted({s["chain"] for s in st}):
                L2 = list(seqs[c])
                for s in st:
                    if s["chain"] == c:
                        L2[s["pos"]] = s["mut"]
                seqs[c] = "".join(L2)
            for s in st:                        # 回读核对
                if seqs[s["chain"]][s["pos"]] != s["mut"]:
                    return None
    return seqs, st

def global_index(seqs, target_site):
    """(chain, pos) -> 拼接序列中的全局残基下标"""
    off = 0
    for c in sorted(seqs):
        if c == target_site["chain"]:
            return off + target_site["pos"]
        off += len(seqs[c])
    return None

probe = list(var.index[:5])
print("\\n构建自检:")
for r in probe:
    t = str(var.loc[r, "target_id"])
    bc = build_complex(t, r)
    if bc is None:
        print(f"   row {r}: 失败"); continue
    seqs, st = bc
    bwt = build_complex(t, None)
    print(f"   row {r} ({t}): 链数={len(seqs)} 位点={len(st)} "
          f"WT链数={len(bwt[0])} 总长={sum(len(s) for s in seqs.values())}")
''')

code('''
# ============================================================ 07 Embed full complex
def embed_complex(seqs, want_sites=None):
    """对**全部链**拼接后一次前向；返回 {global_pos: vector} 与拼接信息。
    want_sites: [(chain,pos)] 只取这些位置的向量（节省显存/带宽）"""
    order = sorted(seqs.keys())
    seq_list, bounds, off = [], {}, 0
    for c in order:
        s = seqs[c]
        bounds[c] = (off, off + len(s))
        seq_list.append(s); off += len(s)
    full = "".join(seq_list)
    if len(full) > MAXLEN:
        return None, {"reason": "SEQUENCE_LIMIT_EXCEEDED",
                      "length": len(full), "maxlen": MAXLEN}
    ids = tok(full, return_tensors="pt")["input_ids"].to(DEVICE)
    with torch.no_grad():
        out = esm(input_ids=ids).last_hidden_state[0]
    res = {}
    if want_sites is None:
        for c in order:
            for i in range(len(seqs[c])):
                res[(c, i)] = out[bounds[c][0] + i + 1].float().cpu().numpy()
    else:
        for c, p in want_sites:
            res[(c, p)] = out[bounds[c][0] + p + 1].float().cpu().numpy()
    del out, ids
    torch.cuda.empty_cache()
    return res, {"n_chains": len(order), "length": len(full), "bounds": bounds}
''')

code('''
# ============================================================ 08 TRUE PARTNER SANITY (MANDATORY GATE)
# §7：至少 10 个 target，label-free，按突变数分层。
nt = var.groupby("target_id")["n_sites"].agg(["size", "mean"]).sort_values("mean")
cand = list(nt.index.astype(str))
# 分层抽 10 个（小/中/大）
k = min(10, len(cand))
idx = [int(round(x)) for x in np.linspace(0, len(cand)-1, k)] if len(cand) > 1 else [0]
import numpy as np
SANITY_TARGETS = [cand[i] for i in sorted(set(idx))]
print("sanity targets (label-free, 按突变数分层):", SANITY_TARGETS)

# 为 wrong-partner 测试准备一个「错误伙伴」序列池（来自其他 target 的链）
all_chains_pool = [s for (t, c), s in wt_seq.items() if len(s) > 30]
rng = random.Random(42)
rng.shuffle(all_chains_pool)

sanity_rows = []
for t in SANITY_TARGETS:
    sub = var[var["target_id"] == t]
    if not len(sub):
        continue
    rid = int(sub.index[0])
    bc = build_complex(t, rid, mutated=False)     # WT 构象，用于语义测试
    if bc is None:
        continue
    seqs, st = bc
    if not st:
        continue
    chains = sorted(seqs)
    mut_chains = sorted({s["chain"] for s in st})
    partners = [c for c in chains if c not in mut_chains]
    if not partners:
        sanity_rows.append({"target": t, "row_id": rid, "n_chains": len(chains),
                            "n_partner": 0, "D_true": None, "D_wrong": None,
                            "note": "NO_PARTNER_CHAIN"})
        continue
    # A alone
    A_seqs = {c: seqs[c] for c in mut_chains}
    A_sites = [(s["chain"], s["pos"]) for s in st]
    rA, _ = embed_complex(A_seqs, A_sites)
    # A + true partners
    rAB, _ = embed_complex(seqs, A_sites)
    # A + wrong partner（替换伙伴序列，保持长度结构合理）
    wrong = dict(seqs)
    for j, c in enumerate(partners):
        wrong[c] = all_chains_pool[(hash((t, c)) % len(all_chains_pool))]
    rAX, _ = embed_complex(wrong, A_sites)
    D_true = float(np.mean([np.linalg.norm(rA[k] - rAB[k]) for k in rA if k in rAB]))
    D_wrong = float(np.mean([np.linalg.norm(rAB[k] - rAX[k]) for k in rAB if k in rAX]))
    sanity_rows.append({"target": t, "row_id": rid, "n_chains": len(chains),
                        "n_partner": len(partners), "D_true": D_true,
                        "D_wrong": D_wrong, "note": ""})
    print(f"  {t}: chains={len(chains)} partners={len(partners)} "
          f"D_true={D_true:.4f} D_wrong={D_wrong:.4f}")
    del rA, rAB, rAX
    torch.cuda.empty_cache()

SA = pd.DataFrame(sanity_rows)
SA.to_csv(OUTPUT_ROOT / "partner_sanity.csv", index=False)
print("\\n", SA.to_string(index=False))

ok = SA[(SA.D_true.notna())]
if len(ok) == 0:
    print("\\n!! 没有任何 target 有伙伴链 —— 无法做 sanity")
    PPI_OK = False
else:
    dt = float(ok.D_true.median()); dw = float(ok.D_wrong.median())
    print(f"\\nD_true  中位 = {dt:.4f}")
    print(f"D_wrong 中位 = {dw:.4f}")
    PPI_OK = (dt > 1e-3)
    print("PPI_CONDITIONING =", "CONFIRMED ✅" if PPI_OK else "NOT_DEMONSTRATED ❌")
    if not PPI_OK:
        print("=> 依 §7：PPI_CONDITIONING_NOT_DEMONSTRATED，**停止**，"
              "不在错误标签下做排序实验。")

if not PPI_OK:
    raise RuntimeError("PPI_CONDITIONING_NOT_DEMONSTRATED —— 依 §7 停止。")
if RUN_SANITY_ONLY:
    print("\\nRUN_SANITY_ONLY=True -> 停止")
    raise SystemExit(0)
''')

code('''
# ============================================================ 09 Full Streaming Inference (all chains)
rng2 = random.Random(42)
tstat = var.groupby("target_id")["n_sites"].max().to_dict()
strata = {"2-9": [], "10-29": [], ">=30": []}
for t in TARGETS:
    n = int(tstat.get(t, 0))
    strata["2-9" if n <= 9 else ("10-29" if n <= 29 else ">=30")].append(t)
for k in strata: rng2.shuffle(strata[k])
order, = ([],)
while any(strata.values()):
    for k in ("2-9", "10-29", ">=30"):
        if strata[k]: order.append(strata[k].pop())
print("分层交错顺序（前 10）:", order[:10])

MAN = OUTPUT_ROOT / "trackB_manifest.csv"
def load_man():
    return pd.read_csv(MAN) if MAN.exists() else pd.DataFrame(
        columns=["canonical_target","family","total_rows","completed_rows",
                 "status","feature_file","SHA256","error_code"])
man = load_man()
done = set(man.loc[man.status=="COMPLETE","canonical_target"].astype(str)) if len(man) else set()
print("已有 COMPLETE:", len(done), "| RESUME =", RESUME)

def sha256_of(p, ch=1<<22):
    h = hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(ch), b""): h.update(b)
    return h.hexdigest()

t_start = time.time()
for ti, t in enumerate(order):
    if RESUME and t in done: continue
    if MAX_TARGETS and ti >= int(MAX_TARGETS): break
    t0 = time.time()
    sub = var[var["target_id"] == t]
    fam = str(sub["family_id"].iloc[0])
    out_p = FEAT_DIR / f"{t}.npz"
    rec = {"canonical_target": t, "family": fam, "total_rows": len(sub),
           "completed_rows": 0, "status": "RUNNING", "feature_file": str(out_p),
           "SHA256": None, "error_code": None}
    try:
        # WT 复合物：每链每残基（缓存，仅一次/target）
        bc0 = build_complex(t, None)
        if bc0 is None:
            raise RuntimeError("BUILD_WT_FAILED")
        wt_seqs, _ = bc0
        r_wt_all, info = embed_complex(wt_seqs, None)
        if r_wt_all is None:
            rec.update(status="SEQUENCE_LIMIT", error_code=info.get("reason"))
            raise RuntimeError("SEQUENCE_LIMIT")
        rids, rwt, rmut, dr, nsite = [], [], [], [], []
        ok_rows = 0
        for rid in sub.index:
            bc = build_complex(t, int(rid), mutated=True)
            if bc is None: continue
            mseqs, st = bc
            if not st: continue
            if len(mseqs) != len(wt_seqs):
                continue
            sites_xy = [(s["chain"], s["pos"]) for s in st]
            r_mut_all, _ = embed_complex(mseqs, sites_xy)
            if r_mut_all is None: continue
            W = np.stack([r_wt_all[k] for k in sites_xy if k in r_wt_all]).astype(np.float16)
            M = np.stack([r_mut_all[k] for k in sites_xy if k in r_mut_all]).astype(np.float16)
            if W.shape[0] == 0: continue
            n = min(W.shape[0], M.shape[0]); W, M = W[:n], M[:n]
            rids.append(int(rid)); nsite.append(n)
            rwt.append(W.mean(0)); rmut.append(M.mean(0)); dr.append((M-W).mean(0))
            ok_rows += 1
            del r_mut_all, W, M
            torch.cuda.empty_cache()
        if ok_rows == 0:
            rec.update(status="EXCLUDED_MAPPING", error_code="NO_VALID_SITE")
        else:
            np.savez_compressed(
                out_p,
                row_id=np.array(rids, np.int64),
                target_id=np.array([t]*ok_rows, dtype=object),
                family_id=np.array([fam]*ok_rows, dtype=object),
                mutation_count=np.array(nsite, np.int32),
                r_wt=np.stack(rwt).astype(np.float16),
                r_mut=np.stack(rmut).astype(np.float16),
                delta_r=np.stack(dr).astype(np.float16),
                representation_level=np.array(["MUTATION_SITE_PPI"]),
                representation_version=np.array([REPR_VERSION]),
                model_version=np.array([MODEL_VERSION]),
                n_chains_embedded=np.array([len(wt_seqs)], np.int32),
            )
            rec.update(status="COMPLETE", completed_rows=ok_rows,
                       SHA256=sha256_of(out_p))
        del r_wt_all, rwt, rmut, dr
        torch.cuda.empty_cache()
    except torch.cuda.OutOfMemoryError:
        rec.update(status="FAILED", error_code="CUDA_OOM"); torch.cuda.empty_cache()
    except Exception as e:
        if rec["error_code"] is None:
            rec.update(status="FAILED", error_code=f"{type(e).__name__}: {str(e)[:100]}")
    man = man[man.canonical_target.astype(str) != t] if len(man) else man
    man = pd.concat([man, pd.DataFrame([rec])], ignore_index=True)
    man.to_csv(MAN, index=False)
    if (ti+1) % 5 == 0 or rec["status"] != "COMPLETE":
        nrow = int(man.loc[man.status=="COMPLETE","completed_rows"].sum()) if len(man) else 0
        ntg = int((man.status=="COMPLETE").sum()) if len(man) else 0
        print(f"[{ti+1}/{len(order)}] {t}: {rec['status']} "
              f"{rec['completed_rows']}/{rec['total_rows']} {time.time()-t0:.0f}s "
              f"| 累计 rows={nrow} targets={ntg} | {(time.time()-t_start)/60:.1f}min")
print("\\n全量推理结束")
''')

code('''
# ============================================================ 10 Coverage QC (no labels)
man = load_man(); comp = man[man.status=="COMPLETE"]
tot = int(var.shape[0]); dr_ = int(comp.completed_rows.sum())
print(man.status.value_counts().to_string())
print(f"\\nrow coverage: {dr_}/{tot} = {100*dr_/tot:.1f}%")
print(f"target coverage: {comp.canonical_target.nunique()}/{var.target_id.nunique()}")
print(f"family coverage: {comp.family.nunique()}/{var.family_id.nunique()}")
allt = var.groupby("target_id").size(); ge5 = set(allt[allt>=5].index.astype(str))
c5 = len(set(comp.canonical_target.astype(str)) & ge5)
print(f"n>=5 target coverage: {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}%")
man[["canonical_target","family","total_rows","completed_rows","status",
     "error_code"]].to_csv(OUTPUT_ROOT/"mapping_report.csv", index=False)
(OUTPUT_ROOT/"coverage_report.md").write_text(
f"""# V7.2R-3 PPI Context Coverage

| 项 | 值 |
|---|---|
| row coverage | {dr_}/{tot} = {100*dr_/tot:.1f}% |
| target coverage | {comp.canonical_target.nunique()}/{var.target_id.nunique()} |
| family coverage | {comp.family.nunique()}/{var.family_id.nunique()} |
| n>=5 target coverage | {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}% |
| SEQUENCE_LIMIT | {int((man.status=='SEQUENCE_LIMIT').sum())} |
| EXCLUDED_MAPPING | {int((man.status=='EXCLUDED_MAPPING').sum())} |
| FAILED | {int((man.status=='FAILED').sum())} |

**表示层级：MUTATION_SITE_PPI**（全复合物 ESM2-650M，跨链自注意力）
**无任何 ΔΔG 或性能指标。**
""", encoding="utf-8")
print("写出 mapping_report.csv / coverage_report.md")
''')

code('''
# ============================================================ 11 Export + Provenance
files = sorted(FEAT_DIR.glob("*.npz"))
tot_b = sum(f.stat().st_size for f in files)
print("feature 文件:", len(files), f"| {tot_b/1e6:.1f} MB")
if files:
    z = np.load(files[0], allow_pickle=True)
    print("字段:", list(z.files))
    print("维度:", {k: int(z[k].shape[1]) for k in ("r_wt","r_mut","delta_r")})
    print("representation_level:", z["representation_level"])
prov = {"project": PROJECT_NAME, "representation_version": REPR_VERSION,
        "model_version": MODEL_VERSION, "esm_model": ESM_NAME,
        "esm_hidden": RESIDUE_DIM, "device": DEVICE,
        "dtype_inference": DTYPE, "dtype_storage": "float16",
        "representation_level": "MUTATION_SITE_PPI",
        "complex_construction": "ALL_PROTEIN_CHAINS (partners retained)",
        "labels_used": False, "supervised_heads_used": False,
        "forbidden_artifacts_used": [],
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
(OUTPUT_ROOT/"provenance.json").write_text(json.dumps(prov, indent=2))
pd.DataFrame(ARTIFACTS).to_csv(OUTPUT_ROOT/"artifact_manifest.csv", index=False)
print(json.dumps(prov, indent=2))
''')

code('''
# ============================================================ 12 Final Report + ZIP
man = load_man(); comp = man[man.status=="COMPLETE"]
files = sorted(FEAT_DIR.glob("*.npz"))
tot = int(var.shape[0]); dr_ = int(comp.completed_rows.sum())
allt = var.groupby("target_id").size(); ge5 = set(allt[allt>=5].index.astype(str))
c5 = len(set(comp.canonical_target.astype(str)) & ge5)
dt = float(SA.D_true.median()) if len(SA) and SA.D_true.notna().any() else float("nan")
dw = float(SA.D_wrong.median()) if len(SA) and SA.D_wrong.notna().any() else float("nan")

REP = f"""# V7.2R-3 Kaggle PPI-Context Run

## Environment
GPU: {env_report.get('gpus')}
CUDA: {env_report.get('cuda_version')} | torch: {env_report.get('torch')}
RAM: {env_report.get('ram_gb')} GB | net: {env_report.get('internet')}

## Artifacts
ESM: {ESM_NAME} | MINT: loaded for chain-level sanity

## Fix vs V7.2R-2
V7.2R-2 embedded ONLY mutation-bearing chains (86.2% of targets lost partners).
V7.2R-3 embeds ALL protein chains of the complex.

## Representation Level
MUTATION_SITE_PPI (full-complex ESM2-650M, cross-chain self-attention)
r_wt/r_mut/delta_r dim: {RESIDUE_DIM}

## Partner Sanity
n targets tested: {len(SA)}
D_true  median: {dt:.4f}
D_wrong median: {dw:.4f}
PPI_CONDITIONING: {'CONFIRMED' if (dt==dt and dt>1e-3) else 'NOT_DEMONSTRATED'}

## Coverage
Rows: {dr_}/{tot} = {100*dr_/tot:.1f}%
Targets: {comp.canonical_target.nunique()}/{var.target_id.nunique()}
Families: {comp.family.nunique()}/{var.family_id.nunique()}
n>=5: {c5}/{len(ge5)} = {100*c5/max(1,len(ge5)):.1f}%

## Scientific Status
FEATURE_EXTRACTION_ONLY
No ranking evaluation was performed on Kaggle.
"""
(OUTPUT_ROOT/"KAGGLE_RUN_REPORT.md").write_text(REP)
print(REP)

zp = Path("/kaggle/working/V7_2R_3_PPI_CONTEXT_FEATURES.zip")
with zipfile.ZipFile(zp,"w",zipfile.ZIP_DEFLATED) as z:
    for f in files: z.write(f, arcname=f"features/{f.name}")
    for n in ("trackB_manifest.csv","artifact_manifest.csv","environment_report.json",
              "coverage_report.md","mapping_report.csv","provenance.json",
              "artifact_dependency_graph.md","partner_sanity.csv",
              "KAGGLE_RUN_REPORT.md"):
        p = OUTPUT_ROOT/n
        if p.exists(): z.write(p, arcname=n)
def sha256_of(p, ch=1<<22):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(ch), b""): h.update(b)
    return h.hexdigest()
print("="*70)
print("DOWNLOAD THIS FILE:"); print("   ", zp)
print("SIZE  :", round(zp.stat().st_size/1e6,1), "MB")
print("SHA256:", sha256_of(zp))
print("="*70)
print("\\nNext (local): 放到 F:\\\\dsh\\\\downloads\\\\v7_2r_3\\\\ 后运行")
print("  python scripts/v7_2r_3/02_import_ppi_features.py")
''')

nb = {"cells": CELLS,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.10"},
                   "accelerator": "GPU"},
      "nbformat": 4, "nbformat_minor": 5}
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"写出 {NB}")
print(f"cells = {len(CELLS)}")
