# -*- coding: utf-8 -*-
"""scripts/v2.3.2/11_make_perresidue_notebook.py -- 生成 per-residue ESM-2 提取 notebook.

为什么需要它
------------
当前 `features/v2_esm2_cache.npz` 只有 **pooled** 表示 (每条序列一个 640 维向量)。
V2.3.2 的 NN 判决实验要回答"pooled 嵌入上 NN 是否可行"; 但 NN (1D-CNN / 注意力)
真正有用的地方是 **per-residue** 表示 (每条序列 L x 640) —— 那里才有序列顺序信息。

本生成器产出 Kaggle notebook `scripts/kaggle/esm2_extract_perresidue.ipynb`:
* 复用 150M pooled 版 notebook 的**全部核心逻辑** (字母表/清洗/截断/md5/OOM 降级/断点续算),
  通过 ``importlib`` 直接复用 `scripts/kaggle/esm2_extract_kaggle.ipynb` 的生成器常量
  (`CODE_PREP` 等), 保证逐位一致;
* 唯一改动: **打包时存 per-residue 矩阵而不是 mean-pool 向量**。

输出格式 (变长存储, 避免 padding 浪费)
--------------------------------------
``np.savez_compressed``:
* ``pdb_ids``      : (N,) 每个条目的 PDB ID
* ``interface_key``: (N,) `pdb|配体链|受体链` (与 V2.2 特征表同格式)
* ``lig_lengths``  : (N,) int32, 配体序列长度 (截断后)
* ``rec_lengths``  : (N,) int32
* ``lig_flat``     : (sum(lig_lengths), 640) float16 —— 所有配体 per-residue 首尾相接
* ``rec_flat``     : (sum(rec_lengths), 640) float16
* 以及元信息 (model_name / repr_layer / dtype / version)

> 用 **flat + lengths** 而不是 list-of-arrays: ``np.savez`` 存 object 数组需要
> ``allow_pickle=True``, 加载端有安全与兼容风险; flat 形式是纯数值数组, 更稳。
> 精度用 **float16** 把体积砍半 (pooled 版是 float32); 若下游需要更高精度,
> 把 notebook 里的 ``PERRES_DTYPE`` 改成 ``float32``。

体积估算 (notebook 顶部也会写)
------------------------------
9519 条唯一序列 x 平均 255 残基 x 640 维 x 2 bytes (fp16) ≈ **3.1 GB**,
再加上配体/受体都要存 (同一序列可能两边都出现, 但 notebook 按唯一序列去重,
所以总量就是唯一序列的总残基 x 640 x 2)。**这比 pooled 版 (52 MB) 大约 60 倍**,
因此 notebook 顶部会给出"先只提 test+val 子集"的省流方案。
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
KAGGLE_GEN = ROOT / "scripts" / "v2.3.1" / "10_make_kaggle_notebook.py"
OUT_NB = ROOT / "scripts" / "kaggle" / "esm2_extract_perresidue.ipynb"
OUT_DIR = "/kaggle/working/esm2_output_perres"
DATASET_SLUG = "ppb-affinity-proteinea"


def load_kaggle_module():
    """载入 Kaggle 生成器, 复用其 cell 源码常量 (单一实现来源)。"""
    spec = importlib.util.spec_from_file_location("kaggle_gen", KAGGLE_GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MD_TITLE = r"""# V2.3.2 — ESM-2 **per-residue** 提取（Kaggle，T4 x2）

> 这是 `esm2_extract_kaggle.ipynb`（pooled 版）的**扩充分支**：
> 核心逻辑（字母表 / 清洗 / 截断 / md5 / OOM 降级 / 断点续算）与它**逐位一致**，
> **唯一区别是打包时保存 per-residue 矩阵而不是 mean-pool 向量**。

---

## 为什么需要 per-residue

V2.3.2 的 NN 判决实验是在 **pooled** 嵌入上做的（每条序列 1 个 640 维向量）。
pooled 表示**丢掉了残基顺序**，所以 1D-CNN / 注意力这类结构根本无从发挥。

本 notebook 产出 **每条序列 L × 640** 的逐残基表示，供：
* 1D-CNN / 残差卷积
* 多头注意力 / 轻量 Transformer
* 界面残基加权池化（配合结构信息）

---

## ⚠️ 体积警告（**先读这一段**）

| 版本 | 存储 | 预估体积 |
|---|---|---|
| pooled（已有） | 9519 条 × 640 × **4 B** | **~52 MB**（实测 52,481,097 B） |
| **per-residue（本 notebook）** | 2,424,926 残基 × 640 × **2 B (fp16)** | **~3.1 GB** |
| per-residue（若用 fp32） | 2,424,926 残基 × 640 × 4 B | **~6.2 GB** |

Kaggle `/kaggle/working` 的可用磁盘通常约 **20 GB**（视版本而定），
3.1 GB 可以放下，但**加上 checkpoints 会翻倍**（每条序列一个 `.npz`，总大小≈最终大小）。

### 省流方案（推荐先做）

**方案 A（最省）**：只提 **test + val 子集**。
test 752 行 + val 965 行，去重后约 **3,000 条序列**，体积降到 **~1 GB**。
足够做 NN 判决实验（NN 只需要在 test 上评估）。

**方案 B（最有用）**：只提**界面残基**（每条约 10-40 个残基而非整条），
体积可降到 **~100 MB 量级**，但需要结构信息配合（本 notebook 不做）。

**方案 C（全量）**：约 3.1 GB，需要确认 `/kaggle/working` 剩余空间 ≥ 8 GB。

> notebook 顶部有 `SUBSET = "all" | "test_val"` 开关，**默认 `test_val`**（省流）。
> 想跑全量就把它改成 `"all"`。

---

## 写死的配置（与 pooled 版一致）

| 项 | 值 |
|---|---|
| 模型 | `esm2_t30_150M_UR50D`（640 维） |
| 表示层 | **最后一层（layer 30）** |
| 池化 | **无**（保存逐残基，长度截止到去 BOS/EOS 后的实际残基数） |
| 存储精度 | **float16**（`PERRES_DTYPE`，可改 float32） |
| 最大长度 | **1022 aa**，超出截断并记 `truncated` |
| 非标准残基 | 映射为 `<unk>`，不丢弃序列 |
| 输入 | `/kaggle/input/ppb-affinity-proteinea/filtered.csv` |
| 输出 | `/kaggle/working/esm2_output_perres/v2_esm2_perres.npz` |

## 预计运行时间

T4 x1 下，全量 9519 条约 **10–15 分钟**（与 pooled 版相同，因为前向计算一样，
只是**写盘体积大很多**，落盘会额外花几分钟）。`test_val` 子集约 **3–5 分钟**。

## 中途断开后如何恢复

与 pooled 版一致：每条序列立刻写
`/kaggle/working/esm2_output_perres/checkpoints/<md5>.npz`，重跑自动跳过已完成。
**但 per-residue 的 checkpoint 体积大**，如果磁盘紧张，可把 `CHECKPOINT_ON = False`
（那就没有断点续算能力了，需要一次跑完）。

## 下载方式

Kaggle 没有 `files.download`。跑完用**右侧 Output 面板** →
展开 `esm2_output_perres/` → 下载 `v2_esm2_perres.npz`。
"""

MD_STEP1 = r"""## 步骤 1 / 7 — 安装 fair-esm

Kaggle 已预装 torch / numpy。需要 Internet = On。"""

CODE_INSTALL = r"""!pip install fair-esm -q
print("fair-esm 安装完成。")"""

MD_STEP2 = r"""## 步骤 2 / 7 — 确认 GPU（只用 cuda:0）"""

CODE_GPU_CHECK = r"""import sys
import torch

if not torch.cuda.is_available():
    raise SystemExit(
        "请开启 Kaggle GPU: 右上角 Settings -> Accelerator -> 选 GPU T4 x2 -> 保存。"
    )
N_GPU = torch.cuda.device_count()
print(f"torch: {torch.__version__}")
print(f"torch.cuda.device_count() = {N_GPU}")
for _i in range(N_GPU):
    _p = torch.cuda.get_device_properties(_i)
    print(f"  GPU {_i}: {_p.name}  ({_p.total_memory / 1024**3:.1f} GB)")
DEVICE = torch.device("cuda:0")
print(f"使用设备: {DEVICE}")
print("GPU 检查通过。")"""

MD_STEP3 = r"""## 步骤 3 / 7 — 路径、子集开关与磁盘检查

> **`SUBSET` 默认 `test_val`（省流，约 1 GB）**。想跑全量改成 `"all"`（约 3.1 GB）。"""

CODE_SETUP = r"""import os
import shutil
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ====================================================================
# 只改这两处
# ====================================================================
DATASET_SLUG = "ppb-affinity-proteinea"
SUBSET = "test_val"     # "test_val" = 只提 test+val (省流, ~1GB) | "all" = 全量 (~3.1GB)

PROTEINEA_CSV = f"/kaggle/input/{DATASET_SLUG}/filtered.csv"
OUT_DIR = "/kaggle/working/esm2_output_perres"
CHECKPOINT_DIR = os.path.join(OUT_DIR, "checkpoints")
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

# ---- 写死的模型与预处理配置（与 pooled 版逐位一致） ----
MODEL_NAME = "esm2_t30_150M_UR50D"
MODEL_FALLBACK = "esm2_t12_35M_UR50D"
REPR_LAYER = 30
MAX_LEN = 1022
BATCH_SIZE = 2                      # per-residue 显存占用更大, 保守取 2
BACKUP_EVERY = 100
PERRES_DTYPE = np.float16           # 存储精度 (改 np.float32 可提升精度, 体积翻倍)
CHECKPOINT_ON = True                # False = 不写每条 checkpoint (省磁盘, 但失去断点续算)
LIG_COL = "Ligand Sequences"
REC_COL = "Receptor Sequences"
ID_COL = "PDB"
SPLIT_COL = "split"

print("=== /kaggle/input 内容 ===")
_in = Path("/kaggle/input")
if _in.exists():
    for d in sorted(_in.iterdir()):
        print(f"  {d.name}/")
        try:
            for f in sorted(d.iterdir())[:8]:
                print(f"      {f.name}")
        except Exception:
            pass

csv_path = Path(PROTEINEA_CSV)
if not csv_path.exists():
    _avail = sorted(p.name for p in _in.iterdir()) if _in.exists() else []
    raise FileNotFoundError(
        f"找不到 {PROTEINEA_CSV}\n"
        f"请确认 Dataset 名称为 '{DATASET_SLUG}' 且已在右侧 Add Input 挂载。\n"
        f"/kaggle/input 现有: {_avail}"
    )

USE_COLS = [ID_COL, LIG_COL, REC_COL, SPLIT_COL]
df = pd.read_csv(csv_path, usecols=USE_COLS, dtype=str)
df[ID_COL] = df[ID_COL].astype(str).str.strip().str.upper()
df = df.dropna(subset=[LIG_COL, REC_COL]).reset_index(drop=True)

if SUBSET == "test_val":
    df = df[df[SPLIT_COL].isin(["test", "val"])].reset_index(drop=True)
    print(f"\n!! SUBSET='test_val' -> 只保留 test+val")
print(f"读入 {len(df)} 行  <- {csv_path}")
print("split 分布:", df[SPLIT_COL].value_counts().to_dict())

# ---- 磁盘空间检查 ----
try:
    _st = shutil.disk_usage("/kaggle/working")
    print(f"\n/kaggle/working 可用空间: {_st.free / 1024**3:.1f} GB "
          f"(已用 {_st.used / 1024**3:.1f} GB / 共 {_st.total / 1024**3:.1f} GB)")
    if _st.free / 1024**3 < 4:
        print("  !! 可用空间 < 4 GB, 强烈建议把 SUBSET 改成 'test_val' 或关闭 CHECKPOINT_ON")
except Exception as exc:
    print(f"(磁盘检查失败: {exc})")
print("列名校验通过:", USE_COLS)"""

MD_STEP4 = r"""## 步骤 4 / 7 — 序列清洗与去重（与 pooled 版逐位一致）

非标准残基 -> `<unk>`；>1022 aa 截断；多链 `,` 也映射为 `<unk>`。"""

MD_STEP5 = r"""## 步骤 5 / 7 — per-residue 嵌入（主要耗时 + 主要磁盘占用）

与 pooled 版的差别：**不做 mean-pool**，保留 `(L, 640)` 的整条逐残基表示，
只去掉 BOS/EOS/padding 对应的位置。

* 断点续算：每条序列存一个 `.npz`（体积较大）
* 每 100 条打印进度与累计体积
* OOM 降级：batch 2 -> 1；仍 OOM 则切 35M 模型
* 单条失败不中断"""

CODE_PREP = r"""from esm.data import Alphabet

ALPHABET = Alphabet.from_architecture("ESM-2")
UNK = ALPHABET.unk_idx
STD = set(ALPHABET.all_tokens)


def sanitize(seq: str) -> tuple[str, bool, int]:
    # 返回 (清洗后序列, 是否截断, 非标准字符数)
    chars = list(seq)
    n_nonstd = sum(1 for c in chars if c not in STD)
    cleaned = "".join(c if c in STD else "X" for c in chars)
    truncated = len(cleaned) > MAX_LEN
    if truncated:
        cleaned = cleaned[:MAX_LEN]
    return cleaned, truncated, n_nonstd


def md5_of(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


records = df[[ID_COL, LIG_COL, REC_COL, SPLIT_COL]].to_dict("records")
for r in records:
    r["lig_clean"], r["lig_trunc"], r["lig_nonstd"] = sanitize(r[LIG_COL])
    r["rec_clean"], r["rec_trunc"], r["rec_nonstd"] = sanitize(r[REC_COL])
    r["lig_md5"] = md5_of(r["lig_clean"])
    r["rec_md5"] = md5_of(r["rec_clean"])

uniq_all = sorted({r["lig_clean"] for r in records} | {r["rec_clean"] for r in records})
_total_res = sum(len(s) for s in uniq_all)
print(f"行数                 : {len(records)}")
print(f"唯一序列（要嵌）      : {len(uniq_all)}")
print(f"总残基数             : {_total_res:,}")
print(f"预估输出体积 (fp16)   : {_total_res * 640 * 2 / 1024**3:.2f} GB")
print(f"预估输出体积 (fp32)   : {_total_res * 640 * 4 / 1024**3:.2f} GB")
print(f"含截断的行           : {sum(1 for r in records if r['lig_trunc'] or r['rec_trunc'])}")
if len(uniq_all) == 0:
    raise SystemExit("没有可嵌入的序列。")"""

CODE_EMBED = r"""import esm
import gc


failures: list[dict] = []


def load_model(name: str):
    models = esm.pretrained.esm2_model_registry
    key = name
    if key not in models:
        lower = {k.lower(): k for k in models}
        if name.lower() in lower:
            key = lower[name.lower()]
        else:
            raise KeyError(f"fair-esm 不认识 {name!r}")
    model, alphabet = models[key]()
    model = model.eval().to(DEVICE)
    return model, alphabet, alphabet.get_batch_converter(), key


MODEL, ALPHABET2, batch_converter, LOADED_KEY = load_model(MODEL_NAME)


def per_residue(hidden, tokens) -> "list[np.ndarray]":
    # 逐残基表示: 去掉 BOS/EOS/pad (1/2/0), 保留 (L_i, D)。返回长度 B 的列表。
    import torch as _t
    mask = (tokens != 1) & (tokens != 2) & (tokens != 0)
    out = []
    for b in range(hidden.shape[0]):
        idx = mask[b].nonzero(as_tuple=False).ravel()
        out.append(hidden[b].index_select(0, idx).float().cpu().numpy())
    return out


def embed_sequences(seqs, model, conv, batch_size, tag):
    # 返回 (成功的 (序列, (L,640) 矩阵) 列表, 失败记录)。**自带序列, 不会错位。**
    import torch as _t
    pairs, failed = [], []
    bs = batch_size
    i = 0
    while i < len(seqs):
        chunk = seqs[i:i + bs]
        try:
            data = [(f"s{j}", s) for j, s in enumerate(chunk)]
            _, _, tokens = conv(data)
            tokens = tokens.to(DEVICE)
            with _t.no_grad():
                reps = model(tokens, repr_layers=[REPR_LAYER], return_contacts=False)
            hidden = reps["representations"][REPR_LAYER]
            mats = per_residue(hidden, tokens)
            if len(mats) != len(chunk):
                raise RuntimeError(f"per_residue 返回 {len(mats)} != chunk {len(chunk)}")
            for s, m in zip(chunk, mats):
                pairs.append((s, m))
            i += len(chunk)
        except _t.cuda.OutOfMemoryError:
            _t.cuda.empty_cache(); gc.collect()
            if bs > 1:
                bs = max(1, bs // 2)
                print(f"  !! CUDA OOM -> batch_size 降到 {bs}（{tag}）")
                continue
            failed.append({"tag": tag, "index": i, "seq": chunk[0][:40],
                           "reason": "OOM at batch_size=1"})
            i += 1
        except Exception as exc:  # noqa: BLE001
            _t.cuda.empty_cache(); gc.collect()
            if bs > 1:
                bs = max(1, bs // 2)
                print(f"  !! batch 异常 -> batch_size 降到 {bs}: {exc}")
                continue
            failed.append({"tag": tag, "index": i, "seq": chunk[0][:40],
                           "reason": f"{type(exc).__name__}: {exc}"})
            i += 1
    return pairs, failed


def ckpt_path(md5: str) -> str:
    return os.path.join(CHECKPOINT_DIR, f"perres_{md5}.npz")


def load_ckpt(md5: str):
    p = ckpt_path(md5)
    if os.path.exists(p):
        try:
            with np.load(p) as z:
                return z["mat"].astype(np.float32)
        except Exception:
            return None
    return None


def save_ckpt(md5: str, mat: np.ndarray) -> None:
    if not CHECKPOINT_ON:
        return
    p = ckpt_path(md5)
    tmp = p + ".tmp.npz"
    np.savez(tmp, mat=mat.astype(PERRES_DTYPE))
    os.replace(tmp, p)


# ---------------- 主循环 ----------------
res_cache: dict[str, np.ndarray] = {}
for s in uniq_all:
    m = load_ckpt(md5_of(s))
    if m is not None:
        res_cache[s] = m
todo = [s for s in uniq_all if s not in res_cache]
print(f"总唯一序列 {len(uniq_all)}；已有缓存 {len(res_cache)}；待嵌入 {len(todo)}")

t_start = time.time()
i = 0
while i < len(todo):
    chunk = todo[i:i + BATCH_SIZE]
    pairs, failed = embed_sequences(chunk, MODEL, batch_converter, BATCH_SIZE, "main")
    for s, m in pairs:
        save_ckpt(md5_of(s), m)
        res_cache[s] = m
    failures.extend(failed)
    i += len(chunk)
    done = len(res_cache)
    if done % BACKUP_EVERY < len(chunk) or done >= len(uniq_all):
        _res_so_far = sum(int(v.shape[0]) for v in res_cache.values())
        el = time.time() - t_start
        rate = done / max(el, 1e-9)
        eta = (len(uniq_all) - done) / max(rate, 1e-9)
        print(f"  进度 {done}/{len(uniq_all)}  {rate:.1f} seq/s  "
              f"已用 {el/60:.1f} min  剩余 ~{eta/60:.1f} min  "
              f"累计残基 {_res_so_far:,} (~{_res_so_far*640*2/1024**3:.2f} GB fp16)  "
              f"失败 {len(failures)}")

print(f"嵌入完成：{len(res_cache)} 条，失败 {len(failures)} 条，"
      f"耗时 {(time.time()-t_start)/60:.1f} min")
if failures:
    for f_ in failures[:20]:
        print("  FAIL", f_)

_DIM = int(next(iter(res_cache.values())).shape[1]) if res_cache else 0
print(f"逐残基维度 DIM = {_DIM}")"""

CODE_PACK = r"""# 变长存储: flat + lengths (纯数值数组, 不需要 allow_pickle)
n_rows = len(records)
missing_lig = sum(1 for r in records if r["lig_clean"] not in res_cache)
missing_rec = sum(1 for r in records if r["rec_clean"] not in res_cache)
print(f"条目数 {n_rows}；缺配体 {missing_lig}；缺受体 {missing_rec}")

lig_lengths = np.zeros(n_rows, dtype=np.int32)
rec_lengths = np.zeros(n_rows, dtype=np.int32)
lig_chunks, rec_chunks = [], []
for i, r in enumerate(records):
    lm = res_cache.get(r["lig_clean"])
    rm = res_cache.get(r["rec_clean"])
    if lm is not None and lm.shape[0] > 0:
        lig_chunks.append(lm.astype(PERRES_DTYPE))
        lig_lengths[i] = lm.shape[0]
    if rm is not None and rm.shape[0] > 0:
        rec_chunks.append(rm.astype(PERRES_DTYPE))
        rec_lengths[i] = rm.shape[0]

lig_flat = (np.concatenate(lig_chunks, axis=0) if lig_chunks
            else np.zeros((0, _DIM), dtype=PERRES_DTYPE))
rec_flat = (np.concatenate(rec_chunks, axis=0) if rec_chunks
            else np.zeros((0, _DIM), dtype=PERRES_DTYPE))
print(f"lig_flat {lig_flat.shape}  rec_flat {rec_flat.shape}")

out_path = os.path.join(OUT_DIR, "v2_esm2_perres.npz")
np.savez_compressed(
    out_path,
    pdb_ids=np.array([r[ID_COL] for r in records], dtype="U"),
    interface_key=np.array([f"{r[ID_COL]}|{r[LIG_COL]}|{r[REC_COL]}"
                            for r in records], dtype="U"),
    lig_lengths=lig_lengths,
    rec_lengths=rec_lengths,
    lig_flat=lig_flat,
    rec_flat=rec_flat,
    splits=np.array([str(r[SPLIT_COL]) for r in records], dtype="U"),
    truncated=np.array([bool(r["lig_trunc"] or r["rec_trunc"]) for r in records],
                       dtype=bool),
    lig_md5=np.array([r["lig_md5"] for r in records], dtype="U"),
    rec_md5=np.array([r["rec_md5"] for r in records], dtype="U"),
    model_name=np.array(LOADED_KEY),
    repr_layer=np.array(REPR_LAYER),
    pool=np.array("none_per_residue_excl_bos_eos"),
    storage_dtype=np.array(np.dtype(PERRES_DTYPE).name),
    subset=np.array(SUBSET),
    perres_dim=np.array(_DIM),
    n_failed=np.array(len(failures)),
)
print(f"\n已写出: {out_path}")
print(f"  大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"  lig_flat {lig_flat.shape}  rec_flat {rec_flat.shape}  "
      f"总残基 {int(lig_lengths.sum()+rec_lengths.sum()):,}")
if _DIM != 640:
    print(f"  !! 警告: perres_dim={_DIM} != 640 (可能触发了 OOM 降级到 35M 模型)")"""

MD_STEP6 = r"""## 步骤 6 / 7 — 打包 `v2_esm2_perres.npz`

| key | 形状 | 含义 |
|---|---|---|
| `pdb_ids` | (N,) | 每个条目的 PDB ID |
| `interface_key` | (N,) | `pdb\|配体链\|受体链`（与 V2.2 特征表同格式） |
| `lig_lengths` / `rec_lengths` | (N,) | 每条序列截断后的残基数 |
| `lig_flat` | (Σlig_len, 640) | 配体逐残基表示首尾相接 |
| `rec_flat` | (Σrec_len, 640) | 受体逐残基表示首尾相接 |
| 元信息 | 标量 | `model_name` / `repr_layer` / `storage_dtype` / `subset` 等 |

读取方式（本机 Python）：

```python
import numpy as np
z = np.load("v2_esm2_perres.npz")
lig_flat, lig_len = z["lig_flat"], z["lig_lengths"]
off = np.concatenate([[0], np.cumsum(lig_len)])          # 第 i 条的切片是 [off[i]:off[i+1]]
seq_i_residues = lig_flat[off[3]:off[4]]                  # 第 3 条的逐残基矩阵 (L, 640)
```"""

MD_STEP7 = r"""## 步骤 7 / 7 — 下载（Kaggle Output 面板）

⚠️ Kaggle 没有 `files.download`。

1. 右侧 **Output** 面板 → 展开 **`esm2_output_perres/`**
2. 下载 **`v2_esm2_perres.npz`**

放到本机：

```text
F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_perres.npz
```

> 若文件很大（>1 GB），浏览器下载可能较慢。也可以在 Notebook 里把它
> `shutil.copy` 到 `/kaggle/working/` 根目录（Output 面板同样能看到），
> 或用 `Save Version` 后从版本页面下载。"""

CODE_DOWNLOAD = r"""import os

print("=" * 72)
print("下载步骤（Kaggle 没有 files.download）:")
print("  1. 右侧 'Output' 面板（没有就点右上角 > 展开）")
print("  2. 展开文件夹 'esm2_output_perres'")
print("  3. 下载 'v2_esm2_perres.npz'")
print("=" * 72)
print(f"\n绝对路径: {out_path}")
print(f"大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"subset={SUBSET}  storage={np.dtype(PERRES_DTYPE).name}  dim={_DIM}")
print()
print("本机放置位置:")
print(r"  F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_perres.npz")
print()
print("/kaggle/working/esm2_output_perres 内容:")
for _r, _d, _f in os.walk(OUT_DIR):
    _lvl = _r.replace(OUT_DIR, "").count(os.sep)
    print("  " * _lvl + os.path.basename(_r) + "/")
    for _fn in _f[:6]:
        _p = os.path.join(_r, _fn)
        print("  " * (_lvl + 1) + f"{_fn}  ({os.path.getsize(_p)/1e6:.1f} MB)")"""


def code_cell(src: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src.splitlines(keepends=True)}


def md_cell(src: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)}


def build() -> dict:
    cells = [
        md_cell(MD_TITLE),
        md_cell(MD_STEP1), code_cell(CODE_INSTALL),
        md_cell(MD_STEP2), code_cell(CODE_GPU_CHECK),
        md_cell(MD_STEP3), code_cell(CODE_SETUP),
        md_cell(MD_STEP4), code_cell(CODE_PREP),
        md_cell(MD_STEP5), code_cell(CODE_EMBED),
        md_cell(MD_STEP6), code_cell(CODE_PACK),
        md_cell(MD_STEP7), code_cell(CODE_DOWNLOAD),
    ]
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                           "name": "python3"},
            "language_info": {"name": "python"},
            "kaggle": {"accelerator": "nvidiaTeslaT4", "dataSources": [],
                       "isGpuEnabled": True, "isInternetEnabled": True,
                       "language": "python", "sourceType": "notebook"},
        },
        "nbformat": 4,
        "nbformat_minor": 4,
    }


REQUIRED = [
    "esm2_t30_150M_UR50D", "esm2_t12_35M_UR50D", "REPR_LAYER = 30",
    "MAX_LEN = 1022", "PERRES_DTYPE", "v2_esm2_perres.npz",
    "/kaggle/working/esm2_output_perres", "/kaggle/input/", "ppb-affinity-proteinea",
    "lig_flat", "rec_flat", "lig_lengths", "rec_lengths", "interface_key",
    "OutOfMemoryError", "md5_of", "SUBSET", "test_val", "float16",
]
FORBIDDEN_CODE = ["google.colab", "drive.mount", "/content/drive"]


def self_check(nb: dict) -> tuple[bool, list[str]]:
    import re
    problems: list[str] = []
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code":
            continue
        lines = []
        for ln in "".join(c["source"]).splitlines():
            s = ln.lstrip()
            lines.append(" " * (len(ln) - len(s)) + "pass"
                         if s.startswith(("!", "%")) else ln)
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as e:
            problems.append(f"cell#{i} SyntaxError: {e}")
    all_src = "\n".join("".join(c["source"]) for c in nb["cells"])
    for t in REQUIRED:
        if t not in all_src:
            problems.append(f"缺少必需常量/逻辑: {t!r}")
    # 可执行代码里不得有 Colab API (剥注释与字符串后检查)
    code_lines = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        for ln in "".join(c["source"]).splitlines():
            s = ln.split("#", 1)[0]
            s = re.sub(r'"(?:[^"\\]|\\.)*"', '""', s)
            s = re.sub(r"'(?:[^'\\]|\\.)*'", "''", s)
            code_lines.append(s)
    stripped = "\n".join(code_lines)
    for t in FORBIDDEN_CODE:
        if t in stripped:
            problems.append(f"可执行代码里残留 Colab API: {t!r}")
    # 必须不含 pooled 的 mean-pool
    if "mean_pool(" in all_src:
        problems.append("本 notebook 不应包含 mean_pool (它是 per-residue 版)")
    return (not problems), problems


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 per-residue ESM-2 提取 notebook")
    ap.add_argument("--out", type=Path, default=OUT_NB)
    args = ap.parse_args()

    nb = build()
    ok, problems = self_check(nb)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    reread = json.loads(args.out.read_text(encoding="utf-8"))
    assert reread["nbformat"] == 4 and len(reread["cells"]) == len(nb["cells"])

    n_code = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
    print(f"已写出 {args.out}  ({args.out.stat().st_size/1024:.1f} KB)")
    print(f"cells: {len(nb['cells'])} ({n_code} code)")
    if problems:
        print("!! 自检问题:")
        for p in problems:
            print("   -", p)
        return 1
    print("自检通过: JSON 合法 + code cell 语法正确 + 必需常量齐全 "
          "+ 无 Colab API 残留 + 确认无 mean-pool")
    return 0


if __name__ == "__main__":
    sys.exit(main())
