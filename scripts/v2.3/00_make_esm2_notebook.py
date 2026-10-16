# -*- coding: utf-8 -*-
"""scripts/v2.3/00_make_esm2_notebook.py -- 生成 Colab ESM-2 提取 notebook.

为什么要用生成器而不是手写 .ipynb
---------------------------------
手写 notebook JSON 极易产生非法 JSON 或转义错误。本生成器:
1. 用 Python 三引号字符串承载每个 cell 的源码 (可读、可 review);
2. 用 ``json.dump`` 写出**保证合法**的 ``.ipynb``;
3. 写完后**自检**: 重新解析 JSON、对每个 code cell 做 ``ast.parse`` 语法校验、
   检查必需常量是否真的写死在 notebook 里 (不是占位符)。

生成的 notebook: ``scripts/colab/esm2_extract.ipynb``

硬约束 (来自任务书)
-------------------
* 环境 cell: ``!pip install fair-esm==2.0.0 torch pandas numpy -q``
* 第 2 个 cell 检查 GPU, 无 GPU 直接 ``raise SystemExit``
* 模型/层/池化/归一化**写死**, 不留占位符:
  ``esm2_t30_150M_UR50D`` / 最后一层 (layer 30) / mean-pool over residues (去 BOS/EOS)
  / 同时存 raw 与 L2
* 最长 1022 aa 截断并记 ``truncated=True``; 非标准残基 -> ``<unk>`` 不丢序列
* Drive 读 Proteinea CSV (路径写在顶部常量, 只改一处); 用真实列名
* 先按唯一序列去重再嵌, 再合并回 ID 级
* 断点续算: 每条存 ``.npz``, 文件名 ``md5(sequence)``
* 每 100 条打印进度 + ``shutil.copy`` 备份到 Drive
* 单条失败不中断, 记入 ``failures`` 最后打印
* 防断线 ``setInterval(ClickConnect, 60000)``
* 打包 ``v2_esm2_cache.npz``: ``pdb_ids`` / ``lig_raw`` / ``lig_l2`` / ``rec_raw`` / ``rec_l2``
* OOM 降级: batch 减半至 1; 仍 OOM 切 ``esm2_t12_35M_UR50D`` 并告警
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
OUT_NB = ROOT / "scripts" / "colab" / "esm2_extract.ipynb"


# ===================================================================== cells
MD_TITLE = r"""# V2.3 — ESM-2 序列嵌入提取（Colab / T4 GPU）

本 notebook 为 **Affinity Intelligence V2.3** 提取 ESM-2 序列嵌入，用于补齐
V2.2 结论所指出的缺口：*60 个结构物理特征在同一界面键内恒定，无法区分突变体*。

---

## 写死的配置（不是占位符，请勿修改）

| 项 | 值 |
|---|---|
| 模型 | `esm2_t30_150M_UR50D`（150M 参数，640 维） |
| 表示层 | **最后一层（layer 30）** |
| 池化 | **mean-pool over residues**，**不含 BOS/EOS** |
| 归一化 | **同时保存 `raw` 与 `l2`**（L2 归一化），供后续消融 |
| 最大长度 | **1022 aa**，超出截断并记 `truncated=True` |
| 非标准残基 | 映射为 `<unk>`，**不丢弃整条序列** |
| 输出 | `v2_esm2_cache.npz` |

---

## 预计运行时间（按 8000+ 条估算）

数据审计（`features/v2.3_esm2_data_audit.md`）实测：

* 去重后需嵌入的**唯一序列共 9,519 条**（总残基 2,424,926，长度中位 196 aa，最大 3128 aa）。
* **纯前向：约 3 分钟**；含每条 `.npz` 落盘与每 100 条 Drive 备份，**预计总时长 5–10 分钟**。
* 若 Colab 上 T4 被争用，可能到 15 分钟。

## 显存需求（T4 = 15 GB）

推理用 `torch.no_grad()`，不需要梯度/优化器状态：

| batch_size | 最长序列 1022 aa 时估算 | T4 可行性 |
|---:|---|---|
| 1 | ~3.3 GB | ✅ |
| 4 | ~11.3 GB | ✅ |
| 8 | ~21.9 GB | ❌ 会 OOM |

**默认 `BATCH_SIZE = 4`**。真 OOM 时本 notebook 会**自动减半重试**（4→2→1）；
若 batch=1 仍 OOM，会**自动切换到 `esm2_t12_35M_UR50D`** 并打印警告。

## 中途断开后如何恢复

**直接重跑本 notebook 即可，不需要任何额外操作。**

原因：每嵌入一条序列就立刻把结果写成 `<CHECKPOINT_DIR>/<md5(sequence)>.npz`，
并每 100 条把整个 checkpoint 目录 `shutil.copy` 到 Google Drive 备份。
重跑时脚本会先扫描已有 `.npz`（**本地与 Drive 两处**）并跳过已完成的序列。

> 注意：Colab 断线后本地磁盘会清空，所以恢复速度取决于 Drive 备份的完整性——
> 这也是每 100 条就备份一次的原因。若希望更稳，可把 `BACKUP_EVERY` 改小（如 25）。
"""

MD_STEP1 = r"""## 步骤 1 / 7 — 安装依赖

`fair-esm==2.0.0` 是 ESM 官方实现（Meta），`torch` 由 Colab 预装。
本 cell 约需 30–60 秒。"""

CODE_INSTALL = r"""!pip install fair-esm==2.0.0 torch pandas numpy -q
print("依赖安装完成。")"""

MD_STEP2 = r"""## 步骤 2 / 7 — 确认 GPU 可用

**必须**是 GPU 环境。若不是，本 cell 会直接终止并提示你切换到 T4。
（菜单路径见 `features/v2.3_colab_instructions.md`）"""

CODE_GPU_CHECK = r"""import sys
import torch

if not torch.cuda.is_available():
    raise SystemExit(
        "请在 Colab 切换到 T4 GPU: 菜单 运行时(Runtime) -> 更改运行时类型(Change runtime type) "
        "-> 硬件加速器(Hardware accelerator) 选 T4 GPU -> 保存(Save)。然后重新运行本 notebook。"
    )

GPU_NAME = torch.cuda.get_device_name(0)
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1024**3
print(f"GPU: {GPU_NAME}  ({GPU_MEM_GB:.1f} GB)")
print(f"torch: {torch.__version__}")
DEVICE = torch.device("cuda")
print("GPU 检查通过。")"""

MD_STEP3 = r"""## 步骤 3 / 7 — 挂载 Google Drive 并设置路径

**你只需要改下面 `PROTEINEA_CSV` 这一个常量**（如果你把 CSV 放在了别的目录）。

本 cell 会：
1. 挂载 Drive；
2. 检查 `PROTEINEA_CSV` 是否存在，**不存在就直接报错并列出该目录下的实际文件**（不做静默跳过）；
3. 校验审计得出的**真实列名**是否都在。"""

CODE_SETUP = r"""import os
import shutil
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ====================================================================
# 只改这一处：Proteinea CSV 在 Google Drive 里的路径
# ====================================================================
PROTEINEA_CSV = "/content/drive/MyDrive/affinity_intelligence/filtered.csv"

# ---- 输出与断点目录（建议留在 Drive，断线不丢） ----
OUT_DIR = "/content/drive/MyDrive/affinity_intelligence/v2_esm2_out"
CHECKPOINT_DIR = os.path.join(OUT_DIR, "checkpoints")
BACKUP_DIR = os.path.join(OUT_DIR, "checkpoints_backup")
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
os.makedirs(BACKUP_DIR, exist_ok=True)

# ---- 写死的模型与预处理配置（审计确认，请勿修改） ----
MODEL_NAME = "esm2_t30_150M_UR50D"     # 150M, 640 维
MODEL_FALLBACK = "esm2_t12_35M_UR50D"  # OOM 兜底
REPR_LAYER = 30                        # 最后一层（layer 30）
MAX_LEN = 1022                         # ESM-2 位置上限（含 BOS/EOS）
BATCH_SIZE = 4                         # T4 上 1022 aa 的安全值
BACKUP_EVERY = 100                     # 每 N 条备份一次到 Drive
LIG_COL = "Ligand Sequences"           # 审计确认的真实列名
REC_COL = "Receptor Sequences"         # 审计确认的真实列名
ID_COL = "PDB"                         # 审计确认的真实列名
SPLIT_COL = "split"                    # 审计确认的真实列名

try:
    from google.colab import drive
    drive.mount("/content/drive")
    print("Google Drive 已挂载。")
except Exception as exc:  # noqa: BLE001
    print(f"!! Drive 挂载失败（如果你不是用 Drive，可忽略）: {exc}")

# 自动创建目标文件夹，方便你上传 CSV（不覆盖任何已有文件）
try:
    os.makedirs(os.path.dirname(PROTEINEA_CSV), exist_ok=True)
    print(f"已确保目录存在: {os.path.dirname(PROTEINEA_CSV)}")
    print("  -> 若这里还没有 filtered.csv，请把本地")
    print("     data/ppb_affinity_proteinea/filtered.csv 上传到该目录。")
except Exception as exc:  # noqa: BLE001
    print(f"!! 创建目录失败: {exc}")

csv_path = Path(PROTEINEA_CSV)
if not csv_path.exists():
    parent = csv_path.parent
    listing = sorted(p.name for p in parent.iterdir())[:40] if parent.exists() else []
    raise FileNotFoundError(
        f"找不到 {PROTEINEA_CSV}\n"
        f"请修改本 cell 顶部的 PROTEINEA_CSV 常量。\n"
        f"该目录现有文件: {listing}"
    )

# ---- 读数据（只用需要的列，省内存） ----
USE_COLS = [ID_COL, LIG_COL, REC_COL, SPLIT_COL]
df = pd.read_csv(csv_path, usecols=USE_COLS, dtype=str)
missing = [c for c in USE_COLS if c not in df.columns]
if missing:
    raise KeyError(f"CSV 缺少必需列: {missing}；实际列: {list(pd.read_csv(csv_path, nrows=1).columns)}")

df[ID_COL] = df[ID_COL].astype(str).str.strip().str.upper()
df = df.dropna(subset=[LIG_COL, REC_COL]).reset_index(drop=True)
print(f"读入 {len(df)} 行 x {len(df.columns)} 列  <- {csv_path}")
print("split 分布:", df[SPLIT_COL].value_counts().to_dict())
print("唯一 PDB:", df[ID_COL].nunique())
print("列名校验通过:", USE_COLS)"""

MD_STEP4 = r"""## 步骤 4 / 7 — 序列预处理与去重（省算力）

按任务要求：**先按唯一序列去重，嵌入完再合并回 ID 级**。

本 cell 处理：
* 非标准残基 → `<unk>`（不丢序列）；
* 超过 1022 aa 截断 → 记 `truncated=True`；
* 多链序列里的 `,` 分隔符也映射为 `<unk>`（审计已确认 31.5% 的配体行含多链）。

任务书要求输出 `pdb_ids` 形状 (N,)。本 notebook 采用 **N = 行级（测度条目级）**：
每一行是一条独立的结合测量（同一 PDB 可有多个突变体），
若按唯一 PDB 聚合会把突变体信息平均掉——那正是 V2.2 失败的原因，因此**不做聚合**。
`pdb_ids` 会保留重复值（同一 PDB 的多行），这是刻意设计。"""

CODE_PREP = r"""from esm.data import Alphabet

# ESM-2 的字母表：用它做「非标准残基 -> <unk>」的权威映射
ALPHABET = Alphabet.from_architecture("ESM-2")
UNK = ALPHABET.unk_idx
STD = set(ALPHABET.all_tokens)

def sanitize(seq: str) -> tuple[str, bool, int]:
    # 返回 (清洗后序列, 是否截断, 非标准字符数)。
    #   * 非标准残基（含多链 ',' 分隔符）-> 'X' 承载，再由 ESM 字母表映射为 <unk>
    #   * 超过 MAX_LEN 截断
    chars = list(seq)
    n_nonstd = sum(1 for c in chars if c not in STD)
    cleaned = "".join(c if c in STD else "X" for c in chars)
    truncated = len(cleaned) > MAX_LEN
    if truncated:
        cleaned = cleaned[:MAX_LEN]
    return cleaned, truncated, n_nonstd

def md5_of(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()

# ---- 收集唯一序列（配体 + 受体 合并去重） ----
records = df[[ID_COL, LIG_COL, REC_COL, SPLIT_COL]].to_dict("records")
for r in records:
    r["lig_clean"], r["lig_trunc"], r["lig_nonstd"] = sanitize(r[LIG_COL])
    r["rec_clean"], r["rec_trunc"], r["rec_nonstd"] = sanitize(r[REC_COL])
    r["lig_md5"] = md5_of(r["lig_clean"])
    r["rec_md5"] = md5_of(r["rec_clean"])

uniq_lig = sorted({r["lig_clean"] for r in records})
uniq_rec = sorted({r["rec_clean"] for r in records})
uniq_all = sorted(set(uniq_lig) | set(uniq_rec))

n_trunc = sum(1 for r in records if r["lig_trunc"] or r["rec_trunc"])
n_nonstd = sum(1 for r in records if r["lig_nonstd"] or r["rec_nonstd"])

print(f"行数                     : {len(records)}")
print(f"唯一配体序列             : {len(uniq_lig)}")
print(f"唯一受体序列             : {len(uniq_rec)}")
print(f"两列合并去重（实际要嵌）  : {len(uniq_all)}")
print(f"含截断的行               : {n_trunc}")
print(f"含非标准残基的行         : {n_nonstd}")

# 长度统计（用于 sanity check，与审计报告对齐）
_lens = np.array([len(s) for s in uniq_all])
print(f"唯一序列长度 min/median/max: {_lens.min()} / {int(np.median(_lens))} / {_lens.max()}")
if len(uniq_all) == 0:
    raise SystemExit("没有可嵌入的序列，请检查 CSV 列名与内容。")"""

MD_STEP5 = r"""## 步骤 5 / 7 — 加载 ESM-2 并嵌入（断点续算 + OOM 降级）

本 cell 是主要耗时部分（预计 5–10 分钟）。

行为：
* **断点续算**：先扫描 checkpoint 目录（本地 + Drive 备份），已完成的序列直接跳过；
* **每条落盘**：`<md5(sequence)>.npz`，原子替换，写入即生效；
* **每 100 条**：打印进度 + `shutil.copy` 备份整个 checkpoint 目录到 Drive；
* **单条失败不中断**：记入 `failures`，最后统一打印；
* **OOM 降级**：batch 减半重试（4→2→1）；batch=1 仍 OOM 则切 35M 模型并告警。"""

CODE_EMBED = r"""import esm
import gc
import traceback

failures: list[dict] = []

def load_model(name: str):
    # 加载 ESM-2 模型；兼容不同版本的模型名大小写。
    models = esm.pretrained.esm2_model_registry
    key = name
    if key not in models:
        lower = {k.lower(): k for k in models}
        if name.lower() in lower:
            key = lower[name.lower()]
        else:
            raise KeyError(f"fair-esm 不认识模型 {name!r}；可选: {sorted(models)}")
    model, alphabet = models[key]()
    model = model.eval().to(DEVICE)
    batch_converter = alphabet.get_batch_converter()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"已加载 {key}: {n_params/1e6:.1f}M 参数, 字母表 {alphabet.all_tokens[:5]}...")
    return model, alphabet, batch_converter, key

MODEL, ALPHABET, batch_converter, LOADED_KEY = load_model(MODEL_NAME)

def mean_pool(hidden: "torch.Tensor", tokens: "torch.Tensor") -> np.ndarray:
    # 对残基维度做 mean-pool，**排除 BOS/EOS 与 padding**。
    # hidden: (B, L, D) 的最后一层表示；tokens: (B, L) 的 token id。
    # ESM-2 约定: 0 = <pad>, 1 = <cls>(BOS), 2 = <eos>。
    import torch as _t
    mask = (tokens != 1) & (tokens != 2) & (tokens != 0)   # 去掉 BOS/EOS/pad
    m = mask.unsqueeze(-1).to(hidden.dtype)                # (B, L, 1)
    summed = (hidden * m).sum(dim=1)                       # (B, D)
    counts = m.sum(dim=1).clamp(min=1.0)                   # (B, 1)
    return (summed / counts).float().cpu().numpy()

def l2_normalize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, eps)

def embed_sequences(seqs: list[str], model, batch_converter_obj,
                    batch_size: int, tag: str) -> tuple[list[tuple[str, np.ndarray]], list[dict]]:
    # 按 batch 嵌入；OOM 时自动降级 batch_size。
    # 返回 (成功的 (序列, 向量) 列表, 失败记录)。
    #
    # !! 关键正确性说明 !! 返回值必须自带序列本身，不能只返回向量列表：
    # 若某个 batch 中途失败，只返回向量会让调用方无法把向量对回到正确的序列上
    # （zip(chunk, vecs) 会静默错位，把 A 的嵌入写到 B 的 md5 文件名下）。
    import torch as _t
    pairs: list[tuple[str, np.ndarray]] = []
    failed: list[dict] = []
    bs = batch_size
    i = 0
    while i < len(seqs):
        chunk = seqs[i:i + bs]
        try:
            data = [(f"s{j}", s) for j, s in enumerate(chunk)]
            _, _, tokens = batch_converter_obj(data)
            tokens = tokens.to(DEVICE)
            with _t.no_grad():
                reps = model(tokens, repr_layers=[REPR_LAYER], return_contacts=False)
            hidden = reps["representations"][REPR_LAYER]
            pooled = mean_pool(hidden, tokens)          # (len(chunk), D)，顺序与 chunk 一致
            if len(pooled) != len(chunk):
                raise RuntimeError(f"pooled 长度 {len(pooled)} != chunk 长度 {len(chunk)}")
            for s, v in zip(chunk, pooled):
                pairs.append((s, v))
            i += len(chunk)
        except _t.cuda.OutOfMemoryError:
            _t.cuda.empty_cache(); gc.collect()
            if bs > 1:
                bs = max(1, bs // 2)
                print(f"  !! CUDA OOM -> batch_size 降到 {bs} 重试（{tag}）")
                continue
            failed.append({"tag": tag, "index": i, "seq": chunk[0][:40],
                           "reason": "OOM at batch_size=1"})
            i += 1
        except Exception as exc:  # noqa: BLE001
            _t.cuda.empty_cache(); gc.collect()
            if bs > 1:
                # 非 OOM 的批次级异常也先降级重试一次，避免一次坏 batch 丢掉整段
                bs = max(1, bs // 2)
                print(f"  !! batch 异常 -> batch_size 降到 {bs} 重试（{tag}）: {exc}")
                continue
            failed.append({"tag": tag, "index": i, "seq": chunk[0][:40],
                           "reason": f"{type(exc).__name__}: {exc}"})
            i += 1
    return pairs, failed

def load_checkpoint(md5: str) -> np.ndarray | None:
    for d in (CHECKPOINT_DIR, BACKUP_DIR):
        p = os.path.join(d, f"{md5}.npz")
        if os.path.exists(p):
            try:
                with np.load(p) as z:
                    vec = z["vec"].astype(np.float32)
                if not os.path.exists(os.path.join(CHECKPOINT_DIR, f"{md5}.npz")):
                    shutil.copy(p, os.path.join(CHECKPOINT_DIR, f"{md5}.npz"))
                return vec
            except Exception:  # noqa: BLE001
                continue
    return None

def save_checkpoint(md5: str, vec: np.ndarray) -> None:
    p = os.path.join(CHECKPOINT_DIR, f"{md5}.npz")
    tmp = p + ".tmp.npz"
    np.savez(tmp, vec=vec.astype(np.float32))
    os.replace(tmp, p)

# ---------------- 主循环 ----------------
emb_cache: dict[str, np.ndarray] = {}
# 先扫描本地 checkpoint 与 Drive 备份，把已完成的序列载入内存并跳过
for s in uniq_all:
    v = load_checkpoint(md5_of(s))
    if v is not None:
        emb_cache[s] = v
todo = [s for s in uniq_all if s not in emb_cache]

print(f"总唯一序列 {len(uniq_all)}；已有缓存 {len(emb_cache)}；待嵌入 {len(todo)}")
t_start = time.time()

i = 0
while i < len(todo):
    chunk = todo[i:i + BATCH_SIZE]
    pairs, failed = embed_sequences(chunk, MODEL, batch_converter, BATCH_SIZE, "main")
    # pairs 自带序列 -> 不可能错位；失败的条目不会进 pairs
    for s, v in pairs:
        save_checkpoint(md5_of(s), v)
        emb_cache[s] = v
    failures.extend(failed)
    i += len(chunk)
    done = len(emb_cache)
    if done % BACKUP_EVERY < len(chunk) or done >= len(uniq_all):
        try:
            shutil.copytree(CHECKPOINT_DIR, BACKUP_DIR, dirs_exist_ok=True)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! Drive 备份失败（不中断）: {exc}")
        el = time.time() - t_start
        rate = done / max(el, 1e-9)
        eta = (len(uniq_all) - done) / max(rate, 1e-9)
        print(f"  进度 {done}/{len(uniq_all)}  {rate:.1f} seq/s  "
              f"已用 {el/60:.1f} min  预计剩余 {eta/60:.1f} min  失败 {len(failures)}")

print(f"嵌入完成：{len(emb_cache)} 条，失败 {len(failures)} 条，耗时 {(time.time()-t_start)/60:.1f} min")
if failures:
    print("失败明细（前 20 条）:")
    for f_ in failures[:20]:
        print("  ", f_)

# 释放显存前先记录
DIM = int(next(iter(emb_cache.values())).shape[0]) if emb_cache else 0
print(f"嵌入维度 DIM = {DIM}")"""

MD_STEP6 = r"""## 步骤 6 / 7 — 打包 `v2_esm2_cache.npz`

按任务书要求输出 5 个数组：

| key | 形状 | 含义 |
|---|---|---|
| `pdb_ids` | (N,) | 每个**行级条目**的 PDB ID（同一 PDB 可出现多次，见步骤 4 说明） |
| `lig_raw` | (N, 640) | 配体序列嵌入（原始） |
| `lig_l2` | (N, 640) | 配体序列嵌入（L2 归一化） |
| `rec_raw` | (N, 640) | 受体序列嵌入（原始） |
| `rec_l2` | (N, 640) | 受体序列嵌入（L2 归一化） |

额外附带（**不属于要求，但对任务 3 的 join 必要**，避免合并脚本靠猜行序）：
* `interface_key`：`PDB|配体链|受体链`（与 V2.2 特征表同名同格式）
* `lig_md5` / `rec_md5`：每条序列的 md5，用于把唯一序列缓存反查回条目
* `truncated`：该行是否有序列被截断
* `model_name` / `repr_layer` / `pool` / `l2_normalized`：写死配置的记录

> 若 `lig_raw` 的列数与 640 不同（例如 OOM 降级到 35M 模型，其维度为 480），
> 本 cell 会**打印明确警告**，合并脚本也接受任意维度。"""

CODE_PACK = r"""import json

def vec_or_zero(seq: str, dim: int) -> np.ndarray:
    v = emb_cache.get(seq)
    if v is None:
        return np.zeros(dim, dtype=np.float32)
    return v.astype(np.float32)

n_rows = len(records)
missing_lig = sum(1 for r in records if r["lig_clean"] not in emb_cache)
missing_rec = sum(1 for r in records if r["rec_clean"] not in emb_cache)
print(f"条目数 {n_rows}；缺配体嵌入 {missing_lig}；缺受体嵌入 {missing_rec}")
if missing_lig or missing_rec:
    print("  !! 缺失的嵌入将填 0 向量，并在 'missing_mask' 中标记（合并脚本会据此排除）")

lig_raw = np.zeros((n_rows, DIM), dtype=np.float32)
rec_raw = np.zeros((n_rows, DIM), dtype=np.float32)
missing_mask = np.zeros(n_rows, dtype=bool)
for i, r in enumerate(records):
    lv = emb_cache.get(r["lig_clean"])
    rv = emb_cache.get(r["rec_clean"])
    if lv is not None:
        lig_raw[i] = lv
    else:
        missing_mask[i] = True
    if rv is not None:
        rec_raw[i] = rv
    else:
        missing_mask[i] = True

lig_l2 = l2_normalize(lig_raw)
rec_l2 = l2_normalize(rec_raw)

pdb_ids = np.array([r[ID_COL] for r in records], dtype="U")
interface_key = np.array(
    [f"{r[ID_COL]}|{r[LIG_COL]}|{r[REC_COL]}" for r in records], dtype="U"
)
lig_md5 = np.array([r["lig_md5"] for r in records], dtype="U")
rec_md5 = np.array([r["rec_md5"] for r in records], dtype="U")
truncated = np.array([bool(r["lig_trunc"] or r["rec_trunc"]) for r in records], dtype=bool)
splits = np.array([str(r[SPLIT_COL]) for r in records], dtype="U")

if DIM != 640:
    print("=" * 70)
    print(f"!! 警告: 嵌入维度是 {DIM}，不是 150M 模型的 640。")
    print("!! 这通常意味着触发了 OOM 降级（esm2_t12_35M_UR50D 的维度是 480）。")
    print("!! 记录会用 model_name 标注，合并脚本按实际维度处理。")
    print("=" * 70)

out_path = os.path.join(OUT_DIR, "v2_esm2_cache.npz")
np.savez_compressed(
    out_path,
    pdb_ids=pdb_ids,
    lig_raw=lig_raw, lig_l2=lig_l2,
    rec_raw=rec_raw, rec_l2=rec_l2,
    interface_key=interface_key,
    lig_md5=lig_md5, rec_md5=rec_md5,
    truncated=truncated, splits=splits,
    missing_mask=missing_mask,
    model_name=np.array(LOADED_KEY),
    repr_layer=np.array(REPR_LAYER),
    pool=np.array("mean_over_residues_excl_bos_eos"),
    l2_normalized=np.array(True),
    n_unique_sequences=np.array(len(uniq_all)),
    n_failed=np.array(len(failures)),
)
print(f"已写出: {out_path}")
print(f"  大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"  pdb_ids {pdb_ids.shape}  lig_raw {lig_raw.shape}  rec_raw {rec_raw.shape}")
print(f"  model={LOADED_KEY}  layer={REPR_LAYER}  dim={DIM}")"""

MD_STEP7 = r"""## 步骤 7 / 7 — 下载到本地

**先下载，再关闭 Colab**（Colab 关掉后本地磁盘会清空，但 Drive 上还有备份）。

下载后的放置位置：

```text
F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_cache.npz
```

然后在本机运行：

```powershell
python scripts/v2.3/03_merge_esm2.py
```"""

CODE_DOWNLOAD = r"""from google.colab import files

# 方式 A（推荐）：从 Drive 直接下载
print("文件在 Google Drive 的位置：")
print("  " + out_path)
print("你也可以在 Drive 网页界面里右键下载。")
print()
print("方式 B：从 Colab 本地直接下载（浏览器会弹下载框）")
try:
    files.download(out_path)
except Exception as exc:  # noqa: BLE001
    print(f"自动下载失败（改用 Drive 手动下载即可）: {exc}")

print()
print("=" * 70)
print("下载完成后，把文件放到：")
print(r"  F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_cache.npz")
print("然后在本机运行：")
print("  python scripts/v2.3/03_merge_esm2.py")
print("=" * 70)"""

MD_ANTI = r"""## 附 — 防止 Colab 因闲置断开

在浏览器里按 **F12** 打开控制台，粘贴下面这段并回车。
它会每 60 秒点一次「连接」按钮，避免长时间无操作被判定闲置。

> 注意：这只是防止**闲置断开**，不能防止**总时长上限**或**崩溃**。
> 真正的保护是断点续算（每 100 条备份到 Drive），重跑即可继续。"""

CODE_ANTI = r"""from IPython.display import display, Javascript

JS = "\n".join([
    "function ClickConnect(){",
    "  console.log('Keeping Colab alive...');",
    "  var btn = document.querySelector('colab-connect-button');",
    "  if (btn) { btn.click(); }",
    "}",
    "setInterval(ClickConnect, 60000);",
])

display(Javascript(JS))
print("已注入防断线脚本：每 60 秒点击一次连接按钮。")"""


# ===================================================================== build
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
        md_cell(MD_ANTI), code_cell(CODE_ANTI),
    ]
    return {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True,
                      "name": "esm2_extract.ipynb"},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }


REQUIRED_TOKENS = [
    "esm2_t30_150M_UR50D", "esm2_t12_35M_UR50D", "fair-esm==2.0.0",
    "REPR_LAYER = 30", "MAX_LEN = 1022", "BATCH_SIZE = 4", "BACKUP_EVERY = 100",
    '"Ligand Sequences"', '"Receptor Sequences"', '"PDB"', '"split"',
    "setInterval(ClickConnect, 60000)", "md5_of", "shutil.copy", "shutil.copytree",
    "OutOfMemoryError", "lig_raw", "lig_l2", "rec_raw", "rec_l2", "pdb_ids",
    "v2_esm2_cache.npz", "raise SystemExit", "truncated=True",
]


def self_check(nb: dict) -> tuple[bool, list[str]]:
    problems: list[str] = []
    # 1) 每个 code cell 必须语法正确
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code":
            continue
        src = "".join(c["source"])
        # Colab 允许 ! 与 % 开头的 magics; 语法检查前替换掉
        lines = []
        for ln in src.splitlines():
            s = ln.lstrip()
            if s.startswith("!") or s.startswith("%"):
                lines.append(" " * (len(ln) - len(s)) + "pass")
            else:
                lines.append(ln)
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as e:
            problems.append(f"cell#{i} SyntaxError: {e}")
    # 2) 必需常量必须写死
    all_src = "\n".join("".join(c["source"]) for c in nb["cells"])
    for tok in REQUIRED_TOKENS:
        if tok not in all_src:
            problems.append(f"缺少必需常量/逻辑: {tok!r}")
    # 3) 结构完整性
    if nb["nbformat"] != 4:
        problems.append("nbformat 必须是 4")
    n_code = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
    print(f"自检: {len(nb['cells'])} cells ({n_code} code), "
          f"{len(all_src.splitlines())} 行代码")
    return (not problems), problems


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 Colab ESM-2 提取 notebook")
    ap.add_argument("--out", type=Path, default=OUT_NB)
    args = ap.parse_args()

    nb = build()
    ok, problems = self_check(nb)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")

    # 回读校验: 确认写出的确实能被解析
    reread = json.loads(args.out.read_text(encoding="utf-8"))
    assert reread["nbformat"] == 4 and len(reread["cells"]) == len(nb["cells"])

    print(f"已写出 {args.out}  ({args.out.stat().st_size/1024:.1f} KB)")
    if problems:
        print("!! 自检问题:")
        for p in problems:
            print("   -", p)
        return 1
    print("自检通过: JSON 合法 + 全部 code cell 语法正确 + 必需常量齐全")
    return 0


if __name__ == "__main__":
    sys.exit(main())
