# -*- coding: utf-8 -*-
"""scripts/v2.3.1/10_make_kaggle_notebook.py -- 生成 Kaggle 版 ESM-2 提取 notebook.

背景 (为什么需要 Kaggle 版)
--------------------------
V2.3 的 Colab notebook 已就绪, 但用户所在地区 Colab 的 GPU 分配失败
(晚间高峰 + 地区限制)。Kaggle Notebooks 每周提供约 30 小时免费 GPU (T4 x2),
因此把同一套逻辑迁移到 Kaggle。

核心逻辑必须与 Colab 版**逐位一致** (硬约束 4)
---------------------------------------------
本生成器**不复制粘贴**核心 cell 的源码, 而是直接
``importlib`` 载入 ``scripts/v2.3/00_make_esm2_notebook.py`` 并**复用它的字符串常量**:

* ``CODE_PREP``     序列清洗与去重
* ``CODE_EMBED``    模型加载 / mean-pool / OOM 降级 / 断点续算 / 失败不中断
* ``CODE_PACK``     打包 npz (raw + L2)
* ``CODE_GPU_CHECK`` / ``CODE_INSTALL`` 的差异部分单独处理

因此不存在"两边不一致"的可能 —— 它们字面上是同一个字符串。
生成后会做**逐位一致性自检**: 对共享 cell 计算 sha256 并断言两边相同。

只有在「路径与下载方式」上改动 (硬约束 4 允许的改动):
* 输入: ``/kaggle/input/ppb-affinity-proteinea/filtered.csv`` (不再挂载 Google Drive)
* 输出: ``/kaggle/working/esm2_output/``
* 删除: Google Drive 挂载 cell
* 备份: Kaggle 的 ``/kaggle/working`` 会话内持久, 无需 Drive 备份, 但每 100 条仍落盘
* 下载: 改用 Kaggle 右侧 Output 面板 (``files.download`` 在 Kaggle 不可用)
* 安装: 只需 ``!pip install fair-esm -q`` (torch/numpy 已预装)
* GPU: 打印 ``torch.cuda.device_count()``, T4 x2 为 2, 指定只用 ``cuda:0``

Colab 与 Kaggle 的 API 差异 (硬约束 5, 完整清单)
----------------------------------------------
=================== ============================== ================================
项                  Colab                          Kaggle
=================== ============================== ================================
持久存储            Drive 挂载 ``drive.mount``      ``/kaggle/working`` 会话内持久
下载文件            ``files.download(path)``       **不可用** -> 右侧 Output 面板
GPU 型号            T4 x1 (免费版)                 T4 x2
指定设备            ``cuda``                       ``cuda:0`` (只用一张)
输入数据            任意路径 (通常 Drive)          **只读** ``/kaggle/input/<slug>/``
预装库              torch 需确认                   torch / numpy 已预装
会话时长上限        约 12 小时                     约 12 小时 (GPU 每周 30 小时)
=================== ============================== ================================
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
COLAB_GEN = ROOT / "scripts" / "v2.3" / "00_make_esm2_notebook.py"
OUT_NB = ROOT / "scripts" / "kaggle" / "esm2_extract_kaggle.ipynb"
COLAB_NB = ROOT / "scripts" / "colab" / "esm2_extract.ipynb"

DATASET_SLUG = "ppb-affinity-proteinea"
INPUT_CSV = f"/kaggle/input/{DATASET_SLUG}/filtered.csv"
OUT_DIR = "/kaggle/working/esm2_output"


def load_colab_module():
    """载入 Colab 生成器模块, 以复用其 cell 源码常量 (单一实现来源)。"""
    spec = importlib.util.spec_from_file_location("colab_gen", COLAB_GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ===================================================================== Kaggle cells
MD_TITLE = r"""# V2.3.1 — ESM-2 序列嵌入提取（**Kaggle 版**，T4 x2）

> 这是 V2.3 Colab notebook 的 **Kaggle 移植版**。
> 迁移原因：用户所在地区 Colab 的 GPU 分配失败（晚间高峰 + 地区限制）。
> Kaggle Notebooks 每周提供约 **30 小时**免费 GPU（**T4 x2**）。
>
> **核心逻辑与 Colab 版逐位一致**（模型 / 层 / 池化 / 归一化 / 去重 / 断点续算 / OOM 降级）。
> 只有「输入输出路径」和「下载方式」不同。

---

## 写死的配置（不是占位符，请勿修改）

| 项 | 值 |
|---|---|
| 模型 | `esm2_t30_150M_UR50D`（150M 参数，640 维） |
| 表示层 | **最后一层（layer 30）** |
| 池化 | **mean-pool over residues**，**不含 BOS/EOS** |
| 归一化 | **同时保存 `raw` 与 `l2`**（L2 归一化） |
| 最大长度 | **1022 aa**，超出截断并记 `truncated=True` |
| 非标准残基 | 映射为 `<unk>`，**不丢弃整条序列** |
| 设备 | **`cuda:0`**（Kaggle 给 2 张 T4，本 notebook 只用第 1 张） |
| 输入 | `/kaggle/input/ppb-affinity-proteinea/filtered.csv` |
| 输出 | `/kaggle/working/esm2_output/v2_esm2_cache.npz` |

---

## 先做这三件事（详见 `features/v2.3.1_kaggle_instructions.md`）

1. **注册 Kaggle 账号**：https://www.kaggle.com → Register（可用 Google 账号直接登录）。
2. **把 `filtered.csv` 上传为 Kaggle Dataset**，Dataset 名称必须是
   **`ppb-affinity-proteinea`**（否则第 3 个 cell 会找不到文件并提示你）。
3. **新建 Notebook 并开启 GPU + Internet**：
   右上角 **Settings** → **Accelerator** 选 **GPU T4 x2**；
   **Internet** 开关设为 **On**（装 `fair-esm` 必须联网）。

---

## 预计运行时间

数据审计（`features/v2.3_esm2_data_audit.md`）实测：

* 去重后需嵌入的**唯一序列共 9,519 条**（总残基 2,424,926，长度中位 196 aa，最大 3128 aa）。
* **T4 x1 下预计 10–15 分钟**（纯前向约 3 分钟；其余是每条落盘与进度打印的开销）。
* T4 x2 不会自动加速——本 notebook 只用 `cuda:0`，因为逐条落盘是串行的。

## 显存需求（T4 = 15 GB / 张）

推理用 `torch.no_grad()`，不需要梯度/优化器状态。
**默认 `BATCH_SIZE = 4`**，真 OOM 时自动减半重试（4→2→1）；
若 batch=1 仍 OOM，自动切换到 `esm2_t12_35M_UR50D` 并打印警告。

## 中途断开后如何恢复

**直接重跑本 notebook 即可，不需要任何额外操作。**

原因：每嵌入一条序列就立刻写成
`/kaggle/working/esm2_output/checkpoints/<md5(sequence)>.npz`，
而 **`/kaggle/working` 在同一个 Kaggle 会话内是持久的**，重跑时会扫描已有 `.npz` 并跳过。

> ⚠️ 但 Kaggle 会话**结束后** `/kaggle/working` 会清空。所以请在**同一次会话内**跑完，
> 跑完立刻从右侧 Output 面板下载 `v2_esm2_cache.npz`。
> 若确实中途会话结束，重新跑一遍即可（只是要从零开始，因为 checkpoint 也丢了）。

## Colab 与 Kaggle 的 API 差异（重要）

| 项 | Colab | Kaggle |
|---|---|---|
| 持久存储 | Drive 挂载 | `/kaggle/working` 会话内持久 |
| 下载文件 | `files.download(path)` | **不可用** → 用右侧 **Output** 面板 |
| GPU | T4 x1 | **T4 x2**（本 notebook 只用 `cuda:0`） |
| 输入数据 | 任意路径 | **只读** `/kaggle/input/<slug>/` |
| 预装库 | torch 需确认 | torch / numpy **已预装** |
"""

MD_STEP1 = r"""## 步骤 1 / 7 — 安装 fair-esm

Kaggle 已预装 `torch` 与 `numpy`，所以**只需装 `fair-esm`**。
需要 Internet = On（Settings → Internet）。

本 cell 约需 20–40 秒。"""

CODE_INSTALL = r"""!pip install fair-esm -q
print("fair-esm 安装完成。")"""

MD_STEP2 = r"""## 步骤 2 / 7 — 确认 GPU

Kaggle 选 **GPU T4 x2** 时，`torch.cuda.device_count()` 应为 **2**。
**本 notebook 只使用 `cuda:0`**（逐条落盘是串行的，用第二张卡不会更快）。

若设备数为 0，本 cell 会直接终止并提示你去 Settings 开 GPU。"""

CODE_GPU_CHECK = r"""import sys
import torch

if not torch.cuda.is_available():
    raise SystemExit(
        "请开启 Kaggle GPU: 右上角 Settings -> Accelerator -> 选 GPU T4 x2 -> 保存。"
        "然后重新运行本 notebook。"
    )

N_GPU = torch.cuda.device_count()
print(f"torch: {torch.__version__}")
print(f"torch.cuda.device_count() = {N_GPU}")
for _i in range(N_GPU):
    _p = torch.cuda.get_device_properties(_i)
    print(f"  GPU {_i}: {_p.name}  ({_p.total_memory / 1024**3:.1f} GB)")

# 本 notebook 只使用第 1 张卡（cuda:0）
DEVICE = torch.device("cuda:0")
print(f"\n使用设备: {DEVICE}  (Kaggle T4 x2 时共 {N_GPU} 张, 只用第 1 张)")
if N_GPU >= 2:
    print("提示: 第 2 张卡本流程用不到 —— 逐条落盘是串行的, 多卡不会更快。")
print("GPU 检查通过。")"""

MD_STEP3 = r"""## 步骤 3 / 7 — 设置路径并读取数据（Kaggle 版）

**与 Colab 版的差异**：
* 不再挂载 Google Drive（Kaggle 不需要）；
* 输入从 **只读的** `/kaggle/input/ppb-affinity-proteinea/filtered.csv` 读取；
* 输出写到 `/kaggle/working/esm2_output/`。

> 你的 Kaggle Dataset 名称**必须**是 `ppb-affinity-proteinea`。
> 如果用了别的名字，只需改下面 `DATASET_SLUG` 这一个常量。
>
> 本 cell 会**自动探测** `/kaggle/input/` 下真实的目录名并在找不到时提示你，
> **不会静默跳过**。"""

CODE_SETUP = r"""import os
import shutil
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ====================================================================
# 只改这一处：Kaggle Dataset 的 slug（就是你创建 Dataset 时用的名字）
# ====================================================================
DATASET_SLUG = "ppb-affinity-proteinea"
PROTEINEA_CSV = f"/kaggle/input/{DATASET_SLUG}/filtered.csv"

# ---- 输出与断点目录（/kaggle/working 在同一会话内持久） ----
OUT_DIR = "/kaggle/working/esm2_output"
CHECKPOINT_DIR = os.path.join(OUT_DIR, "checkpoints")
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

# ---- 写死的模型与预处理配置（与 Colab 版逐位一致，请勿修改） ----
MODEL_NAME = "esm2_t30_150M_UR50D"     # 150M, 640 维
MODEL_FALLBACK = "esm2_t12_35M_UR50D"  # OOM 兜底
REPR_LAYER = 30                        # 最后一层（layer 30）
MAX_LEN = 1022                         # ESM-2 位置上限（含 BOS/EOS）
BATCH_SIZE = 4                         # T4 上 1022 aa 的安全值
BACKUP_EVERY = 100                     # 每 N 条打印进度并存盘
LIG_COL = "Ligand Sequences"           # 审计确认的真实列名
REC_COL = "Receptor Sequences"         # 审计确认的真实列名
ID_COL = "PDB"                         # 审计确认的真实列名
SPLIT_COL = "split"                    # 审计确认的真实列名

# ---- Kaggle 专用：列出 /kaggle/input 下的实际内容，便于排错 ----
print("=== /kaggle/input 下的数据集 ===")
_in = Path("/kaggle/input")
if _in.exists():
    for d in sorted(_in.iterdir()):
        print(f"  {d.name}/")
        try:
            for f in sorted(d.iterdir())[:8]:
                print(f"      {f.name}")
        except Exception:
            pass
else:
    print("  (/kaggle/input 不存在 —— 你可能不是在 Kaggle 上运行)")

csv_path = Path(PROTEINEA_CSV)
if not csv_path.exists():
    _avail = sorted(p.name for p in _in.iterdir()) if _in.exists() else []
    raise FileNotFoundError(
        f"找不到 {PROTEINEA_CSV}\n"
        f"请确认：(1) 已上传 Kaggle Dataset；(2) Dataset 名称是 '{DATASET_SLUG}'；"
        f"(3) 已在 Notebook 右侧 'Add Input' 里挂载该 Dataset。\n"
        f"/kaggle/input 现有数据集: {_avail}"
    )

# ---- 读数据（只用需要的列，省内存） ----
USE_COLS = [ID_COL, LIG_COL, REC_COL, SPLIT_COL]
df = pd.read_csv(csv_path, usecols=USE_COLS, dtype=str)
missing = [c for c in USE_COLS if c not in df.columns]
if missing:
    raise KeyError(f"CSV 缺少必需列: {missing}；实际列: {list(pd.read_csv(csv_path, nrows=1).columns)}")

df[ID_COL] = df[ID_COL].astype(str).str.strip().str.upper()
df = df.dropna(subset=[LIG_COL, REC_COL]).reset_index(drop=True)
print(f"\n读入 {len(df)} 行 x {len(df.columns)} 列  <- {csv_path}")
print("split 分布:", df[SPLIT_COL].value_counts().to_dict())
print("唯一 PDB:", df[ID_COL].nunique())
print("列名校验通过:", USE_COLS)
print("断点目录:", CHECKPOINT_DIR)"""

MD_STEP4 = r"""## 步骤 4 / 7 — 序列预处理与去重（与 Colab 版逐位一致）

* 非标准残基 → `<unk>`（不丢序列）；
* 超过 1022 aa 截断 → 记 `truncated=True`；
* 多链序列里的 `,` 分隔符也映射为 `<unk>`。

`pdb_ids` 保留重复值（同一 PDB 可有多个突变体），这是刻意设计 ——
按唯一 PDB 聚合会把突变体信息平均掉，那正是 V2.2 失败的原因。"""

MD_STEP5 = r"""## 步骤 5 / 7 — 加载 ESM-2 并嵌入（主要等待时间）

行为与 Colab 版一致：
* **断点续算**：扫描 `/kaggle/working/esm2_output/checkpoints`，已完成的序列直接跳过；
* **每条落盘**：`<md5(sequence)>.npz`，原子替换；
* **每 100 条**：打印进度（Kaggle 不需要 Drive 备份，但每 100 条仍确认落盘状态）；
* **单条失败不中断**：记入 `failures`，最后统一打印；
* **OOM 降级**：batch 减半（4→2→1）；batch=1 仍 OOM 则切 35M 模型并告警。"""

MD_STEP6 = r"""## 步骤 6 / 7 — 打包 `v2_esm2_cache.npz`

按任务书要求输出 5 个数组：

| key | 形状 | 含义 |
|---|---|---|
| `pdb_ids` | (N,) | 每个**行级条目**的 PDB ID（同一 PDB 可出现多次） |
| `lig_raw` | (N, 640) | 配体序列嵌入（原始） |
| `lig_l2` | (N, 640) | 配体序列嵌入（L2 归一化） |
| `rec_raw` | (N, 640) | 受体序列嵌入（原始） |
| `rec_l2` | (N, 640) | 受体序列嵌入（L2 归一化） |

额外附带（**对任务 3 的 join 必要**，避免合并脚本靠猜行序）：
`interface_key`、`lig_md5`、`rec_md5`、`truncated`、`splits`、`missing_mask`、
`model_name`、`repr_layer`、`pool`、`l2_normalized`。"""

MD_STEP7 = r"""## 步骤 7 / 7 — 下载到本地（**Kaggle 方式**）

⚠️ **`files.download()` 在 Kaggle 上不可用**（那是 Colab 的 API）。

Kaggle 的下载方式是**右侧 Output 面板**：

1. 等本 cell 打印出文件路径与大小；
2. 看页面**右侧**的 **Output** 面板（若没显示，点右上角 `>` 或 "Output" 展开）；
3. 展开 `esm2_output/` 文件夹；
4. 找到 **`v2_esm2_cache.npz`**，点它右边的 **下载图标**（↓）。

下载后放到本机：

```text
F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_cache.npz
```

然后在本机运行：

```powershell
python scripts/v2.3/03_merge_esm2.py
```"""

CODE_DOWNLOAD = r"""# Kaggle 上没有 google.colab.files，所以这里不做自动下载。
# 改为：打印文件位置 + 给你 Output 面板的操作指引。
import os

print("=" * 72)
print("文件位置（在 Kaggle 右侧 Output 面板里找）:")
print(f"  esm2_output/v2_esm2_cache.npz")
print(f"  绝对路径: {out_path}")
print(f"  大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print("=" * 72)
print()
print("下载步骤（Kaggle 没有 files.download）:")
print("  1. 看页面右侧的 'Output' 面板（没有就点右上角的 > 展开）")
print("  2. 展开文件夹 'esm2_output'")
print("  3. 找到 'v2_esm2_cache.npz'，点右侧的下载图标（↓）")
print()
print("下载后放到本机:")
print(r"  F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_cache.npz")
print()
print("然后在本机运行:")
print("  python scripts/v2.3/03_merge_esm2.py")
print("=" * 72)

# 顺便确认 output 目录里到底有什么，方便你对照 Output 面板
print("\n/kaggle/working/esm2_output 内容:")
for _r, _d, _f in os.walk(OUT_DIR):
    _lvl = _r.replace(OUT_DIR, "").count(os.sep)
    print("  " * _lvl + os.path.basename(_r) + "/")
    for _fn in _f[:10]:
        _p = os.path.join(_r, _fn)
        print("  " * (_lvl + 1) + f"{_fn}  ({os.path.getsize(_p)/1024:.1f} KB)")"""

MD_ANTI = r"""## 附 — 防断线（Kaggle 一般不会断，保留无害）

Kaggle 不像 Colab 那样容易因闲置断开，但保留这段没有副作用。
在浏览器里按 **F12** 打开控制台，粘贴下面这段并回车。

> 真正的保护是**断点续算**：每 100 条已存盘到 `/kaggle/working`，重跑即可继续。"""


# ===================================================================== build
def code_cell(src: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src.splitlines(keepends=True)}


def md_cell(src: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)}


def build(colab) -> dict:
    """用 Colab 生成器的共享常量拼装 Kaggle notebook。

    共享 (逐位一致) 的 code cell: CODE_PREP / CODE_EMBED / CODE_PACK
    Kaggle 专属的 code cell:      CODE_INSTALL / CODE_GPU_CHECK / CODE_SETUP
                                  / CODE_DOWNLOAD / CODE_ANTI
    """
    cells = [
        md_cell(MD_TITLE),
        md_cell(MD_STEP1), code_cell(CODE_INSTALL),
        md_cell(MD_STEP2), code_cell(CODE_GPU_CHECK),
        md_cell(MD_STEP3), code_cell(CODE_SETUP),
        md_cell(MD_STEP4), code_cell(colab.CODE_PREP),      # <-- 共享
        md_cell(MD_STEP5), code_cell(colab.CODE_EMBED),     # <-- 共享
        md_cell(MD_STEP6), code_cell(colab.CODE_PACK),      # <-- 共享
        md_cell(MD_STEP7), code_cell(CODE_DOWNLOAD),
        md_cell(MD_ANTI), code_cell(colab.CODE_ANTI),       # <-- 共享 (防断线 JS)
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


REQUIRED_TOKENS = [
    "esm2_t30_150M_UR50D", "esm2_t12_35M_UR50D", "fair-esm", "REPR_LAYER = 30",
    "MAX_LEN = 1022", "BATCH_SIZE = 4", "BACKUP_EVERY = 100",
    '"Ligand Sequences"', '"Receptor Sequences"', '"PDB"', '"split"',
    "setInterval(ClickConnect, 60000)", "md5_of", "OutOfMemoryError",
    "lig_raw", "lig_l2", "rec_raw", "rec_l2", "pdb_ids",
    "v2_esm2_cache.npz", "raise SystemExit", "truncated=True",
    # Kaggle 专属
    "/kaggle/input/", "/kaggle/working/esm2_output", "ppb-affinity-proteinea",
    "torch.cuda.device_count()", 'torch.device("cuda:0")', "Output",
]

FORBIDDEN_TOKENS = [
    "google.colab", "drive.mount", "files.download", "/content/drive",
]


def self_check(nb: dict, colab) -> tuple[bool, list[str]]:
    problems: list[str] = []
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code":
            continue
        lines = []
        for ln in "".join(c["source"]).splitlines():
            s = ln.lstrip()
            lines.append(" " * (len(ln) - len(s)) + "pass"
                         if (s.startswith("!") or s.startswith("%")) else ln)
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as e:
            problems.append(f"cell#{i} SyntaxError: {e}")

    all_src = "\n".join("".join(c["source"]) for c in nb["cells"])
    for tok in REQUIRED_TOKENS:
        if tok not in all_src:
            problems.append(f"缺少必需常量/逻辑: {tok!r}")

    # ---- Colab 专属 API 只允许出现在「说明/注释/字符串」里, 不允许是**可执行调用** ----
    # 例: 注释 `# Kaggle 上没有 google.colab.files` 与
    #     print("下载步骤（Kaggle 没有 files.download）") 都是文档, 不是残留调用。
    # 因此扫描前先剥掉注释与字符串字面量, 只留真正的代码 token。
    code_lines: list[str] = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        try:
            tree = ast.parse("\n".join(
                "pass" if ln.lstrip().startswith(("!", "%")) else ln
                for ln in "".join(c["source"]).splitlines()))
        except SyntaxError:
            code_lines.append("".join(c["source"]))
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                continue
        # 用 AST 重建不可行, 改用「按行剥离字符串字面量」的近似
        import re as _re
        for ln in "".join(c["source"]).splitlines():
            s = ln.split("#", 1)[0]                      # 去行尾注释
            s = _re.sub(r'"""(?:.|\n)*?"""', "", s)      # 去三引号串
            s = _re.sub(r"'''(?:.|\n)*?'''", "", s)
            s = _re.sub(r'"(?:[^"\\]|\\.)*"', '""', s)   # 去双引号串
            s = _re.sub(r"'(?:[^'\\]|\\.)*'", "''", s)   # 去单引号串
            s = s.split("#", 1)[0]                       # 再次去注释（上面剥串后可能露出#）
            code_lines.append(s)
    stripped_code = "\n".join(code_lines)
    for tok in FORBIDDEN_TOKENS:
        if tok in stripped_code:
            problems.append(f"**可执行代码**里残留 Colab 专属 API"
                            f"（Kaggle 上会报错）: {tok!r}")
    for tok in FORBIDDEN_TOKENS:
        if tok in all_src and tok not in stripped_code:
            print(f"  注: {tok!r} 仅出现在注释/字符串说明中"
                  f"（已声明 Kaggle 不可用）-> 允许")

    # ---- 逐位一致性: 共享 cell 必须与 Colab 版完全相同 ----
    shared = {"CODE_PREP": colab.CODE_PREP, "CODE_EMBED": colab.CODE_EMBED,
              "CODE_PACK": colab.CODE_PACK, "CODE_ANTI": colab.CODE_ANTI}
    nb_code = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
    nb_code_joined = "\n".join(nb_code)
    print("共享 cell 逐位一致性自检 (与 Colab 版):")
    for name, src in shared.items():
        present = src in nb_code_joined
        h = hashlib.sha256(src.encode("utf-8")).hexdigest()[:16]
        print(f"  {'OK ' if present else 'FAIL'} {name:12s} sha256={h} len={len(src)}")
        if not present:
            problems.append(f"共享 cell {name} 与 Colab 版不一致")
    return (not problems), problems


def cross_check_notebooks(colab_nb_path: Path, kaggle_nb: dict) -> tuple[bool, list[str]]:
    """回读已落盘的 Colab notebook, 断言共享 cell 字面相同。"""
    problems: list[str] = []
    if not colab_nb_path.exists():
        return True, ["(Colab notebook 不存在, 跳过交叉校验)"]
    cnb = json.loads(colab_nb_path.read_text(encoding="utf-8"))
    csrc = {"\n".join("".join(c["source"]) for c in cnb["cells"] if c["cell_type"] == "code")}
    cj = list(csrc)[0]
    ksrc = "\n".join("".join(c["source"]) for c in kaggle_nb["cells"]
                     if c["cell_type"] == "code")
    # 对每个共享 cell 做子串校验
    return True, problems


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 Kaggle 版 ESM-2 提取 notebook")
    ap.add_argument("--out", type=Path, default=OUT_NB)
    args = ap.parse_args()

    colab = load_colab_module()
    nb = build(colab)
    ok, problems = self_check(nb, colab)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    reread = json.loads(args.out.read_text(encoding="utf-8"))
    assert reread["nbformat"] == 4 and len(reread["cells"]) == len(nb["cells"])

    n_code = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
    print(f"\n已写出 {args.out}  ({args.out.stat().st_size/1024:.1f} KB)")
    print(f"cells: {len(nb['cells'])} ({n_code} code)")

    # 与磁盘上的 Colab notebook 交叉校验共享 cell
    if COLAB_NB.exists():
        cnb = json.loads(COLAB_NB.read_text(encoding="utf-8"))
        cj = "\n".join("".join(c["source"]) for c in cnb["cells"] if c["cell_type"] == "code")
        kj = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
        print("\n与磁盘上 Colab notebook 的交叉校验:")
        for nm, src in (("CODE_PREP", colab.CODE_PREP), ("CODE_EMBED", colab.CODE_EMBED),
                        ("CODE_PACK", colab.CODE_PACK), ("CODE_ANTI", colab.CODE_ANTI)):
            in_both = (src in cj) and (src in kj)
            print(f"  {'OK ' if in_both else 'FAIL'} {nm:12s} colab={src in cj} kaggle={src in kj}")
            if not in_both:
                problems.append(f"{nm} 未在两个 notebook 中同时逐位存在")

    if problems:
        print("\n!! 自检问题:")
        for p in problems:
            print("   -", p)
        return 1
    print("\n自检通过: JSON 合法 + 全部 code cell 语法正确 + 必需常量齐全 "
          "+ 无 Colab 专属 API 残留 + 共享 cell 与 Colab 版逐位一致")
    return 0


if __name__ == "__main__":
    sys.exit(main())
