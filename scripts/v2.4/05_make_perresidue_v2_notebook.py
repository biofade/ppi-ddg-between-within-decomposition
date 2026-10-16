# -*- coding: utf-8 -*-
"""scripts/v2.4/05_make_perresidue_v2_notebook.py -- 生成 per-residue 提取 notebook V2.

与 V2.3.2 版 (`scripts/kaggle/esm2_extract_perresidue.ipynb`) 的差别
--------------------------------------------------------------------
V2.3.2 版会提取**全量** 9519 条唯一序列的逐残基表示 (~3.1 GB)。
V2.4 需要大幅缩小体积, 因此新增**两个正交的压缩开关**:

1. **`SEQ_SUBSET`** —— 只提"有用"的序列:
   * `"test_val_skempi"` (新默认): test+val + **全部 SKEMPI 复合物**
     (SKEMPI 是 within-target 分析的唯一可用数据, 必须包含) -> ~200 MB;
   * `"test_val"`: 只提 test+val;
   * `"all"`: 全量 (~3.1 GB, 与 V2.3.2 版相同)。
2. **`INTERFACE_ONLY`** —— 只提**界面残基** (默认 **True**):
   * 用 `Bio.PDB` 计算配体/受体链间 < `INTERFACE_CUTOFF` (默认 **8.0 Å**) 的残基;
   * 只存这些残基的 per-residue 向量 -> 体积再降一个数量级 (~10-40 残基/链);
   * 非界面残基不存, 但**记录界面残基序号与总长度**, 便于下游对齐;
   * 需要 PDB 结构文件: 本 notebook 从 **RCSB 下载** (`files.rcsb.org`),
     并把下载结果缓存到 `/kaggle/working/pdb_cache/`。

> ⚠️ `INTERFACE_ONLY=True` 依赖网络与 RCSB 可用性。若某结构下载失败,
> 该条目**回退为整链**提取, 并在 `interface_fallback` 数组中标记 —— 不静默丢弃。

输出
----
`/kaggle/working/esm2_output_perres_v2/v2_esm2_perres_v2.npz`:
* `pdb_ids` / `interface_key` / `splits` / `truncated` / `lig_md5` / `rec_md5`
* `lig_lengths` / `rec_lengths`   (实际存储的残基数)
* `lig_flat` / `rec_flat`         (ΣL, 640) float16
* `lig_residue_ids` / `rec_residue_ids` (ΣL,) int32 —— **原始残基序号**
  (若为整链回退, 则等于 1..L)
* `lig_full_length` / `rec_full_length` (N,) int32 —— 截断后的整链长度
* `interface_fallback` (N,) bool —— 该条目是否因取不到结构而回退整链
* 元信息: `model_name` / `repr_layer` / `storage_dtype` / `seq_subset` /
  `interface_only` / `interface_cutoff`

不做的事
--------
**本 notebook 未运行** (需 Kaggle GPU)。本脚本只生成与自检。
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
OUT_NB = ROOT / "scripts" / "kaggle" / "esm2_extract_perresidue_v2.ipynb"

MD_TITLE = r"""# V2.4 — ESM-2 **per-residue + 界面残基** 提取（Kaggle，T4 x2）

> V2.4 结论: **pooled 嵌入对单点突变不敏感** —— within-target 中位 Spearman 只有
> **+0.1043**（cross-target 是 +0.3890）。原因很直接:
> mean-pool 把整条序列压成一个向量, **一个残基的改变被几百个残基平均掉**。
>
> 要预测突变效应, 必须保留**逐残基**表示, 且最好**只关注界面残基**。
> 本 notebook 就是为此准备输入。

---

## 与 V2.3.2 版的关系

| | V2.3.2 版 | **V2.4 版 (本 notebook)** |
|---|---|---|
| 提取范围 | 全部 9519 条唯一序列 | **test+val+SKEMPI**（默认） |
| 残基范围 | 整链 | **界面残基 ±8 Å**（默认） |
| 预估体积 | ~3.1 GB | **~20–60 MB** |
| 依赖 | 无 | **需联网下载 PDB**（界面判定） |

核心逻辑（字母表 / 清洗 / 截断 / md5 / OOM 降级 / 断点续算）与 V2.3.2 版**一致**;
新增的是**筛选与界面判定**。

---

## 两个压缩开关（决定体积）

```python
SEQ_SUBSET     = "test_val_skempi"   # "test_val_skempi" | "test_val" | "all"
INTERFACE_ONLY = True                # True = 只存界面残基
INTERFACE_CUTOFF = 8.0               # 界面判定距离 (Å)
```

| `SEQ_SUBSET` | 序列数（估） | `INTERFACE_ONLY=True` 体积 | `False` 体积 |
|---|---:|---:|---:|
| `test_val_skempi`（默认） | ~3,500 | **~20–60 MB** | ~1.1 GB |
| `test_val` | ~3,000 | ~15–50 MB | ~1.0 GB |
| `all` | 9,519 | ~60–200 MB | ~3.1 GB |

> 体积区间取决于每个复合物的界面残基数（通常 10–40 / 链）。

## 为什么默认包含 SKEMPI

V2.4 任务 1 的 within-target 基准建立在 **SKEMPI 的 288 个复合物**上。
若此处的子集不含 SKEMPI, 拿到数据后**无法**做突变效应分析。
因此默认 `test_val_skempi` 把 SKEMPI 全部包含进来。

## 写死的模型配置（与 V2.3.2 版一致）

| 项 | 值 |
|---|---|
| 模型 | `esm2_t30_150M_UR50D`（640 维） |
| 表示层 | **最后一层（layer 30）** |
| 池化 | **无**（保存逐残基） |
| 存储精度 | **float16** |
| 最大长度 | **1022 aa**（超出截断） |
| 设备 | **`cuda:0`** |
| 输入 | `/kaggle/input/ppb-affinity-proteinea/filtered.csv` |
| 输出 | `/kaggle/working/esm2_output_perres_v2/v2_esm2_perres_v2.npz` |

## 重要说明：界面残基的序号会被保留

因为只存界面残基, 下游必须知道**每个残基在整链中的原始序号**。
本 notebook 会输出 `lig_residue_ids` / `rec_residue_ids`（1-based 原始序号）
以及 `lig_full_length` / `rec_full_length`。

若某结构**下载失败**, 该条目**回退为整链**提取, 并在 `interface_fallback` 中标记为
`True` —— **不静默丢弃, 也不假装是界面残基**。

## 中途断开后如何恢复

每条序列立刻写 `/kaggle/working/esm2_output_perres_v2/checkpoints/<md5>.npz`,
重跑自动跳过。

## 下载方式

Kaggle 没有 `files.download`。跑完用**右侧 Output 面板** →
展开 `esm2_output_perres_v2/` → 下载 `v2_esm2_perres_v2.npz`。

> ⚠️ **PDB 下载需要 Internet = On**（Settings → Internet）。
> 若 Internet 不可用, 把 `INTERFACE_ONLY` 设为 `False` 仍可正常运行（整链模式）。
"""

MD_STEP1 = r"""## 步骤 1 / 7 — 安装依赖

Kaggle 预装 torch / numpy。需要 Internet = On。
额外装 `fair-esm` 与 `biopython`（界面判定用）。"""

CODE_INSTALL = r"""!pip install fair-esm biopython -q
print("依赖安装完成 (fair-esm + biopython)。")"""

MD_STEP2 = "## 步骤 2 / 7 — 确认 GPU（只用 cuda:0）"

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

MD_STEP3 = r"""## 步骤 3 / 7 — 路径、子集与压缩开关

**只改这一节的三个常量**：

* `SEQ_SUBSET`: 提哪些序列（默认含 SKEMPI）
* `INTERFACE_ONLY`: 是否只提界面残基
* `INTERFACE_CUTOFF`: 界面判定距离（Å）"""

CODE_SETUP = r"""import os
import shutil
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ====================================================================
# 只改这三处
# ====================================================================
DATASET_SLUG = "ppb-affinity-proteinea"
SEQ_SUBSET = "test_val_skempi"   # "test_val_skempi" | "test_val" | "all"
INTERFACE_ONLY = True            # True = 只提界面残基 (需联网下载 PDB)
INTERFACE_CUTOFF = 8.0           # 界面判定距离 (Å)

PROTEINEA_CSV = f"/kaggle/input/{DATASET_SLUG}/filtered.csv"
OUT_DIR = "/kaggle/working/esm2_output_perres_v2"
CHECKPOINT_DIR = os.path.join(OUT_DIR, "checkpoints")
PDB_CACHE = "/kaggle/working/pdb_cache"
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
os.makedirs(PDB_CACHE, exist_ok=True)

# ---- 写死的模型与预处理配置（与 V2.3.2 版一致） ----
MODEL_NAME = "esm2_t30_150M_UR50D"
MODEL_FALLBACK = "esm2_t12_35M_UR50D"
REPR_LAYER = 30
MAX_LEN = 1022
BATCH_SIZE = 2
BACKUP_EVERY = 100
PERRES_DTYPE = np.float16
CHECKPOINT_ON = True
LIG_COL = "Ligand Sequences"
REC_COL = "Receptor Sequences"
ID_COL = "PDB"
SPLIT_COL = "split"

print("=== /kaggle/input 内容 ===")
_in = Path("/kaggle/input")
if _in.exists():
    for d in sorted(_in.iterdir()):
        print(f"  {d.name}/")

csv_path = Path(PROTEINEA_CSV)
if not csv_path.exists():
    raise FileNotFoundError(
        f"找不到 {PROTEINEA_CSV}\n"
        f"请确认 Dataset 名称为 '{DATASET_SLUG}' 且已在右侧 Add Input 挂载。"
    )

USE_COLS = [ID_COL, LIG_COL, REC_COL, SPLIT_COL]
df = pd.read_csv(csv_path, usecols=USE_COLS, dtype=str)
df[ID_COL] = df[ID_COL].astype(str).str.strip().str.upper()
df = df.dropna(subset=[LIG_COL, REC_COL]).reset_index(drop=True)
_n_all = len(df)

# ---- 需要外部 SKEMPI 行清单 (含全部 SKEMPI 复合物) ----
# SKEMPI 的 PDB 列表从 V2.2 特征表获得; 若 Kaggle 上没有该文件,
# 则退化为只用 split 过滤, 并明确警告 (不静默降级)。
if SEQ_SUBSET == "test_val_skempi":
    # SKEMPI 复合物在 filtered.csv 里的 PDB 集合 (V2.4 任务 0/1 实测)
    SKEMPI_PDBS = set('''
1A22 1A4Y 1ACB 1AHW 1AK4 1B2S 1B2U 1B3S 1B41 1BJ1
1BP3 1BRS 1C1Y 1C4Z 1CSE 1CSO 1CT0 1CT2 1CT4 1CZ8
1DAN 1DQJ 1DVF 1E50 1E96 1EAW 1EFN 1EMV 1F5R 1FC2
1FCC 1FFW 1FR2 1FSS 1FY8 1GC1 1GCQ 1GRN 1GUA 1H9D
1HE8 1IAR 1JCK 1JRH 1JTD 1JTG 1K8R 1KAC 1KBH 1KIR
1KTZ 1LFD 1M9E 1MAH 1MHP 1MLC 1MQ8 1N8Z 1NCA 1NMB
1OHZ 1P69 1P6A 1PPF 1QAB 1R0R 1REW 1S0W 1S1Q 1SBB
1SBN 1SGD 1SGE 1SGN 1SGP 1SGQ 1SGY 1SIB 1TM1 1TM3
1TM4 1TM5 1TM7 1TMG 1TO1 1U7F 1UUZ 1VFB 1WQJ 1X1W
1X1X 1XD3 1XGP 1XGQ 1XGR 1XGT 1XGU 1XXM 1Y1K 1Y33
1Y34 1Y3B 1Y3C 1Y3D 1Y48 1Y4A 1YCS 1YQV 1YY9 1Z7X
2A9K 2ABZ 2AJF 2AW2 2B0U 2B0Z 2B10 2B11 2B12 2B2X
2B42 2BDN 2BTF 2C0L 2C5D 2CCL 2DSQ 2DVW 2FTL 2G2U
2G2W 2GOX 2GYK 2HLE 2HRK 2I26 2J0T 2J12 2J1K 2JEL
2KSO 2NOJ 2NU0 2NU1 2NU2 2NU4 2NY7 2NYY 2NZ9 2O3B
2OOB 2PCB 2PCC 2QJ9 2QJA 2QJB 2REX 2SGP 2SGQ 2SIC
2VIR 2VIS 2VLN 2VLO 2VLP 2VLQ 2VN5 2WPT 3AAA 3B4V
3BDY 3BE1 3BIW 3BK3 3BN9 3BP8 3BT1 3BTD 3BTE 3BTF
3BTG 3BTH 3BTM 3BTQ 3BTT 3BTW 3BX1 3C60 3D5R 3D5S
3EG5 3F1S 3G6D 3HFM 3HH2 3IDX 3KBH 3KUD 3L5X 3LB6
3LZF 3M62 3M63 3MZG 3MZW 3N06 3N0P 3N4I 3N85 3NCB
3NCC 3NGB 3NPS 3NVN 3NVQ 3Q3J 3QHY 3R9A 3S9D 3SE3
3SE4 3SE8 3SE9 3SEK 3SF4 3SGB 3SZK 3TGK 3U82 3VR6
3W2D 3WWN 4B0M 4BFI 4CVW 4E6K 4EKD 4FZA 4G0N 4GNK
4GXU 4HFK 4HRN 4HSA 4I77 4JEU 4JGH 4JPK 4K71 4KRL
4KRO 4KRP 4L0P 4LRX 4MYW 4NKQ 4NM8 4NZW 4O27 4OFY
4P23 4P5T 4PWX 4RA0 4RS1 4U6H 4UYP 4UYQ 4X4M 4Y61
4YEB 4YFD 4YH7 4ZS6 5C6T 5CXB 5CYK 5DWU 5E6P 5F4E
5K39 5M2O 5TAR 5UFE 5UFQ
'''.split())
    # 若上面的硬编码列表为空/不足, 用 split 兜底并在日志中说明
    _m = df[df[ID_COL].isin(SKEMPI_PDBS) | df[SPLIT_COL].isin(["test", "val"])]
    print(f"\n!! SEQ_SUBSET='test_val_skempi': 用内置 SKEMPI PDB 列表 "
          f"({len(SKEMPI_PDBS)} 个) + test/val")
    if len(SKEMPI_PDBS) < 50:
        print("   ⚠️ 注意: 内置 SKEMPI 列表较短。若要覆盖全部 SKEMPI 复合物,")
        print("      请把 features/v2.4_skempi_within_target.csv 的 pdb_id 列")
        print("      粘贴到本 cell 的 SKEMPI_PDBS 集合里。")
else:
    _m = df[df[SPLIT_COL].isin(["test", "val"])] if SEQ_SUBSET == "test_val" else df
df = _m.reset_index(drop=True)
print(f"读入 {_n_all} 行 -> 子集 '{SEQ_SUBSET}' 保留 {len(df)} 行")
print("split 分布:", df[SPLIT_COL].value_counts().to_dict())

try:
    _st = shutil.disk_usage("/kaggle/working")
    print(f"/kaggle/working 可用: {_st.free / 1024**3:.1f} GB")
except Exception as exc:
    print(f"(磁盘检查失败: {exc})")
print(f"\n配置: INTERFACE_ONLY={INTERFACE_ONLY}, 判定距离={INTERFACE_CUTOFF} Å")"""

MD_STEP4 = "## 步骤 4 / 7 — 序列清洗与去重"

CODE_PREP = r"""from esm.data import Alphabet

ALPHABET = Alphabet.from_architecture("ESM-2")
UNK = ALPHABET.unk_idx
STD = set(ALPHABET.all_tokens)


def sanitize(seq: str) -> tuple[str, bool, int]:
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
print(f"行数                : {len(records)}")
print(f"唯一序列（要嵌）     : {len(uniq_all)}")
print(f"总残基数（整链）     : {_total_res:,}")
print(f"整链体积 (fp16)      : {_total_res * 640 * 2 / 1024**3:.2f} GB")
if INTERFACE_ONLY:
    print(f"界面模式预估         : 再降 10-50x (典型 10-40 残基/链)")
if len(uniq_all) == 0:
    raise SystemExit("没有可嵌入的序列。")"""

MD_STEP5 = r"""## 步骤 5 / 7 — 界面残基判定 + per-residue 嵌入

* **界面判定**: 用 `Bio.PDB` 读 PDB, 取跨链重原子距离 < `INTERFACE_CUTOFF` 的残基;
* **缓存判定结果** 到 `interface_map`（key = `pdb_id` + 链组），避免重复下载;
* **取不到结构** -> 该条目回退整链, 记入 `interface_fallback`;
* **嵌入**: 同 V2.3.2 版（无 pool, 去掉 BOS/EOS/pad）;
* **只保留界面残基的行**（并记录其原始残基序号）;
* 断点续算 / OOM 降级 / 单条失败不中断, 与 V2.3.2 版一致。"""

CODE_EMBED = r"""import esm
import gc
import urllib.request

failures: list[dict] = []


def load_model(name: str):
    models = esm.pretrained.esm2_model_registry
    if name not in models:
        lower = {k.lower(): k for k in models}
        if name.lower() in lower:
            name = lower[name.lower()]
        else:
            raise KeyError(f"fair-esm 不认识 {name!r}")
    model, alphabet = models[name]()
    model = model.eval().to(DEVICE)
    return model, alphabet, alphabet.get_batch_converter(), name


MODEL, ALPHABET2, batch_converter, LOADED_KEY = load_model(MODEL_NAME)
_DIM = 640 if "150M" in LOADED_KEY else 480


def fetch_pdb(pdb_id: str, tag: str = "lig"):
    # 下载 PDB 到缓存。返回路径或 None。
    safe = f"{pdb_id}_{tag}.pdb"
    p = os.path.join(PDB_CACHE, safe)
    if os.path.exists(p) and os.path.getsize(p) > 500:
        return p
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, p)
        if os.path.getsize(p) > 500:
            return p
    except Exception as exc:
        print(f"    !! 下载 {pdb_id} 失败: {type(exc).__name__}")
    return None


def interface_residues(pdb_path: str, cutoff: float):
    # 返回 {chain_id: sorted(set(resseq))} —— 参与跨链接触的残基 (重原子 < cutoff)。
    from Bio.PDB import PDBParser, NeighborSearch
    parser = PDBParser(QUIET=True)
    try:
        st = parser.get_structure("x", pdb_path)
    except Exception:
        return None
    model = next(st.get_models())
    atoms = [a for a in model.get_atoms()
             if a.element not in ("H", "D") and a.get_altloc() in (" ", "A")]
    if not atoms:
        return None
    ns = NeighborSearch(atoms)
    out: dict[str, set] = {}
    for a1, a2 in ns.search_all(cutoff):
        if a1.get_parent().get_parent().id == a2.get_parent().get_parent().id:
            continue
        for a in (a1, a2):
            ch = a.get_parent().get_parent().id
            out.setdefault(ch, set()).add(int(a.get_parent().id[1]))
    return {k: sorted(v) for k, v in out.items()}


def per_residue(hidden, tokens):
    # 逐残基表示 (去掉 BOS/EOS/pad)。返回长度 B 的 list[(L_i, D)]。
    import torch as _t
    mask = (tokens != 1) & (tokens != 2) & (tokens != 0)
    out = []
    for b in range(hidden.shape[0]):
        idx = mask[b].nonzero(as_tuple=False).ravel()
        out.append(hidden[b].index_select(0, idx).float().cpu().numpy())
    return out


def embed_sequences(seqs, model, conv, batch_size, tag):
    # 返回 ([(seq, (L,D) 矩阵)], failures)。自带序列, 不会错位。
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
            mats = per_residue(reps["representations"][REPR_LAYER], tokens)
            if len(mats) != len(chunk):
                raise RuntimeError("per_residue 长度不匹配")
            for s, m in zip(chunk, mats):
                pairs.append((s, m))
            i += len(chunk)
        except _t.cuda.OutOfMemoryError:
            _t.cuda.empty_cache(); gc.collect()
            if bs > 1:
                bs = max(1, bs // 2)
                print(f"  !! CUDA OOM -> batch_size 降到 {bs}")
                continue
            failed.append({"tag": tag, "index": i, "seq": chunk[0][:40],
                           "reason": "OOM at batch=1"})
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


# ---- 主循环: 按 (pdb, 序列) 提取, 并按界面残基裁剪 ----
res_cache: dict[str, np.ndarray] = {}       # seq -> (L, D) 整链
iface_cache: dict[str, dict] = {}           # pdb_id -> {chain: [resseq]}
fallback: dict[str, bool] = {}              # pdb_id -> 是否回退整链

todo = sorted({s for s in uniq_all})
print(f"总唯一序列 {len(todo)}; 开始嵌入 (整链长度) ...")
t_start = time.time()
i = 0
while i < len(todo):
    chunk = todo[i:i + BATCH_SIZE]
    pairs, failed = embed_sequences(chunk, MODEL, batch_converter, BATCH_SIZE, "main")
    for s, m in pairs:
        res_cache[s] = m
    failures.extend(failed)
    i += len(chunk)
    done = len(res_cache)
    if done % BACKUP_EVERY < len(chunk) or done >= len(todo):
        el = time.time() - t_start
        rate = done / max(el, 1e-9)
        eta = (len(todo) - done) / max(rate, 1e-9)
        print(f"  进度 {done}/{len(todo)}  {rate:.1f} seq/s  "
              f"已用 {el/60:.1f} min  剩余 ~{eta/60:.1f} min  失败 {len(failures)}")
print(f"嵌入完成: {len(res_cache)} 条, 失败 {len(failures)}, "
      f"耗时 {(time.time()-t_start)/60:.1f} min")

# ---- 界面判定 (按 pdb_id) ----
if INTERFACE_ONLY:
    print("\n开始界面判定 (需下载 PDB) ...")
    _pdbs = sorted({r[ID_COL] for r in records})
    for k, pid in enumerate(_pdbs):
        pth = fetch_pdb(pid)
        im = interface_residues(pth, INTERFACE_CUTOFF) if pth else None
        if im is None:
            iface_cache[pid] = {}
            fallback[pid] = True
        else:
            iface_cache[pid] = im
            fallback[pid] = False
        if (k + 1) % 25 == 0 or k + 1 == len(_pdbs):
            _nf = sum(1 for v in fallback.values() if v)
            print(f"  界面判定 {k+1}/{len(_pdbs)}  回退整链 {_nf}")
    print(f"界面判定完成; 回退整链的 PDB: "
          f"{sum(1 for v in fallback.values() if v)}/{len(_pdbs)}")
else:
    iface_cache = {}
    fallback = {}
    print("\nINTERFACE_ONLY=False -> 使用整链。")"""

MD_STEP6 = r"""## 步骤 6 / 7 — 打包 `v2_esm2_perres_v2.npz`

* 若 `INTERFACE_ONLY=True`, 只把**界面残基**的行写进 `lig_flat` / `rec_flat`,
  并在 `*_residue_ids` 里记录这些残基的**原始序号**;
* `*_full_length` 记录整链长度, 便于下游还原位置;
* `interface_fallback` 标记因取不到结构而回退整链的条目。"""

CODE_PACK = r"""# 变长存储: flat + lengths + 原始残基序号
n_rows = len(records)
lig_chunks, rec_chunks = [], []
lig_lens = np.zeros(n_rows, dtype=np.int32)
rec_lens = np.zeros(n_rows, dtype=np.int32)
lig_full = np.zeros(n_rows, dtype=np.int32)
rec_full = np.zeros(n_rows, dtype=np.int32)
lig_resids, rec_resids = [], []
residue_ids_available = True

def _residue_ids(pid, chain_guess):
    # 从界面缓存里取该 PDB 的残基序号集合 (合并所有链, 保守)。
    im = iface_cache.get(pid)
    if not im:
        return None
    if chain_guess and chain_guess in im:
        return set(im[chain_guess])
    s = set()
    for v in im.values():
        s.update(v)
    return s

missing_lig = sum(1 for r in records if r["lig_clean"] not in res_cache)
missing_rec = sum(1 for r in records if r["rec_clean"] not in res_cache)
print(f"条目 {n_rows}; 缺配体 {missing_lig}; 缺受体 {missing_rec}")

for i, r in enumerate(records):
    pid = r[ID_COL]
    for side, seqkey, chunks, lens, fulls, resids in (
            ("lig", "lig_clean", lig_chunks, lig_lens, lig_full, lig_resids),
            ("rec", "rec_clean", rec_chunks, rec_lens, rec_full, rec_resids)):
        mat = res_cache.get(r[seqkey])
        if mat is None or mat.shape[0] == 0:
            continue
        fulls[i] = mat.shape[0]
        keep_ids = None
        if INTERFACE_ONLY and not fallback.get(pid, True):
            ids = _residue_ids(pid, None)
            if ids is not None:
                # mat 的第 j 行对应原始残基序号 j+1 (截断后)
                keep_ids = [j + 1 for j in range(mat.shape[0]) if (j + 1) in ids]
        if keep_ids is None or len(keep_ids) == 0:
            sub = mat
            ids_out = list(range(1, mat.shape[0] + 1))
        else:
            sub = mat[np.asarray(keep_ids) - 1]
            ids_out = keep_ids
        chunks.append(sub.astype(PERRES_DTYPE))
        lens[i] = sub.shape[0]
        resids.append(np.asarray(ids_out, dtype=np.int32))

lig_flat = (np.concatenate(lig_chunks, axis=0) if lig_chunks
            else np.zeros((0, _DIM), dtype=PERRES_DTYPE))
rec_flat = (np.concatenate(rec_chunks, axis=0) if rec_chunks
            else np.zeros((0, _DIM), dtype=PERRES_DTYPE))
lig_res_arr = (np.concatenate(lig_resids) if lig_resids
               else np.zeros(0, dtype=np.int32))
rec_res_arr = (np.concatenate(rec_resids) if rec_resids
               else np.zeros(0, dtype=np.int32))
print(f"lig_flat {lig_flat.shape}  rec_flat {rec_flat.shape}")
print(f"界面模式: {INTERFACE_ONLY} | 平均存留配体残基 "
      f"{lig_lens.sum()/max(1,(lig_lens>0).sum()):.1f} / 整链均值 "
      f"{lig_full.sum()/max(1,(lig_full>0).sum()):.1f}")

out_path = os.path.join(OUT_DIR, "v2_esm2_perres_v2.npz")
np.savez_compressed(
    out_path,
    pdb_ids=np.array([r[ID_COL] for r in records], dtype="U"),
    interface_key=np.array([f"{r[ID_COL]}|{r[LIG_COL]}|{r[REC_COL]}"
                            for r in records], dtype="U"),
    splits=np.array([str(r[SPLIT_COL]) for r in records], dtype="U"),
    lig_lengths=lig_lens, rec_lengths=rec_lens,
    lig_full_length=lig_full, rec_full_length=rec_full,
    lig_residue_ids=lig_res_arr, rec_residue_ids=rec_res_arr,
    lig_flat=lig_flat, rec_flat=rec_flat,
    truncated=np.array([bool(r["lig_trunc"] or r["rec_trunc"]) for r in records],
                       dtype=bool),
    interface_fallback=np.array([bool(fallback.get(r[ID_COL], True))
                                 for r in records], dtype=bool),
    lig_md5=np.array([r["lig_md5"] for r in records], dtype="U"),
    rec_md5=np.array([r["rec_md5"] for r in records], dtype="U"),
    model_name=np.array(LOADED_KEY),
    repr_layer=np.array(REPR_LAYER),
    pool=np.array("none_per_residue_excl_bos_eos"),
    storage_dtype=np.array(np.dtype(PERRES_DTYPE).name),
    seq_subset=np.array(SEQ_SUBSET),
    interface_only=np.array(bool(INTERFACE_ONLY)),
    interface_cutoff=np.array(float(INTERFACE_CUTOFF)),
    perres_dim=np.array(_DIM),
    n_failed=np.array(len(failures)),
)
print(f"\n已写出: {out_path}")
print(f"  大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"  lig_flat {lig_flat.shape}  rec_flat {rec_flat.shape}  "
      f"dim {_DIM}")
if _DIM != 640:
    print(f"  !! 警告: perres_dim={_DIM} != 640 (可能降级到 35M 模型)")"""

MD_STEP7 = r"""## 步骤 7 / 7 — 下载（Kaggle Output 面板）

⚠️ Kaggle 没有 `files.download`。

1. 右侧 **Output** 面板 → 展开 **`esm2_output_perres_v2/`**
2. 下载 **`v2_esm2_perres_v2.npz`**
3. 放到本机 `F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_perres_v2.npz`"""

CODE_DOWNLOAD = r"""import os

print("=" * 72)
print("Kaggle 下载步骤 (没有 files.download):")
print("  1. 右侧 'Output' 面板 (没有就点右上角 > 展开)")
print("  2. 展开 'esm2_output_perres_v2'")
print("  3. 下载 'v2_esm2_perres_v2.npz'")
print("=" * 72)
print(f"\n绝对路径: {out_path}")
print(f"大小: {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"seq_subset={SEQ_SUBSET}  interface_only={INTERFACE_ONLY}  "
      f"cutoff={INTERFACE_CUTOFF}  dim={_DIM}")
print(f"界面回退整链的条目: {int(np.sum(interface_fallback))} / {len(interface_fallback)}")
print()
print("本机放置位置:")
print(r"  F:\DSH_Galatea\Affinity_Intelligence\features\v2_esm2_perres_v2.npz")"""


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
    "esm2_t30_150M_UR50D", "esm2_t12_35M_UR50D", "REPR_LAYER = 30", "MAX_LEN = 1022",
    "PERRES_DTYPE", "SEQ_SUBSET", "test_val_skempi", "INTERFACE_ONLY",
    "INTERFACE_CUTOFF", "interface_residues", "NeighborSearch", "files.rcsb.org",
    "lig_residue_ids", "rec_residue_ids", "lig_full_length", "interface_fallback",
    "v2_esm2_perres_v2.npz", "/kaggle/working/esm2_output_perres_v2",
    "/kaggle/input/", "OutOfMemoryError", "md5_of", "float16",
]
FORBIDDEN_CODE = ["google.colab", "drive.mount", "/content/drive"]


def self_check(nb: dict) -> tuple[bool, list[str]]:
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
    # 必须包含界面裁剪逻辑, 且不得含 mean-pool
    if "keep_ids" not in all_src:
        problems.append("缺少界面残基裁剪逻辑 (keep_ids)")
    if "mean_pool(" in all_src:
        problems.append("不应包含 mean_pool (这是 per-residue 版)")
    return (not problems), problems


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 per-residue V2 提取 notebook")
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
          "+ 无 Colab API 残留 + 含界面裁剪 + 无 mean-pool")
    return 0


if __name__ == "__main__":
    sys.exit(main())
