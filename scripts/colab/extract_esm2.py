# -*- coding: utf-8 -*-
"""scripts/colab/extract_esm2.py -- 在免费 GPU (Colab / Kaggle) 上提取 ESM-2 嵌入.

**本脚本不属主流程依赖。** 主流程 (阶段 0-4) 全程纯 CPU、不 import torch。
这是刻意的设计: 提示词红线第 5 条要求「禁止引入 GPU 硬依赖」。

为什么要单独放这里
------------------
``esm2_t30_150M_UR50D`` 前向一次在 CPU 上可行但 3000+ 条序列很慢; 免费 GPU 上几分钟即可。
因此把嵌入**算一次、缓存下来**, 训练时只读 ``.npz``。

用法 (Colab / Kaggle)
---------------------
1. 把本仓库的 ``features/v2_features_enhanced.csv`` (或 ``features_v1.1.csv``)
   上传到 Colab, 或从 Google Drive 挂载。
2. 运行::

       !pip install -q transformers torch
       !python extract_esm2.py \
           --features /content/v2_features_enhanced.csv \
           --out /content/v2_esm2_cache.npz \
           --model facebook/esm2_t30_150M_UR50D \
           --pool both

   CPU 上也能跑通 (加 ``--limit 50`` 先冒烟)。
3. 把 ``v2_esm2_cache.npz`` 拷回本仓库的 ``features/`` 目录。
4. 本机执行 ``python scripts/18_extra_features.py``, 它会自动发现该缓存并按
   **训练集拟合** 的 PCA 降到 32 维拼进特征表 (不会泄漏测试集)。

输出格式 (与 18_extra_features.py 的约定)
----------------------------------------
``np.savez_compressed`` 保存:

* ``interface_key`` : ``<U`` 字符串数组, 长度 N
* ``embeddings``    : ``float32`` 矩阵 ``(N, D)``, **mean-pool 后**的链级嵌入
* ``model_name``    : 标量字符串
* ``pool``          : 标量字符串 (``mean`` / ``interface`` / ``both``)
* ``layer``         : 标量整数 (取第几层 hidden state; ``-1`` = 最后一层)

当 ``--pool both`` 时, 额外保存 ``embeddings_iface`` (界面残基 mean-pool), 由调用方选择。

诚实说明
--------
* 本脚本**未在本机运行** (本机无 GPU 且刻意不装 ESM-2 权重); 报告里标注为
  「Colab 手动步骤」。主流程缺该缓存时一切照常运行, 只是没有 ESM-2 维度。
* 嵌入的 PCA **必须**在训练集上 fit —— 本脚本不做任何 PCA, 只存原始嵌入,
  PCA 在 ``18_extra_features.py`` 里按划分正确执行。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="ESM-2 嵌入提取 (Colab/Kaggle GPU)")
    ap.add_argument("--features", required=True, type=Path,
                    help="含 interface_key / ligand_chains / receptor_chains 的特征表")
    ap.add_argument("--seq-cache", type=Path, default=None,
                    help="可选: data/pdb_seqs.json (链序列缓存)")
    ap.add_argument("--pdb-dir", type=Path, default=None,
                    help="可选: data/pdb 目录 (序列缓存链 ID 不匹配时从 ATOM 兜底)")
    ap.add_argument("--out", required=True, type=Path, help="输出 .npz 路径")
    ap.add_argument("--model", default="facebook/esm2_t30_150M_UR50D",
                    help="HF 模型名 (150M 是提示词首选; 8M 更快, 650M 更准)")
    ap.add_argument("--layer", type=int, default=-1, help="取第几层 hidden state, -1=最后")
    ap.add_argument("--pool", default="both", choices=("mean", "interface", "both"),
                    help="链级 mean-pool / 界面残基 mean-pool / 两者都存")
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--max-len", type=int, default=1022, help="ESM-2 上限 1022")
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 条 (冒烟)")
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}


def seq_from_structure(pdb_path: Path, chains: list[str]) -> dict[str, str]:
    """从 ATOM 记录提链序列 (缓存链 ID 不匹配时的兜底, 与 18/13 同逻辑)。"""
    want = set(chains)
    per: dict[str, list[tuple[int, str]]] = {c: [] for c in chains}
    seen: dict[str, set] = {c: set() for c in chains}
    if not pdb_path.exists():
        return {c: "" for c in chains}
    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6] != "ATOM  ":
                continue
            ch = line[21]
            if ch not in want:
                continue
            aa = AA3_TO_1.get(line[17:20].strip().upper())
            if aa is None:
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            key = (resseq, line[26])
            if key in seen[ch]:
                continue
            seen[ch].add(key)
            per[ch].append((resseq, aa))
    return {c: "".join(a for _r, a in sorted(v, key=lambda t: t[0]))
            for c, v in per.items()}


def build_chain_map(args: argparse.Namespace, df: pd.DataFrame) -> dict[str, str]:
    """为每个 interface_key 构造 "ligand_seq|receptor_seq" 拼接序列。"""
    cache: dict = {}
    if args.seq_cache and args.seq_cache.exists():
        cache = json.loads(args.seq_cache.read_text(encoding="utf-8"))
    out: dict[str, str] = {}
    for r in df.to_dict("records"):
        pid = str(r["pdb_id"]).upper()
        lig = [c.strip() for c in str(r["ligand_chains"]).split(",") if c.strip()]
        rec = [c.strip() for c in str(r["receptor_chains"]).split(",") if c.strip()]
        cs = cache.get(pid, {}) or {}
        parts = [cs.get(c, {}).get("seq", "") for c in lig + rec]
        if all(parts):
            out[r["interface_key"]] = "".join(parts)
            continue
        if args.pdb_dir:
            p = args.pdb_dir / f"{pid}.pdb"
            if not p.exists():
                p = args.pdb_dir / f"{pid}.ent"
            fb = seq_from_structure(p, lig + rec)
            merged = [a if len(a) >= len(fb.get(c, "")) else fb.get(c, "")
                      for a, c in zip(parts, lig + rec)]
            out[r["interface_key"]] = "".join(merged)
        else:
            out[r["interface_key"]] = "".join(parts)
    return out


def main() -> int:
    args = parse_args()
    t0 = time.time()

    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except Exception as exc:  # noqa: BLE001
        print(f"!! 需要 torch + transformers: {exc}")
        print("   pip install -q torch transformers")
        return 2

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("=" * 78)
    print("ESM-2 嵌入提取")
    print("=" * 78)
    print(f"模型   : {args.model}")
    print(f"设备   : {device}"
          + (f" ({torch.cuda.get_device_name(0)})" if device == "cuda" else " (CPU, 会很慢)"))
    print(f"层     : {args.layer}   pooling: {args.pool}")

    df = pd.read_csv(args.features, low_memory=False)
    if args.limit:
        df = df.head(args.limit).copy()
    print(f"输入   : {args.features.name} — {len(df)} 行")

    seq_map = build_chain_map(args, df)
    keys, seqs = [], []
    n_no_seq = 0
    for k, s in seq_map.items():
        if not s or len(s) < 5:
            n_no_seq += 1
            continue
        keys.append(k)
        seqs.append(s[: args.max_len])
    print(f"可用序列: {len(seqs)} / {len(df)} (无序列 {n_no_seq})")
    if not seqs:
        print("!! 没有任何可用序列")
        return 3

    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModel.from_pretrained(args.model).to(device).eval()

    embs: list[np.ndarray] = []
    t1 = time.time()
    with torch.no_grad():
        for i in range(0, len(seqs), args.batch_size):
            batch = seqs[i:i + args.batch_size]
            enc = tok(batch, return_tensors="pt", padding=True,
                      truncation=True, max_length=args.max_len).to(device)
            out = model(**enc, output_hidden_states=True)
            h = out.hidden_states[args.layer]            # (B, L, D)
            mask = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
            pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1e-6)
            embs.append(pooled.float().cpu().numpy())
            done = min(i + args.batch_size, len(seqs))
            if (i // args.batch_size) % 20 == 0 or done == len(seqs):
                el = time.time() - t1
                rate = done / max(el, 1e-6)
                eta = (len(seqs) - done) / max(rate, 1e-6)
                print(f"  {done}/{len(seqs)}  {rate:.1f} seq/s  ETA {eta / 60:.1f} min")

    E = np.vstack(embs).astype(np.float32)
    print(f"嵌入矩阵: {E.shape}")

    payload = {
        "interface_key": np.asarray(keys, dtype="U"),
        "embeddings": E,
        "model_name": np.asarray(args.model),
        "pool": np.asarray(args.pool),
        "layer": np.asarray(args.layer, dtype=np.int32),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **payload)
    print(f"写出: {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    print(f"耗时: {time.time() - t0:.1f}s")
    print("")
    print("下一步: 把该 npz 拷到本仓库 features/v2_esm2_cache.npz,")
    print("然后在本机运行: python scripts/18_extra_features.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
