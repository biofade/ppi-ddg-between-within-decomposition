# -*- coding: utf-8 -*-
"""scripts/v2.5/01_build_perres_cache.py -- V2.5 任务 1.1: 构建可逐行对齐的紧凑缓存.

为什么需要这一步
----------------
per-residue npz 是 **1.9 GB** (lig_flat 867502x640 fp16 + rec_flat 749687x640 fp16),
而本机只有 7.5 GB 内存且已被系统占用大半 (实测空闲 ~1.5 GB)。
因此**不能**把整个 per-residue 矩阵载入内存来训练 CNN。

同时必须解决**逐行对齐**问题 (任务 1 判决实验的前提):

* per-residue 缓存 6885 行 (train 5163 / val 965 / test 757);
* V2.2 特征表 8207 行 (train 6485 / val 965 / test 757);
* 两者 `interface_key` **格式不同** (sequence-based vs chain-based), 交集 0;
* 但用 `(pdb_id, 组内出现次序)` 可 100% 定位, 且 split 一致率 100%
  (`logs/v2.5/join_check.log`)。

本脚本做三件事
---------------
1. **建立映射**: perres row -> V2.2 row -> pooled(ESM-2 缓存) row,
   并用 **`lig_md5`/`rec_md5` 断言**映射正确 (两个缓存由同一份 filtered.csv 生成,
   同一行的 md5 必须相同)。**md5 不一致就报错停止** —— 绝不用错位的嵌入。
2. **构建残基下采样缓存**: 每个条目保留 **每链最多 `--topk` 个残基**
   (默认 48), 选择依据是 **ESM-2 上一条目内的残基显著性分数**
   = `||h_i - mean_j h_j||_2` (离该链平均表示越远的残基越可能是功能位点)。
   这是**确定性**规则, 不使用任何标签。
3. 落盘紧凑 npz (`v2.5_perres_cache.npz`), 体积从 1.9 GB 降到 ~几十 MB,
   使 CNN 可在内存内训练。

> ⚠️ 下采样是**为算力所迫的近似**, 会在报告中明确标注:
> 它丢掉了部分界面残基, 因此 CNN 的结果是 per-residue 能力的**保守下界**。

输出
----
* `features/v2.5_perres_cache.npz`  (紧凑 per-residue, 逐行对齐 V2.2)
* `logs/v2.5/build_cache.log`
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
LOG = ROOT / "logs" / "v2.5"
LOG.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG / "failures.log"

PERRES = FEAT / "v2_esm2_perres_v2.npz"
POOLED = FEAT / "v2_esm2_cache.npz"
V22 = FEAT / "v2.2_proteinea_features.csv"
OUT = FEAT / "v2.5_perres_cache.npz"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 构建紧凑 per-residue 缓存")
    ap.add_argument("--topk", type=int, default=48,
                    help="每链保留的残基数上限 (按显著性分数降序)")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    t0 = time.time()
    emit("=== V2.5 任务 1.1: 构建逐行对齐的紧凑 per-residue 缓存 ===")
    emit("")

    z = np.load(PERRES, allow_pickle=True, mmap_mode="r")
    pr_pdb = np.asarray(z["pdb_ids"]).astype(str)
    pr_sp = np.asarray(z["splits"]).astype(str)
    pr_lig_len = np.asarray(z["lig_lengths"])
    pr_rec_len = np.asarray(z["rec_lengths"])
    pr_lig_md5 = np.asarray(z["lig_md5"]).astype(str)
    pr_rec_md5 = np.asarray(z["rec_md5"]).astype(str)
    n_pr = len(pr_pdb)
    emit(f"per-residue 缓存: {n_pr} 行 (mmap, 不载入内存)")

    d = pd.read_csv(V22, low_memory=False)
    emit(f"V2.2 特征表: {len(d)} 行")

    # ---------- 1. 映射 perres -> V2.2 ----------
    emit("")
    emit("--- 1) perres -> V2.2 逐行映射 `(pdb_id, 组内次序)` ---")
    a = d.copy()
    a["_seq"] = a.groupby("pdb_id").cumcount()
    b = pd.DataFrame({"pdb_id": pr_pdb})
    b["_seq"] = b.groupby("pdb_id").cumcount()
    pos = a.set_index(["pdb_id", "_seq"]).index
    loc22 = pos.get_indexer(pd.MultiIndex.from_arrays([b["pdb_id"], b["_seq"]]))
    if (loc22 < 0).any():
        emit(f"!! 有 {int((loc22 < 0).sum())} 行无法定位, 停止")
        (LOG / "build_cache.log").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3
    sp22 = d["split"].astype(str).to_numpy()[loc22]
    agree = float((sp22 == pr_sp).mean())
    emit(f"定位成功 {len(loc22)}/{n_pr}; split 一致率 **{agree * 100:.2f}%**")
    if agree < 1.0:
        emit("!! split 不一致, 停止 (映射不可信)")
        (LOG / "build_cache.log").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 4

    # ---------- 2. 映射 V2.2 -> pooled 缓存, 用 md5 断言 ----------
    emit("")
    emit("--- 2) V2.2 -> pooled(ESM-2) 逐行映射, 用 md5 断言 ---")
    po = np.load(POOLED, allow_pickle=True)
    po_pdb = po["pdb_ids"].astype(str)
    po_lig_md5 = po["lig_md5"].astype(str)
    po_rec_md5 = po["rec_md5"].astype(str)
    po_lig = po["lig_l2"].astype(np.float32)      # 用 l2 (V2.3.2 最优变体)
    po_rec = po["rec_l2"].astype(np.float32)
    emit(f"pooled 缓存: {len(po_pdb)} 行")
    # pooled 与 V2.2 全表同序 (V2.3.2 已验证), 但仍用 md5 复核
    pa = pd.DataFrame({"pdb_id": po_pdb})
    pa["_seq"] = pa.groupby("pdb_id").cumcount()
    pos2 = pa.set_index(["pdb_id", "_seq"]).index
    loc_po = pos2.get_indexer(pd.MultiIndex.from_arrays([b["pdb_id"], b["_seq"]]))
    if (loc_po < 0).any():
        emit(f"!! pooled 定位失败 {int((loc_po < 0).sum())} 行")
        (LOG / "build_cache.log").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 5
    okm = (po_lig_md5[loc_po] == pr_lig_md5) & (po_rec_md5[loc_po] == pr_rec_md5)
    emit(f"md5 一致的行: **{int(okm.sum())} / {n_pr}** ({okm.mean() * 100:.2f}%)")
    if not okm.all():
        nbad = int((~okm).sum())
        emit(f"!! **{nbad} 行 md5 不一致 -> 映射错位, 停止** (绝不用错位嵌入)")
        bad = np.where(~okm)[0][:5]
        for i in bad:
            emit(f"   pdb={pr_pdb[i]} lig_md5 perres={pr_lig_md5[i][:8]} "
                 f"pooled={po_lig_md5[loc_po[i]][:8]}")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[01_build_perres_cache] md5 mismatch {nbad} rows\n")
        (LOG / "build_cache.log").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 6
    emit("-> ✅ 映射经 md5 逐行确认, 无错位")
    # 额外健康检查: md5 是否真的区分不同突变体
    emit(f"   (健康检查) perres 内唯一 lig_md5: {len(np.unique(pr_lig_md5))} / {n_pr}; "
         f"唯一 rec_md5: {len(np.unique(pr_rec_md5))} / {n_pr}")

    # ---------- 3. 残基下采样 ----------
    K = args.topk
    emit("")
    emit(f"--- 3) 残基下采样 (每链最多 {K} 个, 按残基显著性) ---")
    emit("    显著性 = ||h_i - mean_j h_j||_2  (同一链内, 确定性, 不用标签)")
    lig_flat = z["lig_flat"]          # mmap, 分片读取
    rec_flat = z["rec_flat"]
    lig_off = np.concatenate([[0], np.cumsum(pr_lig_len)])
    rec_off = np.concatenate([[0], np.cumsum(pr_rec_len)])

    Lk = np.zeros((n_pr, K, 640), dtype=np.float16)
    Rk = np.zeros((n_pr, K, 640), dtype=np.float16)
    Lmask = np.zeros((n_pr, K), dtype=bool)
    Rmask = np.zeros((n_pr, K), dtype=bool)
    t1 = time.time()
    for i in range(n_pr):
        for flat, off, lens, out, mask in (
                (lig_flat, lig_off, pr_lig_len, Lk, Lmask),
                (rec_flat, rec_off, pr_rec_len, Rk, Rmask)):
            n_i = int(lens[i])
            if n_i <= 0:
                continue
            block = np.asarray(flat[off[i]:off[i + 1]], dtype=np.float32)
            if n_i > K:
                mu = block.mean(axis=0, keepdims=True)
                score = np.linalg.norm(block - mu, axis=1)
                sel = np.argsort(-score, kind="mergesort")[:K]
            else:
                sel = np.arange(n_i)
            out[i, :len(sel)] = block[sel].astype(np.float16)
            mask[i, :len(sel)] = True
        if (i + 1) % 1000 == 0:
            el = time.time() - t1
            emit(f"    下采样 {i+1}/{n_pr}  ({el:.0f}s, "
                 f"预计剩余 {el / (i + 1) * (n_pr - i - 1):.0f}s)")
    emit(f"    下采样完成, 耗时 {time.time() - t1:.0f}s")
    emit(f"    保留残基数: 配体 均值 {Lmask.sum(1).mean():.1f} (原 {pr_lig_len.mean():.1f}); "
         f"受体 均值 {Rmask.sum(1).mean():.1f} (原 {pr_rec_len.mean():.1f})")

    # ---------- 4. 落盘 ----------
    #
    # ⚠️ **不能用 `np.savez_compressed`**: 实测在本机 (7.5 GB 内存, 空闲 ~2 GB)
    # 会抛 `MemoryError: Unable to allocate output buffer`
    # —— 压缩需要一次性申请整个数组的输出缓冲。
    # 改为**不压缩** + 全部 fp16 (体积仍从 1.9 GB 降到 ~0.4 GB)。
    tmp = args.out.with_suffix(args.out.suffix + ".tmp.npz")
    np.savez(
        tmp,
        pdb_ids=pr_pdb, splits=pr_sp,
        v22_row=loc22.astype(np.int64),
        pooled_row=loc_po.astype(np.int64),
        lig_topk=Lk, rec_topk=Rk, lig_mask=Lmask, rec_mask=Rmask,
        lig_pooled=po_lig[loc_po].astype(np.float16),
        rec_pooled=po_rec[loc_po].astype(np.float16),
        lig_len_full=pr_lig_len, rec_len_full=pr_rec_len,
        topk=np.array(K), source=np.array("v2_esm2_perres_v2.npz"),
    )
    tmp.replace(args.out)          # 原子替换, 避免半写文件被误读
    mb = args.out.stat().st_size / 1e6
    emit("")
    emit(f"写出 `{args.out.name}`: **{mb:.1f} MB** "
         f"(原 per-residue npz 是 {PERRES.stat().st_size / 1e6:.0f} MB; "
         f"未压缩, 因压缩会 OOM)")
    emit(f"内容: lig_topk/rec_topk ({n_pr}x{K}x640 **fp16**), 掩码, "
         f"pooled 嵌入 (fp16), 以及 V2.2/pooled 的行映射")
    emit(f"- 总耗时 {time.time() - t0:.0f}s")
    (LOG / "build_cache.log").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {LOG / 'build_cache.log'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
