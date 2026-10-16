# -*- coding: utf-8 -*-
"""scripts/v2.5/00_join_check.py -- V2.5 前置: per-residue 缓存与 V2.2 固定划分的 join 校验.

**为什么必须先做这一步**
------------------------
任务 1 的判决实验要求「用 Proteinea 固定划分, 对比 pooled vs per-residue」。
若 per-residue 缓存与 V2.2 特征表**不能逐行对齐**, 那么:
* 训练集不同 -> 与 V2.3.2 的 0.3890 **不可比**;
* 测试集不同 -> 判决实验本身失效。

实测已知: per-residue 缓存的 splits 是 test 757 / val 965 / **train 5163**,
合计 6885; 而 V2.2 是 test 757 / val 965 / **train 6485** (可用 6482), 合计 8207。
**train 少了 1319 行** -> 必须查清是哪些行, 以及能否建立精确映射。

本脚本输出 `logs/v2.5/join_check.log` 与结论。
"""

from __future__ import annotations

import sys
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
OUT = LOG / "join_check.log"

PERRES = FEAT / "v2_esm2_perres_v2.npz"
V22 = FEAT / "v2.2_proteinea_features.csv"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    emit("=== V2.5 前置: per-residue 与 V2.2 固定划分的 join 校验 ===")
    emit("")

    z = np.load(PERRES, allow_pickle=True, mmap_mode="r")
    pr_pdb = np.asarray(z["pdb_ids"]).astype(str)
    pr_ik = np.asarray(z["interface_key"]).astype(str)
    pr_sp = np.asarray(z["splits"]).astype(str)
    n_pr = len(pr_pdb)
    emit(f"per-residue 缓存: {n_pr} 行")
    emit(f"  splits: {dict(zip(*np.unique(pr_sp, return_counts=True)))}")

    d = pd.read_csv(V22, low_memory=False)
    n_v22 = len(d)
    emit(f"V2.2 特征表: {n_v22} 行")
    emit(f"  splits: {d['split'].value_counts().to_dict()}")

    # ---------- 1. interface_key 能否直接 join ----------
    emit("")
    emit("--- 1) interface_key 直接 join ---")
    v22_ik = set(d["interface_key"].astype(str))
    pr_ik_set = set(pr_ik)
    inter = v22_ik & pr_ik_set
    emit(f"  V2.2 唯一 interface_key : {len(v22_ik)}")
    emit(f"  perres 唯一 interface_key: {len(pr_ik_set)}")
    emit(f"  交集                    : {len(inter)}")
    if len(inter) == 0:
        emit("  -> interface_key **格式不同** (与 V2.3/V2.4 同一个坑: perres 是 "
             "`pdb|配体序列|受体序列`)")
    emit(f"  perres 样例: {pr_ik[0][:70]}...")
    emit(f"  V2.2  样例: {str(d['interface_key'].iloc[0])}")

    # ---------- 2. pdb_id join 能否覆盖 ----------
    emit("")
    emit("--- 2) pdb_id 覆盖 ---")
    v22_pdb = set(d["pdb_id"].astype(str))
    pr_pdb_set = set(pr_pdb)
    emit(f"  V2.2 唯一 pdb: {len(v22_pdb)} | perres 唯一 pdb: {len(pr_pdb_set)} | "
         f"交集: {len(v22_pdb & pr_pdb_set)}")
    missing_pdb = v22_pdb - pr_pdb_set
    emit(f"  V2.2 有但 perres 没有的 pdb: {len(missing_pdb)} 个 "
         f"-> 覆盖 {int(d['pdb_id'].astype(str).isin(missing_pdb).sum())} 行")

    # ---------- 3. 逐行精确映射 (pdb_id 序列 + 组内次序) ----------
    emit("")
    emit("--- 3) 逐行精确映射尝试 ---")
    emit(f"  行数: V2.2 {n_v22} vs perres {n_pr} -> "
         f"{'相等' if n_v22 == n_pr else '**不相等**'}")
    if n_v22 == n_pr:
        same = bool((d["pdb_id"].astype(str).to_numpy() == pr_pdb).all())
        emit(f"  pdb_id 逐位置一致: {same}")
        same_sp = bool((d["split"].astype(str).to_numpy() == pr_sp).all())
        emit(f"  split  逐位置一致: {same_sp}")
    # 用 (interface_key 不行) -> 试 (pdb_id, 组内次序)
    emit("")
    emit("  --- 用 (pdb_id, 组内出现次序) 做键 ---")
    a = d.copy()
    a["_seq"] = a.groupby("pdb_id").cumcount()
    b = pd.DataFrame({"pdb_id": pr_pdb})
    b["_seq"] = b.groupby("pdb_id").cumcount()
    pos = a.set_index(["pdb_id", "_seq"]).index
    loc = pos.get_indexer(pd.MultiIndex.from_arrays([b["pdb_id"], b["_seq"]]))
    hit = int((loc >= 0).sum())
    emit(f"  perres 的 {n_pr} 行中, 能在 V2.2 里定位的: **{hit}** "
         f"({hit / n_pr * 100:.1f}%)")
    if hit < n_pr:
        miss = np.where(loc < 0)[0]
        emit(f"  未命中 {len(miss)} 行; 前 10 个:")
        for i in miss[:10]:
            emit(f"     pdb={b['pdb_id'].iloc[i]} seq={b['_seq'].iloc[i]} "
                 f"split={pr_sp[i]}")
        msp = pd.Series(pr_sp[miss]).value_counts().to_dict()
        emit(f"  未命中行的 split 分布: {msp}")
    if hit:
        matched = loc >= 0
        v22_sp = d["split"].astype(str).to_numpy()[loc[matched]]
        agree = float((v22_sp == pr_sp[matched]).mean())
        emit(f"  命中行里 split 一致率: **{agree * 100:.2f}%**")
        if agree < 1.0:
            bad = np.where(matched)[0][v22_sp != pr_sp[matched]]
            emit(f"  !! split 不一致 {len(bad)} 行; 前 5:")
            for i in bad[:5]:
                emit(f"     pdb={pr_pdb[i]} perres_split={pr_sp[i]} "
                     f"v22_split={d['split'].astype(str).iloc[loc[i]]}")

    # ---------- 4. 结论 ----------
    emit("")
    emit("=== 结论 ===")
    if hit == n_pr:
        emit(f"* **per-residue 缓存全部 {n_pr} 行都能映射到 V2.2 固定划分** -> "
             f"可直接用于判决实验。")
        n_tr = int((pr_sp == "train").sum())
        emit(f"* 但 perres 的 train 只有 {n_tr} 行 (V2.2 可用 6482 行) ->")
        emit(f"  训练集**小于** V2.3.2 基线。因此判决实验必须:")
        emit(f"  (a) pooled 基线**在同一批行上重跑** (而不是直接引用 0.3890), 或")
        emit(f"  (b) 明确标注训练量差异。")
        emit(f"  -> 采用 (a): 用 perres 覆盖的训练行重训 pooled 基线, 保证同训练集同测试集。")
    else:
        emit(f"* 有 {n_pr - hit} 行无法映射 -> 判决实验需在这些行上受限, 必须标注。")

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
