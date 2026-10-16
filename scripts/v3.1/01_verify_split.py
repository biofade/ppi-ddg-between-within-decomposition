# -*- coding: utf-8 -*-
"""scripts/v3.1/01_verify_split.py -- V3.1 前置: 复核"family-level split"是否成立.

为什么必须先做这一步
--------------------
硬约束 12 要求**评估口径统一用 family-level split**。所有 V3.1 结论都建立在这个
划分上，因此必须**独立复核**它，而不是假设它成立:

1. 训练/验证/测试之间是否存在**序列同一性泄漏**（>30%）？
   注意本项目是**蛋白-蛋白复合物**，有**两侧**序列（配体 + 受体），
   必须**两侧都查**，否则一侧泄漏就会虚高指标。
2. 同一 `pdb_id` 是否跨划分？（结构级泄漏）
3. 数据集自带的 split 是否就是 <=30% 序列同一性的？（Proteinea 声称如此）

方法
----
用 3-mer Jaccard 近似序列同一性（与 V2.0 家族隔离同口径，k=3, 阈值 0.30）:
* 对 test 的每条序列，找 train 中 Jaccard 最高的序列；
* 若存在 >0.30 的配对 -> 该行存在潜在泄漏。

> 这是**近似**（3-mer Jaccard 不是真实比对同一性）。因此报告里给的是
> **上界式检查**: 若连近似条件下都没有 >30% 的配对, 则真实同一性泄漏更不可能。
> 反之若发现 >30% 配对, 需要进一步用真实比对确认 —— 这一点会在报告中说明。

输出: `logs/v3.1/split_verification.log`
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
LOG_DIR = ROOT / "logs" / "v3.1"
LOG_DIR.mkdir(parents=True, exist_ok=True)
OUT = LOG_DIR / "split_verification.log"

V22 = FEAT / "v2.2_proteinea_features.csv"
RAW = ROOT / "data" / "ppb_affinity_proteinea" / "filtered.csv"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def kmers(s: str, k: int = 3) -> set:
    return {s[i:i + k] for i in range(max(0, len(s) - k + 1))}


def main() -> int:
    emit("# V3.1 前置 — family-level split 复核")
    emit("")
    if not RAW.exists():
        emit(f"!! 缺少 {RAW}")
        return 2
    raw = pd.read_csv(RAW, low_memory=False)
    emit(f"- `filtered.csv`: {raw.shape[0]} 行 x {raw.shape[1]} 列")
    emit(f"- 列名含: {[c for c in raw.columns if 'eq' in c or 'split' in c.lower()][:6]}")
    LIG, REC, SPD = "Ligand Sequences", "Receptor Sequences", "split"
    for c in (LIG, REC, SPD):
        if c not in raw.columns:
            emit(f"!! 缺少列 `{c}`")
            return 3

    emit("")
    emit("## 1. 数据自带 split 的规模")
    emit("")
    emit(f"- split 分布: {raw[SPD].value_counts().to_dict()}")

    tr = raw[raw[SPD] == "train"]
    te = raw[raw[SPD] == "test"]
    va = raw[raw[SPD] == "val"]
    emit(f"- train {len(tr)} / val {len(va)} / test {len(te)}")

    # ---------- 2. 结构级泄漏 ----------
    emit("")
    emit("## 2. 结构级（pdb_id）泄漏检查")
    emit("")
    pdb_tr = set(tr["PDB"].astype(str).str.upper())
    pdb_te = set(te["PDB"].astype(str).str.upper())
    pdb_va = set(va["PDB"].astype(str).str.upper())
    emit(f"- train 唯一 pdb: {len(pdb_tr)} | val {len(pdb_va)} | test {len(pdb_te)}")
    emit(f"- **train ∩ test = {len(pdb_tr & pdb_te)}**")
    emit(f"- **train ∩ val  = {len(pdb_tr & pdb_va)}**")
    emit(f"- **val   ∩ test = {len(pdb_va & pdb_te)}**")
    struct_leak = len(pdb_tr & pdb_te) + len(pdb_tr & pdb_va) + len(pdb_va & pdb_te)
    emit(f"-> 结构级跨划分交集合计: **{struct_leak}** "
         + ("（无结构级泄漏 ✅）" if struct_leak == 0 else "（**存在泄漏** ❌）"))

    # ---------- 3. 序列级泄漏 (两侧都查) ----------
    emit("")
    emit("## 3. 序列同一性泄漏检查（3-mer Jaccard 近似，阈值 0.30）")
    emit("")
    emit("> 对 **test 与 val** 的每条序列，找 **train** 中 Jaccard 最高的序列。")
    emit("> 配体侧与受体侧**分别检查**。")
    emit("")

    rep = {}

    def max_jaccard_to_train(query_seqs, train_seqs, k=3):
        """返回每条 query 序列与 train 的最大 Jaccard。O(n*m)，用 k-mer 集合。"""
        tr_sets = [kmers(s, k) for s in train_seqs]
        out = np.zeros(len(query_seqs))
        for i, q in enumerate(query_seqs):
            qs = kmers(q, k)
            if not qs:
                out[i] = np.nan
                continue
            best = 0.0
            for ts in tr_sets:
                if not ts:
                    continue
                inter = len(qs & ts)
                if inter == 0:
                    continue
                j = inter / len(qs | ts)
                if j > best:
                    best = j
                    if best == 1.0:
                        break
            out[i] = best
        return out

    for tag, part in (("test", te), ("val", va)):
        if not len(part):
            continue
        # 去重后再比对，避免重复计算
        for side, col in (("配体", LIG), ("受体", REC)):
            q = sorted(set(part[col].dropna().astype(str)))
            t = sorted(set(tr[col].dropna().astype(str)))
            emit(f"### {tag} / {side}侧: {len(q)} 条唯一 query vs train {len(t)} 条唯一")
            # 控制规模: query 超过 800 条时抽样 (记录抽样)
            if len(q) > 800:
                rs = np.random.RandomState(42)
                qi = rs.choice(len(q), 800, replace=False)
                qs_sub = [q[i] for i in qi]
                emit(f"  （query 过多，随机抽 800 条，seed=42）")
            else:
                qs_sub = q
            mj = max_jaccard_to_train(qs_sub, t)
            mj = mj[np.isfinite(mj)]
            if not len(mj):
                continue
            n30 = int((mj > 0.30).sum())
            n50 = int((mj > 0.50).sum())
            n90 = int((mj > 0.90).sum())
            emit(f"  - 最大 Jaccard: 中位 {np.median(mj):.3f}, "
                 f"均值 {mj.mean():.3f}, 最大 {mj.max():.3f}")
            emit(f"  - **>0.30 的条数: {n30} / {len(mj)} ({n30/len(mj)*100:.2f}%)**")
            emit(f"  - >0.50: {n50} | >0.90: {n90}")
            rep[f"{tag}_{side}_gt30"] = n30
            rep[f"{tag}_{side}_n"] = len(mj)
            rep[f"{tag}_{side}_max"] = float(mj.max())
            emit("")

    # ---------- 4. 结论 ----------
    emit("## 4. 结论")
    emit("")
    tot_gt30 = sum(v for k, v in rep.items() if k.endswith("_gt30"))
    emit(f"- 结构级泄漏: **{struct_leak}** 个跨划分 pdb")
    emit(f"- 序列级 >0.30 的 query 总数（四项相加, 去重前）: **{tot_gt30}**")
    emit("")
    if struct_leak == 0 and tot_gt30 == 0:
        emit("> ✅ **family-level split 成立**：")
        emit("> 训练/验证/测试在**结构层面完全不相交**，")
        emit("> 且在两侧序列的 3-mer Jaccard 近似下**没有任何 >0.30 的配对**。")
        emit("> -> V3.1 所有评估可以直接使用该固定划分，无需重新划分。")
    elif struct_leak == 0:
        emit("> ⚠️ 结构层面无泄漏，但**近似序列检查发现少量 >0.30 配对**。")
        emit("> 需要注意: 3-mer Jaccard 会**高估**真实同一性（共享局部基序的序列），")
        emit("> 报告中标注该不确定性; 如需严格结论需用真实比对工具复核。")
    else:
        emit("> ❌ **发现结构级泄漏**，必须先修划分再谈评估。")
    emit("")
    emit("### 方法学限制（诚实标注）")
    emit("")
    emit("1. **3-mer Jaccard 不是真实比对同一性**。它对相似但不同的序列会给出偏低值，")
    emit("   对共享局部基序的序列会给出偏高值。因此本检查是**筛查**，不是定论。")
    emit("2. query 抽样的那几项（若发生）只覆盖部分序列。")
    emit("3. 未检查**交叉侧**组合泄漏（如 test 的配体 vs train 的受体）。")
    emit("   本项目预测的是复合物，两侧都来自同一 PDB，结构级检查已覆盖主要风险。")
    emit("")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
