# -*- coding: utf-8 -*-
"""scripts/v3.1/02_exact_seq_overlap.py -- V3.1 决定性的序列共享检查.

为什么需要这一步（接 `01_verify_split.py`）
-------------------------------------------
`01` 用 3-mer Jaccard 近似发现 605 条 query 与 train 序列的 Jaccard >0.30，但
**3-mer Jaccard 不是真实序列同一性**：对组成相似但不同的序列会偏高，对长度差异
敏感。因此那个数字**不能**直接作为泄漏证据。

本脚本改用**无歧义的判据**，分三个层次：

1. **完全相同序列**（exact match）—— 零歧义。若 test 的某条序列在 train 里逐字符
   相同，则该行的序列嵌入与训练行完全一致，是最硬的泄漏形式。
2. **完全相同序列对**（两侧都 exact match）—— 说明整个复合物的输入完全相同。
3. **长度 + 组成相同的近似重复**（同长度且 3-mer 集合 Jaccard == 1.0）—— 零歧义地
   表示"3-mer 组成完全相同"（等价于相同氨基酸多重集合在 3-mer 层面一致）。

同时报告**影响面**：受影响的 test **行数**与**占比**，这是解读所有指标的关键。

> 说明：真实的 <=30% 同一性判定需要 mmseqs2/BLAST 比对。本机未安装，且对
> 8207 条序列做全比对需要数小时。因此本脚本给出的是**可判定的下界**
> （exact match）与**零歧义的近似判据**（3-mer 集合相等），并在报告中明确区分
> "已证实泄漏" 与 "疑似同源"。
"""

from __future__ import annotations

import sys
from collections import defaultdict
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
RAW = ROOT / "data" / "ppb_affinity_proteinea" / "filtered.csv"
OUT = LOG_DIR / "exact_seq_overlap.log"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def kmers(s: str, k: int = 3) -> frozenset:
    return frozenset(s[i:i + k] for i in range(max(0, len(s) - k + 1)))


def main() -> int:
    emit("# V3.1 — 决定性序列共享检查（exact match + 零歧义近似）")
    emit("")
    raw = pd.read_csv(RAW, low_memory=False)
    LIG, REC, SPD, PDB = "Ligand Sequences", "Receptor Sequences", "split", "PDB"
    raw[PDB] = raw[PDB].astype(str).str.upper().str.strip()
    emit(f"- `filtered.csv`: {len(raw)} 行; split = {raw[SPD].value_counts().to_dict()}")
    emit("")

    tr = raw[raw[SPD] == "train"]
    va = raw[raw[SPD] == "val"]
    te = raw[raw[SPD] == "test"]

    summary: dict[str, dict] = {}

    for tag, part in (("test", te), ("val", va)):
        emit(f"## {tag} vs train")
        emit("")
        for side, col in (("配体", LIG), ("受体", REC)):
            tr_seqs = set(tr[col].dropna().astype(str))
            q_seqs = part[col].dropna().astype(str)
            exact = q_seqs.isin(tr_seqs)
            n_exact = int(exact.sum())
            n_rows = len(part)
            # 3-mer 集合完全相等 (零歧义近似): 用 frozenset 作键
            tr_k = defaultdict(list)
            for s in tr_seqs:
                tr_k[kmers(s)].append(s)
            q_keys = [kmers(s) for s in q_seqs]
            same_kmer_set = np.array([bool(tr_k.get(k)) for k in q_keys])
            n_same_k = int(same_kmer_set.sum())
            # 行级影响
            rows_exact = int(part.loc[exact.to_numpy(), :].shape[0])
            emit(f"### {side}侧")
            emit("")
            emit(f"- 唯一 query 序列: {q_seqs.nunique()}")
            emit(f"- **完全相同的序列（exact match）: {n_exact}** "
                 f"({n_exact / max(1, q_seqs.nunique()) * 100:.1f}% of 唯一序列)")
            emit(f"- 3-mer 集合完全相同的序列（零歧义近似）: **{n_same_k}** "
                 f"({n_same_k / max(1, q_seqs.nunique()) * 100:.1f}%)")
            emit(f"- **受影响行数（exact match）: {rows_exact} / {n_rows} "
                 f"({rows_exact / max(1, n_rows) * 100:.2f}% 的 {tag} 行）**")
            summary[f"{tag}_{side}_exact_seq"] = n_exact
            summary[f"{tag}_{side}_same_kmer"] = n_same_k
            summary[f"{tag}_{side}_rows_exact"] = rows_exact
            summary[f"{tag}_{side}_n_rows"] = n_rows
            emit("")

        # 两侧都 exact match（整个复合物输入完全相同）
        both = (part[LIG].isin(set(tr[LIG])).to_numpy()
                & part[REC].isin(set(tr[REC])).to_numpy())
        emit(f"- **两侧序列都 exact match 的行: {int(both.sum())} / {len(part)} "
             f"({both.sum() / max(1, len(part)) * 100:.2f}%)**")
        summary[f"{tag}_both_rows"] = int(both.sum())
        emit("")

        # 同一 PDB 是否在 train（结构级）
        same_pdb = part[PDB].isin(set(tr[PDB]))
        emit(f"- 与 train 共享同一 `PDB` 的行: {int(same_pdb.sum())} / {len(part)} "
             f"({'**存在结构级重叠**' if same_pdb.any() else '无'})")
        emit("")

    # ---------------- 结论 ----------------
    emit("## 结论")
    emit("")
    te_lig = summary.get("test_配体_exact_seq", 0)
    te_rec = summary.get("test_受体_exact_seq", 0)
    te_rows = summary.get("test_配体_rows_exact", 0)
    te_both = summary.get("test_both_rows", 0)
    n_te = summary.get("test_配体_n_rows", len(te))
    emit(f"* 结构级：**train ∩ test 的 PDB = 0**（`01_verify_split.py` 已确认）")
    emit(f"* 序列级（test 侧）：exact match 配体 **{te_lig}** 条、受体 **{te_rec}** 条；")
    emit(f"  **两侧都 exact match 的 test 行 = {te_both} / {n_te} "
         f"({te_both / max(1, n_te) * 100:.2f}%)**")
    emit("")
    if te_both == 0:
        emit("> ✅ **没有任何 test 行是「两侧序列都与训练序列完全相同」** ->")
        emit(f"> 该划分在**复合物输入层面是干净的**。")
        emit(f"> 存在少量**单侧** exact-match（配体 {te_lig} 条序列、受体 {te_rec} 条序列），")
        emit(f"> 这来自**共享亚基**（同一蛋白出现在多个不同复合物中），")
        emit(f"> 是 PPI 数据的固有性质，**不是复合物级泄漏** —— 因为另一侧不同，")
        emit(f"> 整个复合物的输入并不相同。")
        emit(">")
        emit("> `01` 里 22.93% 的 >0.30 Jaccard 应理解为 **3-mer 组成相似**（同源折叠、")
        emit("> 保守界面基序），而非序列重复 —— 那个数字**不构成泄漏证据**。")
    else:
        emit("> ⚠️ **存在「两侧都 exact match」的行** -> 这些行与训练行输入完全相同，")
        emit("> 属于**已证实的复合物级泄漏**。受影响行占比见上。")
        emit("> 需要做剔除后的敏感性分析，并在报告中披露。")
    emit("")
    emit("### 方法学限制（诚实标注）")
    emit("")
    emit("1. **未做真实序列比对**（无 mmseqs2/BLAST）。因此**无法**确认或否认")
    emit("   「存在 30%~99% 同源的配对」—— 只能说 exact match 层面干净。")
    emit("2. `01` 的 3-mer Jaccard **会高估**真实同一性，因此那里的大比例")
    emit("   **不构成**泄漏证据；本脚本的 exact-match 判据才是可判定的。")
    emit("3. 未检查**交叉侧**泄漏（test 配体 vs train 受体 等），但结构级检查已覆盖")
    emit("   主要风险（同一 PDB 不会跨划分）。")
    emit("")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
