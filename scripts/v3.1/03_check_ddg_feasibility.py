# -*- coding: utf-8 -*-
"""scripts/v3.1/03_check_ddg_feasibility.py -- V3.1 前置: ΔΔG 双头是否可行?

为什么必须先查
--------------
文献里最可迁移的两个方向都**依赖"同一复合物内的多突变体"结构**:

* **Twin Peaks 双头**：ΔΔG 头需要 (WT, mutant) 配对或同复合物内的相对变化 [Twin Peaks p7]
* **Pythia-PPI 自蒸馏**：需要界面残基 + 饱和突变 [Pythia-PPI p10]

如果我的 `ppb_affinity` 里**同一复合物只有一条记录**，这两个方向**根本无法构造目标**，
必须在最小验证前就知道，而不是写完代码才发现。

本脚本统计（只用结构/标签，不使用测试标签做选择）：
1. 每个 `interface_key` / `pdb_id` 有多少条记录（突变体数）；
2. 同一复合物内 ΔG 的**离散程度**（有变化才能定义 ΔΔG）；
3. 有多少组可以构造**成对的相对 ΔΔG**（组内两两差值）；
4. 突变计数列（`Mutations`）的可解析性。
"""

from __future__ import annotations

import re
import sys
from collections import Counter
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
OUT = LOG_DIR / "ddg_feasibility.log"
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
    emit("# V3.1 前置 — ΔΔG 双头可行性检查")
    emit("")
    d = pd.read_csv(V22, low_memory=False)
    d = d[d["row_structure_usable"] == 1].reset_index(drop=True)
    emit(f"- 可用行: **{len(d)}**; split = {d['split'].value_counts().to_dict()}")
    emit("")

    # ---------- 1. 分组规模 ----------
    emit("## 1. 分组规模")
    emit("")
    emit("| group by | 组数 | 每组行数 中位 | 均值 | 最大 | 组数(>=2) | 覆盖行(>=2) |")
    emit("|---|---:|---:|---:|---:|---:|---:|")
    for gcol in ("interface_key", "pdb_id"):
        vc = d.groupby(gcol).size()
        emit(f"| `{gcol}` | {len(vc)} | {vc.median():.0f} | {vc.mean():.2f} | "
             f"{vc.max()} | {int((vc >= 2).sum())} | {int(vc[vc >= 2].sum())} |")
    emit("")

    # ---------- 2. 组内 ΔG 变化 ----------
    emit("## 2. 组内 ΔG 是否变化（决定 ΔΔG 能否定义）")
    emit("")
    emit("| group by | 组数(>=2) | 组内 ΔG 有变化 | 组内 ΔG 全相同 | 有变化占比 |")
    emit("|---|---:|---:|---:|---:|")
    for gcol in ("interface_key", "pdb_id"):
        sub = d.groupby(gcol).filter(lambda x: len(x) >= 2)
        if not len(sub):
            emit(f"| `{gcol}` | 0 | 0 | 0 | — |")
            continue
        g = sub.groupby(gcol)["real_dG"]
        n_var = int(g.nunique().gt(1).sum())
        n_tot = int(sub[gcol].nunique())
        emit(f"| `{gcol}` | {n_tot} | {n_var} | {n_tot - n_var} | "
             f"{n_var / max(1, n_tot) * 100:.1f}% |")
    emit("")

    # 训练集内（用于构造 ΔΔG 训练目标）
    tr = d[d["split"] == "train"]
    emit(f"### 训练集内（{len(tr)} 行）")
    emit("")
    for gcol in ("interface_key", "pdb_id"):
        sub = tr.groupby(gcol).filter(lambda x: len(x) >= 2)
        if not len(sub):
            emit(f"- `{gcol}`: **没有 >=2 行的组** -> 无法构造组内 ΔΔG")
            continue
        g = sub.groupby(gcol)["real_dG"]
        vc = sub.groupby(gcol).size()
        n_var = int(g.nunique().gt(1).sum())
        n_tot = int(sub[gcol].nunique())
        # 可构造成对差值数
        n_pairs = int((vc * (vc - 1) / 2).sum())
        emit(f"- `{gcol}`: 组数(>=2) **{n_tot}**, 覆盖 **{len(sub)}** 行; "
             f"组内 ΔG 有变化 **{n_var}** ({n_var / max(1, n_tot) * 100:.1f}%); "
             f"**可构造成对 ΔΔG 的对数 = {n_pairs}**")
    emit("")

    # ---------- 3. 突变列可解析性 ----------
    emit("## 3. `Mutations` 列可解析性")
    emit("")
    if "Mutations" in d.columns:
        mu = d["Mutations"].astype(str)
        nonnull = d["Mutations"].notna()
        emit(f"- 非空: {int(nonnull.sum())} / {len(d)} ({nonnull.mean() * 100:.1f}%)")
        pat = re.compile(r"^[A-Za-z]_?[A-Z]\d+[A-Za-z]$")
        ok = mu.str.match(r"^[A-Za-z]{1,3}_?[A-Z]\d+[A-Za-z]$", na=False)
        emit(f"- 匹配 `XXX?_?X###X` 单点突变模式: {int(ok.sum())} "
             f"({ok.mean() * 100:.1f}%)")
        emit(f"- 样例: {mu[nonnull].head(6).tolist()}")
        # 多点突变（含 ; 或 ,）
        multi = mu.str.contains(r"[;,]", na=False)
        emit(f"- 含多突变分隔符 (`;` 或 `,`): {int(multi.sum())}")
    else:
        emit("- **没有 `Mutations` 列**")
    emit("")

    # ---------- 4. 结论 ----------
    emit("## 4. 结论")
    emit("")
    tr_ik = tr.groupby("interface_key").filter(lambda x: len(x) >= 2)
    tr_pd = tr.groupby("pdb_id").filter(lambda x: len(x) >= 2)
    n_pairs_ik = int((tr_ik.groupby("interface_key").size().pipe(
        lambda s: s * (s - 1) / 2)).sum()) if len(tr_ik) else 0
    n_pairs_pd = int((tr_pd.groupby("pdb_id").size().pipe(
        lambda s: s * (s - 1) / 2)).sum()) if len(tr_pd) else 0
    emit(f"* 训练集可构造成对 ΔΔG 的对数：`interface_key` **{n_pairs_ik}**、"
         f"`pdb_id` **{n_pairs_pd}**")
    emit("")
    if n_pairs_ik >= 200:
        emit("> ✅ **ΔΔG 双头可行**：训练集内有足够的同复合物多突变体可构造相对目标。")
        emit("> （Twin Peaks 的双头设计 [Twin Peaks p7] 可以在本数据上落地。）")
    elif n_pairs_ik > 0:
        emit(f"> ⚠️ **勉强可行但样本很少**（{n_pairs_ik} 对）。")
        emit("> 双头训练目标会很稀疏 -> 预期收益有限，必须如实报告这一点。")
    else:
        emit("> ❌ **ΔΔG 双头不可行**：训练集内没有同复合物多突变体。")
        emit("> 这与 Twin Peaks [p7] / Pythia-PPI [p10] 所需的数据形态**不匹配**。")
    emit("")
    emit("### 与 V2.4 的联系")
    emit("")
    emit("V2.4 曾用 **200 个复合物 / 4987 行**做 within-target 基准（SKEMPI 子集）。")
    emit("那些行**只存在于 SKEMPI 子集**，而本脚本统计的是**整个 V2.2 表**。")
    emit("若上表显示 SKEMPI 之外的分组多为单行，则说明：")
    emit("* 我的 8207 行里**大部分复合物只有一条记录**（跨靶数据集的常态）；")
    emit("* \"同复合物多突变体\"的样本主要集中在 SKEMPI，且 **98% 落在 train**（V2.4 已实测）。")
    emit("")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
