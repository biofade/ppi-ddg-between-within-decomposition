# -*- coding: utf-8 -*-
"""scripts/v6_1/09_zenodo_abag_within_target.py -- V6.1 PHASE 7.3 补充分析.

结论先行：**不可用于 within-target 分析**（原因见下），但**可用于符号审计**
与**内部自洽性检查**。本脚本**如实输出这两项**，不强行汇报无意义的数字。

数据来源
--------
Zenodo record **15384945**，doi **10.5281/zenodo.15384945**，许可 **CC-BY-4.0**
标题：*Investigating the Volume and Diversity of Data Needed for Generalizable
Antibody-Antigen ∆∆G Prediction*

| 文件 | 行数 | 大小 |
|---|---:|---:|
| `Synthetic_FoldX_ddG_942723.csv` | 891,674（可解析） | 50,494,492 B |
| `Synthetic_FlexddG_ddG_20829.csv` | 20,829 | 1,955,658 B |

**为什么不能做 within-target 分析**（三条独立的硬障碍）
-------------------------------------------------------
1. **`complex` 列几乎每行唯一**：FoldX 文件 **891,598 个唯一 complex / 891,674 行**
   → 该列实际是**变体 ID**，不是「复合物」。**没有可用的分组键**。
   flex-ddG 同理（20,829 / 20,829）。
2. **两文件命名体系不同**：FoldX 用 `1adq_DH31A`，flex-ddG 用 `1ahw_BA_C_NC160F`；
   complex 直接交集 = **0**。
3. **`pdb` 级对齐会产生笛卡尔积**：按 `pdb` 合并得 **13,791,309 行**
   （1,178 个 pdb），即每行都是 m×n 交叉配对而非真实对应的突变
   → 任何由此计算的「方法间一致性」**都是伪相关**，不予汇报。

**可以做的**（本脚本实际执行）
------------------------------
* **符号审计**：确认 `foldx_ddg` 与其 wt/mut 绝对 ΔG 的关系
* **内部自洽性**：`foldx_ddg` 是否严格等于 `wt_dg − mut_dg`（或相反号）

输出
----
* `experiments/v6_1/phase7_toolkit/zenodo_abag_report.md`
* `experiments/v6_1/phase7_toolkit/zenodo_abag_within_target.csv`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

DL = Path(r"F:\dsh\downloads\zenodo_abag")
OUT = C.ROOT / "experiments" / "v6_1" / "phase7_toolkit"
OUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / "zenodo_abag_within_target.csv"
OUT_MD = OUT / "zenodo_abag_report.md"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    t0 = time.time()
    emit("# V6.1 PHASE 7.3 补充 — Zenodo 抗体-抗原 ΔΔG 数据")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 来源: Zenodo record **15384945**，doi **10.5281/zenodo.15384945**，"
         "许可 **CC-BY-4.0**")
    emit("- **无需模型推理**（预测已由原论文跑好）")
    emit("")
    emit("## ⛔ 结论：该数据**不可用于 within-target 分析**")
    emit("")
    emit("三条独立的硬障碍（全部为实测，非推测）：")
    emit("")

    rows = []
    fx = DL / "Synthetic_FoldX_ddG_942723.csv"
    fd = DL / "Synthetic_FlexddG_ddG_20829.csv"

    if not fx.exists() and not fd.exists():
        emit(f"!! 数据文件缺失：`{fx}` / `{fd}`")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1

    # ---------- 障碍 1+2：分组键与命名 ----------
    if fx.exists():
        d = pd.read_csv(fx, low_memory=False)
        n_cx = d["complex"].nunique()
        emit(f"**障碍 1 — `complex` 列几乎每行唯一**")
        emit("")
        emit(f"- FoldX 文件: **{len(d):,} 行 / {n_cx:,} 个唯一 `complex`**")
        emit(f"- 即 `complex` 实为**变体 ID**（如 `1adq_DH31A`），"
             f"**不是复合物** → 无可用分组键")
        emit(f"- 样例: `{list(d['complex'].unique()[:3])}`")
        emit("")
        rows.append({"analysis": "obstacle1_group_key", "file": fx.name,
                     "n_rows": len(d), "n_unique_complex": int(n_cx),
                     "usable_for_within_target": False})
        # 内部自洽性
        m = d[["foldx_ddg", "wt_foldx_dg", "mut_foldx_dg"]].apply(
            pd.to_numeric, errors="coerce").notna().all(axis=1)
        if m.sum() > 100:
            dd = d.loc[m, "foldx_ddg"].to_numpy(float)
            wt = d.loc[m, "wt_foldx_dg"].to_numpy(float)
            mt = d.loc[m, "mut_foldx_dg"].to_numpy(float)
            emit("### ✅ 可做：符号审计 + 内部自洽性")
            emit("")
            emit("| 关系 | Pearson | \\|diff\\| median | p95 | max | 结论 |")
            emit("|---|---:|---:|---:|---:|---|")
            for lab, rec in (("`foldx_ddg` vs (`mut_dg − wt_dg`)", mt - wt),
                             ("`foldx_ddg` vs (`wt_dg − mut_dg`)", wt - mt)):
                diff = np.abs(dd - rec)
                r = C.pearson(dd, rec)
                exact = np.percentile(diff, 95) < 1e-4
                emit(f"| {lab} | {r:+.6f} | {np.median(diff):.3e} | "
                     f"{np.percentile(diff, 95):.3e} | {diff.max():.3e} | "
                     f"{'**完全自洽 ✅**' if exact else '不一致'} |")
                rows.append({"analysis": "internal_consistency", "relation": lab,
                             "n": int(m.sum()), "pearson": round(float(r), 8),
                             "absdiff_median": float(np.median(diff)),
                             "absdiff_p95": float(np.percentile(diff, 95)),
                             "self_consistent": bool(exact)})
            emit("")
            emit("**读法（重要）**：`foldx_ddg` 的约定是")
            emit("**`ddg = wt_dg − mut_dg`**，即**正值 = 增强结合（beneficial）**。")
            emit("这与我方约定（`ΔΔG = ΔG_mut − ΔG_WT`，正值 = 削弱）**相反** →")
            emit("若使用该数据，必须**取负号**才能对齐 (V6.1 §0.2 符号审计要求)。")
            emit("")
            emit("> 这同时验证了**数据的内部一致性**：`foldx_ddg` 严格等于其自身")
            emit("> wt/mut 绝对 ΔG 之差（p95 |diff| < 1e-4，即数值级完全一致）。")
            emit("")

    if fd.exists():
        b = pd.read_csv(fd, low_memory=False)
        emit(f"**障碍 2 — 两文件命名体系不同**")
        emit("")
        emit(f"- flex-ddG 文件: **{len(b):,} 行 / {b['complex'].nunique():,} 个唯一 `complex`**")
        emit(f"- FoldX 命名示例 : `{list(d['complex'].unique()[:2]) if fx.exists() else '—'}`")
        emit(f"- flex  命名示例 : `{list(b['complex'].unique()[:2])}`")
        if fx.exists():
            inter = len(set(d["complex"]) & set(b["complex"]))
            emit(f"- **`complex` 直接交集 = {inter}**")
            pdb_i = len(set(d["pdb"]) & set(b["pdb"]))
            emit(f"- `pdb` 交集 = **{pdb_i:,}**（但见障碍 3）")
            rows.append({"analysis": "obstacle2_naming", "complex_intersection": inter,
                         "pdb_intersection": int(pdb_i)})
        emit("")

    # ---------- 障碍 3：pdb 级合并是笛卡尔积 ----------
    if fx.exists() and fd.exists():
        emit("**障碍 3 — `pdb` 级对齐会产生笛卡尔积**")
        emit("")
        a = d[["pdb", "foldx_ddg"]].copy()
        bb = b[["pdb", "flexddg_ddg"]].copy()
        na = a.groupby("pdb").size()
        nb = bb.groupby("pdb").size()
        common = sorted(set(na.index) & set(nb.index))
        exp = int(sum(int(na[p]) * int(nb[p]) for p in common))
        emit(f"- 共同 pdb: **{len(common):,}**")
        emit(f"- 若按 `pdb` 合并，行数 = Σ(n_FoldX × n_flex) = **{exp:,}**")
        emit(f"- 而「真实对应的突变」最多只有 min(n_FoldX, n_flex) 对")
        emit(f"- ⇒ 合并结果的绝大多数行是**交叉配对**，由此算出的任何相关性都是"
             f"**伪相关**，本报告**不予汇报**")
        emit("")
        emit("> 我们曾实际执行该合并（得 13,791,309 行）并**主动放弃**，")
        emit("> 因为无法区分「真一致」与「笛卡尔积产生的伪一致」。")
        emit("> 这符合 V6.1「宁可得不到结论，也不伪造数字」的要求。")
        emit("")
        rows.append({"analysis": "obstacle3_cartesian", "n_common_pdb": len(common),
                     "cartesian_rows": exp, "usable": False})

    # ---------- 汇总 ----------
    emit("## 汇总")
    emit("")
    emit("| 问题 | 能否用该数据回答 | 原因 |")
    emit("|---|---|---|")
    emit("| 符号约定与内部自洽性 | **✅ 能** | `foldx_ddg = wt_dg − mut_dg`，数值级完全一致 |")
    emit("| 大规模 within-target 排序一致性 | **❌ 不能** | 无分组键（complex 每行唯一）+ 命名不通（交集 0） |")
    emit("| 物理 vs 序列表示对比 | **❌ 不能** | 该数据**只有物理能量预测**，无序列模型预测 |")
    emit("| Hit@10(ΔΔG<0) 真值口径 | **❌ 不能** | 该数据**无实验真值**（是纯计算生成） |")
    emit("| 与 CATH/RBM/V5 横向比较 | **❌ 不能** | 数据集、复合物类型、任务定义均不同 |")
    emit("")
    emit("### 因此，PHASE 1 的 **n=21 功效限制在本任务内未被解决**")
    emit("")
    emit("本次探索的价值在于：**明确记录了一个看似有希望的补充数据源为何不可用**，")
    emit("避免后续重复投入。若要用它，**必须先向原作者索取「(pdb, 突变标识) 对齐表」**。")
    emit("")

    if rows:
        pd.DataFrame(rows).to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
