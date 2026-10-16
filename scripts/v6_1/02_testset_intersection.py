# -*- coding: utf-8 -*-
"""scripts/v6_1/02_testset_intersection.py -- V6.1 PHASE 1.2: 测试集交集分析.

问：CATH held-out test 与我们的 V5 test **是不是同一批复合物**？
若不是，则 PHASE 4 只能做**分布对比**，不能做同集直接对比。

步骤
----
A. **PDB 级别（宽松）**：CATH test 唯一 PDB vs V5 test 唯一 PDB（小写化、剥离链）
B. **复合物级别（严格）**：用 (PDB, ligand_chain, receptor_chain) 或 interface_key

输出
----
* `experiments/v6_1/phase1_ci/intersection_report.md`
* `experiments/v6_1/phase1_ci/intersection_pdb.csv`
"""

from __future__ import annotations

import json
import re
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

PH1 = C.ROOT / "experiments" / "v6_1" / "phase1_ci"
PH1.mkdir(parents=True, exist_ok=True)
OUT_MD = PH1 / "intersection_report.md"
OUT_CSV = PH1 / "intersection_pdb.csv"
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def norm_pdb(x) -> str:
    """标准化 PDB ID：转小写、剥离链 ID 与非字母数字。"""
    s = str(x).strip().lower()
    s = re.split(r"[_|,;:\s]", s)[0]
    return re.sub(r"[^a-z0-9]", "", s)


def main() -> int:
    t0 = time.time()
    emit("# V6.1 PHASE 1.2 — 测试集交集分析")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------- CATH test 的 PDB ----------
    cath_pdb: set[str] = set()
    cath_cx: set[str] = set()
    std = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
    if std.exists():
        d = pd.read_csv(std, low_memory=False)
        ca = d[d["dataset_name"] == "CATH"]
        cath_cx = set(ca["complex_id"].astype(str))
        cath_pdb = {norm_pdb(x) for x in cath_cx}
    emit(f"- **CATH held-out test**（`standardized_preds.csv`, N≥5 过滤后）："
         f"**{len(cath_cx)} complex / {len(cath_pdb)} 唯一 PDB**")
    emit("")

    # ---------- V5 test 的 PDB ----------
    v5_pdb: set[str] = set()
    v5_cx: set[str] = set()
    v5_ik: set[str] = set()
    z = C.V5 / "p4_conformal" / "p4_predictions.npz"
    if z.exists():
        zz = np.load(z, allow_pickle=True)
        v5_cx = set(str(x) for x in zz["test_group"])
        v5_pdb = {norm_pdb(x) for x in v5_cx}
    emit(f"- **V5 test**（`p4_predictions.npz` test_group）："
         f"**{len(v5_cx)} complex / {len(v5_pdb)} 唯一 PDB**")
    emit("")

    # V5 的 interface_key（更严格的复合物定义）
    v22 = C.FEAT / "v2.2_proteinea_features.csv"
    if v22.exists() and v5_cx:
        keep = {"interface_key", "pdb_id"}
        try:
            v = pd.read_csv(v22, usecols=lambda c: c in keep, low_memory=False)
            if "pdb_id" in v.columns and "interface_key" in v.columns:
                v["pdb_id"] = v["pdb_id"].astype(str)
                m = v[v["pdb_id"].isin(v5_cx)]
                v5_ik = set(m["interface_key"].astype(str))
        except Exception as exc:  # noqa: BLE001
            emit(f"  （interface_key 读取失败：{type(exc).__name__}）")
    emit(f"- V5 test 涉及的 `interface_key`（严格复合物定义）：**{len(v5_ik)}**")
    emit("")

    # ---------- A. PDB 级别交集 ----------
    emit("## A. PDB 级别交集（宽松）")
    emit("")
    inter = sorted(cath_pdb & v5_pdb)
    emit(f"- CATH 唯一 PDB: **{len(cath_pdb)}**")
    emit(f"- V5   唯一 PDB: **{len(v5_pdb)}**")
    emit(f"- **交集: {len(inter)}**")
    emit(f"- Jaccard = {len(inter) / max(1, len(cath_pdb | v5_pdb)):.4f}")
    emit("")
    if inter:
        emit("交集 PDB（小写，已剥离链）：")
        emit("")
        emit("| # | PDB |")
        emit("|---:|---|")
        for i, x in enumerate(inter[:60], 1):
            emit(f"| {i} | `{x}` |")
        if len(inter) > 60:
            emit(f"| … | （共 {len(inter)} 个，全部见 CSV） |")
        emit("")
    else:
        emit("> **交集为 0** → CATH held-out test 与 V5 test **无任何共同 PDB**。")
        emit("")

    # ---------- B. 复合物级别交集 ----------
    emit("## B. 复合物级别交集（严格）")
    emit("")
    cath_cx_norm = {re.sub(r"[^a-z0-9]", "", str(x).lower()) for x in cath_cx}
    v5_cx_norm = {re.sub(r"[^a-z0-9]", "", str(x).lower()) for x in v5_cx}
    inter_cx = sorted(cath_cx_norm & v5_cx_norm)
    emit(f"- CATH complex_id 数: **{len(cath_cx_norm)}**")
    emit(f"- V5   complex_id 数: **{len(v5_cx_norm)}**")
    emit(f"- **交集: {len(inter_cx)}**")
    emit("")
    if v5_ik:
        ik_norm = {re.sub(r"[^a-z0-9]", "", x.lower()) for x in v5_ik}
        emit(f"- V5 `interface_key`（含链信息）与 CATH complex_id 的交集: "
             f"**{len(ik_norm & cath_cx_norm)}**")
        emit("")

    # ---------- 结论 ----------
    emit("## C. 结论与对后续 Phase 的影响")
    emit("")
    if len(inter) == 0:
        emit("> ### ❌ **两个测试集无可比交集**")
        emit("")
        emit("CATH held-out test（SKEMPI v2.0 的 CATH 同源超家族 hold-out）")
        emit("与我们的 V5 test（PPB-Affinity 三重隔离划分）**完全不重叠**。")
        emit("")
        emit("**这意味着**：")
        emit("")
        emit("1. **禁止**把 CATH 的数值与 V5 的数值做「同一测试集上的胜负比较」。")
        emit("2. PHASE 4 只能做**分布对比**（per-complex Spearman 的分布形状），")
        emit("   且必须标注「跨数据集比较，非同集比较」。")
        emit("3. 任何「我们的模型 vs Flex-ddG」的结论**只能**表述为")
        emit("   「在**各自最严格协议下**的组内排序能力量级对比」，")
        emit("   符合 V6.0 报告已采用的 `NOT_DIRECTLY_COMPARABLE` 口径。")
    elif len(inter) >= 10:
        emit(f"> ### ✅ 交集充足（{len(inter)} 个 PDB）→ 可以做**同集直接对比**")
    else:
        emit(f"> ### ⚠️ 交集偏少（{len(inter)} 个 PDB）→ 只做探索性对比，结论需谨慎")
    emit("")

    pd.DataFrame({"pdb_in_intersection": inter}).to_csv(OUT_CSV, index=False,
                                                       encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    print(f"写出 {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
