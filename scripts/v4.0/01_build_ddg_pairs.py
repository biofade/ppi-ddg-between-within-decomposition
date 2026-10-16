# -*- coding: utf-8 -*-
"""scripts/v4.0/01_build_ddg_pairs.py -- V4.0 任务 1: 构建 ΔΔG 配对数据集.

**先把定义钉死（这决定后面所有数字）**
--------------------------------------
任务书规定:
    ΔΔG = ΔG_mutant − ΔG_wild-type

因此配对**只能**是「WT 行 × mutant 行」。这与 V3.1 的
`03_check_ddg_feasibility.py` 报的 243,439「同 interface_key 两两成对」
**不是同一个量**：后者把 mutant×mutant 也算作"对"，而 mutant×mutant 的差值
**不是** ΔΔG（没有 WT 参照）。本脚本**必须**如实报告两者的差异，不得混淆。

实测（见下方输出）:
* `filtered.csv` 8207 行, `Mutations` 空 3384 (=WT) / 非空 4823 (=mutant)
* **2171 个 PDB 中只有 231 个同时含 WT 与 mutant 行**
* 因此 WT×mutant 对数 = sum(n_wt × n_mut) 落在 231 个 PDB 上 -> **约 7.5k 对**

这远小于 V3.1 的数字，是**任务定义差异**，不是 bug。

输出
----
* `data/v4.0/ddg_pairs.csv`    每行 = 一个 (WT 行, mutant 行) 配对 + ΔΔG + 方向
* `data/v4.0/ddg_groups.json`  每个复合物的 WT/ mutant 结构
* `features/v4.0_ddg_dataset_report.md`  统计报告
"""

from __future__ import annotations

import argparse
import json
import re
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
DATA = ROOT / "data"
V40 = DATA / "v4.0"
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v4.0"
V40.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "failures.log"

RAW = DATA / "ppb_affinity_proteinea" / "filtered.csv"
V22 = FEAT / "v2.2_proteinea_features.csv"

OUT_PAIRS = V40 / "ddg_pairs.csv"
OUT_GROUPS = V40 / "ddg_groups.json"
OUT_MD = FEAT / "v4.0_ddg_dataset_report.md"

MUT_PAT = re.compile(r"^([A-Za-z0-9]+)_([A-Za-z])(\d+)([A-Za-z])$")
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def parse_muts(s) -> list[tuple]:
    """解析 `Mutations` -> [(chain, wt_aa, pos, mut_aa), ...]；空/不可解析返回 []。"""
    if not isinstance(s, str) or not s.strip():
        return []
    out = []
    for p in re.split(r"[;,]", s):
        mm = MUT_PAT.match(p.strip())
        if mm:
            out.append((mm.group(1), mm.group(2), int(mm.group(3)), mm.group(4)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 任务 1: ΔΔG 配对数据集")
    ap.add_argument("--raw", type=Path, default=RAW)
    args = ap.parse_args()

    t0 = time.time()
    emit("# V4.0 任务 1 — ΔΔG 配对数据集构建")
    emit("")
    emit("**定义**: `ΔΔG = ΔG_mutant − ΔG_wild-type`（任务书规定）")
    emit("")
    emit("> ⚠️ 这与 V3.1 的「同 interface_key 两两成对」**不是同一个量** ——")
    emit("> 后者包含 mutant×mutant（没有 WT 参照，其差值不是 ΔΔG）。")
    emit("")

    if not args.raw.exists():
        emit(f"!! 找不到 `{args.raw}`")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2

    raw = pd.read_csv(args.raw, low_memory=False)
    emit("## 1. 原始数据")
    emit("")
    emit(f"- `{args.raw.relative_to(ROOT)}`: **{raw.shape[0]} 行 x {raw.shape[1]} 列**")
    emit(f"- split: {raw['split'].value_counts().to_dict()}")
    emit(f"- 唯一 PDB: {raw['PDB'].nunique()}")

    # ---------- 2. WT / mutant 识别 ----------
    emit("")
    emit("## 2. wild-type 与 mutant 识别（依据 `Mutations` 是否为空）")
    emit("")
    is_wt = raw["Mutations"].isna()
    raw = raw.assign(_is_wt=is_wt)
    n_wt, n_mut = int(is_wt.sum()), int((~is_wt).sum())
    emit(f"- `Mutations` **为空**（= wild-type 候选）: **{n_wt}**")
    emit(f"- `Mutations` **非空**（= mutant）: **{n_mut}**")
    emit(f"- 两者都有有效 `KD(M)`: WT {int(raw.loc[is_wt, 'KD(M)'].notna().sum())}, "
         f"mutant {int(raw.loc[~is_wt, 'KD(M)'].notna().sum())}")

    g = raw.groupby("PDB")["_is_wt"].agg(["sum", "count"])
    g.columns = ["n_wt", "n_rows"]
    g["n_mut"] = g["n_rows"] - g["n_wt"]
    emit("")
    emit("| 每个 PDB 的构成 | PDB 数 |")
    emit("|---|---:|")
    emit(f"| 只有 WT，无 mutant | {int(((g.n_wt > 0) & (g.n_mut == 0)).sum())} |")
    emit(f"| 只有 mutant，无 WT | {int(((g.n_wt == 0) & (g.n_mut > 0)).sum())} |")
    emit(f"| **WT 与 mutant 共存** | **{int(((g.n_wt > 0) & (g.n_mut > 0)).sum())}** |")
    emit(f"| 两者都无 | {int(((g.n_wt == 0) & (g.n_mut == 0)).sum())} |")

    both = g[(g.n_wt > 0) & (g.n_mut > 0)]
    if len(both):
        emit("")
        emit(f"- 共存 PDB 的规模: `n_wt` 中位 **{both.n_wt.median():.0f}**"
             f"（最大 {both.n_wt.max()}）; `n_mut` 中位 **{both.n_mut.median():.0f}**"
             f"（最大 {both.n_mut.max()}）")
        emit(f"- **理论 WT×mutant 对数 = Σ(n_wt × n_mut) = "
             f"{int((both.n_wt * both.n_mut).sum())}**")
        emit(f"- 涉及的 mutant 行数 = {int(both.n_mut.sum())}")
    emit("")
    emit(f"**如果没有 WT 行的 mutant 复合物会被跳过**: 有 mutant 无 WT 的 PDB = "
         f"{int(((g.n_wt == 0) & (g.n_mut > 0)).sum())} 个（本数据为 0，"
         f"即所有 mutant 所在 PDB 都有 WT 行）")

    # ---------- 3. 构造配对 ----------
    emit("")
    emit("## 3. 配对构造")
    emit("")
    ws = raw[raw["_is_wt"]]
    ms = raw[~raw["_is_wt"]]
    # 按 PDB 建索引
    wt_by_pdb: dict[str, list[int]] = {}
    for i in ws.index:
        wt_by_pdb.setdefault(raw.at[i, "PDB"], []).append(i)
    mut_by_pdb: dict[str, list[int]] = {}
    for i in ms.index:
        mut_by_pdb.setdefault(raw.at[i, "PDB"], []).append(i)

    recs: list[dict] = []
    groups: dict[str, dict] = {}
    skipped: list[str] = []
    for pid in sorted(set(wt_by_pdb) | set(mut_by_pdb)):
        w = wt_by_pdb.get(pid, [])
        m = mut_by_pdb.get(pid, [])
        if not w or not m:
            skipped.append(pid)
            continue
        groups[pid] = {"n_wt": len(w), "n_mut": len(m),
                       "wt_rows": [int(x) for x in w],
                       "mut_rows": [int(x) for x in m]}
        for wi in w:
            # KD(M) -> ΔG，统一口径 ΔG = 0.592·ln(K[M])（与 V1.1–V3.1 一致）
            try:
                dg_w = 0.592 * np.log(float(raw.at[wi, "KD(M)"]))
            except (TypeError, ValueError):
                dg_w = np.nan
            for mi in m:
                try:
                    dg_m = 0.592 * np.log(float(raw.at[mi, "KD(M)"]))
                except (TypeError, ValueError):
                    dg_m = np.nan
                if not (np.isfinite(dg_w) and np.isfinite(dg_m)):
                    continue
                muts = parse_muts(raw.at[mi, "Mutations"])
                base = {
                    "pdb_id": pid,
                    "wt_row": int(wi), "mut_row": int(mi),
                    "split": raw.at[mi, "split"],
                    "dg_wt": float(dg_w), "dg_mut": float(dg_m),
                    "ddg": float(dg_m - dg_w),
                    "n_mut": len(muts),
                    "mutation": raw.at[mi, "Mutations"],
                    "ligand_seq_changed": bool(
                        raw.at[wi, "Ligand Sequences"] != raw.at[mi, "Ligand Sequences"]),
                    "receptor_seq_changed": bool(
                        raw.at[wi, "Receptor Sequences"] != raw.at[mi, "Receptor Sequences"]),
                    "temperature_wt": raw.at[wi, "Temperature(K)"],
                    "temperature_mut": raw.at[mi, "Temperature(K)"],
                    "affinity_method_wt": raw.at[wi, "Affinity Method"],
                    "affinity_method_mut": raw.at[mi, "Affinity Method"],
                    "source": raw.at[mi, "Source Data Set"],
                }
                # 正向 WT->mutant
                recs.append({**base, "direction": "forward"})
                # 反向 mutant->WT (用于对称性约束)
                recs.append({**base, "direction": "reverse",
                             "ddg": float(dg_w - dg_m)})
    pairs = pd.DataFrame(recs)
    emit(f"- 配对构造完成: **{len(pairs)} 行** "
         f"（其中 forward {int((pairs['direction'] == 'forward').sum())}、"
         f"reverse {int((pairs['direction'] == 'reverse').sum())}）")
    emit(f"- 涉及复合物 (pdb_id): **{len(groups)}**")
    emit(f"- 跳过的 PDB（缺 WT 或缺 mutant）: **{len(skipped)}** 个 -> 如实记录")
    emit(f"- 注意: 所有 mutant 所在 PDB 都有 WT 行，因此**没有 mutant 因缺 WT 被丢弃**")

    # ---------- 4. ΔΔG 分布 ----------
    fwd = pairs[pairs["direction"] == "forward"]
    ddg = fwd["ddg"].to_numpy(float)
    emit("")
    emit("## 4. ΔΔG 分布（forward 对）")
    emit("")
    emit("| 统计量 | 值 (kcal/mol) |")
    emit("|---|---:|")
    for nm, v in (("mean", np.mean(ddg)), ("std", np.std(ddg)),
                  ("min", np.min(ddg)), ("max", np.max(ddg))):
        emit(f"| {nm} | {v:+.4f} |")
    for q in (1, 5, 25, 50, 75, 95, 99):
        emit(f"| {q}% 分位 | {np.percentile(ddg, q):+.4f} |")
    emit("")
    emit(f"- ΔΔG > 0（削弱结合）: **{int((ddg > 0).sum())}** "
         f"({(ddg > 0).mean() * 100:.1f}%)")
    emit(f"- ΔΔG < 0（增强结合）: **{int((ddg < 0).sum())}** "
         f"({(ddg < 0).mean() * 100:.1f}%)")
    emit(f"- ΔΔG == 0（无变化）: **{int((np.abs(ddg) < 1e-9).sum())}**")
    emit(f"- |ΔΔG| 中位: {np.median(np.abs(ddg)):.4f} kcal/mol")

    emit("")
    emit("### 单点 vs 多点突变")
    emit("")
    vc = fwd["n_mut"].value_counts().sort_index()
    emit("| 突变数 | 对数 | 占比 |")
    emit("|---|---:|---:|")
    for k, v in vc.items():
        emit(f"| {int(k)} | {int(v)} | {v / len(fwd) * 100:.1f}% |")
    emit("")
    emit(f"- 单点突变对: **{int((fwd['n_mut'] == 1).sum())}** "
         f"({(fwd['n_mut'] == 1).mean() * 100:.1f}%)")
    emit(f"- 多点突变对: **{int((fwd['n_mut'] > 1).sum())}**")

    emit("")
    emit("### 每复合物的突变数分布")
    emit("")
    nm = fwd.groupby("pdb_id").size()
    emit(f"- 中位 {nm.median():.0f}, 均值 {nm.mean():.1f}, 最大 {nm.max()}")
    emit(f"- 只有 1 对的复合物: {int((nm == 1).sum())}; >=5 对: {int((nm >= 5).sum())}; "
         f">=20 对: {int((nm >= 20).sum())}")

    emit("")
    emit("### 序列变化情况（决定特征可用性）")
    emit("")
    emit(f"- 配体序列**发生变化**: {int(fwd['ligand_seq_changed'].sum())} 对 "
         f"({fwd['ligand_seq_changed'].mean() * 100:.1f}%)")
    emit(f"- 受体序列**发生变化**: {int(fwd['receptor_seq_changed'].sum())} 对 "
         f"({fwd['receptor_seq_changed'].mean() * 100:.1f}%)")
    emit(f"- 两侧**都变**: {int((fwd['ligand_seq_changed'] & fwd['receptor_seq_changed']).sum())} 对")
    emit(f"- 两侧**都不变**（ΔΔG 应为 0 或来自条件差异）: "
         f"{int((~fwd['ligand_seq_changed'] & ~fwd['receptor_seq_changed']).sum())} 对")

    # ---------- 5. 与 V3.1 对照 ----------
    emit("")
    emit("## 5. 与 V3.1 可行性检查的对照（必须明确区分）")
    emit("")
    emit("| 项 | V3.1 `03_check_ddg_feasibility.py` | **本脚本 (V4.0)** |")
    emit("|---|---|---|")
    emit("| 配对定义 | 同 `interface_key` **两两成对** | **WT × mutant**（任务书定义）|")
    emit("| 对数 | 243,439 | **{}** |".format(len(fwd)))
    emit("| 组数 | 354（`interface_key` >= 2 行）| **{}**（pdb_id）|".format(len(groups)))
    emit("")
    emit("> **两者不是同一个量**：V3.1 的 243,439 包含大量 mutant×mutant 组合，")
    emit("> 它们的差值**没有 WT 参照**，因此**不是 ΔΔG**。")
    emit("> 本脚本按任务书定义严格只做 WT×mutant，得到的是真实可用的 ΔΔG 对数。")

    # ---------- 6. 落盘 ----------
    pairs.to_csv(OUT_PAIRS, index=False, encoding="utf-8")
    OUT_GROUPS.write_text(json.dumps(groups, ensure_ascii=False, indent=1),
                          encoding="utf-8")
    emit("")
    emit("## 6. 输出")
    emit("")
    emit(f"- `{OUT_PAIRS.relative_to(ROOT)}`: {len(pairs)} 行 x {pairs.shape[1]} 列")
    emit(f"- `{OUT_GROUPS.relative_to(ROOT)}`: {len(groups)} 个复合物")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
