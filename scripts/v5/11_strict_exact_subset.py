# -*- coding: utf-8 -*-
"""scripts/v5/11_strict_exact_subset.py -- V5.0: 从 residue_map 抽出**严格 EXACT** 子集.

背景（承接 GATE B 的两次迭代，全部为 PROJECT_EVIDENCE）
-------------------------------------------------------
* `04` v1: 自报覆盖 74.06%，但 `05` 用真实 PDB 残基号独立验证 → E1 通过率 **7.96%**
  ⇒ v1 `INVALID_DUE_TO_BUG`。
* `06` v2（序列引导）：残基一致率升到 **76.82%**，但 offset 分布仍显示假匹配
  （`|offset|` 最大 1745）。说明「残基自洽」不足以排除偶然匹配。

因此本脚本按**最严条件**筛：
  `validated_status ∈ {EXACT, ALIGNMENT_CONFIDENT}`
  **且** `pdb_resseq_real > 0`（成功定位到真实 PDB 残基）
  **且** `real_aa == wt_aa`（该位置真实氨基酸 == annotation 声明的 WT 氨基酸）
  **且** `ann_pos1 == pdb_resseq_real`（**编号直接一致，offset = 0**）

→ 只保留「位置编号 + 氨基酸身份」双重精确匹配的位点。
这是**最保守**的口径；宁可覆盖低，也不要让 MIPA 建立在错位坐标上（V5 §23）。

输出
----
* `data/v5_0/strict_exact_sites.csv`
* `data/v5_0/strict_exact_index.json`（含可用于实验的 pair 列表）
* `experiments/v5_0/manifests/strict_exact_report.md`
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.DATA / "v5_0" / "strict_exact_sites.csv"
OUT_JSON = C.DATA / "v5_0" / "strict_exact_index.json"
OUT_MD = C.V5 / "manifests" / "strict_exact_report.md"
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
    emit("# V5.0 — 严格 EXACT 位点子集（GATE B 的最保守口径）")
    emit("")
    emit("> 筛选条件: `validated_status ∈ {EXACT, ALIGNMENT_CONFIDENT}` **且**")
    emit("> `pdb_resseq_real > 0` **且** `real_aa == wt_aa` **且** `ann_pos1 == pdb_resseq_real`。")
    emit("> 目的: 排除 `|offset|` 极大的偶然匹配（GATE B 已实测存在该现象）。")
    emit("")

    rm = pd.read_csv(C.DATA / "v5_0" / "residue_map.csv", low_memory=False)
    n0 = len(rm)
    emit(f"- 输入位点: **{n0}**")
    emit("")

    rm["ann_pos1"] = pd.to_numeric(rm["ann_pos1"], errors="coerce")
    rm["pdb_resseq_real"] = pd.to_numeric(rm["pdb_resseq_real"], errors="coerce")
    c1 = rm["validated_status"].isin(["EXACT", "ALIGNMENT_CONFIDENT"])
    c2 = rm["pdb_resseq_real"] > 0
    c3 = rm["real_aa"].astype(str) == rm["wt_aa"].astype(str)
    c4 = rm["ann_pos1"] == rm["pdb_resseq_real"]
    emit("| 条件 | 剩余位点 |")
    emit("|---|---:|")
    m = pd.Series(True, index=rm.index)
    for name, c in (("validated_status ∈ {EXACT, ALIGNMENT_CONFIDENT}", c1),
                    ("pdb_resseq_real > 0", c2),
                    ("real_aa == wt_aa（氨基酸身份一致）", c3),
                    ("ann_pos1 == pdb_resseq_real（offset=0）", c4)):
        m &= c
        emit(f"| + {name} | {int(m.sum())} |")
    emit("")
    sub = rm[m].copy()
    emit(f"- **严格 EXACT 位点: {len(sub)} / {n0} = {len(sub) / max(1, n0) * 100:.2f}%**")
    emit(f"- 覆盖的唯一 pair: **{sub['pair_row'].nunique()}**")
    emit(f"- 覆盖的唯一 PDB: {sub['pdb_id'].nunique()}")
    emit(f"- 覆盖的 side: {sub['side'].value_counts().to_dict()}")
    emit("")

    # 每条 pair 的注释数与严格位点数是否相等（要求「全注释都严格映射」）
    pmap = pd.read_csv(C.DATA / "v5_0" / "mutation_pair_map.csv", low_memory=False)
    annn = dict(zip(pmap["pair_row"].astype(int), pmap["n_ann"].astype(int)))
    cnt = sub.groupby("pair_row").size().to_dict()
    full = [int(p) for p, c in cnt.items() if c == annn.get(int(p), -1)]
    # 该 pair 只能有一条链变化（MIPA 需要区分哪条链是突变链）
    d = C.build_pair_features(verbose=False)
    df = d["df"]
    lig_ch = df["ligand_seq_changed"].astype(bool).to_numpy()
    rec_ch = df["receptor_seq_changed"].astype(bool).to_numpy()
    n_usable_rows = len(df)
    # pair_row 来自**未映射**的 7526 行表；build_pair_features 丢弃 3 对后只剩
    # n_usable_rows 行 → 必须滤掉越界索引（否则 IndexError）
    dropped = [p for p in full if p >= n_usable_rows]
    if dropped:
        emit(f"> ⚠️ 有 {len(dropped)} 条 pair 的 `pair_row` 超出映射后的行数"
             f"（{n_usable_rows}）—— 这些属于 V4 映射阶段被丢弃的 3 对，已排除。")
        emit("")
    full = [p for p in full if p < n_usable_rows]
    assert all(0 <= p < n_usable_rows for p in full), "pair_row 越界"
    one_side = [p for p in full if lig_ch[p] ^ rec_ch[p]]
    pdb_of = dict(zip(df.index.astype(int), df["pdb_id"].astype(str)))
    n_pdb_full = len({pdb_of.get(p, "") for p in full})
    n_pdb_one = len({pdb_of.get(p, "") for p in one_side})
    emit("### 收紧到「可做 MIPA」的条件")
    emit("")
    emit("| 条件 | pair 数 | 涉及 PDB |")
    emit("|---|---:|---:|")
    emit(f"| 该 pair 的**全部** annotation 都严格映射 | **{len(full)}** | {n_pdb_full} |")
    emit(f"| + 仅**一侧**序列变化（可判定突变链） | **{len(one_side)}** | {n_pdb_one} |")
    emit("")
    emit("> 「两侧都变化」的 pair 无法确定某条 annotation 属于哪条链，因此不入 MIPA 主实验。")
    emit("")

    # 与 split 交叉
    spl = json.loads((C.DATA / "v4.0" / "ddg_split.json").read_text(encoding="utf-8"))
    in_tr = set(int(x) for x in spl["pair_index"]["train"])
    in_te = set(int(x) for x in spl["pair_index"]["test"])
    in_va = set(int(x) for x in spl["pair_index"]["val"])
    cnt_split = {"train": len(set(one_side) & in_tr),
                 "val": len(set(one_side) & in_va),
                 "test": len(set(one_side) & in_te)}
    emit("### 严格 EXACT 单侧突变 pair 的 split 分布")
    emit("")
    emit("| split | pair 数 |")
    emit("|---|---:|")
    for k, v in cnt_split.items():
        emit(f"| {k} | {v} |")
    emit("")
    usable_test = sorted(set(one_side) & in_te)
    emit(f"- **test 上可用于 MIPA 的 pair: {len(usable_test)}** "
         f"（占 test {len(in_te)} 的 {len(usable_test) / max(1, len(in_te)) * 100:.1f}%）")
    emit("")

    # GATE B 最终判定
    cover_all = len(one_side) / max(1, n0 / max(1, 1))
    n_pairs_total = len(df)
    cov_pairs = len(one_side) / max(1, n_pairs_total)
    emit("## GATE B 最终判定（严格口径）")
    emit("")
    emit(f"- 可用于 MIPA 的 pair 覆盖率 = **{len(one_side)} / {n_pairs_total} = "
         f"{cov_pairs * 100:.2f}%**")
    emit("")
    if cov_pairs >= 0.90:
        v = "FULL_MIPA"
    elif cov_pairs >= 0.70:
        v = "SUBSET_PILOT"
    else:
        v = "SOLVE_MAPPING_FIRST"
    emit(f"**GATE B（严格 EXACT）= `{v}`**")
    emit("")
    emit("> 结论：按 V5 §124，覆盖率不足 70% ⇒ **不做 full MIPA**；")
    emit("> 只做**严格 EXACT 子集 pilot**，并且**必须**在相同 intersection rows 上")
    emit("> 重跑 baseline 做 coverage fairness（V5 §34）。")
    emit("")

    sub.to_csv(OUT_CSV, index=False, encoding="utf-8")
    OUT_JSON.write_text(json.dumps({
        "rule": ("validated_status in {EXACT,ALIGNMENT_CONFIDENT} AND pdb_resseq_real>0 "
                 "AND real_aa==wt_aa AND ann_pos1==pdb_resseq_real"),
        "n_sites_input": int(n0), "n_sites_strict": int(len(sub)),
        "pairs_all_annotations_mapped": sorted(int(x) for x in full),
        "pairs_single_side_changed": sorted(int(x) for x in one_side),
        "pairs_split": {k: sorted(set(one_side) & s)
                        for k, s in (("train", in_tr), ("val", in_va), ("test", in_te))},
        "pair_coverage": cov_pairs, "verdict": v,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    for p in (OUT_CSV, OUT_JSON, OUT_MD):
        print(f"写出 {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
