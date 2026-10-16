# -*- coding: utf-8 -*-
"""scripts/v4.0/02_build_ddg_split.py -- V4.0 任务 2: ΔΔG 抗泄漏 split.

为什么必须**重建** split（而不是沿用 Proteinea 划分）
-----------------------------------------------------
ΔΔG 配对来自 SKEMPI v2.0，且集中在 231 个复合物上。实测继承划分完全退化：

| split | 配对数 | 复合物数 |
|---|---:|---:|
| train | 7,437 | 208 |
| val | 62 | 13 |
| **test** | **27** | **10** |

**27 个 test 对 / 10 个复合物完全无法支撑评估** —— 既算不出稳定的
per-structure Pearson（多数复合物只有 1–3 对），也做不了配对 bootstrap。
因此本任务重建划分，并**如实报告**继承划分的退化。

三层隔离（按任务书要求）
------------------------
1. **复合物级**：同一 `pdb_id` 的全部配对必须在同一侧（硬保证）。
2. **序列同一性**：以 40% 为阈值做序列聚类。
   * **精确匹配**（exact sequence equality）是**零歧义**的强泄漏 -> 必查；
   * 3-mer Jaccard 作为**同源代理**（V3.1 已证明它会**高估**真实同一性），
     取 **保守方向**：只要 Jaccard > 0.40 就当成同源，强制同侧。
     > 这会把真实同一性 <40% 的序列也判为同源（过度保守）。
     > 代价是 split 更小，收益是**不会**低估泄漏。报告中明确标注该取舍。
3. **家族级**：用 (2) 得到的 cluster 作为家族，确保不跨侧。

划分算法（确定性）
------------------
* 先做序列聚类得到 cluster；把复合物映射到 cluster（任一链落入即算同族）。
* cluster 按总配对数降序排列，**贪心地**分配到当前"最欠额"的 split，
  目标比例 70/15/15（按配对数）。降序 + 最欠额可避免大 cluster 集中到一侧。
* 全程无随机数 -> **完全可复现**。

输出
----
* `data/v4.0/ddg_split.json`   三划分的 `pair_index`（指向 ddg_pairs.csv 的行号）
                                + cluster 映射 + 规则说明
* `features/v4.0_ddg_split_report.md`
"""

from __future__ import annotations

import argparse
import json
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
V40 = ROOT / "data" / "v4.0"
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v4.0"
LOG_DIR.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "failures.log"

PAIRS = V40 / "ddg_pairs.csv"
OUT_JSON = V40 / "ddg_split.json"
OUT_MD = FEAT / "v4.0_ddg_split_report.md"

SEQ_THRESH = 0.40          # 任务书要求的序列同一性阈值
TARGET = {"train": 0.70, "val": 0.15, "test": 0.15}
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


def build_sequence_clusters(seqs: set[str], thresh: float) -> dict[str, int]:
    """单链接聚类：Jaccard > thresh 或序列完全相同 -> 同簇。返回 seq -> cluster。"""
    uniq = sorted(seqs)
    km = {s: kmers(s) for s in uniq}
    parent = {s: s for s in uniq}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, a in enumerate(uniq):
        ka = km[a]
        for b in uniq[i + 1:]:
            kb = km[b]
            if a == b:
                union(a, b)
                continue
            inter = len(ka & kb)
            if inter == 0:
                continue
            if inter / len(ka | kb) > thresh:
                union(a, b)
    roots: dict[str, int] = {}
    out: dict[str, int] = {}
    for s in uniq:
        r = find(s)
        if r not in roots:
            roots[r] = len(roots)
        out[s] = roots[r]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 任务 2: ΔΔG 抗泄漏 split")
    ap.add_argument("--pairs", type=Path, default=PAIRS)
    ap.add_argument("--thresh", type=float, default=SEQ_THRESH)
    args = ap.parse_args()

    t0 = time.time()
    emit("# V4.0 任务 2 — ΔΔG 抗泄漏 split")
    emit("")
    emit(f"**序列同一性阈值**: {args.thresh:.2f}（任务书要求 40%）")
    emit("**目标比例**: 70 / 15 / 15（按配对数）")
    emit("")

    pairs = pd.read_csv(args.pairs, low_memory=False)
    fwd = pairs[pairs["direction"] == "forward"].reset_index(drop=True)
    emit("## 1. 输入")
    emit("")
    emit(f"- `ddg_pairs.csv`: {len(pairs)} 行；forward **{len(fwd)}** 对；"
         f"复合物 **{fwd['pdb_id'].nunique()}**")
    emit("")

    # ---------- 2. 继承划分的退化（如实报告）----------
    emit("## 2. 继承划分（Proteinea `split`）的退化 — 必须记录")
    emit("")
    emit("| split | 配对数 | 复合物数 |")
    emit("|---|---:|---:|")
    for s in ("train", "val", "test"):
        sub = fwd[fwd["split"] == s]
        emit(f"| {s} | {len(sub)} | {sub['pdb_id'].nunique()} |")
    n_test = int((fwd["split"] == "test").sum())
    emit("")
    emit(f"> ⚠️ 继承划分的 test 只有 **{n_test} 对 / "
         f"{fwd[fwd['split'] == 'test']['pdb_id'].nunique()} 个复合物** ——")
    emit("> 无法支撑 per-structure Pearson 或配对 bootstrap。**因此重建划分。**")
    emit(f"> （此划分本身是**复合物级隔离**的：跨 split 的 pdb_id = "
         f"{int((fwd.groupby('pdb_id')['split'].nunique() > 1).sum())}，无泄漏，只是太小。）")
    emit("")

    # 序列（用 mutant 侧的序列，因为要预测的就是 mutant 复合物）
    lig: dict[str, str] = {}
    rec: dict[str, str] = {}
    for pid, sub in fwd.groupby("pdb_id"):
        lig[pid] = str(sub.iloc[0]["pdb_id"])
        # 序列不在 pairs 里 -> 从原始表取（pairs 只存了是否变化）
    # 需要从 v2.2 特征表取序列
    v22 = pd.read_csv(FEAT / "v2.2_proteinea_features.csv", low_memory=False)
    v22["pdb_id"] = v22["pdb_id"].astype(str)
    seqmap = (v22.groupby("pdb_id")
              .agg(lig=("Ligand Sequences", "first"),
                   rec=("Receptor Sequences", "first")))
    miss = [p for p in fwd["pdb_id"].unique() if p not in seqmap.index]
    emit(f"## 3. 序列来源")
    emit("")
    emit(f"- 从 `v2.2_proteinea_features.csv` 取每个 `pdb_id` 的代表序列")
    emit(f"- 覆盖: {len(fwd['pdb_id'].unique()) - len(miss)} / "
         f"{fwd['pdb_id'].nunique()} 个复合物"
         + (f"（缺失 {len(miss)} 个，如实报告）" if miss else ""))
    emit("")

    # ---------- 4. 序列聚类 ----------
    emit("## 4. 序列聚类（40% 阈值，两侧分别聚类后合并）")
    emit("")
    pids = sorted(p for p in fwd["pdb_id"].unique() if p in seqmap.index)
    lig_seqs = {p: str(seqmap.at[p, "lig"]) for p in pids}
    rec_seqs = {p: str(seqmap.at[p, "rec"]) for p in pids}
    t1 = time.time()
    lig_cl = build_sequence_clusters(set(lig_seqs.values()), args.thresh)
    rec_cl = build_sequence_clusters(set(rec_seqs.values()), args.thresh)
    emit(f"- 配体侧: {len(set(lig_seqs.values()))} 条唯一序列 -> "
         f"**{len(set(lig_cl.values()))} 簇**")
    emit(f"- 受体侧: {len(set(rec_seqs.values()))} 条唯一序列 -> "
         f"**{len(set(rec_cl.values()))} 簇**")
    emit(f"- 聚类耗时 {time.time() - t1:.0f}s")

    # 复合物 -> 家族。**必须用并查集合并"共享任一链簇"的复合物**：
    # 若简单取 (lig_cl, rec_cl) 组合作为家族 id，则 L8_R67 与 L34_R67 会被当成
    # 两个家族 —— 但它们**共享受体簇 R67**（同一溶菌酶），是同一同源家族。
    # 实测该 bug 会让 6 条完全相同的受体序列跨 train/test（零歧义泄漏）。
    fam = {}
    _parent: dict[str, str] = {}

    def _find(x: str) -> str:
        _parent.setdefault(x, x)
        while _parent[x] != x:
            _parent[x] = _parent[_parent[x]]
            x = _parent[x]
        return x

    def _union(a: str, b: str) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            _parent[rb] = ra

    for p in pids:
        _union(f"L{lig_cl[lig_seqs[p]]}", f"R{rec_cl[rec_seqs[p]]}")
    for p in pids:
        fam[p] = _find(f"L{lig_cl[lig_seqs[p]]}")
    n_fam = len(set(fam.values()))
    emit(f"- 复合物家族（**共享任一链簇即同族**，并查集连通分量）: **{n_fam}** 个")
    # 检查精确重复
    dup_lig = len(set(lig_seqs.values())) - len(set(lig_cl.values()))
    emit("")
    emit("> 注意: 3-mer Jaccard 在 V3.1 已证明会**高估**真实同一性（同源折叠、")
    emit("> 保守基序），因此把它当同源判据是**保守方向** —— 可能把真实同一性")
    emit("> <40% 的序列也判为同源，使 split 更小但**不会低估泄漏**。")
    emit("")

    # ---------- 5. 贪心分配（确定性）----------
    emit("## 5. 划分分配（确定性贪心）")
    emit("")
    fam_size: dict[str, int] = {}
    fam_pids: dict[str, list[str]] = {}
    for p, f_ in fam.items():
        fam_size[f_] = fam_size.get(f_, 0) + int((fwd["pdb_id"] == p).sum())
        fam_pids.setdefault(f_, []).append(p)
    order = sorted(fam_size, key=lambda k: (-fam_size[k], k))   # 大族优先
    total = sum(fam_size.values())
    quota = {s: TARGET[s] * total for s in TARGET}
    got = {s: 0 for s in TARGET}
    assign: dict[str, str] = {}
    for f_ in order:
        # 选"最欠额"的 split（got/target 最小）；并列时按 train>val>test
        best = min(("train", "val", "test"),
                   key=lambda s: (got[s] - quota[s], s))
        assign[f_] = best
        got[best] += fam_size[f_]
    emit("| split | 目标对数 | 实际对数 | 实际占比 | 家族数 | 复合物数 |")
    emit("|---|---:|---:|---:|---:|---:|")
    for s in ("train", "val", "test"):
        fs = [f_ for f_, a in assign.items() if a == s]
        np_ = len([p for f_ in fs for p in fam_pids[f_]])
        emit(f"| {s} | {quota[s]:.0f} | {got[s]} | {got[s] / total * 100:.1f}% | "
             f"{len(fs)} | {np_} |")
    emit("")

    # 行索引
    split_rows: dict[str, list[int]] = {}
    pid2split = {p: assign[fam[p]] for p in pids}
    for s in ("train", "val", "test"):
        idx = fwd.index[fwd["pdb_id"].map(pid2split) == s].tolist()
        split_rows[s] = [int(i) for i in idx]

    # ---------- 6. 验证 ----------
    emit("## 6. Split 验证（全部用精确判据）")
    emit("")
    sets = {s: set(fwd.loc[split_rows[s], "pdb_id"]) for s in split_rows}
    emit("### 6.1 复合物级隔离（精确）")
    emit("")
    emit("| 交集 | pdb_id 数 |")
    emit("|---|---:|")
    for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
        emit(f"| {a} ∩ {b} | **{len(sets[a] & sets[b])}** |")
    ok_comp = all(len(sets[a] & sets[b]) == 0
                  for a, b in (("train", "val"), ("train", "test"), ("val", "test")))
    emit("")
    emit(f"-> {'✅ 完全隔离' if ok_comp else '❌ **存在复合物级泄漏**'}")
    emit("")

    emit("### 6.2 序列同一性隔离（40% 阈值）")
    emit("")
    fam_sets = {s: {fam[p] for p in sets[s]} for s in sets}
    emit("| 交集 | 家族数 |")
    emit("|---|---:|")
    for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
        emit(f"| {a} ∩ {b} | **{len(fam_sets[a] & fam_sets[b])}** |")
    ok_fam = all(len(fam_sets[a] & fam_sets[b]) == 0
                 for a, b in (("train", "val"), ("train", "test"), ("val", "test")))
    emit("")
    emit(f"-> {'✅ 完全隔离（无 >40% 同源跨侧）' if ok_fam else '❌ **存在同源跨侧**'}")

    # 精确重复序列的跨界（零歧义检查）
    emit("")
    emit("### 6.3 精确重复序列跨侧检查（零歧义）")
    emit("")
    exact_leak = 0
    for side, seqmap_ in (("配体", lig_seqs), ("受体", rec_seqs)):
        bad = {}
        for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
            sa = {seqmap_[p] for p in sets[a]}
            sb = {seqmap_[p] for p in sets[b]}
            bad[f"{a}-{b}"] = len(sa & sb)
            exact_leak += len(sa & sb)
        emit(f"- {side}侧精确相同序列跨侧数: {bad}")
    emit("")
    # 硬断言: 完全相同序列**不可能**合法跨侧（它是 100% 同一性，必 >40%）
    if exact_leak > 0:
        emit(f"> ❌ **硬断言失败**: 有 {exact_leak} 条完全相同的序列跨侧。")
        emit("> 完全相同序列的同一性 = 100% > 40%，必然属于同一家族，")
        emit("> 因此这**一定**是家族定义或聚类的 bug，而不是可接受的现象。")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[02_build_ddg_split] exact-duplicate cross-side leak = "
                     f"{exact_leak}\n")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 5
    emit("> ✅ 无任何完全相同的序列跨侧（最高强度的泄漏形式已排除）。")
    emit("")

    emit("### 6.4 各 split 的 ΔΔG 分布")
    emit("")
    emit("| split | n | mean | std | 中位 | \\|ΔΔG\\|中位 | >0 占比 | 单点占比 |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    for s in ("train", "val", "test"):
        sub = fwd.loc[split_rows[s]]
        d = sub["ddg"].to_numpy(float)
        emit(f"| {s} | {len(sub)} | {d.mean():+.3f} | {d.std():.3f} | "
             f"{np.median(d):+.3f} | {np.median(np.abs(d)):.3f} | "
             f"{(d > 0).mean() * 100:.0f}% | "
             f"{(sub['n_mut'] == 1).mean() * 100:.0f}% |")
    emit("")

    # ---------- 7. 落盘 ----------
    out = {
        "pairs_file": "data/v4.0/ddg_pairs.csv",
        "definition": "ddG = dG_mutant - dG_wildtype ; forward pairs only",
        "rules": {
            "complex_level": "all pairs of one pdb_id on the same side (verified disjoint)",
            "sequence_identity_threshold": args.thresh,
            "sequence_identity_metric": (
                "single-linkage on 3-mer Jaccard; EXACT sequence equality always merged. "
                "NOTE: 3-mer Jaccard overestimates true identity (V3.1 finding), so this "
                "is the conservative direction (may over-merge, never under-merges)."),
            "family_level": "family = (ligand cluster, receptor cluster); no family crosses sides",
            "assignment": "deterministic greedy, families sorted by pair count desc, "
                          "assigned to most-under-quota split; target 70/15/15",
            "no_randomness": True,
        },
        "counts": {s: len(split_rows[s]) for s in split_rows},
        "complexes": {s: sorted(sets[s]) for s in sets},
        "family_of_complex": {p: fam[p] for p in pids},
        "pair_index": split_rows,
        "inherited_split_counts": {
            s: int((fwd["split"] == s).sum()) for s in ("train", "val", "test")},
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    emit("## 7. 输出")
    emit("")
    emit(f"- `{OUT_JSON.relative_to(ROOT)}`")
    emit(f"- `{OUT_MD.relative_to(ROOT)}`")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
