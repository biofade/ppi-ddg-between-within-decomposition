# -*- coding: utf-8 -*-
"""scripts/v5/03b_case_study.py -- V5.0: mutation 可映射性案例研究（只读诊断）.

问题
----
`03a` 显示: 76% 的 pair 至少一侧序列相同；但仍有 22.96% 的 pair 存在
identity <0.70 的一侧。且 annotation（如 `A_P2A`）声明的是**链 ID**，
而序列列是**多链拼接**（`Ligand Sequences` 可能含 `,`）。

本脚本回答三个具体问题（不做映射决定，只报事实）:
Q1. 若只看**序列确实变化的那一侧**，annotation 是否更可能被 recovery？
    (即：把「另一侧相同」的 pair 单独看，其变化侧 identity 如何)
Q2. annotation 的位置（如 P2）在变化侧序列里是否真的是该残基？
    测试 0-based 与 1-based 两种解释。
Q3. 多链拼接（含 `,`）是否解释了失败？分桶比较。

输出: experiments/v5_0/manifests/case_study.md
"""

from __future__ import annotations

import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT = C.V5 / "manifests" / "case_study.md"
MUT_PAT = re.compile(r"^([A-Za-z0-9]+)_([A-Za-z])(\d+)([A-Za-z])$")
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
    emit("# V5.0 — mutation 可映射性案例研究")
    emit("")

    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score = 2.0
    al.mismatch_score = -1.0
    al.open_gap_score = -5.0
    al.extend_gap_score = -1.0

    raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    fwd = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = fwd[fwd["direction"] == "forward"].reset_index(drop=True)

    # ---------- 解析 annotation 的链 ----------
    def parse_ann(s):
        if not isinstance(s, str) or not s.strip():
            return []
        out = []
        for p in re.split(r"[;,]", s):
            mm = MUT_PAT.match(p.strip())
            if mm:
                out.append((mm.group(1), mm.group(2), int(mm.group(3)), mm.group(4)))
        return out

    emit("## Q1/Q3: 分类统计")
    emit("")
    rows = []
    for i, r in fwd.iterrows():
        ann = parse_ann(r["mutation"])
        if not ann:
            continue
        try:
            wi, mi = int(r["wt_row"]), int(r["mut_row"])
        except Exception:  # noqa: BLE001
            continue
        ls_w, ls_m = raw.at[wi, "Ligand Sequences"], raw.at[mi, "Ligand Sequences"]
        rs_w, rs_m = raw.at[wi, "Receptor Sequences"], raw.at[mi, "Receptor Sequences"]
        lig_changed = int(ls_w != ls_m)
        rec_changed = int(rs_w != rs_m)
        lig_multi = int("," in str(ls_w) or "," in str(ls_m))
        rec_multi = int("," in str(rs_w) or "," in str(rs_m))
        ann_chains = sorted({a[0] for a in ann})
        rows.append({
            "pair_row": int(i), "pdb_id": r["pdb_id"], "n_ann": len(ann),
            "lig_changed": lig_changed, "rec_changed": rec_changed,
            "lig_multi": lig_multi, "rec_multi": rec_multi,
            "ann_chains": "|".join(ann_chains),
            "ann_first": f"{ann[0][1]}{ann[0][2]}{ann[0][3]}",
            "n_chains_ann": len(ann_chains),
        })
        if (i + 1) % 2000 == 0:
            emit(f"  ... {i + 1}/{len(fwd)}")
    d = pd.DataFrame(rows)
    emit("")
    emit(f"- 有 annotation 的 pair: **{len(d)} / {len(fwd)}**")
    emit("")
    emit("| 变化模式 | pair 数 | 占比 |")
    emit("|---|---:|---:|")
    pat = d.apply(lambda x: ("lig" if x.lig_changed else "") + ("+rec" if x.rec_changed else "")
                  or "none", axis=1)
    for k, v in pat.value_counts().items():
        emit(f"| {k} | {v} | {v / len(d) * 100:.2f}% |")
    emit("")
    emit("| 多链拼接（含 `,`） | pair 数 | 占比 |")
    emit("|---|---:|---:|")
    emit(f"| 配体侧 | {int(d.lig_multi.sum())} | {d.lig_multi.mean() * 100:.2f}% |")
    emit(f"| 受体侧 | {int(d.rec_multi.sum())} | {d.rec_multi.mean() * 100:.2f}% |")
    emit("")

    # ---------- Q2: annotation 位置能否在变化侧被找到 ----------
    emit("## Q2: annotation 位置在原残基上是否匹配（只测「单侧变化且单突变」子集）")
    emit("")
    sub = d[(d.n_ann == 1) & (d.lig_changed + d.rec_changed == 1)].copy()
    emit(f"- 子集规模: **{len(sub)}** pair（单点突变 + 仅一侧序列变化）")
    emit("")
    hit_1based = hit_0based = 0
    tested = 0
    examples = []
    for _, r in sub.iterrows():
        i = r["pair_row"]
        fr = fwd.iloc[i]
        wi, mi = int(fr["wt_row"]), int(fr["mut_row"])
        side_col = "Ligand Sequences" if r["lig_changed"] else "Receptor Sequences"
        w, m = str(raw.at[wi, side_col]), str(raw.at[mi, side_col])
        if len(w) > 1500 or len(m) > 1500:
            continue
        try:
            a = al.align(w, m)[0]
            aw, am = str(a[0]), str(a[1])
        except Exception:  # noqa: BLE001
            continue
        # 建立 alignment 列 → wt 序列下标映射
        col2wi = {}
        k = 0
        for ci, cw in enumerate(aw):
            if cw != "-":
                col2wi[ci] = k
                k += 1
        # 收集所有差异列
        diffs = [(ci, aw[ci], am[ci]) for ci in range(len(aw))
                 if aw[ci] != "-" and am[ci] != "-" and aw[ci] != am[ci]]
        if len(diffs) != 1:
            continue
        ci, cw, cm = diffs[0]
        wi_idx = col2wi.get(ci)
        if wi_idx is None:
            continue
        tested += 1
        # annotation 的 pos 与 wt_aa
        ann = parse_ann(fr["mutation"])[0]
        pos, wt_aa, mut_aa = ann[2], ann[1].upper(), ann[3].upper()
        if cw.upper() == wt_aa and cm.upper() == mut_aa:
            # alignment 与 annotation 自洽；再定位是 0-based 还是 1-based
            if pos - 1 == wi_idx:
                hit_1based += 1
            elif pos == wi_idx:
                hit_0based += 1
            elif len(examples) < 6:
                examples.append((r["pdb_id"], ann, cw, cm, wi_idx, pos, len(w), len(m)))
    emit(f"- 可测试 pair（alignment 恰好 1 个替换）: **{tested}**")
    emit(f"- 其中 annotation 与该唯一替换**自洽**（原/新残基都对）:")
    emit(f"  * 位置按 **1-based** 解释匹配: **{hit_1based}** （{hit_1based / max(1, tested) * 100:.2f}%）")
    emit(f"  * 位置按 **0-based** 解释匹配: **{hit_0based}** （{hit_0based / max(1, tested) * 100:.2f}%）")
    emit(f"  * 两者都匹配（=第一个残基位置附近）或都不匹配: "
         f"{tested - hit_1based - hit_0based}")
    emit("")
    if examples:
        emit("### 位置仍不匹配的样例（供人工判断，不是映射决定）")
        emit("")
        emit("| pdb | annotation | aln_wt | aln_mut | aln_wt_idx0 | ann_pos | len_wt | len_mut |")
        emit("|---|---|---|---|---:|---:|---:|---:|")
        for e in examples:
            emit(f"| {e[0]} | {e[1][1]}{e[1][2]}{e[1][3]} | {e[2]} | {e[3]} | {e[4]} | "
                 f"{e[5]} | {e[6]} | {e[7]} |")
        emit("")

    emit("## 结论（事实）")
    emit("")
    emit(f"- 「单点突变 + 仅一侧变化」子集规模 **{len(sub)}**（占 {len(sub) / len(d) * 100:.1f}%）")
    emit(f"- 其中 alignment 恰好 1 个替换的 **{tested}** 条里，annotation 残基自洽的 "
         f"**{hit_1based + hit_0based}** 条")
    emit("")
    emit("> 判读提示: 若「自洽」比例高，则**位置编号**只是 offset 问题（可解）；")
    emit("> 若自洽比例低，则两侧序列不是同一 construct，alignment 映射的前提不成立。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    sub.to_csv(C.V5 / "manifests" / "case_study_subset.csv", index=False, encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
