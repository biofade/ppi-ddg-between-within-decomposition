# -*- coding: utf-8 -*-
"""scripts/v7_4r2/01_audit_structure_preparation.py -- V7.4R2 PHASE 1：纯 Python 结构审计（§6）.

**不需要 OpenMM**。在提交 Colab 之前，本地先回答：

1. 每个 target 有哪些**链内缺口**（残基编号断层）？PDBFixer 会在这里插入缺失残基
   —— V7.4R2 **禁止**该行为（§5），因此缺口附近的突变位点必须**排除**。
2. 突变残基的**骨架原子是否齐全**（N/CA/C/O）？
3. 突变残基是否有**插入码**（§6 排除）？
4. 同一 target 的**观测序列**是否与 `mutation_sites.wt_aa` 一致？

产出 `CANONICAL_WT_STRICT`：只保留**缺口安全 + 骨架完整 + 无插入码 + WT 一致**的行。

这是纯静态检查，成本秒级，**把上一轮要烧两小时才暴露的问题前移到本地**。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
PH1 = R / "experiments" / "v7_4r2" / "phase1_audit"
PH1.mkdir(parents=True, exist_ok=True)
FEAT = R / "features" / "v7_4r2"
FEAT.mkdir(parents=True, exist_ok=True)
PDBDIR = R / "data" / "pdb"
STD3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
BACKBONE = ("N", "CA", "C", "O")
L: list[str] = []

# 缺口安全半径：位点必须在距最近缺口 >= GAP_SAFE 个残基处
GAP_SAFE = 3


def emit(m=""):
    print(m); L.append(str(m))


def parse_pdb(p: Path):
    """返回 chain -> list[((resseq,icode,resname), {atom:xyz})]（保持文件顺序）。

    注意：必须按**上一次追加的残基键**判断是否新残基；此前误用
    `res[c][-1][0]` 与 `key` 比较（两者都是 tuple），比较恒为 False，
    导致每个原子都新建一个"残基" —— 这是上一版的 bug。
    """
    res = defaultdict(list)
    for ln in open(p, "r", errors="replace"):
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        rn = ln[17:20].strip().upper()
        if rn not in STD3:
            continue
        c = ln[21]
        try:
            rs = int(ln[22:26])
        except ValueError:
            continue
        ic = ln[26].strip()
        key = (rs, ic, rn)
        if not res[c] or res[c][-1][0] != key:      # 两者都是 tuple -> 正确比较
            res[c].append((key, {}))
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el != "H":
            res[c][-1][1][ln[12:16].strip()] = (
                float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return res


def chain_gaps(entries):
    """返回 (gap_positions, break_indices)。

    gap = 相邻残基的**编号差 != 1**（含不连续编号）或链内序号不连续。
    返回每个"链内下标"是否靠近缺口（距离 < GAP_SAFE）。
    """
    n = len(entries)
    near = [False] * n
    gaps = []
    for i in range(1, n):
        prev = entries[i - 1][0][0]
        cur = entries[i][0][0]
        step = cur - prev
        if step != 1:
            gaps.append({"after_index": i - 1, "prev_resseq": prev,
                         "next_resseq": cur, "gap_size": step - 1})
            for j in range(max(0, i - GAP_SAFE), min(n, i + GAP_SAFE)):
                near[j] = True
    return gaps, near


def main() -> int:
    sp = importlib.util.spec_from_file_location(
        "v7f", R / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    T2O = v7f.THREE2ONE            # 三字母 -> 一字母
    ONE2THREE = {v: k for k, v in T2O.items()}   # 一字母 -> 三字母

    emit("# V7.4R2 PHASE 1 — 结构准备审计（纯 Python，无需 OpenMM）")
    emit("")
    emit(f"- 缺口安全半径 `GAP_SAFE = {GAP_SAFE}` 残基（位点距最近缺口需 ≥ {GAP_SAFE}）")
    emit("")

    man = pd.read_csv(R / "features" / "v7_4" /
                      "u_relax_strict_manifest.csv", low_memory=False)
    site = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    emit(f"- 起点 U_RELAX_STRICT：**{len(man)}** 行 / "
         f"{man.target_id_canonical.nunique()} targets")

    cache = {}
    gap_done = set()
    exc = defaultdict(int)
    exc_rows, keep_rows = [], []
    tgt_stats = []
    for r in man.itertuples():
        t = str(r.target_id_canonical)
        p = PDBDIR / f"{t.upper()}.pdb"
        if t not in cache:
            cache[t] = parse_pdb(p) if p.exists() else None
        P = cache[t]
        if P is None:
            exc["missing_pdb"] += 1
            exc_rows.append({"row_id": int(r.row_id), "target": t,
                             "reason": "missing_pdb"}); continue
        # 逐链缺口（**每个 target 只算一次**）
        if t not in gap_done:
            gap_done.add(t)
            gapinfo = {c: chain_gaps(P[c]) for c in P}
            tgt_stats.append({"target": t, "n_chains": len(P),
                              "n_gaps": sum(len(gg) for gg, _ in gapinfo.values()),
                              "n_res": sum(len(P[c]) for c in P)})
        gapinfo = {c: chain_gaps(P[c]) for c in P}
        g = s_by_row.get(int(r.row_id))
        if g is None or not len(g):
            exc["no_site"] += 1
            exc_rows.append({"row_id": int(r.row_id), "target": t,
                             "reason": "no_site"}); continue
        bad = None
        for x in g.itertuples():
            c = str(x.ann_chain)
            if c not in P:
                bad = "missing_chain"; break
            # **用 seq_index0 寻址**（V7.3 的成功做法）
            # 实测：seq_index0 -> 100.0% WT 匹配；pdb_resseq -> 仅 78.4%
            # （pdb_resseq 带每链 numbering_offset，见 mutation_sites.numbering_offset）
            idx = int(x.seq_index0)
            if not (0 <= idx < len(P[c])):
                bad = "site_index_out_of_range"; break
            key, atoms = P[c][idx]
            if key[2] != ONE2THREE.get(str(x.wt_aa), "?"):
                bad = "wt_identity_mismatch"; break
            if key[1] != "":
                bad = "insertion_code_at_site"; break
            if not all(a in atoms for a in BACKBONE):
                bad = "site_backbone_incomplete"; break
            _, near = gapinfo[c]
            if near[idx]:
                bad = "site_near_chain_gap"; break
        if bad:
            exc[bad] += 1
            exc_rows.append({"row_id": int(r.row_id), "target": t, "reason": bad})
            continue
        keep_rows.append({"row_id": int(r.row_id),
                          "target_id_canonical": t,
                          "family_id_canonical": str(r.family_id_canonical),
                          "n_sites": int(r.n_sites),
                          "mutated_chains": str(r.mutated_chains),
                          "mutation_bodies": str(r.mutation_bodies)})

    K = pd.DataFrame(keep_rows)
    emit("")
    emit("## 1.1 CANONICAL_WT_STRICT（§6 严格结构排除）")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| 起点 rows | {len(man)} |")
    emit(f"| **CANONICAL_WT_STRICT rows** | **{len(K)}** |")
    emit(f"| targets | **{K.target_id_canonical.nunique() if len(K) else 0}** |")
    emit(f"| families | **{K.family_id_canonical.nunique() if len(K) else 0}** |")
    for k in (1, 2, 3, 5, 10):
        v = int((K.groupby('target_id_canonical').size() >= k).sum()) if len(K) else 0
        emit(f"| n>={k} targets | **{v}** |")
    emit("")
    emit("### 排除原因（§6）")
    emit("")
    emit("| reason | rows |")
    emit("|---|---:|")
    for k, v in sorted(exc.items(), key=lambda x: -x[1]):
        emit(f"| `{k}` | {v} |")
    emit(f"| **合计排除** | **{sum(exc.values())}** |")
    emit("")

    TS = pd.DataFrame(tgt_stats)
    if len(TS):
        emit("## 1.2 逐 target 缺口统计（PDBFixer 会在此插入缺失残基）")
        emit("")
        emit(f"- target 数：**{len(TS)}**")
        emit(f"- **含 ≥1 链内缺口的 target：{int((TS.n_gaps > 0).sum())} / {len(TS)}**")
        emit(f"- 缺口总数：**{int(TS.n_gaps.sum())}**")
        emit(f"- 缺口分布：median={TS.n_gaps.median():.0f} "
             f"max={TS.n_gaps.max():.0f}")
        emit("")
        emit("> V7.4R2 **禁止**填补缺失残基（§5）。缺口附近的位点已按 "
             f"`GAP_SAFE={GAP_SAFE}` 安全半径排除。")
    emit("")

    K.to_csv(FEAT / "canonical_wt_strict_manifest.csv", index=False, encoding="utf-8")
    if exc_rows:
        pd.DataFrame(exc_rows).to_csv(PH1 / "structural_exclusions.csv",
                                      index=False, encoding="utf-8")
    TS.to_csv(PH1 / "target_gap_stats.csv", index=False, encoding="utf-8")
    json.dump({"n_rows": len(K),
               "n_targets": int(K.target_id_canonical.nunique()) if len(K) else 0,
               "n_families": int(K.family_id_canonical.nunique()) if len(K) else 0,
               "gap_safe": GAP_SAFE,
               "exclusions": dict(exc),
               "n_targets_with_gaps": int((TS.n_gaps > 0).sum()) if len(TS) else 0},
              open(PH1 / "audit_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    (PH1 / "audit_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n>>> CANONICAL_WT_STRICT rows={len(K)} "
          f"targets={K.target_id_canonical.nunique() if len(K) else 0}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
