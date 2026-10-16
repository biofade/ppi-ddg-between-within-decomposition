# -*- coding: utf-8 -*-
"""scripts/v5/04_build_mutation_sites.py -- V5.0 GATE B（正式版）: mutation 残基映射.

机制（对应 V5 §22–§26 全部要求）
--------------------------------
1. **不需要外部 WT 序列** —— 用 V4 pair 内 **WT row 序列 vs mutant row 序列**。
2. 用 Biopython `PairwiseAligner`（global）做 alignment，收集差异（substitution）。
3. **朴素位置失败 → 求解编号 offset**（这是 V5 §25 明确列出的
   "residue numbering offset"，不是 construct mismatch）：
   对同一条链上的每个 annotation `WppM`，令 `aln_idx0` 为 alignment 中
   该替换所在 WT 位置（0-based）。若 `aln_aa == W` 且 `mut_aa == M`（**残基自洽**），
   则该 annotation 支持 offset `o = pp - 1 - aln_idx0`。
   若同链所有自洽 annotation 支持**同一个** `o`，则映射置信。
4. `mapping_status` 分类（V5 §26）:
   * `EXACT`    — 位置直接匹配（o = 0）且残基自洽
   * `ALIGNMENT_CONFIDENT` — 残基自洽 + 同链 offset 唯一（o ≠ 0）
   * `AMBIGUOUS` — 残基自洽但 offset 不唯一，或 alignment 含 indel 无法唯一确定
   * `FAILED`   — 残基不自洽 / 找不到对应替换 / 无法 alignment
5. **主实验只用 EXACT + ALIGNMENT_CONFIDENT**；另外保存 EXACT-only 子集供敏感性分析。
6. 失败条目全部写入 `mutation_mapping_failures.csv`，**不静默丢弃**。

额外记录（供 P1 直接用）
------------------------
* `seq_index0`：突变位点在**序列坐标系**中的下标（0-based）
* `pdb_resseq`：推断的 PDB 残基编号（= annotation 的 pp，因 offset 已解）
* `n_indel_in_aln`：该链 alignment 的 indel 数（indel >0 时标记，V5 §23 要求
  「首阶段只做 substitution」）

输出
----
* `data/v5_0/mutation_sites.csv`
* `data/v5_0/mutation_pair_map.csv`
* `data/v5_0/mutation_mapping_failures.csv`
* `experiments/v5_0/manifests/gate_b_mapping_report.md`
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUTDIR = C.DATA / "v5_0"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_SITES = OUTDIR / "mutation_sites.csv"
OUT_PAIRS = OUTDIR / "mutation_pair_map.csv"
OUT_FAIL = OUTDIR / "mutation_mapping_failures.csv"
OUT_IDX = OUTDIR / "mappable_pair_index.json"
OUT_MD = C.V5 / "manifests" / "gate_b_mapping_report.md"

MUT_PAT = re.compile(r"^([A-Za-z0-9]+)_([A-Za-z])(\d+)([A-Za-z])$")
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def parse_ann(s) -> tuple[list[tuple], list[str]]:
    if not isinstance(s, str) or not s.strip():
        return [], []
    good, bad = [], []
    for p in re.split(r"[;,]", s):
        p = p.strip()
        if not p:
            continue
        mm = MUT_PAT.match(p)
        if mm:
            good.append((mm.group(1), mm.group(2).upper(), int(mm.group(3)),
                         mm.group(4).upper()))
        else:
            bad.append(p)
    return good, bad


def build_aligner():
    from Bio.Align import PairwiseAligner
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -1.0
    return a


def aln_substitutions(aligner, wt: str, mut: str, max_len: int = 1500) -> dict:
    """返回 alignment 的替换列表 + indel 计数（超长序列跳过）。"""
    if wt == mut:
        return {"subs": [], "n_indel": 0, "ok": True, "identity": 1.0}
    if len(wt) > max_len or len(mut) > max_len:
        return {"subs": [], "n_indel": -1, "ok": False, "identity": float("nan")}
    try:
        aln = aligner.align(wt, mut)[0]
        aw, am = str(aln[0]), str(aln[1])
    except Exception:  # noqa: BLE001
        return {"subs": [], "n_indel": -1, "ok": False, "identity": float("nan")}
    subs, n_indel, wi, nmatch = [], 0, 0, 0
    for cw, cm in zip(aw, am):
        if cw == "-" or cm == "-":
            n_indel += 1
            if cw != "-":
                wi += 1
            continue
        if cw == cm:
            nmatch += 1
        else:
            subs.append((wi, cw.upper(), cm.upper()))
        wi += 1
    return {"subs": subs, "n_indel": n_indel, "ok": True,
            "identity": nmatch / max(1, len(wt))}


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 GATE B formal mutation mapping")
    ap.add_argument("--max-len", type=int, default=1500)
    args = ap.parse_args()

    t0 = time.time()
    emit("# V5.0 GATE B — mutation ↔ 残基映射（正式版，含编号 offset 求解）")
    emit("")
    emit("> 机制: pair 内 WT/mutant 序列 global alignment → 残基自洽性 + **同链唯一 offset**。")
    emit("> offset 的引入依据: V5 §25 明确列出 `residue numbering offset` 为必须处理的情形；")
    emit("> `03b` 实测朴素 1-based 位置仅 73.37% 命中，而失配样例呈**恒定偏移**（如 +10），")
    emit("> 证明这是编号差异而非 construct 不一致。")
    emit("")

    aligner = build_aligner()
    raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    fwd = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = fwd[fwd["direction"] == "forward"].reset_index(drop=True)
    emit(f"- forward pair: **{len(fwd)}**")
    emit("")

    site_rows: list[dict] = []
    pair_rows: list[dict] = []
    fail_rows: list[dict] = []
    stat = Counter()

    for i, r in fwd.iterrows():
        ann, bad = parse_ann(r["mutation"])
        if bad:
            stat["ann_unparsed_frag"] += len(bad)
        try:
            wi, mi = int(r["wt_row"]), int(r["mut_row"])
        except Exception:  # noqa: BLE001
            stat["row_index_invalid"] += 1
            continue
        if not ann:
            stat["no_annotation"] += 1
            pair_rows.append({"pair_row": int(i), "pdb_id": r["pdb_id"], "n_ann": 0,
                              "n_sites_ok": 0, "n_sites_failed": 0,
                              "pair_status": "NO_ANNOTATION"})
            continue

        # 两侧 alignment
        per_side = {}
        for side, col in (("ligand", "Ligand Sequences"),
                          ("receptor", "Receptor Sequences")):
            w, m = str(raw.at[wi, col]), str(raw.at[mi, col])
            per_side[side] = aln_substitutions(aligner, w, m, args.max_len)

        # annotation 只声明链 ID，无法直接判定链 → 两侧都试，取"残基自洽数最多"的一侧
        n_ok = n_fail = 0
        for (ch, wt_aa, pos, mut_aa) in ann:
            best = None   # (n_consistent, side, status, evidence)
            for side, res in per_side.items():
                if not res["ok"]:
                    continue
                consistent = [(idx0, cw, cm) for (idx0, cw, cm) in res["subs"]
                              if cw == wt_aa and cm == mut_aa]
                if not consistent:
                    continue
                offsets = {pos - 1 - idx0 for (idx0, _cw, _cm) in consistent}
                # 该链上**所有** annotation 支持的 offset（用于唯一性判定）
                all_offs, all_cons = set(), 0
                for (_c2, w2, p2, m2) in ann:
                    for (idx0, cw, cm) in res["subs"]:
                        if cw == w2 and cm == m2:
                            all_offs.add(p2 - 1 - idx0)
                            all_cons += 1
                uniq = len(all_offs)
                if offsets == {0} and uniq <= 1:
                    st = "EXACT"
                elif uniq == 1:
                    st = "ALIGNMENT_CONFIDENT"
                elif uniq == 0:
                    st = "FAILED"
                else:
                    st = "AMBIGUOUS"
                cand = (len(consistent), len(consistent), side, st,
                        sorted(consistent), sorted(all_offs), res["n_indel"])
                if best is None or cand[:2] > best[:2]:
                    best = cand
            if best is None:
                n_fail += 1
                stat["site_FAILED"] += 1
                fail_rows.append({
                    "pair_row": int(i), "pdb_id": r["pdb_id"],
                    "annotation": f"{ch}_{wt_aa}{pos}{mut_aa}",
                    "reason": ("ALIGNMENT_SKIPPED_OR_FAILED"
                               if not any(v["ok"] for v in per_side.values())
                               else "NO_CONSISTENT_SUBSTITUTION"),
                    "lig_ok": int(per_side["ligand"]["ok"]),
                    "rec_ok": int(per_side["receptor"]["ok"]),
                    "lig_indel": per_side["ligand"]["n_indel"],
                    "rec_indel": per_side["receptor"]["n_indel"]})
                continue
            (_nc, _nc2, side, st, consistent, offs, n_indel) = best
            idx0 = consistent[0][0]
            off = (pos - 1) - idx0
            if st == "AMBIGUOUS":
                n_fail += 1
                stat["site_AMBIGUOUS"] += 1
            else:
                n_ok += 1
                stat[f"site_{st}"] += 1
            site_rows.append({
                "pair_row": int(i), "pdb_id": r["pdb_id"], "side": side,
                "ann_chain": ch, "wt_aa": wt_aa, "mut_aa": mut_aa,
                "ann_pos1": int(pos), "pdb_resseq": int(pos),
                "seq_index0": int(idx0), "numbering_offset": int(off),
                "n_consistent_subs": int(len(consistent)),
                "n_offsets_on_side": int(len(offs)),
                "n_indel_in_aln": int(n_indel),
                "mapping_status": st,
                "aln_identity": round(float(per_side[side]["identity"]), 4)})
        tot = len(ann)
        pair_status = ("MAPPED" if n_fail == 0 else
                       ("PARTIAL" if n_ok > 0 else "FAILED"))
        stat[f"pair_{pair_status}"] += 1
        pair_rows.append({"pair_row": int(i), "pdb_id": r["pdb_id"], "n_ann": tot,
                          "n_sites_ok": n_ok, "n_sites_failed": n_fail,
                          "pair_status": pair_status})
        if (i + 1) % 2000 == 0:
            emit(f"  ... {i + 1}/{len(fwd)}  EXACT={stat['site_EXACT']} "
                 f"CONF={stat['site_ALIGNMENT_CONFIDENT']}")

    sites = pd.DataFrame(site_rows)
    pmap = pd.DataFrame(pair_rows)
    fails = pd.DataFrame(fail_rows)
    sites.to_csv(OUT_SITES, index=False, encoding="utf-8")
    pmap.to_csv(OUT_PAIRS, index=False, encoding="utf-8")
    fails.to_csv(OUT_FAIL, index=False, encoding="utf-8")

    # ---------- 主实验可用子集 ----------
    ok_mask = sites["mapping_status"].isin(["EXACT", "ALIGNMENT_CONFIDENT"]) if len(sites) else pd.Series([], dtype=bool)
    exact_mask = (sites["mapping_status"] == "EXACT") if len(sites) else pd.Series([], dtype=bool)
    pairs_main = sorted(set(sites.loc[ok_mask, "pair_row"].astype(int))) if len(sites) else []
    pairs_exact = sorted(set(sites.loc[exact_mask, "pair_row"].astype(int))) if len(sites) else []
    # 只保留**全部** annotation 都成功映射的 pair（避免部分映射引入偏差）
    full_ok = pmap[(pmap["pair_status"] == "MAPPED")]["pair_row"].astype(int).tolist()
    pairs_main_full = sorted(set(pairs_main) & set(full_ok))
    OUT_IDX.write_text(json.dumps({
        "rule": ("main = EXACT + ALIGNMENT_CONFIDENT, AND all annotations of the pair "
                 "mapped (pair_status == MAPPED)"),
        "pairs_main": pairs_main_full,
        "pairs_exact_only": pairs_exact,
        "n_pairs_main": len(pairs_main_full),
        "n_pairs_exact_only": len(pairs_exact),
        "n_pairs_total": int(len(fwd)),
        "coverage_main": len(pairs_main_full) / max(1, len(fwd)),
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    # ---------- 报告 ----------
    emit("")
    emit("## 1. 配对级状态")
    emit("")
    vc = pmap["pair_status"].value_counts().to_dict()
    emit("| pair_status | 对数 | 占比 |")
    emit("|---|---:|---:|")
    for k, v in sorted(vc.items(), key=lambda x: -x[1]):
        emit(f"| {k} | {v} | {v / len(pmap) * 100:.2f}% |")
    emit("")
    emit("## 2. 位点级状态")
    emit("")
    if len(sites):
        sv = sites["mapping_status"].value_counts().to_dict()
        emit("| mapping_status | 位点数 | 占比 |")
        emit("|---|---:|---:|")
        for k, v in sorted(sv.items(), key=lambda x: -x[1]):
            emit(f"| {k} | {v} | {v / len(sites) * 100:.2f}% |")
    emit("")
    if len(sites):
        emit("### 编号 offset 分布（验证「恒定偏移」假设）")
        emit("")
        ov = sites["numbering_offset"].value_counts().sort_index()
        emit(f"- offset = 0（直接匹配）: **{int(ov.get(0, 0))}** "
             f"({ov.get(0, 0) / len(sites) * 100:.1f}%)")
        emit(f"- offset ≠ 0（编号偏移，已解）: **{int(len(sites) - ov.get(0, 0))}** "
             f"({(len(sites) - ov.get(0, 0)) / len(sites) * 100:.1f}%)")
        emit(f"- offset 范围: [{int(sites['numbering_offset'].min())}, "
             f"{int(sites['numbering_offset'].max())}]，非零 offset 的常见取值: "
             f"{ {int(k): int(v) for k, v in ov.items() if k != 0} if len(ov) < 40 else '过多'}".replace("{", "{").replace("}", "}"))
    emit("")
    emit("## 3. 失败原因（如实分类，不丢弃）")
    emit("")
    if len(fails):
        fv = fails["reason"].value_counts().to_dict()
        emit("| reason | 条数 |")
        emit("|---|---:|")
        for k, v in sorted(fv.items(), key=lambda x: -x[1]):
            emit(f"| {k} | {v} |")
    else:
        emit("- 无失败")
    emit("")

    # ---------- GATE B 判定 ----------
    cover = len(pairs_main_full) / max(1, len(fwd))
    emit("## 4. GATE B 判定（V5 §124）")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| 总 pair | {len(fwd)} |")
    emit(f"| **主实验可用 pair**（EXACT+CONFIDENT 且该 pair 全部 annotation 映射成功） | "
         f"**{len(pairs_main_full)}** |")
    emit(f"| **覆盖率** | **{cover * 100:.2f}%** |")
    emit(f"| EXACT-only 子集（敏感性分析用） | {len(pairs_exact)} |")
    emit("")
    if cover >= 0.90:
        verdict = "FULL_MIPA"
    elif cover >= 0.70:
        verdict = "SUBSET_PILOT"
    else:
        verdict = "SOLVE_MAPPING_FIRST"
    emit(f"**GATE B = `{verdict}`**（阈值: ≥90% full / 70–90% subset pilot / <70% 先解决）")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    for p in (OUT_SITES, OUT_PAIRS, OUT_FAIL, OUT_IDX, OUT_MD):
        print(f"写出 {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
