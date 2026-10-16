# -*- coding: utf-8 -*-
"""scripts/v5/03_mutation_mapping.py -- V5.0 GATE B: mutation ↔ 残基映射.

目的（V5 §22–§26）
------------------
V4 已有 WT×mutant pair。因此**不需要**外部 WT 序列：
pair 内的 **WT row 序列** 与 **mutant row 序列** 可直接做 global alignment
（Biopython `PairwiseAligner`，V5 §23 明确要求不要手写 alignment）。

对每条 pair 的每条链（配体/受体）：
1. 取 WT row 与 mutant row 的序列；若完全相同 → 该链无突变（`NO_CHANGE`）。
2. 否则做 global alignment，定位差异位点。
3. 用 `Mutations` annotation（格式 `CHAIN_ORIGppMUT`，如 `A_R157E`；分隔符 `;`/`,`）
   独立核对：annotation 声明的链/原残基/位置/新残基是否与 alignment 一致。
4. 分类 `mapping_status` ∈ {EXACT, ALIGNMENT_CONFIDENT, AMBIGUOUS, FAILED}。

**主实验只允许 EXACT + ALIGNMENT_CONFIDENT**（V5 §26），并做 EXACT-only 敏感性分析。
失败条目写入 `mutation_mapping_failures.csv`，**不静默丢弃**（V5 §24）。

输出
----
* `data/v5_0/mutation_sites.csv`     每条 pair 的每个突变位点（含映射证据）
* `data/v5_0/mutation_pair_map.csv`  每条 pair 的汇总（覆盖/状态）
* `data/v5_0/mutation_mapping_failures.csv`
* `experiments/v5_0/manifests/gate_b_mapping_report.md`
"""

from __future__ import annotations

import argparse
import csv
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

OUTDIR = C.DATA / "v5_0"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_SITES = OUTDIR / "mutation_sites.csv"
OUT_PAIRS = OUTDIR / "mutation_pair_map.csv"
OUT_FAIL = OUTDIR / "mutation_mapping_failures.csv"
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


def parse_muts(s) -> tuple[list[tuple], list[str]]:
    """解析 Mutations → ([(chain, wt, pos, mut)], [无法解析的原始片段])。"""
    if not isinstance(s, str) or not s.strip():
        return [], []
    good, bad = [], []
    for p in re.split(r"[;,]", s):
        p = p.strip()
        if not p:
            continue
        mm = MUT_PAT.match(p)
        if mm:
            good.append((mm.group(1), mm.group(2), int(mm.group(3)), mm.group(4)))
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


def align_diffs(aligner, wt_seq: str, mut_seq: str) -> dict:
    """global alignment → 差异位点列表（含 gap 情况）。

    返回 {"subs": [(wt_idx0, wt_aa, mut_aa)], "n_indel": int, "identity": float,
          "aln_ok": bool}
    """
    if wt_seq == mut_seq:
        return {"subs": [], "n_indel": 0, "identity": 1.0, "aln_ok": True}
    try:
        alns = aligner.align(wt_seq, mut_seq)
        aln = alns[0]
    except Exception:  # noqa: BLE001
        return {"subs": [], "n_indel": -1, "identity": float("nan"), "aln_ok": False}
    a_wt, a_mut = str(aln[0]), str(aln[1])
    subs, n_indel = [], 0
    wi = 0
    n_match = 0
    for cw, cm in zip(a_wt, a_mut):
        if cw == "-" or cm == "-":
            n_indel += 1
            if cw != "-":
                wi += 1
            continue
        if cw == cm:
            n_match += 1
        else:
            subs.append((wi, cw, cm))
        wi += 1
    return {"subs": subs, "n_indel": n_indel,
            "identity": n_match / max(1, len(a_wt.replace("-", ""))),
            "aln_ok": True}


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 GATE B: mutation mapping")
    ap.add_argument("--limit", type=int, default=0, help=">0 时只处理前 N 条（pilot）")
    args = ap.parse_args()

    t0 = time.time()
    emit("# V5.0 GATE B — mutation ↔ 残基映射")
    emit("")
    emit("> 机制: pair 内 **WT row 序列 vs mutant row 序列** 做 Biopython global")
    emit("> alignment，再用 `Mutations` annotation **独立核对**（V5 §22–§24）。")
    emit("> **不需要外部 WT 序列字段** —— 这正是 V4 的 pair 结构提供的。")
    emit("")

    pairs = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = pairs[pairs["direction"] == "forward"].reset_index(drop=True)
    if args.limit:
        fwd = fwd.head(args.limit).reset_index(drop=True)
    emit(f"- forward pair: **{len(fwd)}**（来自 `data/v4.0/ddg_pairs.csv`）")
    emit("")

    aligner = build_aligner()
    site_rows: list[dict] = []
    pair_rows: list[dict] = []
    fail_rows: list[dict] = []

    stat = Counter()
    for i, r in fwd.iterrows():
        ann, bad = parse_muts(r["Mutations"])
        if bad:
            stat["annotation_unparsed_fragment"] += len(bad)
        if not ann:
            stat["no_annotation"] += 1
            pair_rows.append({"pair_row": int(i), "pdb_id": r["pdb_id"],
                              "n_ann": 0, "n_sites_mapped": 0, "n_sites_failed": 0,
                              "pair_status": "NO_ANNOTATION"})
            continue
        # 按链聚合 annotation
        by_chain: dict[str, list[tuple]] = {}
        for (ch, wt, pos, mu) in ann:
            by_chain.setdefault(ch, []).append((wt, pos, mu))

        lig_changed = (r["Ligand Sequences"] != r["wt_lig_seq"]) if "wt_lig_seq" in r else None
        n_ok = n_fail = 0
        for side, wtcol, mutcol, key_wt, key_mut in (
                ("ligand", "wt_lig_seq", "mut_lig_seq", "wt_lig_seq", "mut_lig_seq"),
                ("receptor", "wt_rec_seq", "mut_rec_seq", "wt_rec_seq", "mut_rec_seq")):
            wt_seq = str(r.get(key_wt, "")) if key_wt in r else ""
            mut_seq = str(r.get(key_mut, "")) if key_mut in r else ""
            if not wt_seq or not mut_seq:
                continue
            res = align_diffs(aligner, wt_seq, mut_seq)
            if not res["aln_ok"]:
                stat["alignment_failed"] += 1
                n_fail += len(by_chain.get(side, []))
                fail_rows.append({"pair_row": int(i), "pdb_id": r["pdb_id"],
                                  "side": side, "reason": "ALIGNMENT_FAILED",
                                  "annotation": r["Mutations"]})
                continue
            submap = {wi: (cw, cm) for wi, cw, cm in res["subs"]}
            # 与该链的 annotation 交叉核对
            for (wt_aa, pos1, mut_aa) in by_chain.get(side, []):
                idx0 = pos1 - 1                       # annotation 通常是 1-based
                match = None
                if idx0 in submap:
                    cw, cm = submap[idx0]
                    match = (cw == wt_aa.upper() and cm == mut_aa.upper())
                status = ("EXACT" if match else
                          ("AMBIGUOUS" if idx0 in submap else "FAILED"))
                if match:
                    n_ok += 1
                    stat["site_EXACT"] += 1
                else:
                    n_fail += 1
                    stat[f"site_{status}"] += 1
                    fail_rows.append({
                        "pair_row": int(i), "pdb_id": r["pdb_id"], "side": side,
                        "reason": "ANNOTATION_SEQUENCE_MISMATCH" if idx0 in submap
                                  else "POSITION_NOT_DIFFERENT",
                        "annotation": f"{wt_aa}{pos1}{mut_aa}",
                        "aln_wt_aa": submap.get(idx0, ("?", "?"))[0],
                        "aln_mut_aa": submap.get(idx0, ("?", "?"))[1],
                        "n_indel_in_aln": res["n_indel"],
                        "identity": round(float(res["identity"]), 4)})
                site_rows.append({
                    "pair_row": int(i), "pdb_id": r["pdb_id"], "side": side,
                    "chain_ann": r["Mutations"][:80] if isinstance(r["Mutations"], str) else "",
                    "seq_index0": idx0, "annotation_pos1": int(pos1),
                    "wt_aa_ann": wt_aa, "mut_aa_ann": mut_aa,
                    "aln_wt_aa": submap.get(idx0, ("?", "?"))[0],
                    "aln_mut_aa": submap.get(idx0, ("?", "?"))[1],
                    "mapping_status": status,
                    "n_indel_in_aln": int(res["n_indel"]),
                    "alignment_identity": round(float(res["identity"]), 4)})
        tot = len(by_chain.get("ligand", [])) + len(by_chain.get("receptor", []))
        pair_status = ("MAPPED" if n_fail == 0 and n_ok > 0 else
                       ("PARTIAL" if n_ok > 0 else "FAILED"))
        stat[f"pair_{pair_status}"] += 1
        pair_rows.append({"pair_row": int(i), "pdb_id": r["pdb_id"],
                          "n_ann": tot, "n_sites_mapped": n_ok,
                          "n_sites_failed": n_fail, "pair_status": pair_status})
        if (i + 1) % 1000 == 0:
            emit(f"  ... 已处理 {i + 1}/{len(fwd)}（EXACT 位点 {stat['site_EXACT']}）")

    sites = pd.DataFrame(site_rows)
    pmap = pd.DataFrame(pair_rows)
    fails = pd.DataFrame(fail_rows)
    sites.to_csv(OUT_SITES, index=False, encoding="utf-8")
    pmap.to_csv(OUT_PAIRS, index=False, encoding="utf-8")
    fails.to_csv(OUT_FAIL, index=False, encoding="utf-8")

    emit("")
    emit("## 1. 配对级覆盖")
    emit("")
    vc = pmap["pair_status"].value_counts().to_dict()
    emit("| pair_status | 对数 | 占比 |")
    emit("|---|---:|---:|")
    for k, v in sorted(vc.items(), key=lambda x: -x[1]):
        emit(f"| {k} | {v} | {v / len(pmap) * 100:.2f}% |")
    mapped = int((pmap["pair_status"] == "MAPPED").sum())
    partial = int((pmap["pair_status"] == "PARTIAL").sum())
    coverage = (mapped + partial) / max(1, len(pmap))
    emit("")
    emit(f"- **可用 pair（MAPPED + PARTIAL）= {mapped + partial} / {len(pmap)} = "
         f"{coverage * 100:.2f}%**")
    emit(f"- 完全映射（MAPPED）= {mapped}（{mapped / max(1, len(pmap)) * 100:.2f}%）")
    emit("")

    emit("## 2. 位点级状态")
    emit("")
    if len(sites):
        sv = sites["mapping_status"].value_counts().to_dict()
        emit("| mapping_status | 位点数 |")
        emit("|---|---:|")
        for k, v in sorted(sv.items(), key=lambda x: -x[1]):
            emit(f"| {k} | {v} |")
    else:
        emit("- 无位点")
    emit("")

    emit("## 3. 失败原因（如实分类，不丢弃）")
    emit("")
    if len(fails):
        fv = fails["reason"].value_counts().to_dict()
        emit("| reason | 条数 |")
        emit("|---|---:|")
        for k, v in sorted(fv.items(), key=lambda x: -x[1]):
            emit(f"| {k} | {v} |")
        if "n_indel_in_aln" in fails.columns:
            sub = fails[fails["n_indel_in_aln"] > 0]
            emit(f"- 其中 alignment 含 indel 的失败条目: **{len(sub)}**")
    else:
        emit("- 无失败")
    emit("")

    # ---------- GATE B 判定 ----------
    emit("## 4. GATE B 判定（V5 §124）")
    emit("")
    emit("| 覆盖率 | 动作（V5 §124 规定） |")
    emit("|---|---|")
    emit("| ≥90% | full MIPA |")
    emit("| 70–90% | subset pilot + matched baseline |")
    emit("| <70% | 先解决 mapping |")
    emit("")
    if coverage >= 0.90:
        verdict = "FULL_MIPA"
    elif coverage >= 0.70:
        verdict = "SUBSET_PILOT"
    else:
        verdict = "SOLVE_MAPPING_FIRST"
    emit(f"**覆盖率 = {coverage * 100:.2f}% → GATE B = `{verdict}`**")
    emit("")
    emit("> ⚠️ 注意：主实验只允许 `EXACT` + `ALIGNMENT_CONFIDENT`（V5 §26）。"
         "本脚本目前只产出 `EXACT`（annotation 与 alignment 完全一致）；"
         "未标 `EXACT` 的位点不进入主实验。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_SITES}")
    print(f"写出 {OUT_PAIRS}")
    print(f"写出 {OUT_FAIL}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
