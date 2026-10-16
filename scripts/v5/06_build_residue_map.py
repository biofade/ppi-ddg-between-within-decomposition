# -*- coding: utf-8 -*-
"""scripts/v5/06_build_residue_map.py -- V5.0 GATE B（修订版）: 序列引导的 PDB 残基映射.

为什么重做（V5 §123 Bug protocol）
---------------------------------
`04`/`05` 的失败根因（已定位，PROJECT_EVIDENCE）:
  * `04` 只从 dataset 序列对齐解出「序列内下标 seq_index0」，**从未**确定
    dataset 序列 ↔ **具体 PDB 链**的对应关系；
  * `05` 用「取最长链」这个**错误假设**做验证，得到 E1 通过率仅 **7.96%**，
    且 99/226 个 pdb 的 offset 自相矛盾。
  ⇒ `04` 产出的 `numbering_offset` / `pdb_resseq` 列 **INVALID_DUE_TO_BUG**，
    不得用于任何科学结论。本脚本取代它。

修订后的正确流程（严格遵照 V5 §25「优先: PDB sequence ↔ dataset sequence global alignment」）
-------------------------------------------------------------------------------------------
对每个 pdb_id：
1. 读入 PDB，对**每条蛋白链**抽取有序残基（resseq, icode, 1-letter），
   同时处理 `SEQRES` 缺失情况（只用 ATOM 亦可，但需记录）。
2. 把 dataset 的 WT 序列（`Ligand Sequences` / `Receptor Sequences`）按 `,`
   拆成**子链**（多链拼接情形，V5 §23「chain-specific」）。
3. 对每个 (dataset 子链, PDB 链) 组合做 global alignment；
   选 **identity 最高且 ≥0.95** 的组合作为该 dataset 子链的映射。
4. 由 alignment 建立 **dataset 下标 ↔ PDB 残基**的**逐位置**映射
   （含 gap：dataset 位置插入 → PDB 无对应残基）。
5. 用该逐位置映射，把 `04` 解出的 `seq_index0` 转成真实 `pdb_resseq` / `icode`，
   并复核该位置的实际氨基酸 = annotation 的 WT 残基。
6. `mapping_status` 重新定义:
   * `EXACT`        — 链 identity = 1.0 且逐位置映射无 gap，残基一致
   * `ALIGNMENT_CONFIDENT` — 链 identity ≥0.95，残基一致（允许 gap）
   * `AMBIGUOUS`    — 多个链/子链组合并列最优
   * `FAILED`       — 无 ≥0.95 的组合，或残基不一致

输出（**取代 04 的对应列**）
---------------------------
* `data/v5_0/residue_map.csv`     逐位点: pair_row, side, dataset_subchain,
                                   pdb_chain, seq_index0, pdb_resseq, icode,
                                   wt_aa, real_aa, mapping_status, chain_identity
* `data/v5_0/residue_map_summary.json`
* `experiments/v5_0/manifests/gate_b_residue_map_report.md`
"""

from __future__ import annotations

import argparse
import json
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

OUT_CSV = C.DATA / "v5_0" / "residue_map.csv"
OUT_JSON = C.DATA / "v5_0" / "residue_map_summary.json"
OUT_MD = C.V5 / "manifests" / "gate_b_residue_map_report.md"

L: list[str] = []
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O"}
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def read_pdb_chains(path: Path) -> tuple[dict[str, list[tuple]], dict[str, str]]:
    """返回 ({chain: [(resseq, icode, aa)]}, {chain: seqres_string})。"""
    atom: dict[str, list[tuple]] = defaultdict(list)
    seen: set[tuple] = set()
    seqres: dict[str, list[str]] = defaultdict(list)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("SEQRES"):
                    ch = line[11]
                    for i in range(19, min(len(line.rstrip("\n")), 70), 3):
                        res = line[i:i + 3].strip()
                        if res:
                            seqres[ch].append(AA3.get(res, "X"))
                    continue
                if not line.startswith("ATOM"):
                    continue
                res = line[17:20].strip()
                aa = AA3.get(res)
                if aa is None:
                    continue
                ch = line[21]
                try:
                    seq = int(line[22:26])
                except ValueError:
                    continue
                ic = line[26].strip()
                key = (ch, seq, ic)
                if key in seen:
                    continue
                seen.add(key)
                atom[ch].append((seq, ic, aa))
    except Exception:  # noqa: BLE001
        return {}, {}
    return dict(atom), {k: "".join(v) for k, v in seqres.items()}


def build_aligner():
    from Bio.Align import PairwiseAligner
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -6.0
    a.extend_gap_score = -0.5
    return a


def align_posmap(aligner, ds_seq: str, pdb_seq: str) -> dict | None:
    """global alignment → {ds_idx0: (pdb_idx0, ds_aa, pdb_aa)} + metrics。

    pdb_idx0 是 pdb_seq 内的下标（0-based）；调用方再映射到 (resseq, icode)。
    """
    if not ds_seq or not pdb_seq:
        return None
    try:
        aln = aligner.align(ds_seq, pdb_seq)[0]
        a_ds, a_pd = str(aln[0]), str(aln[1])
    except Exception:  # noqa: BLE001
        return None
    posmap: dict[int, tuple] = {}
    di = pi = 0
    nmatch = ngap = 0
    for cd, cp in zip(a_ds, a_pd):
        if cd == "-":
            pi += 1
            continue
        if cp == "-":
            posmap[di] = (-1, cd, "-")
            di += 1
            ngap += 1
            continue
        posmap[di] = (pi, cd, cp)
        if cd == cp:
            nmatch += 1
        di += 1
        pi += 1
    aln_len = sum(1 for cd, cp in zip(a_ds, a_pd) if cd != "-")
    return {"posmap": posmap, "identity": nmatch / max(1, aln_len),
            "n_gap": ngap, "aln_len": aln_len}


def main() -> int:
    ap = argparse.ArgumentParser(description="V5 GATE B revised: residue map")
    ap.add_argument("--min-identity", type=float, default=0.95)
    ap.add_argument("--max-pairs", type=int, default=0)
    args = ap.parse_args()

    t0 = time.time()
    emit("# V5.0 GATE B（修订版）— 序列引导的 PDB 残基映射")
    emit("")
    emit("> **前置负结果**: `04`+`05` 的映射方案未通过独立验证（E1 通过率 7.96%），")
    emit("> 根因是**从未确定 dataset 序列 ↔ 具体 PDB 链**的对应（错取最长链）。")
    emit("> 本脚本按 V5 §25「PDB sequence ↔ dataset sequence global alignment」重做。")
    emit("")

    aligner = build_aligner()
    sites = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    pidx = json.loads((C.DATA / "v5_0" / "mappable_pair_index.json")
                      .read_text(encoding="utf-8"))
    main_pairs = set(int(x) for x in pidx["pairs_main"])
    sites = sites[sites["pair_row"].isin(main_pairs)].copy()
    emit(f"- 待重映射位点: **{len(sites)}**（来自 {len(main_pairs)} 个主实验 pair）")

    raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    fwd = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = fwd[fwd["direction"] == "forward"].reset_index(drop=True)
    pdbdir = C.DATA / "pdb"
    emit(f"- 本地 PDB: {len(list(pdbdir.glob('*.pdb')))} 个")
    emit("")

    # 预解析涉及的 PDB
    need_pdb = sorted(set(str(s).upper() for s in sites["pdb_id"]))
    emit(f"- 涉及 PDB: **{len(need_pdb)}** → 解析中 ...")
    pdbcache: dict[str, tuple] = {}
    for i, pid in enumerate(need_pdb):
        p = pdbdir / f"{pid}.pdb"
        if not p.exists():
            p = pdbdir / f"{pid.lower()}.pdb"
        pdbcache[pid] = read_pdb_chains(p) if p.exists() else ({}, {})
        if (i + 1) % 50 == 0:
            emit(f"   ... {i + 1}/{len(need_pdb)}")
    emit(f"- 解析完成；有 ATOM 蛋白链的 PDB: "
         f"**{sum(1 for a, _s in pdbcache.values() if a)} / {len(need_pdb)}**")
    emit("")

    # 预计算 dataset 序列（WT 行）的子链
    def subchains(seq: str) -> list[str]:
        return [s.strip() for s in str(seq).split(",") if s.strip()]

    rows: list[dict] = []
    stat = Counter()
    # 缓存 (pid, side, ds_sub_index) -> (chain, identity, posmap, chain_reslist)
    maps: dict[tuple, tuple | None] = {}

    for n, (_, r) in enumerate(sites.iterrows()):
        pid = str(r["pdb_id"]).upper()
        atom, seqres = pdbcache.get(pid, ({}, {}))
        if not atom:
            stat["pdb_no_protein_chain"] += 1
            rows.append({**r.to_dict(), "pdb_resseq_real": -1, "icode": "",
                         "real_aa": "?", "pdb_chain": "", "chain_identity": np.nan,
                         "validated_status": "PDB_MISSING"})
            continue
        pr = fwd.iloc[int(r["pair_row"])]
        try:
            wi = int(pr["wt_row"])
        except Exception:  # noqa: BLE001
            stat["wt_row_invalid"] += 1
            continue
        col = "Ligand Sequences" if r["side"] == "ligand" else "Receptor Sequences"
        subs = subchains(raw.at[wi, col])
        idx0 = int(r["seq_index0"])

        # 判断 idx0 落在哪个子链
        off = 0
        sub_i, local_idx = None, None
        for si, s in enumerate(subs):
            if off <= idx0 < off + len(s):
                sub_i, local_idx = si, idx0 - off
                break
            off += len(s)
        if sub_i is None:
            stat["idx_outside_any_subchain"] += 1
            rows.append({**r.to_dict(), "pdb_resseq_real": -1, "icode": "",
                         "real_aa": "?", "pdb_chain": "", "chain_identity": np.nan,
                         "validated_status": "SUBCHAIN_NOT_FOUND"})
            continue
        ds_seq = subs[sub_i]

        key = (pid, r["side"], sub_i)
        if key not in maps:
            best = None
            for ch, reslist in atom.items():
                pdb_seq = "".join(a for _s, _i, a in reslist)
                res = align_posmap(aligner, ds_seq, pdb_seq)
                if res is None:
                    continue
                if best is None or res["identity"] > best[2]["identity"]:
                    best = (ch, reslist, res)
            maps[key] = best
        best = maps[key]
        if best is None:
            stat["no_alignment"] += 1
            rows.append({**r.to_dict(), "pdb_resseq_real": -1, "icode": "",
                         "real_aa": "?", "pdb_chain": "", "chain_identity": np.nan,
                         "validated_status": "ALIGNMENT_FAILED"})
            continue
        ch, reslist, res = best
        ident = res["identity"]
        pm = res["posmap"].get(local_idx)
        if pm is None:
            stat["idx_not_in_aln"] += 1
            rows.append({**r.to_dict(), "pdb_resseq_real": -1, "icode": "",
                         "real_aa": "?", "pdb_chain": ch,
                         "chain_identity": round(float(ident), 4),
                         "validated_status": "INDEX_NOT_IN_ALIGNMENT"})
            continue
        pdb_i, ds_aa, pdb_aa = pm
        if pdb_i < 0 or pdb_aa == "-":
            stat["gap_at_position"] += 1
            rows.append({**r.to_dict(), "pdb_resseq_real": -1, "icode": "",
                         "real_aa": "-", "pdb_chain": ch,
                         "chain_identity": round(float(ident), 4),
                         "validated_status": "GAP_AT_MUTATION_SITE"})
            continue
        real_seq, real_ic, real_aa = reslist[pdb_i]
        aa_ok = (real_aa == str(r["wt_aa"]).upper())
        if ident >= 0.999 and aa_ok and res["n_gap"] == 0:
            st = "EXACT"
        elif ident >= args.min_identity and aa_ok:
            st = "ALIGNMENT_CONFIDENT"
        elif aa_ok:
            st = "LOW_IDENTITY_CHAIN"
        else:
            st = "AA_MISMATCH"
        stat[f"val_{st}"] += 1
        rows.append({**r.to_dict(), "pdb_resseq_real": int(real_seq),
                     "icode": real_ic, "real_aa": real_aa, "pdb_chain": ch,
                     "ds_subchain_index": sub_i,
                     "chain_identity": round(float(ident), 4),
                     "chain_n_gap": int(res["n_gap"]),
                     "validated_status": st})
        if (n + 1) % 2000 == 0:
            emit(f"  ... {n + 1}/{len(sites)}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False, encoding="utf-8")

    emit("")
    emit("## 1. 修订后逐位点状态")
    emit("")
    vc = out["validated_status"].value_counts().to_dict()
    tot = len(out)
    emit("| validated_status | 位点数 | 占比 |")
    emit("|---|---:|---:|")
    for k, v in sorted(vc.items(), key=lambda x: -x[1]):
        emit(f"| `{k}` | {v} | {v / max(1, tot) * 100:.2f}% |")
    emit("")

    ok = out[out["validated_status"].isin(["EXACT", "ALIGNMENT_CONFIDENT"])]
    emit("## 2. 链 identity 分布（验证「链对应」是否找对）")
    emit("")
    idt = out["chain_identity"].dropna()
    if len(idt):
        for lo, hi, name in ((0.999, 1.01, "=1.000"), (0.95, 0.999, "0.95–0.999"),
                             (0.80, 0.95, "0.80–0.95"), (0.0, 0.80, "<0.80")):
            m = (idt >= lo) & (idt < hi if hi <= 1.0 else idt <= 1.0)
            emit(f"- identity {name}: {int(m.sum())} ({m.mean() * 100:.2f}%)")
    emit("")
    emit("## 3. 独立验证：`pdb_resseq_real` 处真实残基是否 == annotation WT 残基")
    emit("")
    chk = out[out["real_aa"] != "?"]
    if len(chk):
        okk = int((chk["real_aa"] == chk["wt_aa"]).sum())
        emit(f"- 可检查位点: **{len(chk)}**；残基一致: **{okk}** "
             f"（**{okk / len(chk) * 100:.2f}%**）")
    emit("")
    emit("## 4. GATE B 判定（修订版）")
    emit("")
    n_ok = len(ok)
    cover = n_ok / max(1, tot)
    emit(f"- 可用位点（EXACT + ALIGNMENT_CONFIDENT）: **{n_ok} / {tot} = "
         f"{cover * 100:.2f}%**")
    emit(f"- 覆盖的主实验 pair: **{out[out['validated_status'].isin(['EXACT','ALIGNMENT_CONFIDENT'])]['pair_row'].nunique()}**")
    emit("")
    verdict = ("FULL_MIPA" if cover >= 0.90
               else ("SUBSET_PILOT" if cover >= 0.70 else "SOLVE_MAPPING_FIRST"))
    emit(f"**GATE B（修订版）= `{verdict}`**")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    OUT_JSON.write_text(json.dumps({
        "n_sites_input": int(tot), "n_sites_ok": int(n_ok),
        "site_coverage": cover,
        "verdict": verdict,
        "validated_status_counts": {k: int(v) for k, v in vc.items()},
        "min_identity": args.min_identity,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    for p in (OUT_CSV, OUT_JSON, OUT_MD):
        print(f"写出 {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
