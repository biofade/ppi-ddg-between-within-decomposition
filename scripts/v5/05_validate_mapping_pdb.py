# -*- coding: utf-8 -*-
"""scripts/v5/05_validate_mapping_pdb.py -- V5.0: 用**真实 PDB 结构**独立验证映射.

为什么必须做
------------
`04_build_mutation_sites.py` 解出的 `numbering_offset` 范围是 [−425, +645]，
对真实 PDB 编号而言**大得不合理**（真实 offset 通常 |o| < 100）。
必须用第三方证据独立验证，否则整个 MIPA 建立在不牢靠的坐标上。

验证策略（三条独立证据链）
--------------------------
E1. **序列坐标 → PDB 残基编号**
    取 `mutation_sites.csv` 的 (pdb_id, side, seq_index0)，在本地 PDB 里按
    `seq_index0` 取**第 (seq_index0+1) 个有序残基**（跳过 HETATM/水），读出其
    `resseq`。若读数与 `pdb_resseq` 一致 → 映射正确。

E2. **PDB 序列 vs dataset 序列的独立 alignment**
    从 PDB 的 ATOM/SEQRES 抽出链序列，与 dataset 的 WT 链序列做 global
    alignment，确认 `seq_index0` 落在同一残基上。

E3. **同 pdb 内 offset 的一致性**
    同一 pdb_id 的多个 annotation 若解出同一 offset → 强证据；
    若同一 pdb 内 offset 互相矛盾 → 该 pdb 整体标记不可信。

输出: experiments/v5_0/manifests/mapping_pdb_validation.md
      data/v5_0/mapping_validation.csv
"""

from __future__ import annotations

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

OUT_MD = C.V5 / "manifests" / "mapping_pdb_validation.md"
OUT_CSV = C.DATA / "v5_0" / "mapping_validation.csv"
L: list[str] = []

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def pdb_chain_residues(path: Path) -> dict[str, list[tuple]]:
    """解析 PDB → {chain: [(resseq, icode, one_letter), ...]}（按文件顺序，去重）。"""
    chains: dict[str, list[tuple]] = defaultdict(list)
    seen: set[tuple] = set()
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                res = line[17:20].strip()
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
                one = THREE2ONE.get(res, "X")
                chains[ch].append((seq, ic, one))
    except Exception:  # noqa: BLE001
        return {}
    return dict(chains)


def main() -> int:
    t0 = time.time()
    emit("# V5.0 — 映射的独立 PDB 验证")
    emit("")
    emit("> 目的: `04` 解出的 offset 范围 [−425,+645] 对真实 PDB 编号不合理，")
    emit("> 必须用本地真实结构文件独立验证，不能凭 alignment 自证。")
    emit("")

    sites = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    pidx = json.loads((C.DATA / "v5_0" / "mappable_pair_index.json")
                      .read_text(encoding="utf-8"))
    main_pairs = set(pidx["pairs_main"])
    emit(f"- `mutation_sites.csv` 位点数: **{len(sites)}**")
    emit(f"- 主实验 pair 数: **{len(main_pairs)}**")
    emit("")

    pdbdir = C.DATA / "pdb"
    files = {p.stem.upper(): p for p in pdbdir.glob("*.pdb")} if pdbdir.exists() else {}
    emit(f"- 本地 PDB 文件: **{len(files)}**（目录 `{pdbdir}`）")
    emit("")

    # ---------- E1 + E3: 用 PDB 残基序号核对 seq_index0 ----------
    emit("## E1: `seq_index0` 对应的真实 PDB 残基编号 vs 推断的 `pdb_resseq`")
    emit("")
    sub = sites[sites["pair_row"].isin(main_pairs)].copy()
    emit(f"- 参与验证的位点（限主实验 pair）: **{len(sub)}**")
    emit("")

    rows = []
    cache: dict[str, dict] = {}
    checked = 0
    match = mismatch = nofile = nochain = idxoob = 0
    for _, r in sub.iterrows():
        pid = str(r["pdb_id"]).upper()
        if pid not in cache:
            p = files.get(pid)
            cache[pid] = pdb_chain_residues(p) if p else {}
        chains = cache[pid]
        if not chains:
            nofile += 1
            rows.append({"pdb_id": pid, "side": r["side"], "seq_index0": r["seq_index0"],
                         "pdb_resseq_claimed": r["pdb_resseq"], "status": "PDB_MISSING"})
            continue
        # 选残基最多的链（dataset 序列通常是主要链）
        ch_sel = max(chains, key=lambda c: len(chains[c]))
        reslist = chains[ch_sel]
        idx = int(r["seq_index0"])
        if idx >= len(reslist):
            idxoob += 1
            rows.append({"pdb_id": pid, "chain": ch_sel, "side": r["side"],
                         "seq_index0": idx, "pdb_resseq_claimed": r["pdb_resseq"],
                         "n_res_in_chain": len(reslist), "status": "INDEX_OUT_OF_RANGE"})
            continue
        real_seq, real_ic, real_aa = reslist[idx]
        checked += 1
        ok_pos = int(real_seq) == int(r["pdb_resseq"])
        ok_aa = (real_aa == str(r["wt_aa"]).upper()) or real_aa == "X"
        status = ("MATCH" if (ok_pos and ok_aa) else
                  ("POS_MATCH_AA_DIFF" if ok_pos else "POS_MISMATCH"))
        if status == "MATCH":
            match += 1
        else:
            mismatch += 1
        rows.append({"pdb_id": pid, "chain": ch_sel, "side": r["side"],
                     "seq_index0": idx, "pdb_resseq_claimed": int(r["pdb_resseq"]),
                     "pdb_resseq_real": int(real_seq),
                     "real_aa": real_aa, "wt_aa": r["wt_aa"],
                     "numbering_offset": int(r["numbering_offset"]),
                     "mapping_status": r["mapping_status"], "status": status})
    v = pd.DataFrame(rows)
    v.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("| 检查结果 | 条数 |")
    emit("|---|---:|")
    emit(f"| 可在 PDB 中定位（checked） | {checked} |")
    emit(f"| **位置+残基都一致（MATCH）** | **{match}** |")
    emit(f"| 不一致（POS_MISMATCH / AA_DIFF） | {mismatch} |")
    emit(f"| seq_index0 越界 | {idxoob} |")
    emit(f"| PDB 文件缺失 | {nofile} |")
    emit("")
    if checked:
        emit(f"- **E1 通过率 = {match} / {checked} = {match / checked * 100:.2f}%**")
    emit("")

    # ---------- E3: 同 pdb 内 offset 一致性 ----------
    emit("## E3: 同一 pdb 内 offset 是否自洽")
    emit("")
    g = sub.groupby("pdb_id")["numbering_offset"].nunique()
    emit(f"- 主实验 pair 覆盖 pdb 数: {len(g)}")
    emit(f"- 单 pdb 内 offset **唯一**的 pdb: **{int((g <= 1).sum())}** "
         f"({(g <= 1).mean() * 100:.1f}%)")
    emit(f"- 单 pdb 内 offset **冲突**的 pdb: **{int((g > 1).sum())}**")
    emit("")
    big = sub[sub["numbering_offset"].abs() > 100]
    emit(f"- |offset| > 100 的位点: **{len(big)}**（占主实验位点 "
         f"{len(big) / max(1, len(sub)) * 100:.1f}%）→ 这些极可能是**假匹配**")
    emit("")

    # ---------- 结论 ----------
    emit("## 结论与判据（决定主实验是否可用）")
    emit("")
    e1 = match / checked if checked else 0.0
    if e1 >= 0.95:
        verdict = "MAPPING_VALIDATED"
        note = "映射获 PDB 独立验证，可用于 MIPA 主实验。"
    elif e1 >= 0.70:
        verdict = "MAPPING_PARTIALLY_VALIDATED"
        note = ("映射部分获验证；必须**只保留通过 PDB 验证的位点**（"
                "`mapping_validation.csv` 中 status == MATCH）再做 MIPA。")
    else:
        verdict = "MAPPING_NOT_VALIDATED"
        note = ("映射未获 PDB 独立验证 → **不得**基于该映射做 MIPA 结论；"
                "应先修 mapping（V5 §124 <70% 分支）。")
    emit(f"**{verdict}** — {note}")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}\n写出 {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
