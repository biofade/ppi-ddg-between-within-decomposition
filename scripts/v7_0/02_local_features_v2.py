# -*- coding: utf-8 -*-
"""scripts/v7_0/02_local_features_v2.py -- V7.0 PHASE 1b: 修正后的 BLOCK L.

映射约定修正（关键）
--------------------
V7.0 PHASE 1 曾用 `(ann_chain, pdb_resseq)` 定位残基 → **WT 一致率仅 46.8%**，
误判为「映射不可靠」。PHASE 1b 系统排查四种约定后确定**正确约定**：

```
链 = ann_chain
链内第 seq_index0 个残基（0-based）
```

实测 **WT 一致率 93.1%**（EXACT 档，5123/5503），
而 `pdb_resseq` / `ann_pos1` / `numbering_offset` 三种写法均只有 46.9%
（≈ 随机水平）→ **修正的是本项目的调用方式，不是上游数据**。

本脚本据此重建 BLOCK L（局部结构环境），并只使用 **WT 校验通过**的位点。
"""

from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec)
sys.modules["v7f"] = v7f
_spec.loader.exec_module(v7f)

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH1 = C.ROOT / "experiments" / "v7_0" / "phase1_features"
PH1.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
T2O = v7f.THREE2ONE
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def chain_residues(path: Path):
    """{chain: [(resseq, icode, one_letter), ...]}（按 PDB 文件顺序）。"""
    ch = defaultdict(list)
    for ln in open(path, "r", errors="replace"):
        if not ln.startswith("ATOM") or ln[12:16].strip() != "CA":
            continue
        rn = ln[17:20].strip().upper()
        if rn not in T2O:
            continue
        try:
            seq = int(ln[22:26])
        except ValueError:
            continue
        ch[ln[21]].append((seq, ln[26].strip(), T2O[rn]))
    return dict(ch)


def main() -> int:
    t0 = time.time()
    emit("# V7.0 PHASE 1b — 修正后的 BLOCK L（局部结构环境，93.1% WT 校验）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **正确映射约定**：`链 = ann_chain`，`链内第 seq_index0 个残基`（0-based）")
    emit("- 实测 WT 一致率 **93.1%**（EXACT 档）；对照 `pdb_resseq`/`ann_pos1`/"
         "`numbering_offset` 均 46.9%")
    emit("- 只使用 **WT 校验通过** 的位点")
    emit("")

    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    s = site[site.pair_row.isin(can.row.values)].copy()
    emit(f"- canonical rows: **{len(can)}**；有映射记录: "
         f"**{s.pair_row.nunique()}**")
    emit(f"- mapping_status: "
         f"{s.mapping_status.value_counts().to_dict()}")
    emit("")

    cache: dict = {}
    def get_chains(pdb):
        if pdb.upper() not in cache:
            p = PDBDIR / f"{pdb.upper()}.pdb"
            cache[pdb.upper()] = chain_residues(p) if p.exists() else {}
        return cache[pdb.upper()]

    recs = []
    stat = {"n_total": 0, "n_nopdb": 0, "n_chain_missing": 0,
            "n_idx_oob": 0, "n_wt_mismatch": 0, "n_ok": 0}
    # 按 pair_row 聚合（一个 variant 可能有多条记录）
    for pr, g in s.groupby("pair_row"):
        stat["n_total"] += 1
        pdb = str(g.pdb_id.iloc[0])
        ch = get_chains(pdb)
        if not ch:
            stat["n_nopdb"] += 1
            continue
        locs = []
        for r in g.itertuples():
            c = str(r.ann_chain)
            if c not in ch:
                stat["n_chain_missing"] += 1
                continue
            ent = ch[c]
            i = int(r.seq_index0)
            if not (0 <= i < len(ent)):
                stat["n_idx_oob"] += 1
                continue
            resseq, icode, aa = ent[i]
            if aa != str(r.wt_aa):
                stat["n_wt_mismatch"] += 1
                continue
            locs.append({"chain": c, "resseq": resseq, "icode": icode,
                         "wt": str(r.wt_aa), "mut": str(r.mut_aa),
                         "side": str(r.side)})
        if not locs:
            continue
        lf = v7f.build_local_features(
            pdb, [{"pdb_chain": x["chain"], "pdb_resseq": x["resseq"],
                   "icode": x["icode"]} for x in locs], cache) if False else None
        # 直接调用（build_local_features 内部自行解析 PDB 并缓存）
        got = v7f.build_local_features(
            pdb, [{"pdb_chain": x["chain"], "pdb_resseq": x["resseq"],
                   "icode": x["icode"]} for x in locs],
            cache.setdefault("_v7f_cache", {}))
        if got is None:
            continue
        ok = [h for h in got if h is not None]
        if not ok:
            continue
        agg = {}
        for key in set(ok[0].keys()):
            vals = [h[key] for h in ok if key in h]
            agg[key] = float(np.mean(vals))
            agg[f"{key}__max"] = float(np.max(vals))
        idv = [v7f.mut_identity(x["wt"], x["mut"]) for x in locs
               if x["wt"] in v7f.AA1 and x["mut"] in v7f.AA1]
        if idv:
            for key in idv[0]:
                agg[key] = float(np.mean([d0[key] for d0 in idv]))
        agg["pair_row"] = int(pr)
        agg["pdb_id"] = pdb
        agg["n_sites_verified"] = len(locs)
        agg["n_sites_total"] = len(g)
        recs.append(agg)
        stat["n_ok"] += 1

    emit("## 映射统计（逐 pair_row）")
    emit("")
    emit("| 项 | 数 |")
    emit("|---|---:|")
    for k, v in stat.items():
        emit(f"| {k} | {v} |")
    emit("")
    LF = pd.DataFrame(recs)
    emit(f"### BLOCK L 覆盖")
    emit(f"- **rows: {len(LF)} / {len(can)} = {len(LF) / len(can) * 100:.1f}%**")
    emit(f"- **targets: {LF.pdb_id.nunique()} / {can.pdb_id.nunique()}**")
    emit(f"- 特征列（除元数据）: **{LF.shape[1] - 3}**")
    emit("")
    LF.to_csv(PH1 / "block_L_local.csv", index=False, encoding="utf-8")

    # 只保留 BLOCK L 可用的 canonical 子集索引（供 PHASE 2 建立统一 row universe）
    pd.DataFrame({"pair_row": LF.pair_row.astype(int),
                  "pdb_id": LF.pdb_id}).to_csv(
        PH1 / "matched_coverage_rows.csv", index=False, encoding="utf-8")
    emit(f"- 写出 `matched_coverage_rows.csv`（PHASE 2 的统一 row universe）")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "local_feature_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
    print(f"\n写出 {PH1 / 'block_L_local.csv'}  ({len(LF)} rows, "
          f"{LF.pdb_id.nunique()} targets)")
    print(f"写出 {PH1 / 'matched_coverage_rows.csv'}")
    print(f"写出 {PH1 / 'local_feature_report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
