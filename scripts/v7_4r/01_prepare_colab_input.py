# -*- coding: utf-8 -*-
"""scripts/v7_4r/01_prepare_colab_input.py -- V7.4R：准备**无标签** Colab 输入包（§7）.

从 `features/v7_4/u_relax_strict_manifest.csv`（U_RELAX_STRICT）生成
Colab 端所需的最小输入，并打成 zip。

输出
----
`colab/v7_4r/colab_input_package.zip`（**零标签**）
`colab/v7_4r/input_manifest.csv`
`experiments/v7_4r/phase1_prepare/prepare_report.md`
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
PH1 = R / "experiments" / "v7_4r" / "phase1_prepare"
PH1.mkdir(parents=True, exist_ok=True)
COLAB = R / "colab" / "v7_4r"
COLAB.mkdir(parents=True, exist_ok=True)
STAGE = COLAB / "staging"
STAGE.mkdir(parents=True, exist_ok=True)
PDBDIR = R / "data" / "pdb"
STD3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    sp = importlib.util.spec_from_file_location(
        "v7f", R / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    T2O = v7f.THREE2ONE

    emit("# V7.4R PHASE 1 — Colab 输入包（无标签）")
    emit("")

    man = pd.read_csv(R / "features" / "v7_4" /
                      "u_relax_strict_manifest.csv", low_memory=False)
    emit(f"- U_RELAX_STRICT：**{len(man)}** 行 / "
         f"**{man.target_id_canonical.nunique()}** targets / "
         f"**{man.family_id_canonical.nunique()}** families")

    # 解析每个 target 的链序列（只需突变涉及链）
    ch_cache = {}
    tc = []
    for t in sorted(man.target_id_canonical.astype(str).unique()):
        p = PDBDIR / f"{t.upper()}.pdb"
        if not p.exists():
            continue
        res = defaultdict(list)
        order = defaultdict(list)
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
            if not res[c] or res[c][-1] != key:
                res[c].append(key)
        ch_cache[t] = res
        for c in sorted(res):
            seq = "".join(T2O[a] for (_, _, a) in res[c])
            tc.append({"target_id": t, "chain": c, "wt_sequence": seq,
                       "length": len(seq)})
    TC = pd.DataFrame(tc)
    emit(f"- 解析 (target,chain)：**{len(TC)}** 条 / "
         f"**{TC.target_id.nunique()}** targets")

    # 位点表（seq_index0 由 manifest 的 mutation_bodies 无法反推 -> 用 mutation_sites + 链下标）
    site = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    sit_rows = []
    for r in man.itertuples():
        rid = int(r.row_id)
        g = s_by_row.get(rid)
        if g is None or not len(g):
            continue
        for k, x in enumerate(g.itertuples()):
            sit_rows.append({"row_id": rid,
                             "target_id": str(r.target_id_canonical),
                             "family_id": str(r.family_id_canonical),
                             "site_index": k, "chain": str(x.ann_chain),
                             "seq_index0": int(x.seq_index0),
                             "pdb_resseq": int(x.pdb_resseq),
                             "icode": str(getattr(x, "icode", "") or ""),
                             "wt_aa": str(x.wt_aa), "mut_aa": str(x.mut_aa)})
    SIT = pd.DataFrame(sit_rows)
    emit(f"- 位点行：**{len(SIT)}**")
    VAR = man[["row_id", "target_id_canonical", "family_id_canonical",
               "n_sites", "mutated_chains", "mutation_bodies"]].rename(
        columns={"target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})
    CANON = man[["row_id", "target_id_canonical", "family_id_canonical"]].rename(
        columns={"target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})

    outp = {"variants.csv": VAR, "mutation_sites.csv": SIT,
            "target_chains.csv": TC, "canonical_rows.csv": CANON}
    for n, df in outp.items():
        df.to_csv(STAGE / n, index=False, encoding="utf-8")
    (STAGE / "PACKAGE_README.md").write_text(
        "# V7.4R Colab Input Package（OpenMM 松弛，**无标签**）\n\n"
        "| file | 说明 |\n|---|---|\n"
        "| `canonical_rows.csv` | `row_id` / `target_id` / `family_id`（权威命名空间）|\n"
        "| `variants.csv` | 每变体位点数、涉及链、突变体串 |\n"
        "| `mutation_sites.csv` | 链、`seq_index0`、`pdb_resseq`、`icode`、WT、MUT |\n"
        "| `target_chains.csv` | 每 (target, chain) 的 WT 链序列 |\n\n"
        "结构：`target_id` 即 PDB ID，从 "
        "`https://files.rcsb.org/download/<ID>.pdb` 下载。\n\n"
        "**无 ΔΔG / KD / 任何标签。**\n", encoding="utf-8")

    man_rows = []
    for n in list(outp) + ["PACKAGE_README.md"]:
        p = STAGE / n
        man_rows.append({"file": n, "size_bytes": p.stat().st_size,
                         "sha256_prefix16": sha(p)[:16],
                         "contains_labels": "NO"})
    M = pd.DataFrame(man_rows)
    M.to_csv(COLAB / "input_manifest.csv", index=False, encoding="utf-8")
    zp = COLAB / "colab_input_package.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in list(outp) + ["PACKAGE_README.md"]:
            zf.write(STAGE / n, arcname=n)
        zf.write(COLAB / "input_manifest.csv", arcname="input_manifest.csv")

    LABTOK = ["ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred"]
    viol = [(n, c) for n, df in outp.items() for c in df.columns
            if any(t in c.lower() for t in LABTOK)]
    emit("")
    emit("| file | bytes | 含标签 |")
    emit("|---|---:|---|")
    for r in man_rows:
        emit(f"| `{r['file']}` | {r['size_bytes']:,} | **{r['contains_labels']}** |")
    emit("")
    emit(f"- ZIP：`colab_input_package.zip`（{zp.stat().st_size / 1e6:.2f} MB）")
    emit(f"- SHA256：`{sha(zp)}`")
    emit(f"- **无标签检查：{len(viol)} 处违规 {viol if viol else '→ ✅ 无'}**")
    (PH1 / "prepare_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    import json
    json.dump({"zip": str(zp), "sha256": sha(zp), "n_rows": len(man),
               "n_targets": int(man.target_id_canonical.nunique()),
               "n_families": int(man.family_id_canonical.nunique()),
               "n_chains": len(TC), "n_sites": len(SIT),
               "zip_bytes": zp.stat().st_size},
              open(PH1 / "package_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n>>> rows={len(man)} targets={man.target_id_canonical.nunique()} "
          f"chains={len(TC)} sites={len(SIT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
