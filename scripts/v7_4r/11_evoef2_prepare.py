# -*- coding: utf-8 -*-
"""scripts/v7_4r/11_evoef2_prepare.py -- V7.4R Track B：准备 EvoEF2 输入（§34–§37）.

从 V7.3 `U_PHYS_STRICT`（或 V7.4 `U_RELAX_STRICT`）生成 EvoEF2 所需：

1. **蛋白链专用 PDB**（去水/配体/离子/氢）——EvoEF2 只认氨基酸
2. **`individual_list.txt`**：每行一个变体，语法
   `<WT><chain><pdb_resseq><MUT>`（逗号分隔、无空格、以 `;` 结尾）
3. **变体索引 CSV**：把 EvoEF2 的输出文件名映射回 `row_id`

**不做**任何标签关联（零 ΔΔG）。
"""

from __future__ import annotations

import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "evoef2" / "v7_4r"
OUT.mkdir(parents=True, exist_ok=True)
STRUCT = OUT / "structures"
STRUCT.mkdir(parents=True, exist_ok=True)
PH3 = R / "experiments" / "v7_4r" / "phase3_evoef2"
PH3.mkdir(parents=True, exist_ok=True)
PDBDIR = R / "data" / "pdb"
STD3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))


def protein_only(src: Path, dst: Path) -> int:
    """写出去水/去配体/去氢的蛋白 PDB；返回保留原子数。"""
    n = 0
    with open(src, "r", errors="replace") as fi, open(dst, "w") as fo:
        for ln in fi:
            if ln[:6] != "ATOM  ":
                continue
            rn = ln[17:20].strip().upper()
            if rn not in STD3:
                continue
            el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if el == "H":
                continue
            fo.write(ln if ln.endswith("\n") else ln + "\n")
            n += 1
    return n


def main() -> int:
    emit("# V7.4R Track B — EvoEF2 输入准备")
    emit("")
    man = pd.read_csv(R / "features" / "v7_4" /
                      "u_relax_strict_manifest.csv", low_memory=False)
    site = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    emit(f"- 来源 U_RELAX_STRICT：**{len(man)}** 行 / "
         f"**{man.target_id_canonical.nunique()}** targets")

    rows = []
    indiv = []              # individual_list.txt 行
    n_pdb = 0
    n_missing = 0
    for r in man.itertuples():
        rid = int(r.row_id)
        t = str(r.target_id_canonical).upper()
        g = s_by_row.get(rid)
        if g is None or not len(g):
            n_missing += 1
            continue
        src = PDBDIR / f"{t}.pdb"
        if not src.exists():
            n_missing += 1
            continue
        dst = STRUCT / f"{t}_protein.pdb"
        if not dst.exists():
            protein_only(src, dst)
            n_pdb += 1
        muts = []
        ok = True
        for x in g.itertuples():
            wt = str(x.wt_aa); mu = str(x.mut_aa)
            ch = str(x.ann_chain); rs = int(x.pdb_resseq)
            if wt not in "ACDEFGHIKLMNPQRSTVWY" or mu not in "ACDEFGHIKLMNPQRSTVWY":
                ok = False; break
            if not ch or not ch.isalnum():
                ok = False; break
            muts.append(f"{wt}{ch}{rs}{mu}")
        if not ok or not muts:
            n_missing += 1
            continue
        token = ",".join(muts) + ";"
        indiv.append(token)
        rows.append({"row_id": rid, "target_id": t,
                     "family_id": str(r.family_id_canonical),
                     "n_sites": len(muts), "evoef2_mutation": token,
                     "protein_pdb": dst.name})
    E = pd.DataFrame(rows)
    (OUT / "individual_list.txt").write_text("\n".join(indiv) + "\n", encoding="utf-8")
    E.to_csv(OUT / "evoef2_variant_index.csv", index=False, encoding="utf-8")

    emit(f"- 写出蛋白链 PDB：**{n_pdb}** 个（去水/配体/离子/氢）")
    emit(f"- 变体数：**{len(E)}**；跳过 **{n_missing}**（位点/结构缺失）")
    emit(f"- 多点突变（n_sites>=2）：**{int((E.n_sites >= 2).sum())}**")
    emit("")
    emit("## `individual_list.txt` 语法示例（前 5 行）")
    emit("")
    emit("```")
    for s in indiv[:5]:
        emit(s)
    emit("```")
    emit("")
    emit("> 语法：`<WT><chain><pdb_resseq><MUT>`，逗号分隔、**无空格**、以 `;` 结尾")
    emit("")
    emit("## 执行顺序（§35）")
    emit("")
    emit("```bash")
    emit("# 1) 每个 target repair 一次（只做一次，缓存复用）")
    emit("EvoEF2 --command=RepairStructure --pdb=<T>_protein.pdb --num_of_runs=3")
    emit("#    -> <T>_protein_Repair.pdb")
    emit("")
    emit("# 2) 逐个变体 build mutant")
    emit("EvoEF2 --command=BuildMutant --pdb=<T>_protein_Repair.pdb \\")
    emit("       --mutant_file=individual_list.txt --num_of_runs=10")
    emit("#    -> <T>_protein_Repair_Model_1.pdb")
    emit("")
    emit("# 3) WT/MUT 各算一次 binding")
    emit("EvoEF2 --command=ComputeBinding --pdb=<T>_protein_Repair.pdb")
    emit("EvoEF2 --command=ComputeBinding --pdb=<T>_protein_Repair_Model_1.pdb")
    emit("```")
    emit("")
    emit("> ⚠️ `BuildMutant` 一次只产出一个 `Model_1.pdb`；多突变**写在同一条 token 内**")
    emit("> （逗号分隔），不是多行。脚本 `12_evoef2_run.py` 需按此逐个调用。")
    emit("")
    emit("> ⚠️ 多链复合物需显式 `--split=<binder>,<partner>`；见 `12_evoef2_run.py`。")
    (PH3 / "evoef2_prepare.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    json.dump({"n_variants": len(E), "n_targets": int(E.target_id.nunique()),
               "n_pdb": n_pdb, "n_skipped": n_missing,
               "individual_list": str(OUT / "individual_list.txt")},
              open(PH3 / "evoef2_prepare_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n>>> variants={len(E)} targets={E.target_id.nunique()} pdb={n_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
