# -*- coding: utf-8 -*-
"""scripts/v7_4r2/02_prepare_colab.py -- V7.4R2：Colab **无标签**输入包（§26）.

来源：`features/v7_4r2/canonical_wt_strict_manifest.csv`（CANONICAL_WT_STRICT）
输出：`colab/v7_4r2/colab_input_package.zip`

**关键**：位点表以 **`seq_index0`** 为权威寻址（V7.3 做法，实测 100% WT 匹配），
同时保留 `pdb_resseq` 供核对。
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
PH2 = R / "experiments" / "v7_4r2" / "phase2_prepare"
PH2.mkdir(parents=True, exist_ok=True)
COLAB = R / "colab" / "v7_4r2"
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
    O2T = {v: k for k, v in T2O.items()}

    emit("# V7.4R2 PHASE 2 — Colab 输入包（无标签）")
    emit("")
    man = pd.read_csv(R / "features" / "v7_4r2" /
                      "canonical_wt_strict_manifest.csv", low_memory=False)
    emit(f"- CANONICAL_WT_STRICT：**{len(man)}** 行 / "
         f"**{man.target_id_canonical.nunique()}** targets / "
         f"**{man.family_id_canonical.nunique()}** families")

    # 链序列（只用观测残基；不补缺口）
    tc, ch_cache = [], {}
    for t in sorted(man.target_id_canonical.astype(str).unique()):
        p = PDBDIR / f"{t.upper()}.pdb"
        if not p.exists():
            continue
        res = defaultdict(list)
        for ln in open(p, "r", errors="replace"):
            if ln[:6] != "ATOM  ":
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
            if not res[c] or res[c][-1][0] != key:
                res[c].append((key, {}))
        ch_cache[t] = res
        for c in sorted(res):
            seq = "".join(T2O[k[2]] for k, _ in res[c])
            tc.append({"target_id": t, "chain": c, "wt_sequence": seq,
                       "length": len(seq)})
    TC = pd.DataFrame(tc)
    emit(f"- 解析 (target,chain)：**{len(TC)}** 条 / {TC.target_id.nunique()} targets")

    site = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    sit_rows = []
    for r in man.itertuples():
        rid = int(r.row_id)
        g = s_by_row.get(rid)
        if g is None or not len(g):
            continue
        for k, x in enumerate(g.itertuples()):
            t = str(r.target_id_canonical); c = str(x.ann_chain)
            si = int(x.seq_index0)
            # 记录 PDB 中该 index 的真实编号与插入码（供 canonical WT 克隆使用）
            real_rs, real_ic = None, ""
            if t in ch_cache and c in ch_cache[t] and 0 <= si < len(ch_cache[t][c]):
                real_rs, real_ic = ch_cache[t][c][si][0][0], ch_cache[t][c][si][0][1]
            sit_rows.append({"row_id": rid, "target_id": t,
                             "family_id": str(r.family_id_canonical),
                             "site_index": k, "chain": c,
                             "seq_index0": si,
                             "pdb_resseq": real_rs,
                             "icode": real_ic,
                             "wt_aa": str(x.wt_aa), "mut_aa": str(x.mut_aa)})
    SIT = pd.DataFrame(sit_rows)
    emit(f"- 位点行：**{len(SIT)}**（以 `seq_index0` 为权威寻址）")

    VAR = man[["row_id", "target_id_canonical", "family_id_canonical", "n_sites",
               "mutated_chains", "mutation_bodies"]].rename(
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
        "# V7.4R2 Colab Input Package（Canonical-WT 局部松弛，**无标签**）\n\n"
        "| file | 说明 |\n|---|---|\n"
        "| `canonical_rows.csv` | `row_id` / `target_id` / `family_id` |\n"
        "| `variants.csv` | 每变体位点数、涉及链、突变体串 |\n"
        "| `mutation_sites.csv` | **`seq_index0` 为权威寻址**；另附 PDB 真实 "
        "`pdb_resseq` 与 `icode` |\n"
        "| `target_chains.csv` | 每 (target, chain) 的**观测** WT 序列（不补缺口）|\n\n"
        "结构：`target_id` 即 PDB ID，从 "
        "`https://files.rcsb.org/download/<ID>.pdb` 下载。\n\n"
        "**无 ΔΔG / KD / 任何标签。**\n", encoding="utf-8")

    man_rows = []
    for n in list(outp) + ["PACKAGE_README.md"]:
        p = STAGE / n
        man_rows.append({"file": n, "size_bytes": p.stat().st_size,
                         "sha256_prefix16": sha(p)[:16],
                         "contains_labels": "NO"})
    pd.DataFrame(man_rows).to_csv(COLAB / "input_manifest.csv", index=False,
                                  encoding="utf-8")
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
    (PH2 / "prepare_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    import json
    json.dump({"zip": str(zp), "sha256": sha(zp), "n_rows": len(man),
               "n_targets": int(man.target_id_canonical.nunique()),
               "n_families": int(man.family_id_canonical.nunique()),
               "n_chains": len(TC), "n_sites": len(SIT),
               "zip_bytes": zp.stat().st_size},
              open(PH2 / "package_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n>>> rows={len(man)} targets={man.target_id_canonical.nunique()} "
          f"chains={len(TC)} sites={len(SIT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
