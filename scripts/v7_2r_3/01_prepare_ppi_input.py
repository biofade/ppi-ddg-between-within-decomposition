# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/01_prepare_ppi_input.py -- V7.2R-3 PHASE 0-1.

**修复 V7.2R-2 的链丢弃缺陷**（§10/§29）
----------------------------------------
V7.2R-2 的错误：
    `target_chains.csv` 由 `used = sit[["target_id","chain"]].drop_duplicates()`
    写出 → **只含携带突变位点的链**；notebook 的 `build_variant()` 也只返回
    `sites_by_chain` 中的链 → **伙伴链从未进入前向**。
实测后果：**150/174 target（86.2%）缺失伙伴链，共丢 421 条链**。

本脚本改为：**每个 target 的全部蛋白链都写入包**（`build_complex()` 语义），
突变**只施加**到受影响链，伙伴链保持 WT（§10/§11）。

输出
----
`kaggle/v7_2r_3/kaggle_input_package.zip`（无标签）
`kaggle/v7_2r_3/kaggle_input_manifest.csv`
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

KAG = C.ROOT / "kaggle" / "v7_2r_3"
V723 = C.ROOT / "experiments" / "v7_2r_3"
PH1 = V723 / "phase1_prepare"
for d in (KAG, PH1):
    d.mkdir(parents=True, exist_ok=True)
STAGE = KAG / "staging"
STAGE.mkdir(parents=True, exist_ok=True)
V711 = C.ROOT / "experiments" / "v7_1_1"
PDBDIR = C.DATA / "pdb"
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    sp = importlib.util.spec_from_file_location(
        "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    T2O = v7f.THREE2ONE

    emit("# V7.2R-3 PHASE 0–1 — 修复链丢弃缺陷 + 构建含全部伙伴链的无标签包")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)

    # 只做 V7.2R-2 已成功覆盖的 174 个 target（保持 U_AB 可比）
    bdir = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
    b_rows = set()
    for p in sorted(bdir.glob("*.npz")):
        try:
            b_rows |= set(np.load(p)["row_id"].tolist())
        except Exception:  # noqa: BLE001
            continue
    lab = lab[lab.pair_row.isin(b_rows)].copy()
    emit(f"- 目标宇宙：**{len(lab)}** 行 / **{lab.target_id_canonical.nunique()}** "
         f"targets / **{lab.family_id_canonical.nunique()}** families"
         f"（= V7.2R-2 已覆盖集，保证 U_ABC 可比）")
    emit("")

    # ---------- 全部蛋白链 ----------
    chain_cache: dict = {}

    def all_chains(t):
        k = t.upper()
        if k not in chain_cache:
            p = PDBDIR / f"{k}.pdb"
            if not p.exists():
                chain_cache[k] = {}
            else:
                ch = defaultdict(list)
                for ln in open(p, "r", errors="replace"):
                    if not ln.startswith("ATOM") or ln[12:16].strip() != "CA":
                        continue
                    rn = ln[17:20].strip().upper()
                    if rn not in T2O:
                        continue
                    try:
                        sq = int(ln[22:26])
                    except ValueError:
                        continue
                    ch[ln[21]].append((sq, ln[26].strip(), T2O[rn]))
                chain_cache[k] = dict(ch)
        return chain_cache[k]

    # ---------- 位点（按 row 索引）----------
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}

    targets = sorted(lab.target_id_canonical.astype(str).unique())
    tc_rows, var_rows, site_rows = [], [], []
    ok_rows = 0
    for t in targets:
        ch = all_chains(t)
        if not ch:
            continue
        # **全部链**（修复点）
        for c, res in ch.items():
            tc_rows.append({"target_id": t, "chain": c,
                            "wt_sequence": "".join(x[2] for x in res),
                            "length": len(res),
                            "pdb_resseq_first": res[0][0],
                            "pdb_resseq_last": res[-1][0],
                            "is_protein_chain": 1})
        sub = lab[lab.target_id_canonical.astype(str) == t]
        for r in sub.itertuples():
            pr = int(r.pair_row)
            g = s_by_row.get(pr)
            if g is None or not len(g):
                continue
            acc = []
            for x in g.itertuples():
                c = str(x.ann_chain); i = int(x.seq_index0)
                if c not in ch or not (0 <= i < len(ch[c])):
                    continue
                if ch[c][i][2] != str(x.wt_aa):      # §19 绝不接受 WT mismatch
                    continue
                acc.append({"chain": c, "seq_index0": i,
                            "pdb_resseq": ch[c][i][0], "icode": ch[c][i][1],
                            "wt_aa": str(x.wt_aa), "mut_aa": str(x.mut_aa)})
            if not acc:
                continue
            ok_rows += 1
            for k, a in enumerate(acc):
                site_rows.append({"row_id": pr, "target_id": t,
                                  "family_id": str(r.family_id_canonical),
                                  "site_index": k, "chain": a["chain"],
                                  "seq_index0": a["seq_index0"],
                                  "pdb_resseq": a["pdb_resseq"],
                                  "icode": a["icode"], "wt_aa": a["wt_aa"],
                                  "mut_aa": a["mut_aa"]})
            var_rows.append({"row_id": pr, "target_id": t,
                             "family_id": str(r.family_id_canonical),
                             "n_sites": len(acc),
                             "mutated_chains": "|".join(sorted({a["chain"]
                                                                for a in acc})),
                             "mutation_bodies": ";".join(
                                 f"{a['wt_aa']}{a['pdb_resseq']}{a['mut_aa']}"
                                 for a in acc)})

    TC = pd.DataFrame(tc_rows); VAR = pd.DataFrame(var_rows)
    SIT = pd.DataFrame(site_rows)
    emit("## 1.1 链覆盖（修复后）")
    emit("")
    emit("| 项 | V7.2R-2（错） | V7.2R-3（修复） |")
    emit("|---|---:|---:|")
    emit(f"| (target,chain) 记录 | 233 | **{len(TC)}** |")
    emit(f"| 每 target 链数中位 | 1 | **{TC.groupby('target_id').size().median():.0f}** |")
    emit(f"| 行覆盖 | 4219 | **{ok_rows}** |")
    emit("")
    lab_out = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
        columns={"pair_row": "row_id", "target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})
    outp = {"canonical_rows.csv": lab_out,
            "variants.csv": VAR,
            "mutation_sites.csv": SIT,
            "target_chains.csv": TC,
            "variant_reconstruction.csv":
                VAR[["row_id", "target_id", "n_sites", "mutated_chains",
                     "mutation_bodies"]]}
    for n, df in outp.items():
        df.to_csv(STAGE / n, index=False, encoding="utf-8")

    # 无标签检查
    LABTOK = ["ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred",
              "score", "ipsae"]
    viol = [(n, c) for n, df in outp.items() for c in df.columns
            if any(t in c.lower() for t in LABTOK)]
    emit(f"- 无标签检查：**{len(viol)}** 处违规 "
         f"{viol if viol else '→ ✅ 无'}")
    emit("")
    readme = STAGE / "PACKAGE_README.md"
    readme.write_text(
        "# V7.2R-3 Kaggle Input Package（含**全部伙伴链**，无标签）\n\n"
        "## ⚠️ 与 V7.2R-2 的关键差别\n\n"
        "V7.2R-2 的 `target_chains.csv` **只含携带突变位点的链**，"
        "导致 86.2% 的 target 缺失伙伴链。\n"
        "本包的 `target_chains.csv` **含该 PDB 的每一条蛋白链**。\n\n"
        "## 复合物构建（强制）\n\n```\n"
        "WT  复合物 = 全部链的 WT 序列\n"
        "MUT 复合物 = 受影响链施加该 row_id 的**全部**替换，伙伴链保持 WT\n```\n\n"
        "**伙伴链必须保留在前向输入中**，否则不是 PPI-conditioned。\n\n"
        "## 位点定位\n\n`seq_index0` 为 0-based 链内下标；"
        "`target_chains.csv[wt_sequence][seq_index0]` 必须等于 `mutation_sites.csv[wt_aa]`。\n\n"
        "## 无标签\n\n不含 ΔΔG / KD(M) / 任何 V7.x 预测。\n",
        encoding="utf-8")
    man = []
    for n in list(outp) + ["PACKAGE_README.md"]:
        p = STAGE / n
        man.append({"file": n, "size_bytes": p.stat().st_size,
                    "sha256_prefix16": sha(p)[:16], "contains_labels": "NO"})
    pd.DataFrame(man).to_csv(KAG / "kaggle_input_manifest.csv", index=False,
                             encoding="utf-8")
    zp = KAG / "kaggle_input_package.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for n in list(outp) + ["PACKAGE_README.md"]:
            z.write(STAGE / n, arcname=n)
        z.write(KAG / "kaggle_input_manifest.csv",
                arcname="kaggle_input_manifest.csv")
    emit("## 1.2 包内容")
    emit("")
    emit("| file | bytes | sha256(16) | 含标签 |")
    emit("|---|---:|---|---|")
    for r in man:
        emit(f"| `{r['file']}` | {r['size_bytes']:,} | `{r['sha256_prefix16']}` | "
             f"**{r['contains_labels']}** |")
    emit("")
    emit(f"- **ZIP**: `kaggle_input_package.zip`（{zp.stat().st_size / 1e6:.2f} MB）")
    emit(f"- **SHA256**: `{sha(zp)}`")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "prepare_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {zp} ({zp.stat().st_size / 1e6:.2f} MB)")
    print(f">>> rows={ok_rows} targets={VAR.target_id.nunique()} "
          f"(target,chain)={len(TC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
