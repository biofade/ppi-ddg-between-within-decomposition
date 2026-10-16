# -*- coding: utf-8 -*-
"""scripts/v7_4/01_prepare_relax_input.py -- V7.4 PHASE 1：U_RELAX_STRICT + Kaggle 输入包.

从 **U_PHYS_STRICT** 出发（§9），追加**松弛安全**要求：
    * 突变不落在插入码残基上（§8 已知风险）
    * 突变残基的**骨架原子 N/CA/C/O 齐全**
    * 突变是 20 种标准氨基酸之一
    * 位点周围（8 Å 内）**没有链断裂**（局部骨架可用）
    * 仅蛋白链（§10）
多突变**全部**位点须通过。

**不做**任何静默修复（§8/§11）。

输出（§28/§62）
--------------
`kaggle/v7_4/kaggle_input_package.zip`（**无标签**）
`kaggle/v7_4/kaggle_input_manifest.csv`
`features/v7_4/u_relax_strict_manifest.csv`
`experiments/v7_4/phase1_prepare/prepare_report.md`
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

V711 = C.ROOT / "experiments" / "v7_1_1"
V74 = C.ROOT / "experiments" / "v7_4"
PH1 = V74 / "phase1_prepare"
PH1.mkdir(parents=True, exist_ok=True)
FEAT = C.ROOT / "features" / "v7_4"
FEAT.mkdir(parents=True, exist_ok=True)
KAG = C.ROOT / "kaggle" / "v7_4"
KAG.mkdir(parents=True, exist_ok=True)
STAGE = KAG / "staging"
STAGE.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
L: list[str] = []
# 20 标准氨基酸（三字母）
STD3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


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
    ONE2THREE = {v: k for k, v in T2O.items()}

    emit("# V7.4 PHASE 1 — U_RELAX_STRICT + Kaggle 输入包")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    inc = pd.read_csv(C.ROOT / "features" / "v7_3" /
                      "u_phys_strict_manifest.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    s_by_row = {int(k): g for k, g in site.groupby("pair_row")}
    inc_rows = set(inc.row_id.astype(int))
    lab = lab[lab.pair_row.isin(inc_rows)].copy()
    emit(f"- 起点 U_PHYS_STRICT：**{len(lab)}** 行 / "
         f"{lab.target_id_canonical.nunique()} targets")
    emit("")

    cache = {}          # target -> 解析结果
    exc = defaultdict(int)
    exc_rows, keep_rows = [], []

    def parse(t):
        if t in cache:
            return cache[t]
        p = PDBDIR / f"{t.upper()}.pdb"
        if not p.exists():
            cache[t] = None
            return None
        ch_res = defaultdict(list)      # chain -> [(resseq,icode,aa3)]
        ch_atoms = defaultdict(list)    # chain -> [ {atom:xyz} ]
        ch_order = defaultdict(list)    # chain -> [resseq+icode] 顺序
        for ln in open(p, "r", errors="replace"):
            if ln[:6] not in ("ATOM  ", "HETATM"):
                continue
            elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if elem == "H":
                continue
            rn = ln[17:20].strip().upper()
            if rn not in STD3:
                continue
            try:
                x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            except ValueError:
                continue
            c = ln[21]
            try:
                rs = int(ln[22:26])
            except ValueError:
                continue
            ic = ln[26].strip()
            key = (rs, ic, rn)
            if not ch_res[c] or ch_res[c][-1] != key:
                ch_res[c].append(key); ch_atoms[c].append({})
                ch_order[c].append(f"{rs}{ic}")
            an = ln[12:16].strip()
            ch_atoms[c][-1][an] = (x, y, z)
        out = {"res": ch_res, "atoms": ch_atoms, "order": ch_order}
        cache[t] = out
        return out

    for ti, r in enumerate(lab.itertuples()):
        t = str(r.target_id_canonical)
        pr = int(r.pair_row)
        P = parse(t)
        if P is None:
            exc["missing_pdb"] += 1
            exc_rows.append({"row_id": pr, "target": t, "reason": "missing_pdb"})
            continue
        g = s_by_row.get(pr)
        if g is None or not len(g):
            exc["no_site"] += 1
            exc_rows.append({"row_id": pr, "target": t, "reason": "no_site"})
            continue
        bad, sites = None, []
        for x in g.itertuples():
            c = str(x.ann_chain); i = int(x.seq_index0)
            if c not in P["res"]:
                bad = "missing_chain"; break
            if not (0 <= i < len(P["res"][c])):
                bad = "idx_oob"; break
            rs, ic, aa3 = P["res"][c][i]
            if aa3 not in STD3:
                bad = "nonstandard_residue"; break
            if str(x.wt_aa) != T2O[aa3]:
                bad = "wt_mismatch"; break
            if str(x.mut_aa) not in list(T2O.values()):
                bad = "nonstandard_mutant"; break
            if ic != "":
                bad = "insertion_code_at_site"; break
            at = P["atoms"][c][i]
            if not all(k in at for k in ("N", "CA", "C", "O")):
                bad = "missing_backbone"; break
            sites.append({"chain": c, "seq_index0": i, "resseq": rs,
                          "icode": ic, "aa3": aa3, "wt": str(x.wt_aa),
                          "mut": str(x.mut_aa)})
        if bad is None:
            # 局部骨架连续性：8 Å 内同链邻居索引必须连续（无断裂）
            for s in sites:
                c = s["chain"]; i = s["seq_index0"]
                lo, hi = max(0, i - 12), min(len(P["res"][c]) - 1, i + 12)
                for j in range(lo, hi):
                    a = P["atoms"][c][j]; b = P["atoms"][c][j + 1]
                    if "CA" not in a or "CA" not in b:
                        bad = f"gap_near_site_{c}"; break
                    d = float(np.linalg.norm(np.array(a["CA"]) - np.array(b["CA"])))
                    if d > 5.0:          # CA-CA 相邻 > 5 Å ⇒ 断裂
                        bad = f"chain_break_near_site_{c}"; break
                if bad:
                    break
        if bad is not None:
            exc[bad] += 1
            exc_rows.append({"row_id": pr, "target": t, "reason": bad})
            continue
        keep_rows.append({"row_id": pr, "target_id_canonical": t,
                          "family_id_canonical": str(r.family_id_canonical),
                          "n_sites": len(sites),
                          "mutated_chains": "|".join(sorted({s["chain"]
                                                             for s in sites})),
                          "mutation_bodies": ";".join(
                              f"{s['wt']}{s['resseq']}{s['mut']}" for s in sites),
                          "sites": sites})
        if (ti + 1) % 1000 == 0:
            emit(f"  - {ti + 1}/{len(lab)} | keep={len(keep_rows)} | "
                 f"{(time.time() - t0) / 60:.1f} min")

    emit("")
    emit("## 1.1 U_RELAX_STRICT（§9）")
    emit("")
    K = pd.DataFrame([{k: v for k, v in r.items() if k != "sites"}
                      for r in keep_rows])
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| U_PHYS_STRICT rows | {len(lab)} |")
    emit(f"| **U_RELAX_STRICT rows** | **{len(K)}** |")
    emit(f"| targets | **{K.target_id_canonical.nunique()}** |")
    emit(f"| families | **{K.family_id_canonical.nunique()}** |")
    for k in (1, 2, 3, 5, 10):
        emit(f"| n>={k} targets | "
             f"**{int((K.groupby('target_id_canonical').size() >= k).sum())}** |")
    emit("")
    emit("### 追加排除原因（V7.4 §9）")
    emit("")
    emit("| reason | rows |")
    emit("|---|---:|")
    for k, v in sorted(exc.items(), key=lambda x: -x[1]):
        emit(f"| `{k}` | {v} |")
    emit(f"| **额外排除** | **{sum(exc.values())}** |")
    emit("")
    if exc_rows:
        pd.DataFrame(exc_rows).to_csv(PH1 / "relax_exclusions.csv", index=False,
                                      encoding="utf-8")
    K.to_csv(FEAT / "u_relax_strict_manifest.csv", index=False, encoding="utf-8")

    # ---------- Kaggle 输入包（只含推理必需，无标签）----------
    # 需要：pdb_id（供 Kaggle 从 RCSB 下载）、链序列、位点、期望 WT/突变
    tc_rows, sit_rows = [], []
    used = defaultdict(set)
    for r in keep_rows:
        for s in r["sites"]:
            used[r["target_id_canonical"]].add(s["chain"])
    for t, chs in used.items():
        P = parse(t)
        for c in sorted(chs):
            seq = "".join(T2O[a3] for (_, _, a3) in P["res"][c])
            tc_rows.append({"target_id": t, "chain": c, "wt_sequence": seq,
                            "length": len(seq)})
    for r in keep_rows:
        for k, s in enumerate(r["sites"]):
            sit_rows.append({"row_id": r["row_id"], "target_id": r["target_id_canonical"],
                             "family_id": r["family_id_canonical"],
                             "site_index": k, "chain": s["chain"],
                             "seq_index0": s["seq_index0"],
                             "pdb_resseq": s["resseq"], "icode": s["icode"],
                             "wt_aa": s["wt"], "mut_aa": s["mut"]})
    VAR = K[["row_id", "target_id_canonical", "family_id_canonical", "n_sites",
             "mutated_chains", "mutation_bodies"]].rename(
        columns={"target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})
    SIT = pd.DataFrame(sit_rows)
    TC = pd.DataFrame(tc_rows)
    outp = {"variants.csv": VAR, "mutation_sites.csv": SIT,
            "target_chains.csv": TC,
            "canonical_rows.csv": K[["row_id", "target_id_canonical",
                                     "family_id_canonical"]].rename(
                columns={"target_id_canonical": "target_id",
                         "family_id_canonical": "family_id"})}
    for n, df in outp.items():
        df.to_csv(STAGE / n, index=False, encoding="utf-8")
    readme = STAGE / "PACKAGE_README.md"
    readme.write_text(
        "# V7.4 Kaggle Input Package（OpenMM 松弛，**无标签**）\n\n"
        "## 内容\n\n"
        "| file | 说明 |\n|---|---|\n"
        "| `canonical_rows.csv` | `row_id` / `target_id` / `family_id`（权威命名空间）|\n"
        "| `variants.csv` | 每变体的位点数、涉及链、突变体字符串 |\n"
        "| `mutation_sites.csv` | 位点：链、`seq_index0`（0-based 链内下标）、"
        "`pdb_resseq`、`icode`、期望 WT、突变后残基 |\n"
        "| `target_chains.csv` | 每 (target, chain) 的 **WT 链序列** |\n\n"
        "## 结构获取\n\n"
        "`target_id` 即 **PDB ID**（如 `1A22`）⇒ 从 `https://files.rcsb.org/download/"
        "<ID>.pdb` 下载。**不要**把完整 PDB 打进包里。\n\n"
        "## 突变体构建（强制）\n\n"
        "```\nWT  体系 = 全部蛋白链的 WT 序列（去掉水/配体/离子）\n"
        "MUT 体系 = 受影响链施加该 row_id 的**全部**替换，伙伴链保持 WT\n```\n\n"
        "## ⚠️ 无标签\n\n不含 ΔΔG / KD(M) / 任何 V7.x 预测或性能。\n", encoding="utf-8")
    man = []
    for n in list(outp) + ["PACKAGE_README.md"]:
        p = STAGE / n
        man.append({"file": n, "size_bytes": p.stat().st_size,
                    "sha256_prefix16": sha(p)[:16], "contains_labels": "NO"})
    M = pd.DataFrame(man)
    M.to_csv(KAG / "kaggle_input_manifest.csv", index=False, encoding="utf-8")
    zp = KAG / "kaggle_input_package.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in list(outp) + ["PACKAGE_README.md"]:
            zf.write(STAGE / n, arcname=n)
        zf.write(KAG / "kaggle_input_manifest.csv",
                 arcname="kaggle_input_manifest.csv")

    # 无标签检查
    LABTOK = ["ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred"]
    viol = [(n, c) for n, df in outp.items() for c in df.columns
            if any(t in c.lower() for t in LABTOK)]
    emit("## 1.2 Kaggle 输入包")
    emit("")
    emit("| file | bytes | 含标签 |")
    emit("|---|---:|---|")
    for r in man:
        emit(f"| `{r['file']}` | {r['size_bytes']:,} | **{r['contains_labels']}** |")
    emit("")
    emit(f"- ZIP：`kaggle_input_package.zip`（{zp.stat().st_size / 1e6:.2f} MB）")
    emit(f"- SHA256：`{sha(zp)}`")
    emit(f"- **无标签检查：{len(viol)} 处违规 {viol if viol else '→ ✅ 无'}**")
    emit(f"- (target,chain) 记录：**{len(TC)}**；target 数 **{TC.target_id.nunique()}**")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH1 / "prepare_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    import json
    json.dump({"zip": str(zp), "sha256": sha(zp), "n_rows": len(K),
               "n_targets": int(K.target_id_canonical.nunique()),
               "n_families": int(K.family_id_canonical.nunique()),
               "n_chains": len(TC), "zip_bytes": zp.stat().st_size},
              open(PH1 / "package_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n写出 u_relax_strict_manifest.csv / kaggle_input_package.zip")
    print(f">>> U_RELAX_STRICT rows={len(K)} targets="
          f"{K.target_id_canonical.nunique()} chains={len(TC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
