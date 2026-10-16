# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/01_prepare_kaggle_input.py -- V7.2R-2：构建 **无标签** Kaggle 推理包.

原则（§6）
---------
* **只含推理必需字段**：canonical row/target/family id、WT/突变体链序列、
  链分配、突变位点、期望 WT、partner 链序列、row→target 映射、多突变重建元数据
* **绝不包含**：ΔΔG 标签、任何测试标签、V7.1/V7.2 预测、outer 性能、bootstrap 结果
* target identity **只**用 `canonical_target_labels.csv`（§7），按 **row_id** 连接
* 输出 `kaggle_input_package.zip` + `kaggle_input_manifest.csv`

输出
----
`kaggle/v7_2r_2/kaggle_input_package.zip`
`kaggle/v7_2r_2/kaggle_input_manifest.csv`
`experiments/v7_2r_2/phase0_prepare/prepare_report.md`
"""

from __future__ import annotations

import hashlib
import io
import json
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

KAG = C.ROOT / "kaggle" / "v7_2r_2"
V722 = C.ROOT / "experiments" / "v7_2r_2"
PH0 = V722 / "phase0_prepare"
for d in (KAG, PH0):
    d.mkdir(parents=True, exist_ok=True)
STAGE = KAG / "staging"
STAGE.mkdir(parents=True, exist_ok=True)
V711 = C.ROOT / "experiments" / "v7_1_1"
PDBDIR = C.DATA / "pdb"
T2O = None
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
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    T2O = v7f.THREE2ONE

    emit("# V7.2R-2 PHASE 0 — 构建无标签 Kaggle 推理包")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------- 权威标签（只用这个做 target identity）----------
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    fwd = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = fwd[fwd.direction == "forward"].reset_index(drop=True)
    emit(f"- canonical rows: **{len(lab)}** / targets "
         f"**{lab.target_id_canonical.nunique()}** / families "
         f"**{lab.family_id_canonical.nunique()}**")
    emit("")

    # ---------- 只保留 Track A 已覆盖的行（与完整合约可比）----------
    v72per = C.ROOT / "features" / "v7_2" / "per_target"
    have = set()
    for p in sorted(v72per.glob("*.npz")):
        try:
            have |= set(np.load(p)["row_id"].tolist())
        except Exception:  # noqa: BLE001
            continue
    lab = lab[lab.pair_row.isin(have)].copy()
    emit(f"- 限定为 Track A 已覆盖行（保证 U_AB 可比）：**{len(lab)}** 行 / "
         f"**{lab.target_id_canonical.nunique()}** targets / "
         f"**{lab.family_id_canonical.nunique()}** families")
    emit("")

    # ---------- 链序列（WT）----------
    chain_cache: dict = {}

    def chains_of(t):
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

    # ---------- 构建「变体」表 ----------
    s_by_row: dict = {}
    for pr, gg in site.groupby("pair_row"):
        s_by_row[int(pr)] = gg

    sites_rows, variant_rows = [], []
    tgt_chain_rows = []
    ok_rows, no_site = [], 0
    for r in lab.itertuples():
        pr = int(r.pair_row); t = str(r.target_id_canonical)
        ch = chains_of(t)
        if not ch:
            no_site += 1
            continue
        g = s_by_row.get(pr)
        if g is None or not len(g):
            no_site += 1
            continue
        acc = []
        for x in g.itertuples():
            c = str(x.ann_chain); i = int(x.seq_index0)
            if c not in ch or not (0 <= i < len(ch[c])):
                continue
            resseq, icode, aa = ch[c][i]
            if aa != str(x.wt_aa):
                continue
            acc.append({"chain": c, "seq_index0": i, "pdb_resseq": resseq,
                        "icode": icode, "wt_aa": aa, "mut_aa": str(x.mut_aa)})
        if not acc:
            no_site += 1
            continue
        ok_rows.append(pr)
        # 每 variant 一条记录 + 每 site 一条记录
        for k, a in enumerate(acc):
            sites_rows.append({
                "row_id": pr, "target_id": t,
                "family_id": str(r.family_id_canonical),
                "site_index": k, "chain": a["chain"],
                "seq_index0": a["seq_index0"], "pdb_resseq": a["pdb_resseq"],
                "icode": a["icode"], "wt_aa": a["wt_aa"],
                "mut_aa": a["mut_aa"]})
        variant_rows.append({
            "row_id": pr, "target_id": t,
            "family_id": str(r.family_id_canonical),
            "n_sites": len(acc),
            "chains": "|".join(sorted({a["chain"] for a in acc})),
            "mutation_bodies": ";".join(
                f"{a['wt_aa']}{a['pdb_resseq']}{a['mut_aa']}" for a in acc)})
    emit(f"- 可用变体（WT 校验通过）: **{len(variant_rows)}** 行；"
         f"跳过（无可用位点）: {no_site}")
    emit("")

    # ---------- WT / 突变体链序列 ----------
    var = pd.DataFrame(variant_rows)
    sit = pd.DataFrame(sites_rows)
    # 逐 (target, chain) 写 WT 序列
    used = sit[["target_id", "chain"]].drop_duplicates()
    for u in used.itertuples():
        ch = chains_of(str(u.target_id))
        if str(u.chain) not in ch:
            continue
        seq = "".join(x[2] for x in ch[str(u.chain)])
        tgt_chain_rows.append({"target_id": str(u.target_id),
                               "chain": str(u.chain),
                               "wt_sequence": seq,
                               "length": len(seq),
                               "pdb_resseq_first": ch[str(u.chain)][0][0],
                               "pdb_resseq_last": ch[str(u.chain)][-1][0]})
    tc = pd.DataFrame(tgt_chain_rows)
    emit(f"- 链序列表: **{len(tc)}** 条 (target,chain) 记录")
    emit("")

    # ---------- 写出 staging ----------
    outp = {}
    lab_out = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
        columns={"pair_row": "row_id", "target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})
    outp["canonical_rows.csv"] = lab_out
    outp["variants.csv"] = var
    outp["mutation_sites.csv"] = sit
    outp["target_chains.csv"] = tc
    # 多突变重建元数据（无标签）
    outp["variant_reconstruction.csv"] = var[["row_id", "target_id", "n_sites",
                                              "chains", "mutation_bodies"]]
    for n, df in outp.items():
        df.to_csv(STAGE / n, index=False, encoding="utf-8")

    # ---------- 反向检查：**确认包内无任何标签** ----------
    emit("## 0.1 无标签保证（强制检查）")
    emit("")
    LABEL_TOKENS = ["ddg", "dg", "kd(", "affinity", "label", "y_true", "target_y",
                    "delta_g", "ddelta", "score", "pred", "ipsae", "pae"]
    violations = []
    for n, df in outp.items():
        for c in df.columns:
            cl = c.lower()
            if any(tok in cl for tok in LABEL_TOKENS):
                violations.append((n, c))
    emit("| 检查 | 结果 |")
    emit("|---|---|")
    emit(f"| 包内列名含标签类词 | **{len(violations)}** 处 "
         f"{violations if violations else '→ ✅ 无'} |")
    for tok in ("ddG", "KD(M)", "affinity"):
        emit(f"| 未包含原始来源文件（含 `{tok}`） | ✅ 仅导出剥离后的推理表 |")
    emit("| 未包含 V7.1/V7.2 预测 | ✅ |")
    emit("| 未包含 bootstrap / outer 性能 | ✅ |")
    emit("")

    # ---------- manifest ----------
    man = []
    for n, df in outp.items():
        p = STAGE / n
        man.append({"file": n, "rows": len(df), "cols": df.shape[1],
                    "size_bytes": p.stat().st_size,
                    "sha256_prefix16": sha(p)[:16],
                    "contains_labels": "NO"})
    # 说明文件
    readme = KAG / "staging" / "PACKAGE_README.md"
    readme.write_text(
        "# V7.2R-2 Kaggle Input Package（**无标签推理包**）\n\n"
        "## 内容\n\n"
        "| file | 说明 |\n|---|---|\n"
        "| `canonical_rows.csv` | `row_id` / `target_id` / `family_id`"
        "（**权威命名空间**，按 row_id 连接） |\n"
        "| `variants.csv` | 每个变体的 row_id、target、family、位点数、涉及链 |\n"
        "| `mutation_sites.csv` | 每个位点：链、`seq_index0`（0-based 链内下标）、"
        "`pdb_resseq`、`icode`、**期望 WT**、突变后残基 |\n"
        "| `target_chains.csv` | 每个 (target, chain) 的 **WT 链序列**（由 PDB CA 记录重建） |\n"
        "| `variant_reconstruction.csv` | 多突变重建元数据（`mutation_bodies`） |\n\n"
        "## ⚠️ 本包含有**零个标签**\n\n"
        "- **不含** ΔΔG / KD(M) / 任何亲和力标签\n"
        "- **不含** 任何 V7.x 预测或 outer-test 性能\n"
        "- 突变体序列**必须**由 `target_chains.csv` + `mutation_sites.csv` 现场重建\n\n"
        "## 位点定位约定（关键）\n\n"
        "`seq_index0` 是 **0-based，指向该链在 PDB 中 CA 记录序列的第几个残基**"
        "（链内下标），**不是** PDB 编号。\n\n"
        "⇒ 从 `target_chains.csv` 取 `wt_sequence`，"
        "`wt_sequence[seq_index0]` **必须**等于 `mutation_sites.csv` 的 `wt_aa`。\n"
        "若不等，**该位点必须被拒绝**（不要静默修复）。\n\n"
        "## 突变体序列构建\n\n"
        "对每个 `row_id`：拷贝所有涉及链的 WT 序列，"
        "把该 row_id 的**全部**位点一次性加上去（多突变**同时**施加），\n"
        "得到突变体链序列。\n",
        encoding="utf-8")
    man.append({"file": "PACKAGE_README.md", "rows": 0, "cols": 0,
                "size_bytes": readme.stat().st_size,
                "sha256_prefix16": sha(readme)[:16],
                "contains_labels": "NO"})
    M = pd.DataFrame(man)
    M.to_csv(KAG / "kaggle_input_manifest.csv", index=False, encoding="utf-8")

    # ---------- zip ----------
    zpath = KAG / "kaggle_input_package.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for n in list(outp.keys()) + ["PACKAGE_README.md"]:
            z.write(STAGE / n, arcname=n)
        z.write(KAG / "kaggle_input_manifest.csv",
                arcname="kaggle_input_manifest.csv")
    zsha = sha(zpath)
    emit("## 0.2 包内容")
    emit("")
    emit("| file | rows | cols | bytes | sha256(16) | 含标签 |")
    emit("|---|---:|---:|---:|---|---|")
    for r in man:
        emit(f"| `{r['file']}` | {r['rows']} | {r['cols']} | {r['size_bytes']:,} | "
             f"`{r['sha256_prefix16']}` | **{r['contains_labels']}** |")
    emit("")
    emit(f"- **ZIP**: `kaggle_input_package.zip`（{zpath.stat().st_size / 1e6:.2f} MB）")
    emit(f"- **SHA256**: `{zsha}`")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH0 / "prepare_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    json.dump({"zip": str(zpath), "sha256": zsha, "n_rows": len(var),
               "n_targets": int(var.target_id.nunique()),
               "n_families": int(var.family_id.nunique()),
               "zip_bytes": zpath.stat().st_size},
              open(PH0 / "package_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n写出 {zpath}  ({zpath.stat().st_size / 1e6:.2f} MB)")
    print(f"写出 {KAG / 'kaggle_input_manifest.csv'}")
    print(f"写出 {PH0 / 'prepare_report.md'}")
    print(f">>> rows={len(var)} targets={var.target_id.nunique()} "
          f"families={var.family_id.nunique()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
