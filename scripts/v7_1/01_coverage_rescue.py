# -*- coding: utf-8 -*-
"""scripts/v7_1/01_coverage_rescue.py -- V7.1 PHASE 1: 全量覆盖审计 + 救援.

V7.0 的 BLOCK L 只覆盖 4473/7523（59.5%）。本脚本查明**确切排除机制**并救援。

审计发现的排除机制（两类，量级差异极大）
------------------------------------------
| 类别 | target 数 | 行数 |
|---|---:|---:|
| **A. 完全排除**（该 target 无任何已验证位点） | **16** | **107** |
| **B. 部分排除**（target 有一些位点验证通过，但**部分 variant** 的位点未通过） | 136 | **2943** |

⇒ V7.0 的损失**主要不是来自 16 个坏 target（只占 107 行）**，而是来自
**多突变 variant 中「部分位点未通过 WT 校验」导致整条 variant 被丢弃**。

救援策略（科学上有效，不伪造映射）
----------------------------------
对每个 variant：
* 只要有 **≥1 个位点通过 WT 校验**，就用这些位点构造 BLOCK L（**部分覆盖**）；
* 记录 `n_sites_verified / n_sites_total`、`site_coverage_frac`、`partial_mapping` 标志；
* 若 `n_sites_total > n_sites_verified`，标 `partial_mapping=1`（供敏感性分析剔除）。

**不做**的：不强行使用未通过 WT 校验的位点（Scientific correctness > coverage）。

两个 universe（§1.5）
--------------------
* **U_FULL** = 全部 7523 行（凡不需要 Local 的模型都用它）
* **U_LOCAL** = 救援后有 BLOCK L 的行
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec)
sys.modules["v7f"] = v7f
_spec.loader.exec_module(v7f)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V70 = C.ROOT / "experiments" / "v7_0"
PH0 = C.ROOT / "experiments" / "v7_1" / "phase0_provenance"
PH1 = C.ROOT / "experiments" / "v7_1" / "phase1_coverage"
for d in (PH0, PH1):
    d.mkdir(parents=True, exist_ok=True)
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


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def chain_residues(path: Path):
    from collections import defaultdict
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
    emit("# V7.1 PHASE 1 — 全量覆盖审计 + 救援")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ================= PHASE 0 =================
    emit("## PHASE 0 — 输入冻结")
    emit("")
    FILES = [
        ("v62_matched_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
        ("v70_block_L", "experiments/v7_0/phase1_features/block_L_local.csv"),
        ("v70_block_S", "experiments/v7_0/phase1_features/block_S_set.csv"),
        ("v70_block_G", "experiments/v7_0/phase1_features/block_G_global.csv"),
        ("v70_matched_rows", "experiments/v7_0/phase1_features/matched_coverage_rows.csv"),
        ("v70_ablation_pred", "experiments/v7_0/phase2_ablation/ablation_predictions.npz"),
        ("v70_loss_pred", "experiments/v7_0/phase3_losses/loss_predictions.npz"),
        ("v621_donly_pred", "experiments/v6_2_1/phase3_d_only/d_only_predictions.npz"),
        ("v63_target_table", "experiments/v6_3/phase1_target_table/target_mechanism_table.csv"),
        ("ddg_split", "data/v4.0/ddg_split.json"),
        ("mutation_sites", "data/v5_0/mutation_sites.csv"),
        ("v22_features", "features/v2.2_proteinea_features.csv"),
        ("esm2_cache", "features/v2_esm2_cache.npz"),
        ("v5_common", "scripts/v5/v5_common.py"),
    ]
    rows = []
    for name, rel in FILES:
        p = C.ROOT / rel
        rows.append({"name": name, "path": rel,
                     "size": p.stat().st_size if p.exists() else 0,
                     "sha256_prefix40": sha(p)[:40] if p.exists() else "MISSING"})
    pd.DataFrame(rows).to_csv(PH0 / "input_manifest.csv", index=False,
                              encoding="utf-8")
    emit(f"- 冻结 **{len(rows)}** 个输入 → `phase0_provenance/input_manifest.csv`")
    emit("")

    # ================= PHASE 1 =================
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    lf0 = pd.read_csv(V70 / "phase1_features" / "block_L_local.csv",
                      low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    s = site[site.pair_row.isin(can.row.values)].copy()
    emit(f"- canonical: **{len(can)}** 行 / **{can.pdb_id.nunique()}** targets")
    emit(f"- V7.0 U_LOCAL: **{len(lf0)}** 行 / **{lf0.pdb_id.nunique()}** targets")
    emit("")

    # ---- 1.2 排除机制审计 ----
    have0 = set(lf0.pair_row.astype(int))
    can["inL"] = can.row.isin(have0)
    g = can.groupby("pdb_id").agg(n_rows=("row", "size"), n_inL=("inL", "sum"))
    ex = g[g.n_inL == 0]
    part = g[(g.n_inL > 0) & (g.n_inL < g.n_rows)]
    emit("### 1.2 排除机制分解")
    emit("")
    emit("| 类别 | targets | 行数 |")
    emit("|---|---:|---:|")
    emit(f"| A 完全排除（0 个已验证位点） | **{len(ex)}** | "
         f"**{int(ex.n_rows.sum())}** |")
    emit(f"| B 部分排除（部分 variant 未通过） | **{len(part)}** | "
         f"**{int((part.n_rows - part.n_inL).sum())}** |")
    emit(f"| 合计排除 | {len(ex) + len(part)} | "
         f"**{int(can[~can.inL].shape[0])}** |")
    emit("")
    emit(f"> **关键**：主要损失来自 **B 类（部分排除）**，"
         f"而非 16 个坏 target（仅 {int(ex.n_rows.sum())} 行）。")
    emit("")

    # 逐 target 审计表
    auc = []
    for pdb in sorted(ex.index.tolist()):
        sub = s[s.pdb_id.astype(str).str.upper() == str(pdb).upper()]
        auc.append({"pdb_id": pdb, "family": str(can.loc[can.pdb_id == pdb,
                                                         "family_label"].iloc[0]),
                    "n_rows": int(ex.loc[pdb, "n_rows"]),
                    "n_sites_records": len(sub),
                    "mapping_status": "|".join(sorted(
                        sub.mapping_status.unique())) if len(sub) else "NONE",
                    "has_pdb_file": (PDBDIR / f"{str(pdb).upper()}.pdb").exists(),
                    "failure_reason": ("NO_SITE_RECORDS" if not len(sub)
                                       else "WT_MISMATCH_OR_CHAIN_MISSING")})
    pd.DataFrame(auc).to_csv(PH1 / "excluded_target_audit.csv", index=False,
                             encoding="utf-8")
    emit(f"- 写出 `excluded_target_audit.csv`（{len(auc)} 个完全排除的 target）")
    emit("")

    # ---- 1.3 救援：允许部分位点 ----
    emit("### 1.3 救援（允许部分位点映射，仅用 WT 校验通过的位点）")
    emit("")
    v22 = pd.read_csv(C.V22, low_memory=False)
    d = C.build_pair_features(verbose=False)
    fwd_mutrow = d["df"]["mut_row"].to_numpy()
    raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                      low_memory=False)
    raw["PDB"] = raw["PDB"].astype(str)
    rk = (raw["PDB"] + "|" + raw["KD(M)"].astype(str) + "|"
          + raw["Mutations"].fillna("WT").astype(str))
    vk = (v22["pdb_id"].astype(str) + "|" + v22["KD(M)"].astype(str) + "|"
          + v22["Mutations"].fillna("WT").astype(str))
    k2v: dict = {}
    for i, k in enumerate(vk):
        k2v.setdefault(k, i)
    rows = can["row"].to_numpy(int)
    mut_v22 = np.array([k2v.get(rk.iloc[r], -1) for r in fwd_mutrow[rows]])
    emit(f"- v2.2 行映射覆盖: {int((mut_v22 >= 0).sum())} / {len(rows)}")
    emit("")

    chain_cache: dict = {}   # 本脚本用（chain_residues）
    v7f_cache: dict = {}     # 传给 v7f.build_local_features（须存 (res, ca_map)）
    recs = []
    stat = {"n_total": 0, "no_pdb": 0, "no_chain": 0, "idx_oob": 0,
            "wt_mismatch": 0, "verified": 0, "no_site_ok": 0}
    per_target = {}
    for pr, gg in s.groupby("pair_row"):
        stat["n_total"] += 1
        pdb = str(gg.pdb_id.iloc[0])
        key = pdb.upper()
        if key not in chain_cache:
            p = PDBDIR / f"{key}.pdb"
            chain_cache[key] = chain_residues(p) if p.exists() else {}
        ch = chain_cache[key]
        if not ch:
            stat["no_pdb"] += 1
            continue
        locs = []
        for r in gg.itertuples():
            c = str(r.ann_chain)
            if c not in ch:
                stat["no_chain"] += 1
                continue
            i = int(r.seq_index0)
            if not (0 <= i < len(ch[c])):
                stat["idx_oob"] += 1
                continue
            resseq, icode, aa = ch[c][i]
            if aa != str(r.wt_aa):
                stat["wt_mismatch"] += 1
                continue
            locs.append({"chain": c, "resseq": resseq, "icode": icode,
                         "wt": str(r.wt_aa), "mut": str(r.mut_aa)})
        if not locs:
            stat["no_site_ok"] += 1
            continue
        stat["verified"] += 1
        got = v7f.build_local_features(
            pdb, [{"pdb_chain": x["chain"], "pdb_resseq": x["resseq"],
                   "icode": x["icode"]} for x in locs], v7f_cache)
        if got is None:
            continue
        ok = [h for h in got if h is not None]
        if not ok:
            continue
        agg = {}
        for kk in set(ok[0].keys()):
            vals = [h[kk] for h in ok if kk in h]
            agg[kk] = float(np.mean(vals))
            agg[f"{kk}__max"] = float(np.max(vals))
        idv = [v7f.mut_identity(x["wt"], x["mut"]) for x in locs
               if x["wt"] in v7f.AA1 and x["mut"] in v7f.AA1]
        if idv:
            for kk in idv[0]:
                agg[kk] = float(np.mean([d0[kk] for d0 in idv]))
        agg["pair_row"] = int(pr)
        agg["pdb_id"] = pdb
        agg["n_sites_verified"] = len(locs)
        agg["n_sites_total"] = len(gg)
        agg["site_coverage_frac"] = len(locs) / max(1, len(gg))
        agg["partial_mapping"] = int(len(locs) < len(gg))
        recs.append(agg)
        per_target[pdb] = per_target.get(pdb, 0) + 1

    emit("| 统计 | 值 |")
    emit("|---|---:|")
    for k, v in stat.items():
        emit(f"| {k} | {v} |")
    emit("")
    LF = pd.DataFrame(recs)
    LF.to_csv(PH1 / "block_L_rescued.csv", index=False, encoding="utf-8")
    emit(f"### 救援后 U_LOCAL")
    emit(f"- **rows: {len(LF)} / {len(can)} = {len(LF) / len(can) * 100:.1f}%**"
         f"（V7.0 为 {len(lf0)} = {len(lf0) / len(can) * 100:.1f}%）")
    emit(f"- **targets: {LF.pdb_id.nunique()} / {can.pdb_id.nunique()}**"
         f"（V7.0 为 {lf0.pdb_id.nunique()}）")
    emit(f"- 其中 `partial_mapping=1` 的行: "
         f"**{int(LF.partial_mapping.sum())}**"
         f"（{LF.partial_mapping.mean() * 100:.1f}%）")
    emit("")
    emit(f"- target 覆盖 **{LF.pdb_id.nunique() / can.pdb_id.nunique() * 100:.1f}%**"
         f"（目标 ≥90% 行覆盖：**{'达成' if len(LF) / len(can) >= 0.90 else '未达成'}**）")
    emit("")

    # ---- 1.5 两个 universe ----
    pd.DataFrame({"pair_row": LF.pair_row.astype(int), "pdb_id": LF.pdb_id,
                  "partial_mapping": LF.partial_mapping}).to_csv(
        PH1 / "u_local_rows.csv", index=False, encoding="utf-8")
    pd.DataFrame({"pair_row": can.row.astype(int), "pdb_id": can.pdb_id}).to_csv(
        PH1 / "u_full_rows.csv", index=False, encoding="utf-8")
    emit("### 1.5 两个 universe")
    emit(f"- **U_FULL** = {len(can)} 行（不需要 Local 的模型用这个）")
    emit(f"- **U_LOCAL** = {len(LF)} 行（涉及 Local 的比较用这个）")
    emit("")

    # ---- 1.6 覆盖偏倚 ----
    emit("### 1.6 Coverage-bias 检验（included vs excluded）")
    emit("")
    T63 = pd.read_csv(C.ROOT / "experiments" / "v6_3" / "phase1_target_table" /
                      "target_mechanism_table.csv", low_memory=False)
    tgt_in = set(LF.pdb_id.astype(str))
    T63["included"] = T63.pdb_id.astype(str).isin(tgt_in)
    bias = []
    emit("| 变量 | included 中位 | excluded 中位 | Mann-Whitney p | 类型 |")
    emit("|---|---:|---:|---:|---|")
    from scipy.stats import mannwhitneyu
    for c, kind in (("n_mutations", "DEPLOY"), ("n_unique_mutation_sites", "DEPLOY"),
                    ("site_entropy", "DEPLOY"),
                    ("alanine_fraction", "DEPLOY"),
                    ("single_mut_fraction", "DEPLOY"),
                    ("true_mean_ddg", "LABEL"), ("ddg_std", "LABEL")):
        if c not in T63.columns:
            continue
        a = pd.to_numeric(T63.loc[T63.included, c], errors="coerce").dropna()
        b = pd.to_numeric(T63.loc[~T63.included, c], errors="coerce").dropna()
        if len(a) < 3 or len(b) < 3:
            emit(f"| `{c}` | — | — | — | 样本不足 |")
            continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = abs(1 - 2 * u / (len(a) * len(b)))
        emit(f"| `{c}` | {a.median():.3f} | {b.median():.3f} | {p:.4g} | "
             f"{kind} (rb={rb:.3f}) |")
        bias.append({"variable": c, "kind": kind, "included_median": a.median(),
                     "excluded_median": b.median(), "p": p,
                     "rank_biserial": rb, "n_included": len(a),
                     "n_excluded": len(b)})
    emit("")
    B = pd.DataFrame(bias)
    B.to_csv(PH1 / "coverage_bias.csv", index=False, encoding="utf-8")
    sig = B[B.p < 0.05] if len(B) else B
    verdict = "COVERAGE_BIAS_PRESENT" if len(sig) > 0 else "NO_CLEAR_COVERAGE_BIAS"
    emit(f"### 判定：**`{verdict}`**"
         f"（{len(sig)}/{len(B)} 个变量 p<0.05）")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "coverage_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH0}/input_manifest.csv")
    print(f"写出 {PH1}/block_L_rescued.csv ({len(LF)} rows, "
          f"{LF.pdb_id.nunique()} targets)")
    print(f"写出 {PH1}/coverage_report.md, coverage_bias.csv, "
          f"excluded_target_audit.csv")
    print(f"\n>>> U_LOCAL 救援: {len(lf0)} -> {len(LF)} 行; {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
