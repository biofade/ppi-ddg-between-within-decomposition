# -*- coding: utf-8 -*-
"""scripts/v6_3/01_build_target_table.py -- V6.3 PHASE 1: target-level 机制表.

每个 target（pdb_id）一行，主键 `pdb_id`（含 `family` 与 `L4_fold`）。

**最高优先级原则（V6.3 §1.1）**：严格区分
* `LABEL_DERIVED = YES`：用到 true ΔΔG 的字段 → **只可用于机制解释**，
  **不可**用于未来可部署的 calibration 模型
* `DEPLOYABLE = YES`：预测前可知 → 可用于未来模型

所有字段写入 `feature_dictionary.csv` 并标注归属。

输出
----
* `phase1_target_table/target_mechanism_table.csv`
* `phase1_target_table/feature_dictionary.csv`
* `phase1_target_table/missingness_report.md`
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V621 = C.ROOT / "experiments" / "v6_2_1"
PH1 = C.ROOT / "experiments" / "v6_3" / "phase1_target_table"
PH1.mkdir(parents=True, exist_ok=True)
OUT_T = PH1 / "target_mechanism_table.csv"
OUT_D = PH1 / "feature_dictionary.csv"
OUT_M = PH1 / "missingness_report.md"
L: list[str] = []
DICT: list[dict] = []

# 三字母 → 一字母
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O"}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def shannon(counts) -> float:
    """Shannon entropy（bit）。counts 为非负计数。"""
    c = np.asarray([x for x in counts if x > 0], float)
    if c.sum() <= 0:
        return 0.0
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def gini(vals) -> float:
    """Gini 系数（0 = 完全均匀；→1 = 高度集中）。"""
    x = np.sort(np.asarray([v for v in vals if v >= 0], float))
    n = len(x)
    if n == 0 or x.sum() <= 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2 * (idx * x).sum()) / (n * x.sum()) - (n + 1) / n)


# 侧链性质分类（用于 composition）
CHARGED = set("DEKRH")
POLAR = set("STNQYC")
HYDROPHOBIC = set("AVLIMFWPG")
SPECIAL = set("PG")


def parse_mutation(s: str):
    """解析 mutation 字符串，返回 (wt, mut, pos, chain) 列表。

    支持形态（本项目 ddg_pairs 实际使用）：
    * `A_P2A`      —— 单字母三字符码 `P2A`
    * `AB_Y32A`
    * `A_P2A;B_Q5R`（多突变，分号/逗号分隔）
    * 三字母码 `ALA2GLY`（兼容保留）
    """
    out = []
    if not isinstance(s, str) or not s.strip():
        return out
    for tok in re.split(r"[;,]", s.strip()):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^([A-Za-z0-9]+)_(.+)$", tok)
        chain, body = (m.group(1), m.group(2)) if m else ("", tok)
        body = body.strip().upper()
        # 形态 1：单字母码  X<pos>Y   （如 P2A）
        for mm in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])", body):
            out.append((mm.group(1), mm.group(3), int(mm.group(2)), chain))
        if out:
            continue
        # 形态 2：三字母码  ALA2GLY
        for mm in re.finditer(r"(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|"
                              r"LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|MSE)"
                              r"(\d+)"
                              r"(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|"
                              r"LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|MSE)", body):
            out.append((AA3.get(mm.group(1), "X"), AA3.get(mm.group(3), "X"),
                        int(mm.group(2)), chain))
    return out


def main() -> int:
    t0 = time.time()
    emit("# V6.3 PHASE 1 — Target-level 机制表")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 主键: `pdb_id`（= target）；另含 `family`")
    emit("- **严格区分 `LABEL_DERIVED` 与 `DEPLOYABLE`**（V6.3 §1.1）")
    emit("")

    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    pairs = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = pairs[pairs["direction"] == "forward"].reset_index(drop=True)
    # 与 canonical 对齐
    rows = can["row"].to_numpy(int)
    fwd2 = fwd.iloc[rows].reset_index(drop=True)
    v22 = pd.read_csv(C.V22, low_memory=False)

    emit(f"- matched set: **{len(can)}** 行 / **{can.pdb_id.nunique()}** targets")
    emit(f"- ddg_pairs forward 对齐: **{len(fwd2)}** 行")
    emit("")

    # ---------- L4 OOF ensemble 预测 ----------
    z6 = np.load(V62P2 / "split4_predictions.npz", allow_pickle=True)
    z2 = np.load(V621 / "phase3_d_only" / "d_only_predictions.npz",
                 allow_pickle=True)
    lv = "L4_family"
    y = z6[f"{lv}_y"].astype(float)
    g = np.array([str(x) for x in z6[f"{lv}_pdb"]])
    ks = [f"{lv}_M1_s{s}" for s in (42, 43, 44) if f"{lv}_M1_s{s}" in z6.files]
    pM1 = np.mean([z6[k] for k in ks], axis=0).astype(float)
    ks0 = [f"{lv}_M0_s{s}" for s in (42, 43, 44) if f"{lv}_M0_s{s}" in z6.files]
    pM0 = np.mean([z6[k] for k in ks0], axis=0).astype(float)
    ks2 = [f"{lv}_M2_s{s}" for s in (42, 43, 44) if f"{lv}_M2_s{s}" in z2.files]
    pM2 = np.mean([z2[k] for k in ks2], axis=0).astype(float) if ks2 else None
    m = np.isfinite(y) & np.isfinite(pM1) & np.isfinite(pM0)
    emit(f"- L4 预测可用行: **{int(m.sum())}**；M1 ensemble 来自 {len(ks)} seeds")
    emit("")

    # ---------- 逐 target 聚合 ----------
    # 预计算每行的 mutation 解析结果
    parsed = [parse_mutation(s) for s in fwd2["mutation"].astype(str)]
    # 每行的 v2.2 索引（用于结构/实验字段）
    v22_idx_wt = None
    try:
        raw = pd.read_csv(C.DATA / "ppb_affinity_proteinea" / "filtered.csv",
                          low_memory=False)
        raw["PDB"] = raw["PDB"].astype(str)
        rk = (raw["PDB"] + "|" + raw["KD(M)"].astype(str) + "|"
              + raw["Mutations"].fillna("WT").astype(str))
        vk = (v22["pdb_id"].astype(str) + "|" + v22["KD(M)"].astype(str) + "|"
              + v22["Mutations"].fillna("WT").astype(str))
        k2v: dict[str, int] = {}
        for i, k in enumerate(vk):
            k2v.setdefault(k, i)
        v22_idx_mut = np.array([k2v.get(rk.iloc[r], -1)
                                for r in fwd2["mut_row"]])
    except Exception as exc:  # noqa: BLE001
        emit(f"- ⚠️ v2.2 行映射失败（{type(exc).__name__}）→ 结构字段将用 pdb_id 聚合")
        v22_idx_mut = np.full(len(fwd2), -1)

    have_v22 = (v22_idx_mut >= 0)
    emit(f"- 能映射到 v2.2 行的 mutation: **{int(have_v22.sum())} / {len(fwd2)}**"
         f"（{have_v22.mean() * 100:.1f}%）")
    # 结构字段（每 target 取该 target 内突变的 v2.2 行均值）
    # 只保留**真正数值型**的结构/界面描述符；排除字符串列与 label 列，
    # 否则会产生一堆全 NaN 的 `struct_<字符串列>`（本脚本 v1 的 bug）。
    EXCLUDE = {"real_dG", "KD(M)", "n_missing_chain", "chain_map_ok", "split",
               "structure_usable", "row_structure_usable", "n_features_present",
               "pdb_id", "interface_key"}
    STRUCT_COLS = [c for c in v22.columns
                   if c not in EXCLUDE
                   and pd.api.types.is_numeric_dtype(v22[c])]
    emit(f"- 可用结构描述符（数值型）: **{len(STRUCT_COLS)}** 个")
    emit("")

    recs = []
    for pdb, sub in fwd2.groupby("pdb_id", sort=True):
        idx = sub.index.to_numpy()
        sel = m[idx]
        if sel.sum() < 2:
            continue
        yy = y[idx][sel]
        pm1 = pM1[idx][sel]
        pm0 = pM0[idx][sel]
        # ---- 1.1 calibration outcomes ----
        r = {
            "pdb_id": str(pdb),
            "family": str(can.loc[can.pdb_id == pdb, "family_label"].iloc[0])
            if (can.pdb_id == pdb).any() else "NA",
            "n_rows": int(sel.sum()),
        }
        r["true_mean_ddg"] = float(yy.mean())
        r["pred_mean_ddg"] = float(pm1.mean())
        r["signed_mean_error"] = r["pred_mean_ddg"] - r["true_mean_ddg"]
        r["abs_mean_error"] = abs(r["signed_mean_error"])
        r["squared_mean_error"] = r["signed_mean_error"] ** 2
        r["pred_mean_M0"] = float(pm0.mean())
        r["signed_mean_error_M0"] = r["pred_mean_M0"] - r["true_mean_ddg"]
        # within metrics
        from scipy.stats import spearmanr, kendalltau
        if len(yy) >= 5:
            rho = spearmanr(yy, pm1).statistic
            r["within_spearman"] = float(rho) if rho == rho else np.nan
            try:
                r["within_kendall"] = float(kendalltau(yy, pm1).statistic)
            except Exception:  # noqa: BLE001
                r["within_kendall"] = np.nan
        else:
            r["within_spearman"] = np.nan
            r["within_kendall"] = np.nan

        # ---- 1.2 mutation-library size ----
        muts = [p for i in idx for p in parsed[i]]
        sites = {(p[3], p[2]) for p in muts}          # (chain, pos)
        r["n_mutations"] = int(len(idx))              # = 该 target 的变体数
        r["log_n_mutations"] = float(np.log1p(len(idx)))
        r["n_unique_mutation_sites"] = int(len(sites))
        r["mutations_per_site"] = (len(muts) / len(sites)) if sites else np.nan
        per_site = pd.Series([(p[3], p[2]) for p in muts]).value_counts()
        r["max_mutations_per_site"] = int(per_site.max()) if len(per_site) else 0

        # ---- 1.4 mutation composition ----
        wt_res = [p[0] for p in muts]
        mu_res = [p[1] for p in muts]
        n_ok = len(muts)
        if n_ok:
            r["alanine_fraction"] = float(np.mean([x == "A" for x in mu_res]))
            r["glycine_fraction"] = float(np.mean([x == "G" for x in mu_res]))
            r["charged_to_neutral_fraction"] = float(np.mean(
                [(w in CHARGED) and (mm not in CHARGED) for w, mm in zip(wt_res, mu_res)]))
            r["hydrophobic_to_polar_fraction"] = float(np.mean(
                [(w in HYDROPHOBIC) and (mm in POLAR) for w, mm in zip(wt_res, mu_res)]))
            # conservative：同类别内互换
            def cls(a):
                return ("chg" if a in CHARGED else
                        "pol" if a in POLAR else
                        "hyd" if a in HYDROPHOBIC else "oth")
            r["conservative_fraction"] = float(np.mean(
                [cls(w) == cls(mm) for w, mm in zip(wt_res, mu_res)]))
            r["nonconservative_fraction"] = 1.0 - r["conservative_fraction"]
            r["WT_residue_entropy"] = shannon(pd.Series(wt_res).value_counts())
            r["Mut_residue_entropy"] = shannon(pd.Series(mu_res).value_counts())
            pairs_vc = pd.Series([f"{w}{mm}" for w, mm in zip(wt_res, mu_res)]).value_counts()
            r["mutation_pair_entropy"] = shannon(pairs_vc)
        else:
            for k in ("alanine_fraction", "glycine_fraction",
                      "charged_to_neutral_fraction",
                      "hydrophobic_to_polar_fraction", "conservative_fraction",
                      "nonconservative_fraction", "WT_residue_entropy",
                      "Mut_residue_entropy", "mutation_pair_entropy"):
                r[k] = np.nan

        # ---- 1.5 single / multiple mutation composition ----
        nm = pd.to_numeric(sub["n_mut"], errors="coerce").to_numpy(float)
        r["single_mut_fraction"] = float(np.mean(nm == 1))
        r["multi_mut_fraction"] = float(np.mean(nm > 1))
        r["mean_mutations_per_variant"] = float(np.nanmean(nm))
        r["max_mutations_per_variant"] = float(np.nanmax(nm)) if len(nm) else np.nan

        # ---- 1.6 site diversity ----
        sc = pd.Series([f"{p[3]}_{p[2]}" for p in muts]).value_counts()
        r["site_entropy"] = shannon(sc)
        r["mutation_count_gini"] = gini(nm)
        r["site_count_gini"] = gini(sc.to_numpy())

        # ---- 1.3 coverage + 1.7 structural（取 target 内 v2.2 行均值）----
        vi = v22_idx_mut[idx]
        vi = vi[vi >= 0]
        if len(vi):
            for c in STRUCT_COLS:
                try:
                    r[f"struct_{c}"] = float(
                        pd.to_numeric(v22[c].iloc[vi], errors="coerce").mean())
                except Exception:  # noqa: BLE001
                    r[f"struct_{c}"] = np.nan
            prot_len = float(pd.to_numeric(v22["n_res_lig"].iloc[vi],
                                           errors="coerce").mean()
                             + pd.to_numeric(v22["n_res_rec"].iloc[vi],
                                             errors="coerce").mean())
            iface_cnt = float(pd.to_numeric(v22["n_iface_res_lig"].iloc[vi],
                                            errors="coerce").mean()
                              + pd.to_numeric(v22["n_iface_res_rec"].iloc[vi],
                                              errors="coerce").mean())
            r["protein_length"] = prot_len
            r["interface_residue_count"] = iface_cnt
            r["site_coverage_all"] = (len(sites) / prot_len) if prot_len > 0 else np.nan
            r["site_coverage_interface"] = ((len(sites) / iface_cnt)
                                            if iface_cnt > 0 else np.nan)
            # 1.8 sequence
            r["lig_seq_len"] = float(v22["Ligand Sequences"].iloc[vi]
                                     .astype(str).str.len().mean())
            r["rec_seq_len"] = float(v22["Receptor Sequences"].iloc[vi]
                                     .astype(str).str.len().mean())
            r["lig_seq_nuniq"] = int(v22["Ligand Sequences"].iloc[vi].nunique())
            r["rec_seq_nuniq"] = int(v22["Receptor Sequences"].iloc[vi].nunique())
        else:
            for c in ("protein_length", "interface_residue_count",
                      "site_coverage_all", "site_coverage_interface",
                      "lig_seq_len", "rec_seq_len", "lig_seq_nuniq",
                      "rec_seq_nuniq"):
                r[c] = np.nan

        # ---- 1.9 experimental / source metadata ----
        if len(vi):
            src = v22["Source Data Set"].iloc[vi]
            r["source_nuniq"] = int(src.nunique())
            r["source_top"] = str(src.mode().iloc[0]) if len(src.mode()) else "NA"
            am = v22["Affinity Method"].iloc[vi]
            r["assay_nuniq"] = int(am.nunique())
            r["assay_top"] = str(am.mode().iloc[0]) if len(am.mode()) else "NA"
            tk = pd.to_numeric(
                v22["Temperature(K)"].iloc[vi].astype(str).str.extract(
                    r"(\d+)")[0], errors="coerce")
            r["temp_mean"] = float(tk.mean()) if tk.notna().any() else np.nan
            r["temp_missing_frac"] = float(tk.isna().mean())
        else:
            for c in ("source_nuniq", "assay_nuniq", "temp_mean",
                      "temp_missing_frac"):
                r[c] = np.nan
            r["source_top"] = "NA"
            r["assay_top"] = "NA"

        # ---- 1.10 label-derived ----
        r["ddg_std"] = float(yy.std(ddof=0))
        r["ddg_range"] = float(yy.max() - yy.min())
        r["ddg_IQR"] = float(np.percentile(yy, 75) - np.percentile(yy, 25))
        r["beneficial_fraction"] = float((yy < 0).mean())
        r["deleterious_fraction"] = float((yy > 0).mean())
        r["extreme_fraction"] = float((np.abs(yy) > 2.0).mean())
        recs.append(r)

    T = pd.DataFrame(recs)
    emit(f"## 结果：**{len(T)}** 个 target（原始 {can.pdb_id.nunique()}）")
    emit("")
    T.to_csv(OUT_T, index=False, encoding="utf-8")

    # ---------- feature dictionary ----------
    LABEL_DERIVED = {"true_mean_ddg", "signed_mean_error", "abs_mean_error",
                     "squared_mean_error", "ddg_std", "ddg_range", "ddg_IQR",
                     "beneficial_fraction", "deleterious_fraction",
                     "extreme_fraction", "within_spearman", "within_kendall"}
    for c in T.columns:
        kind = ("LABEL_DERIVED" if c in LABEL_DERIVED else
                "ID" if c in ("pdb_id", "family") else
                "DEPLOYABLE")
        DICT.append({"column": c, "dtype": str(T[c].dtype),
                     "kind": kind, "deployable": kind == "DEPLOYABLE",
                     "label_derived": kind == "LABEL_DERIVED",
                     "non_null_pct": round(float(T[c].notna().mean() * 100), 2),
                     "n_unique": int(T[c].nunique()),
                     "note": ("含 true ΔΔG 信息 → 仅可用于机制解释"
                              if kind == "LABEL_DERIVED" else "")})
    pd.DataFrame(DICT).to_csv(OUT_D, index=False, encoding="utf-8")

    # ---------- missingness ----------
    emit("## 缺失率（非空比例 < 90% 的字段）")
    emit("")
    emit("| 字段 | 非空 % | 类型 |")
    emit("|---|---:|---|")
    for c in T.columns:
        nn = T[c].notna().mean() * 100
        if nn < 90:
            emit(f"| `{c}` | {nn:.1f} | "
                 f"{'LABEL_DERIVED' if c in LABEL_DERIVED else 'DEPLOYABLE'} |")
    emit("")
    n_unavail = [c for c in T.columns if T[c].notna().sum() == 0]
    emit(f"- **完全不可用（UNAVAILABLE）字段: {len(n_unavail)}** "
         f"{n_unavail if n_unavail else ''}")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_M.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_T}  ({len(T)} targets, {T.shape[1]} cols)")
    print(f"写出 {OUT_D}")
    print(f"写出 {OUT_M}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
