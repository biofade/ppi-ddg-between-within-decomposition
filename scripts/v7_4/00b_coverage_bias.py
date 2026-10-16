# -*- coding: utf-8 -*-
"""scripts/v7_4/00b_coverage_bias.py -- V7.4 PHASE 0（§3）：U_PHYS_STRICT 覆盖偏倚审计.

比较 **included**（进去 U_PHYS_STRICT 的行/target）与 **excluded**（被排除的）。

可部署属性（§3）：
    target 规模、single/multi 比例、突变位点数、丙氨酸比例、序列长度、
    界面距离可得性、family 分布、PDB 链数
标签派生量**仅作诊断**（§3 明确）。

**不**用于选择，只用于确认「结论是否受覆盖偏倚威胁」。
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

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V711 = C.ROOT / "experiments" / "v7_1_1"
V73 = C.ROOT / "experiments" / "v7_3"
PH0 = C.ROOT / "experiments" / "v7_4" / "phase0_confirmation"
PH0.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    t0 = time.time()
    emit("# V7.4 PHASE 0b — U_PHYS_STRICT 覆盖偏倚审计（§3）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    inc = pd.read_csv(C.ROOT / "features" / "v7_3" /
                      "u_phys_strict_manifest.csv", low_memory=False)
    exc = pd.read_csv(V73 / "phase1_universe" / "exclusions.csv", low_memory=False)
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)
    from scipy.stats import mannwhitneyu, chi2_contingency

    inc_rows = set(inc.row_id.astype(int))
    lab["included"] = lab.pair_row.isin(inc_rows)
    emit("## 0b.1 行级覆盖")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| canonical rows | {len(lab)} |")
    emit(f"| **included** | **{int(lab.included.sum())}** |")
    emit(f"| excluded | **{int((~lab.included).sum())}** |")
    emit(f"| included targets | **{lab[lab.included].target_id_canonical.nunique()}** |")
    emit(f"| excluded-only targets | "
         f"**{lab[~lab.included].target_id_canonical.nunique()}** |")
    emit("")

    # ---- 排除原因（行级）----
    emit("### 排除原因（行级，§10）")
    emit("")
    if len(exc):
        vc = exc.reason.value_counts()
        emit("| reason | rows | 占比 |")
        emit("|---|---:|---:|")
        for k, v in vc.items():
            emit(f"| `{k}` | {v} | {v / vc.sum() * 100:.1f}% |")
        emit(f"| **合计** | **{int(vc.sum())}** | 100% |")
    emit("")

    # ---- 逐 target 属性 ----
    emit("## 0b.2 逐 target 属性（included vs excluded-only targets）")
    emit("")
    T63 = pd.read_csv(C.ROOT / "experiments" / "v6_3" / "phase1_target_table" /
                      "target_mechanism_table.csv", low_memory=False)
    T63 = T63.rename(columns={"pdb_id": "target_id_canonical"})
    T63["target_id_canonical"] = T63["target_id_canonical"].astype(str)
    per = lab.groupby("target_id_canonical").agg(
        n_rows=("pair_row", "size"),
        n_incl=("included", "sum"),
        fam=("family_id_canonical", "first")).reset_index()
    per["grp"] = np.where(per.n_incl == per.n_rows, "full_included",
                          np.where(per.n_incl == 0, "fully_excluded", "partial"))
    emit("| 组 | targets | rows |")
    emit("|---|---:|---:|")
    for k, g in per.groupby("grp"):
        emit(f"| {k} | {len(g)} | {int(g.n_rows.sum())} |")
    emit("")

    m = T63.merge(per[["target_id_canonical", "grp", "n_rows"]],
                  on="target_id_canonical", how="left")
    m["grp"] = m["grp"].fillna("not_in_canonical")
    m["included_any"] = m.grp.isin(["full_included", "partial"])
    from scipy.stats import mannwhitneyu as mwu
    emit("| 变量 | included 中位 | excluded 中位 | p | 类型 |")
    emit("|---|---:|---:|---:|---|")
    bias = []
    for c, kind in (("n_mutations", "DEPLOY"),
                    ("n_unique_mutation_sites", "DEPLOY"),
                    ("single_mut_fraction", "DEPLOY"),
                    ("alanine_fraction", "DEPLOY"),
                    ("site_entropy", "DEPLOY"),
                    ("true_mean_ddg", "LABEL"),
                    ("ddg_std", "LABEL")):
        if c not in m.columns:
            continue
        A = pd.to_numeric(m.loc[m.included_any, c], errors="coerce").dropna()
        B = pd.to_numeric(m.loc[~m.included_any, c], errors="coerce").dropna()
        if len(A) < 3 or len(B) < 3:
            emit(f"| `{c}` | — | — | — | 样本不足 |")
            continue
        u, p = mwu(A, B, alternative="two-sided")
        rb = abs(1 - 2 * u / (len(A) * len(B)))
        emit(f"| `{c}` | {A.median():.3f} | {B.median():.3f} | {p:.4g} | "
             f"{kind} (rb={rb:.3f}) |")
        bias.append({"variable": c, "kind": kind, "included_median": A.median(),
                     "excluded_median": B.median(), "p": float(p),
                     "rank_biserial": float(rb),
                     "n_included": len(A), "n_excluded": len(B)})
    emit("")

    # ---- PDB 链数 / 序列长度（Deployable）----
    emit("## 0b.3 结构属性（PDB 链数、序列长度）")
    emit("")
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
    v7f = importlib.util.module_from_spec(sp)
    sys.modules["v7f"] = v7f
    sp.loader.exec_module(v7f)
    T2O = v7f.THREE2ONE
    props = {}
    for t in per.target_id_canonical:
        p = PDBDIR / f"{t.upper()}.pdb"
        nch, maxlen = 0, 0
        if p.exists():
            ch = defaultdict(int)
            for ln in open(p, "r", errors="replace"):
                if ln.startswith("ATOM") and ln[12:16].strip() == "CA" \
                        and ln[17:20].strip().upper() in T2O:
                    ch[ln[21]] += 1
            nch = len(ch); maxlen = max(ch.values()) if ch else 0
        props[t] = (nch, maxlen)
    per["n_pdb_chains"] = per.target_id_canonical.map(
        lambda t: props.get(t, (np.nan, np.nan))[0])
    per["max_chain_len"] = per.target_id_canonical.map(
        lambda t: props.get(t, (np.nan, np.nan))[1])
    emit("| 变量 | included 中位 | excluded 中位 | p | 类型 |")
    emit("|---|---:|---:|---:|---|")
    for c in ("n_pdb_chains", "max_chain_len"):
        A = per.loc[per.n_incl > 0, c].dropna()
        B = per.loc[per.n_incl == 0, c].dropna()
        if len(A) < 3 or len(B) < 3:
            emit(f"| `{c}` | — | — | — | 样本不足 |"); continue
        u, p = mwu(A, B, alternative="two-sided")
        rb = abs(1 - 2 * u / (len(A) * len(B)))
        emit(f"| `{c}` | {A.median():.1f} | {B.median():.1f} | {p:.4g} | "
             f"DEPLOY (rb={rb:.3f}) |")
        bias.append({"variable": c, "kind": "DEPLOY",
                     "included_median": float(A.median()),
                     "excluded_median": float(B.median()), "p": float(p),
                     "rank_biserial": float(rb),
                     "n_included": len(A), "n_excluded": len(B)})
    emit("")

    # ---- family 分布 ----
    emit("## 0b.4 Family 分布")
    emit("")
    fam_inc = set(per.loc[per.n_incl > 0, "fam"].astype(str))
    fam_all = set(per.fam.astype(str))
    emit(f"- included 涉及 family：**{len(fam_inc)}** / 全体 **{len(fam_all)}**")
    emit(f"- **完全无覆盖的 family：{len(fam_all - fam_inc)}** "
         f"{sorted(fam_all - fam_inc)[:10] if fam_all - fam_inc else '→ ✅ 无'}")
    emit("")

    # ---- 界面距离可得性（Deployable）----
    emit("## 0b.5 界面距离可得性")
    emit("")
    BL = pd.read_csv(C.ROOT / "experiments" / "v7_0" / "phase1_features" /
                     "block_L_local.csv", low_memory=False)
    bl_rows = set(BL.pair_row.astype(int))
    lab["has_bl"] = lab.pair_row.isin(bl_rows)
    ct = pd.crosstab(lab.included, lab.has_bl)
    emit("| | 有 BLOCK L（可算环境） | 无 BLOCK L |")
    emit("|---|---:|---:|")
    for k in (True, False):
        if k in ct.index:
            emit(f"| {'included' if k else 'excluded'} | {ct.loc[k].get(True, 0)} | "
                 f"{ct.loc[k].get(False, 0)} |")
    if ct.shape == (2, 2):
        chi2, p, _, _ = chi2_contingency(ct)
        emit("")
        emit(f"- χ²={chi2:.2f}, **p={p:.4g}** → "
             f"{'环境可得性与纳入**相关**' if p < 0.05 else '无显著相关'}")
    emit("")

    # ---- 判定 ----
    B = pd.DataFrame(bias)
    B.to_csv(PH0 / "coverage_bias.csv", index=False, encoding="utf-8")
    dep = B[B.kind == "DEPLOY"] if len(B) else B
    sig = dep[dep.p < 0.05] if len(dep) else dep
    verdict = ("COVERAGE_BIAS_PRESENT" if len(sig) > 0
               else "NO_CLEAR_COVERAGE_BIAS")
    emit("## 0b.6 判定")
    emit("")
    emit(f"- DEPLOYABLE 变量中 p<0.05 的个数：**{len(sig)} / {len(dep)}**")
    emit(f"- **`{verdict}`**")
    emit("")
    if len(B):
        lab_only = B[B.kind == "LABEL"]
        if len(lab_only):
            emit("> 标签派生量（**仅诊断，不用于选择**）：")
            for _, r in lab_only.iterrows():
                emit(f">   - `{r.variable}`: included 中位 {r.included_median:.3f} "
                     f"vs excluded {r.excluded_median:.3f}（p={r.p:.4g}）")
        emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH0 / "coverage_bias_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
    print(f"\n写出 coverage_bias.csv 与 coverage_bias_report.md")
    print(f">>> {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
