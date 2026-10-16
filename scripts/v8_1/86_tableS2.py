# -*- coding: utf-8 -*-
"""scripts/v8_1/86_tableS2.py -- 新增 Supplementary Table S2（model family / configuration provenance）。

数据源全部真实：
  v8_1/MODEL_FAMILY_REGISTRY.csv        配置级 provenance
  v8_1/all_models_by_rowset.csv         每配置的 universe / 指标
  v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv  weighting audit 是否包含
  v8_1/MATCHED_UNIVERSE_REGISTRY.csv    fingerprint
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
A = Path("v8_1")
OUT = Path("latex_v6")
mf = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv", low_memory=False)
am = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
wa = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv", low_memory=False)

# 主分析：rowset f8c33c7e (3,351 行) 的 10 个 method family
pu = am[am.n_rows == 3351].copy()
wa_keys = set(wa.model_name.astype(str))
wkeys = set(am.key.astype(str)) & wa_keys

rows = []
for fam, g in pu.groupby("method_family"):
    # 与 MODEL_FAMILY_REGISTRY 对齐
    reg = mf[mf.method_family == fam] if "method_family" in mf.columns else mf.iloc[0:0]
    splits = sorted(set(g.split_regime.astype(str)))
    trains = (sorted(set(reg.training_provenance.dropna().astype(str)))
              if "training_provenance" in reg.columns and len(reg) else [])
    reason = ""
    if "representative_reason" in reg.columns and len(reg):
        rr = reg.representative_reason.dropna().astype(str)
        reason = rr.iloc[0] if len(rr) else ""
    rows.append(dict(
        family=fam,
        n_configs=len(g),
        n_seeds=g.key.astype(str).str.extract(r"_(s\d+)$")[0].nunique(),
        ablation=bool(g.is_ablation.any()),
        primary_rep=bool(g.is_primary_representative.any()),
        rowset=g.rowset_hash.iloc[0],
        n_rows=int(g.n_rows.iloc[0]),
        n_ge5=int(g.n_targets_ge5.iloc[0]),
        split=", ".join(splits),
        in_weighting=bool(set(g.key.astype(str)) & wa_keys),
        train_prov="; ".join(trains)[:42] if trains else "-",
        rep_rule=reason[:46] if reason else "performance-blind, by configuration id",
    ))

df = pd.DataFrame(rows).sort_values("family")
df.to_csv(A / "TABLE_S2_MODEL_FAMILIES.csv", index=False, encoding="utf-8")

LATEX = []
for r in df.itertuples():
    LATEX.append(
        f"{r.family.replace('_', chr(92) + '_')} & {r.n_configs} & {r.n_seeds} & "
        f"{'yes' if r.ablation else 'no'} & {r.rowset[:8]} & {r.n_rows:,} & {r.n_ge5} & "
        f"{'yes' if r.in_weighting else 'no'} & "
        f"{r.rep_rule.replace('_', chr(92) + '_')} " + chr(92) * 2)
body = "\n".join(LATEX)

TBL = r"""\begin{landscape}
\begin{table}[h]
  \centering
  \scriptsize
  \begin{threeparttable}
  \caption{Method families in the primary matched universe: composition and provenance
  (Table~S2).}\label{tab:s2families}
  \begin{tabular}{@{}lrrllrrll@{}}
    \toprule
    Method family & Configs & Seeds & Ablation & Fingerprint & Rows & Targets ($n \geq 5$)
      & In weighting audit & Representative-selection rule \\
    \midrule
__S2BODY__
    \bottomrule
  \end{tabular}
  \begin{tablenotes}\footnotesize
    \item The primary matched universe comprises 30 configurations drawn from these 10 method
    families (three seeds per family). A \emph{method family} groups configurations that share
    an architectural or feature-design lineage; configurations within a family are not
    independent observations, so family-level representatives are used wherever independent
    evidence is required. All families share the same row-set fingerprint, which means they
    were evaluated on identical mutation rows and are therefore directly comparable. The
    \emph{weighting audit} column indicates whether the family's configurations are also
    present in the 109-configuration weighting audit; that audit draws on the wider artifact
    registry and is not restricted to the primary universe.
  \end{tablenotes}
  \end{threeparttable}
\end{table}
\end{landscape}
"""
TBL = TBL.replace("__S2BODY__", body)

p = OUT / "supplementary_v6.tex"
s = p.read_text(encoding="utf-8")
anchor = "\\section*{Supplementary Figures}"
if "tab:s2families" not in s and anchor in s:
    s = s.replace(anchor, TBL + "\n" + anchor, 1)
    p.write_text(s, encoding="utf-8")
    print(f"inserted Table S2 with {len(df)} method families")
else:
    print("Table S2 already present or anchor missing")

print(df.to_string(index=False))
