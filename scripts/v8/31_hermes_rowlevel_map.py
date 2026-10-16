# -*- coding: utf-8 -*-
"""scripts/v8/31_hermes_rowlevel_map.py -- HERMES 官方行级预测 vs V8 canonical 行映射审计。

HERMES 在 `hermes_experiments` 仓库**直接发布了逐突变 ΔΔG 预测 CSV**（§9 reuse-first
的首选：官方 row-level 预测，无需推理、无需 GPU）。

本脚本：
  1. 下载官方 CSV 到 `F:\\dsh\\downloads\\hermes`
  2. 摸清列名与键（`PDB_filename` / `mutant_chain` / `mutant` / `Mutation(s)_cleaned`）
  3. 与 V8 canonical（4473 行, key=`row_id`/`pair_row`）做**可行映射率**评估
  4. 给出 evidence_tier 判定（§11）+ §50 污染状态

**只做映射审计**，不进入记分板（映射率决定能否进 Tier 1/2）。
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
DL = Path("F:/dsh/downloads/hermes")
DL.mkdir(parents=True, exist_ok=True)
OUT = R / "experiments" / "v8" / "external"
OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0"}
BASE = "https://raw.githubusercontent.com/StatPhysBio/hermes_experiments/main"

FILES = {
    "hermes_py000_cdna117k":
        "skempi/results/hermes_py_000_ft_cdna117k_ddg_st/"
        "skempi_v2_cleaned_NO_1KBH-hermes_py_000_ft_cdna117k_ddg_st-use_mt_structure=0.csv",
    "hermes_py050_cdna117k":
        "skempi/results/hermes_py_050_ft_cdna117k_ddg_st/"
        "skempi_v2_cleaned_NO_1KBH-hermes_py_050_ft_cdna117k_ddg_st-use_mt_structure=0.csv",
    "hermes_skempi_targets":
        "skempi/skempi_v2_cleaned_NO_1KBH.csv",
    "hermes_train_test_csv":
        "skempi/skempi_v2_cleaned.csv",
}

local = {}
for name, rel in FILES.items():
    dst = DL / f"{name}.csv"
    if dst.exists() and dst.stat().st_size > 1000:
        print(f"  [cached] {name} ({dst.stat().st_size:,} B)")
        local[name] = dst
        continue
    url = f"{BASE}/{rel}"
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=90) as r:
            data = r.read()
        dst.write_bytes(data)
        print(f"  [downloaded] {name} ({len(data):,} B)  <- {rel}")
        local[name] = dst
    except Exception as e:  # noqa: BLE001
        print(f"  [FAIL] {name}: {type(e).__name__} {e}")

print("\n" + "=" * 74)
print("1) HERMES 官方预测文件的列与样本")
print("=" * 74)
H = {}
for name in ("hermes_py000_cdna117k", "hermes_py050_cdna117k"):
    if name not in local:
        continue
    d = pd.read_csv(local[name], low_memory=False)
    H[name] = d
    print(f"\n[{name}] shape={d.shape}")
    print(f"  列({len(d.columns)}): {list(d.columns)[:16]}")
    print(d.head(5).to_string(index=False)[:2200])

print("\n" + "=" * 74)
print("2) 与 V8 canonical 的行映射可行性")
print("=" * 74)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
print(f"canonical: {len(lab)} 行 / {lab.target_id_canonical.nunique()} targets")
print(f"canonical 列: {list(lab.columns)[:20]}")

if "hermes_py000_cdna117k" in H:
    d = H["hermes_py000_cdna117k"]
    # 找 PDB 列与突变列
    pdb_col = next((c for c in d.columns
                    if c.lower() in ("pdb_filename", "pdb", "pdb_id", "pdbid")), None)
    mut_col = next((c for c in d.columns
                    if "mut" in c.lower() and "chain" not in c.lower()), None)
    chn_col = next((c for c in d.columns if "chain" in c.lower()), None)
    ddg_col = next((c for c in d.columns
                    if "ddg" in c.lower() or "neg_ddg" in c.lower()), None)
    print(f"\n  推定列: pdb={pdb_col!r} mut={mut_col!r} chain={chn_col!r} "
          f"ddg={ddg_col!r}")
    if pdb_col:
        hp = d[pdb_col].astype(str).str.upper().str.strip()
        # 取 PDB 前 4 位（HERMES 可能带 .pdb 或 _A_B 后缀）
        hp4 = hp.str.extract(r"^([0-9A-Z]{4})")[0]
        cp4 = lab.target_id_canonical.astype(str).str.upper().str.strip().str[:4]
        hs, cs = set(hp4.dropna()), set(cp4.dropna())
        inter = hs & cs
        print(f"\n  HERMES 唯一 PDB(前4位) = {len(hs)}")
        print(f"  canonical 唯一 target  = {len(cs)}")
        print(f"  交集 = {len(inter)}")
        print(f"  canonical 覆盖率 = {len(inter)/max(1,len(cs)):.4f}")
        # 行级：按 (PDB4) 粗匹配
        n_h = len(d)
        print(f"\n  HERMES 行数 = {n_h}")
        print(f"  canonical 行数 = {len(lab)}")
        # 逐行：canonical 的 (pdb4) 是否出现在 HERMES
        lab_hit = cp4.isin(hs)
        print(f"  canonical 行中 PDB4 能在 HERMES 找到的 = "
              f"{int(lab_hit.sum())} / {len(lab)} = {lab_hit.mean():.4f}")
    if ddg_col:
        v = pd.to_numeric(d[ddg_col], errors="coerce")
        print(f"\n  HERMES {ddg_col}: n={v.notna().sum()} "
              f"mean={v.mean():.4f} sd={v.std():.4f} "
              f"min={v.min():.3f} max={v.max():.3f}")
        print(f"  ⇒ V8 约定 ΔΔG>0 = 结合**变弱**；HERMES 此列为 "
              f"{'neg_ddg_pred（可能方向相反，需审计）' if 'neg' in ddg_col.lower() else 'ddg_pred'}")

print("\n" + "=" * 74)
print("3) HERMES 发布的 split（easy/medium/hard fold1-3）")
print("=" * 74)
for sub in ("easy", "medium", "hard"):
    for fold in (1, 2, 3):
        u = (f"{BASE}/skempi/supervised_models/{sub}/fold{fold}"
             f"/skempi_v2_cleaned_{sub}_fold{fold}.csv")
        try:
            req = urllib.request.Request(u, headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                b = r.read()
            print(f"  OK  {sub}/fold{fold}  {len(b):,} B")
            (DL / f"split_{sub}_fold{fold}.csv").write_bytes(b)
        except Exception as e:  # noqa: BLE001
            print(f"  --  {sub}/fold{fold}  ({type(e).__name__})")

meta = {
    "source": "https://github.com/StatPhysBio/hermes_experiments",
    "official_row_level_predictions": sorted(list(H.keys())),
    "note": ("HERMES 是 §9 的首选资源类型：**官方发布的 row-level 预测**，"
             "无需推理、无需 GPU、无需 checkpoint。"),
    "evidence_tier_pending": ("须先确认 (a) 行映射率, (b) 其 SKEMPI 微调是否与 "
                              "OUTER_FOLDS_V7_2 家族重叠（§50）"),
    "ddg_sign_convention": "见脚本第 2 节输出；neg_ddg_pred 可能与本项目相反",
}
(OUT / "hermes_audit.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 {OUT / 'hermes_audit.json'}")
