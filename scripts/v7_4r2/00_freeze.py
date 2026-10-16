# -*- coding: utf-8 -*-
"""scripts/v7_4r2/00_freeze.py -- V7.4R2 上游冻结校验（只读）。

    python scripts/v7_4r2/00_freeze.py            # 校验
    python scripts/v7_4r2/00_freeze.py --write    # 重新生成
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import sys
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
PH0 = R / "experiments" / "v7_4r2" / "phase0_freeze"
PH0.mkdir(parents=True, exist_ok=True)

FILES = [
    ("canonical_target_labels", "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv"),
    ("outer_folds_v7_2_protocol", "experiments/v7_1_1/phase3_nested_cv/protocol.json"),
    ("outer_folds_v7_2_manifest", "experiments/v7_1_1/phase3_nested_cv/fold_manifest.csv"),
    ("v7_2_evaluation_contract", "experiments/v7_1_1/V7_2_EVALUATION_CONTRACT.md"),
    ("v73_u_phys_strict_manifest", "features/v7_3/u_phys_strict_manifest.csv"),
    ("v73_physics_features", "features/v7_3/physics_features.parquet"),
    ("v73_family_spec", "experiments/v7_3/phase5_rankers/family_spec.json"),
    ("v73_family_predictions", "experiments/v7_3/phase5_rankers/family_predictions.npz"),
    ("v74_u_relax_strict", "features/v7_4/u_relax_strict_manifest.csv"),
    ("v74r_variant_index", "evoef2/v7_4r/evoef2_variant_index.csv"),
    ("mutation_sites", "data/v5_0/mutation_sites.csv"),
    ("matched_analysis_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
    ("v5_common", "scripts/v5/v5_common.py"),
    ("v6_2_split4", "scripts/v6_2/02_split4_ddg.py"),
    ("v6_1_fastmetrics", "scripts/v6_1/fastmetrics.py"),
    ("v7_0_features", "scripts/v7_0/01_features.py"),
    ("r2_canonical_wt_strict", "features/v7_4r2/canonical_wt_strict_manifest.csv"),
]


def sha(p: Path, ch: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def compute() -> pd.DataFrame:
    rows = []
    for nm, rel in FILES:
        p = R / rel
        rows.append({"name": nm, "path": rel,
                     "size": p.stat().st_size if p.exists() else 0,
                     "sha256_prefix40": sha(p)[:40] if p.exists() else "MISSING"})
    pdbs = sorted(glob.glob(str(R / "data" / "pdb" / "*.pdb")))
    rows.append({"name": "pdb_dir(aggregate)", "path": "data/pdb/*.pdb",
                 "size": len(pdbs),
                 "sha256_prefix40": f"{len(pdbs)}_"
                                    f"{sum(Path(x).stat().st_size for x in pdbs)}"})
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    out = PH0 / "input_hashes.csv"
    cur = compute()
    if a.write or not out.exists():
        cur.to_csv(out, index=False, encoding="utf-8")
        print(f"写出 {out} ({len(cur)} 项)")
        return 0
    ref = pd.read_csv(out).set_index("name")["sha256_prefix40"]
    bad = []
    for _, r in cur.iterrows():
        if r["name"] not in ref.index:
            bad.append((r["name"], "NEW", r["sha256_prefix40"])); continue
        if ref[r["name"]] == "MISSING":
            continue
        if ref[r["name"]] != r["sha256_prefix40"]:
            bad.append((r["name"], ref[r["name"]], r["sha256_prefix40"]))
    print(f"校验 {len(cur)} 项")
    for b in bad:
        print("  CHANGED:", b)
    print(f"UPSTREAM_HASH_CHANGE = {len(bad)}{'  ✅' if not bad else '  ❌'}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
