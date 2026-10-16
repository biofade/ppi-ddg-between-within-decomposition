# -*- coding: utf-8 -*-
"""scripts/v8/00_freeze.py -- V8 PHASE 0：冻结上游 + 制品清单 + 协议/溯源（§8）。

只读。生成：
    experiments/v8/phase0_freeze/input_hashes.csv
    experiments/v8/phase0_freeze/artifact_inventory.csv
    experiments/v8/phase0_freeze/protocol.json
    experiments/v8/phase0_freeze/provenance.md
要求 UPSTREAM_HASH_CHANGE = 0。
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
PH0 = R / "experiments" / "v8" / "phase0_freeze"
PH0.mkdir(parents=True, exist_ok=True)

# ---- 权威上游（只读）----
FILES = [
    ("canonical_target_labels", "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv"),
    ("outer_folds_v7_2_protocol", "experiments/v7_1_1/phase3_nested_cv/protocol.json"),
    ("outer_folds_v7_2_manifest", "experiments/v7_1_1/phase3_nested_cv/fold_manifest.csv"),
    ("v7_2_evaluation_contract", "experiments/v7_1_1/V7_2_EVALUATION_CONTRACT.md"),
    ("v6_2_matched_analysis_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
    ("mutation_sites", "data/v5_0/mutation_sites.csv"),
    ("ddg_pairs", "data/v4.0/ddg_pairs.csv"),
    # V7.3 冻结结果
    ("v73_u_phys_strict_manifest", "features/v7_3/u_phys_strict_manifest.csv"),
    ("v73_physics_features", "features/v7_3/physics_features.parquet"),
    ("v73_family_spec", "experiments/v7_3/phase5_rankers/family_spec.json"),
    ("v73_family_predictions", "experiments/v7_3/phase5_rankers/family_predictions.npz"),
    ("v73_family_ensemble", "experiments/v7_3/phase5_rankers/family_ensemble.csv"),
    ("v73_family_results", "experiments/v7_3/phase5_rankers/family_results.csv"),
    ("v73_bootstrap", "experiments/v7_3/phase6_bootstrap/contrasts.csv"),
    # V7.4R / V7.4（含 5000 bootstrap 确认）
    ("v74_u_relax_strict", "features/v7_4/u_relax_strict_manifest.csv"),
    ("v74_confirm_bootstrap5000", "experiments/v7_4/phase0_confirmation/v73_bootstrap_5000.csv"),
    ("v74_confirm_gate", "experiments/v7_4/phase0_confirmation/confirm_gate.json"),
    # V7.4R2（已暂停，但其审计产物是 V8 的结构口径依据）
    ("v74r2_canonical_wt_strict", "features/v7_4r2/canonical_wt_strict_manifest.csv"),
    ("v74r2_audit_meta", "experiments/v7_4r2/phase1_audit/audit_meta.json"),
    # 代码（口径依赖）
    ("v5_common", "scripts/v5/v5_common.py"),
    ("v6_2_split4", "scripts/v6_2/02_split4_ddg.py"),
    ("v6_1_fastmetrics", "scripts/v6_1/fastmetrics.py"),
    ("v7_0_features", "scripts/v7_0/01_features.py"),
]

# ---- 制品清单扫描根 ----
SCAN_GLOBS = [
    "experiments/*/phase*/**/*.csv",
    "experiments/*/phase*/**/*.npz",
    "experiments/*/*.csv",
    "features/*/*.npz",
    "features/*/*.parquet",
    "logs/*/**/*.npz",
]
# 预测类文件名特征
PRED_TOKENS = ("predict", "pred", "ranking", "results", "ensemble", "contrasts",
               "metrics", "family_results", "seedmean", "bootstrap", "noise")


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
                 "sha256_prefix40": f"{len(pdbs)}_{sum(Path(x).stat().st_size for x in pdbs)}"})
    return pd.DataFrame(rows)


def inventory() -> pd.DataFrame:
    seen, rows = set(), []
    for g in SCAN_GLOBS:
        for f in glob.glob(str(R / g), recursive=True):
            p = Path(f)
            if p.name == "input_hashes.csv" or str(p) in seen:
                continue
            seen.add(str(p))
            low = p.name.lower()
            kind = "prediction_like" if any(t in low for t in PRED_TOKENS) else "other"
            try:
                sz = p.stat().st_size
            except OSError:
                continue
            rows.append({"path": str(p.relative_to(R)).replace("\\", "/"),
                         "ext": p.suffix.lower(), "size_bytes": sz,
                         "kind": kind,
                         "phase": p.parent.name})
    return pd.DataFrame(rows).sort_values(["kind", "size_bytes"],
                                          ascending=[True, False])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    out = PH0 / "input_hashes.csv"
    cur = compute()
    rc = 0
    if a.write or not out.exists():
        cur.to_csv(out, index=False, encoding="utf-8")
        print(f"写出 {out} ({len(cur)} 项)")
    else:
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
        rc = 0 if not bad else 1

    inv = inventory()
    inv.to_csv(PH0 / "artifact_inventory.csv", index=False, encoding="utf-8")
    n_pred = int((inv.kind == "prediction_like").sum())
    print(f"制品清单：{len(inv)} 个文件（prediction_like {n_pred}）")

    proto = {
        "benchmark": "CR-PPI-ddG (Capability-Resolved PPI ddG Benchmark)",
        "status": "working name, not publication-final",
        "axes": ["A_calibration", "B_within_target_ranking", "C_selection_utility"],
        "no_universal_composite_score_yet": True,
        "primary_eligibility_n": 5,
        "sensitivity_eligibility_n": [3, 10],
        "bootstrap": {"requested": 5000, "cluster": "target", "seed": 42,
                      "require_valid_5000": True},
        "weighting_schemes": ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"],
        "k_primary": 10,
        "practical_tie_tol": {"overall_pearson": 0.01, "macro_spearman": 0.01},
        "target_namespace": "canonical_target_labels.csv joined by row_id",
        "frozen_split": "OUTER_FOLDS_V7_2",
        "forbidden_target_identity": ["block_L.pdb_id", "mutation_sites.pdb_id"],
        "reuse_first": True, "retrain_last": True, "gpu_last_resort": True,
        "hard_stops": ["restart V7.4R2", "retrain PLMs", "rerun ESM embeddings",
                       "rerun OpenMM", "new neural architectures",
                       "new handcrafted feature families"],
        "language": "task-specific family-OOD",
        "n_workers": 2,
    }
    (PH0 / "protocol.json").write_text(json.dumps(proto, indent=2, ensure_ascii=False),
                                       encoding="utf-8")
    (PH0 / "provenance.md").write_text(
        "# V8 PHASE 0 provenance\n\n"
        "* 本目录及全部 V8 产物**只读复用**上游；**不修改** V1–V7 任何输出。\n"
        "* 权威 target 命名空间：`canonical_target_labels.csv`，以 **`row_id`** 连接。\n"
        "* **禁止**以 `block_L.pdb_id` / `mutation_sites.pdb_id` 作为 target 身份。\n"
        "* 冻结 split：`OUTER_FOLDS_V7_2`。\n"
        "* 平台：**本地 CPU 优先**；GPU（Kaggle）仅在 §51 全部条件满足时使用。\n"
        "* `V7.4R2_STATUS = PAUSED_RESOURCE_LIMIT`；V8 期间不继续 V7.4R2。\n",
        encoding="utf-8")
    return rc


if __name__ == "__main__":
    sys.exit(main())
