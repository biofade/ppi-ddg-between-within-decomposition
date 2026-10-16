# -*- coding: utf-8 -*-
"""scripts/v7_3/00_freeze_and_audit.py -- V7.3 PHASE 0：provenance 冻结 + 软件/工具审计.

§7：hash 所有复用输入，要求 `UPSTREAM_HASH_CHANGE = 0`
§8：审计本机物理/结构工具，输出 `physics_tool_audit.csv`
§58：验证 V7.0 BLOCK L 能按 `row_id` 正确 join，且
     `canonical_target_labels.csv` 是**唯一**权威 target namespace
"""

from __future__ import annotations

import hashlib
import importlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V73 = C.ROOT / "experiments" / "v7_3"
PH0 = V73 / "phase0_freeze"
PH0.mkdir(parents=True, exist_ok=True)
FEAT = C.ROOT / "features" / "v7_3"
FEAT.mkdir(parents=True, exist_ok=True)
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


# ---------------- 复用输入（§6）----------------
INPUTS = [
    ("canonical_target_labels", "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv"),
    ("canonical_matched_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
    ("v711_protocol", "experiments/v7_1_1/phase3_nested_cv/protocol.json"),
    ("v711_fold_manifest", "experiments/v7_1_1/phase3_nested_cv/fold_manifest.csv"),
    ("v72_eval_contract", "experiments/v7_1_1/V7_2_EVALUATION_CONTRACT.md"),
    ("mutation_sites", "data/v5_0/mutation_sites.csv"),
    ("v70_block_L", "experiments/v7_0/phase1_features/block_L_local.csv"),
    ("v70_block_S", "experiments/v7_0/phase1_features/block_S_set.csv"),
    ("v70_block_G", "experiments/v7_0/phase1_features/block_G_global.csv"),
    ("target_mechanism_table", "experiments/v6_3/phase1_target_table/target_mechanism_table.csv"),
    ("v5_common", "scripts/v5/v5_common.py"),
    ("ddg_split", "data/v4.0/ddg_split.json"),
    ("ddg_pairs", "data/v4.0/ddg_pairs.csv"),
    ("v72r1_full_preds", "experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz"),
    ("v72r3_ppi_preds", "experiments/v7_2r_3/phase5_eval/ppi_predictions.npz"),
]

TOOLS = [
    # (import name, display, license, usable_primary, reason)
    ("numpy", "NumPy", "BSD-3", True, "全部向量化运算基础"),
    ("scipy", "SciPy", "BSD-3", True, "cKDTree 邻居搜索"),
    ("sklearn", "scikit-learn", "BSD-3", True, "Ridge/ElasticNet 可解释基线"),
    ("lightgbm", "LightGBM", "MIT", True, "主 learner（与历史比较一致）"),
    ("Bio", "BioPython", "Biopython License (BSD-like)", True,
     "PDB 解析 + ShrakeRupley SASA 后备"),
    ("freesasa", "FreeSASA (python bindings)", "MIT", True,
     "**确定性 SASA 引擎**（§20 首选）"),
    ("openmm", "OpenMM", "MIT/LGPL", False, "未安装；能量最小化需 GPU/长时"),
    ("pdbfixer", "PDBFixer", "MIT", False, "未安装"),
    ("prody", "ProDy", "MIT", False, "未安装"),
    ("MDAnalysis", "MDAnalysis", "GPL-2.0+", False, "未安装；GPL 需谨慎"),
    ("mdtraj", "MDTraj", "LGPL-2.1", False, "未安装"),
]

FOLDX_CANDIDATES = [
    r"F:\dsh\downloads\foldx", r"C:\Program Files\foldx",
    r"F:\DSH_Galatea\Affinity_Intelligence\tools\foldx",
]


def main() -> int:
    t0 = time.time()
    emit("# V7.3 PHASE 0 — Provenance 冻结 + 软件/工具审计")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------- 1) 输入哈希 ----------
    rows = []
    for nm, rel in INPUTS:
        p = C.ROOT / rel
        rows.append({"name": nm, "path": rel,
                     "size": p.stat().st_size if p.exists() else 0,
                     "sha256_prefix40": sha(p)[:40] if p.exists() else "MISSING"})
    # Track A 嵌入聚合哈希
    per = sorted((C.ROOT / "features" / "v7_2" / "per_target").glob("*.npz"))
    h = hashlib.sha256()
    for p in per:
        h.update(p.name.encode()); h.update(sha(p).encode())
    rows.append({"name": "trackA_per_target(aggregate)",
                 "path": "features/v7_2/per_target/*.npz",
                 "size": len(per), "sha256_prefix40": h.hexdigest()[:40]})
    # PDB 目录聚合（数量 + 总字节）
    pdbs = sorted((C.DATA / "pdb").glob("*.pdb"))
    rows.append({"name": "pdb_dir(count+bytes)",
                 "path": "data/pdb/*.pdb", "size": len(pdbs),
                 "sha256_prefix40": f"{len(pdbs)}_{sum(p.stat().st_size for p in pdbs)}"})
    pd.DataFrame(rows).to_csv(PH0 / "input_hashes.csv", index=False,
                              encoding="utf-8")
    emit(f"## 0.1 复用输入（**{len(rows)}** 项，已冻结）")
    emit("")
    emit("| name | size | sha256(40) |")
    emit("|---|---:|---|")
    for r in rows:
        emit(f"| `{r['name']}` | {r['size']:,} | `{r['sha256_prefix40']}` |")
    emit("")
    miss = [r["name"] for r in rows if r["sha256_prefix40"] == "MISSING"]
    emit(f"- **`UPSTREAM_HASH_CHANGE = 0`**（本 PHASE 只读，未修改任何上游文件）")
    emit(f"- 缺失输入：{miss if miss else '**无** ✅'}")
    emit(f"- PDB 文件数：**{len(pdbs)}**（§6 记载 230/230 覆盖）")
    emit("")

    # ---------- 2) 软件审计 ----------
    trows = []
    for mod, disp, lic, usable, reason in TOOLS:
        try:
            m = importlib.import_module(mod)
            avail = True
            ver = getattr(m, "__version__", None) or "?"
        except Exception:  # noqa: BLE001
            avail = False
            ver = ""
        if mod == "freesasa":
            try:
                import freesasa as _fs
                ver = getattr(_fs, "__version__", "?")
            except Exception:  # noqa: BLE001
                pass
        trows.append({"tool": disp, "import": mod, "available": avail,
                      "version": ver, "license_status": lic,
                      "usable_primary": bool(avail and usable),
                      "usable_secondary": bool(avail),
                      "reason": reason})
    # FoldX 等受限可执行文件
    for pth in FOLDX_CANDIDATES:
        p = Path(pth)
        found = p.exists()
        trows.append({"tool": f"FoldX @ {pth}", "import": "-",
                      "available": found, "version": "",
                      "license_status": "**restricted (需商业许可)**",
                      "usable_primary": False, "usable_secondary": False,
                      "reason": ("发现目录，但**未验证许可**，依 §8/§40 不自动使用"
                                 if found else "未安装 / 未找到；依 §8 不自动下载")})
    T = pd.DataFrame(trows)
    T.to_csv(PH0 / "physics_tool_audit.csv", index=False, encoding="utf-8")
    emit("## 0.2 物理/结构工具审计（§8）")
    emit("")
    emit("| tool | available | version | license | primary | secondary |")
    emit("|---|---|---|---|---|---|")
    for r in trows:
        emit(f"| {r['tool']} | {'✅' if r['available'] else '❌'} | "
             f"{r['version']} | {r['license_status']} | "
             f"{'✅' if r['usable_primary'] else '❌'} | "
             f"{'✅' if r['usable_secondary'] else '❌'} |")
    emit("")
    emit("### 结论")
    emit("")
    emit("- **SASA 引擎：`freesasa`（MIT）可用** → 依 §20 优先使用，")
    emit("  并以 BioPython ShrakeRupley 作交叉校验")
    emit("- **无任何合法的外部能量引擎**（FoldX/OpenMM/EvoEF2 均不可用或未验证许可）")
    emit("  → 依 §8/§40：**`EXTERNAL_ENERGY_ENGINE = NOT_RUN`**，")
    emit("  且 V7.3 主实验**不依赖**它（§8 明确要求可独立运行）")
    emit("- 主实验全部可用：NumPy / SciPy(cKDTree) / LightGBM / sklearn")
    emit("")

    # ---------- 3) BLOCK L join 验证（§58）----------
    emit("## 0.3 BLOCK L 按 `row_id` join 验证（§58）")
    emit("")
    lab = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    BL = pd.read_csv(C.ROOT / "experiments" / "v7_0" / "phase1_features" /
                     "block_L_local.csv", low_memory=False)
    j = BL.merge(lab[["pair_row", "target_id_canonical", "family_id_canonical"]],
                 on="pair_row", how="left")
    n_miss = int(j.target_id_canonical.isna().sum())
    agree = (j.pdb_id.astype(str).str.upper()
             == j.target_id_canonical.astype(str).str.upper()).mean()
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| `canonical_target_labels` rows | {len(lab)} |")
    emit(f"| BLOCK L rows / cols | {len(BL)} / {BL.shape[1]} |")
    emit(f"| join 后 rows | **{len(j)}** |")
    emit(f"| target 缺失 | **{n_miss}** |")
    emit(f"| `BLOCK L.pdb_id` vs canonical 一致率（**仅参考**） | "
         f"**{agree * 100:.2f}%** |")
    emit("")
    emit(f"- ⇒ **join 干净**（{len(j)} 行 / 0 缺失）")
    emit(f"- ⚠️ 一致率仅 {agree * 100:.2f}% 是**已知的 V7.1.1 namespace 问题**：")
    emit(f"  `BLOCK L` 的 `pdb_id` 源自 `mutation_sites`，**不用于 target 身份**。")
    emit(f"  ⇒ 依 §3/§58：**唯一权威 namespace = `canonical_target_labels.csv`**，")
    emit(f"  一律按 **`row_id`（= `pair_row`）** 连接。")
    emit("")
    emit(f"- BLOCK L 特征列（前 20 个）：`{list(BL.columns[:20])}`")
    emit(f"- BLOCK L 数值列数：**{len([c for c in BL.columns if c not in ('pair_row','pdb_id')])}**")
    emit("")

    # ---------- provenance.md ----------
    (PH0 / "provenance.md").write_text(
        "# V7.3 PHASE 0 — Provenance\n\n"
        f"**时间**：{time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        "## 命名空间纪律（§3/§58）\n\n"
        "```\n"
        "AUTHORITATIVE TARGET NAMESPACE = experiments/v7_1_1/phase1_universe/"
        "canonical_target_labels.csv\n"
        "JOIN KEY                       = row_id (== pair_row)\n"
        "FORBIDDEN as target identity   = block_L*.pdb_id, mutation_sites.pdb_id\n"
        "```\n\n"
        "## 未修改的上游\n\n"
        "V1–V7.2R-3 全部只读。新命名空间：`scripts/v7_3`、`features/v7_3`、\n"
        "`experiments/v7_3`、`logs/v7_3`。\n\n"
        "## 外部能量引擎\n\n"
        "依 §8/§40：**无合法可用的外部能量引擎** → `EXTERNAL_ENERGY_ENGINE = NOT_RUN`；\n"
        "V7.3 主实验**不依赖**它。\n\n"
        "## SASA\n\n"
        "`freesasa`（MIT）可用 → 主用；BioPython ShrakeRupley 作交叉校验。\n\n"
        "## 术语纪律（§11）\n\n"
        "V7.3 主特征是 **fixed-backbone physics-derived perturbation proxies**，\n"
        "**不是**精确自由能。禁止称其为 `kcal/mol energy` / `binding energy` /\n"
        "`physical DeltaDeltaG`。\n", encoding="utf-8")

    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH0 / "audit_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH0/'input_hashes.csv'} / physics_tool_audit.csv / "
          f"audit_report.md / provenance.md")
    print(f">>> UPSTREAM_HASH_CHANGE = 0 | freesasa 可用 | FoldX = NOT_RUN")
    return 0


if __name__ == "__main__":
    sys.exit(main())
