# -*- coding: utf-8 -*-
"""scripts/v8/45_build_final_package.py -- 组装 `V8_FINAL_PACKAGE/` 并压成 ZIP。

目录结构（按用户指定）：
    V8_FINAL_PACKAGE/
    ├── README.md                      <- 本包导览 + 状态 + 复核入口
    ├── V8_FINAL_SUMMARY.md
    ├── EXTERNAL_AUDIT.md
    ├── BENCHMARK_NOVELTY_MATRIX.md
    ├── model_registry.csv
    ├── BENCHMARK_SCOREBOARD.csv
    ├── reports/      (关键分析报告)
    ├── external/     (PHASE 1 外部验证)
    ├── tables/       (T1..T8 + 索引)
    ├── figures/      (fig1..fig7)
    └── bootstrap/    (统计验证)

**排除**：模型权重、原始 PDB、venv、缓存、C: 盘写入。
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
V8 = R / "experiments" / "v8"
STAGE = R / "V8_FINAL_PACKAGE"
ZIP = R / "V8_FINAL_PACKAGE.zip"

# ---------- 顶层报告 -> reports/ 的归类 ----------
REPORTS = {
    "RANK_REVERSAL_REPORT.md": "关键分析",
    "WEIGHTING_SENSITIVITY_REPORT.md": "关键分析",
    "SPLIT_INFLATION_REPORT.md": "关键分析",
    "DECISION_ALIGNMENT_REPORT.md": "关键分析",
    "AXIS_C_DIAGNOSIS.md": "关键分析 (AXIS C floor)",
    "BENCHMARK_GUIDED_REMEDY.md": "关键分析 (H5)",
    "SCENARIO_RECOMMENDATIONS.md": "关键分析 (scenarios)",
    "PARETO_MODELS.md": "关键分析",
    # PHASE 13 / 方法学
    "EXTERNAL_REPLICATION_REPORT.md": "PHASE13 外部复现",
    "EXTERNAL_MULTIMODEL.md": "PHASE13 外部多模型",
    "EXTERNAL_HELDOUT_CATH.md": "PHASE13 外部 held-out",
    "EXTERNAL_POOL_EXPANSION.md": "PHASE13 扩池",
    "EXTERNAL_INFERENCE_FEASIBILITY.md": "PHASE13 推理可行性(阻断记录)",
    "CROSS_UNIVERSE_CONSISTENCY.md": "PHASE13 跨宇宙",
    "DISSOCIATION_BOUNDARY_CONDITIONS.md": "方法学 边界条件",
    "BOUNDARY_BOOTSTRAP.md": "方法学 bootstrap",
    "POOL_COMPOSITION_ARTIFACT.md": "方法学 池组成",
    "POOL_PROVENANCE_AUDIT.md": "方法学 池 provenance（★ 可复核）",
    "V8_LOCAL_CHECKPOINT.md": "阶段检查点",
}

EXTERNAL_KEEP = [
    "external_audit.csv", "hermes_audit.json", "hermes_matched_meta.json",
    "hermes_matched_final.csv", "external_multimodel_scoreboard.csv",
    "external_multimodel_meta.json", "external_heldout_cath.csv",
    "external_heldout_cath_meta.json", "external_expanded_pool.csv",
    "external_expanded_pool_meta.json", "cross_universe_external.csv",
    "cross_universe_agreement.csv", "cross_universe_meta.json",
    "external_replication_scoreboard.csv", "external_replication_meta.json",
    "pool_composition_models.csv", "pool_composition_meta.json",
    "pool_provenance_audit.csv", "pool_provenance_audit.json",
    "boundary_bootstrap.csv", "boundary_bootstrap_meta.json",
    "discrepancy_diagnosis.csv", "multipoint_mechanism_meta.json",
    "external_weight_probe.json",
]

# ---------- 清空并创建 ----------
if STAGE.exists():
    shutil.rmtree(STAGE)
if ZIP.exists():
    ZIP.unlink()
for sub in ("reports", "external", "tables", "figures", "bootstrap"):
    (STAGE / sub).mkdir(parents=True, exist_ok=True)
copied = []


def cp(src: Path, dst_dir: Path, name=None):
    if not src.exists():
        return False
    dst = dst_dir / (name or src.name)
    shutil.copy2(src, dst)
    copied.append(str(dst.relative_to(STAGE)))
    return True


# 1) 顶层
for f in ("README.md", "V8_FINAL_SUMMARY.md", "EXTERNAL_AUDIT.md",
          "BENCHMARK_NOVELTY_MATRIX.md", "model_registry.csv",
          "BENCHMARK_SCOREBOARD.csv", "EXTERNAL_MATCHED_SCOREBOARD.csv"):
    src = V8 / f
    if not src.exists():
        src = V8 / "external" / f          # EXTERNAL_MATCHED_SCOREBOARD 在 external/
    cp(src, STAGE, name=f)
    print(f"  [top]      {f}")

# 2) reports/
for f, tag in REPORTS.items():
    if cp(V8 / f, STAGE / "reports"):
        print(f"  [reports]  {f}")

# 3) external/
for f in EXTERNAL_KEEP:
    if cp(V8 / "external" / f, STAGE / "external"):
        pass
print(f"  [external] {len(list((STAGE/'external').iterdir()))} 文件")

# 4) tables/
for f in sorted((V8 / "tables").glob("*")):
    cp(f, STAGE / "tables")
print(f"  [tables]   {len(list((STAGE/'tables').iterdir()))} 文件")

# 5) figures/
for f in sorted((V8 / "figures").glob("*.png")):
    cp(f, STAGE / "figures", name=f.name)
print(f"  [figures]  {len(list((STAGE/'figures').iterdir()))} 文件")

# 6) bootstrap/（统计验证）
for f in sorted((V8 / "bootstrap").glob("*")):
    if f.suffix in (".csv", ".md", ".json"):
        cp(f, STAGE / "bootstrap")
cp(V8 / "metrics" / "decomposition_qc.csv", STAGE / "bootstrap")
cp(V8 / "external" / "boundary_bootstrap_clean.csv", STAGE / "bootstrap")
cp(V8 / "external" / "boundary_bootstrap_final.csv", STAGE / "bootstrap")
# pub 记分板与逐宇宙表
for f in sorted((V8 / "metrics").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="metrics_" + f.name)
for f in sorted((V8 / "rank_reversal").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="rankrev_" + f.name)
for f in sorted((V8 / "weighting").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="weight_" + f.name)
for f in sorted((V8 / "selection").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="axisC_" + f.name)
for f in sorted((V8 / "split_inflation").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="split_" + f.name)
for f in sorted((V8 / "scenarios").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="scen_" + f.name)
for f in sorted((V8 / "decision_alignment").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="align_" + f.name)
for f in sorted((V8 / "remedy").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="remedy_" + f.name)
for f in sorted((V8 / "pareto").glob("*.csv")):
    cp(f, STAGE / "bootstrap", name="pareto_" + f.name)
for f in sorted((V8 / "registry").glob("*")):
    cp(f, STAGE / "bootstrap", name="registry_" + f.name)
for f in sorted((V8 / "phase0_freeze").glob("*")):
    cp(f, STAGE / "bootstrap", name="freeze_" + f.name)
for f in sorted((V8 / "literature").glob("*")):
    cp(f, STAGE / "bootstrap", name="lit_" + f.name)
print(f"  [bootstrap] {len(list((STAGE/'bootstrap').iterdir()))} 文件")

# 7) 能力卡
(STAGE / "cards").mkdir(exist_ok=True)
for f in sorted((V8 / "cards").glob("*.md")):
    cp(f, STAGE / "cards")
print(f"  [cards]    {len(list((STAGE/'cards').iterdir()))} 文件")

# 8) 统计验证自检报告（由 _test_metrics 日志生成）
selftest = R / "logs" / "v8" / "metrics_test.log"
if selftest.exists():
    txt = selftest.read_text(encoding="utf-8", errors="replace")
    (STAGE / "bootstrap" / "METRIC_ENGINE_SELFTEST.md").write_text(
        "# V8 度量引擎自检报告\n\n"
        "来源：`python -W ignore -u scripts/v8/_test_metrics.py`\n"
        "（合成数据：perfect / reversed / noise / const 四个已知答案用例）\n\n"
        "```text\n" + txt.strip() + "\n```\n\n"
        "关键：**§28 分解恒等式 `r_full = C_B + C_W` 残差 ≤ 3.3e-16**（机器精度）。\n",
        encoding="utf-8")
    copied.append("bootstrap/METRIC_ENGINE_SELFTEST.md")
    print("  [selftest] METRIC_ENGINE_SELFTEST.md")

# ---------- 清单 + 哈希 ----------
manifest = []
for p in sorted(STAGE.rglob("*")):
    if p.is_file():
        h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        manifest.append({"path": str(p.relative_to(STAGE)),
                         "bytes": p.stat().st_size, "sha256_16": h})
(STAGE / "MANIFEST.json").write_text(
    json.dumps({"n_files": len(manifest), "files": manifest},
               indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n清单：{len(manifest)} 个文件")

# ---------- ZIP ----------
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in sorted(STAGE.rglob("*")):
        if p.is_file():
            z.write(p, p.relative_to(STAGE.parent))
print(f"ZIP -> {ZIP}  ({ZIP.stat().st_size/1048576:.2f} MB)")
