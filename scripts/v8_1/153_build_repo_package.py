# -*- coding: utf-8 -*-
"""scripts/v8_1/153_build_repo_package.py -- 生成要上传 GitHub 的仓库内容（staging 目录）。

策略：
  * **自产内容**（分析脚本、冻结 artifact、LaTeX 源、矢量图、校验器）→ 直接打包上传。
  * **第三方数据**（SKEMPI 派生表 skempi_v2.csv、CATH-ddG 仓库等）→ **不重新分发**，
    只记录来源 URL 与 SHA-256，附 fetch 脚本让读者自行下载。避免许可证问题。
  * 同时产出 REPO_MANIFEST.csv（路径/字节/SHA-256）与上传说明。

输出目录：FINAL_SUBMISSION/repo_staging/
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
STAGE = R / "FINAL_SUBMISSION/repo_staging"
if STAGE.exists():
    shutil.rmtree(STAGE)
STAGE.mkdir(parents=True)

# ---------------- 1) 自产内容 ----------------
COPY = [
    # (源, 目标, 说明)
    ("scripts/v8_1", "scripts/v8_1", "V8.1 分析/审计/制稿/校验脚本"),
    ("scripts/v6_2", "scripts/v6_2", "V6.2 四档受控划分与分解（含 L3 键级修复）"),
    ("scripts/v7_1", "scripts/v7_1", "V7.1 ranker 与嵌套 CV"),
    ("v8_1", "artifacts/v8_1", "冻结 artifact：所有头条数字的机器可读来源"),
    ("latex_v6", "manuscript/latex_v6", "LaTeX 源（主稿/补充/参考文献/构建脚本）"),
    ("figures_v5", "manuscript/figures", "图（PDF 矢量 + PNG + 生成脚本）"),
    ("validate_final_manuscript.py", "validate_final_manuscript.py", "自动一致性校验器"),
    ("FINAL_SUBMISSION/VALIDATION_REPORT.md", "VALIDATION_REPORT.md", "校验报告"),
    ("FINAL_SUBMISSION/CHANGELOG_FINAL.md", "CHANGELOG_FINAL.md", "变更记录"),
    ("FINAL_SUBMISSION/ROUND8_FIX_REPORT.md", "reports/ROUND8_FIX_REPORT.md", "审计修复记录"),
    ("FINAL_SUBMISSION/RECOMPUTE_FEASIBILITY.md", "reports/RECOMPUTE_FEASIBILITY.md",
     "键级泄漏评估与修复"),
    ("FINAL_SUBMISSION/AUTHOR_METADATA.template.json",
     "metadata/AUTHOR_METADATA.template.json", "metadata 模板"),
]

rows = []
for src, dst, desc in COPY:
    sp = R / src
    if not sp.exists():
        print(f"  [跳过] {src} 不存在"); continue
    dp = STAGE / dst
    if sp.is_dir():
        shutil.copytree(sp, dp, dirs_exist_ok=True)
    else:
        dp.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sp, dp)
    rows.append((dst, desc))

# 冻结输入清单（experiments 下的小体积 npz/csv —— 复现必需）
FROZEN = [
    "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
    "experiments/v7_1_1/phase3_nested_cv/fold_manifest.csv",
    "experiments/v7_1_1/phase3_nested_cv/protocol.json",
    "experiments/v7_3/phase5_rankers/family_predictions.npz",
    "experiments/v7_1/phase2_rank/z_predictions.npz",
    "experiments/v6_2/phase2_split4/split4_predictions.npz",
    "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
    "experiments/v6_2/phase2_split4/family_labels.csv",
    "experiments/v6_2/phase2_split4/split4_results.csv",
    "experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz",
    "experiments/v7_2r_2/phase5_eval/trackB_predictions.npz",
    "experiments/v7_2r_3/phase5_eval/ppi_predictions.npz",
    "experiments/v8/BENCHMARK_SCOREBOARD.csv",
    "data/v5_0/mutation_sites.csv",
]
for f in FROZEN:
    sp = R / f
    if not sp.exists():
        print(f"  [跳过-缺失] {f}"); continue
    dp = STAGE / "artifacts" / f
    dp.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(sp, dp)
    rows.append(("artifacts/" + f, "冻结输入（复现必需）"))

# ---------------- 2) 第三方数据：只记录，不分发 ----------------
THIRD = [
    ("F:/dsh/downloads/cathddg/skempi_v2.csv",
     "https://github.com/<CATH-ddG-repo>/.../data/SKEMPI2/skempi_v2.csv",
     "SKEMPI-derived annotation table distributed with CATH-ddG; source of the five "
     "external prediction columns"),
]
third_rows = []
for path, url, note in THIRD:
    p = Path(path)
    if p.exists():
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        third_rows.append(dict(path=path, bytes=p.stat().st_size, sha256=h,
                               source_url=url, note=note))
    else:
        third_rows.append(dict(path=path, bytes=None, sha256=None,
                               source_url=url, note=note + " (not present locally)"))

# ---------------- 3) REPO_MANIFEST.csv ----------------
man = []
for dp in sorted(STAGE.rglob("*")):
    if dp.is_file():
        man.append(dict(relpath=str(dp.relative_to(STAGE)).replace("\\", "/"),
                        bytes=dp.stat().st_size,
                        sha256=hashlib.sha256(dp.read_bytes()).hexdigest()))
pd.DataFrame(man).to_csv(STAGE / "REPO_MANIFEST.csv", index=False, encoding="utf-8")
pd.DataFrame(third_rows).to_csv(STAGE / "THIRD_PARTY_INPUTS.csv", index=False,
                                encoding="utf-8")

total = sum(m["bytes"] for m in man)
print("=" * 96)
print("仓库 staging 生成完成")
print("=" * 96)
print(f"  文件数 = {len(man):,}    总体积 = {total/1e6:.1f} MB")
print(f"  输出目录 = {STAGE}")
print()
print("  顶层结构:")
for d in sorted({m["relpath"].split("/")[0] for m in man}):
    n = sum(1 for m in man if m["relpath"].startswith(d + "/") or m["relpath"] == d)
    sz = sum(m["bytes"] for m in man if m["relpath"].startswith(d + "/") or m["relpath"] == d)
    print(f"    {d:<24} {n:>4} 文件  {sz/1e6:>7.2f} MB")
print()
print("  第三方输入（不分发，仅记录来源与校验和）:")
for t in third_rows:
    print(f"    {Path(t['path']).name:<24} sha256={(t['sha256'] or 'n/a')[:16]}  {t['source_url'][:52]}")
json.dump(dict(n_files=len(man), total_bytes=total,
               third_party=third_rows),
          open(STAGE / "REPO_PACKAGE_SUMMARY.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print()
print("  已写 REPO_MANIFEST.csv / THIRD_PARTY_INPUTS.csv / REPO_PACKAGE_SUMMARY.json")
