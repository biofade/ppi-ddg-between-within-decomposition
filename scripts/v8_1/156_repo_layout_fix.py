# -*- coding: utf-8 -*-
"""scripts/v8_1/156_repo_layout_fix.py -- 修正 staging 目录布局，使仓库自带校验器可直接运行。

验证器 (validate_final_manuscript.py) 假定如下相对布局：
    <root>/v8_1/...
    <root>/latex_v6/*.pdf
    <root>/figures_v5/*.pdf
    <root>/experiments/...
    <root>/data/v5_0/mutation_sites.csv
    <root>/FINAL_SUBMISSION/AUTHOR_METADATA.json
因此 staging 顶层必须复刻这些路径（而不是放在 artifacts/ 下），
否则作者 clone 下来无法自校验。
"""
import shutil
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
STAGE = R / "FINAL_SUBMISSION/repo_staging"
assert STAGE.exists(), "请先运行 153_build_repo_package.py"

# 需要复制到 staging 顶层的路径（保持相对结构）
NEED = [
    "v8_1",
    "latex_v6",
    "figures_v5",
    "experiments",
    "data/v5_0/mutation_sites.csv",
    "scripts",
    "validate_final_manuscript.py",
]
SKIP_SUFFIX = (".aux", ".log", ".out", ".blg", ".bbl", ".fls", ".fdb_latexmk",
               ".synctex.gz")


def ignore(d, names):
    return [n for n in names if n.endswith(SKIP_SUFFIX) or n == "__pycache__"]


for rel in NEED:
    src = R / rel
    if not src.exists():
        print(f"  [缺失] {rel}"); continue
    dst = STAGE / rel
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True, ignore=ignore)
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    print(f"  copied {rel}")

# 清掉旧的 artifacts/ 里的重复内容（改为顶层布局后不再需要）
old = STAGE / "artifacts"
if old.exists():
    # 只删除与顶层重复的部分，保留 v8_1 的副本无意义 → 直接移除整个 artifacts
    shutil.rmtree(old)
    print("  removed 旧 artifacts/（内容已在顶层）")

# metadata
(STAGE / "FINAL_SUBMISSION").mkdir(exist_ok=True)
shutil.copy2(R / "FINAL_SUBMISSION/AUTHOR_METADATA.json",
             STAGE / "FINAL_SUBMISSION/AUTHOR_METADATA.json")
print("  copied FINAL_SUBMISSION/AUTHOR_METADATA.json")

total = sum(f.stat().st_size for f in STAGE.rglob("*") if f.is_file())
n = sum(1 for f in STAGE.rglob("*") if f.is_file())
print()
print(f"  staging: {n:,} 文件  {total/1e6:.1f} MB")
