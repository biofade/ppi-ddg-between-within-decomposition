# -*- coding: utf-8 -*-
"""scripts/v7_4r/10_evoef2_audit.py -- V7.4R Track B：EvoEF2 溯源审计（§33）.

审计来源：官方 README（`tommyhuangthu/EvoEF2`，master）
    https://raw.githubusercontent.com/tommyhuangthu/EvoEF2/master/README.md
    https://github.com/tommyhuangthu/EvoEF2

**科学定位（§32）**：EvoEF2 是**独立的物理引擎验证**，
**绝不**等于 OpenMM 松弛。EvoEF2 用的是**优化过的经验能量函数**
（其权重是针对 ΔΔG 预测与序列设计**拟合**出来的），
而 OpenMM 用的是**力场 + 显式几何最小化**。两者机制不同。

同时记录一个**必须写进报告**的溯源事实：EvoEF2 的能量权重
是**在实验 ΔΔG 数据上优化过的**（README 原文：
"optimized by maximizing its ability to predict the experimental
thermodynamic change data (ΔΔG)"）。
⇒ 它不是一个「未见过标签」的引擎；这在解释其表现时**必须**说明。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "evoef2" / "v7_4r"
OUT.mkdir(parents=True, exist_ok=True)
PH3 = R / "experiments" / "v7_4r" / "phase3_evoef2"
PH3.mkdir(parents=True, exist_ok=True)
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))


PROV = {
    "engine": "EvoEF2",
    "purpose": "INDEPENDENT PHYSICS-ENGINE VALIDATION",
    "not_equivalent_to": "OpenMM relaxation (different mechanism: fitted empirical energy function vs force-field minimization)",
    "repository": "https://github.com/tommyhuangthu/EvoEF2",
    "readme": "https://raw.githubusercontent.com/tommyhuangthu/EvoEF2/master/README.md",
    "branch_audited": "master",
    "author": "Xiaoqiang Huang",
    "contact": "xiaoqiah@umich.edu / xiaoqiah@outlook.com",
    "license": "Free to academic users (per README 'COPYRIGHT & CONTACT'); not an OSI license",
    "license_caveat": ("Academic-use-only. Must be stated in any publication. "
                       "Not equivalent to a permissive open-source license."),
    "citation": [
        "Huang X, Pearce R, Zhang Y. EvoEF2: accurate and fast energy function "
        "for computational protein design. Bioinformatics (2020) 36:1135-1142",
        "Pearce R, Huang X, Setiawan D, Zhang Y. EvoDesign: Designing "
        "Protein-Protein Binding Interactions Using Evolutionary Interface "
        "Profiles in Conjunction with an Optimized Physical Energy Function. "
        "J. Mol. Biol. (2019) 431:2467-2476",
    ],
    "build": {
        "linux_mac": "./build.sh    OR    g++ -O3 --fast-math -o EvoEF2 src/*.cpp",
        "linux_fallback": "g++ -O3 -o EvoEF2 src/*.cpp   (若 --fast-math 失败)",
        "mac_note": "Mac 用 -fast-math",
        "windows": "可用预编译 EvoEF2.exe；否则 MinGW g++ 自行编译",
        "cpp_standard": "C++ (g++ -O3)",
    },
    "capabilities_used": ["RepairStructure", "BuildMutant", "ComputeBinding"],
    "capabilities_available": ["ComputeStability", "ComputeBinding",
                               "RepairStructure", "BuildMutant",
                               "OptimizeHydrogen", "ProteinDesign"],
    "mutation_file_syntax": {
        "per_line": "CA171A,DB180E;   (每行一个 mutant，以 ';' 结尾)",
        "single_mutation": "<WT_aa><chain><position><MUT_aa>  例: CA171A",
        "multi_mutation": "逗号分隔，无空格：CA171A,DB180E",
        "note": "第一位 = 参考氨基酸；第二位 = 链 ID；随后为链内位置；最后 = 突变后氨基酸",
    },
    "commands": {
        "RepairStructure": ("EvoEF2 --command=RepairStructure --pdb=your.pdb "
                            "--num_of_runs=3   -> your_Repair.pdb"),
        "BuildMutant": ("EvoEF2 --command=BuildMutant --pdb=your_Repair.pdb "
                        "--mutant_file=individual_list.txt --num_of_runs=10 "
                        "-> your_Repair_Model_1.pdb"),
        "ComputeBinding": ("EvoEF2 --command=ComputeBinding --pdb=struct.pdb "
                           "[--split=A,BC]"),
        "ComputeStability": "EvoEF2 --command=ComputeStability --pdb=protein.pdb",
    },
    "known_constraints": [
        "仅支持氨基酸/蛋白；**不**处理核苷酸、DNA、RNA、水或其他分子",
        "⇒ 输入必须先去水/去配体（与 V7.4R 的 protein-chains-only 策略一致）",
        "单次 BuildMutant 输出 your_Model_1.pdb（每个变异体一个文件）",
    ],
    "weight_fitting_provenance": {
        "fact": ("EvoEF2 energy weights were optimized/reoptimized against "
                 "experimental thermodynamic change (ddG) data and sequence-design "
                 "benchmarks (README: 'optimized by maximizing its ability to "
                 "predict the experimental thermodynamic change data (dG)' and "
                 "'reoptimizing energy weights' for EvoEF2)."),
        "implication": ("EvoEF2 is NOT a label-naive engine. Any comparison against "
                        "V7.3 F4 must acknowledge that EvoEF2's function was tuned "
                        "on ddG-like data, so it is an external *fitted* physics-"
                        "inspired score, not a from-first-principles simulator."),
    },
    "scientific_status": "INDEPENDENT PHYSICS ENGINE, NOT OPENMM EQUIVALENT",
    "recommended_execution_environment": ("WSL2 Ubuntu (preferred if available) "
                                          "or Google Colab CPU runtime; GPU not required"),
    "term_discipline": {
        "allowed": ["EvoEF2 binding score change",
                    "EvoEF2 energy-function delta"],
        "forbidden": ["experimental binding free energy",
                      "true ddG", "OpenMM-equivalent relaxation"],
    },
}
(OUT / "provenance.json").write_text(json.dumps(PROV, indent=2, ensure_ascii=False),
                                     encoding="utf-8")

emit("# V7.4R Track B — EvoEF2 溯源审计（§33）")
emit("")
emit("| 项 | 值 |")
emit("|---|---|")
emit(f"| 引擎 | **EvoEF2** |")
emit(f"| 仓库 | `{PROV['repository']}` |")
emit(f"| 审计分支 | `master` |")
emit(f"| 许可 | **{PROV['license']}** |")
emit(f"| 构建 | `{PROV['build']['linux_mac']}` |")
emit(f"| 科学定位 | **{PROV['scientific_status']}** |")
emit("")
emit("## 命令与语法")
emit("")
for k, v in PROV["commands"].items():
    emit(f"- **`{k}`**: `{v}`")
emit("")
emit("**突变文件语法**（`individual_list.txt`）")
emit("")
emit(f"- 每行一个 mutant，以 `;` 结尾")
emit(f"- 单点：`{PROV['mutation_file_syntax']['single_mutation']}`")
emit(f"- 多点：`{PROV['mutation_file_syntax']['multi_mutation']}`（**逗号分隔、无空格**）")
emit("")
emit("## ⚠️ 两项必须写进报告的溯源事实")
emit("")
emit("### 1. 能量权重是在 ΔΔG 数据上**优化**过的")
emit("")
emit(f"> {PROV['weight_fitting_provenance']['fact']}")
emit("")
emit(f"**含义**：{PROV['weight_fitting_provenance']['implication']}")
emit("")
emit("### 2. 许可为**学术免费**，非 OSI 许可")
emit("")
emit(f"> {PROV['license_caveat']}")
emit("")
emit("## 已知约束")
emit("")
for c in PROV["known_constraints"]:
    emit(f"- {c}")
emit("")
emit("## 术语纪律")
emit("")
emit(f"- ✅ 允许：{'、'.join(PROV['term_discipline']['allowed'])}")
emit(f"- ❌ 禁止：{'、'.join(PROV['term_discipline']['forbidden'])}")
emit("")
emit(f"→ `evoef2/v7_4r/provenance.json`")
(PH3 / "evoef2_audit.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print(f"\n写出 provenance.json 与 evoef2_audit.md")
