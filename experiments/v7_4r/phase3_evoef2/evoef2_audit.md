# V7.4R Track B — EvoEF2 溯源审计（§33）

| 项 | 值 |
|---|---|
| 引擎 | **EvoEF2** |
| 仓库 | `https://github.com/tommyhuangthu/EvoEF2` |
| 审计分支 | `master` |
| 许可 | **Free to academic users (per README 'COPYRIGHT & CONTACT'); not an OSI license** |
| 构建 | `./build.sh    OR    g++ -O3 --fast-math -o EvoEF2 src/*.cpp` |
| 科学定位 | **INDEPENDENT PHYSICS ENGINE, NOT OPENMM EQUIVALENT** |

## 命令与语法

- **`RepairStructure`**: `EvoEF2 --command=RepairStructure --pdb=your.pdb --num_of_runs=3   -> your_Repair.pdb`
- **`BuildMutant`**: `EvoEF2 --command=BuildMutant --pdb=your_Repair.pdb --mutant_file=individual_list.txt --num_of_runs=10 -> your_Repair_Model_1.pdb`
- **`ComputeBinding`**: `EvoEF2 --command=ComputeBinding --pdb=struct.pdb [--split=A,BC]`
- **`ComputeStability`**: `EvoEF2 --command=ComputeStability --pdb=protein.pdb`

**突变文件语法**（`individual_list.txt`）

- 每行一个 mutant，以 `;` 结尾
- 单点：`<WT_aa><chain><position><MUT_aa>  例: CA171A`
- 多点：`逗号分隔，无空格：CA171A,DB180E`（**逗号分隔、无空格**）

## ⚠️ 两项必须写进报告的溯源事实

### 1. 能量权重是在 ΔΔG 数据上**优化**过的

> EvoEF2 energy weights were optimized/reoptimized against experimental thermodynamic change (ddG) data and sequence-design benchmarks (README: 'optimized by maximizing its ability to predict the experimental thermodynamic change data (dG)' and 'reoptimizing energy weights' for EvoEF2).

**含义**：EvoEF2 is NOT a label-naive engine. Any comparison against V7.3 F4 must acknowledge that EvoEF2's function was tuned on ddG-like data, so it is an external *fitted* physics-inspired score, not a from-first-principles simulator.

### 2. 许可为**学术免费**，非 OSI 许可

> Academic-use-only. Must be stated in any publication. Not equivalent to a permissive open-source license.

## 已知约束

- 仅支持氨基酸/蛋白；**不**处理核苷酸、DNA、RNA、水或其他分子
- ⇒ 输入必须先去水/去配体（与 V7.4R 的 protein-chains-only 策略一致）
- 单次 BuildMutant 输出 your_Model_1.pdb（每个变异体一个文件）

## 术语纪律

- ✅ 允许：EvoEF2 binding score change、EvoEF2 energy-function delta
- ❌ 禁止：experimental binding free energy、true ddG、OpenMM-equivalent relaxation

→ `evoef2/v7_4r/provenance.json`
