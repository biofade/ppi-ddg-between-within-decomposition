# V7.4R Track B — EvoEF2 输入准备

- 来源 U_RELAX_STRICT：**3253** 行 / **162** targets
- 写出蛋白链 PDB：**162** 个（去水/配体/离子/氢）
- 变体数：**3253**；跳过 **0**（位点/结构缺失）
- 多点突变（n_sites>=2）：**704**

## `individual_list.txt` 语法示例（前 5 行）

```
PA2A;
TA3A;
IA4A;
PA5A;
LA6A;
```

> 语法：`<WT><chain><pdb_resseq><MUT>`，逗号分隔、**无空格**、以 `;` 结尾

## 执行顺序（§35）

```bash
# 1) 每个 target repair 一次（只做一次，缓存复用）
EvoEF2 --command=RepairStructure --pdb=<T>_protein.pdb --num_of_runs=3
#    -> <T>_protein_Repair.pdb

# 2) 逐个变体 build mutant
EvoEF2 --command=BuildMutant --pdb=<T>_protein_Repair.pdb \
       --mutant_file=individual_list.txt --num_of_runs=10
#    -> <T>_protein_Repair_Model_1.pdb

# 3) WT/MUT 各算一次 binding
EvoEF2 --command=ComputeBinding --pdb=<T>_protein_Repair.pdb
EvoEF2 --command=ComputeBinding --pdb=<T>_protein_Repair_Model_1.pdb
```

> ⚠️ `BuildMutant` 一次只产出一个 `Model_1.pdb`；多突变**写在同一条 token 内**
> （逗号分隔），不是多行。脚本 `12_evoef2_run.py` 需按此逐个调用。

> ⚠️ 多链复合物需显式 `--split=<binder>,<partner>`；见 `12_evoef2_run.py`。
