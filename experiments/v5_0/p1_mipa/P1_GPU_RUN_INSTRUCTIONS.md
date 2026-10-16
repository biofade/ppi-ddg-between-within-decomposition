# P1 MIPA — GPU / 外部资源执行说明

> 状态：**`NOT_RUN_RESOURCE_BLOCKED`**（本机无 GPU；且当前映射质量不足以支撑 full MIPA）
> 本文件给出**可直接运行**的完整流程，使该实验在具备 GPU 与官方序列资源时可一次跑通。
> 严格遵循 V5 §36：**禁止伪造 full 结果**；被阻塞部分必须在报告中标记为
> `NOT_RUN_RESOURCE_BLOCKED`。

---

## 1. 为什么当前不能做 full MIPA（PROJECT_EVIDENCE）

| 步骤 | 结果 | 结论 |
|---|---|---|
| `03a_seq_sanity.py` | 配体 36.4% / 受体 24.1% 的 pair **两侧序列相同**；朴素 1-based 位置命中率 73.37%（0-based 0.43%） | 存在编号 offset，且并非每侧都变化 |
| `04_build_mutation_sites.py`（v1） | 自报覆盖 74.06% | 但见下一行 |
| `05_validate_mapping_pdb.py` | 用**真实 PDB 残基序号**独立验证：E1 通过率 **7.96%**（568/7134）；99/226 pdb 的 offset 自相矛盾 | **v1 映射 `INVALID_DUE_TO_BUG`** |
| `06_build_residue_map.py`（v2，序列引导） | 链 identity=1.000 占 63.06%；`pdb_resseq_real` 处真实残基 == annotation WT 残基 = **76.82%**（5930/7719）；可用位点 **4828/7869 = 61.35%** | **GATE B = `SOLVE_MAPPING_FIRST`**（<70%） |

**根因（已定位）**：数据集只提供 `Ligand Sequences` / `Receptor Sequences`
（多链用 `,` 拼接），**没有**「dataset 序列 ↔ 具体 PDB 链」的显式对应表，
也**没有**独立的 WT 序列字段。仅靠 WT/mutant 两条序列互相 alignment，
在 22.07% 的位点上出现 `AA_MISMATCH`（映射到的真实残基与 annotation 的 WT 残基不符），
说明两侧序列并非严格同 construct（可能含不同长度构建体 / 缺失残基差异）。

**按 V5 §124**：覆盖率 61.35% < 70% → **先解决 mapping**，不做 full MIPA。
**按 V5 §40**：不得把该结果解读为「MIPA 假设被证伪」——只能报
「在现有映射质量下无法对 MIPA 下结论」。

---

## 2. 解决映射的两条正路（按优先级）

### 路线 A（推荐，不需要 GPU）：用 PDB 官方 SEQRES + 逐残基比对

本地已有 **3031 个 PDB**（`data/pdb/`）。关键改进点：

1. **不要从 mutant 序列反推 WT**；改为
   `PDB chain SEQRES` ↔ `dataset Ligand/Receptor 子链` 做 alignment，
   用 **SEQRES（完整构建体）** 而不是只含 ATOM 的记录。
2. 处理 PDB 缺失残基：ATOM 中缺的残基在 SEQRES 里存在 → 建立
   「SEQRES 下标 ↔ ATOM 残基号」双表，避免因缺失残基导致的**系统性错位**。
3. 对 `,` 拼接的多链，做**子链级** alignment（已实现于 `06`），
   并对每个子链单独解 offset，**不允许跨子链共享 offset**。
4. 用 **insertion code** 与 **altLoc** 显式处理（`06` 已读 icode）。

预计可把 `AA_MISMATCH` 从 22.07% 降到 <5%，从而把覆盖推到 ≥90%。

### 路线 B（需要 GPU + 外部资源）：重跑 per-residue ESM-2，按**真实残基**索引

现有 `features/v2_esm2_perres_v2.npz`（1.9 GB）是**按 dataset 序列下标**
提取的残基表示，因此它同样继承了上面的映射不确定性。要支持 MIPA 必须
**按 PDB 真实链 + 残基号**重新提取。

---

## 3. 可直接运行的 MIPA 提取流程（GPU）

### 3.1 输入

```
data/v5_0/residue_map.csv        # 逐位点: pair_row, side, pdb_chain, pdb_resseq,
                                 #         icode, seq_index0, wt_aa, real_aa, mapping_status
data/pdb/<PDB>.pdb               # 本地结构
features/v2.2_proteinea_features.csv
data/v4.0/ddg_pairs.csv
data/v4.0/ddg_split.json
```

### 3.2 提取脚本骨架（把下面内容存为 `scripts/v5/gpu_mipa_extract.py` 并在 GPU 机器运行）

```python
# 关键约束（与 V5 §28–§37 一致）
#  * 模型: esm2_t30_150M_UR50D（**不换大模型**，V5 §119）
#  * 层: 30（最后一层）；池化: 不需要（我们要逐残基）
#  * 只缓存 mutation site 与 local partner residues，**不缓存全蛋白残基**
#  * 输出 fp16 embedding + float32 距离 + offset 数组（V5 §37）

import numpy as np, pandas as pd, torch
import esm  # fair-esm
from Bio.PDB import PDBParser, NeighborSearch

R_CUT = 8.0           # V5 §27 主值
SIGMA = 4.0           # V5 §29 主值

rm   = pd.read_csv("data/v5_0/residue_map.csv", low_memory=False)
use  = rm[rm.mapping_status.isin(["EXACT", "ALIGNMENT_CONFIDENT"])].copy()
print("mutation sites usable:", len(use), "pairs:", use.pair_row.nunique())

model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
model = model.eval().cuda()
bc = alphabet.get_batch_converter()

# 1) 对每个 (pdb_id, pdb_chain) 提取一次该链的逐残基表示（缓存复用）
#    注: dataset 序列可能含非 PDB 残基，因此用 PDB 的 SEQRES 序列作为模型输入，
#    再用 SEQRES↔ATOM 映射把表示对齐到真实残基号。
chain_cache = {}
for (pid, ch), g in use.groupby(["pdb_id", "pdb_chain"]):
    seq = seqres_of(pid, ch)                    # 从 PDB SEQRES 读，非 ATOM
    _, _, toks = bc([(f"{pid}_{ch}", seq)])
    with torch.no_grad():
        rep = model(toks.cuda(), repr_layers=[30], return_contacts=False)["representations"][30]
    rep = rep[0, 1:1+len(seq)]                  # 去 BOS/EOS
    chain_cache[(pid, ch)] = rep.half().cpu().numpy()

# 2) 对每个 mutation site 计算 MIPA 事件向量
#    m = mutation site; partner chain = 另一侧的链
#    Δh_m = h_m(mut) - h_m(WT)
#    多链拼接情形: 对每个 dataset 子链分别求 Δh，再按 §32 做 SUM/MEAN
#    c_m  = Σ_j w_j h_j / Σ w_j,  w_j = exp(-d_mj² / (2σ²))   (V5 §29)
#    blocks: Δh, c, Δh⊙c, cos(Δh,c), ||Δh||, ||c||          (V5 §30)
# 3) 每个 640D block 在 **train-only** 上 PCA 到 32 维（V5 §31）
# 4) 落盘: mipa_cache.npz（fp16）
```

### 3.3 必须实现的断言

```python
assert set(pca_fit_ids) <= set(train_ids)          # V5 §87
assert len(train_pdb & test_pdb) == 0              # V5 §86
assert len(train_pdb & val_pdb) == 0
assert oof_row_not_in_its_fold_train               # V5 §88
assert fwd_and_rev_in_same_split                   # V5 §89
```

### 3.4 运行顺序

```bash
python scripts/v5/06_build_residue_map.py            # CPU，先修映射
python scripts/v5/gpu_mipa_extract.py --device cuda  # GPU，提取 MIPA 缓存
python scripts/v5/10_p1_mipa.py --cache mipa_cache.npz   # CPU，跑消融
```

---

## 4. P1 最小消融矩阵（脚本已就绪，等缓存）

必须在**完全同一 test rows** 上比较（V5 §33），并做 coverage fairness（V5 §34）：

| ID | 内容 |
|---|---|
| `P1-0` | 现有 pooled `D` baseline |
| `P1-1` | 仅 mutation-site `Δh` |
| `P1-2` | 仅 partner context `c` |
| `P1-3` | `Δh` + `c` |
| `P1-4` | `Δh` + `c` + elementwise interaction `Δh⊙c` |
| `P1-5` | 仅 interaction |
| `P1-6` | pooled `D` + MIPA |

**覆盖公平性三档**（V5 §34，V2.5 曾在此吃亏）：

```
FULL baseline          （全部 test rows）
SUBSET-matched baseline（仅 MIPA 可覆盖的 test rows，重跑）
MIPA subset            （同 subset）
```

---

## 5. 成功判据（V5 §39）

- `Δ overall Pearson ≥ +0.03` **且** `Δ macro per-structure Spearman ≥ +0.05`
- 3 seeds 中 ≥2 个方向一致
- cluster bootstrap 至少一个核心指标 `q < 0.05`
- breakthrough 线：macro per-structure Spearman 接近 **0.20–0.30+**

**当前 baseline 参考值**（本项目 `PROJECT_EVIDENCE`，cluster bootstrap，n_clusters=61）：
overall Pearson **+0.3861**（3-seed 均）、macro-PDB Spearman **+0.0956**。

---

## 6. 若映射最终无法解决

则按 V5 §123 记录为阻塞，并把 P1 降级为：
**「mutation-centered representation 在当前数据可映射性下不可检验」** —— 这是一条
**有信息量的负结论**（说明数据本身的构建体标注是限制因素），而不是「MIPA 无效」。
