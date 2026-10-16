# V6.2 PHASE 1 — Target Provenance Audit

- 时间: 2026-10-14 08:50:15
- **判定原则**：只依据 provenance / 列名 / 公式 / 行宇宙 / hash，
  **不依据「数值看起来是否合理」**（V6.2 §11）

## 证据 1 — V6.1 PHASE 2 脚本中的 target 赋值语句

- 脚本: `scripts/v6_1/04_split_inflation.py`
- SHA256: `acb1a415e6f7e8fb8c635f518e1498d2…`

相关源码行（target 相关）：

```python
y = pd.to_numeric(v22["real_dG"], errors="coerce").to_numpy(float)
y = d["y"][idx_all]
```

- **精确匹配**：`y = pd.to_numeric(v22["real_dG"]...)`
- 即 target 取自 `features/v2.2_proteinea_features.csv` 的 **`real_dG`** 列


## 证据 2 — 行宇宙对比（绝对值表 vs ΔΔG 配对表）

| 表 | 行数 | 唯一 PDB | 唯一 complex/interface | 时间戳列 |
|---|---:|---:|---:|---|
| `features/v2.2_proteinea_features.csv` | **8207** | 2171 | 3197 | — |
| `data/v4.0/ddg_pairs.csv`（全部） | **15052** | — | — | 含 `direction` |
| `data/v4.0/ddg_pairs.csv`（forward） | **7526** | 231 | — | — |

> **关键行数签名**：绝对值表 **8207** 行；ΔΔG forward 配对 **7526** 行
> （映射到 v2.2 后 **7523**）。V6.1 PHASE 2 用的是 **8207 行**那一批。

- V6.1 PHASE 2 报告的「全量可用行」= **8207** → 与 v2.2 表一致，
  **不是** ΔΔG 配对的 7523。

## 证据 3 — 数值分布（**辅助证据，非判据**）

> ⚠️ 按 V6.2 §11，数值分布**不作为判定依据**，仅列出以便交叉核对。

| 列 | n | mean | std | min | max | 典型 ΔG 量级? |
|---|---:|---:|---:|---:|---:|---|
| `v2.2.real_dG`（V6.1 实际用的） | 8207 | -10.7166 | 2.8211 | -21.3997 | -1.8982 | ✅ 是（−21~−2） |
| `ddg_pairs.ddg`（真正的 ΔΔG 侧） | 7526 | +1.5215 | 2.3957 | -11.4939 | +12.4388 | 否（≈0 附近） |
| `ddg_pairs.dg_wt`（真正的 ΔΔG 侧） | 7526 | -12.3198 | 2.8503 | -20.8573 | -4.3582 | 是 |
| `ddg_pairs.dg_mut`（真正的 ΔΔG 侧） | 7526 | -10.7983 | 2.9079 | -21.3997 | -1.8982 | 是 |

## 判定

### 判据（严格依据 provenance）

1. **列名来源**：V6.1 的 target 取自 `features/v2.2_proteinea_features.csv`
   的 **`real_dG`** 列 —— 该列是 PPB-Affinity 每一条实验记录的
   **绝对结合自由能**（`ΔG = 0.592·ln KD`），**每条突变体一条记录**。
2. **行宇宙**：V6.1 使用 **8207 行**（v2.2 表全量），
   而真正的 ΔΔG 配对表 forward 只有 **7526 行**（映射后 7523）。
   两者是**不同的行宇宙**。
3. **分组键**：V6.1 用 `pdb_id` 分组（**2171 个唯一值**）作为
   「complex-level」划分 —— 这是**绝对值的 PDB 分组**，
   而不是 ΔΔG 配对所依据的 **231 个复合物**。
4. **没有 ΔG_mut − ΔG_WT 步骤**：V6.1 PHASE 2 的 `load_rows()` 与
   `make_folds()` 中**不存在**任何 WT/mutant 配对、差分或正反向增广逻辑。

> # **TASK_IS_ABSOLUTE_DG**

**结论**：V6.1+ PHASE 2 的 4 档划分实验（以及其 Split Inflation = +0.2663）
是在**绝对 ΔG 回归**任务上做的，**不是**在 ΔΔG 配对任务上做的。
因此：

- **Split Inflation = +0.2663 这个现象本身仍然有效**（它是关于
  「划分严格性 → 可见 Pearson」的陈述，与任务是 ΔG 还是 ΔΔG 无关）；
- **但它不能**被表述为「ΔΔG 任务的 split inflation」；
- **宏观结论**（「该模型全局 Pearson 来自跨复合物区分」）在**绝对值任务**上成立，
  需要 V6.2 PHASE 2 在**真 ΔΔG 配对**上重做才能迁移。

> 依据 V6.2 §11：本判定**不**依赖 `real_dG` 的均值是 −10 还是 +1.5，
> 而依赖**列名来源、行宇宙、分组键、缺少差分步骤**这四条 provenance 证据。

## 文件哈希（可复现性）

| 文件 | 大小 (B) | SHA256 |
|---|---:|---|
| `scripts/v6_1/04_split_inflation.py` | 15142 | `acb1a415e6f7e8fb8c635f518e1498d2…` |
| `features/v2.2_proteinea_features.csv` | 11852893 | `134a406175c975541ae40aed37c3f824…` |
| `data/v4.0/ddg_pairs.csv` | 2238640 | `e43b19d69a8b54efe5cc4b3e647e8765…` |
| `experiments/v6_1/phase2_split_inflation/split_inflation_results.csv` | 969 | `2ec87adb29a1614a903f9bbe6bd51623…` |

- 耗时 1.7s
