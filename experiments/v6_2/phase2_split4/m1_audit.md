# V6.2 PHASE 2.0 — M1 历史实现审计

- 时间: 2026-10-14 08:50:51
- **强制约束**：只允许替换 `split generation` 与 `OOF generation` 层；
  **禁止**修改模型定义、特征定义、超参数（V6.2 §2.0）

## 候选实现与哈希

| 名称 | 路径 | 存在 | 大小 (B) | SHA256 |
|---|---|---|---:|---|
| V5 共享机器（特征装配 + 指标） | `scripts/v5/v5_common.py` | ✅ | 27748 | `91bee0814f66737b50edadd937d96764…` |
| V5 Gate A 复现脚本（C+DGaux_fwd+rev） | `scripts/v5/02_baseline_repro.py` | ✅ | 13749 | `f3b333c3fe88a828385eae945795150d…` |
| V4 原始实现（创新2） | `scripts/v4.0/11_dualhead_dg_ddg.py` | ✅ | 12098 | `991068c7a4bf9342a39c400fc0202356…` |
| V4 原始实现（创新1） | `scripts/v4.0/10_dualhead_symmetric.py` | ✅ | 12363 | `9c7443af6c8410a70e4ec1fafdd15e60…` |
| V4 共享机器 | `scripts/v4.0/v40_common.py` | ✅ | 14213 | `74f2dbb20401ea9bc97c338bae6d2a54…` |

## M1（`C+DGaux_fwd+rev`）语义提取

### DGaux（OOF 辅助任务）语义

- V4 原始: `GroupKFold(n_splits=args.folds)`
- V4 原始: OOF 作为额外一列 → `Daux = np.hstack([Dn, oof[:, None]]).astype(np.float32)`
- V4 原始: OOF 作为额外一列 → `Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)`
- V4 原始: 反向语义 → `Xc = np.vstack([Dn[mtr], -Dn[mtr]])`
- V4 原始: 反向语义 → `Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)`
- V4 原始: 反向语义 → `Xc2 = np.vstack([Daux[mtr], Daux_rev[mtr]])`
- V5 复现: `GroupKFold(n_splits=args.folds)`
- V5 复现: OOF 作为额外一列 → `Daux = np.hstack([Dn, oof[:, None]]).astype(np.float32)`
- V5 复现: OOF 作为额外一列 → `Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)   # S/aux 不取负`
- V5 复现: 反向语义 → `Daux_rev = np.hstack([-Dn, oof[:, None]]).astype(np.float32)   # S/aux 不取负`
- V5 复现: 反向语义 → `Xc = np.vstack([Daux[mtr], Daux_rev[mtr]])`

### 特征定义（`build_pair_features`）

```python
blocks = {"phys_mut": P_mut, "phys_wt": P_wt, "phys_diff": P_mut - P_wt,
              "Z_mut": Z_mut, "Z_wt": Z_wt, "Z_diff": Z_mut - Z_wt}
```
- PCA: `PCA(n_components=min(pca_dim, both.shape[1], len(widx)`
- **PCA 拟合范围**: `pca.fit(both[widx])` （须为 train 行并集）

### 超参数（LightGBM）

```python
LGB_PARAMS = dict(objective="regression", n_estimators=800, learning_rate=0.03,
                  num_leaves=15, max_depth=5, min_child_samples=20, subsample=0.8,
                  subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                  n_jobs=2, verbose=-1)
```

## 复用策略（PHASE 2 采用）

**V6.2 PHASE 2 不重写 M1**，而是：

| 层 | 处理 |
|---|---|
| 特征装配 | **复用** `v5_common.build_pair_features`（PCA 只在各折 train 行拟合）|
| 模型定义 | **复用** `v5_common.LGB_PARAMS`（800 树固定，不用早停）|
| 正/反向增广 | **复用** V5 Gate A 的规则：`y_rev = −y`，**辅助列不取负** |
| ΔG 辅助任务 | **复用** `GroupKFold(pdb_id)` 结构，但 **OOF 分组键按 split 级别分级**（V6.2 §2.4 要求）|
| **split generation** | **替换**（V6.2 的四档 L1–L4）|
| **OOF generation** | **替换**（nested-OOF，继承 outer protocol grouping）|

### ⚠️ 与 V5 的一处**有意偏离**（必须记录）

V5 Gate A 的 DGaux OOF 用 `GroupKFold(groups=groups[tr_i])` 即 **pdb_id** 分组。
V6.2 §2.4 要求 **inner OOF 必须继承 outer protocol 的分组级别**：

- L1 inner → 按 biological pair / duplicate 分组
- L2 inner → 按 complex（pdb_id）分组
- L3 inner → 按 seq40 cluster 分组
- L4 inner → 按 family 分组

因此 V6.2 的 M1 在 **OOF generation layer 被替换**，
这**符合** V6.2 §2.0 的允许范围（允许替换 split + OOF 层）。
同时额外强制 **state_id**（`pdb|变异序列` 哈希）与 outer group 的**双重分组**。

### 实施状态

```
M1_implementation_audited            = True
M1_uses_historical_script            = True
M1_only_split_layer_replaced         = True   # split + OOF，符合允许范围
M1_implementation_not_found          = False
```

- 耗时 —
