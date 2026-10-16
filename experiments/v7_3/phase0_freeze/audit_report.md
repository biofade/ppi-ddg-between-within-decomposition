# V7.3 PHASE 0 — Provenance 冻结 + 软件/工具审计

- 时间: 2026-10-15 00:15:24

## 0.1 复用输入（**17** 项，已冻结）

| name | size | sha256(40) |
|---|---:|---|
| `canonical_target_labels` | 66,663 | `b0f43c069c0baf32e979beb4aff4f7a06a69b415` |
| `canonical_matched_set` | 754,492 | `8d209eea08cc6058fb235b059c9e3725fcbc9ba8` |
| `v711_protocol` | 983 | `4f8aa6b5729a8a632b5ffdaca39ae124f242f318` |
| `v711_fold_manifest` | 526 | `818ccd5d792fba6b3e621b571b85d2b9731a90c5` |
| `v72_eval_contract` | 4,394 | `1c3abcbf2a761db69072e75b551482833735898b` |
| `mutation_sites` | 706,966 | `ca23d3b42eadf078dce62b0af67d01895a3d6c57` |
| `v70_block_L` | 3,414,403 | `dfd0886c54cd96f21e1f646a6e74b36741a79472` |
| `v70_block_S` | 54,408 | `35b0bb7df68dbce9585c997d74126d0122b15abe` |
| `v70_block_G` | 164,572 | `9989a4e8e37732d00a4081bb5419ec2c557c6e9f` |
| `target_mechanism_table` | 285,776 | `b2b241cd2afc0dacb1c05692a637699cbec4d96d` |
| `v5_common` | 27,748 | `91bee0814f66737b50edadd937d967644bbf3ec7` |
| `ddg_split` | 82,187 | `92d4e8d78f0778ff82f12896aeedc9607e5e8cf7` |
| `ddg_pairs` | 2,238,640 | `e43b19d69a8b54efe5cc4b3e647e87657106e86c` |
| `v72r1_full_preds` | 393,550 | `fd3b528822aabd7718b0479dee76585f98ca355f` |
| `v72r3_ppi_preds` | 363,992 | `0b3f20e68309c6cc2da3ca462a7c5994efbd0203` |
| `trackA_per_target(aggregate)` | 174 | `abb8ee3744aec1f7b8fffe423093a134e02d6f80` |
| `pdb_dir(count+bytes)` | 3,031 | `3031_3315016935` |

- **`UPSTREAM_HASH_CHANGE = 0`**（本 PHASE 只读，未修改任何上游文件）
- 缺失输入：**无** ✅
- PDB 文件数：**3031**（§6 记载 230/230 覆盖）

## 0.2 物理/结构工具审计（§8）

| tool | available | version | license | primary | secondary |
|---|---|---|---|---|---|
| NumPy | ✅ | 2.4.6 | BSD-3 | ✅ | ✅ |
| SciPy | ✅ | 1.17.1 | BSD-3 | ✅ | ✅ |
| scikit-learn | ✅ | 1.9.1 | BSD-3 | ✅ | ✅ |
| LightGBM | ✅ | 4.7.0 | MIT | ✅ | ✅ |
| BioPython | ✅ | 1.88 | Biopython License (BSD-like) | ✅ | ✅ |
| FreeSASA (python bindings) | ✅ | ? | MIT | ✅ | ✅ |
| OpenMM | ❌ |  | MIT/LGPL | ❌ | ❌ |
| PDBFixer | ❌ |  | MIT | ❌ | ❌ |
| ProDy | ❌ |  | MIT | ❌ | ❌ |
| MDAnalysis | ❌ |  | GPL-2.0+ | ❌ | ❌ |
| MDTraj | ❌ |  | LGPL-2.1 | ❌ | ❌ |
| FoldX @ F:\dsh\downloads\foldx | ❌ |  | **restricted (需商业许可)** | ❌ | ❌ |
| FoldX @ C:\Program Files\foldx | ❌ |  | **restricted (需商业许可)** | ❌ | ❌ |
| FoldX @ F:\DSH_Galatea\Affinity_Intelligence\tools\foldx | ❌ |  | **restricted (需商业许可)** | ❌ | ❌ |

### 结论

- **SASA 引擎：`freesasa`（MIT）可用** → 依 §20 优先使用，
  并以 BioPython ShrakeRupley 作交叉校验
- **无任何合法的外部能量引擎**（FoldX/OpenMM/EvoEF2 均不可用或未验证许可）
  → 依 §8/§40：**`EXTERNAL_ENERGY_ENGINE = NOT_RUN`**，
  且 V7.3 主实验**不依赖**它（§8 明确要求可独立运行）
- 主实验全部可用：NumPy / SciPy(cKDTree) / LightGBM / sklearn

## 0.3 BLOCK L 按 `row_id` join 验证（§58）

| 项 | 值 |
|---|---:|
| `canonical_target_labels` rows | 4473 |
| BLOCK L rows / cols | 4473 / 121 |
| join 后 rows | **4473** |
| target 缺失 | **0** |
| `BLOCK L.pdb_id` vs canonical 一致率（**仅参考**） | **94.21%** |

- ⇒ **join 干净**（4473 行 / 0 缺失）
- ⚠️ 一致率仅 94.21% 是**已知的 V7.1.1 namespace 问题**：
  `BLOCK L` 的 `pdb_id` 源自 `mutation_sites`，**不用于 target 身份**。
  ⇒ 依 §3/§58：**唯一权威 namespace = `canonical_target_labels.csv`**，
  一律按 **`row_id`（= `pair_row`）** 连接。

- BLOCK L 特征列（前 20 个）：`['partner_mean_kd_r6', 'partner_mean_kd_r6__max', 'partner_frac_aromatic_r8', 'partner_frac_aromatic_r8__max', 'partner_frac_hydrophobic_r8', 'partner_frac_hydrophobic_r8__max', 'mut_n_samechain_contacts_r8', 'mut_n_samechain_contacts_r8__max', 'mut_n_contacts_r6', 'mut_n_contacts_r6__max', 'partner_mean_kd_r10', 'partner_mean_kd_r10__max', 'mut_n_partner_contacts_r6', 'mut_n_partner_contacts_r6__max', 'partner_mean_vol_r6', 'partner_mean_vol_r6__max', 'partner_frac_charged_r10', 'partner_frac_charged_r10__max', 'mut_is_iface', 'mut_is_iface__max']`
- BLOCK L 数值列数：**119**

- 耗时 7s
