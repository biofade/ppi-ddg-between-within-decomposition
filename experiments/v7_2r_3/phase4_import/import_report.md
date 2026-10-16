# V7.2R-3 PHASE 4 — PPI-context 特征导入

- ZIP: `F:\dsh\downloads\v7_2r_3\V7_2R_3_PPI_CONTEXT_FEATURES.zip`
- 实际 SHA256: `b7e9c039e9a8fdc29d874193576d9c0bcfcdab317158945c9e86e395393087e8`
- features npz: **132**
  - 40/132
  - 80/132
  - 120/132

## 4.0 Kaggle provenance

```json
{
  "project": "v7_2r_3_ppi",
  "representation_version": "v7_2r_3_ppi_context_v1",
  "model_version": "esm2_t33_650M_UR50D|full_complex",
  "esm_model": "facebook/esm2_t33_650M_UR50D",
  "esm_hidden": 1280,
  "device": "cuda:0",
  "dtype_inference": "float16",
  "dtype_storage": "float16",
  "representation_level": "MUTATION_SITE_PPI",
  "complex_construction": "ALL_PROTEIN_CHAINS (partners retained)",
  "labels_used": false,
  "supervised_heads_used": false,
  "forbidden_artifacts_used": [],
  "generated_utc": "2026-10-14T14:19:58Z"
}
```

## 4.1 导入结果

| 项 | 值 |
|---|---|
| per-target npz | **132** |
| rows | **3452** |
| unique row_id | **3452** |
| targets | **132** |
| families | **70** |
| dim(r_wt) | **1280** |
| dim(r_mut) | **1280** |
| dim(delta_r) | **1280** |
| 嵌入链数中位 | **3** |
| 非预期 row_id | **0** 文件 → ✅ 无 |

## 4.2 伙伴 sanity（Kaggle 侧）

- 测试 target 数：**10**
- `D_true` 中位 = **4.0289**
- `D_wrong` 中位 = **2.6557**
- **PPI_CONDITIONING = CONFIRMED ✅**

- 耗时 5s
