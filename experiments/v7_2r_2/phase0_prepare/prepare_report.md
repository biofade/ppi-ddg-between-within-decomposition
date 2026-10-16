# V7.2R-2 PHASE 0 — 构建无标签 Kaggle 推理包

- 时间: 2026-10-14 20:29:57

- canonical rows: **4473** / targets **214** / families **108**

- 限定为 Track A 已覆盖行（保证 U_AB 可比）：**4233** 行 / **174** targets / **91** families

- 可用变体（WT 校验通过）: **4233** 行；跳过（无可用位点）: 0

- 链序列表: **233** 条 (target,chain) 记录

## 0.1 无标签保证（强制检查）

| 检查 | 结果 |
|---|---|
| 包内列名含标签类词 | **0** 处 → ✅ 无 |
| 未包含原始来源文件（含 `ddG`） | ✅ 仅导出剥离后的推理表 |
| 未包含原始来源文件（含 `KD(M)`） | ✅ 仅导出剥离后的推理表 |
| 未包含原始来源文件（含 `affinity`） | ✅ 仅导出剥离后的推理表 |
| 未包含 V7.1/V7.2 预测 | ✅ |
| 未包含 bootstrap / outer 性能 | ✅ |

## 0.2 包内容

| file | rows | cols | bytes | sha256(16) | 含标签 |
|---|---:|---:|---:|---|---|
| `canonical_rows.csv` | 4233 | 3 | 63,002 | `0416bec39ad8d054` | **NO** |
| `variants.csv` | 4233 | 6 | 118,530 | `b0e01bf426819ac5` | **NO** |
| `mutation_sites.csv` | 7280 | 10 | 219,579 | `f76ef23b8f0018f4` | **NO** |
| `target_chains.csv` | 233 | 6 | 55,553 | `b771bca6ae2033d4` | **NO** |
| `variant_reconstruction.csv` | 4233 | 5 | 101,515 | `3c94daad8246baca` | **NO** |
| `PACKAGE_README.md` | 0 | 0 | 1,451 | `768dabf35aa94910` | **NO** |

- **ZIP**: `kaggle_input_package.zip`（0.13 MB）
- **SHA256**: `5bf0a25005f675db2f7360139bd8de87d14609011a704302565e6f6402bf537c`

- 耗时 23s
