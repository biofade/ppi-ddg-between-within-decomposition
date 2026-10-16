# δ-sensitive Concordance（预注册 δ ∈ {0, 0.2, 0.5}）

**规则（写死，不调参）**：真实 `|y_i − y_j| ≤ δ` → 不比较；
prediction tie → 记 0.5；无 eligible pair 的靶点记 NA，不进 macro mean。

| 档 | 模型 | δ | MacroAcc(δ) | n_eligible |
|---|---|---:|---:|---:|
| L1_random | M0 | 0.0 | 0.5097 | 200 |
| L1_random | M0 | 0.2 | 0.5123 | 199 |
| L1_random | M0 | 0.5 | 0.5223 | 194 |
| L1_random | M1 | 0.0 | 0.5106 | 200 |
| L1_random | M1 | 0.2 | 0.5126 | 199 |
| L1_random | M1 | 0.5 | 0.5182 | 194 |
| L4_family | M0 | 0.0 | 0.5023 | 200 |
| L4_family | M0 | 0.2 | 0.4986 | 199 |
| L4_family | M0 | 0.5 | 0.5021 | 194 |
| L4_family | M1 | 0.0 | 0.5040 | 200 |
| L4_family | M1 | 0.2 | 0.5009 | 199 |
| L4_family | M1 | 0.5 | 0.5014 | 194 |
