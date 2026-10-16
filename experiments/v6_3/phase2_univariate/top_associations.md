# V6.3 PHASE 2 — Top associations

Primary outcome：`abs_mean_error`（target-level，N=203）

## DEPLOYABLE

| 特征 | ρ | p | BH q |
|---|---:|---:|---:|
| `signed_mean_error_M0` | -0.4774 | 5.97e-13 | 1.54e-11 |
| `struct_chainpair_contact_density` | +0.2230 | 0.00138 | 0.0146 |
| `struct_contacts_per_res_total` | +0.1845 | 0.00842 | 0.0723 |
| `struct_n_contact_pairs` | +0.1817 | 0.00946 | 0.075 |
| `pred_mean_ddg` | -0.1789 | 0.0106 | 0.0783 |
| `struct_n_contact_res_pairs` | +0.1757 | 0.0121 | 0.0834 |
| `mutations_per_site` | +0.1691 | 0.0159 | 0.102 |
| `pred_mean_M0` | -0.1669 | 0.0173 | 0.102 |
| `struct_weakest_chainpair_dsasa` | +0.1664 | 0.0188 | 0.102 |
| `struct_interface_residue_ratio` | +0.1654 | 0.0183 | 0.102 |
| `max_mutations_per_site` | +0.1597 | 0.0229 | 0.118 |
| `site_entropy` | -0.1553 | 0.027 | 0.132 |

## LABEL_DERIVED

| 特征 | ρ | p | BH q |
|---|---:|---:|---:|
| `squared_mean_error` | +1.0000 | 0 | 0 |
| `extreme_fraction` | +0.5940 | 9.55e-21 | 4.92e-19 |
| `signed_mean_error` | -0.4813 | 3.6e-13 | 1.24e-11 |
| `true_mean_ddg` | +0.4130 | 9.11e-10 | 1.88e-08 |
| `ddg_std` | +0.3064 | 8.72e-06 | 0.00015 |
| `ddg_IQR` | +0.2668 | 0.000119 | 0.00176 |
| `ddg_range` | +0.2319 | 0.00087 | 0.0112 |
| `deleterious_fraction` | +0.2226 | 0.00142 | 0.0146 |
| `beneficial_fraction` | -0.1890 | 0.00691 | 0.0647 |
| `within_kendall` | -0.0772 | 0.315 | 0.541 |
| `within_spearman` | -0.0717 | 0.352 | 0.57 |

## 关键读法

- `SIZE_HAS_INDEPENDENT_ASSOCIATION` 的前提是 `log_n_mutations` 在**多变量控制后**仍保持强关联（见 PHASE 3）。
- 单变量显著**不等于**因果；本项目**没有** causal design。
