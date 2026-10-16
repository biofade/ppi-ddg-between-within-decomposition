# V8 模型池 Provenance 审计（审稿人可复核）

> `rho(Overall r, Macro rho)` 在本次工作中被报出过多个不同值。
> 根因是**模型池组成不同**。本表把每个池的定义、成员构成与来源脚本全部固定下来。

| pool_id | 定义 | 来源脚本 | n | PERMISSIVE | STRICT | OTHER | **rho(r,macro)** | rho(r,eff) | rho(macro,eff) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| `P1_U_CANONICAL_all` | all U_CANONICAL scoreboard rows | `03_model_registry / 06_capability_metrics` | 160 | 7 | 139 | 14 | **+0.1837** | -0.0071 | +0.1878 |
| `P2_U_CANONICAL_nrows4473` | U_CANONICAL AND n_rows==4473 (V6.2 + V7.1 families only) | `08_rank_reversal` | 67 | 7 | 55 | 5 | **-0.2685** | +0.3029 | +0.6365 |
| `P3_P2_minus_V6_2` | P2 with V6_2_split4 removed | `decomposition of P2` | 16 | 0 | 11 | 5 | **+0.1861** | +0.1987 | +0.9730 |
| `P4_U_CANONICAL_minus_V6_2` | U_CANONICAL with V6_2_split4 removed | `43_fix_fig1_and_scenarios (fig1 panel b)` | 109 | 0 | 95 | 14 | **+0.4088** | +0.3085 | +0.0933 |
| `P5_pool_composition_ALL` | all scoreboards incl. mas-align artifacts, 10 artifacts | `35_pool_composition_artifact` | 178 | 37 | 133 | 8 | **+0.2789** | +0.2722 | +0.4555 |
| `P5a_pool_composition_PERMISSIVE` | PERMISSIVE subset | `35_pool_composition_artifact` | 37 | 37 | 0 | 0 | **+0.6262** | +0.8088 | +0.8521 |
| `P5b_pool_composition_STRICT` | STRICT subset | `35_pool_composition_artifact` | 133 | 0 | 133 | 0 | **+0.2214** | -0.0289 | +0.4038 |
| `P6_HERMES_variants` | 7 HERMES/proteinMPNN variants (2 families) | `39/40/41 external scripts` | 7 | 0 | 0 | 7 | **+0.6786** | +0.7143 | +0.3214 |
| `P7_EXT_CANONICAL` | 5 CATH-ddG baselines on canonical rows | `39/40/41 external scripts` | 5 | 0 | 0 | 5 | **+0.7000** | +0.7000 | +0.3000 |
| `P8_EXT_HELDOUT_CATH` | 5 CATH-ddG baselines on official held-out CATH | `39/40/41 external scripts` | 5 | 0 | 0 | 5 | **+0.7000** | +0.4000 | +0.9000 |

**观测范围：-0.2685 .. +0.7000**

## 使用规则（写入论文方法部分）

1. 引用任何秩相关时**必须**同时给出 `pool_id` 与其 provenance 构成；
2. **不得**把任一单值表述为 benchmark 的性质；
3. 报告时给出**池间范围**而非点估计；
4. 最小池（n<20）的秩相关**不应**作为结论依据（有效样本量受伪重复削弱）。
