# V8 Scenario-Specific Recommendations（§41/§42）

> **不使用通用复合分**（§41）。每个场景只按其**任务相关**指标排序。

| 场景 | 说明 | 主排序指标 | 首选模型 | 该指标值 | 按 Overall r 会选 | 场景内排名 | 在 top-3? |
|---|---|---|---|---:|---|---:|:--:|
| `A_quantitative_ddG` | 定量 ΔΔG 预测（需要绝对数值可信） | `rmse` | `L1_random_A1_base_set` | 1.5464 | `L1_random_A1_base_set` | 1 | ✅ |
| `B_library_ranking` | 突变文库排序（只要序对，不要绝对值） | `macro_spearman` | `R1_Base_F0_s42` | 0.3388 | `L1_random_A1_base_set` | 54 | ❌ |
| `C_wetlab_shortlist` | 湿实验短名单（Top-K 命中） | `U_norm` | `L1_random_A3_base_local` | 0.8746 | `L1_random_A1_base_set` | 10 | ❌ |

## 核心问题（§41）

- 三场景中，按 Overall Pearson 会给出**不同首选模型**的个数 = **2 / 3**
- 按 Overall Pearson 选的模型进入该场景 top-3 的个数 = **1 / 3**

> ⚠️ 场景 C 的效用轴已证明在数据层面缺乏功效（见 `AXIS_C_DIAGNOSIS.md`），
> 因此场景 C 的推荐**不可作为强结论**，只作记录。

