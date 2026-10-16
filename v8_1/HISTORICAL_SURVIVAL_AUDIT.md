# V8.1 历史头条存活审计（§7.1/§8.2/§14）

> 规则：任何历史头条若在**匹配 + 去伪重复 + provenance-clean** 后不存活，
> 必须降级为历史/补充材料，或删除。

## 汇总

| 历史头条 | 判定 | 是否仍可辩护 |
|---|---|:--:|
| S1 | SURVIVES_ONLY_AS_POOL_DEPENDENT_RANGE — 单值正交/-0.26 不成立；跨池范围包含强正值(+0.82)，因此「近乎正交」**不能作为普遍结论** | ❌ |
| S2 | DOES_NOT_SURVIVE_IN_PROVENANCE_CLEAN_POOL | ❌ |
| S3 | DOES_NOT_QUALIFY_AS_CAUSAL_EVIDENCE | ❌ |
| S4 | NOT_PRESENT_IN_PRIMARY_POOL | ❌ |

## S1 详情：「近乎正交 / ρ=−0.26」

跨 10 个匹配 stratum 的 `rho(Overall r, Macro rho)` 范围 = **-0.2685 .. +0.8242**（负值池 3，正值池 7）。

**⇒ 单值不成立**；且范围上界为**强正值 (+0.82)**，
所以「全局与靶内指标近乎正交」**不能作为普遍结论**。
可辩护的表述见 `FINAL_DEFENSIBLE_CLAIMS.md`。

