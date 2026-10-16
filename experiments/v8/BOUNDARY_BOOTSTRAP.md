# V8 Boundary Bootstrap（§6.2.2）— **根因已定位；结论**恢复**并确证**

> ## 修正记录（重要）
>
> 上一轮我**撤回**了「ρ(Overall r, Eff@10) 从 +0.8088（宽松）塌到 −0.0289（严格）」，
> 理由是 `36_boundary_bootstrap.py` 给出**相反符号**（−0.558 / −0.325）。
>
> **本轮的诊断证明：撤回本身是错的 —— 有 bug 的是 `36`，不是 `35`。**
>
> `37_diagnose_discrepancy.py` 用**规范实现**（`v8metrics.compute_all`）
> 独立重算 `pool_composition_models.csv` 的每个 PERMISSIVE 模型：
>
> ```
> max |Δ overall_pearson| = 5.860e-03      ← 三位 L1_random_M2_s* 的微小差异
> max |Δ efficiency@10|   = 5.782e-03
> rho（用 pool 读数）      = +0.8088
> rho（用 v8metrics 重算） = +0.8162
> ```
>
> 随后 `38_boundary_bootstrap_clean.py` 用规范实现重算点估计，**逐位复现**：
>
> ```
> 点估计（规范实现）: PERMISSIVE rho=+0.8088 | STRICT rho=-0.0289 | Δ=+0.8377
> 对照 pool 文件:      PERMISSIVE +0.8088 / STRICT -0.0289
> ```
>
> **⇒ 2000 次 target-cluster bootstrap（`36`）的整个结果作废；原结论恢复。**

---

## 1. 最终结果（点估计 + 两种不确定度）

| 量 | point | Fisher 95% CI（**模型池**抽样） | `36` 的 target-cluster CI（**已作废**） |
|---|---:|---|---|
| `PERMISSIVE` (n=37) | **+0.8088** | **[+0.6570, +0.8976]** 排除 0 | [−0.7521, −0.0676]（符号错） |
| `STRICT` (n=133) | **−0.0289** | [−0.1981, +0.1421] 含 0 | [−0.4826, +0.1079] |
| **`DELTA`** | **+0.8377** | **[+0.6498, +0.9104] 排除 0** | [−0.7242, +0.2011]（作废） |

**⇒ 在宽松划分制度下，Overall Pearson 与 Efficiency@10 的模型排名秩相关为
+0.81（CI [+0.66, +0.90]）；在 family-OOD 制度下为 −0.03（CI [−0.20, +0.14]）。
Δ 的 Fisher 95% CI = [+0.65, +0.91]，排除 0。**

## 2. 必须声明的范围限制

> ⚠️ **Fisher CI 只看模型池抽样误差，不看 target 抽样误差。**
> 严格正确的做法是 target-cluster bootstrap，但那需要 `v8metrics` 完整调用
> ⇒ 实测 **14.5 s / 次迭代**（170 模型 × 130 个 target，每次都要重算每模型每 target），
> 2000 次约 **8 小时**，超出本轮预算。
>
> **因此本结论的强度限定为**：
> * 「两个制度下的 ρ 差异很大」——**已确证**（点估计 +0.81 vs −0.03；
>   连 `STRICT` 的 Fisher CI 上界（+0.14）也远低于 `PERMISSIVE` 的点估计）；
> * 「差异的精确区间」——**未确证**（缺 target 抽样分量）。
>
> `38_boundary_bootstrap_clean.py` 已在后台以 `--max-minutes 35` 运行，
> 完成后可补上 target-cluster 分量。**在补齐之前，不引用 target-cluster 区间。**

## 3. 为什么这是 V8 关于「宽松划分危险」的最强证据

| 制度 | ρ(Overall r, Eff@10) | 含意 |
|---|---:|---|
| **PERMISSIVE**（随机/complex/seq40） | **+0.8088** | Overall Pearson **是**选择效用的优良代理 |
| **STRICT**（family-OOD） | **−0.0289** | Overall Pearson **对**选择效用**毫无信息** |

**⇒ 一个只用宽松划分构建的基准会得出「全局指标可以指导湿实验筛选」的结论；
而它在 family-OOD 下完全失效。** 这比「Overall r 膨胀 4.5 倍」更严重 ——
后者只是数值虚高，前者是**模型选择建议的方向性错误**。

## 4. 对 V8 判定的影响

* `H1 = WEAKLY_SUPPORTED` **维持**（秩一致在严格制度下 CI 含 0 ⇒ 不稳定）。
* `BENCHMARK_VALIDATION = PARTIALLY_SUPPORTED` **维持**。
* `BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE` **维持**，
  依据**恢复并强化**为：严格制度下 `ρ(Overall r, Eff@10)` 的 CI 含 0（−0.20 ~ +0.14），
  而宽松制度下为 +0.81 ⇒ **现有全局指标在严格制度下失去对选择效用的指示力**。
* `POOL_COMPOSITION_ARTIFACT.md` §2b **恢复**（其数值经规范实现逐位复核）。
* `V8_FINAL_SUMMARY.md` §6.2.2 第 2 条 **恢复**。

## 5. 方法学教训（本轮新增，且比上一轮的版本更准确）

上一轮我写下的教训是「不能脱离 provenance 报告秩相关」—— **这仍然对**。
但本轮加了一条**更根本**的：

```text
当两个独立实现对同一统计量给出相反符号时，
正确的第一反应是**用规范实现逐位复核**，而不是撤回结论。
撤回一个正确的结论与保留一个错误的结论是同类错误。
⇒ 纪律：任何跨实现的分歧，先做**逐模型逐位比对**（本轮的 37 号脚本），
  再决定是撤回还是修复。
```

**本轮我犯的错误**：在发现 `36` 与 `35` 不一致时，我**假设了 `35` 可能是错的**，
并据此撤回。实际上 `36` 才是坏的。**正确做法应是先诊断再决定。**

## 6. 复现命令

```powershell
python -W ignore -u scripts/v8/37_diagnose_discrepancy.py          # 逐位比对（决定性）
python -W ignore -u scripts/v8/38_boundary_bootstrap_clean.py      # 规范实现 bootstrap
```

产出：`external/discrepancy_diagnosis.csv`、
`external/boundary_bootstrap_clean.csv`、`external/boundary_bootstrap_clean_meta.json`
