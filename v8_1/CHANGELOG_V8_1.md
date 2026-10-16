# CHANGELOG — V8.1（清理与冻结）

## A. 解析性修正（analytic corrections）

1. **指标方向登记表**（`METRIC_DIRECTION_REGISTRY.csv`，29 指标）首次建立。
2. **方向 bug 修正**：`12_decision_alignment.py` 用未取负的 `rmse` 做秩相关 ⇒
   校准轴汇总作废，V8.1 用 `-RMSE` 重算。受影响：
   `DECISION_ALIGNMENT_REPORT.md` 的相关句、`T3`/`Fig1` 的校准轴。
3. **行集分层**：12 个 artifact 分属 **7 个不同 `rowset_hash`**
   ⇒ 「U_CANONICAL」被确证为**命名空间而非匹配宇宙**，禁止跨行集比较。
4. **去伪重复三层**：artifact（描述性）/ configuration / method-family（LEVEL C）。
   `MODEL_FAMILY_REGISTRY.csv` 记录 202 键、73 configuration、
   63 个代表性 LEVEL-C 模型（**performance-blind** 选取）。
5. **provenance 分类**：GREEN=165、
   RED=37（依 `split_regime`，非仅凭名称）。
6. **泄漏规则改写**：`r>0.98` 由「自动排除」改为 `LEAKAGE_ALERT=TRUE` +
   须有身份/来源证据才排除；共记录 23 条。

## B. 解释性修正（interpretation corrections）

1. `C_B`/`C_W` 统一称 **covariance-signal contributions**，禁止称 correlation / skill / 因果机制。
2. 删除「model memorized PDB identity」类表述，改为「consistent with reliance on
   target-specific information that does not transfer across families」。
3. AXIS C 的「data floor」改为 **candidate-set geometry**（§11.3）。
4. 删除「n ≥ 50 因为区分度随 n 增加」——该推断在数学上不成立（见 §11.2）。

## C. 图表重生成（figure/table regeneration）

* 全部最终图与头条表将由 `v8_1/FROZEN_RESULTS.json` 生成（§16/§17/§21）。
* `fig1` 的 160-artifact 热图降级为 **Supplement**（artifact-level, not independent）。

## D. 上游模型预测：**未做任何改动**

* `UPSTREAM_HASH_CHANGE = 0`（见 `FROZEN_INPUT_MANIFEST.csv`）；
* 未训练、未重跑 PLM、未跑 OpenMM/V7.4R2、未新建架构、未新建物理特征；
* `GPU_USED = NO`。
