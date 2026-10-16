# BUG 修正记录 — `efficiency_at_10` 反号

**时间**: 2026-10-13 20:39:13
**状态**: `INVALID_DUE_TO_BUG`（旧产物已重命名保留为 `*.INVALID_DUE_TO_BUG`）

## 现象（合成自检证据）

| 合成模型 | 修正前 efficiency_at_10 | 修正后 |
|---|---:|---:|
| 完美 oracle | **−0.6667** | **+1.0000** |
| 完全反向 | **+1.0000** | **−0.9816** |
| 含噪中间模型 | −0.6667 | **+0.7644** |

## 根因

`v5_common.topk_metrics` 中：
1. 每个 complex 的 **random 基线只抽 1 次**，在 n 较小的 complex 上方差极大，
   使分母 `v_or − v_rnd` 可正可负，符号不稳定；
2. oracle 与 random 的语义在符号变换中混用。

## 修正

* 随机基线改为 **200 次抽样取均值**；
* efficiency 统一定义为 `(random − model) / (random − oracle)`
  → 对合理模型恒 ≥ 0；oracle → +1；random → 0（与 V5 §108 语义一致）。

## 影响范围

* **受影响**：`efficiency_at_10`（及 `efficiency_at_10_median`、
  `utility_*`、`regret_at_10`）一列。
* **不受影响**：Pearson / Spearman / macro per-PDB / RMSE / MAE / AUROC /
  P@k / EF@x% / Hit@k —— 这些不依赖该函数，因此 **P2/P3/P4 的主要结论不变**。

## 处置（V5 §123）

1. 保存 bug evidence ✅（本文件）
2. 写 reproduction ✅（上表合成自检）
3. 修 bug ✅
4. **从受影响节点之后全部重新运行** → 见下方重算链
5. 旧结果标 `INVALID_DUE_TO_BUG` ✅
6. 不只更新最终表 ✅（重跑全部下游脚本）
