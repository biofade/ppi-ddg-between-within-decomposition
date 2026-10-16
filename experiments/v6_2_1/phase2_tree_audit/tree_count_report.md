# V6.2.1 PHASE 2 — LightGBM 树数 Provenance 审计

- 时间: 2026-10-14 13:17:09
- **判据**：只有真实训练时使用的值才是最终值；不根据文本猜

## 证据 1 — 源码中的树数定义

| 来源 | 文件 | 行内容 | SHA256(前16) |
|---|---|---|---|
| V6.2 split4 训练脚本 | `scripts/v6_2/02_split4_ddg.py` | TREES = 800   # 可被 --trees 覆盖（所有档一致）<br>m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": TREES,<br>m0 = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": TREES,<br>"n_estimators": TREES, | `6a391a9c08f64a37` |
| V5 共享超参（LightGBM 配置源） | `scripts/v5/v5_common.py` | LGB_PARAMS = dict(objective="regression", n_estimators=800, learning_rate=0.03, | `91bee0814f66737b` |
| V5 Gate A 复现 | `scripts/v5/02_baseline_repro.py` | ap.add_argument("--trees", type=int, default=800)<br>m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,<br>m_full = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees,<br>m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": args.trees, | `f3b333c3fe88a828` |

## 证据 2 — 实际运行时使用的命令行（**权威来源**）

V6.2 的 PHASE 2 运行方式（记录于当时的 shell 调用）：

```
python -u scripts/v6_2/02_split4_ddg.py --levels L1_random L2_complex L3_seq40 L4_family \
    --seeds 42 43 44 --trees 300 --max-minutes 120
```

**`--trees 300` 是显式传入的**，且 `main()` 中：

```python
ap.add_argument("--trees", type=int, default=800, ...)
args = ap.parse_args()
global TREES
TREES = int(args.trees)        # ← 用命令行值覆盖默认 800
```

⇒ **实际训练树数 = 300**（所有档、所有 seed、所有模型统一）。

## 证据 3 — 运行时日志

- `- 模型定义/超参与 V5 Gate A **完全一致**（800 树，PCA=128，无早停）`

- 日志 mtime: **2026-10-14 10:49:17**
- 日志 SHA256(前16): `033a85c42e9e93dd`

## 证据 4 — 预测文件批次一致性

- `split4_predictions.npz`：**44** 个数组，**36** 个 (level, model, seed) 预测
- mtime: **2026-10-14 11:50:39**
- SHA256(前16): `ba711ea1b349836a`
- **所有 36 个预测在同一批内产出**（L3 因 fold 修正于稍后单独重跑并合并 —— 见 PHASE 2 记录的合并步骤）

> 关键：**不存在「部分 300 / 部分 800」的情形**：
> 唯一一次全量训练使用 `--trees 300`；L3 的单独重跑同样使用 `--trees 300`。

## 证据 5 — 「800」出现在哪里

| 文件 | 「800」出现处 | 性质 |
|---|---|---|
| `scripts/v6_2/02_split4_ddg.py` | `TREES = 800`（默认值）/ `default=800` | **源码默认值**（被 `--trees 300` 覆盖） |
| `experiments/v6_2/V6_2_FINAL_SUMMARY.md` | PHASE 6 段「V6.2 的 M1 已完成 800 树复现」 | **报告文本错误** |
| `experiments/v6_2/PROGRESS.md` | 「800 树 → 本次统一 300 树」 | 正确说明 |

> **定位结论**：`800` 只出现在 (a) 源码**默认值**、(b) 一处**报告文本笔误**。
> **没有任何证据**表明实际训练用了 800 树。

## 判据与结果

| model | split | seed | reported_value | actual_value | source_of_truth |
|---|---|---:|---:|---:|---|
| M0 | L1_random | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L1_random | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L1_random | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L1_random | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L1_random | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L1_random | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L1_random | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L1_random | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L1_random | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L2_complex | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L2_complex | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L2_complex | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L2_complex | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L2_complex | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L2_complex | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L2_complex | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L2_complex | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L2_complex | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L3_seq40 | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L3_seq40 | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L3_seq40 | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L3_seq40 | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L3_seq40 | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L3_seq40 | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L3_seq40 | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L3_seq40 | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L3_seq40 | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L4_family | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L4_family | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M0 | L4_family | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L4_family | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L4_family | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M1 | L4_family | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L4_family | 42 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L4_family | 43 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|
| M2 | L4_family | 44 | 300 | **300** | 命令行 `--trees 300`（覆盖源码默认 800）|

### 标记：**`ACTUAL_300_CONSISTENT`**

- 全部档 / 模型 / seed 的**实际**树数均为 **300**（统一值）
- 因此 **split inflation 比较仍然公平**（各档同值）
- 「800」仅为 (a) 源码默认值、(b) 一处报告文本笔误
- **无需重训**

### 建议的文本修正（V6.2.1 不改 V6.2 原文件，仅在此记录）

> `V6_2_FINAL_SUMMARY.md` PHASE 6 段中的「**已完成 800 树复现**」
> 应改为「**使用 300 树（统一值）**」。

- 耗时 0s
