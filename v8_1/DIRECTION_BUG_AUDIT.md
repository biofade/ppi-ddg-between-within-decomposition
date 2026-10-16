# V8.1 方向 bug 审计（§2）

## 检查项：旧脚本在做「跨指标秩一致」时是否把 RMSE 当 higher-is-better

`scripts/v8/12_decision_alignment.py` 中：

```python
spearmanr(g[pr], g[ut])   # pr = 'rmse' 未取负
```

**这是真实的符号错误**：RMSE 是 lower-is-better，直接与效用指标做秩相关会**翻转符号**。

| utility | rho(raw RMSE) 旧行为 | rho(-RMSE) 修正 | rho(Overall Pearson) | 符号翻转 |
|---|---:|---:|---:|:--:|
| `efficiency@10` | -0.1685 | +0.1685 | -0.0071 | **是** |
| `ndcg@10` | -0.2451 | +0.2451 | +0.0953 | **是** |
| `hit@10` | +0.1141 | -0.1141 | -0.2578 | **是** |

## 处置

* `12_decision_alignment.py` 的「校准轴」汇总**作废**（`OBSOLETE`）；
* `DECISION_ALIGNMENT_REPORT.md` 中依赖该汇总的句子**作废**；
* `T3` / `Fig1` 的**校准轴**在 V8.1 中用 `-RMSE` 重算；
* 其余指标（Overall Pearson、Spearman、Macro Spearman、C_B、C_W、效用）方向本就一致，未受影响。

## 未被污染的部分

* `08_rank_reversal.py` 的**排名**计算用了 `LOWER_BETTER={'rmse'}`，排名本身正确；
  但其 `rank_agreement` 矩阵（`spearmanr(Rk[a], Rk[b])`）是对**已排名次**求相关 ⇒ 亦正确。
* 受影响的仅是**直接用原始列值**做相关的地方（`12_decision_alignment.py`）。
