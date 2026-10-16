# V5.0 — P4 OOF 异常诊断

- Gate A npz 数组: ['test_y', 'test_group', 'aux_dg_mut_oof', 'aux_dg_mut_test', 'C+DGaux_fwd+rev_s42', 'C+DGaux_fwd+rev_s43', 'C+DGaux_fwd+rev_s44']
- Gate A OOF 向量长度: 7523（当前特征表 7523）
- Gate A OOF(train 行) Pearson = **+0.4404**（覆盖 5265/5265）

### Gate A OOF 逐折（按 pdb_id GroupKFold 复算折归属）

- fold 0: n=1053, complexes=23, Pearson **+0.3216**
- fold 1: n=1053, complexes=23, Pearson **+0.5823**
- fold 2: n=1053, complexes=24, Pearson **+0.2887**
- fold 3: n=1053, complexes=24, Pearson **+0.5347**
- fold 4: n=1053, complexes=24, Pearson **+0.4567**

## 关键检查：P4 的 OOF 是否用 aux 列预测了自己

P4 的 fold 内训练矩阵构造：

```python
Xa = np.vstack([Daux[a], np.hstack([-Dn[a], oof_aux[a][:, None]])])
```

其中 `Daux[a] = [Dn[a], oof_aux[a]]`。因此 **训练特征里包含 `oof_aux`**，
而 `oof_aux` 本身是 GroupKFold 的**折外**对该行 ΔG_mut 的预测。

对**折外行 b** 而言，`oof_aux[b]` 来自「未包含 b 的模型」→ 对 b 是合法的；
但 fold 内**训练行 a** 用的 `oof_aux[a]` 也来自「未包含 a 的模型」→ 也合法。
**因此构成上没有泄漏。** 问题在别处。

## 直接复现 P4 fold 0（与 Gate A 同一 fold 划分与种子）

- 本次复算 aux OOF(train) Pearson = **+0.4404**

### P4 式的第 1 折主模型（含 aux 列）

- fold 0: n=1053, complexes=23, Pearson(正向预测) **-0.4234**, Pearson(反向预测) **+0.4110**
    fold0 y 统计: mean +2.590, std 2.577, ΔΔG>0 占 84.3%
    fold0 pred 统计: mean -1.421, std 2.240
- fold 1: n=1053, complexes=23, Pearson(正向预测) **+0.5212**, Pearson(反向预测) **-0.5213**
- fold 2: n=1053, complexes=24, Pearson(正向预测) **+0.2513**, Pearson(反向预测) **-0.3069**
- fold 3: n=1053, complexes=24, Pearson(正向预测) **+0.0015**, Pearson(反向预测) **+0.0424**
- fold 4: n=1053, complexes=24, Pearson(正向预测) **+0.2514**, Pearson(反向预测) **-0.3441**

## 判读

若各折「正向预测」Pearson 的正负号在不同折间翻转，说明 **OOF 逐折 Pearson
在 n≈1000 且组内标签高度偏斜时方差极大**，不能作为校准残差的可靠性指标；
此时 conformal 的分位数估计会被离群折拉偏。这不属于「代码 bug」，
但属于**方法学风险**，必须在报告中如实标注。

