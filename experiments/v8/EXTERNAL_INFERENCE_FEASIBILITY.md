# V8 外部推理可行性审计（PHASE 13 收尾）

> **结论：在不违反 §9（不静默替换 / 不静默降级）与不动用 GPU 的前提下，
> 外部模型推理在本机**不可行**。这不是"没时间做"，而是**有明确代码级与资源级证据的阻断**。**

---

## 1. PPIformer —— **代码级硬阻断**（且是官方的 TODO）

PPIformer 是三条候选路线中**文档最明确支持 CPU** 的一条
（README 原文：`device = 'cuda' if torch.cuda.is_available() else 'cpu'`）。
因此我先选它做可行性验证。

克隆后逐文件扫描，发现**一处无条件 CUDA 调用**：

```python
# ppiformer/data/transforms.py:467-472
def __call__(self, data: Data) -> Data:
    src, dst = data.edge_index[0], data.edge_index[1]
    rel_pos = data.x[dst] - data.x[src]
    # TODO if not torch.cuda.is_available():          <- 作者自己标注的未实现分支
    data.bases = {k: v.cpu() for k, v in get_basis(rel_pos.cuda(), **self.bases_kwargs).items()}
    return data
```

**第 471 行 `rel_pos.cuda()` 是硬编码的**，违反其 README 承诺。
在纯 CPU 机器上，该 transform 会在任何推理开始**之前**就抛
`AssertionError: Torch not compiled with CUDA enabled`。

**⇒ 这是决定性的**：README 声明的 CPU 路径**实际不可用**，
而作者用 `# TODO` 明确承认了这一点。

（其余 13 处 CUDA 命中为注释、文档串或 `test_cuda_installation.py`，非阻断项。）

---

## 2. 依赖与资源阻断（独立于第 1 条）

| 项 | PPIformer 要求 | 本机 | 判定 |
|---|---|---|---|
| `torch` | **`>=2.7`** | **2.5.1+cpu** | ❌ 不满足 |
| `torch_geometric` | `2.3.1` | **未安装** | ❌ |
| `pytorch_lightning` | `2.0.8` | **未安装** | ❌ |
| `graphein` | `1.7.7` | 未安装 | ❌ |
| `equiformer-pytorch` | **fork @ git** | 未安装 | ❌ |
| `mutils` / `ppiref` | git | 未安装 | ❌ |
| 权重 `weights.zip` | **510.0 MB** | — | 需下载 |
| 数据缓存 `.pyg_dataset_cache.zip` | **516.2 MB** | — | 需下载 |
| 可用 RAM | — | **1.9 GB / 7.5 GB** | ⚠️ 装完依赖后余量极低 |

**装 torch>=2.7 会破坏本机现有的 `torch 2.5.1+cpu`** ——
而该 torch 是 V8 全部内部分析的运行时（虽然 V8 的度量引擎只用 numpy/scipy，
但破坏环境会危及后续任何复算）。

**⇒ 为一次**结果不确定**的外部推理去改动工作环境，风险收益比不成立。

---

## 3. 另两条路线的阻断

| 路线 | CPU 声明 | 阻断原因 |
|---|---|---|
| **RDE-Network** | `--device` 默认 `cuda` 但可配 | 权重 30k step，需 `torch_geometric`；**未发布逐突变预测**，必须自跑 |
| **Prompt-DDG** | 代码级 `cuda if available else cpu` | 权重在 Google Drive；**无任意 (结构,突变) 的评分 CLI**（SARS 路径在代码中被注释掉） |

**⇒ 三条路线全部阻断，原因各不相同，且都不是"我可以绕过"的。**

---

## 4. 遵循的纪律

依 §9 与用户长期要求（**"宁可显式失败，也不要静默换物理模型"**）：

* **不**把 `rel_pos.cuda()` 改成 `.cpu()` 后强行运行 ——
  那是**修改上游代码**，会使结果不代表官方实现，
  且违反「**不修改 V1–V7**／不改造上游」的精神（此处虽非 V1–V7，但同类）；
* **不**静默降级到另一个更容易的模型；
* **不**为了凑够外部模型数而放宽 §68 的门槛。

**⇒ 明确记录为阻断，而非静默跳过。**

---

## 5. 对 §68 的最终影响

```text
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE
   证据:
   (a) 官方发布的逐突变预测只有 HERMES 一家（已用完，构成 U_HERMES_CANONICAL）
   (b) 可复用的多模型预测表只有 CATH-ddG 仓库的 5 个基线
       （已用完，构成 U_EXT_CANONICAL 与 U_EXT_HELDOUT_CATH）
   (c) 其余三条路线（PPIformer / RDE-Network / Prompt-DDG）均需自跑推理，
       而 PPIformer 存在**官方 TODO 承认的 CPU 硬阻断**，另两条依赖与权重
       同样不可得或不适用
   ⇒ 两个外部宇宙各 5 个可评价模型（< 8 门槛）⇒ 无法检验能力分化
   ⇒ §68 STRONGLY_SUPPORTED 条件**不成立**（原因：不可检，非阴性）
```

**⇒ 依 §69，V8 的最终判定维持：**

```text
BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE
```

**这是本项目的真实能力边界，不以"再试一次"来绕过。**

---

## 6. 若要突破此边界，需要什么（供后续）

1. **一台有 CUDA GPU 的机器**（PPIformer 路线即通；本机无）；
2. 或**上游修复** `transforms.py:471` 的 `.cuda()`（可作为 issue 反馈给作者）；
3. 或**接受修改上游代码**并明确标注「非官方实现」——
   本轮因纪律原因**不采取**；
4. 或**更多已发布逐突变预测的模型**出现（这是最省力的路径，
   且 V8 的 `U_EXT_*` 构建脚本可直接复用）。

---

## 7. 复现命令

```powershell
git clone --depth 1 https://github.com/anton-bushuiev/PPIformer.git F:\dsh\cache\tmp\PPIformer
# 查看阻断点
Select-String -Path F:\dsh\cache\tmp\PPIformer\ppiformer\data\transforms.py -Pattern 'rel_pos.cuda'
```
