# V7.2R-1 PHASE 5 — Track B 可行性 + 计算可行性

## 5.1 Track B 的**目标**形态（§10）

```
ReCLIP-style residue-context feature construction
  +  generic leakage-safe ESM / MINT representations
  +  OUR ranking head trained inside nested family CV
```

**明确不是**：`ReCLIP` 预训练的 mutation 分类器。

## 5.2 能否在不使用其监督 mutation head 的前提下抽出所需表示？

| 所需表示 | 公共 pipeline 是否提供 | 依据 |
|---|---|---|
| site-local residue representation | **YES** | 特征构建器 `scripts/*/ReCLIP/*_prediction_save.py` 落 `Feature_cache/`；ReCLIP 本身是**残基为中心**的框架 |
| partner-conditioned representation | **YES** | README：*「which interaction partner residues are most informative」*；特征为 residue×partner 上下文 |
| MINT partner representation | **YES** | MINT `MINTWrapper(chains, chain_ids)` 接受**双链**输入，输出链级 (N,1280) 或 `sep_chains=True` (N,2560) |
| 是否需要其 mutation XGBoost head | **NO** | README：*「Release artifacts are hosted separately on Hugging Face」*，且 *「trained task-specific classifier heads are external artifacts」* |

```
TRACK_B_FEASIBLE = YES
```

**但有两个硬性前置条件**（见 5.3 / 5.4）。

## 5.3 计算可行性（§11 要求的估算）

**实测 artifact 大小（HTTP HEAD，未下载）**

| artifact | 许可 | 大小 | 用途 |
|---|---|---:|---|
| MINT checkpoint（mint.ckpt） | MIT | **3,253.8 MB** | MINT 主干权重（基于 ESM2-650M） |
| ESM2-650M（model.safetensors） | MIT | **2,609.5 MB** | MINT/ReCLIP 特征构建所需的底层 PLM |
| Bernett MLP head（bernett_mlp.pth） | MIT | **26.2 MB** | MINT 提供的**下游任务 head**（预测两序列是否结合） |
| ESM2-8M（model.safetensors） | MIT | **31.4 MB** | Track A 已在用（已缓存于 F 盘） |

- **Track B 必需两件合计 = 3253.8 + 2609.5 = 5,863.3 MB = 5.86 GB**
- §29 的自动允许上限 = **2 GB** → **超出 3.86 GB**

```
DEFERRED_LARGE_ARTIFACT = YES
```

**本机资源核算**

| 项 | 值 |
|---|---|
| OS | Windows 10（ReCLIP 参考环境为 **64-bit Linux**，未在 Windows 上测试） |
| GPU | **无 CUDA**（`torch.cuda.is_available() = False`） |
| 参考 GPU | README 实测于 **A100 40GB / T4** |
| RAM | **7.55 GB 总量** |
| ESM2-650M 权重磁盘 | ≈2.61 GB |
| MINT 权重磁盘 | ≈3.25 GB |
| 峰值 RAM（单模型 fp32 权重 + 激活） | ≈3–6 GB |
| **两模型同时驻留** | **≈6–9 GB → 超过 7.55 GB 物理内存** |

**序列规模（U_CONTRACT）**

| 项 | 估算 |
|---|---|
| 行数 | 4473（+ WT 复本） |
| 唯一 WT/MUT 链序列 | 数千条 |
| 最长链 | 3128 aa（需截断） |
| 特征缓存（650M 级、残基级） | **数十 GB 量级** |

⇒ 结论：

```
TRACK_B_GPU_REQUIRED = YES
TRACK_B_READY_FOR_GPU = YES
```

**理由**：① 必需 artifact 合计 5.86 GB > 2 GB 上限；
② 无 CUDA；③ ESM2-650M + MINT 同时驻留超出 7.55 GB 物理内存；
④ ReCLIP 仅实测于 Linux/A100。
**依 §11「Do not substitute a weaker model silently」，
不以更弱模型替代，而是明确标记需要 GPU。**

## 5.4 Track B 启动清单（GPU 就绪后可直接执行）

1. 在 **Linux + NVIDIA GPU**（≥24 GB 显存建议）环境：
   `git clone https://github.com/SiweiLab/ReCLIP.git`；
   `git clone https://github.com/VarunUllanat/mint.git mint`
2. 下载 `mint.ckpt`（3.25 GB, MIT）与 `facebook/esm2_t33_650M_UR50D`（2.61 GB, MIT）
3. **只运行 feature builder**，落 `Feature_cache/`；
   **绝不加载** `classifier_heads/mutation/reclip_mutation_xgb_*.pkl`
4. 用我们自己的 `z` ranker，在 `OUTER_FOLDS_V7_2` + inner family CV 下训练
5. 依 §9，论文须写明：family-OOD 指**任务监督训练**，
   **不保证**基础模型预训练未见这些家族

- 耗时 0.0 分钟
