# V7.2 PHASE 0 — 残基级模型选择决策记录（**下载前**记录，§16 要求）

**时间**：2026-10-14

## 1. 本地缓存审计结果

| 位置 | 内容 |
|---|---|
| `C:\Users\KELE\.cache\huggingface\hub` | **不存在**（无 HF 缓存） |
| `F:\dsh\cache` | 存在（空） |
| `F:\dsh\downloads` | 存在（仅历史下载） |
| `transformers` | **4.46.3**（已安装） |
| `torch` | **2.5.1+cpu**（CUDA **不可用**） |
| `esm` / `prody` / `fairscale` | **未安装** |

⇒ **本地没有可用的 residue-level 模型** → 按 §3.15 允许**获取一个小型预训练 PLM**。

## 2. 决策（下载前登记）

| 项 | 值 |
|---|---|
| **model name** | `facebook/esm2_t6_8M_UR50D` |
| **类** | ESM-2 **8M**（§3.17 首选的 fallback 类） |
| **参数数** | ~8M |
| **embedding dim** | **320** |
| **层数** | 6 |
| **expected size（权重合计）** | **≈ 31 MB**（`model.safetensors` ~31 MB，float32） |
| **URL/source** | `https://huggingface.co/facebook/esm2_t6_8M_UR50D` |
| **license** | **MIT** |
| **download_required** | **YES** |
| **download_size** | **≈ 31 MB** ← **< 100 MB 上限 ✅** |
| **下载目录** | `F:\dsh\downloads`（HF_HOME 指向 F 盘） |
| **缓存目录** | `F:\dsh\cache` |
| **是否微调** | **NO**（仅推理，§3.17 明确要求） |
| **是否用 GPU** | NO（本机无 CUDA）→ CPU 推理 |

## 3. 内存与规模预估（§14 要求，先估算再推理）

| 项 | 估算 |
|---|---|
| 唯一 WT 链序列 | ≤ 2 × 214 targets（lig/rec 两侧） |
| 唯一 mutant 链序列 | ≤ 4473（每 variant 一条 mutant 链） |
| 单序列长度 | median ~180–224 aa，max 3128 aa（需截断） |
| 单序列嵌入（L × 320, float16） | 224×320×2 ≈ **143 KB** |
| 4473 条序列合计（float16） | ≈ **640 MB** |
| 峰值 RAM（单条前向，batch=1，长度 ≤ 1024 截断） | < **1.5 GB** |

⇒ 采用 **逐序列/小批量流式** 推理 + **按 sequence hash 去重缓存**，
**绝不一次性全部载入内存**（§14）。

## 4. 若下载失败的处置

若 `facebook/esm2_t8_8M_UR50D` 不可获取（网络/许可问题）：

```text
V7_2_BLOCKED_NO_RESIDUE_MODEL
```

并**生成 partial final report**（§PHASE 0 明确要求）。
**不会**用 V7.0 的静态局部几何冒充 residue-level（§17 末句明令禁止）。
