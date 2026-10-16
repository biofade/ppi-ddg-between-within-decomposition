# V7.2 PHASE 1 — 残基级嵌入（ESM-2 8M，流式 + 续跑）

- 时间: 2026-10-14 18:24:14
- model: **facebook/esm2_t6_8M_UR50D**（7.84M 参数，hidden 320，MIT）
- 缓存: `F:\dsh\cache\huggingface`（**未写 C 盘**）
- 仅推理，**不微调**；CPU

- 模型加载完成；hidden = 320

- U_CONTRACT rows **4473** / targets **214**
- 权威标签来源: `canonical_target_labels.csv`（§3.3 强制）

  - 已完成 10 个 target（跳过 1），累计 rows 1803，耗时 36.1 min
  - 已完成 20 个 target（跳过 1），累计 rows 2429，耗时 43.1 min
  - 已完成 30 个 target（跳过 2），累计 rows 2850，耗时 48.0 min
  - 已完成 40 个 target（跳过 3），累计 rows 3110，耗时 50.3 min
  - 已完成 50 个 target（跳过 3），累计 rows 3310，耗时 53.2 min
  - 已完成 60 个 target（跳过 3），累计 rows 3463，耗时 58.0 min
  - 已完成 70 个 target（跳过 4），累计 rows 3573，耗时 60.4 min
  - 已完成 80 个 target（跳过 4），累计 rows 3660，耗时 63.4 min
  - 已完成 90 个 target（跳过 5），累计 rows 3735，耗时 64.4 min
  - 已完成 100 个 target（跳过 5），累计 rows 3790，耗时 65.5 min
  - 已完成 110 个 target（跳过 5），累计 rows 3839，耗时 66.1 min
  - 已完成 120 个 target（跳过 5），累计 rows 3880，耗时 67.4 min
  - 已完成 130 个 target（跳过 5），累计 rows 3917，耗时 67.9 min
  - 已完成 140 个 target（跳过 5），累计 rows 3944，耗时 68.5 min
  - 已完成 150 个 target（跳过 9），累计 rows 3964，耗时 68.9 min
  - 已完成 160 个 target（跳过 16），累计 rows 3977，耗时 69.0 min

## 结果
- 本轮新完成 target: **168**；已存在跳过: **46**
- 本轮 rows: **3988**
- **累计已缓存 target: 174；累计 rows: 4233**
- 覆盖率（对 4473 U_CONTRACT 行）: **94.6%**

- 耗时 69.1 分钟
