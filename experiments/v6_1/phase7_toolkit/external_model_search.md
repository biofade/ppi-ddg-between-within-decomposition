# V6.1 PHASE 7.2 — 外部模型预测文件搜索

- 时间: 2026-10-13 23:29:24
- 原则: **不跑模型推理**（V6.1 约束 3）；只找已跑好的预测文件；只下 <2 MB 的小文件
- 检索: GitHub Repository Search + 各仓库 git tree；Zenodo / Figshare 见文末说明

## SAMAffinity


- license: `None` | size 9299 KB | 文件 121 | 数据文件 50 | pushed 2025-12-15T03:12:17Z
- **无候选预测文件**
- 但有较大的数据文件（可能是真值/数据集，无预测列）：
  - `data/S4169.csv` (4278188 B)
  - `data/M1101.csv` (1561914 B)
  - `data/S645.csv` (1006527 B)
  - `data/S1131.csv` (763063 B)
  - `data/test_fold/S4169/S4169_test_fold3.csv` (446643 B)
  - `data/test_fold/S4169/S4169_test_fold8.csv` (442144 B)

## Pythia-PPI


- license: `Apache-2.0` | size 17071 KB | 文件 32 | 数据文件 0 | pushed 2026-06-13T03:14:31Z
- **无候选预测文件**

## MINT


- license: `MIT` | size 3615 KB | 文件 68 | 数据文件 2 | pushed 2026-06-15T16:56:53Z
- **无候选预测文件**

## DeepPPAPred

- GitHub 搜索 `DeepPPAPred` → 0 个候选：
- **未找到 GitHub 仓库**

## MuToN

- GitHub 搜索 `MuToN` → 5 个候选：
  - `zpliulab/MuToN` (★6, 577 KB, license=MIT, pushed 2025-12-02) — 
  - `mutonaothenton/mutonaothenton.github.io` (★1, 45198 KB, license=None, pushed 2022-07-20) — 
  - `Mindera/muton` (★15, 401 KB, license=MIT, pushed 2018-08-01) — A feature toggle tool with support for feature throttling and multivariance testing.
  - `trailofbits/muton` (★9, 530 KB, license=AGPL-3.0, pushed 2026-06-30) — Mutation testing framework for TON smart contracts
  - `raydaly/muton` (★5, 116 KB, license=None, pushed 2011-10-18) — a couchdb administration console for mobile using jquerymobile
- 取最相关的一个做深入检查：`zpliulab/MuToN`

- license: `MIT` | size 577 KB | 文件 38 | 数据文件 2 | pushed 2025-12-02T03:16:05Z
- **无候选预测文件**
- 但有较大的数据文件（可能是真值/数据集，无预测列）：
  - `lists/skempi_v2.csv` (1602208 B)

## ProBASS

- GitHub 搜索 `ProBASS` → 5 个候选：
  - `sagagugit/ProBASS` (★13, 2394 KB, license=None, pushed 2026-05-24) — 
  - `Gryshchenko/probass` (★0, 11485 KB, license=None, pushed 2018-08-16) — 
  - `chancasTIC/probass` (★0, 0 KB, license=None, pushed 2022-02-08) — 
  - `SayuruUSJ/ProbassLK-Admin-frontend` (★0, 94 KB, license=None, pushed 2026-07-25) — 
  - `SayuruUSJ/ProbassLK-frontend-` (★0, 60232 KB, license=None, pushed 2026-07-26) — 
- 取最相关的一个做深入检查：`sagagugit/ProBASS`

- license: `None` | size 2394 KB | 文件 10 | 数据文件 2 | pushed 2026-05-24T07:44:26Z
- **无候选预测文件**
- 但有较大的数据文件（可能是真值/数据集，无预测列）：
  - `Data_set.xlsx` (520566 B)

## CATH-ddG


- license: `None` | size 211284 KB | 文件 70 | 数据文件 6 | pushed 2025-12-15T14:46:32Z

| 候选预测文件 | size (B) |
|---|---:|
| `data/SKEMPI2/skempi_v2_with_all_results.csv` | 1291260 |

## DSSA-PPI


- license: `MIT` | size 2466 KB | 文件 99 | 数据文件 33 | pushed 2026-01-24T08:04:59Z

| 候选预测文件 | size (B) |
|---|---:|
| `validate_on_RBM/pred_vs_exp.csv` | 43492 |

## Zenodo / Figshare 检索

> V6.1 要求同时查 Zenodo / Figshare。以下为检索尝试结果（如实记录）。

- `Zenodo` 查询 `PPI ddG prediction` → 3 条结果
  - spxuw/PPI-prediction: v1.0
  - Cultural heritage – the first research campus in Hungary
  - Train/test splits for experiments with DDG prediction
- `Zenodo` 查询 `SKEMPI prediction results` → 3 条结果
  - MuToN dataset
  - Learning Functional Properties of Proteins with Language Models
  - PPIformer
- `Figshare` 查询 `PPI ddG prediction` → **检索失败/被阻断** （<html>
<head><title>403 Forbidden</title></head>
<body>
<center><h1>403 Forbidde）
- `Figshare` 查询 `SKEMPI prediction results` → **检索失败/被阻断** （<html>
<head><title>403 Forbidden</title></head>
<body>
<center><h1>403 Forbidde）

## 汇总

| 模型 | 状态 | repo | license | 说明 |
|---|---|---|---|---|
| SAMAffinity | `DATA_ONLY_NO_PRED` | `iobio-zjut/SAMAffinity` | None | V6.0 已下载 58 个 test_fold CSV（仅含 ddG 真值） |
| Pythia-PPI | `CODE_ONLY_NO_PRED` | `Wublab/pythia_ppi` | Apache-2.0 | V6.0 未找到预测文件 |
| MINT | `CODE_ONLY_NO_PRED` | `VarunUllanat/mint` | MIT | V6.0 未找到预测文件 |
| DeepPPAPred | `REPO_NOT_FOUND` | `—` | — | 需搜索 |
| MuToN | `DATA_ONLY_NO_PRED` | `zpliulab/MuToN` | MIT | 需搜索 |
| ProBASS | `DATA_ONLY_NO_PRED` | `sagagugit/ProBASS` | None | 需搜索 |
| CATH-ddG | `PREDICTION_FOUND` | `ak422/CATH-ddG` | None | V6.0 已获取基线预测列 |
| DSSA-PPI | `PREDICTION_FOUND` | `Agitw/DSSA-PPI` | MIT | V6.0 已获取 8 模型预测列 |

