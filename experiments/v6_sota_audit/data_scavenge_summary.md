# V6.0 PHASE 1 — SOTA 预测文件搜刮

- 时间: 2026-10-13 22:42:23
- 原则: **不跑任何模型推理**，只找已跑好的预测文件；只下载 < 2 MB 的文件

| 仓库 | 存在 | license | size KB | 文件数 | 候选预测文件 |
|---|---|---|---:|---:|---:|
| SAMAffinity | ✅ | None | 9299 | 121 | **73** |
| MINT | ✅ | MIT | 3615 | 68 | **2** |
| DSSA-PPI | ✅ | MIT | 2466 | 99 | **3** |
| Pythia-PPI | ✅ | Apache-2.0 | 17071 | 32 | **0** |

## 各仓库的候选预测文件（按大小升序，只列前 15）

### SAMAffinity — `iobio-zjut/SAMAffinity`

- license: `None` | size 9299 KB | 文件 121 | pushed 2025-12-15T03:12:17Z

| path | size (B) |
|---|---:|
| `main/M1101/predict/pytorchtools.py` | 1854 |
| `main/S1131/predict/pytorchtools.py` | 1854 |
| `main/S4169/predict/pytorchtools.py` | 1854 |
| `main/S645/predict/pytorchtools.py` | 1854 |
| `main/M1101/predict/lookahead.py` | 2476 |
| `main/S1131/predict/lookahead.py` | 2476 |
| `main/S4169/predict/lookahead.py` | 2476 |
| `main/S645/predict/lookahead.py` | 2476 |
| `main/S645/predict/__pycache__/lookahead.cpython-37.pyc` | 2554 |
| `main/M1101/predict/__pycache__/lookahead.cpython-37.pyc` | 2555 |
| `main/S1131/predict/__pycache__/lookahead.cpython-37.pyc` | 2555 |
| `main/S4169/predict/__pycache__/lookahead.cpython-37.pyc` | 2555 |
| `main/S645/predict/__pycache__/Radam.cpython-37.pyc` | 5552 |
| `main/M1101/predict/__pycache__/Radam.cpython-37.pyc` | 5553 |
| `main/S1131/predict/__pycache__/Radam.cpython-37.pyc` | 5553 |

- 尝试下载 73 个小文件到 `F:\dsh\downloads\samaffinity`:
  - `main/M1101/predict/pytorchtools.py` → 0 B (empty_or_blocked)
  - `main/S1131/predict/pytorchtools.py` → 1854 B (downloaded)
  - `main/S4169/predict/pytorchtools.py` → 1854 B (cached)
  - `main/S645/predict/pytorchtools.py` → 1854 B (cached)
  - `main/M1101/predict/lookahead.py` → 2476 B (downloaded)
  - `main/S1131/predict/lookahead.py` → 2476 B (cached)
  - `main/S4169/predict/lookahead.py` → 2476 B (cached)
  - `main/S645/predict/lookahead.py` → 2476 B (cached)
  - `main/S645/predict/__pycache__/lookahead.cpython-37.pyc` → 2554 B (downloaded)
  - `main/M1101/predict/__pycache__/lookahead.cpython-37.pyc` → 2554 B (cached)
  - `main/S1131/predict/__pycache__/lookahead.cpython-37.pyc` → 2554 B (cached)
  - `main/S4169/predict/__pycache__/lookahead.cpython-37.pyc` → 2554 B (cached)
  - `main/S645/predict/__pycache__/Radam.cpython-37.pyc` → 5552 B (downloaded)
  - `main/M1101/predict/__pycache__/Radam.cpython-37.pyc` → 5552 B (cached)
  - `main/S1131/predict/__pycache__/Radam.cpython-37.pyc` → 5552 B (cached)
  - `main/S4169/predict/__pycache__/Radam.cpython-37.pyc` → 5552 B (cached)
  - `main/M1101/predict/Radam.py` → 8262 B (downloaded)
  - `main/S1131/predict/Radam.py` → 8262 B (cached)
  - `main/S4169/predict/Radam.py` → 8262 B (cached)
  - `main/S645/predict/Radam.py` → 8262 B (cached)
  - `main/S1131/predict/predict_S1131.py` → 12942 B (downloaded)
  - `main/S645/predict/__pycache__/model_S645.cpython-37.pyc` → 13290 B (downloaded)
  - `main/M1101/predict/__pycache__/model_M1101.cpython-37.pyc` → 13292 B (downloaded)
  - `main/M1101/predict/__pycache__/model_AB1101_v49.cpython-37.pyc` → 13297 B (downloaded)
  - `main/S4169/predict/predict_S4169.py` → 13782 B (downloaded)
  - `main/S1131/predict/__pycache__/model_S1131_v56.cpython-37.pyc` → 15069 B (downloaded)
  - `main/S645/predict/predict_S645.py` → 15433 B (downloaded)
  - `main/M1101/predict/predict_M1101.py` → 15514 B (downloaded)
  - `main/S4169/predict/__pycache__/model_S4169.cpython-37.pyc` → 15832 B (downloaded)
  - `main/M1101/predict/model_M1101.py` → 17102 B (downloaded)
  - `main/S645/predict/model_S645.py` → 17102 B (downloaded)
  - `main/S4169/predict/model_S4169.py` → 19508 B (downloaded)
  - `main/S1131/predict/model_S1131.py` → 20804 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold8.csv` → 72862 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold2.csv` → 74695 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold1.csv` → 74790 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold10.csv` → 75886 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold9.csv` → 75941 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold4.csv` → 76512 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold5.csv` → 76521 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold6.csv` → 78436 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold3.csv` → 78612 B (downloaded)
  - `data/test_fold/S1131/S1131_test_fold7.csv` → 79394 B (downloaded)
  - `data/test_fold/S645/S645_test_fold7.csv` → 85254 B (downloaded)
  - `data/test_fold/S645/S645_test_fold10.csv` → 92538 B (downloaded)
  - `data/test_fold/S645/S645_test_fold1.csv` → 97784 B (downloaded)
  - `data/test_fold/S645/S645_test_fold6.csv` → 98084 B (downloaded)
  - `data/test_fold/S645/S645_test_fold9.csv` → 98603 B (downloaded)
  - `data/test_fold/S645/S645_test_fold5.csv` → 101571 B (downloaded)
  - `data/test_fold/S645/S645_test_fold8.csv` → 102811 B (downloaded)
  - `data/test_fold/S645/S645_test_fold2.csv` → 103468 B (downloaded)
  - `data/test_fold/S645/S645_test_fold3.csv` → 113024 B (downloaded)
  - `data/test_fold/S645/S645_test_fold4.csv` → 115883 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold2.csv` → 148316 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold7.csv` → 148740 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold9.csv` → 150168 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold10.csv` → 152791 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold6.csv` → 153416 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold1.csv` → 157518 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold4.csv` → 159984 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold8.csv` → 160509 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold5.csv` → 164065 B (downloaded)
  - `data/test_fold/M1101/M1101_test_fold3.csv` → 172043 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold9.csv` → 403576 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold10.csv` → 411266 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold7.csv` → 411911 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold6.csv` → 429316 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold5.csv` → 432004 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold2.csv` → 432020 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold1.csv` → 433829 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold4.csv` → 436774 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold8.csv` → 442144 B (downloaded)
  - `data/test_fold/S4169/S4169_test_fold3.csv` → 446643 B (downloaded)

### MINT — `VarunUllanat/mint`

- license: `MIT` | size 3615 KB | 文件 68 | pushed 2026-06-15T16:56:53Z

| path | size (B) |
|---|---:|
| `mint/helpers/predict.py` | 428 |
| `downstream/CovidVariants/calculate_metrics.ipynb` | 3353 |

- 尝试下载 2 个小文件到 `F:\dsh\downloads\mint`:
  - `mint/helpers/predict.py` → 428 B (downloaded)
  - `downstream/CovidVariants/calculate_metrics.ipynb` → 3353 B (downloaded)

### DSSA-PPI — `Agitw/DSSA-PPI`

- license: `MIT` | size 2466 KB | 文件 99 | pushed 2026-01-24T08:04:59Z

| path | size (B) |
|---|---:|
| `validate_on_RBM/pred_vs_exp.csv` | 43492 |
| `validate_on_RBM/pred_heatmap_6e0j_RBM.png` | 106582 |
| `validate_on_RBM/RBM_baselines_comparison_ddG.csv` | 174025 |

- 尝试下载 3 个小文件到 `F:\dsh\downloads\dssa_ppi`:
  - `validate_on_RBM/pred_vs_exp.csv` → 43492 B (downloaded)
  - `validate_on_RBM/pred_heatmap_6e0j_RBM.png` → 106582 B (downloaded)
  - `validate_on_RBM/RBM_baselines_comparison_ddG.csv` → 174025 B (downloaded)

### Pythia-PPI — `Wublab/pythia_ppi`

- license: `Apache-2.0` | size 17071 KB | 文件 32 | pushed 2026-06-13T03:14:31Z
- **无候选预测文件** → `PREDICTION_FILE_NOT_FOUND`


