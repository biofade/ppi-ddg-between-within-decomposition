# 外部证据（EXACT_UNIQUE only，§13.3/§13.5）

**宇宙**：`U_EXT_EXACT_UNIQUE` — 1781 行 / 168 targets / 67 个 n≥5 / 93 families

**排除**：`AMBIGUOUS_KEY` 2048 行、`UNMATCHED` 644 行

> 此前基于完整键匹配（约 85%）的读数**主要由歧义键构成**，已作废。

## 方法级能力画像（首选呈现）

| model | n | n≥5 | Overall r | Macro ρ | C_B | C_W | Eff@10 | RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Flex-ddG` | 1755 | 65 | +0.2028 | +0.1092 | +0.0478 | +0.1550 | +nan | 2.1636 |
| `DDGPred` | 1756 | 65 | +0.2406 | +0.1287 | +0.0715 | +0.1690 | +nan | 2.1262 |
| `ddG_FoldX` | 1776 | 66 | +0.2089 | +0.1088 | +0.0657 | +0.1433 | +nan | 2.3691 |
| `ESM-IF` | 1776 | 66 | -0.2250 | +0.0213 | -0.1149 | -0.1101 | +nan | 2.3587 |
| `proteinMPNN` | 1776 | 66 | +0.1509 | +0.0320 | +0.0699 | +0.0810 | +nan | 2.3199 |

## 外部秩一致（**低功效**）

```
{
  "n_methods": 5,
  "rho_overall_pearson_vs_macro_spearman": 0.8999999999999998,
  "tau_overall_pearson_vs_macro_spearman": 0.7999999999999999,
  "perm_p_overall_pearson_vs_macro_spearman": 0.08333333333333333,
  "rho_overall_pearson_vs_efficiency@10": 0.6,
  "tau_overall_pearson_vs_efficiency@10": 0.39999999999999997,
  "perm_p_overall_pearson_vs_efficiency@10": 0.35,
  "rho_macro_spearman_vs_efficiency@10": 0.7,
  "tau_macro_spearman_vs_efficiency@10": 0.6,
  "perm_p_macro_spearman_vs_efficiency@10": 0.23333333333333334,
  "low_power_warning": "n_methods=5 < 8 ⇒ rank correlation is UNDERPOWERED; report method-level profiles rather than a correlation headline",
  "overall_winner": "DDGPred",
  "ranking_winner": "DDGPred",
  "same_winner": true,
  "comparison_to_ambiguous_mapping": [
    {
      "model": "Flex-ddG",
      "overall_pearson_cleanEU": 0.20283191346540427,
      "macro_spearman_cleanEU": 0.10915022730387325,
      "overall_pearson_inclAmbiguous": 0.2618655601054078,
      "macro_spearman_inclAmbiguous": 0.1319757779785927,
      "delta_overall": -0.05903364664000352,
      "delta_macro": -0.022825550674719447
    },
    {
      "model": "DDGPred",
      "overall_pearson_cleanEU": 0.24057145791062529,
      "macro_spearman_cleanEU": 0.12868933419751688,
      "overall_pearson_inclAmbiguous": 0.3357429476893619,
      "macro_spearman_inclAmbiguous": 0.10472869177707,
      "delta_overall": -0.09517148977873663,
      "delta_macro": 0.023960642420446876
    },
    {
      "model": "ddG_FoldX",
      "overall_pearson_cleanEU": 0.2089130565869603,
      "macro_spearman_cleanEU": 0.10878355125874062,
      "overall_pearson_inclAmbiguous": 0.2727472435637179,
      "macro_spearman_inclAmbiguous": 0.1502382989702551,
      "delta_overall": -0.06383418697675758,
      "delta_macro": -0.04145474771151447
    },
    {
      "model": "ESM-IF",
      "overall_pearson_cleanEU": -0.22498534652089358,
      "macro_spearman_cleanEU": 0.021291001382564356,
      "overall_pearson_inclAmbiguous": -0.1013982345349476,
      "macro_spearman_inclAmbiguous": -0.0346327471193198,
      "delta_overall": -0.12358711198594598,
      "delta_macro": 0.055923748501884155
    },
    {
      "model": "proteinMPNN",
      "overall_pearson_cleanEU": 0.15090845873032877,
      "macro_spearman_cleanEU": 0.031994320163390245,
      "overall_pearson_inclAmbiguous": 0.1673771161714917,
      "macro_spearman_inclAmbiguous": 0.0919437798125134,
      "delta_overall": -0.016468657441162937,
      "delta_macro": -0.05994945964912315
    }
  ]
}
```

> n 过小 ⇒ 秩相关**不宜作为头条**；以上方法级画像为准（§13.5）。

