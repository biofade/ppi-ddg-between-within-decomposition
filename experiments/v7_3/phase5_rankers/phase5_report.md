# V7.3 PHASE 5 — 冻结嵌套 family-CV（F0–F9）

- 时间: 2026-10-15 01:16:12
- learner: **LightGBM**（trees=300）；Ridge = 可解释诊断（不选优）
- 预处理：**train-only** 中位数插补 + 标准化；**无 PCA**（§31）

- **U_PHYS_STRICT: rows 3351 / targets 166 / families 86**
- n>=5 target：**102**

## 5.1 冻结特征族（§26）

| 族 | n_features | 构成 |
|---|---:|---|
| `F0_Base` | 188 | 含 Base |
| `F1_Base_Mut` | 236 | 含 Base |
| `F2_Base_Static` | 307 | 含 Base |
| `F3_Base_EnvPartner` | 275 | 含 Base |
| `F4_Base_PartnerPert` | 263 | 含 Base |
| `F5_Base_IntraPert` | 254 | 含 Base |
| `F6_Partner_Intra` | 281 | 含 Base |
| `F7_F6_Burial` | 298 | 含 Base |
| `F8_Full` | 303 | 含 Base |
| `F9_PhysicsOnly` | 115 | **无 Base** |

> 出现在多个族中的列已去重（见 `family_spec.json`）

## 5.2 嵌套 family-CV

| F0_Base | s42 | Macro ρ=+0.0434 | Eff@10=+0.1933 | C_W=+0.2663 | Kendall=+0.0295 | r=+0.3599 |
| F1_Base_Mut | s42 | Macro ρ=+0.0734 | Eff@10=+0.1189 | C_W=+0.2976 | Kendall=+0.0685 | r=+0.3608 |
| F2_Base_Static | s42 | Macro ρ=+0.0629 | Eff@10=+0.1116 | C_W=+0.2986 | Kendall=+0.0945 | r=+0.3887 |
| F3_Base_EnvPartner | s42 | Macro ρ=+0.0873 | Eff@10=+0.1154 | C_W=+0.2932 | Kendall=+0.0428 | r=+0.4042 |
| F4_Base_PartnerPert | s42 | Macro ρ=+0.1054 | Eff@10=+0.0679 | C_W=+0.3117 | Kendall=+0.1364 | r=+0.4310 |
| F5_Base_IntraPert | s42 | Macro ρ=+0.0459 | Eff@10=+0.1203 | C_W=+0.2942 | Kendall=+0.0476 | r=+0.3558 |
| F6_Partner_Intra | s42 | Macro ρ=+0.1258 | Eff@10=+0.0938 | C_W=+0.3171 | Kendall=+0.1068 | r=+0.4164 |
| F7_F6_Burial | s42 | Macro ρ=+0.1047 | Eff@10=+0.1451 | C_W=+0.3179 | Kendall=+0.0985 | r=+0.4582 |
| F8_Full | s42 | Macro ρ=+0.1394 | Eff@10=+0.1762 | C_W=+0.3245 | Kendall=+0.1408 | r=+0.4662 |
| F9_PhysicsOnly | s42 | Macro ρ=+0.0756 | Eff@10=+0.1230 | C_W=+0.1946 | Kendall=+0.0520 | r=+0.3610 |
| F0_Base | s43 | Macro ρ=+0.0153 | Eff@10=+0.0985 | C_W=+0.2468 | Kendall=-0.0160 | r=+0.3335 |
| F1_Base_Mut | s43 | Macro ρ=+0.0483 | Eff@10=+0.1164 | C_W=+0.2566 | Kendall=+0.0728 | r=+0.3291 |
| F2_Base_Static | s43 | Macro ρ=+0.0770 | Eff@10=+0.1302 | C_W=+0.2716 | Kendall=+0.0374 | r=+0.3741 |
| F3_Base_EnvPartner | s43 | Macro ρ=+0.1191 | Eff@10=+0.1746 | C_W=+0.2885 | Kendall=+0.1074 | r=+0.4311 |
| F4_Base_PartnerPert | s43 | Macro ρ=+0.1198 | Eff@10=+0.1302 | C_W=+0.2847 | Kendall=+0.0785 | r=+0.4275 |
| F5_Base_IntraPert | s43 | Macro ρ=+0.0383 | Eff@10=+0.1614 | C_W=+0.2543 | Kendall=+0.0583 | r=+0.3414 |
| F6_Partner_Intra | s43 | Macro ρ=+0.0810 | Eff@10=+0.1139 | C_W=+0.2960 | Kendall=+0.0698 | r=+0.4273 |
| F7_F6_Burial | s43 | Macro ρ=+0.0989 | Eff@10=+0.1890 | C_W=+0.3047 | Kendall=+0.0845 | r=+0.4461 |
| F8_Full | s43 | Macro ρ=+0.0741 | Eff@10=+0.0535 | C_W=+0.2910 | Kendall=+0.0572 | r=+0.4588 |
| F9_PhysicsOnly | s43 | Macro ρ=+0.0839 | Eff@10=+0.1349 | C_W=+0.2026 | Kendall=+0.0513 | r=+0.3963 |
| F0_Base | s44 | Macro ρ=+0.0364 | Eff@10=+0.1495 | C_W=+0.2501 | Kendall=+0.0419 | r=+0.3621 |
| F1_Base_Mut | s44 | Macro ρ=+0.0340 | Eff@10=+0.0819 | C_W=+0.2766 | Kendall=+0.0307 | r=+0.3434 |
| F2_Base_Static | s44 | Macro ρ=+0.0439 | Eff@10=+0.1326 | C_W=+0.2705 | Kendall=+0.0471 | r=+0.3641 |
| F3_Base_EnvPartner | s44 | Macro ρ=+0.0645 | Eff@10=+0.1393 | C_W=+0.2747 | Kendall=+0.0976 | r=+0.4025 |
| F4_Base_PartnerPert | s44 | Macro ρ=+0.1141 | Eff@10=+0.1947 | C_W=+0.3089 | Kendall=+0.0922 | r=+0.4332 |
| F5_Base_IntraPert | s44 | Macro ρ=+0.0495 | Eff@10=+0.1010 | C_W=+0.2767 | Kendall=+0.0763 | r=+0.3345 |
| F6_Partner_Intra | s44 | Macro ρ=+0.1095 | Eff@10=+0.1991 | C_W=+0.3123 | Kendall=+0.0845 | r=+0.4433 |
| F7_F6_Burial | s44 | Macro ρ=+0.1230 | Eff@10=+0.1572 | C_W=+0.3114 | Kendall=+0.1387 | r=+0.4429 |
| F8_Full | s44 | Macro ρ=+0.1007 | Eff@10=+0.1489 | C_W=+0.3121 | Kendall=+0.0899 | r=+0.4571 |
| F9_PhysicsOnly | s44 | Macro ρ=+0.0567 | Eff@10=+0.1461 | C_W=+0.2173 | Kendall=+0.0371 | r=+0.4158 |

- **泄漏断言**：outer violation = **0**
- **train-only 插补/标准化：PASS**；**无 PCA**

### 3-seed 均值（敏感性口径）

| 族 | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Overall r |
|---|---:|---:|---:|---:|---:|---:|
| `F0_Base` | +0.0317 | +0.1471 | +0.2544 | +0.0185 | +0.6962 | +0.3519 |
| `F1_Base_Mut` | +0.0519 | +0.1058 | +0.2769 | +0.0573 | +0.7117 | +0.3444 |
| `F2_Base_Static` | +0.0613 | +0.1248 | +0.2802 | +0.0597 | +0.7154 | +0.3756 |
| `F3_Base_EnvPartner` | +0.0903 | +0.1431 | +0.2855 | +0.0826 | +0.7169 | +0.4126 |
| `F4_Base_PartnerPert` | +0.1131 | +0.1309 | +0.3018 | +0.1024 | +0.7187 | +0.4306 |
| `F5_Base_IntraPert` | +0.0446 | +0.1276 | +0.2751 | +0.0608 | +0.7104 | +0.3439 |
| `F6_Partner_Intra` | +0.1054 | +0.1356 | +0.3085 | +0.0870 | +0.7202 | +0.4290 |
| `F7_F6_Burial` | +0.1089 | +0.1638 | +0.3113 | +0.1073 | +0.7211 | +0.4491 |
| `F8_Full` | +0.1048 | +0.1262 | +0.3092 | +0.0960 | +0.7262 | +0.4607 |
| `F9_PhysicsOnly` | +0.0720 | +0.1347 | +0.2049 | +0.0468 | +0.6320 | +0.3910 |

### ensemble prediction 口径（主口径）

| 族 | Macro ρ | Eff@10 | C_W |
|---|---:|---:|---:|
| `F0_Base` | +0.0240 | +0.1261 | +0.2638 |
| `F1_Base_Mut` | +0.0528 | +0.1378 | +0.2923 |
| `F2_Base_Static` | +0.0668 | +0.1369 | +0.2939 |
| `F3_Base_EnvPartner` | +0.0937 | +0.1470 | +0.2985 |
| `F4_Base_PartnerPert` | +0.1259 | +0.1573 | +0.3168 |
| `F5_Base_IntraPert` | +0.0531 | +0.1276 | +0.2898 |
| `F6_Partner_Intra` | +0.1202 | +0.1546 | +0.3235 |
| `F7_F6_Burial` | +0.1220 | +0.0975 | +0.3246 |
| `F8_Full` | +0.1147 | +0.0891 | +0.3226 |
| `F9_PhysicsOnly` | +0.0859 | +0.1593 | +0.2123 |

- 耗时 19.0 分钟
