# V7.2R-2 — 负对照（NC1–NC4，§46）

- 时间: 2026-10-14 21:26:32
- U_AB: rows **4219** / targets **174** / families **91**

| 对照 | 构成 | Macro ρ | Eff@10 | C_W | Overall r |
|---|---|---:|---:|---:|---:|
| `NC1_B2_dr` | base+B_dr | +0.0305 | +0.3314 | +0.1215 | +0.1723 |
| `NC1_B1_wt_only` | base+B_wt | +0.0462 | +0.2992 | +0.1262 | +0.1702 |
| `NC2_mut_identity` | base+F0 | +0.0061 | +0.4351 | +0.1163 | +0.1709 |
| `NC3_dr_norm_only` | base+B_drnorm | +0.0240 | +0.3017 | +0.1217 | +0.1768 |
| `NC4_random_partner` | base+B_dr+B_randpartner | +0.0173 | +0.3264 | +0.1236 | +0.1941 |

## 判定要点（§46）

- **NC1**：若 `B1_wt_only`（仅 `r_wt`）已与 `B2_dr` 相当 → **扰动本身不必要**，增益来自结构身份
- **NC2**：若 `mut_identity` 与 B2 相当 → 增益只是**氨基酸替换编码**
- **NC3**：若 `dr_norm_only` 已接近 B2 → 增益只是**幅度**信息，方向/内容无贡献
- **NC4**：若 `random_partner` 与 B3 相当 → **不得**声称 partner-specific 信号

- 耗时 11.1 分钟
