# V6.2 — Scientific Repair + Decoupling Verification

**状态**：8 个 PHASE **全部完成**（PHASE 0–8），约 4 小时（上限 6 小时，未超时）
**上游**：V6.1+（发现物理方法 within-target 优于学习表示 p=0.0038；Split Inflation +0.2663）

## 本任务修复的 6 个关键风险

| # | V6.1 风险 | V6.2 处置 | 结论 |
|---|---|---|---|
| 1 | **Phase 2 任务错位** | PHASE 1 provenance 审计 | **确认错位**：V6.1 实际在**绝对 ΔG**（8207 行 `real_dG`）上做，非 7523 ΔΔG 配对 |
| 2 | **Bootstrap `unique≥15` 阈值** | PHASE 5 无阈值重算 | **确认 bug**：CI 系统性偏窄 **20–26%**；结论方向不变 |
| 3 | **「Global Pearson 来自 between」未验证** | PHASE 3 协方差分解 + 负对照 | **证实**：within-shuffle 后全局 Pearson 从 **+0.699 → −0.019** |
| 4 | **DGaux OOF 可能泄漏** | PHASE 2.4 nested-OOF + PHASE 7 断言 | **无泄漏**（50/50 通过） |
| 5 | **Forward/reverse 顺序** | PHASE 2.2 reverse 只在 outer-train 内 | **无泄漏**（test = canonical forward only） |
| 6 | **M1 重写会与历史不一致** | PHASE 2.0 审计 + **只复用不重写** | **复用** `v5_common`（仅替换 split + OOF 层） |

## PHASE 清单

| PHASE | 内容 | 状态 | 关键产物 |
|---|---|---|---|
| 0 | 环境 + 目录 | DONE | `phase0_environment/environment.md` |
| 1 | Target Provenance Audit | DONE | **判定 `TASK_IS_ABSOLUTE_DG`** |
| 2 | M1 审计 + Canonical + 真 ΔΔG 4 档 | DONE | `phase2_split4/` |
| 3 | Pearson 分解 + Primary Statistics | DONE | `phase3_decomposition/` |
| 4 | 排序指标 + δ 敏感 | DONE | `phase4_ranking_metrics/` |
| 5 | Bootstrap 修正 | DONE | `phase5_bootstrap/` |
| 6 | Leave-one-out + BH-FDR | DONE | **`ROBUST`** |
| 7 | 泄漏检查 + 浮点 sanity | DONE | **50/50 通过** |
| 8 | 汇总 + Decoupling 判定 | DONE | `V6_2_FINAL_SUMMARY.md` |

## 三条核心结论

1. **全局 Pearson 几乎全部来自「跨靶点区分」**
   - L1 下 within-shuffle 把全局 Pearson 从 **+0.6990 打到 −0.0192**
   - L4 下 between-only 相关系数 **−0.0301**、within-only **+0.1923**（与 L1 完全相反）
2. **within 信号在 OOD 后显著下降但未消失**
   - `C_W`: L1 **+0.2658** → L4 **+0.1361**（保留 **51%**）
   - `Δ_OOD = +0.1144 [+0.0642, +0.1898]`，CI 排除 0
   - **与 V6.2 预设的「丢失」不符** → 路线判定为 **`PARTIAL_MIPA`**
3. **Split Inflation（真 ΔΔG，L1−L4）= +0.5590**
   - M1: **+0.6953 → +0.1250**（对照 V6.1 绝对 ΔG 的 +0.2663，任务不同不可直接比）

## 执行命令

```powershell
python -u scripts/v6_2/01_provenance_audit.py
python -u scripts/v6_2/02a_m1_audit.py
python -u scripts/v6_2/02b_canonical_set.py
python -u scripts/v6_2/02_split4_ddg.py --levels L1_random L2_complex L3_seq40 L4_family --seeds 42 43 44 --trees 300
python -u scripts/v6_2/03_decomposition.py
python -u scripts/v6_2/05_bootstrap_fixed.py
python -u scripts/v6_2/06_leave_one_out.py
python -u scripts/v6_2/07_leakage_assertions.py
python -u scripts/v6_2/04_ranking_metrics.py
python -u scripts/v6_2/10_final_verify.py
```

## 约束遵守

- ✅ 未修改 V1.1–V6.1 任何文件（自查脚本已核验）
- ✅ `v6_2` 命名空间
- ✅ CPU only；**全程无下载**
- ✅ **未伪造任何数字**；与预设不符的结果已显式报告
- ✅ 每完成一个 PHASE 立即更新 `PROGRESS.md`
- ✅ **判定只依据 provenance / IDs / formula / sign convention / hashes**

## 三处必须注意的如实记录

1. **四档有效层级数为 3**：本数据集中每个 PDB 恰属一个 family 且
   `family == seq40_cluster` → **L3 ≡ L2**（数据性质，非实现缺陷）。
2. **M2 ≡ M0**：故 `Gain_W(M2−M0) = 0` 表示**该对照未提供信息**，
   不能解读为「D 差分无增益」。
3. **PHASE 6 的 nuance**：去掉 ESM-IF 后物理 vs 序列差异**不再显著**（p=0.111）
   → 结论主要由 ESM-IF 这一个弱模型驱动；
   论文措辞须限定为 *"among the evaluated explicit-physics and
   structure-conditioned inverse-folding models"*。
