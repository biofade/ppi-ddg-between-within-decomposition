# P0 — 外部 benchmark 来源溯源（V5 §120）

生成时间: 2026-10-13 20:17:23

## 一手来源

| 项 | 值 |
|---|---|
| 论文 | CATH-ddG: towards robust mutation effect prediction on protein–protein interactions out of CATH homologous superfamily |
| 作者 | Guanglei Yu, Xuehua Bi, Teng Ma, Yaohang Li, Jianxin Wang |
| 期刊/卷期 | Bioinformatics, 2025, 41(Suppl 1): i362–i372 |
| DOI | 10.1093/bioinformatics/btaf228 |
| PMID / PMCID | 40662779 / PMC12261453 |
| 许可（论文） | PMC Copyright Notice（开放获取） |
| 官方代码 | https://github.com/ak422/CATH-ddG |
| **代码许可** | **`license` 字段为 `None` —— 仓库未声明任何许可** |
| 代码大小 | 211,284 KB |
| 下载日期 | 2026-10-13 |

## 已下载并落盘的外部文件

| 文件 | 大小 (B) | SHA256 | 说明 |
|---|---:|---|---|
| `F:\dsh\downloads\cathddg_oa.bin` | 246,881 | - | Europe PMC 全文 XML（用于提取表 2–5 数值） |
| `F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv` | 1,291,260 | - | 含官方 `cath_fold` split 标签 + 各 baseline 预测 |
| `F:\dsh\downloads\cathddg\CATH.yml` | 2,373 | - | 官方训练/测试配置（3 个 checkpoint） |
| `F:\dsh\downloads\cathddg\zero_domain_foldseek.csv` | 1,960 | - | 未注释链的 Foldseek 域注释补全 |
| `F:\dsh\downloads\cathddg\README.md` | 3,160 | - | 官方运行说明 |

> 注：`configs/` 与部分 `data/` 文件通过 GitHub contents API（base64）获取；
> `raw.githubusercontent.com` 直连返回 0 字节（已如实记录）。
> 按 V5 §121：未执行任何来源不明的 shell script；未写入任何 API key。

## 官方 split 的实际取值（从 `cath_fold` 列读出，PROJECT_EVIDENCE）

| cath_fold | mutations | 唯一 complex | 唯一 PPI |
|---|---:|---:|---:|
| `train` | 5,893 | 287 | 290 |
| `val` | **813** | **53** | **53** |
| `NaN`（被过滤的不可用条目） | 379 | — | — |
| 合计 | 7,085 | 345 | — |

**关键口径澄清**：论文正文所说的「held-out CATH blind test set ... 813 mutations」
对应本表的 **`val` 分区**（不是 `test` 分区，该列无 `test` 取值）。
论文 Table 2 的数字即在该 813 条上评估。

## sign convention 核实（V5 §17）

* CATH-ddG 表内 `ddG` 列取值范围查证：mean 约 **+1.5 kcal/mol**，
  正值为**削弱结合**（destabilizing），与 SKEMPI v2.0 的 `ddG` 定义一致。
* 我方定义 `ΔΔG = ΔG_mut − ΔG_WT`，其中 `ΔG = 0.592·ln(KD[M])`
  → ΔΔG > 0 同样是削弱结合。**两侧 sign convention 一致**，无需反号。
* 我方 ΔΔG 分布（PROJECT_EVIDENCE, `data/v4.0/ddg_pairs.csv`）：
  mean **+1.5215**、ΔΔG>0 占 **74.5%** —— 与 CATH `ddG` 均值同号同量级。
