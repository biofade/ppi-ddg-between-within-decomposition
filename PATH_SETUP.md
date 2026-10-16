# PATH_SETUP

**为什么需要这一步。** 部分分析脚本（约 48 个，集中在 `scripts/v8_1/`、`scripts/v6_*/`、
`scripts/v7_*`）硬编码了本研究的原始工作路径，另有 12 个脚本硬编码第三方数据路径
（`F:/dsh/downloads/...`）。因此把仓库 clone 到别的位置后，**部分脚本**会找不到文件。

有两种解决方式，任选其一。**校验器不需要任何设置**——它用自身位置解析全部输入。

---

## 方式 A（推荐）：建立路径链接，零代码改动

```bash
python setup_original_paths.py          # 建立链接
python setup_original_paths.py --check  # 确认就绪
```

它会在**原路径不存在**时，创建指向你 clone 位置的 junction（Windows，无需管理员权限）
或符号链接（macOS/Linux）：

```
F:\DSH_Galatea\Affinity_Intelligence   ->   <你 clone 的位置>
F:\dsh\downloads                       ->   你放置第三方数据的目录
```

已经是真实目录的路径**不会被改动**（脚本会拒绝并提示）。
用 `--remove` 可撤销。

若你的平台不允许创建链接（Windows 未开开发者模式、且非 NTFS），改用方式 B。

可以覆盖默认路径：

```bash
# Windows PowerShell
$env:PPIDDG_ORIG_ROOT = "D:\work\Affinity_Intelligence"
$env:PPIDDG_EXT_DIR   = "D:\data\ppi_ddg"
python setup_original_paths.py
```

---

## 方式 B：把仓库放回原路径

最省事，不需要任何链接：

```
F:\DSH_Galatea\Affinity_Intelligence
```

把仓库内容放到这个位置即可（Windows 上即 `F:` 盘的该目录）。

---

## 哪些内容开箱即用，哪些需要上面任一步

| 功能 | 是否需要 PATH_SETUP |
|---|---|
| `python validate_final_manuscript.py` —— 全部 67 项一致性检查 | **不需要** |
| 编译 LaTeX 重建主稿/补充 PDF | **不需要**（`latex_v6/` 自包含） |
| 重绘全部图（`figures_v5/_generator/make_all_figures.py`） | 一般**不需要**（部分输入在仓库内） |
| 重跑四档受控划分 / 分解（`scripts/v6_2/`） | **需要** |
| 重跑 V8.1 全部分析链（`scripts/v8_1/01…09`） | **需要** |
| 重跑 ranker 与嵌套 CV（`scripts/v7_1/`） | **需要** |

---

## 第三方数据

`F:\dsh\downloads` 下的**第三方数据不随仓库分发**（许可证原因）。下列文件需自行下载：

| 文件 | 用途 | 获取 |
|---|---|---|
| `cathddg/skempi_v2.csv` | 五个外部预测列的来源表（SKEMPI 派生，CATH-ddG 发布） | 见 `THIRD_PARTY_INPUTS.csv`（含 SHA-256） |
| `hermes_variants/…` | HERMES 语义审计（可选） | 见 `scripts/v8_1/04_external_audits.py` 中的官方来源 URL |

`fetch_third_party.py` 会按 `THIRD_PARTY_INPUTS.csv` 记录的 URL 下载并校验 SHA-256。
**注意**：该 CSV 里的 URL 目前是占位符，需要填成 CATH-ddG 公开仓库的确切地址。

---

## 如果我们后续把路径改成完全可移植

这些硬编码来自研究过程，不是设计选择。更彻底的做法是引入 `scripts/_paths.py`
统一解析根目录——我尝试过自动改写，但 110 个文件中有多个含
`from __future__ import`，自动注入会破坏语法位置，风险高于收益，因此改为本方案。
若你希望我逐个文件手工改成可移植形式（需要逐个验证可运行性），告诉我即可。
