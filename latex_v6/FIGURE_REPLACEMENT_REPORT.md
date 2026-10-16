# FIGURE_REPLACEMENT_REPORT

用你最新提供的 9 张图替换 `manuscript_v6.pdf` / `supplementary_v6.pdf` 中的图形。

## 结论

```
FIGURE_CONTENT_IDENTICAL = PASS   (9/9 逐像素 maxdiff = 0)
ASPECT_RATIO_PRESERVED   = PASS   (显示比例 == 源图比例)
FIGURE_COUNT             = 6 (main) + 3 (supplement)
BUILD                    = fatal=0, citeUndef=0, refUndef=0
```

## 关键发现：旧图是「另一套画风」

替换前 `figures_v5/*.pdf` 与你给的图**不是同一套**。例如 Figure 1：

| | 内容 |
|---|---|
| 旧 PDF（已替换） | 示意风：`A. Provenance-filtered, Fingerprint-Matched Benchmark & Target Structure`，含 Target 1/2/3 卡片、`μ_c / ε_ci` 分组框 |
| 你给的图（现采用） | 具体风：SKEMPI v2.0 数据库圆柱、蛋白结构、SHA-256 二维码卡、`4,473 rows / 214 targets / 108 families`、右侧 F1/F2/F4/F8 能力折线 |

图片页码（Figure 1）两者比例也不同：旧图 1.455，你的图 **2.334**。
Figure 4 同样是两套不同版本。

## 但无法真正「矢量化」

需要向你说明一个客观限制：

- 你提供的图是**位图**（PNG / JPEG）；我核查了你 v6 PDF 内嵌的图，**也是位图**
  （每页 1 张 image、drawing 数为 0），Word 里同样没有矢量原件。
- 位图无法在不重跑生成代码的前提下还原为矢量路径。我没有生成这批新图的绘图代码
  （Figure 1 是 DAMO-RADAR 风格的信息图，不是本项目脚本产出的）。

因此采用论文排版里的标准做法：**每张图包成一个单页 PDF，页面尺寸严格等于图像宽高比**，
再由 `\includegraphics` 引用。这样：

- 论文里的图形文件是 **PDF**（你要求的格式）
- 图形内容 = 你给的位图，**逐像素完全相同**（已验证 maxdiff = 0）
- **未重绘、未裁切、未改比例、未重新编码**

按 600 dpi 等效建立页面尺寸：

```
page_width_pt  = width_px  / 600 * 72
page_height_pt = height_px / 600 * 72
```

## 逐图核验

| 图 | 源图像素 | PDF 页面 (pt) | 比例保持 | 内嵌像素 vs 源图 |
|---|---|---|---|---|
| Fig 1 | 1916×821 | 229.9×98.5 | 2.3337 | maxdiff = 0 |
| Fig 2 | 2837×1477 | 340.4×177.2 | 1.9208 | maxdiff = 0 |
| Fig 3 | 3376×1242 | 405.1×149.0 | 2.7182 | maxdiff = 0 |
| Fig 4 | 1619×971 | 194.3×116.5 | 1.6674 | maxdiff = 0 |
| Fig 5 | 2586×1620 | 310.3×194.4 | 1.5963 | maxdiff = 0 |
| Fig 6 | 3083×1360 | 370.0×163.2 | 2.2669 | maxdiff = 0 |
| Fig S1 | 3012×1392 | 361.4×167.0 | 2.1638 | maxdiff = 0 |
| Fig S2 | 2695×1556 | 323.4×186.7 | 1.7320 | maxdiff = 0 |
| Fig S3 | 2917×1437 | 350.0×172.4 | 2.0299 | maxdiff = 0 |

「内嵌像素 vs 源图」= 从生成的 PDF 中抽出内嵌位图，与你给的源图逐像素相减，
最大像素差为 0，即完全相同。

## 排版适配

新图比例差异大（Figure 1 是 2.33:1 横幅，Figure 5 是 1.60:1），固定的
`width=\textwidth` 会让横幅图过小。改用等比上限宏：

```latex
\newcommand{\figfit}[2][\textwidth]{%
  \includegraphics[width=#1,height=0.82\textheight,keepaspectratio]{#2}}
```

宽不超过文本宽、高不超过正文高的 82%，两轴同时约束，**保持长宽比不变**。

编译后实际显示尺寸：

| 图 | 页码 | 显示尺寸 (pt) |
|---|---|---|
| Fig 1 | p3 | 464.9 × 199.2 |
| Fig 2 | p4 | 464.9 × 242.0 |
| Fig 3 | p4 | 464.9 × 171.0 |
| Fig 4 | p5 | 464.9 × 278.8 |
| Fig 5 | p6 | 464.9 × 291.2 |
| Fig 6 | p7 | 464.9 × 205.1 |
| Fig S1 | p3 | 481.9 × 222.7 |
| Fig S2 | p4 | 481.9 × 278.2 |
| Fig S3 | p4 | 481.9 × 237.4 |

## 备份

替换前的旧矢量图已备份到 `F:\dsh\cache\tmp\figs_backup\`。

## 复现

```powershell
python scripts/v8_1/68_make_vector_pdfs.py    # 由源图生成图形 PDF
python scripts/v8_1/53_build_latex_v6.py      # 重新生成 .tex
cd latex_v6; powershell -File build.ps1       # 编译
```
