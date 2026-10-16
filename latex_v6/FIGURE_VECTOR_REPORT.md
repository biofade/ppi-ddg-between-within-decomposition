# FIGURE_VECTOR_REPORT

用你提供的 `make_all_figures.py` 重新生成图形，并把**真矢量 PDF** 装入论文。

## 结论

```
VECTOR_FIGURES      = 9/9 真矢量（Figure 6 / S2 各含 2 个位图元素，见下）
FIGURE_COUNT        = 6 (main) + 3 (supplement)
BUILD               = fatal=0, citation_undef=0, reference_undef=0
MAIN_PAGES          = 10   SUPPLEMENT_PAGES = 4
MAIN_PDF_SIZE       = 384 KB（此前的位图包装版为 8.0 MB）
```

## 之前的判断是错的，这里更正

我在上一轮说「你给的是位图，无法矢量化」。**这个结论错了** —— 我当时只有位图，
没有生成代码。你这次给了 `make_all_figures.py`，直接运行即可得到真矢量。

同时发现一个更重要的事实：**你之前发我的 PNG 是脚本输出的裁切/压缩版**。
对比 Figure 1：

| | 尺寸 | 内容 |
|---|---|---|
| 脚本输出（现采用） | 1879×1190 | 完整：`True ΔΔG` / `Predicted ΔΔG` 框、`14 families`、精确恒等式、所有标注齐全 |
| 你此前发的 PNG | 1916×821 | 上下边缘被裁掉，`14 families`、`True ΔΔG` 标题等缺失 |

这也解释了我先前"同一张图却对不上"的困惑：它们同源，但前者完整、后者被裁切。

## 生成与安装

`make_all_figures.py` 已审查：仅本地文件写入（svg/pdf/png），无网络、无 subprocess。

```powershell
python make_all_figures.py --figure all     # 输出 outputs/*.{svg,pdf,png}
python scripts/v8_1/72_install_true_vectors.py   # 装入 figures_v5/
```

名称对应：

| 脚本输出 | 论文用名 |
|---|---|
| `Figure1_large_manuscript.pdf` | `Fig1_overview_EN.pdf` |
| `Figure2.pdf` | `Fig2_controlled_split_collapse_EN.pdf` |
| `Figure3.pdf` | `Fig3_target_weighting_sensitivity_EN.pdf` |
| `Figure4.pdf` | `Fig4_fixed_K_geometry_EN.pdf` |
| `Figure5.pdf` | `Fig5_retrospective_remedy_EN.pdf` |
| `Figure6.pdf` | `Fig6_external_methods_EN.pdf` |
| `FigureS1.pdf` | `FigS1_matched_universe_strata_EN.pdf` |
| `FigureS2.pdf` | `FigS2_capability_heatmap_provenance_clean_EN.pdf` |
| `FigureS3.pdf` | `FigS3_pool_sensitivity_rank_agreement_EN.pdf` |

## 矢量性核验

从绘图脚本产出的 PDF 中读取内容类型：

| 图 | PDF 页面 (pt) | 矢量对象数 | 位图元素 |
|---|---|---|---|
| Fig 1 | 1020.1×646.2 | 139 | 0 |
| Fig 2 | 861.6×442.6 | 85 | 0 |
| Fig 3 | 980.1×323.1 | 123 | 0 |
| Fig 4 | 802.9×427.3 | 69 | 0 |
| Fig 5 | 623.2×392.7 | 61 | 0 |
| Fig 6 | 774.9×331.4 | 65 | 2 |
| Fig S1 | 727.9×398.8 | 58 | 0 |
| Fig S2 | 622.8×360.0 | 29 | 2 |
| Fig S3 | 610.9×373.5 | 37 | 0 |

编译后的论文 PDF 内：

```
manuscript_v6.pdf    总矢量对象 558   位图 2   （p3–p7 含图）
supplementary_v6.pdf 总矢量对象 129   位图 2   （p3–p4 含图）
```

**Figure 6 与 Figure S2 各含 2 个位图元素**：这两个图带有 `imshow` 生成的色标
（colorbar），matplotlib 对 `imshow` 数据栅格化是正常行为。主体（热图单元格描边、
文字、坐标轴、色标刻度）均为矢量。这是绘图脚本本身的输出方式，我没有改动它。

## 排版

新增等比上限宏，宽 ≤ 文本宽、高 ≤ 正文高 82%，`keepaspectratio` 保持长宽比：

```latex
\newcommand{\figfit}[2][\textwidth]{%
  \includegraphics[width=#1,height=0.82\textheight,keepaspectratio]{#2}}
```

图仍位于正文 p3–p7，紧邻首次引用。

## 备份

- 脚本输出原件：`F:\dsh\cache\tmp\figuregen\outputs\`
- 上一轮位图包装版：`F:\dsh\cache\tmp\figs_raster_backup\`
- 最早那套旧矢量图：`F:\dsh\cache\tmp\figs_backup\`

## 复现

```powershell
cd F:\dsh\cache\tmp\figuregen
python make_all_figures.py --figure all
cd F:\DSH_Galatea\Affinity_Intelligence
python scripts/v8_1/72_install_true_vectors.py
python scripts/v8_1/53_build_latex_v6.py
cd latex_v6; powershell -File build.ps1
```
