# IOBR 分析示例项目（基于 GSE289743）

目录结构说明见项目根目录。

## 先决条件
- 已在系统安装 R（>=4.1 推荐）以及 VSCode（可选）。
- 如果在 Windows，请安装 Rtools（若需要本地编译包）。
- 确保电脑可联网以下载安装包与 GEO 数据。

## 推荐执行顺序（在项目根目录）
1. 在 R 终端或者 PowerShell 中运行（仅需运行一次）：
   `Rscript install_packages.R`

2. 检查包是否正确安装：
   `Rscript check_packages.R`

3. 下载并保存 GEO 数据（GSE289743）：
   `Rscript download_data.R`

4. 预处理（注释、去低表达、TPM）：
   `Rscript preprocess.R`

5. IOBR 分析与可视化（会把图表输出到 output/）：
   `Rscript analysis_iobr.R`

也可以一步执行：
`Rscript run_all.R`

## 输出
- `data/`：原始与中间数据（csv/rds）
- `output/`：signature scores、deconvolution 结果、热图、PCA、箱线图等

## 参考论文（IOBR 2.0）
请参阅本目录下上传的 PDF：`/mnt/data/3_Cell Reports - 2024 - Enhancing immuno-oncology investigations through multidimensional decoding of tumor microenvironment with IOBR 2.0.pdf`
