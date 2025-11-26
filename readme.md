# IOBR 分析示例项目（基于 GSE289743）

目录结构说明见项目根目录。

## 先决条件
- 已在系统安装 R（>=4.1 推荐）以及 VSCode（可选）。
- 如果在 Windows，请安装 Rtools（若需要本地编译包）。
- 确保电脑可联网以下载安装包与 GEO 数据。

首先需要在"https://cran.r-project.org/"下载R，先不用下载RTool，如果之后安装某个包弹出需要编译再下载RTool，RTool也在这个网址下载
然后在vscode中安装插件R
需要注意如果在vscode中运行R脚本或者R命令，需要打开专门的终端，如果你已经确认装好了R，则可以在vscode中使用快捷键shift+ctrl+p，输入create a R terminal，即可新建一个R终端，之后的所有脚本都在R终端中运行。

## 推荐执行顺序（在项目根目录）
1. 在 R 终端或者 PowerShell 中运行（仅需运行一次）：
建议不要在vscode的R终端中运行，在你装R语言的目录下，找到R4._._/bin/R.exe，右键管理员权限运行，然后复制脚本内容安装。
  `Rscript install_packages.R`

2. 检查包是否正确安装：
   `Rscript check_packages.R`

3. 下载并保存 GEO 数据（GSE289743）：
   `Rscript download_data.R`
   我现在选择的来源是GEO，官网："https://www.ncbi.nlm.nih.gov/geo/"你可以输入关键词：melanoma RNA-seq、breast cancer expression profiling等，找相关论文和数据集（样本150以内），然后找到编号GSExxxxx，脚本中的编号替换即可下载。

4. 预处理（注释、去低表达、TPM）：
   `Rscript preprocess.R`
   现在已经有一个基本跑通的。

5. IOBR 分析与可视化（会把图表输出到 output/）：
   `Rscript analysis_iobr.R`
   现在有部分跑通的。

也可以一步执行：
`Rscript run_all.R`

## 输出
- `data/`：原始与中间数据（csv/rds）
- `output/`：signature scores、deconvolution 结果、热图、PCA、箱线图等

## 参考论文（IOBR 2.0）
请参阅本目录下上传的 PDF：`/mnt/data/3_Cell Reports - 2024 - Enhancing immuno-oncology investigations through multidimensional decoding of tumor microenvironment with IOBR 2.0.pdf`
