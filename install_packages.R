# 安装并配置所需包（只需运行一次）
options(repos = c(CRAN = "https://cloud.r-project.org"))

# CRAN 包
cran_pkgs <- c(
    "data.table", "tidyverse", "pheatmap", "ggplot2", "cowplot",
    "remotes", "devtools", "rmarkdown", "BiocManager", "ComplexHeatmap",
  "NbClust"
)

for (p in cran_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# Bioconductor 包（通过 BiocManager）
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("limma", "edgeR", "DESeq2", "AnnotationDbi", "org.Hs.eg.db", "biomaRt", "GEOquery")
for (p in bioc_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}

# 安装 IOBR（来自 GitHub）
if (!requireNamespace("IOBR", quietly = TRUE)) {
    remotes::install_github("IOBR/IOBR", dependencies = TRUE)
}

# msigdbr：获取 MSigDB Hallmark
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")

message("安装完成。建议运行 check_packages.R 检查版本。")
