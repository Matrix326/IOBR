# 检查关键包是否安装并显示版本
pkgs <- c("IOBR", "GEOquery", "data.table", "msigdbr", "biomaRt", "pheatmap", "ggplot2", "ComplexHeatmap")

for (p in pkgs) {
    cat("----", p, "----\n")
    if (!requireNamespace(p, quietly = TRUE)) {
        cat(p, "未安装\n")
    } else {
        v <- as.character(utils::packageVersion(p))
        cat(p, "version:", v, "\n")
    }
}
