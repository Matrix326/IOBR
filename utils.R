# 通用工具函数
library(data.table)

# 将 counts data.table-> matrix 的辅助函数
dt_to_matrix <- function(dt) {
    # 假设第一列为 gene id, 其余为样本列
    gene_col <- names(dt)[1]
    m <- as.matrix(dt[, -1, with = FALSE])
    rownames(m) <- dt[[gene_col]]
    return(m)
}

# 去掉 ensembl 版本号（例如 ENSG000001234.5 -> ENSG000001234）
remove_ensembl_version <- function(ens_ids) {
    sapply(strsplit(ens_ids, split = "\\."), `[`, 1)
}

# 简单 counts -> tpm 转换（若 IOBR::count2tpm 可用则优先使用）
counts_to_tpm_simple <- function(counts, gene_length) {
    # counts: genes x samples matrix (numeric)
    # gene_length: vector of gene lengths (bp), names match rownames(counts)
    if (is.null(gene_length)) stop("需要基因长度向量")
    if (length(gene_length) != nrow(counts)) stop("基因长度与基因数不匹配")
    rpk <- counts / (gene_length / 1000)
    tpm <- t(t(rpk) / colSums(rpk, na.rm = TRUE) * 1e6)
    return(tpm)
}
