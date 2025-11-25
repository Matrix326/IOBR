# debug_fix_counts.R — 诊断并修复 counts -> tpm 的问题
setwd("d:/VscodePrograms/IOBR")
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# 加载 utils（如果有）
if (file.exists("utils.R")) {
    source("utils.R")
} else {
    dt_to_matrix <- function(dt) {
        gene_col <- names(dt)[1]
        m <- as.matrix(dt[, -1, with = FALSE])
        rownames(m) <- dt[[gene_col]]
        return(m)
    }
    remove_ensembl_version <- function(ens_ids) sapply(strsplit(ens_ids, split = "\\."), `[`, 1)
    counts_to_tpm_simple <- function(counts, gene_length) {
        rpk <- counts / (gene_length / 1000)
        tpm <- t(t(rpk) / colSums(rpk, na.rm = TRUE) * 1e6)
        return(tpm)
    }
}

gse_id <- "GSE289743"
raw_rds <- file.path("data", paste0(gse_id, "_raw_counts.rds"))
if (!file.exists(raw_rds)) stop("找不到 raw counts rds: ", raw_rds)

# 1. 读取原始 data.table 并显示前几行（原始格式检查）
counts_dt <- readRDS(raw_rds)
cat("counts_dt class:", class(counts_dt), "\n")
cat("counts_dt 列名（前20）：\n")
print(names(counts_dt)[1:min(20, length(names(counts_dt)))])
cat("counts_dt 前 6 行（显示 1~6 列）:\n")
print(head(counts_dt[, 1:min(6, ncol(counts_dt)), with = FALSE], 10))

# 2. 若文件是 data.table，检查第一列名与前几行是否是 header 重复（例如第一行再次包含列名）
first_col <- names(counts_dt)[1]
cat("第一列名:", first_col, "\n")
# 查看第一列若包含 "gene" 等关键字
cat("第一列前 6 值示例：", paste(head(counts_dt[[first_col]], 6), collapse = " | "), "\n")

# 3. 如果发现第一行是列名的重复（例如第一行值等于列名的样式），需要删掉该行
# 检查是否有一整行等于列名（常见情况）
maybe_header_row <- which(apply(counts_dt, 1, function(r) all(as.character(r[1]) == names(counts_dt)[1])))
# 简单检测：若第一行第1列等于 "Gene" 或 "gene", 或第一行与列名相似
if (tolower(as.character(counts_dt[[first_col]][1])) %in% c("gene", "genes", "ensembl", "gene_id")) {
    cat("检测到第一行可能是重复 header（比如 'Gene'），将删除第一行。\n")
    counts_dt <- counts_dt[-1, ]
}

# 4. 进一步检查每列是否含非数字（排查导致 as.numeric -> NA 的原因）
# 先查看每列中非数字（含字母）的比例
col_nonnum_rate <- sapply(names(counts_dt)[-1], function(cn) {
    vals <- as.character(counts_dt[[cn]])
    # 非数字视作不能被转换为数值（排除 NA）
    nonnum <- sum(is.na(suppressWarnings(as.numeric(vals))))
    return(nonnum / length(vals))
})
# 报告样本列非数字比例高的前 10
ord <- order(col_nonnum_rate, decreasing = TRUE)
high_nonnum <- data.frame(sample = names(col_nonnum_rate)[ord], nonnum_rate = col_nonnum_rate[ord])
print(head(high_nonnum, 10))
cat("如果某些样本列非数字率 ==1，说明该列被误解析为文本或包含多余注释。\n")

# 5. 如果出现某列非数字率接近 1（比如第一个样本列被当作列名），把该列移除或修正
problem_cols <- names(col_nonnum_rate)[which(col_nonnum_rate > 0.9)]
if (length(problem_cols) > 0) {
    cat("检测到以下列非数字率>0.9，可能是误解析，打印样例查看：\n")
    print(problem_cols)
    # 打印这些列的前 6 个值看看
    for (cn in problem_cols[1:min(length(problem_cols), 5)]) {
        cat("列", cn, "前 6 值：", paste(head(as.character(counts_dt[[cn]]), 6), collapse = " | "), "\n")
    }
    # 如果 problem_cols 包含名为 "Gene" 的列（或首列当作 sample），提示用户人工检查
    if ("Gene" %in% problem_cols) {
        cat("样本列包含 'Gene'，可能是 header 混入列，请手动检查 counts 文件或告诉我输出以便我帮你处理。\n")
    }
}

# 6. 强制清洗并转换：保证第一列为 gene id（字符），其余列均转换 numeric（转换失败会变 NA）
gene_col <- names(counts_dt)[1]
expr_mat <- as.matrix(counts_dt[, -1, with = FALSE])
rownames_expr <- as.character(counts_dt[[gene_col]])
# If any duplicate gene ids, preserve as is (will aggregate later if needed)
colnames(expr_mat) <- names(counts_dt)[-1]

# 强制 numeric 转换（逐列）
expr_num <- apply(expr_mat, 2, function(col) suppressWarnings(as.numeric(as.character(col))))
rownames(expr_num) <- rownames_expr

# 报告转换后 NA 情况（逐列与逐行）
na_rows <- sum(rowSums(is.na(expr_num)) > 0)
na_cols <- sum(colSums(is.na(expr_num)) > 0)
cat("转换 numeric 后：含 NA 的基因行数（至少一个 NA）：", na_rows, "\n")
cat("转换 numeric 后：含 NA 的样本列数（至少一个 NA）：", na_cols, "\n")
if (na_cols > 0) {
    cat("含 NA 的样本（示例前 10）:\n")
    print(names(which(colSums(is.na(expr_num)) > 0))[1:10])
}

# 7. 若某些样本列全部为 NA（被误解析为文本），将其移除以免影响后续
all_na_cols <- names(which(colSums(!is.na(expr_num)) == 0))
if (length(all_na_cols) > 0) {
    cat("移除全为 NA 的样本列：", paste(all_na_cols, collapse = ", "), "\n")
    expr_num <- expr_num[, setdiff(colnames(expr_num), all_na_cols), drop = FALSE]
}

# 8. 去掉行名可能的空值或 NA
valid_genes <- !is.na(rownames(expr_num)) & rownames(expr_num) != ""
expr_num <- expr_num[valid_genes, , drop = FALSE]

cat("当前矩阵维度（清洗后）:", dim(expr_num), "\n")
cat("前 5 基因 id:", paste(head(rownames(expr_num), 5), collapse = ", "), "\n")
cat("前 5 样本名:", paste(head(colnames(expr_num), 5), collapse = ", "), "\n")

# 9. 去掉 Ensembl 版本号（如 ENSG00000000003.14 -> ENSG00000000003）
if (any(grepl("^ENSG", rownames(expr_num)))) {
    rownames(expr_num) <- remove_ensembl_version(rownames(expr_num))
    cat("已去除 Ensembl 版本号，检查前 5 ids:", paste(head(rownames(expr_num), 5), collapse = ", "), "\n")
}

# 10. 聚合重复基因名（若同一 gene symbol/ensembl 出现多次），用行均值合并
if (any(duplicated(rownames(expr_num)))) {
    cat("检测到重复基因名，使用行均值合并重复基因（aggregate）...\n")
    expr_df2 <- as.data.frame(expr_num)
    expr_df2$gene <- rownames(expr_num)
    agg <- aggregate(. ~ gene, data = expr_df2, FUN = mean)
    rownames(agg) <- agg$gene
    expr_num <- as.matrix(agg[, -ncol(agg)])
    cat("合并后矩阵维度:", dim(expr_num), "\n")
}

# 11. 保存一个中间文件，便于回顾
saveRDS(expr_num, file = file.path("data", paste0(gse_id, "_expr_num_intermediate.rds")))
cat("已保存中间矩阵：data/", paste0(gse_id, "_expr_num_intermediate.rds\n"))

# 12. 使用 biomaRt 查询 gene length（优先用 ensembl id，否则用 symbol）
library(biomaRt)
mart <- tryCatch(
    {
        useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    },
    error = function(e) useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

genes_now <- rownames(expr_num)
cat("biomaRt 查询基因数:", length(genes_now), "\n")

# 选择过滤器
if (all(grepl("^ENSG", genes_now))) {
    filter_name <- "ensembl_gene_id"
    attrs <- c("ensembl_gene_id", "hgnc_symbol", "transcript_length")
    query_vals <- genes_now
} else {
    filter_name <- "hgnc_symbol"
    attrs <- c("hgnc_symbol", "ensembl_gene_id", "transcript_length")
    query_vals <- genes_now
}

cat("执行 biomaRt 查询（可能较慢）...\n")
bm <- getBM(attributes = attrs, filters = filter_name, values = query_vals, mart = mart)

if (nrow(bm) == 0) {
    cat("biomaRt 返回 0 行：可能是 ID 类型不匹配。请确认 rownames 是否为 Ensembl (无版本) 或 HGNC symbol。\n")
    cat("示例前 10 gene ids:\n")
    print(head(genes_now, 10))
    stop("biomaRt 未返回注释，停止以便人工检查。你可以把上面的示例 gene ids 发给我，我帮你诊断。")
}

# 聚合 transcript_length 取最大
if (filter_name == "ensembl_gene_id") {
    bm2 <- bm %>%
        group_by(ensembl_gene_id) %>%
        summarize(gene_symbol = first(hgnc_symbol), gene_length = max(transcript_length, na.rm = TRUE))
    idx <- match(genes_now, bm2$ensembl_gene_id)
    gene_lengths <- bm2$gene_length[idx]
    names(gene_lengths) <- genes_now
} else {
    bm2 <- bm %>%
        group_by(hgnc_symbol) %>%
        summarize(ensembl_gene_id = first(ensembl_gene_id), gene_length = max(transcript_length, na.rm = TRUE))
    idx <- match(genes_now, bm2$hgnc_symbol)
    gene_lengths <- bm2$gene_length[idx]
    names(gene_lengths) <- genes_now
}

cat("查询到有 gene_length 的基因数:", sum(!is.na(gene_lengths) & gene_lengths > 0), "\n")
# 移除没有 gene length 的基因
sel_len <- !is.na(gene_lengths) & gene_lengths > 0
expr_num2 <- expr_num[sel_len, , drop = FALSE]
gene_lengths2 <- gene_lengths[sel_len]
cat("保留下来用于 TPM 转换的基因数:", nrow(expr_num2), "\n")

if (nrow(expr_num2) < 50) {
    cat("警告：可用于 TPM 转换的基因数非常少（<50），请检查 biomaRt 匹配是否正确或原始 counts 文件格式是否异常。\n")
}

# 13. 过滤低表达基因（默认至少在 5% 样本中 counts >=1，且最少2个样本）
min_pct <- 0.05
min_samps_min <- 2
min_samps_pct <- max(ceiling(min_pct * ncol(expr_num2)), min_samps_min)
keep_gene <- rowSums(expr_num2 >= 1, na.rm = TRUE) >= min_samps_pct
cat("过滤阈值：至少在", min_samps_pct, "个样本中 counts >= 1. 过滤后基因数:", sum(keep_gene), "\n")
if (sum(keep_gene) < 100) {
    cat("过滤后基因数 < 100，尝试更宽松阈值（至少2个样本）。\n")
    keep_gene <- rowSums(expr_num2 >= 1, na.rm = TRUE) >= min_samps_min
    cat("重新过滤后基因数:", sum(keep_gene), "\n")
}
expr_final <- expr_num2[keep_gene, , drop = FALSE]
gene_lengths_final <- gene_lengths2[keep_gene]

# 14. 转 TPM 并保存
if (nrow(expr_final) == 0) stop("最终没有基因保留（expr_final 行数为0），请回顾上面诊断信息并把输出发给我以便进一步诊断。")
tpm_mat <- counts_to_tpm_simple(as.matrix(expr_final), gene_lengths_final)

saveRDS(expr_final, file = file.path("data", paste0(gse_id, "_counts_filtered_fixed2.rds")))
saveRDS(tpm_mat, file = file.path("data", paste0(gse_id, "_tpm_fixed2.rds")))
write.csv(tpm_mat, file = file.path("data", paste0(gse_id, "_tpm_fixed2.csv")), row.names = TRUE)
saveRDS(gene_lengths_final, file = file.path("data", paste0(gse_id, "_gene_lengths_final.rds")))

cat("修复完成并保存：\n")
cat(" - counts_filtered_fixed2:", file.path("data", paste0(gse_id, "_counts_filtered_fixed2.rds")), "\n")
cat(" - tpm_fixed2:", file.path("data", paste0(gse_id, "_tpm_fixed2.rds")), "\n")
cat(" - tpm csv:", file.path("data", paste0(gse_id, "_tpm_fixed2.csv")), "\n")
