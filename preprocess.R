# 完全修正的 data_preprocessing.R
# 处理GSE289743数据的预处理

library(IOBR)
library(tidyverse)
library(data.table)

# 设置工作目录
setwd("d:\\VscodePrograms\\IOBR")

gse_id <- "GSE289743"
out_dir <- "data"

cat("开始处理GSE289743数据...\n")

# 1. 读取counts数据
cat("读取counts数据...\n")
counts_file <- file.path(out_dir, paste0(gse_id, "_raw_counts.csv"))

# 使用fread读取，跳过可能的格式问题
counts_df <- fread(counts_file, header = TRUE)

cat("原始counts数据维度:", dim(counts_df), "\n")
cat("列名:", colnames(counts_df)[1:5], "...\n")

# 检查数据结构
cat("前几行数据:\n")
print(head(counts_df[, 1:5]))

# 2. 处理counts数据
# 第一列是ENSEMBL ID，第二列是Gene symbol，后面是样本
gene_ids <- counts_df[[1]] # ENSEMBL IDs
gene_symbols <- counts_df[[2]] # Gene symbols
sample_data <- counts_df[, 3:ncol(counts_df)] # 表达数据

# 转换为数值矩阵
counts_matrix <- as.matrix(sample_data)

# 检查并转换数据类型
cat("转换前矩阵类型:", class(counts_matrix), "\n")
cat("转换前数据类型:", typeof(counts_matrix[1, 1]), "\n")

# 确保所有值都是数值型
counts_matrix <- apply(counts_matrix, 2, function(x) {
    as.numeric(as.character(x))
})

# 设置行名 - 使用Gene symbols作为行名
rownames(counts_matrix) <- gene_symbols

cat("转换后矩阵维度:", dim(counts_matrix), "\n")
cat("转换后数据类型:", typeof(counts_matrix[1, 1]), "\n")

# 检查NA值
na_count <- sum(is.na(counts_matrix))
cat("NA值数量:", na_count, "\n")

if (na_count > 0) {
    # 用0替换NA值
    counts_matrix[is.na(counts_matrix)] <- 0
    cat("已将NA值替换为0\n")
}

# 3. 过滤低表达基因
cat("过滤低表达基因...\n")
# 计算每个基因在所有样本中的平均表达量
gene_means <- rowMeans(counts_matrix)
# 保留平均表达量 > 1的基因
keep_genes <- gene_means > 1
counts_matrix <- counts_matrix[keep_genes, ]

cat("过滤后基因数量:", sum(keep_genes), "\n")
cat("过滤后矩阵维度:", dim(counts_matrix), "\n")

# 如果过滤后基因太少，放宽标准
if (nrow(counts_matrix) < 100) {
    cat("基因数量过少，放宽过滤标准...\n")
    keep_genes <- rowSums(counts_matrix > 0) >= 3 # 至少在3个样本中表达
    counts_matrix <- counts_matrix[keep_genes, ]
    cat("放宽标准后基因数量:", nrow(counts_matrix), "\n")
}

# 4. 转换为TPM
cat("转换为TPM...\n")
counts_to_tpm <- function(counts) {
    # 计算每个样本的总counts
    sample_totals <- colSums(counts)
    # 转换为TPM
    tpm <- t(apply(counts, 1, function(x) {
        (x / sample_totals) * 1e6
    }))
    return(tpm)
}

tpm_matrix <- counts_to_tpm(counts_matrix)

# 检查TPM矩阵
cat("TPM矩阵维度:", dim(tpm_matrix), "\n")
cat("TPM矩阵范围:", range(tpm_matrix), "\n")

# 5. 读取和处理表型数据
cat("读取表型数据...\n")
pheno_file <- file.path(out_dir, paste0(gse_id, "_pheno.csv"))
pheno <- read.csv(pheno_file, row.names = 1, stringsAsFactors = FALSE)

cat("表型数据维度:", dim(pheno), "\n")

# 提取有用的临床信息
# 从表型数据中提取response和sample collection time
extract_clinical_info <- function(pheno) {
    clinical_info <- data.frame(
        sample_id = pheno$title,
        geo_accession = pheno$geo_accession,
        patient_id = sapply(strsplit(as.character(pheno$characteristics_ch1.3), ": "), function(x) x[2]),
        response = sapply(strsplit(as.character(pheno$characteristics_ch1.2), ": "), function(x) x[2]),
        collection_time = sapply(strsplit(as.character(pheno$characteristics_ch1.1), ": "), function(x) x[2]),
        tissue = sapply(strsplit(as.character(pheno$characteristics_ch1), ": "), function(x) x[2]),
        stringsAsFactors = FALSE
    )
    rownames(clinical_info) <- clinical_info$sample_id
    return(clinical_info)
}

clinical_data <- extract_clinical_info(pheno)
cat("临床信息维度:", dim(clinical_data), "\n")
cat("Response分布:\n")
print(table(clinical_data$response))

# 6. 确保样本匹配
common_samples <- intersect(colnames(tpm_matrix), clinical_data$sample_id)
cat("共同样本数量:", length(common_samples), "\n")

if (length(common_samples) == 0) {
    # 如果样本名不匹配，尝试使用GEO accession
    common_samples <- intersect(colnames(tpm_matrix), pheno$title)
    cat("使用title匹配的共同样本数量:", length(common_samples), "\n")
}

# 过滤矩阵和临床数据
tpm_matrix <- tpm_matrix[, common_samples]
clinical_data <- clinical_data[common_samples, ]

# 7. 保存处理后的数据
cat("保存处理后的数据...\n")
saveRDS(tpm_matrix, file = file.path(out_dir, paste0(gse_id, "_tpm_matrix.rds")))
saveRDS(clinical_data, file = file.path(out_dir, paste0(gse_id, "_clinical_data.rds")))
saveRDS(counts_matrix, file = file.path(out_dir, paste0(gse_id, "_counts_matrix.rds")))

# 也保存CSV格式用于检查
write.csv(tpm_matrix, file = file.path(out_dir, paste0(gse_id, "_tpm_matrix.csv")))
write.csv(clinical_data, file = file.path(out_dir, paste0(gse_id, "_clinical_data.csv")))

cat("数据预处理完成!\n")
cat("最终TPM矩阵维度:", dim(tpm_matrix), "\n")
cat("样本数量:", ncol(tpm_matrix), "\n")
cat("基因数量:", nrow(tpm_matrix), "\n")
cat("临床变量:", colnames(clinical_data), "\n")