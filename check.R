# 将下面整段复制到 R（工作目录应为 d:/VscodePrograms/IOBR）
gse_id <- "GSE289743"
data_dir <- "data"

# 读取中间产物（你之前脚本保存的）
counts_filt_path <- file.path(data_dir, paste0(gse_id, "_counts_filtered_fixed2.rds"))
expr_inter_path <- file.path(data_dir, paste0(gse_id, "_expr_num_intermediate.rds"))
gene_len_file <- file.path(data_dir, paste0(gse_id, "_gene_lengths.rds"))
gene_len_final_file <- file.path(data_dir, paste0(gse_id, "_gene_lengths_final.rds"))

# 优先读取现有对象（替换为你实际保存的路径）
if (file.exists(file.path(data_dir, paste0(gse_id, "_counts_filtered.rds")))) {
    counts_mat <- readRDS(file.path(data_dir, paste0(gse_id, "_counts_filtered.rds")))
} else if (file.exists(counts_filt_path)) {
    counts_mat <- readRDS(counts_filt_path)
} else if (file.exists(expr_inter_path)) {
    # 如果只有中间矩阵，使用之
    tmp <- readRDS(expr_inter_path)
    # tmp 可能是带行名的 matrix 或 data.frame
    counts_mat <- as.matrix(tmp)
} else {
    stop("未找到 counts_filtered 或 intermediate 文件，请确认 data/ 目录下文件名。")
}

# 尝试读取 gene_lengths（不同阶段有不同文件）
if (file.exists(gene_len_final_file)) {
    gene_lengths <- readRDS(gene_len_final_file)
} else if (file.exists(gene_len_file)) {
    gene_lengths <- readRDS(gene_len_file)
} else {
    # 如果没有保存的 gene_lengths，请提示并停止（否则需重新跑 biomaRt）
    stop("找不到 gene_lengths 的 rds 文件，请先运行 preprocess 的 biomaRt 部分或告诉我让脚本重新查询。")
}

cat("counts 矩阵维度:", dim(counts_mat), "\n")
cat("gene_lengths 长度:", length(gene_lengths), " 名称示例: ", paste(head(names(gene_lengths), 5), collapse = ", "), "\n")
cat("前 5 gene_lengths 值:", paste(head(gene_lengths, 5), collapse = ", "), "\n")
cat("typeof(gene_lengths):", typeof(gene_lengths), "\n")
cat("class(gene_lengths):", class(gene_lengths), "\n")

# 1) 强制把 gene_lengths 转为 numeric（先保存原样以便调试）
gene_lengths_orig <- gene_lengths
# 如果 gene_lengths 是 data.frame 或 matrix，尝试提取合适列
if (is.data.frame(gene_lengths) || is.matrix(gene_lengths)) {
    # 取第一列为长度（根据你的保存格式调整）
    gene_lengths <- as.numeric(as.character(gene_lengths[, 1]))
    # 如果 names 丢失，尝试恢复
    if (is.null(names(gene_lengths)) && !is.null(rownames(gene_lengths_orig))) {
        names(gene_lengths) <- rownames(gene_lengths_orig)
    }
} else {
    # 直接转 numeric（会把非数值变为 NA）
    gene_lengths <- as.numeric(as.character(gene_lengths))
}

# 2) 确保 names(gene_lengths) 与 counts_mat 的行名对齐
if (is.null(names(gene_lengths))) {
    warning("gene_lengths 没有 names，尝试用 counts_mat 的行名顺序直接映射 gene_lengths（谨慎）")
    # 若长度一致，直接设置 names
    if (length(gene_lengths) == nrow(counts_mat)) {
        names(gene_lengths) <- rownames(counts_mat)
    }
}

# 对齐顺序：按 counts_mat 行顺序取 gene_lengths
gene_lengths_aligned <- gene_lengths[rownames(counts_mat)]
# 检查 NA 数量与零值
na_cnt <- sum(is.na(gene_lengths_aligned))
zero_cnt <- sum(gene_lengths_aligned == 0, na.rm = TRUE)
cat("对齐后 gene_lengths NA 数量:", na_cnt, "  零值数量:", zero_cnt, "\n")

# 3) 如果 NA/0 比例很高，则列出部分示例供排查
if (na_cnt > 0) {
    cat("示例：前 20 NA 的基因名：\n")
    print(head(rownames(counts_mat)[which(is.na(gene_lengths_aligned))], 20))
}

# 4) 去掉没有长度或长度为 0 的基因（必须）
keep_len <- !is.na(gene_lengths_aligned) & gene_lengths_aligned > 0
cat("保留有长度基因数:", sum(keep_len), "（占原基因数", nrow(counts_mat), "）\n")
counts_mat2 <- counts_mat[keep_len, , drop = FALSE]
gene_lengths2 <- gene_lengths_aligned[keep_len]

if (nrow(counts_mat2) == 0) stop("无基因保留用于 TPM 转换（gene_lengths 全为 NA 或 0）。请把上面 NA 的基因名贴给我，我来判断应该用哪种 ID 查询策略。")

# 5) 如果 counts_mat2 中仍为 character（防御），强制转换 numeric
if (!is.numeric(counts_mat2)) {
    counts_mat2 <- apply(counts_mat2, 2, function(x) as.numeric(as.character(x)))
    rownames(counts_mat2) <- rownames(counts_mat)[keep_len]
}

# 6) 过滤低表达（可调整阈值）
min_pct <- 0.05
min_samps_min <- 2
min_keep <- max(ceiling(min_pct * ncol(counts_mat2)), min_samps_min)
keep_expr <- rowSums(counts_mat2 >= 1, na.rm = TRUE) >= min_keep
cat("表达过滤阈值：至少在", min_keep, "个样本中 counts >=1. 过滤后基因数:", sum(keep_expr), "\n")
counts_final <- counts_mat2[keep_expr, , drop = FALSE]
gene_lengths_final <- gene_lengths2[keep_expr]

if (nrow(counts_final) < 50) {
    warning("用于 TPM 转换的基因数较少（<50），请确认是否合适。")
}

# 7) 计算 TPM（使用简单函数）
rpk <- counts_final / (gene_lengths_final / 1000) # 如果这里出错说明 gene_lengths_final 仍非数值
tpm <- t(t(rpk) / colSums(rpk, na.rm = TRUE) * 1e6)

# 8) 保存结果
tpm_rds <- file.path(data_dir, paste0(gse_id, "_tpm_fixed3.rds"))
tpm_csv <- file.path(data_dir, paste0(gse_id, "_tpm_fixed3.csv"))
saveRDS(tpm, tpm_rds)
write.csv(tpm, tpm_csv, row.names = TRUE)
cat("已保存 TPM：", tpm_rds, "\n")
cat("TPM 矩阵维度：", dim(tpm), "\n")
