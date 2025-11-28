########################################################################
# preprocess_full_cn.R
# 预处理管线（中文注释与中文输出）
# - 读取 counts 与 pheno
# - 去除 Ensembl 版本号 / 去除 _PAR_Y
# - 可选反 log2(x+1)
# - 低表达过滤（edgeR::filterByExpr 或 简单阈值）
# - 尝试用 IOBR::count2tpm（多种调用方式）
# - 回退到使用 biomaRt 的手动 TPM（若需要）
# - Ensembl -> gene symbol 映射并合并重复（按均值）
# - 保存结果为 rds/csv，并输出摘要
########################################################################

# -------------------- 可由用户修改的参数 -----------------------------
counts_file <- "data/GSE289743_raw_counts.csv" # 你的 counts 文件路径
pheno_file <- "data/GSE289743_pheno.csv" # 你的 pheno 文件路径
outdir <- "results_preprocess" # 输出目录
subset_if_more_than <- 150 # 如果样本数 > 此值，可选择子集
subset_strategy <- c("post_treatment", "responders", "random")[3]
# 子集策略： "post_treatment" / "responders" / "random"
set.seed(123)

# -----------------------------------------------------------------------
# 0. 安装 / 载入需要的包（如缺失将尝试安装）
pkgs_cran <- c("data.table", "ggplot2", "matrixStats")
pkgs_bioc <- c("edgeR", "biomaRt", "IOBR", "GSVA")
need_install_bioc <- pkgs_bioc[!(pkgs_bioc %in% installed.packages()[, "Package"])]
need_install_cran <- pkgs_cran[!(pkgs_cran %in% installed.packages()[, "Package"])]

if (length(need_install_cran) > 0) {
    message("正在安装缺失的 CRAN 包: ", paste(need_install_cran, collapse = ", "))
    install.packages(need_install_cran)
}
if (length(need_install_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    message("正在安装缺失的 Bioconductor 包: ", paste(need_install_bioc, collapse = ", "))
    BiocManager::install(need_install_bioc, ask = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(biomaRt)
    library(IOBR)
    library(matrixStats)
    library(GSVA)
})

# 创建输出目录（若不存在）
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -------------------- 1. 读取输入文件 ------------------------------
message("开始：读取 counts 与 pheno 文件 ...")
counts_dt <- data.table::fread(counts_file)
pheno_dt <- data.table::fread(pheno_file)

message("读取完成：counts 行数:", nrow(counts_dt), " 列数:", ncol(counts_dt))
message("读取完成：pheno 行数:", nrow(pheno_dt))

# 期望 counts_dt 的列为：ID, Gene (可选), sample1, sample2, ...
# 稳健检测 ID/Gene 列，假设第一列为 Ensembl ID
colnames(counts_dt)[1:3] <- colnames(counts_dt)[1:3] # 保持原样，避免校验提示
id_col <- colnames(counts_dt)[1]
gene_col <- ifelse(ncol(counts_dt) >= 2, colnames(counts_dt)[2], NA)
sample_cols <- setdiff(colnames(counts_dt), c(id_col, gene_col))

# 构建表达矩阵
expr_raw <- as.matrix(counts_dt[, ..sample_cols])
rownames(expr_raw) <- as.character(counts_dt[[id_col]])
mode(expr_raw) <- "numeric"

# -------------------- 2. 去除 Ensembl 版本号 与 _PAR_Y ----------------
# 去掉点及其后的版本号（如 ENSG000001234.5 -> ENSG000001234）
rownames(expr_raw) <- sub("\\..*$", "", rownames(expr_raw))

# 如果存在以 _PAR_Y 结尾的基因，移除它们（可选）
par_y_idx <- grepl("_PAR_Y$", rownames(expr_raw), ignore.case = TRUE)
if (any(par_y_idx)) {
    message("检测到并移除后缀为 _PAR_Y 的基因，共计: ", sum(par_y_idx), " 个")
    expr_raw <- expr_raw[!par_y_idx, , drop = FALSE]
} else {
    message("未检测到 _PAR_Y 后缀基因，跳过该步骤。")
}

# -------------------- 3. 样本与 pheno 对齐 -----------------------
message("正在尝试将表达矩阵样本列与 pheno 中的样本 ID 对齐 ...")
possible_sample_cols <- c("title", "geo_accession", "sample", "Sample", "sample_id", "ID")
pheno_samp_col <- intersect(possible_sample_cols, colnames(pheno_dt))
if (length(pheno_samp_col) >= 1) {
    pheno_samp_col <- pheno_samp_col[1]
} else {
    # 回退方法：找一个在 counts 列名中有匹配的 pheno 列
    match_idx <- which(sapply(colnames(pheno_dt), function(x) sum(colnames(expr_raw) %in% pheno_dt[[x]])))
    if (length(match_idx) > 0) pheno_samp_col <- colnames(pheno_dt)[match_idx[1]] else pheno_samp_col <- NA
}
if (is.na(pheno_samp_col)) {
    stop("无法在 pheno 中自动检测样本 ID 列。请检查 pheno 文件并手动设置 pheno_samp_col。")
}
message("使用 pheno 中的样本列: ", pheno_samp_col)

# 对齐样本
common_samples <- intersect(colnames(expr_raw), as.character(pheno_dt[[pheno_samp_col]]))
if (length(common_samples) == 0) {
    # 可能 pheno 的该列值本身即为 counts 的列名
    if (all(as.character(pheno_dt[[pheno_samp_col]]) %in% colnames(expr_raw))) {
        common_samples <- as.character(pheno_dt[[pheno_samp_col]])
    } else {
        stop("counts 的列名与 pheno 中识别的样本名无重叠，请检查命名一致性。")
    }
}
expr_raw <- expr_raw[, common_samples, drop = FALSE]
pheno_dt <- pheno_dt[match(common_samples, pheno_dt[[pheno_samp_col]]), , drop = FALSE]
rownames(pheno_dt) <- pheno_dt[[pheno_samp_col]]
message("样本对齐完成。样本数: ", ncol(expr_raw))

# -------------------- 4. 如果样本数过多，可选子集 -----------------------
if (ncol(expr_raw) > subset_if_more_than) {
    message("样本数量 (", ncol(expr_raw), ") 超过阈值 (", subset_if_more_than, ")，应用子集策略：", subset_strategy)
    if (subset_strategy == "post_treatment") {
        colnm <- grep("sample.collection.time|collection.time|time", tolower(colnames(pheno_dt)), value = TRUE)
        if (length(colnm) > 0) {
            idx <- which(tolower(pheno_dt[[colnm[1]]]) == "post-treatment")
            if (length(idx) > 0) {
                chosen <- rownames(pheno_dt)[idx]
            } else {
                chosen <- sample(rownames(pheno_dt), subset_if_more_than)
            }
        } else {
            chosen <- sample(rownames(pheno_dt), subset_if_more_than)
        }
    } else if (subset_strategy == "responders") {
        colnm <- grep("response", tolower(colnames(pheno_dt)), value = TRUE)
        if (length(colnm) > 0) {
            idx <- which(tolower(pheno_dt[[colnm[1]]]) == "responders")
            if (length(idx) >= subset_if_more_than) {
                chosen <- rownames(pheno_dt)[idx]
            } else {
                chosen <- c(rownames(pheno_dt)[idx], sample(setdiff(rownames(pheno_dt), rownames(pheno_dt)[idx]), subset_if_more_than - length(idx)))
            }
        } else {
            chosen <- sample(rownames(pheno_dt), subset_if_more_than)
        }
    } else {
        chosen <- sample(rownames(pheno_dt), subset_if_more_than)
    }
    expr_raw <- expr_raw[, chosen, drop = FALSE]
    pheno_dt <- pheno_dt[chosen, , drop = FALSE]
    message("子集完成。新的样本数: ", ncol(expr_raw))
}

# -------------------- 5. 判断是否为 log2(x+1) 并反变换 ----------------
message("检测数据是否需要反 log2(x+1) 变换 ...")
max_val <- max(expr_raw, na.rm = TRUE)
are_integers <- all((expr_raw == floor(expr_raw)) | is.na(expr_raw))
message("矩阵最大值: ", round(max_val, 3), "；是否全部为整数: ", are_integers)
expr_counts <- expr_raw
if (!are_integers && max_val < 100) {
    message("数据看起来可能已做 log2(x+1) 变换，正在用 2^x - 1 恢复为 counts ...")
    expr_counts <- (2^expr_raw) - 1
    expr_counts[expr_counts < 0] <- 0
} else {
    message("不做反变换，假定数据为原始 counts。")
}

# -------------------- 6. 低表达基因过滤 -----------------------
message("进行低表达基因过滤 ...")
group_candidates <- grep("response|group|condition|treatment", tolower(colnames(pheno_dt)), value = TRUE)
if (length(group_candidates) > 0) {
    group_col <- group_candidates[1]
    group_vec <- as.factor(pheno_dt[[group_col]])
    message("使用 pheno 列 '", group_col, "' 作为分组用于 filterByExpr。")
    dge <- DGEList(counts = expr_counts)
    keep <- filterByExpr(dge, group = group_vec)
} else {
    message("未检测到明显的分组列，使用简单阈值：保留至少在 2 个样本中 counts >= 10 的基因。")
    keep <- rowSums(expr_counts >= 10, na.rm = TRUE) >= 2
}
message("过滤后保留基因: ", sum(keep), " / ", nrow(expr_counts))
expr_filt <- expr_counts[keep, , drop = FALSE]

# 保存过滤后的 counts
saveRDS(expr_filt, file = file.path(outdir, "counts_filtered.rds"))
write.csv(as.data.frame(expr_filt), file = file.path(outdir, "counts_filtered.csv"), row.names = TRUE)
message("已保存过滤后的 counts 到目录：", normalizePath(outdir))

# -------------------- 7. 尝试用 IOBR::count2tpm（多种调用方式） ----------------
message("尝试使用 IOBR::count2tpm 进行 TPM 转换（尝试多种调用方式）...")
tpm <- NULL
attempts <- list(
    function(x) IOBR::count2tpm(x),
    function(x) IOBR::count2tpm(countMat = x, idType = "Ensembl"),
    function(x) IOBR::count2tpm(counts = x),
    function(x) IOBR::count2tpm(counts = x, gene_id = "ensembl"),
    function(x) IOBR::count2tpm(countMat = x, idType = "Ensembl", org = "hsa", source = "local")
)
for (f in attempts) {
    try(
        {
            tpm_try <- f(expr_filt)
            if (!is.null(tpm_try)) {
                tpm <- tpm_try
                message("IOBR::count2tpm 成功（使用某种调用方式）。")
                break
            }
        },
        silent = TRUE
    )
}
if (is.null(tpm)) {
    message("IOBR::count2tpm 尝试均失败或返回 NULL，接下来将回退到使用 biomaRt 的手动 TPM 计算。")
} else {
    saveRDS(tpm, file = file.path(outdir, "tpm_ensembl_iobr.rds"))
    write.csv(as.data.frame(tpm), file = file.path(outdir, "tpm_ensembl_iobr.csv"), row.names = TRUE)
    message("已保存 IOBR 生成的 TPM 文件：", file.path(outdir, "tpm_ensembl_iobr.rds"))
}

# -------------------- 8. 若 tpm 为 NULL，使用 biomaRt 手动计算 TPM ----------------
if (is.null(tpm)) {
    message("开始手动 TPM 计算：使用 biomaRt 查询基因长度 ...（可能需要联网）")
    genes_need <- rownames(expr_filt)
    mart <- NULL
    try(
        {
            mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
        },
        silent = TRUE
    )
    if (is.null(mart)) {
        try(
            {
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            },
            silent = TRUE
        )
    }
    if (is.null(mart)) stop("无法通过 biomaRt 连接 Ensembl，请检查网络或提供本地基因长度表。")
    bm <- getBM(attributes = c("ensembl_gene_id", "transcript_length"), filters = "ensembl_gene_id", values = genes_need, mart = mart)
    if (nrow(bm) == 0) stop("biomaRt 对提供的 Ensembl IDs 未返回任何结果，请检查 ID 格式。")
    bm_dt <- as.data.table(bm)
    len_dt <- bm_dt[, .(gene_len = median(transcript_length, na.rm = TRUE)), by = ensembl_gene_id]
    gene_len_vec <- len_dt$gene_len[match(genes_need, len_dt$ensembl_gene_id)]
    valid_len <- !is.na(gene_len_vec) & gene_len_vec > 0
    message("在 Ensembl 中找到长度信息的基因数量：", sum(valid_len), " / ", length(valid_len))
    if (sum(valid_len) < 10) stop("可用基因长度数量过少 (<10)，无法可靠计算 TPM。")
    expr_for_tpm <- expr_filt[valid_len, , drop = FALSE]
    gene_len_kb <- gene_len_vec[valid_len] / 1000
    names(gene_len_kb) <- rownames(expr_for_tpm)
    # 计算 RPK（按行除以基因长度 kb）
    rpk <- sweep(expr_for_tpm, 1, gene_len_kb, "/")
    # 计算 TPM（按列归一化）
    tpm <- sweep(rpk, 2, colSums(rpk, na.rm = TRUE), "/") * 1e6
    message("手动计算 TPM 完成，尺寸（基因 x 样本）：", dim(tpm)[1], " x ", dim(tpm)[2])
    saveRDS(tpm, file = file.path(outdir, "tpm_ensembl_manual.rds"))
    write.csv(as.data.frame(tpm), file = file.path(outdir, "tpm_ensembl_manual.csv"), row.names = TRUE)
}

# -------------------- 9. Ensembl -> HGNC symbol 映射并合并重复 ----------------
message("开始将 Ensembl ID 转换为 HGNC gene symbol ...")
mart2 <- NULL
try(
    {
        mart2 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    },
    silent = TRUE
)
if (is.null(mart2)) {
    try(
        {
            mart2 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        },
        silent = TRUE
    )
}
if (is.null(mart2)) stop("无法连接 Ensembl 以获取基因符号映射，请检查网络。")

ensembl_ids <- rownames(tpm)
map_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl_ids, mart = mart2)
map_dt <- as.data.table(map_df)
map_dt <- map_dt[hgnc_symbol != "" & !is.na(hgnc_symbol)]
symbol_vec <- map_dt$hgnc_symbol[match(ensembl_ids, map_dt$ensembl_gene_id)]
keep_with_symbol <- !is.na(symbol_vec) & symbol_vec != ""
message("映射成功（有 symbol 的行数）：", sum(keep_with_symbol), " / ", length(symbol_vec))
tpm_sym <- tpm[keep_with_symbol, , drop = FALSE]
rownames(tpm_sym) <- symbol_vec[keep_with_symbol]

# 若存在重复 symbol，则按行均值合并
if (any(duplicated(rownames(tpm_sym)))) {
    message("检测到重复的基因符号，开始按行均值合并重复项 ...")
    uniq <- unique(rownames(tpm_sym))
    agg_mat <- matrix(nrow = length(uniq), ncol = ncol(tpm_sym))
    rownames(agg_mat) <- uniq
    colnames(agg_mat) <- colnames(tpm_sym)
    for (i in seq_along(uniq)) {
        rows_i <- which(rownames(tpm_sym) == uniq[i])
        if (length(rows_i) == 1) {
            agg_mat[i, ] <- tpm_sym[rows_i, ]
        } else {
            agg_mat[i, ] <- colMeans(tpm_sym[rows_i, , drop = FALSE], na.rm = TRUE)
        }
    }
    tpm_sym <- as.matrix(agg_mat)
}
message("最终用于 signature 分析的 TPM (symbol) 尺寸：", dim(tpm_sym)[1], " 基因 x ", dim(tpm_sym)[2], " 样本")

# -------------------- 10. 保存最终输出 ----------------
saveRDS(tpm, file = file.path(outdir, "tpm_ensembl_final.rds"))
write.csv(as.data.frame(tpm), file = file.path(outdir, "tpm_ensembl_final.csv"), row.names = TRUE)
saveRDS(tpm_sym, file = file.path(outdir, "tpm_symbol_final.rds"))
write.csv(as.data.frame(tpm_sym), file = file.path(outdir, "tpm_symbol_final.csv"), row.names = TRUE)
saveRDS(pheno_dt, file = file.path(outdir, "pheno_aligned.rds"))
write.csv(pheno_dt, file = file.path(outdir, "pheno_aligned.csv"), row.names = TRUE)
message("已将所有主要输出保存到目录：", normalizePath(outdir))

# -------------------- 11. 快速诊断 / 摘要 ----------------
cat("\n=== 预处理摘要 ===\n")
cat("已保存：过滤后的 counts ->", file.path(outdir, "counts_filtered.rds"), "\n")
cat("已保存：TPM (ensembl) 文件 ->", list.files(outdir, pattern = "tpm_ensembl", full.names = TRUE), "\n")
cat("已保存：TPM (symbol) 文件 ->", file.path(outdir, "tpm_symbol_final.rds"), "\n")
cat("已保存：对齐后的 pheno ->", file.path(outdir, "pheno_aligned.rds"), "\n")
cat("样本数量（过滤后）:", ncol(expr_filt), "\n")
cat("基因数量（最终 symbol TPM）:", nrow(tpm_sym), "\n")
cat("输出目录下的文件列表：\n")
print(list.files(outdir))

# 脚本结束
message("预处理流程完成。你现在可以使用 tpm_symbol_final.rds (或 csv) 作为 IOBR 下游分析的输入。")
########################################################################
# 结束
########################################################################
