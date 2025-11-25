# 文件2: iobr_analysis.R
# IOBR免疫浸润分析和signature评分

library(IOBR)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

cat("开始IOBR分析...\n")

# 读取预处理数据
gse_id <- "GSE289743"
out_dir <- "data"

tpm_matrix <- readRDS(file.path(out_dir, paste0(gse_id, "_tpm_matrix.rds")))
clinical_data <- readRDS(file.path(out_dir, paste0(gse_id, "_clinical_data.rds")))

cat("TPM矩阵维度:", dim(tpm_matrix), "\n")
cat("临床数据样本数:", nrow(clinical_data), "\n")

# 确保样本匹配
common_samples <- intersect(colnames(tpm_matrix), rownames(clinical_data))
tpm_matrix <- tpm_matrix[, common_samples]
clinical_data <- clinical_data[common_samples, ]

cat("共同样本数量:", length(common_samples), "\n")

# 1. 免疫细胞浸润分析 - 使用多种方法
cat("进行免疫细胞浸润分析...\n")

# 方法1: CIBERSORT
tryCatch(
    {
        cibersort_result <- deconvo_tme(eset = tpm_matrix, method = "cibersort", arrays = FALSE)
        cat("CIBERSORT分析完成\n")
    },
    error = function(e) {
        cat("CIBERSORT分析失败:", e$message, "\n")
        cibersort_result <- NULL
    }
)

# 方法2: EPIC
tryCatch(
    {
        epic_result <- deconvo_tme(eset = tpm_matrix, method = "epic", arrays = FALSE)
        cat("EPIC分析完成\n")
    },
    error = function(e) {
        cat("EPIC分析失败:", e$message, "\n")
        epic_result <- NULL
    }
)

# 方法3: xCell
tryCatch(
    {
        xcell_result <- deconvo_tme(eset = tpm_matrix, method = "xcell", arrays = FALSE)
        cat("xCell分析完成\n")
    },
    error = function(e) {
        cat("xCell分析失败:", e$message, "\n")
        xcell_result <- NULL
    }
)

# 方法4: MCPcounter
tryCatch(
    {
        mcp_result <- deconvo_tme(eset = tpm_matrix, method = "mcpcounter", arrays = FALSE)
        cat("MCPcounter分析完成\n")
    },
    error = function(e) {
        cat("MCPcounter分析失败:", e$message, "\n")
        mcp_result <- NULL
    }
)

# 2. 计算signature评分
cat("计算signature评分...\n")

# 获取IOBR内置signature列表
signature_list <- signature_collection
cat("可用signature数量:", length(signature_list), "\n")

# 选择一些关键的免疫相关signature进行计算
selected_signatures <- c(
    "CD8_T_cells",
    "T_cell_inflamed_GEP",
    "Immunogenicity",
    "Cytolytic_activity",
    "NK_cells",
    "B_cells",
    "Macrophages",
    "Dendritic_cells"
)

# 计算signature评分
signature_scores <- calculate_sig_score(
    pdata = NULL,
    eset = tpm_matrix,
    signature = selected_signatures,
    method = "pca" # 使用PCA方法
)

cat("signature评分计算完成\n")
cat("signature评分维度:", dim(signature_scores), "\n")

# 3. 添加自定义signature
cat("添加自定义signature...\n")

# 自定义signature 1: HALLMARK_APOPTOSIS (细胞凋亡)
apoptosis_genes <- c(
    "BAX", "BAK1", "BCL2", "BCL2L1", "CASP3", "CASP8", "CASP9",
    "FAS", "FADD", "TNF", "TNFSF10", "TP53", "PIDD", "APAF1"
)

# 自定义signature 2: HALLMARK_ANGIOGENESIS (血管生成)
angiogenesis_genes <- c(
    "VEGFA", "VEGFC", "FLT1", "KDR", "FLT4", "PGF", "ANGPT1",
    "TEK", "PDGFA", "PDGFB", "FGFR1", "FGFR2", "EFNA1", "EPHB4"
)

# 自定义signature 3: HALLMARK_INTERFERON_GAMMA_RESPONSE (干扰素γ反应)
ifn_gamma_genes <- c(
    "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "CIITA",
    "HLA-DRA", "HLA-DRB1", "PSMB8", "PSMB9", "TAP1", "TAP2"
)

# 计算自定义signature评分
calculate_custom_signature <- function(expr_matrix, gene_list, method = "pca") {
    # 只保留在表达矩阵中存在的基因
    available_genes <- gene_list[gene_list %in% rownames(expr_matrix)]

    if (length(available_genes) < 3) {
        cat("警告: 可用基因数量不足(", length(available_genes), ")\n")
        return(rep(NA, ncol(expr_matrix)))
    }

    if (method == "pca") {
        # PCA方法
        sig_expr <- expr_matrix[available_genes, , drop = FALSE]
        # 移除方差为0的基因
        gene_vars <- apply(sig_expr, 1, var)
        sig_expr <- sig_expr[gene_vars > 0, , drop = FALSE]

        if (nrow(sig_expr) < 2) {
            return(rep(NA, ncol(expr_matrix)))
        }

        pca_result <- prcomp(t(sig_expr), scale. = TRUE)
        scores <- pca_result$x[, 1]
        return(scores)
    } else if (method == "mean") {
        # 简单均值方法
        sig_expr <- expr_matrix[available_genes, ]
        scores <- colMeans(sig_expr)
        return(scores)
    }
}

# 计算自定义signature
custom_scores <- data.frame(
    HALLMARK_APOPTOSIS = calculate_custom_signature(tpm_matrix, apoptosis_genes, "pca"),
    HALLMARK_ANGIOGENESIS = calculate_custom_signature(tpm_matrix, angiogenesis_genes, "pca"),
    HALLMARK_INTERFERON_GAMMA_RESPONSE = calculate_custom_signature(tpm_matrix, ifn_gamma_genes, "pca")
)

rownames(custom_scores) <- colnames(tpm_matrix)

# 4. 合并所有结果
cat("合并分析结果...\n")

# 创建综合结果数据框
full_results <- data.frame(
    Sample = colnames(tpm_matrix)
)

# 添加临床信息
full_results <- cbind(full_results, clinical_data[full_results$Sample, ])

# 添加免疫浸润结果
if (!is.null(cibersort_result)) {
    full_results <- cbind(full_results, cibersort_result$score[full_results$Sample, ])
}

# 添加signature评分
full_results <- cbind(full_results, signature_scores[full_results$Sample, ])

# 添加自定义signature
full_results <- cbind(full_results, custom_scores[full_results$Sample, ])

rownames(full_results) <- full_results$Sample

# 5. 保存结果
cat("保存分析结果...\n")

# 保存完整结果
saveRDS(full_results, file.path(out_dir, paste0(gse_id, "_full_iobr_results.rds")))
write.csv(full_results, file.path(out_dir, paste0(gse_id, "_full_iobr_results.csv")))

# 保存各个单独结果
if (!is.null(cibersort_result)) {
    saveRDS(cibersort_result, file.path(out_dir, paste0(gse_id, "_cibersort.rds")))
}
if (!is.null(epic_result)) {
    saveRDS(epic_result, file.path(out_dir, paste0(gse_id, "_epic.rds")))
}
if (!is.null(xcell_result)) {
    saveRDS(xcell_result, file.path(out_dir, paste0(gse_id, "_xcell.rds")))
}
if (!is.null(mcp_result)) {
    saveRDS(mcp_result, file.path(out_dir, paste0(gse_id, "_mcpcounter.rds")))
}

saveRDS(signature_scores, file.path(out_dir, paste0(gse_id, "_signature_scores.rds")))
saveRDS(custom_scores, file.path(out_dir, paste0(gse_id, "_custom_signatures.rds")))

cat("IOBR分析完成!\n")
cat("最终结果维度:", dim(full_results), "\n")
cat("包含的特征:", colnames(full_results), "\n")