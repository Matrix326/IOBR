# 修正版 analysis_iobr.R - 检查可用signature
# IOBR免疫浸润分析和signature评分 - 修正signature名称

library(IOBR)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

cat("开始IOBR分析...\n")

# 读取预处理数据
gse_id <- "GSE289743"
out_dir <- "data"

# 确保 signature_collection 已定义
if (!exists("signature_collection")) {
    stop("Error: 'signature_collection' is not defined. Please load or define it before running the script.")
}

tpm_matrix <- readRDS(file.path(out_dir, paste0(gse_id, "_tpm_matrix.rds")))
clinical_data <- readRDS(file.path(out_dir, paste0(gse_id, "_clinical_data.rds")))

cat("TPM矩阵维度:", dim(tpm_matrix), "\n")
cat("临床数据样本数:", nrow(clinical_data), "\n")

# 确保样本匹配
common_samples <- intersect(colnames(tpm_matrix), rownames(clinical_data))
tpm_matrix <- tpm_matrix[, common_samples]
clinical_data <- clinical_data[common_samples, ]

cat("共同样本数量:", length(common_samples), "\n")

# 1. 检查可用的signature
cat("检查可用signature...\n")
available_signatures <- names(signature_collection)
cat("总signature数量:", length(available_signatures), "\n")

# 打印前20个signature名称
cat("前20个可用signature:\n")
print(head(available_signatures, 20))

# 查找包含特定关键词的signature
immune_related <- available_signatures[grepl("CD8|T_cell|immune|cytolytic|NK|B_cell|macrophage|dendritic",
    available_signatures,
    ignore.case = TRUE
)]
cat("免疫相关signature:\n")
print(immune_related)

# 2. 免疫细胞浸润分析
cat("进行免疫细胞浸润分析...\n")

# 方法1: CIBERSORT
cat(">>> Running CIBERSORT\n")
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
cat(">>> Running EPIC\n")
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

# 方法3: MCPcounter
cat(">>> Running MCP-counter\n")
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

# 3. 计算signature评分 - 使用找到的正确signature名称
cat("计算signature评分...\n")

# 选择一些确定存在的signature
# 根据检查结果选择可用的signature
safe_signatures <- c()
if (length(immune_related) > 0) {
    safe_signatures <- head(immune_related, 5) # 使用前5个免疫相关signature
} else {
    # 如果找不到免疫相关signature，使用一些通用signature
    safe_signatures <- head(available_signatures, 5)
}

cat("将计算以下signature:\n")
print(safe_signatures)

# 稳健的PCA计算方法
calculate_signature_pca <- function(eset, signature_name) {
    tryCatch(
        {
            # 检查signature是否存在
            if (!signature_name %in% names(signature_collection)) {
                cat("警告: signature", signature_name, "不存在\n")
                return(NULL)
            }

            sig_genes <- signature_collection[[signature_name]]
            available_genes <- sig_genes[sig_genes %in% rownames(eset)]

            if (length(available_genes) < 3) {
                cat("警告: signature", signature_name, "可用基因不足(", length(available_genes), ")\n")
                return(NULL)
            }

            # 提取signature基因表达矩阵
            sig_expr <- eset[available_genes, , drop = FALSE]

            # 检查并处理零方差基因
            gene_vars <- apply(sig_expr, 1, var)
            zero_var_genes <- gene_vars == 0

            if (sum(zero_var_genes) > 0) {
                cat("移除", sum(zero_var_genes), "个零方差基因\n")
                sig_expr <- sig_expr[!zero_var_genes, , drop = FALSE]
            }

            # 再次检查剩余基因数量
            if (nrow(sig_expr) < 2) {
                cat("警告: 剩余基因数量不足进行PCA\n")
                # 使用均值作为备选
                scores <- colMeans(sig_expr)
                return(scores)
            }

            # 执行PCA
            pca_result <- prcomp(t(sig_expr), scale. = TRUE, center = TRUE)

            # 使用第一主成分作为signature score
            scores <- pca_result$x[, 1]

            # 解释符号：确保与基因表达正相关
            loading_cor <- cor(pca_result$rotation[, 1], rowMeans(sig_expr))
            if (loading_cor < 0) {
                scores <- -scores
            }

            return(scores)
        },
        error = function(e) {
            cat("计算signature", signature_name, "时PCA出错:", e$message, "\n")
            # 出错时使用均值方法
            sig_genes <- signature_collection[[signature_name]]
            available_genes <- sig_genes[sig_genes %in% rownames(eset)]
            if (length(available_genes) > 0) {
                sig_expr <- eset[available_genes, , drop = FALSE]
                return(colMeans(sig_expr))
            } else {
                return(rep(NA, ncol(eset)))
            }
        }
    )
}

# 计算每个signature
signature_results <- list()
for (sig in safe_signatures) {
    cat("计算signature:", sig, "\n")
    scores <- calculate_signature_pca(tpm_matrix, sig)
    if (!is.null(scores) && !all(is.na(scores))) {
        signature_results[[sig]] <- scores
        cat("  - 成功计算", sig, "\n")
    } else {
        cat("  - 失败计算", sig, "\n")
    }
}

# 转换为数据框
if (length(signature_results) > 0) {
    signature_scores <- as.data.frame(do.call(cbind, signature_results))
    rownames(signature_scores) <- colnames(tpm_matrix)
    cat("signature评分计算完成，维度:", dim(signature_scores), "\n")
} else {
    cat("警告: 没有成功计算任何signature评分\n")
    signature_scores <- data.frame(matrix(nrow = ncol(tpm_matrix), ncol = 0))
    rownames(signature_scores) <- colnames(tpm_matrix)
}

# 4. 添加自定义signature - 使用PCA方法
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

# 稳健的自定义signature PCA计算
calculate_custom_signature_pca <- function(expr_matrix, gene_list) {
    tryCatch(
        {
            # 只保留在表达矩阵中存在的基因
            available_genes <- gene_list[gene_list %in% rownames(expr_matrix)]

            if (length(available_genes) < 3) {
                cat("警告: 可用基因数量不足(", length(available_genes), ")\n")
                return(rep(NA, ncol(expr_matrix)))
            }

            # 提取基因表达
            sig_expr <- expr_matrix[available_genes, , drop = FALSE]

            # 检查并处理零方差基因
            gene_vars <- apply(sig_expr, 1, var)
            zero_var_genes <- gene_vars == 0

            if (sum(zero_var_genes) > 0) {
                cat("移除", sum(zero_var_genes), "个零方差基因\n")
                sig_expr <- sig_expr[!zero_var_genes, , drop = FALSE]
            }

            # 再次检查剩余基因数量
            if (nrow(sig_expr) < 2) {
                cat("警告: 剩余基因数量不足进行PCA\n")
                # 使用均值作为备选
                return(colMeans(sig_expr))
            }

            # 执行PCA
            pca_result <- prcomp(t(sig_expr), scale. = TRUE, center = TRUE)

            # 使用第一主成分作为signature score
            scores <- pca_result$x[, 1]

            # 解释符号：确保与基因表达正相关
            loading_cor <- cor(pca_result$rotation[, 1], rowMeans(sig_expr))
            if (loading_cor < 0) {
                scores <- -scores
            }

            return(scores)
        },
        error = function(e) {
            cat("自定义signature PCA计算出错:", e$message, "\n")
            # 出错时使用均值方法
            available_genes <- gene_list[gene_list %in% rownames(expr_matrix)]
            if (length(available_genes) > 0) {
                sig_expr <- expr_matrix[available_genes, , drop = FALSE]
                return(colMeans(sig_expr))
            } else {
                return(rep(NA, ncol(expr_matrix)))
            }
        }
    )
}

# 计算自定义signature
cat("计算HALLMARK_APOPTOSIS...\n")
apoptosis_scores <- calculate_custom_signature_pca(tpm_matrix, apoptosis_genes)

cat("计算HALLMARK_ANGIOGENESIS...\n")
angiogenesis_scores <- calculate_custom_signature_pca(tpm_matrix, angiogenesis_genes)

cat("计算HALLMARK_INTERFERON_GAMMA_RESPONSE...\n")
ifn_gamma_scores <- calculate_custom_signature_pca(tpm_matrix, ifn_gamma_genes)

custom_scores <- data.frame(
    HALLMARK_APOPTOSIS = apoptosis_scores,
    HALLMARK_ANGIOGENESIS = angiogenesis_scores,
    HALLMARK_INTERFERON_GAMMA_RESPONSE = ifn_gamma_scores
)

rownames(custom_scores) <- colnames(tpm_matrix)
cat("自定义signature计算完成，维度:", dim(custom_scores), "\n")

# 5. 合并所有结果 - 修复score列问题
cat("合并分析结果...\n")

# 创建基础结果数据框
full_results <- data.frame(
    Sample = colnames(tpm_matrix),
    stringsAsFactors = FALSE
)
rownames(full_results) <- full_results$Sample

# 添加临床信息
clinical_subset <- clinical_data[full_results$Sample, , drop = FALSE]
full_results <- cbind(full_results, clinical_subset)

# 修复：正确提取免疫浸润结果
# CIBERSORT
if (!is.null(cibersort_result)) {
    # 检查结果结构
    if (is.list(cibersort_result) && "score" %in% names(cibersort_result)) {
        cibersort_data <- cibersort_result$score[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, cibersort_data)
        cat("已添加CIBERSORT结果\n")
    } else if (is.data.frame(cibersort_result)) {
        # 如果直接是数据框
        cibersort_data <- cibersort_result[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, cibersort_data)
        cat("已添加CIBERSORT结果（直接数据框）\n")
    } else {
        cat("CIBERSORT结果结构未知\n")
    }
}

# EPIC
if (!is.null(epic_result)) {
    if (is.list(epic_result) && "score" %in% names(epic_result)) {
        epic_data <- epic_result$score[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, epic_data)
        cat("已添加EPIC结果\n")
    } else if (is.data.frame(epic_result)) {
        epic_data <- epic_result[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, epic_data)
        cat("已添加EPIC结果（直接数据框）\n")
    } else {
        cat("EPIC结果结构未知\n")
    }
}

# MCPcounter
if (!is.null(mcp_result)) {
    if (is.list(mcp_result) && "score" %in% names(mcp_result)) {
        mcp_data <- mcp_result$score[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, mcp_data)
        cat("已添加MCPcounter结果\n")
    } else if (is.data.frame(mcp_result)) {
        mcp_data <- mcp_result[full_results$Sample, , drop = FALSE]
        full_results <- cbind(full_results, mcp_data)
        cat("已添加MCPcounter结果（直接数据框）\n")
    } else {
        cat("MCPcounter结果结构未知\n")
    }
}

# 添加signature评分
if (ncol(signature_scores) > 0) {
    sig_data <- signature_scores[full_results$Sample, , drop = FALSE]
    full_results <- cbind(full_results, sig_data)
    cat("已添加signature评分\n")
}

# 添加自定义signature
custom_data <- custom_scores[full_results$Sample, , drop = FALSE]
full_results <- cbind(full_results, custom_data)
cat("已添加自定义signature\n")

# 检查最终结果
cat("最终结果维度:", dim(full_results), "\n")
cat("包含的特征数量:", ncol(full_results), "\n")

# 6. 保存结果
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
if (!is.null(mcp_result)) {
    saveRDS(mcp_result, file.path(out_dir, paste0(gse_id, "_mcpcounter.rds")))
}

saveRDS(signature_scores, file.path(out_dir, paste0(gse_id, "_signature_scores.rds")))
saveRDS(custom_scores, file.path(out_dir, paste0(gse_id, "_custom_signatures.rds")))

# 生成结果摘要
cat("\n=== 分析完成摘要 ===\n")
cat("样本数量:", nrow(full_results), "\n")
cat("特征数量:", ncol(full_results), "\n")
cat("成功运行的方法:\n")
if (!is.null(cibersort_result)) cat("- CIBERSORT\n")
if (!is.null(epic_result)) cat("- EPIC\n")
if (!is.null(mcp_result)) cat("- MCPcounter\n")
cat("使用PCA计算的IOBR signature数量:", ncol(signature_scores), "\n")
cat("使用PCA计算的自定义signature数量:", ncol(custom_scores), "\n")
cat("最终数据框包含的列:", colnames(full_results), "\n")

cat("IOBR分析完成!\n")
