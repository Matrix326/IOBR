########################################################################
# analysis_iobr_only.R
# 完全基于 IOBR 内置函数实现（不使用 GSVA 回退）
# 作用：读取预处理产物 -> 用 IOBR 计算 signature scores (IOBR 内置集 + 自定义)
#       -> 使用 IOBR deconvo_tme(Builtin methods, e.g., cibersort) 进行免疫浸润估计
#       -> 保存所有结果到 results_preprocess/
########################################################################

# 输出目录（与预处理一致）

library(GSVA)
outdir <- "results_preprocess"
dir.create(outdir, showWarnings = FALSE)

# 载入包（IOBR 必须已安装）
if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("请先安装 IOBR 包（remotes::install_github('IOBR/IOBR') 或 BiocManager::install）。")
}
suppressPackageStartupMessages({
    library(IOBR)
    library(data.table)
})

# 1. 读取预处理输出
tpm_file <- file.path(outdir, "tpm_symbol_final.rds")
pheno_file <- file.path(outdir, "pheno_aligned.rds")
if (!file.exists(tpm_file) || !file.exists(pheno_file)) stop("请先运行预处理脚本，确保存在 tpm_symbol_final.rds 与 pheno_aligned.rds")

tpm <- readRDS(tpm_file) # 行 = gene symbol, 列 = sample
pheno_dt <- readRDS(pheno_file) # 行 = sample (rownames = sample IDs)

cat("[INFO] 读取 TPM 与 pheno，TPM dims:", dim(tpm), "; pheno dims:", dim(pheno_dt), "\n")

# 2. 样本对齐（必须：IOBR 后续函数依赖样本对应）
if (!all(colnames(tpm) == rownames(pheno_dt))) {
    common <- intersect(colnames(tpm), rownames(pheno_dt))
    if (length(common) == 0) stop("TPM 与 pheno 无公共样本名，请检查样本命名。")
    # 采用交集样本并按 pheno 顺序重排
    tpm <- tpm[, common, drop = FALSE]
    pheno_dt <- pheno_dt[common, , drop = FALSE]
    cat("[WARN] 仅保留 TPM 与 pheno 的交集样本，样本数:", length(common), "\n")
}

# 3. 过滤低表达基因（基于 TPM）以加快计算 - IOBR 官方示例通常会先过滤
keep_gene <- rowSums(tpm >= 1) >= 2 # 基因在至少 2 个样本中 TPM>=1
cat("[INFO] 过滤前基因数:", nrow(tpm), "，过滤后保留:", sum(keep_gene), "\n")
tpm_filt <- tpm[keep_gene, , drop = FALSE]

# 4. 使用 IOBR 内置 signature 打分（ssGSEA）
sig_rds <- file.path(outdir, "iobr_sig_scores_ssgsea.rds")
if (file.exists(sig_rds)) {
    sig_scores <- readRDS(sig_rds)
    cat("[INFO] 已存在 IOBR signature 得分，载入：", sig_rds, "\n")
} else {
    cat("[INFO] 使用 IOBR::calculate_sig_score 计算内置 signature（method = ssGSEA）...\n")
    # IOBR::calculate_sig_score 默认会使用 IOBR 内置 signature 集（如果不传 signature 参数）
    sig_scores <- IOBR::calculate_sig_score(eset = tpm_filt, method = "ssgsea")
    # IOBR 有不同版本的返回格式：确保保存原始返回（一般为 samples x signatures）
    saveRDS(sig_scores, sig_rds)
    write.csv(as.data.frame(sig_scores), file.path(outdir, "iobr_sig_scores_ssgsea.csv"), row.names = TRUE)
    cat("[INFO] signature 打分完成并已保存。\n")
}

# 5. 计算内置 deconvolution（使用 IOBR::deconvo_tme）
# 默认尝试方法列表（你可以在此加入更多 IOBR 支持的方法）
methods_to_try <- c("cibersort", "mcpcounter", "timer", "epic", "xcell") # IOBR 支持的方法集合，实际可用依赖本地 IOBR 配置
deconv_list <- list()
for (m in methods_to_try) {
    out_rds <- file.path(outdir, paste0("iobr_deconv_", m, ".rds"))
    out_csv <- file.path(outdir, paste0("iobr_deconv_", m, ".csv"))
    # 若文件已存在则跳过计算
    if (file.exists(out_rds)) {
        cat("[INFO] 已存在 deconv 结果，载入：", out_rds, "\n")
        deconv_list[[m]] <- readRDS(out_rds)
        next
    }
    # 尝试调用 IOBR::deconvo_tme（捕获错误）
    cat("[INFO] 尝试 deconvo_tme(method='", m, "') ...\n", sep = "")
    res <- tryCatch(
        {
            IOBR::deconvo_tme(eset = tpm_filt, method = m)
        },
        error = function(e) {
            cat("[WARN] deconvo_tme 方法", m, "失败：", e$message, "\n")
            return(NULL)
        }
    )
    if (!is.null(res)) {
        deconv_list[[m]] <- res
        saveRDS(res, out_rds)
        # 若返回是数据框或矩阵则保存csv，方便后续可视化
        if (is.data.frame(res) || is.matrix(res)) write.csv(as.data.frame(res), out_csv, row.names = TRUE)
        cat("[INFO] deconv (", m, ") 完成并保存。\n", sep = "")
    } else {
        cat("[INFO] 跳过方法：", m, "\n")
    }
}

# 6. 计算自定义 signature（仅用 IOBR 内置接口）
# 例子：两个自定义 signature（你可以替换为 MSigDB hallmark 或任意基因集）
custom_sig <- list(
    MY_CYTOTOXIC = c("CD8A", "GZMB", "PRF1"),
    MY_IFNG = c("IFNG", "CXCL9", "CXCL10")
)
# IOBR 提供 format_signatures() 来处理自定义签名（接受 data.frame 或 list）
# 先把 custom_sig 转为 IOBR 期望格式（data.frame）
sig_df <- format_signatures(as.data.frame(sapply(custom_sig, function(x) {
    length(x) <- max(sapply(custom_sig, length))
    x
})))

custom_rds <- file.path(outdir, "iobr_custom_sig_scores.rds")
if (file.exists(custom_rds)) {
    custom_scores <- readRDS(custom_rds)
    cat("[INFO] 自定义 signature 已存在，载入：", custom_rds, "\n")
} else {
    cat("[INFO] 使用 IOBR::calculate_sig_score 对自定义 signature 打分...\n")
    custom_scores <- IOBR::calculate_sig_score(eset = tpm_filt, signature = sig_df, method = "ssgsea")
    saveRDS(custom_scores, custom_rds)
    write.csv(as.data.frame(custom_scores), file.path(outdir, "iobr_custom_sig_scores.csv"), row.names = TRUE)
    cat("[INFO] 自定义 signature 打分已保存。\n")
}

# 7. 合并 built-in 与 custom（若需要）
# sig_scores 与 custom_scores 的行/列方向可能依 IOBR 版本不同：确保得到 samples x signatures 的数据.frames
to_df <- function(x) {
    # 如果是 matrix 或 data.frame 就转 data.frame
    if (is.matrix(x) || is.data.frame(x)) {
        return(as.data.frame(x))
    }
    # 否则尝试转置
    return(as.data.frame(t(as.matrix(x))))
}
sig_df_all <- to_df(sig_scores)
custom_df <- to_df(custom_scores)
# 若样本是列名则转置：以样本名相等为准
if (!all(rownames(sig_df_all) %in% rownames(pheno_dt)) && all(colnames(sig_df_all) %in% rownames(pheno_dt))) {
    sig_df_all <- as.data.frame(t(sig_df_all))
}
if (!all(rownames(custom_df) %in% rownames(pheno_dt)) && all(colnames(custom_df) %in% rownames(pheno_dt))) {
    custom_df <- as.data.frame(t(custom_df))
}
# 仅合并共同样本
common_samples <- intersect(rownames(sig_df_all), rownames(custom_df))
if (length(common_samples) == 0) {
    # 若没有共同样本，尝试按 pheno 行合并
    common_samples <- intersect(rownames(sig_df_all), rownames(pheno_dt))
}
sig_df_all <- sig_df_all[common_samples, , drop = FALSE]
custom_df <- custom_df[common_samples, , drop = FALSE]
sig_combined <- cbind(sig_df_all, custom_df)
saveRDS(sig_combined, file.path(outdir, "iobr_sig_combined_samples_x_signatures.rds"))
write.csv(as.data.frame(sig_combined), file.path(outdir, "iobr_sig_combined_samples_x_signatures.csv"), row.names = TRUE)
cat("[INFO] 已保存合并后的 signature 矩阵（samples x signatures）。\n")

# 8. 计算 signature 相关性矩阵并保存（Spearman）
cor_mat <- tryCatch(
    {
        cor(sig_combined, method = "spearman", use = "pairwise.complete.obs")
    },
    error = function(e) {
        cat("[WARN] 计算相关矩阵失败：", e$message, "\n")
        return(NULL)
    }
)
if (!is.null(cor_mat)) write.csv(cor_mat, file.path(outdir, "iobr_signature_corr_spearman.csv"), row.names = TRUE)

# 9. 保存清单与结束
saveRDS(deconv_list, file.path(outdir, "iobr_deconv_results_list.rds"))
saveRDS(sig_scores, file.path(outdir, "iobr_sig_scores_raw.rds"))
cat("[DONE] 分析（IOBR-only）完成。所有输出保存在：", normalizePath(outdir), "\n")
########################################################################
# END analysis_iobr_only.R
########################################################################
