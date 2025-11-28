########################################################################
# visualization_iobr_only.R
# 基于 analysis_iobr_only.R 的输出用 IOBR 结果绘制论文级图（热图、箱线/小提琴、相关性图、PCA）
# 完全使用 R 的绘图包（ComplexHeatmap / ggplot2），但输入数据均来源 IOBR 输出。
########################################################################

outdir <- "results_preprocess"
dir.create(outdir, showWarnings = FALSE)

# 载包
required <- c("ComplexHeatmap", "circlize", "RColorBrewer", "pheatmap", "ggplot2", "ggpubr", "data.table", "dplyr")
for (p in required) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(pheatmap)
    library(ggplot2)
    library(ggpubr)
    library(data.table)
    library(dplyr)
})

# 1. 加载 IOBR 生成的对象
sig_comb_file <- file.path(outdir, "iobr_sig_combined_samples_x_signatures.rds")
sig_raw_file <- file.path(outdir, "iobr_sig_scores_raw.rds")
deconv_file <- file.path(outdir, "iobr_deconv_results_list.rds")
pheno_file <- file.path(outdir, "pheno_aligned.rds")

if (!file.exists(sig_comb_file) || !file.exists(pheno_file)) stop("缺少 signature 合并文件或 pheno，请先运行 analysis_iobr_only.R")

sig_comb <- readRDS(sig_comb_file) # samples x signatures
pheno_dt <- readRDS(pheno_file)
deconv_list <- if (file.exists(deconv_file)) readRDS(deconv_file) else NULL

cat("[INFO] 读取 sig_comb dims:", dim(sig_comb), "\n")

# 2. 构建样本注释（annotation_col）并生成 annotation_colors
grp_candidates <- intersect(c("response:ch1", "response", "Response", "group"), colnames(pheno_dt))
grp_col <- ifelse(length(grp_candidates) >= 1, grp_candidates[1], NA)
annotation_col <- data.frame(row.names = rownames(pheno_dt))
if (!is.na(grp_col)) {
    annotation_col[[grp_col]] <- as.character(pheno_dt[[grp_col]])
    annotation_col[[grp_col]][is.na(annotation_col[[grp_col]])] <- "NA"
    annotation_col[[grp_col]] <- factor(annotation_col[[grp_col]])
}
# 生成 colors
annotation_colors <- list()
if (ncol(annotation_col) > 0) {
    for (coln in colnames(annotation_col)) {
        levs <- levels(annotation_col[[coln]])
        k <- length(levs)
        if (k <= 12) pal <- brewer.pal(max(3, k), "Set3") else pal <- colorRampPalette(brewer.pal(8, "Set3"))(k)
        names(pal) <- levs
        annotation_colors[[coln]] <- pal
    }
}

# 3. Heatmap：取方差最大的 top 40 signatures
sig_sd <- apply(sig_comb, 2, sd, na.rm = TRUE)
topn <- min(40, length(sig_sd))
top_sigs <- names(sort(sig_sd, decreasing = TRUE))[1:topn]
mat <- t(scale(t(as.matrix(sig_comb[, top_sigs, drop = FALSE])))) # rows signatures, cols samples (zscore)
# ensure sample names match annotation
if (ncol(mat) != nrow(annotation_col)) {
    common <- intersect(colnames(mat), rownames(annotation_col))
    mat <- mat[, common, drop = FALSE]
    annotation_col <- annotation_col[common, , drop = FALSE]
}
png(file.path(outdir, "iobr_heatmap_top40_signatures.png"), width = 2200, height = 1400, res = 200)
ht <- Heatmap(mat,
    name = "zscore", top_annotation = HeatmapAnnotation(df = annotation_col, col = annotation_colors),
    show_row_names = TRUE, show_column_names = FALSE, cluster_rows = TRUE, cluster_columns = TRUE
)
draw(ht)
dev.off()
cat("[INFO] 已输出 heatmap: iobr_heatmap_top40_signatures.png\n")

# 4. 箱线/小提琴图（示例：top sig 在 response 分组的差异）
if (!is.na(grp_col)) {
    sig1 <- top_sigs[1]
    df <- data.frame(sample = rownames(sig_comb), group = annotation_col[[grp_col]], score = sig_comb[, sig1])
    p <- ggplot(df, aes(x = group, y = score, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.12, outlier.size = 0.6) +
        geom_jitter(width = 0.15, size = 0.7) +
        stat_compare_means(aes(label = ..p.format..), method = "wilcox.test") +
        theme_minimal() +
        ggtitle(paste(sig1, "by", grp_col))
    ggsave(file.path(outdir, paste0("iobr_violin_box_", sig1, ".png")), p, width = 5, height = 4, dpi = 300)
    cat("[INFO] 已输出箱线/小提琴图: ", paste0("iobr_violin_box_", sig1, ".png"), "\n")
}

# 5. signature vs cell type 相关性（以 CIBERSORT 为例，若存在）
if (!is.null(deconv_list) && "cibersort" %in% names(deconv_list)) {
    cell_mat <- deconv_list[["cibersort"]]
    # 确保行样本匹配
    common <- intersect(rownames(cell_mat), rownames(sig_comb))
    if (length(common) > 2) {
        cell_mat2 <- cell_mat[common, , drop = FALSE]
        sig2 <- sig_comb[common, , drop = FALSE]
        # 选第一个 signature 与第一个 cell type 做散点图
        cell_type <- colnames(cell_mat2)[1]
        sig_name <- top_sigs[1]
        df_sc <- data.frame(
            sample = common, cell = as.numeric(cell_mat2[, cell_type]), sig = as.numeric(sig2[, sig_name]),
            group = annotation_col[common, grp_col, drop = TRUE]
        )
        p_sc <- ggplot(df_sc, aes(x = sig, y = cell, color = group)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", se = TRUE) +
            stat_cor(method = "spearman") +
            theme_minimal() +
            labs(x = sig_name, y = paste0(cell_type, " fraction")) +
            ggtitle(paste(sig_name, "vs", cell_type))
        ggsave(file.path(outdir, paste0("iobr_scatter_", sig_name, "_vs_", cell_type, ".png")), p_sc, width = 6, height = 5, dpi = 300)
        cat("[INFO] 已输出 signature vs celltype 散点图。\n")
    } else {
        cat("[WARN] CIBERSORT 结果或样本数不足，跳过相关性图。\n")
    }
}

# 6. PCA（基于 sig_comb）
if (ncol(sig_comb) >= 3) {
    pca <- prcomp(sig_comb, center = TRUE, scale. = TRUE)
    pc_df <- data.frame(
        PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = rownames(pca$x),
        group = if (!is.na(grp_col)) annotation_col[rownames(pca$x), grp_col] else NA
    )
    p_pca <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 2) +
        theme_minimal() +
        ggtitle("PCA on IOBR signatures")
    ggsave(file.path(outdir, "iobr_pca_signatures.png"), p_pca, width = 6, height = 5, dpi = 300)
    cat("[INFO] 已输出 PCA 图: iobr_pca_signatures.png\n")
}

# 7. 保存用于审阅的表格（sig_comb, deconv if available）
write.csv(sig_comb, file.path(outdir, "iobr_sig_combined_for_review.csv"), row.names = TRUE)
if (!is.null(deconv_list) && "cibersort" %in% names(deconv_list)) write.csv(deconv_list[["cibersort"]], file.path(outdir, "iobr_cibersort_for_review.csv"), row.names = TRUE)

cat("[DONE] 可视化（IOBR-only）完成，图与表保存在：", normalizePath(outdir), "\n")
########################################################################
# END visualization_iobr_only.R
########################################################################
