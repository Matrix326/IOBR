# 文件3: clustering_visualization.R
# 聚类分析和可视化

library(IOBR)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)

# 读取IOBR分析结果
gse_id <- "GSE289743"
out_dir <- "data"

full_results <- readRDS(file.path(out_dir, paste0(gse_id, "_iobr_results.rds")))
pheno <- readRDS(file.path(out_dir, paste0(gse_id, "_clean_pheno.rds")))

# 1. 选择用于聚类的特征
# 选择免疫细胞比例和关键signature
clustering_features <- full_results[, c(
    # 免疫细胞
    "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting",
    "NK_cells_resting", "Macrophages_M0", "Macrophages_M1", "Macrophages_M2",
    "Dendritic_cells_resting", "Mast_cells_resting", "Neutrophils",
    # 关键signature
    "CD8_T_cells", "T_cell_inflamed_GEP", "Immunogenicity", "Cytolytic_activity",
    # 自定义signature
    "HALLMARK_APOPTOSIS", "HALLMARK_ANGIOGENESIS"
)]

# 移除包含NA值的行
clustering_features <- na.omit(clustering_features)

# 2. 无监督聚类
set.seed(123)
# 使用IOBR的tme_cluster函数进行聚类
cluster_result <- tme_cluster(
    eset = t(clustering_features),
    cluster_method = "kmeans",
    k = 3 # 假设分为3个亚型
)

# 3. 热图可视化
# 准备热图数据
heatmap_data <- as.matrix(clustering_features)
heatmap_data <- t(scale(t(heatmap_data))) # 行标准化

# 创建注释
cluster_annotation <- data.frame(
    Cluster = factor(cluster_result$cluster),
    row.names = rownames(clustering_features)
)

# 颜色设置
cluster_colors <- brewer.pal(3, "Set1")
names(cluster_colors) <- unique(cluster_annotation$Cluster)

ha <- HeatmapAnnotation(
    Cluster = cluster_annotation$Cluster,
    col = list(Cluster = cluster_colors),
    annotation_name_side = "left"
)

# 创建热图
heatmap_plot <- Heatmap(
    heatmap_data,
    name = "Z-score",
    top_annotation = ha,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_title = "Features",
    column_title = "Samples"
)

# 保存热图
png(file.path(out_dir, paste0(gse_id, "_heatmap.png")),
    width = 12, height = 8, units = "in", res = 300
)
draw(heatmap_plot)
dev.off()

# 4. PCA可视化
pca_result <- prcomp(t(clustering_features), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Cluster <- factor(cluster_result$cluster)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = cluster_colors) +
    theme_minimal() +
    labs(
        title = "PCA Plot of TME Features",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(out_dir, paste0(gse_id, "_pca.png")),
    pca_plot,
    width = 8, height = 6, dpi = 300
)

# 5. 箱线图比较不同cluster的特征差异
# 选择几个关键特征进行比较
key_features <- c(
    "T_cells_CD8", "Macrophages_M1", "T_cell_inflamed_GEP",
    "HALLMARK_APOPTOSIS", "HALLMARK_ANGIOGENESIS"
)

plot_list <- list()
for (feature in key_features) {
    plot_df <- data.frame(
        Value = full_results[rownames(cluster_annotation), feature],
        Cluster = cluster_annotation$Cluster
    )

    p <- ggplot(plot_df, aes(x = Cluster, y = Value, fill = Cluster)) +
        geom_boxplot(alpha = 0.7) +
        scale_fill_manual(values = cluster_colors) +
        theme_minimal() +
        labs(title = feature, y = "Score") +
        theme(plot.title = element_text(hjust = 0.5))

    plot_list[[feature]] <- p
}

# 合并箱线图
boxplot_combined <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave(file.path(out_dir, paste0(gse_id, "_boxplots.png")),
    boxplot_combined,
    width = 12, height = 8, dpi = 300
)

cat("聚类分析和可视化完成!\n")