# 修正版 clustering_visualization.R
# 修复PCA图中的连续值错误

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

cat("开始聚类分析和可视化...\n")

# 读取IOBR分析结果
gse_id <- "GSE289743"
out_dir <- "data"

full_results <- readRDS(file.path(out_dir, paste0(gse_id, "_full_iobr_results.rds")))

cat("分析结果维度:", dim(full_results), "\n")

# 1. 识别可用的数值型特征用于聚类
numeric_columns <- sapply(full_results, is.numeric)
clustering_data <- full_results[, numeric_columns, drop = FALSE]

# 移除方差为0或接近0的列
col_vars <- apply(clustering_data, 2, var, na.rm = TRUE)
zero_var_cols <- col_vars == 0 | is.na(col_vars)
clustering_data <- clustering_data[, !zero_var_cols, drop = FALSE]

cat("可用于聚类的数值特征:", colnames(clustering_data), "\n")

# 2. 选择关键特征进行聚类
# 使用所有可用的数值特征
clustering_features <- clustering_data
cat("最终用于聚类的特征:\n")
print(colnames(clustering_features))

# 移除包含NA值的行
clustering_features <- na.omit(clustering_features)
cat("移除NA后用于聚类的样本数量:", nrow(clustering_features), "\n")

# 3. 无监督聚类
set.seed(123)

# 确定最佳聚类数（使用肘部法则）
max_k <- min(10, nrow(clustering_features) - 1)
wss <- sapply(1:max_k, function(k) {
    kmeans(clustering_features, centers = k, nstart = 25)$tot.withinss
})

# 绘制肘部法则图
elbow_data <- data.frame(k = 1:length(wss), wss = wss)
elbow_plot <- ggplot(elbow_data, aes(x = k, y = wss)) +
    geom_line() +
    geom_point() +
    labs(
        title = "Elbow Method for Optimal k",
        x = "Number of Clusters",
        y = "Total Within-Cluster Sum of Squares"
    ) +
    theme_minimal()

# 选择聚类数
optimal_k <- 3

# 执行K-means聚类
kmeans_result <- kmeans(clustering_features, centers = optimal_k, nstart = 25)
cluster_assignments <- kmeans_result$cluster

# 将聚类结果添加到原始数据中
full_results$Cluster <- NA
full_results[rownames(clustering_features), "Cluster"] <- as.factor(cluster_assignments)

cat("聚类完成，各簇样本数:\n")
print(table(full_results$Cluster, useNA = "always"))

# 4. 热图可视化
# 准备热图数据
heatmap_data <- as.matrix(clustering_features)
heatmap_data_scaled <- t(scale(t(heatmap_data))) # 行标准化

cat("热图数据维度:", dim(heatmap_data_scaled), "\n")
cat("热图数据列名数量:", length(colnames(heatmap_data_scaled)), "\n")

# 创建简单的聚类注释 - 确保与热图数据列名完全匹配
cluster_samples <- colnames(heatmap_data_scaled)
cluster_annotation <- data.frame(
    Cluster = full_results[cluster_samples, "Cluster"]
)
rownames(cluster_annotation) <- cluster_samples

cat("注释数据行数:", nrow(cluster_annotation), "\n")
cat("注释数据与热图数据列数是否一致:", nrow(cluster_annotation) == ncol(heatmap_data_scaled), "\n")

# 使用简单的颜色设置
cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")[1:optimal_k]

# 创建热图注释 - 确保使用正确的样本顺序
ha <- HeatmapAnnotation(
    df = cluster_annotation,
    col = list(Cluster = cluster_colors),
    annotation_name_side = "left",
    simple_anno_size = unit(0.5, "cm")
)

# 创建热图
heatmap_plot <- Heatmap(
    heatmap_data_scaled,
    name = "Z-score",
    top_annotation = ha,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_title = "Features",
    column_title = "Samples",
    row_names_max_width = unit(12, "cm")
)

# 保存热图
png(file.path(out_dir, paste0(gse_id, "_clustering_heatmap.png")),
    width = 14, height = 10, units = "in", res = 300
)
draw(heatmap_plot)
dev.off()

cat("热图保存成功!\n")

# 保存肘部法则图
ggsave(file.path(out_dir, paste0(gse_id, "_elbow_plot.png")),
    elbow_plot,
    width = 8, height = 6, dpi = 300
)

# 5. PCA可视化 - 修复连续值错误
pca_result <- prcomp(clustering_features, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Cluster <- full_results[rownames(pca_df), "Cluster"]

# 确保Cluster是因子类型
pca_df$Cluster <- as.factor(pca_df$Cluster)

cat("PCA数据中Cluster的类型:", class(pca_df$Cluster), "\n")
cat("PCA数据中Cluster的水平:", levels(pca_df$Cluster), "\n")

# PCA图 - 按聚类着色
pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = cluster_colors) +
    theme_minimal() +
    labs(
        title = "PCA Plot - Colored by Cluster",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5))

# 保存PCA图
ggsave(file.path(out_dir, paste0(gse_id, "_pca_cluster.png")),
    pca_cluster_plot,
    width = 8, height = 6, dpi = 300
)

# 6. 特征比较箱线图
# 选择前6个特征进行比较
key_features <- head(colnames(clustering_features), 6)

plot_list <- list()
for (feature in key_features) {
    plot_df <- data.frame(
        Value = full_results[rownames(clustering_features), feature],
        Cluster = full_results[rownames(clustering_features), "Cluster"]
    )
    plot_df <- na.omit(plot_df)

    # 确保Cluster是因子类型
    plot_df$Cluster <- as.factor(plot_df$Cluster)

    p <- ggplot(plot_df, aes(x = Cluster, y = Value, fill = Cluster)) +
        geom_boxplot(alpha = 0.7) +
        scale_fill_manual(values = cluster_colors) +
        theme_minimal() +
        labs(title = feature, y = "Score") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    plot_list[[feature]] <- p
}

# 合并箱线图
if (length(plot_list) > 0) {
    boxplot_combined <- wrap_plots(plot_list, ncol = 3) +
        plot_layout(guides = "collect")

    ggsave(file.path(out_dir, paste0(gse_id, "_feature_boxplots.png")),
        boxplot_combined,
        width = 12, height = 8, dpi = 300
    )
}

# 7. 单独的可视化 - 响应状态分布（如果可用）
if ("response" %in% colnames(full_results)) {
    response_data <- full_results[!is.na(full_results$response) & !is.na(full_results$Cluster), ]

    # 确保Cluster是因子类型
    response_data$Cluster <- as.factor(response_data$Cluster)

    response_dist_plot <- ggplot(
        response_data,
        aes(x = Cluster, fill = response)
    ) +
        geom_bar(position = "fill") +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        labs(
            title = "Response Distribution by Cluster",
            y = "Proportion"
        ) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(out_dir, paste0(gse_id, "_response_distribution.png")),
        response_dist_plot,
        width = 8, height = 6, dpi = 300
    )

    # 响应状态与聚类的交叉表
    response_table <- table(response_data$Cluster, response_data$response)
    cat("响应状态与聚类的交叉表:\n")
    print(response_table)
}

# 8. 单独的可视化 - 采样时间分布（如果可用）
if ("collection_time" %in% colnames(full_results)) {
    time_data <- full_results[!is.na(full_results$collection_time) & !is.na(full_results$Cluster), ]

    # 确保Cluster是因子类型
    time_data$Cluster <- as.factor(time_data$Cluster)

    time_dist_plot <- ggplot(
        time_data,
        aes(x = Cluster, fill = collection_time)
    ) +
        geom_bar(position = "fill") +
        scale_fill_brewer(palette = "Set3") +
        theme_minimal() +
        labs(
            title = "Collection Time Distribution by Cluster",
            y = "Proportion"
        ) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(out_dir, paste0(gse_id, "_time_distribution.png")),
        time_dist_plot,
        width = 8, height = 6, dpi = 300
    )
}

# 9. 聚类特征均值热图
# 计算每个簇的特征均值
cluster_means <- full_results %>%
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    column_to_rownames("Cluster")

# 选择用于展示的特征
display_features <- colnames(clustering_features)
if (length(display_features) > 15) {
    display_features <- head(display_features, 15)
}

cluster_means_display <- as.matrix(cluster_means[, display_features, drop = FALSE])
cluster_means_scaled <- t(scale(t(cluster_means_display)))

# 创建聚类均值热图
mean_heatmap <- Heatmap(
    cluster_means_scaled,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_title = "Clusters",
    column_title = "Features"
)

png(file.path(out_dir, paste0(gse_id, "_cluster_means_heatmap.png")),
    width = 12, height = 6, units = "in", res = 300
)
draw(mean_heatmap)
dev.off()

# 10. 保存聚类结果
saveRDS(full_results, file.path(out_dir, paste0(gse_id, "_clustered_results.rds")))
write.csv(full_results, file.path(out_dir, paste0(gse_id, "_clustered_results.csv")))

# 生成分析摘要
cat("\n=== 聚类分析完成摘要 ===\n")
cat("样本数量:", nrow(full_results), "\n")
cat("聚类数量:", optimal_k, "\n")
cat("各簇样本分布:\n")
print(table(full_results$Cluster, useNA = "always"))
cat("使用的特征数量:", ncol(clustering_features), "\n")
cat("主要特征:", paste(colnames(clustering_features), collapse = ", "), "\n")

cat("\n生成的可视化文件:\n")
cat("- 聚类热图:", file.path(out_dir, paste0(gse_id, "_clustering_heatmap.png")), "\n")
cat("- 肘部法则图:", file.path(out_dir, paste0(gse_id, "_elbow_plot.png")), "\n")
cat("- PCA图:", file.path(out_dir, paste0(gse_id, "_pca_cluster.png")), "\n")
cat("- 特征箱线图:", file.path(out_dir, paste0(gse_id, "_feature_boxplots.png")), "\n")
cat("- 聚类均值热图:", file.path(out_dir, paste0(gse_id, "_cluster_means_heatmap.png")), "\n")
if ("response" %in% colnames(full_results)) {
    cat("- 响应分布图:", file.path(out_dir, paste0(gse_id, "_response_distribution.png")), "\n")
}
if ("collection_time" %in% colnames(full_results)) {
    cat("- 时间分布图:", file.path(out_dir, paste0(gse_id, "_time_distribution.png")), "\n")
}

cat("聚类分析和可视化完成!\n")