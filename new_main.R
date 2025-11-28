##############################################################################
# iobr_visualize_and_complementarity.R
# 功能：
# 1) 统一对齐 signature scores 与 deconv 结果（自动处理行列方向）
# 2) 生成热图（signature）、箱线/violin（按分组若有 pheno）、signature-vs-cell 相关性图
# 3) 下载 MSigDB Hallmark（若可联网）作为自定义 signature，计算并加入分析
# 4) 计算“互补性”度量（自定义 signature 与内建签名相关性分布、聚类）
# 输出：results_iobr_analysis 下的多个图片和 csv
##############################################################################

# 载入包
pkgs <- c("data.table", "dplyr", "ggplot2", "pheatmap", "gridExtra", "corrplot", "msigdbr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
    library(gridExtra)
    library(corrplot)
    library(msigdbr)
})

outdir <- "results_iobr_analysis"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 读取 signature scores 与 deconv 结果
sig_path <- file.path(outdir, "signature_scores_all.rds")
if (file.exists(sig_path)) {
    sig_all <- readRDS(sig_path) # 期望为 matrix: rows = signatures, cols = samples
} else if (file.exists(file.path(outdir, "signature_scores_all.csv"))) {
    sig_all <- as.matrix(fread(file.path(outdir, "signature_scores_all.csv"), data.table = FALSE))
    # 如果 CSV 第一列为行名问题，需要调整；尝试自动处理：
    if (colnames(sig_all)[1] == "ID") {
        rownames(sig_all) <- sig_all[, 1]
        sig_all <- sig_all[, -1]
    }
    sig_all <- as.matrix(sig_all)
} else {
    stop("找不到 signature_scores_all，请先运行 signature 打分并保存。")
}

# 找到第一个可用 deconv 结果（rds）
deconv_files <- list.files(outdir, pattern = "^deconv_.*\\.rds$", full.names = TRUE)
deconv_list <- list()
for (f in deconv_files) {
    obj <- readRDS(f)
    # 标准化：若返回为 data.frame/matrix samples x celltypes 或反向，均转换为 samples x celltypes
    if (is.matrix(obj) || is.data.frame(obj)) {
        mat <- as.matrix(obj)
        # 如果行名是 celltypes 而列是 samples（检测 colnames 是否为样本名）
        # 我们判断样本是否在 sig_all 的列中
        if (!is.null(colnames(mat)) && any(colnames(mat) %in% colnames(sig_all))) {
            # colnames 是样本，假设 mat: celltypes x samples? Actually celltypes x samples would have samples as colnames
            # we want samples x celltypes => transpose
            if (all(colnames(mat) %in% colnames(sig_all))) {
                mat2 <- t(mat) # now rows = samples, cols = celltypes
            } else {
                mat2 <- mat
            }
        } else if (!is.null(rownames(mat)) && any(rownames(mat) %in% colnames(sig_all))) {
            # rows are samples; good
            if (all(rownames(mat) %in% colnames(sig_all))) mat2 <- mat else mat2 <- mat
        } else {
            # cannot infer; keep as-is
            mat2 <- mat
        }
        deconv_list[[basename(f)]] <- mat2
    }
}
if (length(deconv_list) == 0) {
    warning("未读取到任何 deconv 结果文件（deconv_*.rds）。")
}

# 选择首个可对齐的方法用于示例相关性分析
chosen_method <- NULL
for (nm in names(deconv_list)) {
    mat <- deconv_list[[nm]]
    # mat rows should be samples
    samples_in_common <- intersect(rownames(mat), colnames(sig_all))
    if (length(samples_in_common) > 0) {
        chosen_method <- nm
        chosen_mat <- mat
        break
    }
}
if (is.null(chosen_method)) {
    # 作为回退，尝试列名匹配
    for (nm in names(deconv_list)) {
        mat <- deconv_list[[nm]]
        samples_in_common <- intersect(colnames(mat), colnames(sig_all))
        if (length(samples_in_common) > 0) {
            chosen_method <- nm
            chosen_mat <- t(mat)
            break
        }
    }
}
if (is.null(chosen_method)) {
    message("无法对齐任何 deconv 结果与 signature 分数（样本名不匹配）。相关性分析将被跳过。")
} else {
    message("用于相关性分析的 deconv 文件：", chosen_method)
}

# -------------------- 绘制 signature 热图（行归一化） --------------------
# 取方差前 40 个 signature
sig_var <- apply(sig_all, 1, var, na.rm = TRUE)
top_sigs <- names(sort(sig_var, decreasing = TRUE))[1:min(40, length(sig_var))]
hm_mat <- sig_all[top_sigs, , drop = FALSE]
# 若有 pheno 可做 annotation
pheno_path <- "results_preprocess/pheno_aligned.rds"
anno_col <- NULL
if (file.exists(pheno_path)) {
    pheno <- readRDS(pheno_path)
    # 选择部分列显示（自动）
    cand <- grep("response|treatment|group|time|tissue", colnames(pheno), value = TRUE, ignore.case = TRUE)
    if (length(cand) > 0) {
        anno_col <- pheno[, cand, drop = FALSE]
        rownames(anno_col) <- rownames(pheno)
    }
}
png(file.path(outdir, "signature_heatmap_top40.png"), width = 2000, height = 1200, res = 200)
pheatmap::pheatmap(hm_mat,
    scale = "row", annotation_col = anno_col,
    clustering_method = "ward.D2", show_rownames = TRUE, show_colnames = FALSE,
    main = "Top 40 signatures (row-scaled)"
)
dev.off()

# -------------------- 若有 pheno：分组箱线/violin 图 --------------------
if (file.exists(pheno_path)) {
    pheno <- readRDS(pheno_path)
    # 尝试找到分组列（response 优先）
    grp_cols <- grep("response|Response|group|treatment", colnames(pheno), value = TRUE, ignore.case = TRUE)
    if (length(grp_cols) > 0) {
        grp_col <- grp_cols[1]
        # 画若干 signature（选择用户关注的或方差最高的）
        sigs_to_plot <- c("cytotoxic", "ifn_gamma")
        sigs_to_plot <- sigs_to_plot[sigs_to_plot %in% rownames(sig_all)]
        if (length(sigs_to_plot) == 0) sigs_to_plot <- rownames(sig_all)[1:2]
        plot_list <- list()
        for (sg in sigs_to_plot) {
            df <- data.frame(sample = colnames(sig_all), score = as.numeric(sig_all[sg, colnames(sig_all)]))
            df$group <- as.factor(pheno[colnames(sig_all), grp_col])
            p <- ggplot(df, aes(x = group, y = score, fill = group)) +
                geom_violin(trim = FALSE) +
                geom_boxplot(width = 0.2, outlier.size = 0.6) +
                ggtitle(paste("Signature:", sg)) +
                theme_bw()
            plot_list[[sg]] <- p
        }
        ggsave(file.path(outdir, "sig_violin_plots.png"), plot = gridExtra::grid.arrange(grobs = plot_list, ncol = 2), width = 10, height = 5, dpi = 300)
    } else {
        message("pheno 存在但未检测到分组列，跳过分组箱线图。")
    }
} else {
    message("未找到 pheno_aligned.rds，跳过分组箱线图（如需分组比较请提供 pheno）。")
}

# -------------------- 若能对齐 deconv，计算 signature vs cell correlation --------------------
if (exists("chosen_mat") && !is.null(chosen_mat)) {
    # 确保 chosen_mat 行 = samples, cols = celltypes
    cellmat <- chosen_mat
    if (!is.null(colnames(cellmat)) && all(colnames(cellmat) %in% colnames(sig_all))) {
        cellmat <- t(cellmat) # adjust if needed (some earlier code may have transposed)
    }
    # ensure samples align: use intersection
    samples_common <- intersect(rownames(cellmat), colnames(sig_all))
    sig_sub <- sig_all[, samples_common, drop = FALSE]
    cell_sub <- as.matrix(cellmat[samples_common, , drop = FALSE]) # rows = samples
    # compute Spearman correlation between signatures (rows) and cell types (columns)
    cor_mat <- cor(t(sig_sub), t(cell_sub), use = "pairwise.complete.obs", method = "spearman")
    write.csv(cor_mat, file = file.path(outdir, paste0("sig_vs_", gsub("\\.rds", "", chosen_method), "_spearman.csv")), row.names = TRUE)
    png(file.path(outdir, paste0("sig_vs_", gsub("\\.rds", "", chosen_method), "_spearman.png")), width = 1600, height = 1200, res = 200)
    corrplot::corrplot(cor_mat, method = "color", tl.cex = 0.6, cl.cex = 0.7, number.cex = 0.6, title = paste("Signature vs", gsub("\\.rds", "", chosen_method)))
    dev.off()
} else {
    message("没有可对齐的 deconv 结果，跳过 signature-vs-cell 相关性绘图。")
}

# -------------------- 载入 MSigDB Hallmark（如可联网）作为自定义 signature --------------------
# 需要 msigdbr 包（上面已尝试安装）；若无法联网会失败
message("尝试获取 MSigDB Hallmark 作为自定义 signature（需要联网）。")
hm <- try(msigdbr(species = "Homo sapiens", category = "H"), silent = TRUE)
if (inherits(hm, "try-error")) {
    message("无法从 MSigDB 获取 Hallmark 名单，可能没有网络或 msigdbr 无法工作。")
    custom_hallmark <- NULL
} else {
    # 以 "HALLMARK_X" 分组为 signature list
    custom_hallmark <- split(hm$gene_symbol, hm$gs_name)
    # 选取 1-2 个示例 signature（如 IFN_GAMMA_RESPONSE / INTERFERON_ALPHA_RESPONSE）
    want <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
    want <- want[want %in% names(custom_hallmark)]
    if (length(want) == 0) want <- names(custom_hallmark)[1:2]
    custom_subset <- custom_hallmark[want]
    # 计算这两个 signature 分数（ssGSEA）
    library(GSVA)
    message("为 Hallmark 计算 ssGSEA 分数（可能较慢）...")
    hm_scores <- gsva::gsva(expr = as.matrix(tpm), gset.idx.list = custom_subset, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
    # gsva 返回 matrix: rows=signature, cols=samples
    # 将它加入现有 sig_all（对样本取交集并合并）
    common_samples2 <- intersect(colnames(sig_all), colnames(hm_scores))
    if (length(common_samples2) > 0) {
        sig_all <- cbind(sig_all[, common_samples2, drop = FALSE], hm_scores[, common_samples2, drop = FALSE])
        write.csv(as.data.frame(hm_scores), file = file.path(outdir, "hallmark_custom_scores.csv"), row.names = TRUE)
        message("已计算并保存 Hallmark 自定义 signature 分数。")
    } else {
        message("Hallmark 分数与现有 signature 无样本交集，未合并。")
    }
}

# -------------------- 计算自定义 signature 与 原 signature 的互补性 --------------------
# 互补性度量（示例）：计算自定义 sig 与每个原 signature 的 Spearman 绝对相关性分布，
# 若相关性低（接近 0）则说明互补性高；我们以 |rho| 的中位数或最大值评估
if (!is.null(custom_subset) && exists("sig_all")) {
    # get names of custom
    custom_names <- names(custom_subset)
    custom_names <- custom_names[custom_names %in% rownames(sig_all)]
    if (length(custom_names) == 0) {
        message("没有在 sig_all 中找到自定义 Hallmark 的名称，跳过互补性计算。")
    } else {
        # compute correlations (spearman) between custom signatures and all other signatures
        other_sigs <- setdiff(rownames(sig_all), custom_names)
        comp_res <- matrix(NA, nrow = length(custom_names), ncol = length(other_sigs), dimnames = list(custom_names, other_sigs))
        for (cn in custom_names) {
            for (os in other_sigs) {
                comp_res[cn, os] <- cor(as.numeric(sig_all[cn, ]), as.numeric(sig_all[os, ]), use = "pairwise.complete.obs", method = "spearman")
            }
        }
        # 统计|rho|的分布：中位数、最大值
        comp_summary <- data.frame(
            custom = custom_names,
            median_abs_rho = apply(abs(comp_res), 1, median, na.rm = TRUE),
            max_abs_rho = apply(abs(comp_res), 1, max, na.rm = TRUE)
        )
        write.csv(comp_summary, file = file.path(outdir, "custom_signature_complementarity_summary.csv"), row.names = FALSE)
        message("已保存自定义 signature 与原 signature 的互补性汇总：", file.path(outdir, "custom_signature_complementarity_summary.csv"))
        # 画出热图（custom x top_N most correlated original signatures）
        topN <- 30
        for (cn in custom_names) {
            top_idx <- order(abs(comp_res[cn, ]), decreasing = TRUE)[1:min(topN, ncol(comp_res))]
            subm <- comp_res[cn, colnames(comp_res)[top_idx], drop = FALSE]
            png(file.path(outdir, paste0("comp_heat_", cn, ".png")), width = 1400, height = 600, res = 200)
            pheatmap::pheatmap(t(subm), main = paste("Correlations with", cn), cluster_rows = TRUE, cluster_cols = TRUE)
            dev.off()
        }
    }
} else {
    message("没有自定义 Hallmark 可用或 sig_all 不存在，跳过互补性分析。")
}

message("可视化与互补性分析完成。主要文件保存在：", normalizePath(outdir))
