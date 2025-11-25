# 下载 GSE289743 的 Series Matrix 与 Supplementary counts 文件，并保存到 data/
library(GEOquery)
library(data.table)

gse_id <- "GSE289743"
out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE)

# 1. 下载 Series Matrix（pheno）
cat("下载 Series Matrix ......\n")
gse_list <- getGEO(gse_id, GSEMatrix = TRUE)
gse <- gse_list[[1]]
pheno <- pData(gse)
saveRDS(pheno, file = file.path(out_dir, paste0(gse_id, "_pheno.rds")))
write.csv(pheno, file = file.path(out_dir, paste0(gse_id, "_pheno.csv")), row.names = TRUE)
cat("Series Matrix 下载并保存到", out_dir, "\n")

# 2. 下载 supplementary files（counts）
cat("下载 Supplementary files ......\n")
getGEOSuppFiles(gse_id, makeDirectory = TRUE, baseDir = out_dir)
# 找到可能的 counts 文件（按常见命名）
supp_dir <- file.path(out_dir, gse_id)
files <- list.files(supp_dir, full.names = TRUE)
cat("Supplement files:\n")
print(files)

# 常见: 包含 Raw_counts 或 counts 的 CSV/GZ 文件，挑选第一个匹配项
counts_file <- files[grepl("counts|raw_counts|Raw_counts|csv", files, ignore.case = TRUE)][1]
if (is.na(counts_file)) {
    stop("未在 supplementary files 中找到 counts 文件，请手动检查:", supp_dir)
}
cat("使用 counts 文件：", counts_file, "\n")

# 读取 counts（支持 .gz）
if (grepl("\\.gz$", counts_file)) {
    counts_df <- fread(cmd = paste("gzip -dc", shQuote(counts_file)))
} else {
    counts_df <- fread(counts_file)
}
# 保存原始 counts 为 rds 以便后续处理
saveRDS(counts_df, file = file.path(out_dir, paste0(gse_id, "_raw_counts.rds")))
write.csv(counts_df, file = file.path(out_dir, paste0(gse_id, "_raw_counts.csv")), row.names = FALSE)
cat("原始 counts 保存到", out_dir, "\n")
