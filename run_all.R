# 一键运行脚本（按顺序执行）
scripts <- c("install_packages.R", "check_packages.R", "download_data.R", "preprocess.R", "analysis_iobr.R")
for (s in scripts) {
    cat("===== Running", s, "=====\n")
    source(s)
}
cat("全部脚本已执行完毕\n")
