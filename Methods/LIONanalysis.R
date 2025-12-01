# =================================================================

# =================================================================

rm(list = ls())
library(dplyr)
set.seed(42)

researchAim = '2String9'

available_datasets <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                       "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

source('/proj/c.zihao/work1/1NT/function/LIONfunction.R')

for(dataset_name in 1:length(available_datasets)) {

codestart_time <- Sys.time()
diease_name <- available_datasets[dataset_name]
print(diease_name)

base_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/', diease_name, "/data")
save_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/', diease_name, "/LION")
ppiFile <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/links.csv')

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

exprSetFile <- paste0(base_path, '/', "exprSet_filtered.csv")
print(exprSetFile)

resultDF <- LIONcal(exprSetFile, ppiFile)

print(save_path)
setwd(save_path)
write.csv(resultDF, file = "result.csv", row.names = FALSE)
print("All sample for LION completed")

codeend_time <- Sys.time()
elapsed_time <- as.numeric(difftime(codeend_time, codestart_time, units = "secs"))
cat("Execution time:", elapsed_time, "\n")
write(
  paste("Start:", codestart_time,
        "| End:", codeend_time,
        "| Duration:", sprintf("%.3f seconds", elapsed_time)),
  file = "runtime_log.txt",
  append = TRUE
)
cat("============================================")

 }