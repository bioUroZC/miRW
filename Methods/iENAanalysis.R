# =================================================================

# =================================================================

rm(list = ls())
library(dplyr)
set.seed(42)

researchAim = '2String9'

dataloop <- data.frame(diease = c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
                                  "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD"),
                       organ = c("Bladder", "Breast", "Colon", "Esophagus","Salivary Gland", 
                                 "Kidney", "Liver", "Lung", "Lung", "Prostate", "Stomach")
                       )
source('/proj/c.zihao/work1/1NT/function/iENAfunction.R')


for(lo in 1:dim(dataloop)[1]) {

codestart_time <- Sys.time()

diease_name <- dataloop$diease[lo]
organ_name  <- dataloop$organ[lo]

print(diease_name)
print(organ_name)

base_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/', diease_name, "/data")
save_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/', diease_name, "/iENA")
ppiFile <- paste0("/proj/c.zihao/work1/1NT/", researchAim, '/links.csv')

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

exprSetFile <- paste0(base_path, '/', "exprSet_filtered.csv")
print(exprSetFile)

NomrlaFile <- "/proj/c.zihao/work1/0ref/GTEx/combined_expr_df.csv"

resultDF <- iENAcal(exprSetFile, NomrlaFile, ppiFile, organ_name)
print(save_path)
setwd(save_path)

write.csv(resultDF, "result_matrix.csv", row.names = FALSE) 
print("All sample for iENA completed")
cat("============================================")


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