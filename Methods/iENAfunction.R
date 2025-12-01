iENAcal <- function (exprSetFile, NomrlaFile, ppiFile, organ){

cancer <- read.csv(exprSetFile, header = TRUE, row.names = 1)
cancer <- cancer[apply(cancer, 1, sd) > 0, ]

# ======================================================================================

# ======================================================================================

Normal <- read.csv(NomrlaFile, row.names = 1)
Normal <- Normal[which(Normal$organ == organ), ]
Normal$organ <- NULL
Normal <- as.data.frame(t(Normal))
Normal <- Normal[apply(Normal, 1, sd) > 0, ]

interGenes <- intersect(rownames(cancer), rownames(Normal))

# ======================================================================================
# Step 2: Load and filter PPI ----------------------
# ======================================================================================

ppi <- read.csv(ppiFile, header = TRUE)
ppi <- ppi[ppi$protein1 %in% interGenes & ppi$protein2 %in% interGenes, ]
ppi <- unique(ppi[, c("protein1", "protein2")])
head(ppi)
colnames(ppi) <- c("Gene1", "Gene2")
rownames(ppi) <- NULL

# Step 3: Compute reference mean and variance ------

Normal_use <- Normal[interGenes, ]
mean_ref <- apply(Normal_use, 1, mean)
var_ref <- apply(Normal_use, 1, var)

# Step 4: Initialize result matrix -----------------

result_df <- ppi
sample_names <- colnames(cancer)
for (s in sample_names) {
  result_df[[s]] <- NA
}

# Step 5: Compute sPCC for each sample -------------

for (sample_idx in 1:length(sample_names)) {
  sname <- sample_names[sample_idx]
  percent_done <- round(100 * sample_idx / length(sample_names), 1)
  cat(sprintf("Processing sample %s (%d of %d, %.1f%% done)\n",
              sname, sample_idx, length(sample_names), percent_done))
  
  sample_vec <- cancer[, sname]
  names(sample_vec) <- rownames(cancer)
  
  for (row_id in 1:nrow(ppi)) {
    g1 <- ppi$Gene1[row_id]
    g2 <- ppi$Gene2[row_id]
    
    xi <- sample_vec[g1]
    xj <- sample_vec[g2]
    mi <- mean_ref[g1]
    mj <- mean_ref[g2]
    vi <- var_ref[g1]
    vj <- var_ref[g2]
    
    cov_ij <- (xi - mi) * (xj - mj)
    denom <- sqrt(vi * vj) + 1e-8
    sPCC <- cov_ij / denom
    
    result_df[row_id, sname] <- sPCC
  }
}

return(result_df)

}