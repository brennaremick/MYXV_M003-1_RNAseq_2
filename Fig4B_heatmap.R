library(org.Hs.eg.db)
library(pheatmap)

# Make PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))

# Extract the vst transformed expression matrix
vst_matrix <- assay(vsd)

# Identify the top 30 genes by variance
gene_variances <- apply(vst_matrix, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:30]
top_genes_df <- data.frame(
  Gene = top_genes,
  Variance = gene_variances[top_genes]
)

# include NF-kB target genes TNF and IL6
# ENSG00000232810.5 = TNF
# ENSG00000136244.12 = IL6
top_genes_TNF_IL6 <- top_genes_df$Gene
top_genes_TNF_IL6 <- c(top_genes_TNF_IL6, "ENSG00000232810.5", "ENSG00000136244.12")
top_30_vst <- vst_matrix[top_genes_TNF_IL6, ]

# Compute Z-scores for each gene (row-wise standardization)
z_scores <- as.data.frame(t(scale(t(top_30_vst))))

# replace ENSG gene names with gene symbol
ens.str <- substr((rownames(z_scores)), 1, 15)
z_scores$symbol <- mapIds(org.Hs.eg.db,
                          keys=ens.str,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
rownames(z_scores) <- z_scores$symbol
z_scores <- z_scores[, -ncol(z_scores)]
colnames(z_scores) <- sample_info[colnames(z_scores), "condition"]

# subset the columns of interest to make figure & specify desired column order
colnames(z_scores) <- make.unique(colnames(z_scores))
desired_order <- c("NTC_cherry","NTC_cherry.1","NTC_cherry.2","NTC_M003","NTC_M003.1","NTC_M003.2", "ZC3H12A_N4BP1_cherry","ZC3H12A_N4BP1_cherry.1","ZC3H12A_N4BP1_cherry.2", "ZC3H12A_N4BP1_M003","ZC3H12A_N4BP1_M003.1", "ZC3H12A_N4BP1_M003.2", "ZC3H12A_cherry","ZC3H12A_cherry.1", "ZC3H12A_cherry.2", "ZC3H12A_M003","ZC3H12A_M003.1","ZC3H12A_M003.2","N4BP1_cherry","N4BP1_cherry.1","N4BP1_cherry.2","N4BP1_M003","N4BP1_M003.1","N4BP1_M003.2" ) 

# Reorder the columns of the matrix based on the desired order
z_scores_ordered <- z_scores[, desired_order]

# Create a heatmap
pheatmap(
  z_scores_ordered,
  cluster_rows = TRUE,          # Cluster genes
  cluster_cols = FALSE,          
  scale = "none",               # No additional scaling (already Z-scores)
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  show_rownames = TRUE,         
  show_colnames = TRUE,        
  main = "Heatmap of Z-scores for Top 30 Genes by Variance + TNF + IL6",
)

