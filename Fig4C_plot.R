library(ggplot2)
library(ggrepel)

# Compare differential gene expression in gRNA:NTC +M003.1 (vs gRNA:NTC +mCherry) and gRNA:ZC3H12A/N4BP1 +mCherry (vs gRNA:NTC +mCherry)

NT_M003 <- results(dds, contrast = c("condition","NTC_M003", "NTC_cherry"))
ZCN4_cherry <- results(dds, contrast = c("condition","ZC3H12A_N4BP1_cherry", "NTC_cherry"))

# Create a combined data frame
df <- data.frame(
  gene = rownames(NT_M003),
  log2FC_NT_M003 = NT_M003$log2FoldChange,
  log2FC_ZCN4_cherry = ZCN4_cherry$log2FoldChange,
  padj_NT_M003 = NT_M003$padj,
  padj_ZCN4_cherry = ZCN4_cherry$padj
)

# Remove NA values
df <- na.omit(df)

# Categorize genes based on log2FC and adjusted p-value
# Use log2FC cutoff of 1 and -1
# Use adjusted p-value cutoff of 0.05

df$category <- "Non-significant"
df$category[
  (df$log2FC_NT_M003 > 1 & df$padj_NT_M003 < 0.05 & 
     df$log2FC_ZCN4_cherry > 1 & df$padj_ZCN4_cherry < 0.05) |
    (df$log2FC_NT_M003 < -1 & df$padj_NT_M003 < 0.05 & 
       df$log2FC_ZCN4_cherry < -1 & df$padj_ZCN4_cherry < 0.05)
] <- "Up or down in both"
df$category[
  (df$log2FC_NT_M003 > 1 & df$padj_NT_M003 < 0.05 & 
     (is.na(df$log2FC_ZCN4_cherry) | df$log2FC_ZCN4_cherry <= 1 | df$padj_ZCN4_cherry >= 0.05)) |
    (df$log2FC_NT_M003 < -1 & df$padj_NT_M003 < 0.05 & 
       (is.na(df$log2FC_ZCN4_cherry) | df$log2FC_ZCN4_cherry >= -1 | df$padj_ZCN4_cherry >= 0.05))
] <- "Up or down in NT_M003 only"
df$category[
  (df$log2FC_ZCN4_cherry > 1 & df$padj_ZCN4_cherry < 0.05 & 
     (is.na(df$log2FC_NT_M003) | df$log2FC_NT_M003 <= 1 | df$padj_NT_M003 >= 0.05)) |
    (df$log2FC_ZCN4_cherry < -1 & df$padj_ZCN4_cherry < 0.05 & 
       (is.na(df$log2FC_NT_M003) | df$log2FC_NT_M003 >= -1 | df$padj_NT_M003 >= 0.05))
] <- "Up or down in ZCN4 only"

# Create four-way plot comparing differential gene expression
plot <- ggplot(df, aes(x = log2FC_NT_M003, y = log2FC_ZCN4_cherry, color = category)) +
  geom_point(alpha = 0.8, size = 3, stroke=0, shape = 16) +
  scale_color_manual(values = c(
    "Non-significant" = "gray",
    "Up or down in both" = "#4DAC26",
    "Up or down in NT_M003 only" = "#0571B0",
    "Up or down in ZCN4 only" = "#FDB863"
  )) +
  geom_hline(yintercept = c(0, 5, -5, 10, -10), linetype = "solid", color = "black") +
  geom_vline(xintercept = c(0,5, -5, 10, -10), linetype = "solid", color = "black") +
  geom_hline(yintercept = c(1, -1), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +
  labs(
    title = "Comparison of Differential Gene Expression",
    x = "Log2FC (NTC_M003.1 vs NTC_cherry)",
    y = "Log2FC (ZC3H12A_N4BP1_cherry vs NTC_cherry)",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  ) +
  coord_fixed(ratio = 1) +
  xlim(-13, 13) + 
  ylim(-13, 13) 

print(plot)