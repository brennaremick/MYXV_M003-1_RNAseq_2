library(DESeq2)

# load in the data
sample_info <- read.table("RVBR004_sample_info.txt", row.names=1)
counts <- read.table("raw_feature_counts2.txt", row.names = 1, stringsAsFactors=FALSE)

# perform differential gene expression with DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

# filter to remove low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
