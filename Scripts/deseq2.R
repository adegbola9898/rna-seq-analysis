
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("ReactomePA")

BiocManager::install("pathview")

# Install pheatmap package if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Load pheatmap library
library(pheatmap)

library(clusterProfiler)
library(ReactomePA)

library(pathview)

BiocManager::install("DESeq2")

BiocManager::install("tximport")

BiocManager::install("apeglm")

BiocManager::install("EnhancedVolcano")

BiocManager::install("org.Hs.eg.db")

install.packages("tidyverse")

library(DESeq2)

library(tidyverse)

library(tximport)
library(apeglm)
library(EnhancedVolcano)
library(org.Hs.eg.db)

count_data <- read.csv("/content/count_data.csv", row.names = 1)
head(count_data)

# Create metadata
metadata <- data.frame(
    row.names = colnames(count_data),
    condition = c(rep("Normal", 3), rep("HER2", 5), rep("NonTNBC", 6), rep("TNBC", 5))
)

# Inspect the metadata
metadata

# Convert non-integer values to integers by rounding
count_data <- round(count_data)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ condition)

# Pre-filtering (Optional) Remove rows with low counts to improve computation
dds <- dds[rowSums(counts(dds)) > 1, ]

# Normalize the data and run DESeq2 analysis
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Inspect the normalized counts
head(normalized_counts)

# Perform differential expression analysis
res_HER2_vs_Normal <- results(dds, contrast = c("condition", "HER2", "Normal"))
res_NonTNBC_vs_Normal <- results(dds, contrast = c("condition", "NonTNBC", "Normal"))
res_TNBC_vs_Normal <- results(dds, contrast = c("condition", "TNBC", "Normal"))

# Inspect results
summary(res_HER2_vs_Normal)
summary(res_NonTNBC_vs_Normal)
summary(res_TNBC_vs_Normal)

# Adjust p-values
res_HER2_vs_Normal <- res_HER2_vs_Normal[order(res_HER2_vs_Normal$padj), ]
res_NonTNBC_vs_Normal <- res_NonTNBC_vs_Normal[order(res_NonTNBC_vs_Normal$padj), ]
res_TNBC_vs_Normal <- res_TNBC_vs_Normal[order(res_TNBC_vs_Normal$padj), ]

# Define significance thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1

# HER2 vs Normal
sig_HER2_vs_Normal <- res_HER2_vs_Normal[which(res_HER2_vs_Normal$padj < padj_threshold & abs(res_HER2_vs_Normal$log2FoldChange) > log2fc_threshold), ]

# NonTNBC vs Normal
sig_NonTNBC_vs_Normal <- res_NonTNBC_vs_Normal[which(res_NonTNBC_vs_Normal$padj < padj_threshold & abs(res_NonTNBC_vs_Normal$log2FoldChange) > log2fc_threshold), ]

# TNBC vs Normal
sig_TNBC_vs_Normal <- res_TNBC_vs_Normal[which(res_TNBC_vs_Normal$padj < padj_threshold & abs(res_TNBC_vs_Normal$log2FoldChange) > log2fc_threshold), ]

# Inspect significant DEGs
head(sig_HER2_vs_Normal)
head(sig_NonTNBC_vs_Normal)
head(sig_TNBC_vs_Normal)

# Save HER2 vs Normal significant DEGs
write.csv(sig_HER2_vs_Normal, file = "sig_HER2_vs_Normal.csv", row.names = TRUE)

# Save NonTNBC vs Normal significant DEGs
write.csv(sig_NonTNBC_vs_Normal, file = "sig_NonTNBC_vs_Normal.csv", row.names = TRUE)

# Save TNBC vs Normal significant DEGs
write.csv(sig_TNBC_vs_Normal, file = "sig_TNBC_vs_Normal.csv", row.names = TRUE)

# MA Plot for HER2 vs Normal
plotMA(res_HER2_vs_Normal, ylim = c(-2, 2), main = "MA Plot: HER2 vs Normal")

# MA Plot for NonTNBC vs Normal
plotMA(res_NonTNBC_vs_Normal, ylim = c(-2, 2), main = "MA Plot: NonTNBC vs Normal")

# MA Plot for TNBC vs Normal
plotMA(res_TNBC_vs_Normal, ylim = c(-2, 2), main = "MA Plot: TNBC vs Normal")

# MA Plot for HER2 vs Normal
head(sig_HER2_vs_Normal)
head(res_HER2_vs_Normal)
plotMA(sig_HER2_vs_Normal, ylim = c(-2, 2), main = "MA Plot: HER2 vs Normal")
plotMA(res_HER2_vs_Normal, ylim = c(-2, 2), main = "MA Plot: HER2 vs Normal")
write.csv(res_HER2_vs_Normal, file = "res_HER2_vs_Normal.csv", row.names = TRUE)

# Extract the list of significant gene names
gene_list_HER2 <- rownames(sig_HER2_vs_Normal)
gene_list_NonTNBC <- rownames(sig_NonTNBC_vs_Normal)
gene_list_TNBC <- rownames(sig_TNBC_vs_Normal)

# List the valid keys for the SYMBOL keytype
valid_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
head(valid_keys)

# Print the head of each gene list
print("Head of gene_list_HER2:")
head(gene_list_HER2)

print("Head of gene_list_NonTNBC:")
head(gene_list_NonTNBC)

print("Head of gene_list_TNBC:")
head(gene_list_TNBC)

# Extract the list of significant gene names
gene_list_HER2 <- rownames(sig_HER2_vs_Normal)
gene_list_NonTNBC <- rownames(sig_NonTNBC_vs_Normal)
gene_list_TNBC <- rownames(sig_TNBC_vs_Normal)

# Convert Ensembl gene IDs to Entrez IDs
gene_list_HER2_entrez <- bitr(gene_list_HER2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_NonTNBC_entrez <- bitr(gene_list_NonTNBC, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_TNBC_entrez <- bitr(gene_list_TNBC, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check the conversion results
print(head(gene_list_HER2_entrez))
print(head(gene_list_NonTNBC_entrez))
print(head(gene_list_TNBC_entrez))

# Save original Ensembl gene lists to CSV files
write.csv(data.frame(Ensembl_ID = gene_list_HER2), file = "gene_list_HER2_ensembl.csv", row.names = FALSE)
write.csv(data.frame(Ensembl_ID = gene_list_NonTNBC), file = "gene_list_NonTNBC_ensembl.csv", row.names = FALSE)
write.csv(data.frame(Ensembl_ID = gene_list_TNBC), file = "gene_list_TNBC_ensembl.csv", row.names = FALSE)

# Save converted Entrez gene lists to CSV files
write.csv(gene_list_HER2_entrez, file = "gene_list_HER2_entrez.csv", row.names = FALSE)
write.csv(gene_list_NonTNBC_entrez, file = "gene_list_NonTNBC_entrez.csv", row.names = FALSE)
write.csv(gene_list_TNBC_entrez, file = "gene_list_TNBC_entrez.csv", row.names = FALSE)

# Calculate sample distances
sampleDists <- dist(t(assay(rld)))

# Convert to matrix for heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds)
colnames(sampleDistMatrix) <- colnames(dds)

# Plot heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# Load required libraries
library(DESeq2)
library(pheatmap)
library(tidyverse)

# Load count data (replace the file path with your actual path)
count_data <- read.csv("/content/count_data.csv", row.names = 1)

# Create metadata for the full dataset
metadata <- data.frame(
  row.names = colnames(count_data),
  condition = c(rep("Normal", 3), rep("HER2", 5), rep("NonTNBC", 6), rep("TNBC", 5))
)

# Convert non-integer values to integers by rounding
count_data <- round(count_data)

# Create DESeq2 dataset for the full dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)

# Pre-filtering (Optional) Remove rows with low counts to improve computation
dds <- dds[rowSums(counts(dds)) > 1, ]

# Normalize the data and run DESeq2 analysis
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Perform differential expression analysis for each condition vs Normal
res_HER2_vs_Normal <- results(dds, contrast = c("condition", "HER2", "Normal"))
res_NonTNBC_vs_Normal <- results(dds, contrast = c("condition", "NonTNBC", "Normal"))
res_TNBC_vs_Normal <- results(dds, contrast = c("condition", "TNBC", "Normal"))

# Define significance thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1

# Extract significant results
sig_HER2_vs_Normal <- subset(res_HER2_vs_Normal, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
sig_NonTNBC_vs_Normal <- subset(res_NonTNBC_vs_Normal, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
sig_TNBC_vs_Normal <- subset(res_TNBC_vs_Normal, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)

# Top genes for each comparison
top_genes_HER2 <- rownames(head(sig_HER2_vs_Normal[order(sig_HER2_vs_Normal$padj), ], 50))
top_genes_NonTNBC <- rownames(head(sig_NonTNBC_vs_Normal[order(sig_NonTNBC_vs_Normal$padj), ], 50))
top_genes_TNBC <- rownames(head(sig_TNBC_vs_Normal[order(sig_TNBC_vs_Normal$padj), ], 50))

# Extract relevant samples for each comparison
her2_samples <- c("NBS1_read", "NBS2", "NBS3", "X26", "IP2.53", "X56_s", "X83", "X171")
nontnbc_samples <- c("NBS1_read", "NBS2", "NBS3", "IP2.42", "IP2.48", "IP2.49", "IP2.65", "IP2.66", "IP2.71")
tnbc_samples <- c("NBS1_read", "NBS2", "NBS3", "IP2.50", "IP2.76", "IP2.78", "IP2.83", "IP2.90")

# Filter the normalized counts for top genes and relevant samples
normalized_counts_HER2 <- normalized_counts[top_genes_HER2, her2_samples]
normalized_counts_NonTNBC <- normalized_counts[top_genes_NonTNBC, nontnbc_samples]
normalized_counts_TNBC <- normalized_counts[top_genes_TNBC, tnbc_samples]

# Scale the counts
scaled_counts_HER2 <- t(scale(t(normalized_counts_HER2)))
scaled_counts_NonTNBC <- t(scale(t(normalized_counts_NonTNBC)))
scaled_counts_TNBC <- t(scale(t(normalized_counts_TNBC)))

# Create metadata for relevant samples
metadata_HER2 <- metadata[her2_samples, , drop=FALSE]
metadata_NonTNBC <- metadata[nontnbc_samples, , drop=FALSE]
metadata_TNBC <- metadata[tnbc_samples, , drop=FALSE]

# Create heatmaps
pheatmap(scaled_counts_HER2, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         annotation_col = metadata_HER2,
         main = "Heatmap: Top 50 Differentially Expressed Genes (HER2 vs Normal)")

pheatmap(scaled_counts_NonTNBC, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         annotation_col = metadata_NonTNBC,
         main = "Heatmap: Top 50 Differentially Expressed Genes (NonTNBC vs Normal)")

pheatmap(scaled_counts_TNBC, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         annotation_col = metadata_TNBC,
         main = "Heatmap: Top 50 Differentially Expressed Genes (TNBC vs Normal)")



# Extract the list of significant gene names
gene_list_HER2 <- rownames(sig_HER2_vs_Normal)
gene_list_NonTNBC <- rownames(sig_NonTNBC_vs_Normal)
gene_list_TNBC <- rownames(sig_TNBC_vs_Normal)


# Print the head of each gene list
print("Head of gene_list_HER2:")
head(gene_list_HER2)

print("Head of gene_list_NonTNBC:")
head(gene_list_NonTNBC)

print("Head of gene_list_TNBC:")
head(gene_list_TNBC)


# Extract the list of significant gene names
gene_list_HER2 <- rownames(sig_HER2_vs_Normal)
gene_list_NonTNBC <- rownames(sig_NonTNBC_vs_Normal)
gene_list_TNBC <- rownames(sig_TNBC_vs_Normal)

# Convert Ensembl gene IDs to Entrez IDs
gene_list_HER2_entrez <- bitr(gene_list_HER2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_NonTNBC_entrez <- bitr(gene_list_NonTNBC, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_TNBC_entrez <- bitr(gene_list_TNBC, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check the conversion results
print(head(gene_list_HER2_entrez))
print(head(gene_list_NonTNBC_entrez))
print(head(gene_list_TNBC_entrez))

# MA Plot
plotMA(res_HER2_vs_Normal, main = "MA Plot: HER2 vs Normal", ylim = c(-5, 5))
plotMA(res_NonTNBC_vs_Normal, main = "MA Plot: NonTNBC vs Normal", ylim = c(-5, 5))
plotMA(res_TNBC_vs_Normal, main = "MA Plot: TNBC vs Normal", ylim = c(-5, 5))

# Volcano Plot
EnhancedVolcano(res_HER2_vs_Normal,
                lab = rownames(res_HER2_vs_Normal),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: HER2 vs Normal',
                pCutoff = padj_threshold,
                FCcutoff = log2fc_threshold)

EnhancedVolcano(res_NonTNBC_vs_Normal,
                lab = rownames(res_NonTNBC_vs_Normal),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: NonTNBC vs Normal',
                pCutoff = padj_threshold,
                FCcutoff = log2fc_threshold)

EnhancedVolcano(res_TNBC_vs_Normal,
                lab = rownames(res_TNBC_vs_Normal),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: TNBC vs Normal',
                pCutoff = padj_threshold,
                FCcutoff = log2fc_threshold)

# PCA Plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot of Normalized Counts")

# Perform GO enrichment analysis
go_enrich_HER2 <- enrichGO(gene = gene_list_HER2_entrez$ENTREZID,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05,
                           readable = TRUE)

# View results
head(go_enrich_HER2)

# Perform KEGG enrichment analysis
kegg_enrich_HER2 <- enrichKEGG(gene = gene_list_HER2_entrez$ENTREZID,
                               organism = 'hsa',
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)

# View results
head(kegg_enrich_HER2)

# Perform Reactome pathway enrichment analysis
reactome_enrich_HER2 <- enrichPathway(gene = gene_list_HER2_entrez$ENTREZID,
                                      organism = "human",
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 0.05,
                                      readable = TRUE)

# View results
head(reactome_enrich_HER2)

# Dotplot for GO enrichment results
dotplot(go_enrich_HER2, showCategory = 20) + ggtitle("GO Enrichment Analysis for HER2 vs Normal")

# Pathway view for KEGG enrichment results
pathview(gene.data = sig_HER2_vs_Normal$log2FoldChange,
         pathway.id = "hsa05200",
         species = "hsa")

# Barplot for Reactome enrichment results
barplot(reactome_enrich_HER2, showCategory = 20) + ggtitle("Reactome Enrichment Analysis for HER2 vs Normal")
