
# WGCNA Analysis for Breast Cancer Subtype Classification

# This script performs Weighted Gene Co-expression Network Analysis (WGCNA)
# to identify gene modules associated with breast cancer subtypes
 


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages and WGCNA (as BiocManager handles its dependencies, explicitly listing WGCNA's Bioc dependencies for robustness)
pkgs_bioc_and_wgcna <- c("DESeq2", "org.Hs.eg.db", "impute", "preprocessCore", "GO.db", "AnnotationDbi", "WGCNA")
for (p in pkgs_bioc_and_wgcna) {
    if (!requireNamespace(p, quietly = TRUE)) {
        BiocManager::install(p, update = FALSE, ask = FALSE)
    }
}

# Other CRAN packages
pkgs_cran_others <- c("tidyverse", "pheatmap", "reshape2", "ggplot2") # Removed flashClust
for (p in pkgs_cran_others) {
    if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p)
    }
}

library(WGCNA)
library(DESeq2)
library(tidyverse)
library(pheatmap)
# library(flashClust) # Removed flashClust
library(reshape2)
library(ggplot2)

# Important for WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ====  Loading and Preparing Data ====

# Load normalized count data
count_data <- read.csv("filtered_count_data.csv", row.names = 1)


# Load sample metadata
col_data <- read.csv("col_data.csv")
# Ensure Sample column matches count_data column names
rownames(col_data) <- col_data$Sample

# Transpose the data: WGCNA expects samples as rows, genes as columns
datExpr <- t(count_data)


# Checking for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
    # Removing offending genes and samples
    if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples) > 0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ==== Sample Clustering to Detect Outliers ====

# Cluster samples to see if there are obvious outliers
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample clustering
pdf("1_sample_clustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# ==== Preparing Trait Data ====

# Convert conditions to binary trait matrix for correlation analysis
# Create binary variables for each condition
trait_data <- data.frame(
    Normal = as.numeric(col_data$Conditions == "Normal Human Breast Organoids"),
    HER2 = as.numeric(col_data$Conditions == "HER2 Positive Breast Tumor"),
    NonTNBC = as.numeric(col_data$Conditions == "Non-TNBC Breast Tumor"),
    TNBC = as.numeric(col_data$Conditions == "TNBC Breast Tumor"),
    row.names = col_data$Sample
)

# Ensure trait data matches datExpr sample order
trait_data <- trait_data[rownames(datExpr), ]

# ==== Choosing Soft-Thresholding Power (β) ====

# Choosing a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
pdf("2_scale_independence.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")  # Threshold line at R^2 = 0.80

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Choose power based on scale-free topology criterion
# Typically choose the lowest power where R^2 > 0.80
softPower <- sft$powerEstimate
if (is.na(softPower)) {
    # If automatic selection fails, choose manually (e.g., 6 or 8)
    softPower <- 6
    cat("Automatic power selection failed. Using softPower =", softPower, "\n")
} else {
    cat("Chosen soft-thresholding power:", softPower, "\n")
}

# ==== Construct Network and Identify Modules ====

net <- blockwiseModules(
  datExpr,
  power = softPower,
  corType = "bicor",
  maxPOutliers = 0.1,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)

# Module assignments
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
table(moduleColors)



MEs <- orderMEs(net$MEs)

moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

sig_threshold <- 0.05
get_sig_modules <- function(trait) {
  p <- moduleTraitPvalue[, trait]
  hit <- !is.na(p) & p < sig_threshold
  gsub("^ME", "", rownames(moduleTraitPvalue)[hit])
}

her2_modules    <- get_sig_modules("HER2")
tnbc_modules    <- get_sig_modules("TNBC")
nontnbc_modules <- get_sig_modules("NonTNBC")

her2_modules
tnbc_modules
nontnbc_modules

# ---- Align metadata samples to expression samples (DO THIS ONCE) ----
expr_ids_clean <- sub("^X", "", rownames(datExpr))
expr_ids_clean <- gsub("\\.", "-", expr_ids_clean)

meta_ids_clean <- sub("^X", "", as.character(col_data$Sample))
meta_ids_clean <- gsub("\\.", "-", meta_ids_clean)

match_idx <- match(expr_ids_clean, meta_ids_clean)
stopifnot(sum(!is.na(match_idx)) == nrow(datExpr))

col_data2 <- col_data[match_idx, , drop = FALSE]
rownames(col_data2) <- rownames(datExpr)

# ---- Build trait matrix using *aligned* metadata ----
trait_data <- data.frame(
  HER2    = as.numeric(col_data2$Conditions == "HER2 Positive Breast Tumor"),
  TNBC    = as.numeric(col_data2$Conditions == "TNBC Breast Tumor"),
  NonTNBC = as.numeric(col_data2$Conditions == "Non-TNBC Breast Tumor"),
  Normal  = as.numeric(col_data2$Conditions == "Normal Human Breast Organoids"),
  row.names = rownames(datExpr)
)

print(colSums(trait_data))          # should be 5,5,6,3
stopifnot(all(colSums(trait_data) > 0))

# ==== Visualising Gene Dendrogram and Module Colors ====

# Plot the dendrogram and module colors
pdf("3_gene_dendrogram_modules.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

# ==== Relating Modules to Traits ====

# Calculate module eigengenes (ME)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# Order MEs by hierarchical clustering
MEs <- orderMEs(MEs)

# Calculate correlation between MEs and traits
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Save results
write.csv(moduleTraitCor, "module_trait_correlation.csv")
write.csv(moduleTraitPvalue, "module_trait_pvalue.csv")

# ==== Visualize Module-Trait Relationships ====

# Create a heatmap of module-trait correlations
pdf("4_module_trait_heatmap.pdf", width = 8, height = 10)
# Display correlation values and p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_data),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()

# ==== Identifying Modules of Interest ====

# Define significance threshold
sig_threshold <- 0.05

# Find modules significantly correlated with HER2
her2_modules <- rownames(moduleTraitPvalue)[moduleTraitPvalue[, "HER2"] < sig_threshold]
cat("\nModules significantly correlated with HER2:\n")
print(her2_modules)

# Find modules significantly correlated with TNBC
tnbc_modules <- rownames(moduleTraitPvalue)[moduleTraitPvalue[, "TNBC"] < sig_threshold]
cat("\nModules significantly correlated with TNBC:\n")
print(tnbc_modules)

# Find modules significantly correlated with NonTNBC
nontnbc_modules <- rownames(moduleTraitPvalue)[moduleTraitPvalue[, "NonTNBC"] < sig_threshold]
cat("\nModules significantly correlated with NonTNBC:\n")
print(nontnbc_modules)

 ====
### Exporting Gene Lists for Significant Modules
 ====

# Get gene names for each module
genes <- colnames(datExpr)
module_gene_list <- data.frame(
    gene_id = genes,
    module_color = moduleColors,
    stringsAsFactors = FALSE
)

# Save complete module assignments
write.csv(module_gene_list, "module_gene_assignments.csv", row.names = FALSE)

# Export genes for each significant module
all_sig_modules <- unique(c(her2_modules, tnbc_modules, nontnbc_modules))

for (module in all_sig_modules) {
    module_name <- gsub("ME", "", module)
    module_genes <- genes[moduleColors == module_name]
    write.table(module_genes,
                file = paste0("genes_module_", module_name, ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("Module", module_name, "contains", length(module_genes), "genes\n")
}

 ====
### HUB Gene Extraction ====

# Collecting ALL significant modules across traits
sig_MEs <- unique(c(her2_modules, tnbc_modules, nontnbc_modules))
sig_MEs <- unique(na.omit(sig_MEs))

# Remove grey if present
sig_MEs <- sig_MEs[sig_MEs != "ME0"]

print(sig_MEs)  # sanity check: should be like "ME1" "ME4" "ME6" ...

#kME matrix (genes x modules)
kME <- as.data.frame(cor(datExpr, net$MEs, use = "p"))
colnames(kME) <- colnames(net$MEs)

# Mapping numeric labels -> colours (for filenames / readability)
label_to_color <- setNames(
  WGCNA::labels2colors(sort(unique(moduleLabels))),
  sort(unique(moduleLabels))
)

# Looping through significant modules
for (MEname in sig_MEs) {

  # Extracting numeric module label (ME1 -> 1)
  lab <- suppressWarnings(as.integer(sub("^ME", "", MEname)))
  if (is.na(lab) || lab == 0) next

  inModule <- moduleLabels == lab
  if (!any(inModule)) next

  # kME for genes in this module
  kME_mod <- kME[inModule, MEname, drop = FALSE]
  kME_mod <- kME_mod[order(-abs(kME_mod[[MEname]])), , drop = FALSE]

  top <- head(kME_mod, 20)

  module_color <- label_to_color[as.character(lab)]

  write.csv(
    data.frame(
      gene = rownames(top),
      kME  = top[[MEname]]
    ),
    file = paste0("hub_genes_", MEname, "_", module_color, ".csv"),
    row.names = FALSE
  )

  cat("\nTop 10 hub genes for", MEname, "(", module_color, "):\n")
  print(head(rownames(top), 10))
}

# ==== Converting Colour-based significant modules to Numeric ME labels ====

# Build colour -> numeric label mapping from the network
color_to_label <- tapply(moduleLabels, moduleColors, function(x) {
  as.integer(names(sort(table(x), decreasing = TRUE))[1])
})

cat("\nColour → numeric label mapping:\n")
print(color_to_label)

convert_color_MEs_to_numeric <- function(MEs_color) {
  cols <- gsub("^ME", "", MEs_color)
  labs <- color_to_label[cols]
  labs <- labs[!is.na(labs) & labs != 0]   # drop NA and grey
  paste0("ME", labs)
}

# Convert significant module lists
her2_modules_num    <- convert_color_MEs_to_numeric(her2_modules)
tnbc_modules_num    <- convert_color_MEs_to_numeric(tnbc_modules)
nontnbc_modules_num <- convert_color_MEs_to_numeric(nontnbc_modules)

# Combine all significant numeric MEs
sig_MEs_numeric <- unique(c(her2_modules_num,
                            tnbc_modules_num,
                            nontnbc_modules_num))

cat("\nSignificant numeric ME labels:\n")
print(sig_MEs_numeric)

# ==== Hub Gene Extraction Using Numeric ME Labels ====

# Use the numeric ME labels produced in the last Step
sig_MEs <- sig_MEs_numeric
print(sig_MEs)

# kME matrix: genes x module eigengenes
kME <- as.data.frame(cor(datExpr, net$MEs, use = "p"))
colnames(kME) <- colnames(net$MEs)

# Map numeric labels to colours (for filenames / readability)
label_to_color <- setNames(
  WGCNA::labels2colors(sort(unique(moduleLabels))),
  sort(unique(moduleLabels))
)

for (MEname in sig_MEs) {

  # Extract numeric label (ME4 -> 4)
  lab <- suppressWarnings(as.integer(sub("^ME", "", MEname)))
  if (is.na(lab) || lab == 0) next

  inModule <- moduleLabels == lab
  if (!any(inModule)) next

  # kME for genes in this module
  kME_mod <- kME[inModule, MEname, drop = FALSE]
  kME_mod <- kME_mod[order(-abs(kME_mod[[MEname]])), , drop = FALSE]

  top <- head(kME_mod, 20)

  module_color <- label_to_color[as.character(lab)]

  out_file <- paste0("hub_genes_", MEname, "_", module_color, ".csv")

  write.csv(
    data.frame(
      gene = rownames(top),
      kME  = top[[MEname]]
    ),
    file = out_file,
    row.names = FALSE
  )

  cat("\nTop 10 hub genes for", MEname, "(", module_color, "):\n", sep = "")
  print(head(rownames(top), 10))
  cat("Saved:", out_file, "\n")
}

# ==== # Reassign canonical module lists to numeric ME labels ====

her2_modules    <- her2_MEs_numeric
tnbc_modules    <- tnbc_MEs_numeric
nontnbc_modules <- nontnbc_MEs_numeric

# Sanity check
cat("HER2 numeric modules:\n")
print(her2_modules)

cat("TNBC numeric modules:\n")
print(tnbc_modules)

cat("NonTNBC numeric modules:\n")
print(nontnbc_modules)

# ==== Summary Statistics ====

# ===== FINAL WGCNA SUMMARY (REPORTING ONLY) =====

label_to_color <- setNames(
  WGCNA::labels2colors(sort(unique(moduleLabels))),
  sort(unique(moduleLabels))
)

clean_modules <- function(x) {
  x <- unique(na.omit(x))
  x <- sub("^ME", "", x)
  x <- x[x != "0"]
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) "None" else paste(label_to_color[x], collapse = ", ")
}

cat("Total genes analyzed:", ncol(datExpr), "\n")
cat("Total samples:", nrow(datExpr), "\n")
cat("Soft-thresholding power (β):", softPower, "\n")
cat("Number of modules identified (excluding grey):",
    length(setdiff(unique(moduleLabels), 0)), "\n")

cat("\nModule sizes:\n")
print(table(moduleColors))

cat("\nModules significantly associated with traits (p < 0.05):\n")
cat("  HER2:", clean_modules(her2_modules), "\n")
cat("  TNBC:", clean_modules(tnbc_modules), "\n")
cat("  NonTNBC:", clean_modules(nontnbc_modules), "\n")

# ==== Saving R workspace for future use ====

save.image("WGCNA_analysis_final.RData")
