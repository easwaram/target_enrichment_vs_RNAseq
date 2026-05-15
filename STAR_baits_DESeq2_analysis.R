# STAR_baits_DESeq2_analysis.R
# DESeq2 differential expression: male vs female (liver)
# Input: STAR + featureCounts output for library and baits-enriched samples
# Organism: Danio rerio (zebrafish)

library(DESeq2)
library(edgeR)
library(ggplot2)
library(apeglm)

# --- Load STAR featureCounts output ---

counts_lib_STAR  <- read.delim("~/Documents/STAR_featureCounts/STAR_library_liver_gene_counts.txt",
                                comment.char = "#", check.names = FALSE)
counts_baits_STAR <- read.delim("~/Documents/STAR_featureCounts/STAR_baits_liver_gene_counts.txt",
                                 comment.char = "#", check.names = FALSE)

# Extract count matrix (columns 7+ are samples; first 6 are featureCounts annotation)
rownames(counts_lib_STAR)  <- counts_lib_STAR$Geneid
rownames(counts_baits_STAR) <- counts_baits_STAR$Geneid

lib_matrix_STAR  <- as.matrix(counts_lib_STAR[,  7:ncol(counts_lib_STAR)])
baits_matrix_STAR <- as.matrix(counts_baits_STAR[, 7:ncol(counts_baits_STAR)])

# Clean .bam suffix from column names
colnames(lib_matrix_STAR)  <- gsub("\\.bam$", "", colnames(lib_matrix_STAR))
colnames(baits_matrix_STAR) <- gsub("\\.bam$", "", colnames(baits_matrix_STAR))

# --- CPM filtering ---

dge_lib_STAR  <- DGEList(counts = lib_matrix_STAR)
dge_baits_STAR <- DGEList(counts = baits_matrix_STAR)

cpm_lib_STAR  <- cpm(dge_lib_STAR)
cpm_baits_STAR <- cpm(dge_baits_STAR)

# Library: CPM > 0.33 in at least 3 samples (~10 reads at 30M depth)
keep_lib  <- rowSums(cpm_lib_STAR  > 0.33) >= 3
# Baits: slightly higher CPM threshold due to lower sequencing depth
keep_baits <- rowSums(cpm_baits_STAR > 0.50) >= 3

counts_filtered_lib_STAR  <- dge_lib_STAR$counts[keep_lib, ]
counts_filtered_baits_STAR <- dge_baits_STAR$counts[keep_baits, ]

cat("Genes retained (library):", nrow(counts_filtered_lib_STAR), "\n")
cat("Genes retained (baits):",   nrow(counts_filtered_baits_STAR), "\n")

# --- Defensome gene subsetting ---

defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv",
                               header = FALSE, stringsAsFactors = FALSE)[, 1]

common_lib  <- intersect(rownames(counts_filtered_lib_STAR),  defensome_genes)
common_baits <- intersect(rownames(counts_filtered_baits_STAR), defensome_genes)

baits_defensome_counts <- counts_filtered_baits_STAR[common_baits, , drop = FALSE]
lib_defensome_counts   <- counts_filtered_lib_STAR[common_lib,   , drop = FALSE]

cat("Defensome genes in library dataset:", length(common_lib),  "\n")
cat("Defensome genes in baits dataset:",   length(common_baits), "\n")

# --- Sample metadata ---

samples_lib  <- c("L1_S1_L001", "L2_S2_L001", "L3_S3_L001",
                  "L4_S4_L001", "L5_S5_L001", "L6_S6_L001")
samples_baits <- c("E1_S1_L001", "E2_S2_L001", "E3_S3_L001",
                   "E4_S4_L001", "E5_S5_L001", "E6_S6_L001")

sex_design_lib <- data.frame(
  sample     = samples_lib,
  background = c("male", "male", "male", "female", "female", "female")
)

sex_design_baits <- data.frame(
  sample     = samples_baits,
  background = c("male", "male", "male", "female", "female", "female")
)

# --- DESeq2 objects ---

dds_lib_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered_lib_STAR,
  colData   = sex_design_lib,
  design    = ~background
)

dds_baits_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered_baits_STAR,
  colData   = sex_design_baits,
  design    = ~background
)

# Defensome-only objects (for PCA visualization)
dds_lib_defensome <- DESeqDataSetFromMatrix(
  countData = lib_defensome_counts,
  colData   = sex_design_lib,
  design    = ~background
)

dds_baits_defensome <- DESeqDataSetFromMatrix(
  countData = baits_defensome_counts,
  colData   = sex_design_baits,
  design    = ~background
)

# --- PCA ---

# Full transcriptome PCA
plotPCA(rlog(dds_lib_all,  blind = TRUE), intgroup = "background", ntop = 10000)
plotPCA(rlog(dds_baits_all, blind = TRUE), intgroup = "background", ntop = 1000)

# Defensome subset PCA
plotPCA(rlog(dds_lib_defensome,  blind = TRUE), intgroup = "background", ntop = 400)
plotPCA(rlog(dds_baits_defensome, blind = TRUE), intgroup = "background", ntop = 400)

# --- DESeq2 and LFC shrinkage ---

dds_lib2  <- DESeq(dds_lib_all)
dds_baits2 <- DESeq(dds_baits_all)

plotDispEsts(dds_lib2,  legend = FALSE)
plotDispEsts(dds_baits2, legend = TRUE)

# LFC shrinkage with apeglm
res_lib_shrunk  <- lfcShrink(dds_lib2,  coef = "background_male_vs_female", type = "apeglm")
res_baits_shrunk <- lfcShrink(dds_baits2, coef = "background_male_vs_female", type = "apeglm")

# Subset to defensome genes
res_lib_defensome  <- na.omit(res_lib_shrunk[rownames(res_lib_shrunk)   %in% defensome_genes, ])
res_baits_defensome <- na.omit(res_baits_shrunk[rownames(res_baits_shrunk) %in% defensome_genes, ])

summary(res_lib_defensome)
summary(res_baits_defensome)

# --- Compare LFC between library and baits (Spearman correlation) ---

res_lib_df  <- as.data.frame(res_lib_shrunk)
res_baits_df <- as.data.frame(res_baits_shrunk)

res_lib_df$gene  <- rownames(res_lib_df)
res_baits_df$gene <- rownames(res_baits_df)

res_lib_df  <- res_lib_df[res_lib_df$gene   %in% defensome_genes, ]
res_baits_df <- res_baits_df[res_baits_df$gene %in% defensome_genes, ]

merged_res <- merge(
  res_lib_df[,  c("gene", "log2FoldChange")],
  res_baits_df[, c("gene", "log2FoldChange")],
  by = "gene", suffixes = c("_lib", "_baits")
)

# Scatter plot: library LFC vs baits LFC
ggplot(merged_res, aes(x = log2FoldChange_lib, y = log2FoldChange_baits, label = gene)) +
  geom_point(color = "darkgreen", size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  labs(
    x = "Library log2FC (male vs female, shrunk)",
    y = "Baits log2FC (male vs female, shrunk)",
    title = "Defensome Genes: LFC Correlation - Library vs Baits"
  ) +
  theme_minimal()

# Spearman correlation
cor_test <- cor.test(merged_res$log2FoldChange_lib, merged_res$log2FoldChange_baits,
                     method = "spearman")
cat("Spearman r:", round(cor_test$estimate, 3),
    " | p-value:", signif(cor_test$p.value, 3), "\n")
