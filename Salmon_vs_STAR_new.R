BiocManager::install("tximport")
install.packages("tidyverse")

# Load packages
library(tximport)
library(tidyverse)

# Locate and name your quant.sf files
lib_dir <- list.dirs("~/Documents/library_salmon_counts_redo", recursive = FALSE)
baits_dir <- list.dirs("~/Documents/enriched_salmon_counts_redo", recursive = FALSE)

lib_quant_files <- file.path(lib_dir, "quant.sf")
baits_quant_files <- file.path(baits_dir, "quant.sf")

# use folder names as sample names
names(lib_quant_files) <- basename(lib_dir)  
names(baits_quant_files) <- basename(baits_dir)

file.exists(lib_quant_files)
file.exists(baits_quant_files)

### Create tx2gene file using gtf file
library(rtracklayer)
library(dplyr)
library(stringr)

# Import the GTF
gtf <- import("~/Documents/GRCz12tu_genomic.gtf")

# Extract metadata as a data frame
gtf_df <- as.data.frame(mcols(gtf))

# Filter to rows where transcript_id is present and not empty
tx2gene <- gtf_df %>%
  filter(!is.na(transcript_id) & transcript_id != "" & !is.na(gene_id)) %>%
  select(transcript_id, gene_id) %>%
  distinct() %>%
  rename(TXNAME = transcript_id, GENEID = gene_id)

# Strip version numbers (e.g., XR_012404707.1 → XR_012404707) to match Salmon
tx2gene$TXNAME <- sub("\\.\\d+$", "", tx2gene$TXNAME)

# Inspect result
head(tx2gene)

# Save as a TSV (recommended for bioinformatics workflows)
#write.table(tx2gene, file = "tx2gene_GRCz12tu.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Load tx2gene file
tx2gene <- read.delim("tx2gene_GRCz12tu.tsv", stringsAsFactors = FALSE)

### Import counts from quant.sf using tximport
txi_lib <- tximport(lib_quant_files, 
                    type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_baits <- tximport(baits_quant_files, 
                      type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

### Filter out zero and low counts (~10 reads)

# Load edgeR
library(edgeR)

# Create a DGEList object (for easy normalization)
# Step 1: Start with raw counts
dge_lib_salmon <- DGEList(counts = txi_lib$counts)
dge_baits_salmon <- DGEList(counts = txi_baits$counts)

# Calculate CPM
cpm_matrix_lib <- cpm(dge_lib_salmon)
cpm_matrix_baits <- cpm(dge_baits_salmon)

head(cpm_matrix_lib)
head(cpm_matrix_baits)

# This filters out those with a CPM of greater than 0.33 in atleast 3 samples (so about 10 reads according to library size of 30M reads)
# Changed the CPM cut-off to 0.5 for baits since much lower depth
# This then filters those out and keeps the counts_filtered as raw counts to use for DESeq or rarefaction (both needs raw counts; rarefaction needs rounded)
keep_lib_salmon <- rowSums(cpm_matrix_lib > 0.33) >= 3
keep_baits_salmon <- rowSums(cpm_matrix_baits > 0.5) >= 3

# Filtered counts
counts_filtered_lib <- dge_lib_salmon$counts[keep_lib_salmon, ]
counts_filtered_baits <- dge_baits_salmon$counts[keep_baits_salmon, ]

colSums(counts_filtered_lib)
colSums(counts_filtered_baits)

# Re-create DGEList objects with filtered counts
dge_lib_filtered <- DGEList(counts = counts_filtered_lib)
dge_baits_filtered <- DGEList(counts = counts_filtered_baits)

# Calculate CPMs after filtering
cpm_filtered_lib <- cpm(dge_lib_filtered)
cpm_filtered_baits <- cpm(dge_baits_filtered)


# Load STAR featurecounts results
# Read STAR featureCounts output (skip first row of comments)
counts_baits_STAR <- read.delim("~/Documents/STAR_featureCounts/STAR_baits_liver_gene_counts.txt", comment.char = "#", check.names = FALSE)
counts_lib_STAR <- read.delim("~/Documents/STAR_featureCounts/STAR_library_liver_gene_counts.txt", comment.char = "#", check.names = FALSE)

# Drop annotation columns (first 6 columns) and keep count matrix only
rownames(counts_baits_STAR) <- counts_baits_STAR$Geneid
baits_matrix_STAR <- as.matrix(counts_baits_STAR[, 7:ncol(counts_baits_STAR)])

rownames(counts_lib_STAR) <- counts_lib_STAR$Geneid
lib_matrix_STAR <- as.matrix(counts_lib_STAR[, 7:ncol(counts_lib_STAR)])

# Convert STAR counts to CPM
library(edgeR)

dge_baits_STAR <- DGEList(counts = baits_matrix_STAR)
cpm_baits_STAR <- cpm(dge_baits_STAR)

dge_lib_STAR   <- DGEList(counts = lib_matrix_STAR)
cpm_lib_STAR <- cpm(dge_lib_STAR)

head (cpm_baits_STAR)
head (cpm_lib_STAR)

# Filter for baits: keep genes with CPM > 0.5 in at least 3 samples
keep_baits_STAR <- rowSums(cpm_baits_STAR > 0.25) >= 3
counts_filtered_baits_STAR <- dge_baits_STAR$counts[keep_baits_STAR, ]
dge_baits_filtered_STAR <- DGEList(counts = counts_filtered_baits_STAR) # Create new filtered DGEList objects
cpm_filtered_baits_STAR <- cpm(dge_baits_filtered_STAR) # Calculate CPMs on filtered data

# Filter for libraries: keep genes with CPM > 0.33 in at least 3 samples
keep_lib_STAR <- rowSums(cpm_lib_STAR > 0.165) >= 3
counts_filtered_lib_STAR   <- dge_lib_STAR$counts[keep_lib_STAR, ]
dge_lib_filtered_STAR   <- DGEList(counts = counts_filtered_lib_STAR)
cpm_filtered_lib_STAR   <- cpm(dge_lib_filtered_STAR)


# Now I want to compare the CPM between Salmon counts and STAR counts

# Read defensome gene list (assuming it's one gene symbol per line)
defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv", header = FALSE, stringsAsFactors = FALSE)[,1]

# Assuming cpm_filtered_baits (Salmon) and cpm_filtered_baits_STAR (STAR) are already loaded

# Keep genes that exist in both matrices
common_defensome_genes_cpm <- intersect(
  intersect(rownames(cpm_filtered_baits), rownames(cpm_filtered_baits_STAR)),
  defensome_genes
)

# Subset CPM matrices
salmon_defensome_cpm <- cpm_filtered_baits[common_defensome_genes_cpm, ]
star_defensome_cpm   <- cpm_filtered_baits_STAR[common_defensome_genes_cpm, ]

library(tidyr)
library(dplyr)

# Clean column names in Salmon CPM
colnames(salmon_defensome_cpm) <- gsub("_quant$", "", colnames(salmon_defensome_cpm))

# Clean column names in STAR CPM
colnames(star_defensome_cpm) <- gsub("\\.bam$", "", colnames(star_defensome_cpm))

male_samples   <- c("E1_S1_L001", "E2_S2_L001", "E3_S3_L001")
female_samples <- c("E4_S4_L001", "E5_S5_L001", "E6_S6_L001")

# Calculate average CPM across males/females
salmon_male_avg   <- rowMeans(salmon_defensome_cpm[, male_samples, drop = FALSE])
salmon_female_avg <- rowMeans(salmon_defensome_cpm[, female_samples, drop = FALSE])

star_male_avg     <- rowMeans(star_defensome_cpm[, male_samples, drop = FALSE])
star_female_avg   <- rowMeans(star_defensome_cpm[, female_samples, drop = FALSE])

# Step 2: Log-transform
log_salmon_male   <- log10(salmon_male_avg + 1)
log_star_male     <- log10(star_male_avg + 1)

log_salmon_female <- log10(salmon_female_avg + 1)
log_star_female   <- log10(star_female_avg + 1)

# Step 3: Define outliers
threshold <- 1
outlier_male   <- abs(log_star_male - log_salmon_male) > threshold
outlier_female <- abs(log_star_female - log_salmon_female) > threshold

# Step 4: Plot with color coding
par(mfrow = c(1, 2))

# Male
plot(log_star_male, log_salmon_male,
     xlab = "STAR Male CPM (log10)", ylab = "Salmon Male CPM (log10)",
     main = "Male: STAR vs Salmon (baits)",
     col = ifelse(outlier_male, "black", "blue"),
     pch = 16)
abline(a = 0, b = 1, lty = 2, col = "gray")

# Female
plot(log_star_female, log_salmon_female,
     xlab = "STAR Female CPM (log10)", ylab = "Salmon Female CPM (log10)",
     main = "Female: STAR vs Salmon (baits)",
     col = ifelse(outlier_female, "black", "red"),
     pch = 16)
abline(a = 0, b = 1, lty = 2, col = "gray")

par(mfrow = c(1, 1))  # Reset

# Determine outliers from the Salmon vs STAR comparison
outlier_genes_male   <- names(log_star_male)[outlier_male]
outlier_genes_female <- names(log_star_female)[outlier_female]

cat("Outlier genes in males:", length(outlier_genes_male), "\n")
print(outlier_genes_male)

cat("Outlier genes in females:", length(outlier_genes_female), "\n")
print(outlier_genes_female)


### Now Im going to redo the above but now for the library samples

# Read defensome gene list (assuming one gene symbol per line)
defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv", 
                              header = FALSE, stringsAsFactors = FALSE)[, 1]

# Get common defensome genes present in both STAR and Salmon CPM matrices
common_defensome_genes_cpm <- intersect(
  intersect(rownames(cpm_filtered_lib), rownames(cpm_filtered_lib_STAR)),
  defensome_genes
)

# Subset CPM matrices for common defensome genes
salmon_defensome_cpm_lib <- cpm_filtered_lib[common_defensome_genes_cpm, ]
star_defensome_cpm_lib   <- cpm_filtered_lib_STAR[common_defensome_genes_cpm, ]

# Clean column names to remove _quant and .bam suffixes for alignment
colnames(salmon_defensome_cpm_lib) <- gsub("_quant$", "", colnames(salmon_defensome_cpm_lib))
colnames(star_defensome_cpm_lib)   <- gsub("\\.bam$", "", colnames(star_defensome_cpm_lib))

# Define sample groupings
male_samples   <- c("L1_S1_L001", "L2_S2_L001", "L3_S3_L001")
female_samples <- c("L4_S4_L001", "L5_S5_L001", "L6_S6_L001")

# Calculate average CPM across male and female samples
salmon_male_avg   <- rowMeans(salmon_defensome_cpm_lib[, male_samples, drop = FALSE])
salmon_female_avg <- rowMeans(salmon_defensome_cpm_lib[, female_samples, drop = FALSE])

star_male_avg     <- rowMeans(star_defensome_cpm_lib[, male_samples, drop = FALSE])
star_female_avg   <- rowMeans(star_defensome_cpm_lib[, female_samples, drop = FALSE])

# Log-transform to avoid skew from large counts
log_salmon_male   <- log10(salmon_male_avg + 1)
log_star_male     <- log10(star_male_avg + 1)

log_salmon_female <- log10(salmon_female_avg + 1)
log_star_female   <- log10(star_female_avg + 1)

# Identify outliers where absolute log-difference > threshold
threshold <- 1
outlier_male   <- abs(log_star_male - log_salmon_male) > threshold
outlier_female <- abs(log_star_female - log_salmon_female) > threshold

# Plot male vs female separately, coloring outliers black
par(mfrow = c(1, 2))

# Male plot
plot(log_star_male, log_salmon_male,
     xlab = "STAR Male CPM (log10)", ylab = "Salmon Male CPM (log10)",
     main = "Library: STAR vs Salmon (Males)",
     col = ifelse(outlier_male, "black", "blue"),
     pch = 16)
abline(a = 0, b = 1, col = "gray", lty = 2)

# Female plot
plot(log_star_female, log_salmon_female,
     xlab = "STAR Female CPM (log10)", ylab = "Salmon Female CPM (log10)",
     main = "Library: STAR vs Salmon (Females)",
     col = ifelse(outlier_female, "black", "red"),
     pch = 16)
abline(a = 0, b = 1, col = "gray", lty = 2)

par(mfrow = c(1, 1))  # Reset

# Report outlier genes
outlier_genes_male   <- names(log_star_male)[outlier_male]
outlier_genes_female <- names(log_star_female)[outlier_female]

cat("Outlier genes in library (males):", length(outlier_genes_male), "\n")
print(outlier_genes_male)

cat("\nOutlier genes in library (females):", length(outlier_genes_female), "\n")
print(outlier_genes_female)


#### Now comparing CPM of library vs baits for STAR analysis
# Read defensome gene list
defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv", 
                              header = FALSE, stringsAsFactors = FALSE)[, 1]

# Get common defensome genes across both datasets
common_defensome_genes <- Reduce(intersect, list(
  rownames(cpm_filtered_baits_STAR),
  rownames(cpm_filtered_lib_STAR),
  defensome_genes
))

# Subset CPM matrices
baits_defensome_cpm <- cpm_filtered_baits_STAR[common_defensome_genes, ]
lib_defensome_cpm   <- cpm_filtered_lib_STAR[common_defensome_genes, ]

# Clean column names
colnames(baits_defensome_cpm) <- gsub("\\.bam$", "", colnames(baits_defensome_cpm))
colnames(lib_defensome_cpm)   <- gsub("\\.bam$", "", colnames(lib_defensome_cpm))

# Define sample groups
baits_male_samples   <- c("E1_S1_L001", "E2_S2_L001", "E3_S3_L001")
baits_female_samples <- c("E4_S4_L001", "E5_S5_L001", "E6_S6_L001")

lib_male_samples     <- c("L1_S1_L001", "L2_S2_L001", "L3_S3_L001")
lib_female_samples   <- c("L4_S4_L001", "L5_S5_L001", "L6_S6_L001")

# Calculate average CPM for each group
baits_male_avg   <- rowMeans(baits_defensome_cpm[, baits_male_samples, drop = FALSE])
baits_female_avg <- rowMeans(baits_defensome_cpm[, baits_female_samples, drop = FALSE])
lib_male_avg     <- rowMeans(lib_defensome_cpm[, lib_male_samples, drop = FALSE])
lib_female_avg   <- rowMeans(lib_defensome_cpm[, lib_female_samples, drop = FALSE])

# Plot comparison of male and female expression: enriched (Y) vs library (X)
par(mfrow = c(1, 2))

# Male plot
plot(lib_male_avg, baits_male_avg,
     log = "xy", pch = 16, col = "blue",
     xlab = "Library Male CPM (avg)", ylab = "Baits Male CPM (avg)",
     main = "Male: Baits vs Library (CPM)")
abline(a = 0, b = 1, col = "gray", lty = 2)

# Female plot
plot(lib_female_avg, baits_female_avg,
     log = "xy", pch = 16, col = "red",
     xlab = "Library Female CPM (avg)", ylab = "Baits Female CPM (avg)",
     main = "Female: Baits vs Library (CPM)")
abline(a = 0, b = 1, col = "gray", lty = 2)

par(mfrow = c(1, 1))  # Reset layout

# Log10-transform (after averaging CPMs)
#log_lib_male    <- log10(lib_male_avg + 1)
#log_baits_male  <- log10(baits_male_avg + 1)

#log_lib_female  <- log10(lib_female_avg + 1)
#log_baits_female <- log10(baits_female_avg + 1)

# Linear model: Enriched ~ Library
#lm_male    <- lm(log_baits_male ~ log_lib_male)
#lm_female  <- lm(log_baits_female ~ log_lib_female)

# R² values
#r2_male   <- summary(lm_male)$r.squared
#r2_female <- summary(lm_female)$r.squared

#cat("R² (males):", round(r2_male, 3), "\n")
#cat("R² (females):", round(r2_female, 3), "\n")

# Spearman correlation
cor_male   <- cor.test(lib_male_avg, baits_male_avg, method = "spearman")
cor_female <- cor.test(lib_female_avg, baits_female_avg, method = "spearman")

cat("Spearman correlation (males):", round(cor_male$estimate, 3), 
    "p-value:", signif(cor_male$p.value, 3), "\n")

cat("Spearman correlation (females):", round(cor_female$estimate, 3), 
    "p-value:", signif(cor_female$p.value, 3), "\n")

cor.test(1:100, 1:100, method = "spearman")

# Shapiro- Wilk Normality test
# Apply Shapiro-Wilk normality test on each log-transformed vector
#shapiro_lib_male    <- shapiro.test(log_lib_male)
#shapiro_baits_male  <- shapiro.test(log_baits_male)

#shapiro_lib_female  <- shapiro.test(log_lib_female)
#shapiro_baits_female <- shapiro.test(log_baits_female)

# Print results
#cat("Shapiro-Wilk for log_lib_male:\n"); print(shapiro_lib_male)
#cat("\nShapiro-Wilk for log_baits_male:\n"); print(shapiro_baits_male)

#cat("\nShapiro-Wilk for log_lib_female:\n"); print(shapiro_lib_female)
#cat("\nShapiro-Wilk for log_baits_female:\n"); print(shapiro_baits_female)

# If I wanna end up going with Pearson's correlation test
# Perform Pearson correlation tests
#pearson_male <- cor.test(log_lib_male, log_baits_male, method = "pearson")
#pearson_female <- cor.test(log_lib_female, log_baits_female, method = "pearson")

# Print the results
#cat("Pearson Correlation (Males):\n")
#print(pearson_male)

#cat("\nPearson Correlation (Females):\n")
#print(pearson_female)

