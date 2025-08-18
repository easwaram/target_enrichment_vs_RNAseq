BiocManager::install("tximport")
install.packages("tidyverse")

# Load packages
library(tximport)
library(tidyverse)

# Locate and name your quant.sf files
lib_dir <- list.dirs("~/Documents/library_salmon_counts", recursive = FALSE)
baits_dir <- list.dirs("~/Documents/enriched_salmon_counts", recursive = FALSE)

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

# Load the NCBI GTF
#gtf <- import("~/Downloads/genomic.gtf")

# Keep only features with transcript_id
#gtf_tx <- gtf[!is.na(mcols(gtf)$transcript_id)]

# Extract transcript_id and db_xref field
#txinfo <- mcols(gtf_tx)[, c("transcript_id", "db_xref")]

# Flatten list-columns into text
#txinfo <- as.data.frame(txinfo)

# Extract GeneID from db_xref field
#txinfo$GeneID <- str_extract(txinfo$db_xref, "GeneID:\\d+") %>%
#  str_remove("GeneID:")

# Keep only valid rows not missing data; removing duplicates
#tx2gene <- txinfo %>%
#  filter(!is.na(GeneID)) %>%
#  distinct(transcript_id, GeneID) %>%
#  rename(TXNAME = transcript_id, GENEID = GeneID)

# Strip version numbers to match Salmon
#tx2gene$TXNAME <- sub("\\.\\d+$", "", tx2gene$TXNAME)
#head(tx2gene)

# Save as a TSV (recommended for bioinformatics workflows)
#write.table(tx2gene, file = "tx2gene.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

tx2gene <- read.delim("tx2gene.tsv", stringsAsFactors = FALSE)

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
dge_lib <- DGEList(counts = txi_lib$counts)
dge_baits <- DGEList(counts = txi_baits$counts)

# Calculate CPM
cpm_matrix_lib <- cpm(dge_lib)
cpm_matrix_baits <- cpm(dge_baits)

head(cpm_matrix_lib)
head(cpm_matrix_baits)

# This filters out those with a CPM of greater than 0.33 in atleast 3 samples (so about 10 reads according to library size of 30M reads)
# Changed the CPM cut-off to 0.5 for baits since much lower depth
# This then filters those out and keeps the counts_filtered as raw counts to use for DESeq or rarefaction (both needs raw counts; rarefaction needs rounded)
keep_lib <- rowSums(cpm_matrix_lib > 0.33) >= 3
keep_baits <- rowSums(cpm_matrix_baits > 0.5) >= 3

counts_filtered_lib <- dge_lib$counts[keep_lib, ]
counts_filtered_baits <- dge_baits$counts[keep_baits, ]

colSums(counts_filtered_lib)
colSums(counts_filtered_baits)

### Round counts for rarefaction curve
counts_filtered_rounded_baits <- round(counts_filtered_baits)
head(counts_filtered_rounded_baits)

counts_filtered_rounded_lib <- round(counts_filtered_lib)
head(counts_filtered_rounded_lib)

####################################################################### TESTING
# Want to see if CYP counts match with expected rank-order based on Kubota et al paper

# --- 1. Load gene ID list from file ---
cyp_gene_ids <- read.table("~/Documents/cyp_gene_ids_only.tsv", header = FALSE, stringsAsFactors = FALSE)$V1

# --- 2. Define male/female sample columns ---
baits_male_samples   <- c("E1_S1_L001_quant", "E2_S2_L001_quant", "E3_S3_L001_quant")
baits_female_samples <- c("E4_S4_L001_quant", "E5_S5_L001_quant", "E6_S6_L001_quant")

# --- 3. Keep only genes present in the matrix ---
present_cyp_genes <- intersect(cyp_gene_ids, rownames(counts_filtered_rounded_baits))

# --- 4. Subset count matrix ---
cyp_counts <- counts_filtered_rounded_baits[present_cyp_genes, , drop = FALSE]

# --- 5. Sum counts across male/female samples ---
cyp_counts_male   <- rowSums(cyp_counts[, baits_male_samples, drop = FALSE])
cyp_counts_female <- rowSums(cyp_counts[, baits_female_samples, drop = FALSE])

# --- 6. Convert to percent of total counts from cyp list ---
percent_male   <- 100 * cyp_counts_male / sum(cyp_counts_male)
percent_female <- 100 * cyp_counts_female / sum(cyp_counts_female)

# --- 7. Assemble and sort ---
ranked_cyp <- data.frame(
  GeneID = names(percent_male),
  MalePercent = percent_male,
  FemalePercent = percent_female
)

ranked_cyp_male   <- ranked_cyp[order(-ranked_cyp$MalePercent), ]
ranked_cyp_female <- ranked_cyp[order(-ranked_cyp$FemalePercent), ]

# --- 8. Optional: Show top 10 ---
cat("Top 15 CYP genes by % in MALES:\n")
print(head(ranked_cyp_male, 15), row.names = FALSE)

cat("\nTop 15 CYP genes by % in FEMALES:\n")
print(head(ranked_cyp_female, 15), row.names = FALSE)





library(ggplot2)
library(dplyr)

# ---- Function to prepare pie chart data ----
prepare_pie_data <- function(percent_vector, threshold = 4) {
  df <- data.frame(
    GeneID = names(percent_vector),
    Percent = as.numeric(percent_vector)
  )
  
  df <- df %>%
    mutate(GeneID = ifelse(Percent < threshold, "Others", GeneID)) %>%
    group_by(GeneID) %>%
    summarise(Percent = sum(Percent)) %>%
    arrange(desc(Percent))
  
  df$Label <- paste0(df$GeneID, " (", round(df$Percent, 1), "%)")
  return(df)
}

# ---- Prepare data for each sex ----
male_pie_data   <- prepare_pie_data(percent_male)
female_pie_data <- prepare_pie_data(percent_female)

# ---- Function to plot pie chart ----
plot_pie_chart <- function(pie_data, title) {
  ggplot(pie_data, aes(x = "", y = Percent, fill = Label)) +
    geom_col(width = 1, color = "white") +
    coord_polar("y") +
    theme_void() +
    labs(title = title, fill = "Gene") +
    theme(legend.position = "right")
}

# ---- Plotting both pie charts ----
plot_pie_chart(male_pie_data, "Male CYP Gene Expression")
plot_pie_chart(female_pie_data, "Female CYP Gene Expression")



###TESTING

cat("Male 'Others':\n")
print(male_result$others)

cat("\nFemale 'Others':\n")
print(female_result$Others)





# Transpose: rarecurve expects samples as rows
abundance_matrix_baits <- t(counts_filtered_rounded_baits)
abundance_matrix_lib <- t(counts_filtered_rounded_lib)

# Transpose: rarecurve expects samples as rows- can directly transpose from filtered_rounded_counts instead
#abundance_matrix_t <- t(abundance_matrix)

# Plot the rarefaction curves
rarecurve(abundance_matrix_baits,
          step = 100,
          sample = min(rowSums(abundance_matrix_baits)),  # use smallest library size
          col = rainbow(nrow(abundance_matrix_baits)),    # unique color per sample
          label = TRUE,
          xlab = "Sequencing Depth (Subsampled Reads)",
          ylab = "Number of Transcripts Detected",
          main = "Rarefaction Curves per Sample")

rareslope(abundance_matrix_baits,
          sample = 1000000)  # use to find slope at specific depth

rareslope(abundance_matrix_lib,
          sample = 4000000)  # use to find slope at specific depth
####THE ABOVE RARECURVE WORKS!!- approx 4M reads is where 96.9%-100% is captured for all samples (enriched)

rarecurve(abundance_matrix_lib,
          step = 1000,
          sample = min(rowSums(abundance_matrix_lib)),  # use smallest library size
          col = rainbow(nrow(abundance_matrix_lib)),    # unique color per sample
          label = TRUE,
          xlab = "Sequencing Depth (Subsampled Reads)",
          ylab = "Number of Genes Detected",
          main = "Rarefaction Curves per Sample")

# Cutoff samples at 4M depth; need to use the counts_filtered_x as input instead of rounded counts
# We need to use rounded counts here since rrarefy requires rounded (integers) 
library(vegan)

# Step 1: Set the target depth
target_depth <- 4e6

# Step 2: Calculate total reads per sample
sample_totals_lib <- rowSums(abundance_matrix_lib)
sample_totals_baits <- rowSums(abundance_matrix_baits)

# Step 3: Filter out samples that don't meet the depth requirement
abundance_matrix_rarefy_lib <- abundance_matrix_lib[sample_totals_lib >= target_depth, ]
abundance_matrix_rarefy_baits <- abundance_matrix_baits[sample_totals_baits >= target_depth, ]

# Step 4: Rarefy to exactly 4 million reads
rarefied_counts_lib <- rrarefy(abundance_matrix_rarefy_lib, sample = target_depth) 
rarefied_counts_baits <- rrarefy(abundance_matrix_rarefy_baits, sample = target_depth) 

# Step 5: Check output
head(rarefied_counts_lib)
head(rarefied_counts_baits)
rowSums(rarefied_counts_lib)  # Should all equal 4e6
rowSums(rarefied_counts_baits)

rarecurve(abundance_matrix_rarefy_lib,
                  step = 1000,  # Smoother curves
                  sample = target_depth,
                  col = rainbow(nrow(abundance_matrix_rarefy_lib)),
                  label = TRUE,
                  xlab = "Sequencing Depth (Subsampled Reads)",
                  ylab = "Number of Genes Detected",
                  main = "Rarefaction Curves per Sample (to 4M)")

rarecurve(abundance_matrix_rarefy_baits,
          step = 1000,  # Smoother curves
          sample = target_depth,
          col = rainbow(nrow(abundance_matrix_rarefy_baits)),
          label = TRUE,
          xlab = "Sequencing Depth (Subsampled Reads)",
          ylab = "Number of Genes Detected",
          main = "Rarefaction Curves per Sample (to 4M)")


# Now within these genes, I want to filter only for chemical defensome genes and see how many on target reads we get; need to filter low counts again

library(edgeR)

# Convert to DGEList for CPM calculation
dge_lib2 <- DGEList(counts = t(rarefied_counts_lib))      # transpose to samples in columns
dge_baits2 <- DGEList(counts = t(rarefied_counts_baits))  # same here

# Calculate CPM
cpm_lib2 <- cpm(dge_lib2)
cpm_baits2 <- cpm(dge_baits2)

# Apply CPM filter: keep genes with CPM > 0.5 in at least 3 samples
keep_lib2 <- rowSums(cpm_lib2 > 0.5) >= 3
keep_baits2 <- rowSums(cpm_baits2 > 0.5) >= 3

# Filter rarefied count matrices (transpose back so samples are rows)
filtered_rarefied_counts_lib <- t(dge_lib2$counts[keep_lib2, ])
filtered_rarefied_counts_baits <- t(dge_baits2$counts[keep_baits2, ])

# Read and clean your list of target gene IDs
defense.genes.list <- readLines("~/Documents/complete_ncbi_gene_ids(redo595).txt")
defense.genes.list<- sub("\\.\\d+$", "", defense.genes.list)  # strip version numbers
defense.genes.list <- trimws(defense.genes.list)              # remove whitespace

# Subset columns in rarefied_counts
defensome_rarefied_counts_lib <- filtered_rarefied_counts_lib[, colnames(filtered_rarefied_counts_lib) %in% defense.genes.list]
defensome_rarefied_counts_baits <- filtered_rarefied_counts_baits[, colnames(filtered_rarefied_counts_baits) %in% defense.genes.list]

# Drop genes with all-zero counts (optional)- There were none so skipped this
# defensome_rarefied_counts_lib <- defensome_rarefied_counts_lib[, colSums(defensome_rarefied_counts_lib) > 0]
# defensome_rarefied_counts_baits <- defensome_rarefied_counts_baits[, colSums(defensome_rarefied_counts_baits) > 0]

# Check result
cat("Matched genes:", ncol(defensome_rarefied_counts_lib), "of", length(defense.genes.list), "\n")
head(defensome_rarefied_counts_lib)
cat("Matched genes:", ncol(defensome_rarefied_counts_baits), "of", length(defense.genes.list), "\n")
head(defensome_rarefied_counts_baits)

# Check if any genes are missing 
missing_genes_lib <- setdiff(defense.genes.list, colnames(defensome_rarefied_counts_lib))
cat("Unmatched gene IDs:", length(missing_genes_lib), "\n")
cat(missing_genes_lib)

missing_genes_baits <- setdiff(defense.genes.list, colnames(defensome_rarefied_counts_baits))
cat("Unmatched gene IDs:", length(missing_genes_baits), "\n")
cat(missing_genes_baits)

# Are any genes missing between the baited and unbaited
# Genes missing from both lib and baits
missing_in_both <- intersect(missing_genes_lib, missing_genes_baits)

# Genes only missing from lib
only_missing_in_lib <- setdiff(missing_genes_lib, missing_genes_baits)

# Genes only missing from baits
only_missing_in_baits <- setdiff(missing_genes_baits, missing_genes_lib)

# Genes present in both filtered lib and baits datasets
genes_in_both <- intersect(
  colnames(defensome_rarefied_counts_lib),
  colnames(defensome_rarefied_counts_baits)
)

# Report counts
cat("Genes found in BOTH lib and baits:", length(genes_in_both), "\n")
#print(genes_in_both)
cat("Missing in BOTH  :", length(missing_in_both), "\n")
cat("Only missing in LIB       :", length(only_missing_in_lib), "\n")
cat("Only missing in BAITS     :", length(only_missing_in_baits), "\n")

cat(only_missing_in_baits)

### Any genes in the baited at high levels that were not actually part of defensome target list?

#Get genes in baits not in defense list
extra_genes_baits <- setdiff(colnames(rarefied_counts_baits), defense.genes.list)

# Step 2: Subset those genes
extra_counts_baits <- rarefied_counts_baits[, extra_genes_baits]

# Step 3: Compute total or average counts per gene
total_counts <- colSums(extra_counts_baits)
# or: avg_counts <- colMeans(extra_counts_baits)

# Step 4: Sort to get top-expressed genes
top_extra_genes <- sort(total_counts, decreasing = TRUE)

# View top ones
head(top_extra_genes, 20)


# View a specific gene's counts; To put it into a data frame for plotting or export- 
df_lib_792610  <- data.frame(
  Sample = rownames(rarefied_counts_lib),
  Counts = rarefied_counts_lib[, "792610"]
)
head(df_lib_792610)

df_baits_792610 <- data.frame(
  Sample = rownames(rarefied_counts_baits),
  Counts = rarefied_counts_baits[, "792610"]
)
head(df_baits_792610)

# Number of reads that are attributed to the defensome gene list at 4M reads in either set
# Step 1: Get gene IDs that are present in both rarefied counts and your defense gene list
present_genes_lib   <- intersect(defense.genes.list, colnames(rarefied_counts_lib))
present_genes_baits <- intersect(defense.genes.list, colnames(rarefied_counts_baits))

# Step 2: Subset the rarefied matrices to defense genes only
defense_counts_lib   <- rarefied_counts_lib[, present_genes_lib, drop = FALSE]
defense_counts_baits <- rarefied_counts_baits[, present_genes_baits, drop = FALSE]

# Step 3: Sum counts per sample
defense_reads_per_sample_lib   <- rowSums(defense_counts_lib)
defense_reads_per_sample_baits <- rowSums(defense_counts_baits)

# Step 4: Output
cat("Defense gene reads per sample (RNA-seq):\n")
print(defense_reads_per_sample_lib)

cat("\nDefense gene reads per sample (Target enrichment):\n")
print(defense_reads_per_sample_baits)


################################################################################
# To determine if baits and library counts are linear for quatification purposes

# Define male and female sample names
lib_male_samples   <- c("L1_S1_L001_quant", "L2_S2_L001_quant", "L3_S3_L001_quant")
lib_female_samples <- c("L4_S4_L001_quant", "L5_S5_L001_quant", "L6_S6_L001_quant")

baits_male_samples   <- c("E1_S1_L001_quant", "E2_S2_L001_quant", "E3_S3_L001_quant")
baits_female_samples <- c("E4_S4_L001_quant", "E5_S5_L001_quant", "E6_S6_L001_quant")

# Get only the defense genes that are actually present in the count matrices
present_genes_lib <- intersect(defense.genes.list, rownames(counts_filtered_rounded_lib))
present_genes_baits <- intersect(defense.genes.list, rownames(counts_filtered_rounded_baits))

# Subset counts using only those genes
lib_defense_counts <- counts_filtered_rounded_lib[present_genes_lib, , drop = FALSE]
baits_defense_counts <- counts_filtered_rounded_baits[present_genes_baits, , drop = FALSE]

# Sum counts across male and female samples
lib_male   <- rowSums(lib_defense_counts[, lib_male_samples, drop = FALSE])
lib_female <- rowSums(lib_defense_counts[, lib_female_samples, drop = FALSE])

baits_male   <- rowSums(baits_defense_counts[, baits_male_samples, drop = FALSE])
baits_female <- rowSums(baits_defense_counts[, baits_female_samples, drop = FALSE])

# Create a data frame for plotting — only include genes found in both lib and baits
common_genes_male <- intersect(names(lib_male), names(baits_male))
common_genes_female <- intersect(names(lib_female), names(baits_female))

# Prepare data for male
df_male <- data.frame(
  Gene = common_genes_male,
  Lib_Male = lib_male[common_genes_male],
  Baits_Male = baits_male[common_genes_male]
)

# Prepare data for female
df_female <- data.frame(
  Gene = common_genes_female,
  Lib_Female = lib_female[common_genes_female],
  Baits_Female = baits_female[common_genes_female]
)

# Plotting
par(mfrow = c(1, 2))  # Two plots side by side

# Male plot
plot(log10(df_male$Lib_Male + 1), log10(df_male$Baits_Male + 1),
     main = "Male Samples", xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)",
     pch = 16, col = "blue")
abline(0, 1, col = "gray", lty = 2)

# Female plot
plot(log10(df_female$Lib_Female + 1), log10(df_female$Baits_Female + 1),
     main = "Female Samples", xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)",
     pch = 16, col = "red")
abline(0, 1, col = "gray", lty = 2)

par(mfrow = c(1,1))  # Reset layout


# To compute R2
# Log-transform counts (+1 to avoid log(0))
log_lib_male   <- log10(df_male$Lib_Male + 1)
log_baits_male <- log10(df_male$Baits_Male + 1)

log_lib_female   <- log10(df_female$Lib_Female + 1)
log_baits_female <- log10(df_female$Baits_Female + 1)

# To identify outliers, define a threshold (adjust as needed)
threshold <- 1

# Outliers = where |log(baits) - log(lib)| > threshold
male_outliers_idx <- abs(log_baits_male - log_lib_male) > threshold
female_outliers_idx <- abs(log_baits_female - log_lib_female) > threshold

par(mfrow = c(1, 2))

# Male
plot(log_lib_male, log_baits_male, pch = 16, col = ifelse(male_outliers_idx, "black", "blue"),
     main = "Male Samples", xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)")
abline(a = 0, b = 1, lty = 2, col = "grey")

# Female
plot(log_lib_female, log_baits_female, pch = 16, col = ifelse(female_outliers_idx, "black", "red"),
     main = "Female Samples", xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)")
abline(a = 0, b = 1, lty = 2, col = "grey")



# Fit linear models
fit_male   <- lm(log_baits_male ~ log_lib_male)
fit_female <- lm(log_baits_female ~ log_lib_female)

# Extract R²
r2_male   <- summary(fit_male)$r.squared
r2_female <- summary(fit_female)$r.squared

# Print results
cat("R² for Male samples:", round(r2_male, 4), "\n")
cat("R² for Female samples:", round(r2_female, 4), "\n")

# Replot without outliers
# Remove outliers
log_lib_male_no_outliers   <- log_lib_male[!male_outliers_idx]
log_baits_male_no_outliers <- log_baits_male[!male_outliers_idx]

log_lib_female_no_outliers   <- log_lib_female[!female_outliers_idx]
log_baits_female_no_outliers <- log_baits_female[!female_outliers_idx]

# Re-plot without outliers
par(mfrow = c(1, 2))  # Side-by-side plots

# Male
plot(log_lib_male_no_outliers, log_baits_male_no_outliers,
     col = "blue", pch = 16, main = "Male Samples (No Outliers)",
     xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)")
abline(0, 1, lty = 2, col = "gray")

# Female
plot(log_lib_female_no_outliers, log_baits_female_no_outliers,
     col = "red", pch = 16, main = "Female Samples (No Outliers)",
     xlab = "Library (log10 counts)", ylab = "Baits (log10 counts)")
abline(0, 1, lty = 2, col = "gray")

# Recalculate R² values
r2_male <- summary(lm(log_baits_male_no_outliers ~ log_lib_male_no_outliers))$r.squared
r2_female <- summary(lm(log_baits_female_no_outliers ~ log_lib_female_no_outliers))$r.squared

cat("R² (Male, no outliers):", round(r2_male, 3), "\n")
cat("R² (Female, no outliers):", round(r2_female, 3), "\n")


# To identify which were the outlier genes that were removed
gene_names_male   <- names(log_lib_male)
gene_names_female <- names(log_lib_female)

# Use logical indices to extract outlier gene names
male_outlier_genes   <- gene_names_male[male_outliers_idx]
female_outlier_genes <- gene_names_female[female_outliers_idx]

# Output results
cat("Number of male outlier genes:", length(male_outlier_genes), "\n")
print(male_outlier_genes)

cat("\nNumber of female outlier genes:", length(female_outlier_genes), "\n")
print(female_outlier_genes)





########################################################################################
# Then filter for how many DEGs between MvF- DESeq2 on raw imported counts from tximport
install.packages("DESeq2")
library(DESeq2)

samples_baits <- c("E1_S1_L001_quant",
             "E2_S2_L001_quant",
             "E3_S3_L001_quant",
             "E4_S4_L001_quant",
             "E5_S5_L001_quant",
             "E6_S6_L001_quant")

samples_lib <- c("L1_S1_L001_quant",
             "L2_S2_L001_quant",
             "L3_S3_L001_quant",
             "L4_S4_L001_quant",
             "L5_S5_L001_quant",
             "L6_S6_L001_quant")

print(samples_lib)
print(samples_baits)

sex_design_lib <- data.frame(
  sample= samples_lib,
  background=c("male","male","male","female","female","female")
)
print(sex_design_lib)

sex_design_baits <- data.frame(
  sample= samples_baits,
  background=c("male","male","male","female","female","female")
)
print(sex_design_baits)

storage.mode(counts_filtered_rounded_lib) <- "integer"
storage.mode(counts_filtered_rounded_baits) <- "integer"

dds_lib <- DESeqDataSetFromMatrix(
  countData = counts_filtered_rounded_lib,  
  colData = sex_design_lib,
  design = ~background
)


dds_baits <- DESeqDataSetFromMatrix(
  countData = counts_filtered_rounded_baits,  
  colData = sex_design_baits,
  design = ~background
)


#This makes PCA for conditions (male vs female)
lib_pca <- rlog(dds_lib, 
                blind = TRUE)

baits_pca <- rlog(dds_baits, 
                blind = TRUE)

dim(lib_pca)
dim(baits_pca)

plotPCA(lib_pca, 
        intgroup=c("background"),
        ntop = 10000)

plotPCA(baits_pca, 
        intgroup=c("background"),
        ntop = 1000)

#Let's test out DESeq2 on the entire imported dataset first
dds_lib2 <- DESeq(dds_lib)

dds_baits2 <- DESeq(dds_baits)

plotDispEsts(dds_lib2,
             legend = F)

plotDispEsts(dds_baits2,
             legend = T)

txi_res_lib <- results(dds_lib2)
txi_res_lib
resultsNames(dds_lib2)
summary(dds_lib2)

txi_res_baits <- results(dds_baits2)
txi_res_baits
resultsNames(dds_baits2)


# LFC shrinkage of estimates; using apeglm since this is newer and is "said" to have better application
BiocManager::install("apeglm")
txi_res_lib_shrunk <- lfcShrink(dds_lib2, coef="background_male_vs_female", type="apeglm")
txi_res_baits_shrunk <- lfcShrink(dds_baits2, coef="background_male_vs_female", type="apeglm")



