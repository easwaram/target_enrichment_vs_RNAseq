# Load STAR featurecounts results
# Read STAR featureCounts output (skip first row of comments)
counts_lib_STAR <- read.delim("~/Documents/STAR_featureCounts/STAR_library_liver_gene_counts.txt", comment.char = "#", check.names = FALSE)
counts_baits_STAR <- read.delim("~/Documents/STAR_featureCounts/STAR_baits_liver_gene_counts.txt", comment.char = "#", check.names = FALSE)

# Drop annotation columns (first 6 columns) and keep count matrix only
rownames(counts_lib_STAR) <- counts_lib_STAR$Geneid
lib_matrix_STAR <- as.matrix(counts_lib_STAR[, 7:ncol(counts_lib_STAR)])

rownames(counts_baits_STAR) <- counts_baits_STAR$Geneid
baits_matrix_STAR <- as.matrix(counts_baits_STAR[, 7:ncol(counts_baits_STAR)])

# Clean column names in STAR CPM
colnames(lib_matrix_STAR) <- gsub("\\.bam$", "", colnames(lib_matrix_STAR))
colnames(baits_matrix_STAR) <- gsub("\\.bam$", "", colnames(baits_matrix_STAR))
                                                                      

# Convert STAR counts to CPM
library(edgeR)

# Create a DGEList object (for easy normalization)
dge_lib_STAR   <- DGEList(counts = lib_matrix_STAR)
dge_baits_STAR <- DGEList(counts = baits_matrix_STAR)

# Calculate CPM
cpm_lib_STAR <- cpm(dge_lib_STAR)
cpm_baits_STAR <- cpm(dge_baits_STAR)

head (cpm_lib_STAR)
head (cpm_baits_STAR)

# Below filters out those with a CPM of greater than 0.33 in atleast 3 samples (so about 10 reads according to library size of 30M reads)
# Changed the CPM cut-off to 0.5 for baits since much lower depth
# This then filters those out and keeps the counts_filtered as raw counts to use for DESeq or rarefaction (both needs raw counts; rarefaction needs rounded)

# Filter for libraries: keep genes with CPM > 0.33 in at least 3 samples
keep_lib_STAR <- rowSums(cpm_lib_STAR > 0.33) >= 3
counts_filtered_lib_STAR   <- dge_lib_STAR$counts[keep_lib_STAR, ]

# Filter for baits: keep genes with CPM > 0.5 in at least 3 samples
keep_baits_STAR <- rowSums(cpm_baits_STAR > 0.5) >= 3
counts_filtered_baits_STAR <- dge_baits_STAR$counts[keep_baits_STAR, ]

# Create new filtered DGEList objects
dge_lib_filtered_STAR   <- DGEList(counts = counts_filtered_lib_STAR)
dge_baits_filtered_STAR <- DGEList(counts = counts_filtered_baits_STAR) 

# Calculate CPMs on filtered data
cpm_filtered_lib_STAR   <- cpm(dge_lib_filtered_STAR)
cpm_filtered_baits_STAR <- cpm(dge_baits_filtered_STAR) 

colSums(counts_filtered_lib_STAR)
colSums(counts_filtered_baits_STAR)

colSums(baits_matrix_STAR)
colSums(lib_matrix_STAR)

# Subset raw filtered counts with those in common with defensome gene list (for downstream purposes)

# Load defensome gene list (assumes one gene symbol per line, no header)
defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv", 
                              header = FALSE, stringsAsFactors = FALSE)[,1]

# Subset defensome genes that are present in the filtered baits and lib counts
common_defensome_baits <- intersect(rownames(counts_filtered_baits_STAR), defensome_genes)
common_defensome_lib   <- intersect(rownames(counts_filtered_lib_STAR), defensome_genes)

# Subset the counts to only defensome genes
baits_defensome_filtered_counts <- counts_filtered_baits_STAR[common_defensome_baits, , drop = FALSE]
lib_defensome_filtered_counts   <- counts_filtered_lib_STAR[common_defensome_lib, , drop = FALSE]

# Print how many defensome genes were found in each dataset
cat("Number of defensome genes found in baits dataset:", length(common_defensome_baits), "\n")
cat("Number of defensome genes found in lib dataset:", length(common_defensome_lib), "\n")



### Round counts for rarefaction curve; for STAR, all counts are already rounded in the raw count matrix
# Therefore, I skipped this step (refer to old_salmon_baits_analysis for old code)

# Transpose: rarecurve expects samples as rows
abundance_matrix_baits <- t(counts_filtered_rounded_baits)
abundance_matrix_lib <- t(counts_filtered_rounded_lib)


###############################################################################################


# Then filter for how many DEGs between MvF- DESeq2 on raw imported counts from tximport
install.packages("DESeq2")
library(DESeq2)

                                                                    

samples_baits <- c("E1_S1_L001",
                   "E2_S2_L001",
                   "E3_S3_L001",
                   "E4_S4_L001",
                   "E5_S5_L001",
                   "E6_S6_L001")

samples_lib <- c("L1_S1_L001",
                 "L2_S2_L001",
                 "L3_S3_L001",
                 "L4_S4_L001",
                 "L5_S5_L001",
                 "L6_S6_L001")

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

# storage.mode(counts_filtered_rounded_lib) <- "integer"
# storage.mode(counts_filtered_rounded_baits) <- "integer"

dds_lib_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered_lib_STAR,  
  colData = sex_design_lib,
  design = ~background
)


dds_baits_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered_baits_STAR,  
  colData = sex_design_baits,
  design = ~background
)

dds_lib_defensome <- DESeqDataSetFromMatrix( # made this just to run pca on subset
  countData = lib_defensome_filtered_counts,  # USE ONLY dds_all for DESeq2
  colData = sex_design_lib,
  design = ~background
)


dds_baits_defensome <- DESeqDataSetFromMatrix(
  countData = baits_defensome_filtered_counts,  
  colData = sex_design_baits,
  design = ~background
)


#This makes PCA for conditions (male vs female)
lib_pca_all <- rlog(dds_lib, 
                blind = TRUE)

baits_pca_all <- rlog(dds_baits, 
                  blind = TRUE)

dim(lib_pca_all)
dim(baits_pca_all)

plotPCA(lib_pca_all, 
        intgroup=c("background"),
        ntop = 10000)

plotPCA(baits_pca_all, 
        intgroup=c("background"),
        ntop = 1000)

lib_pca_defensome <- rlog(dds_lib_defensome, 
                    blind = TRUE)

baits_pca_defensome <- rlog(dds_baits_defensome, 
                      blind = TRUE)

dim(lib_pca_defensome)
dim(baits_pca_defensome)

plotPCA(lib_pca_defensome, 
        intgroup=c("background"),
        ntop = 400)

plotPCA(baits_pca_defensome, 
        intgroup=c("background"),
        ntop = 400)


#First run DESeq2 on the entire dataset - we filter for only defensome genes after
# shrinking estimates.
dds_lib2 <- DESeq(dds_lib_all)

dds_baits2 <- DESeq(dds_baits_all)

plotDispEsts(dds_lib2,
             legend = F)

plotDispEsts(dds_baits2,
             legend = T)

res_lib <- results(dds_lib2)
res_lib

res_baits <- results(dds_baits2)
res_baits

resultsNames(dds_lib2)
summary(dds_lib2)


###################################################################

# Now can compre LFC between lib and baits after DESeq2- will compare before and after 
# shrinking estimates

defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv", 
                              header = FALSE, stringsAsFactors = FALSE)[,1]

# Convert to data.frame and add gene names
res_lib_df   <- as.data.frame(res_lib)
res_baits_df <- as.data.frame(res_baits)

res_lib_df   <- as.data.frame(res_lib_shrunk)
res_baits_df <- as.data.frame(res_baits_shrunk)

res_lib_df$gene   <- rownames(res_lib_df)
res_baits_df$gene <- rownames(res_baits_df)

# Filter to keep only defensome genes
res_lib_df   <- res_lib_df[res_lib_df$gene %in% defensome_genes, ]
res_baits_df <- res_baits_df[res_baits_df$gene %in% defensome_genes, ]

# Merge by gene name (common defensome genes)
merged_res <- merge(res_lib_df[, c("gene", "log2FoldChange")],
                    res_baits_df[, c("gene", "log2FoldChange")],
                    by = "gene",
                    suffixes = c("_lib", "_baits"))

# Plot defensome gene log2FCs: lib vs baits
library(ggplot2)

ggplot(merged_res, aes(x = log2FoldChange_lib, 
                       y = log2FoldChange_baits, 
                       label = gene)) +
  geom_point(color = "darkgreen", size = 3) +
  #geom_text(vjust = -1.1, size = 3.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +  # 1:1 line in black
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray60", size = 0.8) +  # regression line lighter
  xlab("Library log2 Fold Change shrunk (male vs female)") +
  ylab("Baits log2 Fold Change shrunk (male vs female)") +
  ggtitle("Defensome Genes: log2FC shrunk Correlation Between Library and Baits") +
  theme_minimal()

# Spearman correlation
cor_test <- cor.test(merged_res$log2FoldChange_lib, merged_res$log2FoldChange_baits, method = "spearman")
cat("Spearman correlation:", round(cor_test$estimate, 3), 
    "p-value:", signif(cor_test$p.value, 3), "\n")


### LFC shrinkage of estimates; using apeglm since this is newer and is "said" to have better application
BiocManager::install("apeglm")
res_lib_shrunk <- lfcShrink(dds_lib2, coef="background_male_vs_female", type="apeglm")
res_baits_shrunk <- lfcShrink(dds_baits2, coef="background_male_vs_female", type="apeglm")


# Load defensome gene list
defensome_genes <- read.table("~/Documents/defensome_genes_symbol_only.tsv",
                              header = FALSE, stringsAsFactors = FALSE)[, 1]

# Subset results to defensome genes
res_lib_defensome <- res_lib_shrunk[rownames(res_lib_shrunk) %in% defensome_genes, ]
res_baits_defensome <- res_baits_shrunk[rownames(res_baits_shrunk) %in% defensome_genes, ]

# Cleanup subset of defensome results to make sure NA p-values are removed
res_lib_defensome <- na.omit(res_lib_defensome)
res_baits_defensome <- na.omit(res_baits_defensome)

summary(res_lib_defensome)
summary(res_baits_defensome)
