#!/bin/bash
# 04_featurecounts.sh
# Read counting per gene using featureCounts (v2.0.6)
# -s 2: reversely stranded (NEBNext directional library prep)
# -p: paired-end mode; -B: both reads must map; -C: exclude chimeric pairs

set -euo pipefail

gtf="/home/mellissa/baits_liver/ncbi_dataset/data/GCF_049306965.1/genomic.gtf"
threads=8

echo "featureCounts: liver libraries"
featureCounts -T "$threads" \
    -p -B -C -s 2 \
    -a "$gtf" \
    -o "/home/mellissa/baits_liver/STAR_baits_liver_gene_counts.txt" \
    /home/mellissa/baits_liver/STAR_alignment/*.bam

echo "featureCounts: tissue libraries"
featureCounts -T "$threads" \
    -p -B -C -s 2 \
    -a "$gtf" \
    -o "/home/mellissa/baits_tissues/STAR_baits_tissues_gene_counts.txt" \
    /home/mellissa/baits_tissues/STAR_alignment/*.bam

echo "featureCounts complete"
# Enriched reads: ~98-99% successfully assigned
# Unenriched reads: ~95-97% successfully assigned
# Output: tab-separated count matrix with featureCounts annotation columns (1-6)
# Columns 7+ are per-sample counts; import into R with read.delim(comment.char="#")
