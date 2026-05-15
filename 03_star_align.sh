#!/bin/bash
# 03_star_align.sh
# Genome indexing and alignment using STAR (v2.7.11b)
# Reference: Danio rerio GRCz12tu (GCF_049306965.1)
# --outFilterMultimapNmax 1: reads mapping to multiple loci are excluded
# --sjdbOverhang 50: set to read_length - 1 for 2x50 bp configuration
# Run on both liver (run 1) and tissue (run 2) libraries

set -euo pipefail

star="/usr/lib/rna-star/bin/STAR-avx"
threads=16

genome_fa="/home/mellissa/baits_liver/ncbi_dataset/data/GCF_049306965.1/GCF_049306965.1_GRCz12tu_genomic.fna"
gtf="/home/mellissa/baits_liver/ncbi_dataset/data/GCF_049306965.1/genomic.gtf"
index_dir="/home/mellissa/baits_liver/STAR_GRCz12tu_index"
echo "Building STAR index"
STAR --runMode genomeGenerate \
    --genomeDir "$index_dir" \
    --genomeFastaFiles "$genome_fa" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang 50

run_star() {
    local trimmed_dir=$1
    local out_dir=$2

    mkdir -p "$out_dir"

    for infile in "${trimmed_dir}"/*_R1_PE.fastq; do
        base=$(basename "${infile}" _R1_PE.fastq)
        echo "Aligning: $base"

        "$star" \
            --runThreadN "$threads" \
            --genomeDir "$index_dir" \
            --readFilesIn \
                "${trimmed_dir}/${base}_R1_PE.fastq" \
                "${trimmed_dir}/${base}_R2_PE.fastq" \
            --outFileNamePrefix "${out_dir}/${base}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterMatchNmin 40 \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterMultimapNmax 1 \
            --alignEndsType EndToEnd \
            --twopassMode Basic \
            --outSAMmapqUnique 60 \
            --limitBAMsortRAM 16000000000
    done

    # Rename BAM files to remove the verbose STAR suffix
    echo "Renaming BAM files"
    for f in "${out_dir}"/*_Aligned.sortedByCoord.out.bam; do
        base=$(basename "$f" _Aligned.sortedByCoord.out.bam)
        mv "$f" "${out_dir}/${base}.bam"
    done
}
echo "Aligning liver libraries"
run_star \
    "/home/mellissa/baits_liver/trimmed_redo_reads" \
    "/home/mellissa/baits_liver/STAR_alignment"
echo "Aligning tissue libraries"
run_star \
    "/home/mellissa/baits_tissues/trimmed_reads" \
    "/home/mellissa/baits_tissues/STAR_alignment"

echo "STAR alignment complete"
# ~98-99% of enriched reads and ~96-97% of unenriched reads mapped to GRCz12tu
