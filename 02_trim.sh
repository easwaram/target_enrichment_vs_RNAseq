#!/bin/bash
# 02_trim.sh
# Adapter and quality trimming with Trimmomatic (paired-end mode)
# Adapters: NEB kit (NEB_adaptors.fa)
# MINLEN set to 40 (not 50) since sequencing was 2x50 bp configuration

set -euo pipefail

trimmomatic="/opt/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters="/home/mellissa/NEB_adaptors.fa"
threads=16

echo "Trimming liver libraries"
for infile in /home/mellissa/baits_liver/all_raw_fastq/*_R1_001.fastq.gz
do
  base=$(basename ${infile} _R1_001.fastq.gz)
  java -jar "$trimmomatic" PE -threads "$threads" \
    -trimlog /home/mellissa/baits_liver/trimmed_redo_reads/${base}.log \
    /home/mellissa/baits_liver/all_raw_fastq/${base}_R1_001.fastq.gz \
    /home/mellissa/baits_liver/all_raw_fastq/${base}_R2_001.fastq.gz \
    /home/mellissa/baits_liver/trimmed_redo_reads/${base}_R1_PE.fastq \
    /home/mellissa/baits_liver/trimmed_redo_reads/${base}_R1_SE.fastq \
    /home/mellissa/baits_liver/trimmed_redo_reads/${base}_R2_PE.fastq \
    /home/mellissa/baits_liver/trimmed_redo_reads/${base}_R2_SE.fastq \
    ILLUMINACLIP:${adapters}:2:30:10 \
    LEADING:3 TRAILING:3 MAXINFO:40:0.4 MINLEN:40
done

echo "Trimming tissue libraries"
for infile in /home/mellissa/baits_tissues/all_raw_fastq/*_R1_001.fastq.gz
do
  base=$(basename ${infile} _R1_001.fastq.gz)
  java -jar "$trimmomatic" PE -threads "$threads" \
    -trimlog /home/mellissa/baits_tissues/trimmed_reads/${base}.log \
    /home/mellissa/baits_tissues/all_raw_fastq/${base}_R1_001.fastq.gz \
    /home/mellissa/baits_tissues/all_raw_fastq/${base}_R2_001.fastq.gz \
    /home/mellissa/baits_tissues/trimmed_reads/${base}_R1_PE.fastq \
    /home/mellissa/baits_tissues/trimmed_reads/${base}_R1_SE.fastq \
    /home/mellissa/baits_tissues/trimmed_reads/${base}_R2_PE.fastq \
    /home/mellissa/baits_tissues/trimmed_reads/${base}_R2_SE.fastq \
    ILLUMINACLIP:${adapters}:2:30:10 \
    LEADING:3 TRAILING:3 MAXINFO:40:0.4 MINLEN:40
done

echo "Trimming complete"
# Note: >99.6% of paired-end reads survived trimming across all samples
