#!/bin/bash
# 01_fastqc.sh
# FastQC on raw and trimmed reads for both sequencing runs
# Run 1: enriched and unenriched liver libraries
# Run 2: enriched gill, intestine, and kidney libraries

set -euo pipefail

threads=18
fastqc="/opt/local/fastqc/fastqc"

raw_liver="baits_liver/all_raw_fastq"
trimmed_liver="baits_liver/trimmed_redo_reads"
raw_qc_liver="baits_liver/fastqc_all_raw_reads"
trimmed_qc_liver="baits_liver/fastqc_trimmed_redo_paired_reads"

mkdir -p "$raw_qc_liver" "$trimmed_qc_liver"

echo "FastQC: raw liver reads"
"$fastqc" "${raw_liver}"/*.fastq.gz -o "$raw_qc_liver" -t "$threads"

echo "FastQC: trimmed liver reads"
"$fastqc" "${trimmed_liver}"/*_PE.fastq -o "$trimmed_qc_liver" -t "$threads"

raw_tissues="baits_tissues/all_raw_fastq"
trimmed_tissues="baits_tissues/trimmed_reads"
raw_qc_tissues="baits_tissues/fastqc_all_raw_reads"
trimmed_qc_tissues="baits_tissues/fastqc_trimmed_paired_reads"

mkdir -p "$raw_qc_tissues" "$trimmed_qc_tissues"

echo "FastQC: raw tissue reads"
"$fastqc" "${raw_tissues}"/*.fastq.gz -o "$raw_qc_tissues" -t "$threads"

echo "FastQC: trimmed tissue reads"
"$fastqc" "${trimmed_tissues}"/*_PE.fastq -o "$trimmed_qc_tissues" -t "$threads"

echo "FastQC complete"
