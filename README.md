# HTS_2024_Single_cell_RNA_seq_data_analysis
Quality_check_and_trimming.sh
# Overview
This script performs the initial steps of single-cell RNA sequencing data preprocessing, including quality control and read trimming. It uses FastQC for quality control and Trimmomatic for read trimming.
# Prerequisites
Ensure the following tools are installed and accessible in your PATH:
    FastQC
    Trimmomatic
Script Usage
## Usage
./scRNA_seq_step1.sh <input_fastq_R1> <input_fastq_R2> <adapters_fasta>
    <input_fastq_R1>: Path to the forward read FASTQ file.
    <input_fastq_R2>: Path to the reverse read FASTQ file (for paired-end reads) or an empty string (for single-end reads).
    <adapters_fasta>: Path to the FASTA file containing adapter sequences.
