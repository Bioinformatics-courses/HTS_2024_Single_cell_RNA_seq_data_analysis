#!/bin/bash

# Check if necessary arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_fastq_R1> <input_fastq_R2> <adapters_fasta>"
    exit 1
fi

# Input arguments
INPUT_FASTQ_R1=$1
INPUT_FASTQ_R2=$2
ADAPTERS=$3

# Output directories
OUTPUT_DIR="output_step1"
QC_OUTPUT_DIR="${OUTPUT_DIR}/qc_results"
TRIMMED_OUTPUT_DIR="${OUTPUT_DIR}/trimmed_reads"

# Create output directories
mkdir -p $QC_OUTPUT_DIR
mkdir -p $TRIMMED_OUTPUT_DIR

# Full paths to tools (modify these if needed)
FASTQC_PATH="fastqc"
TRIMMOMATIC_PATH="trimmomatic"

# Quality Control using FastQC
echo "Running FastQC on forward reads..."
$FASTQC_PATH $INPUT_FASTQ_R1 -o $QC_OUTPUT_DIR

# Determine if paired-end or single-end
if [ -z "$INPUT_FASTQ_R2" ]; then
    # Single-end reads
    echo "Detected single-end reads."

    # Read Trimming using Trimmomatic
    echo "Running Trimmomatic for single-end reads..."
    TRIMMED_FASTQ_R1="${TRIMMED_OUTPUT_DIR}/$(basename ${INPUT_FASTQ_R1%.fastq}_trimmed.fastq)"
    $TRIMMOMATIC_PATH SE -phred33 \
        $INPUT_FASTQ_R1 $TRIMMED_FASTQ_R1 \
        ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
else
    # Paired-end reads
    echo "Detected paired-end reads."
    
    echo "Running FastQC on reverse reads..."
    $FASTQC_PATH $INPUT_FASTQ_R2 -o $QC_OUTPUT_DIR

    # Read Trimming using Trimmomatic
    echo "Running Trimmomatic for paired-end reads..."
    TRIMMED_FASTQ_R1="${TRIMMED_OUTPUT_DIR}/$(basename ${INPUT_FASTQ_R1%.fastq}_trimmed.fastq)"
    TRIMMED_FASTQ_R2="${TRIMMED_OUTPUT_DIR}/$(basename ${INPUT_FASTQ_R2%.fastq}_trimmed.fastq)"
    TRIMMED_FASTQ_R1_UNPAIRED="${TRIMMED_OUTPUT_DIR}/$(basename ${INPUT_FASTQ_R1%.fastq}_unpaired.fastq)"
    TRIMMED_FASTQ_R2_UNPAIRED="${TRIMMED_OUTPUT_DIR}/$(basename ${INPUT_FASTQ_R2%.fastq}_unpaired.fastq)"

    $TRIMMOMATIC_PATH PE -phred33 \
        $INPUT_FASTQ_R1 $INPUT_FASTQ_R2 \
        $TRIMMED_FASTQ_R1 $TRIMMED_FASTQ_R1_UNPAIRED \
        $TRIMMED_FASTQ_R2 $TRIMMED_FASTQ_R2_UNPAIRED \
        ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
fi

echo "Quality control and trimming completed. Results are in the $OUTPUT_DIR directory."
