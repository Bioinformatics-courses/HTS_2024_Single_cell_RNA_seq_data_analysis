#!/bin/bash

# Check if necessary arguments are provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <reference_genome_index_base> <trimmed_fastq_R1> <trimmed_fastq_R2> <num_threads>"
    exit 1
fi

# Input arguments
REFERENCE_GENOME_INDEX_BASE=$1
TRIMMED_FASTQ_R1=$2
TRIMMED_FASTQ_R2=$3
NUM_THREADS=$4

# Output directory
OUTPUT_DIR="output_alignment"
ALIGNMENT_OUTPUT_DIR="${OUTPUT_DIR}/alignment_results"

# Create output directory
mkdir -p $ALIGNMENT_OUTPUT_DIR

# Full path to BWA (modify this if needed)
BWA_PATH="bwa"

# Align reads using BWA
echo "Aligning reads..."
if [ -z "$TRIMMED_FASTQ_R2" ]; then
    # Single-end alignment
    echo "Aligning single-end reads..."
    /usr/bin/time -v $BWA_PATH mem -t $NUM_THREADS $REFERENCE_GENOME_INDEX_BASE $TRIMMED_FASTQ_R1 > ${ALIGNMENT_OUTPUT_DIR}/aligned_reads.sam
else
    # Paired-end alignment
    echo "Aligning paired-end reads..."
    /usr/bin/time -v $BWA_PATH mem -t $NUM_THREADS $REFERENCE_GENOME_INDEX_BASE $TRIMMED_FASTQ_R1 $TRIMMED_FASTQ_R2 > ${ALIGNMENT_OUTPUT_DIR}/aligned_reads.sam
fi

echo "Alignment completed. Results are in the $ALIGNMENT_OUTPUT_DIR directory."
