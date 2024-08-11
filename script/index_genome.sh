#!/bin/bash

# Check if necessary arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <reference_genome_fasta>"
    exit 1
fi

# Input argument
REFERENCE_GENOME=$1

# Output directory
OUTPUT_DIR="output_index"
INDEX_DIR="${OUTPUT_DIR}/genome_index"

# Create output directories
mkdir -p $INDEX_DIR

# Full path to BWA (modify this if needed)
BWA_PATH="bwa"

# Check if the reference genome is already indexed
if [ ! -f "${REFERENCE_GENOME}.bwt" ]; then
    echo "Indexing reference genome..."
    $BWA_PATH index -p $INDEX_DIR/genome_index $REFERENCE_GENOME
    echo "Reference genome indexing completed. Indexed files are in $INDEX_DIR."
else
    echo "Reference genome is already indexed."
fi
