#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <GTF_FILE> <BAM_FILES> <NUM_THREADS>"
    exit 1
fi

# Arguments
GTF_FILE=$1
BAM_FILES=$2
NUM_THREADS=$3

# Output directory
OUTPUT_DIR="output_counts"
mkdir -p $OUTPUT_DIR

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/count_matrix.txt"

# Run featureCounts
echo "Generating count matrix..."
featureCounts -a $GTF_FILE -o $OUTPUT_FILE -T $NUM_THREADS -p $BAM_FILES

echo "Count matrix generated. Results are in the $OUTPUT_FILE file."

