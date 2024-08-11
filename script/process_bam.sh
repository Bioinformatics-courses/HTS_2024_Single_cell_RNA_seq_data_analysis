#!/bin/bash

# Usage: ./process_bam.sh <input_sam> <output_prefix>

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_sam> <output_prefix>"
    exit 1
fi

# Assign arguments to variables
INPUT_SAM=$1
OUTPUT_PREFIX=$2
OUTPUT_DIR="output_bam"
SAMTOOLS_PATH="samtools"  # Modify this if samtools is not in your PATH

# Create output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

# Convert SAM to BAM
echo "Converting SAM to BAM..."
${SAMTOOLS_PATH} view -Sb ${INPUT_SAM} > ${OUTPUT_DIR}/${OUTPUT_PREFIX}.bam

# Sort BAM file
echo "Sorting BAM file..."
${SAMTOOLS_PATH} sort ${OUTPUT_DIR}/${OUTPUT_PREFIX}.bam -o ${OUTPUT_DIR}/${OUTPUT_PREFIX}_sorted.bam

# Index BAM file
echo "Indexing BAM file..."
${SAMTOOLS_PATH} index ${OUTPUT_DIR}/${OUTPUT_PREFIX}_sorted.bam

echo "Processing completed. Sorted BAM file and index are in the ${OUTPUT_DIR} directory."
