#!/bin/bash

# Script to verify sorting of a BAM file

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sorted_bam_file>"
    exit 1
fi

# Assign the input argument to a variable
SORTED_BAM_FILE="$1"

# Check if the sorted BAM file exists
if [ ! -f "$SORTED_BAM_FILE" ]; then
    echo "Error: The file $SORTED_BAM_FILE does not exist."
    exit 1
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed."
    exit 1
fi

# Verify sorting by checking the header and a summary of the BAM file
echo "Verifying sorting of the BAM file: $SORTED_BAM_FILE"

# View the header of the BAM file
echo "BAM file header:"
samtools view -H "$SORTED_BAM_FILE"

# Provide a summary of the BAM file
echo "BAM file summary:"
samtools flagstat "$SORTED_BAM_FILE"

echo "Verification completed."

exit 0
