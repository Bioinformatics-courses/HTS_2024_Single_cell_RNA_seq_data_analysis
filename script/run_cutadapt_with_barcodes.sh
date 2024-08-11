#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_fastq1> <path_to_fastq2> <path_to_fasta_barcode_file>"
    exit 1
fi

# Paths to FastQ files and FASTA barcode file
fastq1=$1
fastq2=$2
barcode_file=$3

# Output directory for trimmed FASTQs
output_dir="cutadapt_output"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Run cutadapt to trim barcodes using the FASTA barcode file
cutadapt -g file:$barcode_file \
         -o $output_dir/trimmed_1.fastq \
         -p $output_dir/trimmed_2.fastq \
         $fastq1 $fastq2

# Verify the output
echo "Verifying barcode trimming..."
head -n 4 $output_dir/trimmed_1.fastq
head -n 4 $output_dir/trimmed_2.fastq

echo "Barcode trimming complete."

