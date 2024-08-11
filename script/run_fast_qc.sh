#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_fastq1> <path_to_fastq2>"
    exit 1
fi

# Paths to FastQ files
fastq1=$1
fastq2=$2

# Output directory for unzipped files
output_dir="gunzip_data"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Unzip the FastQ files
echo "Unzipping $fastq1 to $output_dir/$(basename ${fastq1%.gz})..."
gunzip -c $fastq1 > $output_dir/$(basename ${fastq1%.gz})

echo "Unzipping $fastq2 to $output_dir/$(basename ${fastq2%.gz})..."
gunzip -c $fastq2 > $output_dir/$(basename ${fastq2%.gz})

# Run FastQC on the unzipped files
echo "Running FastQC on $output_dir/$(basename ${fastq1%.gz})..."
fastqc $output_dir/$(basename ${fastq1%.gz})

echo "Running FastQC on $output_dir/$(basename ${fastq2%.gz})..."
fastqc $output_dir/$(basename ${fastq2%.gz})

echo "FastQC analysis complete."
