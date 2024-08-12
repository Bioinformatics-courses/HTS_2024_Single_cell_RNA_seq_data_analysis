# Single-Cell RNA Sequencing Analysis of Plasmodium Species

#Introduction
This repository contains the code and scripts used to replicate and extend the single-cell RNA sequencing (scRNA-seq) analysis of Plasmodium species, specifically [Plasmodium falciparum](https://en.wikipedia.org/wiki/Plasmodium_falciparum) and [Plasmodium berghei](https://en.wikipedia.org/wiki/Plasmodium_berghei). These species are responsible for malaria, a devastating disease that affects millions worldwide.
The primary goal of this project is to analyze gene expression at the single-cell level in malaria parasites, providing insights into the biology of these pathogens and their interactions with the host. By employing scRNA-seq, we aim to identify species-specific gene expression patterns and potential targets for therapeutic intervention.

# Setting up Our Working Environment
NOTE!
This project was demonstrated on a Linux system.
Throughout this project, we'll be using several tools that will be installed via [conda](https://conda.io/projects/conda/en/latest/).

cd
git clone [https://github.com/Bioinformatics-courses/HTS_2024_Single_cell_RNA_seq_data_analysis]
mv HTS_2024_Single_cell_RNA_seq_data_analysis
cd mv HTS_2024_Single_cell_RNA_seq_data_analysis
mkdir -p data

CODE BREAKDOWN
cd - Moves to your home directory.
git clone - Clones the GitHub repository to your local machine.
cd HTS_2024_Single_cell_RNA_seq_data_analysis - Navigates into the cloned directory.
mkdir -p data - Creates a data directory inside the project.

Next we create our conda environment from the yml file provided in the env folder and activate the environment

conda activate bif
conda env export --name bif --file scRNA_seq.yml

# The Data
This project utilizes [single-cell RNA-sequencing (scRNA-seq) data](https://www.ebi.ac.uk/ena/browser/view/PRJEB19245) from the study titled "Single-cell transcriptomics of malaria parasites" (Project: PRJEB19245), provided by the Wellcome Sanger Institute. This data explores transcriptional variation within populations of Plasmodium parasites, specifically Plasmodium berghei and Plasmodium falciparum, revealing gene families that play a role in immune evasion during transmission to mosquitoes. The study represents one of the first instances of scRNA-seq being used to uncover transcriptional variation in unicellular eukaryotic organisms. For detailed information on data usage and publication policies, refer to the [Wellcome Trust Sanger Institute Data Sharing Policy](http://www.sanger.ac.uk/datasharing/).

## Quality check
The first script we're going to run is the Quality_check_and_trimming.sh. This script runs fastqc on the fastq files in data/raw and generates a fastqc report. It also run trimming, and is used only when the data requires trimming. This script performs the initial steps of single-cell RNA sequencing data preprocessing, including quality control and read trimming. It uses FastQC for quality control and Trimmomatic for read trimming.
# Prerequisites
Ensure the following tools are installed and accessible in your PATH:
    FastQC
    Trimmomatic
Script Usage
## Usage
./Quality_check_and_trimming.sh <input_fastq_R1> <input_fastq_R2> <adapters_fasta>
    <input_fastq_R1>: Path to the forward read FASTQ file.
    <input_fastq_R2>: Path to the reverse read FASTQ file (for paired-end reads) or an empty string (for single-end reads).
    <adapters_fasta>: Path to the FASTA file containing adapter sequences.
NOTE:
If the script does not require trimming use this ./run_fast_qc.sh <path_to_fastq1> <path_to_fastq2>

cd script
sudo chmod u+x Quality_check_and_trimming.sh
./Quality_check_and_trimming.sh

CODE BREAKDOWN
cd script - Move to the script directory.
sudo chmod u+x Quality_check_and_trimming.sh - Make the script executable.
./Quality_check_and_trimming.sh - Run the script.


<img width="501" alt="image" src="https://github.com/user-attachments/assets/9208619d-7391-4215-9348-d435fc82353f">

# Indexing the Genome
Next, Before aligning sequencing reads to the reference genome, it is essential to index the genome. Indexing creates a data structure that allows the aligner to efficiently map the reads to the genome. To perform genome indexing in this project, we use this bash script provided in the script directory:

./index_genome.sh <reference_genome_fasta>

The following files would be generated after indexing:
<img width="782" alt="image" src="https://github.com/user-attachments/assets/0517d0d2-d8eb-4d2a-b52d-5c3242addfcb">

# Barcode Trimming Script
This script is designed to trim barcodes from paired-end FASTQ files using a provided FASTA barcode file. It leverages the cutadapt tool to remove the barcode sequences, which is an essential preprocessing step in single-cell RNA sequencing analysis. To perform Barcode Trimming in this project, we use this bash script provided in the script directory:

 ./run_cutadapt_with_barcodes.sh <path_to_fastq1> <path_to_fastq2> <path_to_fasta_barcode_file>

Parameters

    <path_to_fastq1>: The path to the forward read FASTQ file containing the raw sequencing reads.
    <path_to_fastq2>: The path to the reverse read FASTQ file (for paired-end reads). If you are working with single-end reads, this should be an empty string.
    <path_to_fasta_barcode_file>: The path to the FASTA file containing the barcode sequences that need to be trimmed from the FASTQ files.

Expected Output

After running the script, the trimmed FASTQ files will be saved in the same directory as the input files, with the following naming convention:

    <input_fastq1>_trimmed.fastq: Trimmed forward read FASTQ file.
    <input_fastq2>_trimmed.fastq: Trimmed reverse read FASTQ file (for paired-end reads).

These trimmed files will be ready for alignment and further analysis.

This process ensures that the barcodes, which are used to differentiate between individual cells, are correctly removed from the reads before alignment. This process directly impacts the accuracy of downstream analysis, including cell identification and gene expression quantification.

# Read Alignment
After indexing the reference genome, the next crucial step in the single-cell RNA sequencing analysis pipeline is to align the raw sequencing reads to the indexed genome. This step maps the reads to their respective locations in the reference genome, which is necessary for downstream analysis such as gene expression quantification. To perform Read Alignment in this project, we use this bash script provided in the script directory:

./align_reads_modified_for_thread.sh <reference_genome_index_base> <trimmed_fastq_R1> <trimmed_fastq_R2> <num_threads>

Parameters

    <reference_genome_index_base>: The base name of the indexed reference genome files. This is the prefix used for the Bowtie2 index files generated during the genome indexing step.
    <trimmed_fastq_R1>: The path to the trimmed forward read FASTQ file. This file should be the result of the barcode trimming step.
    <trimmed_fastq_R2>: The path to the trimmed reverse read FASTQ file (for paired-end reads). This should also be the result of the barcode trimming step. If you are working with single-end reads, this should be an empty string.
    <num_threads>: The number of threads to use for parallel processing. This helps speed up the alignment process by utilizing multiple CPU cores.

Expected Output

After running the script, you will obtain the following files:

    <reference_genome_index_base>.sam: The SAM file containing the aligned reads. This file will be used for further processing and analysis.

The SAM file will need to be processed into BAM format in the next step of the pipeline.

# Processing BAM Files
After aligning the sequencing reads to the reference genome, the resulting SAM files must be processed into BAM files, which are more efficient for storage and analysis. This step includes conversion, sorting, and indexing the aligned reads.
Here's the script for it:

./process_bam.sh <input_sam> <output_prefix>

Parameters

    <input_sam>: The path to the SAM file containing the aligned reads. This file is the output from the read alignment step.
    <output_prefix>: The prefix for the output BAM files. The script will generate a sorted BAM file and an index file with this prefix.

Expected Output

After running the script, you will obtain the following files:

    <output_prefix>.bam: The sorted BAM file containing the aligned reads. This file is ready for downstream analysis.
    <output_prefix>.bam.bai: The index file for the sorted BAM file, which is required for efficient data retrieval and visualization.

Importance

Processing BAM files is a crucial step in single-cell RNA-seq analysis. It ensures that the data is organized and accessible for downstream analyses, such as identifying gene expression patterns and visualizing alignments in genome browsers.

Properly sorted and indexed BAM files allow for:

    Efficient data retrieval during analysis.
    Accurate visualization of alignment data.
    Seamless integration with tools that require indexed BAM files.

By automating these steps, the process_bam.sh script simplifies the preparation of your alignment data, ensuring that it's ready for the next phase of your analysis.

# Generating Count Matrix

Once the BAM files are sorted and indexed, the next critical step in single-cell RNA-seq analysis is to quantify gene expression levels. This is done by generating a count matrix, where each cell (or sample) is represented by the number of reads mapped to each gene. To perform Generating Count Matrix in this project, we use this bash script provided in the script directory:

./generate_count_matrix.sh <GTF_FILE> <BAM_FILES> <NUM_THREADS>

Parameters

    <GTF_FILE>: The path to the GTF file that provides gene annotations. This file is used to map reads to gene features and generate the count matrix.
    <BAM_FILES>: A comma-separated list of paths to the BAM files that contain the aligned reads. Ensure that all BAM files are sorted and indexed before running this script.
    <NUM_THREADS>: The number of threads to use for parallel processing. This helps speed up the count matrix generation by utilizing multiple CPU cores.

Expected Output

After running the script, you will obtain:

    <output_prefix>_counts.txt: A tab-separated count matrix file where rows represent genes, columns represent samples (cells), and values represent the number of reads mapped to each gene.

This count matrix is a crucial input for subsequent analysis steps, such as normalization and dimensionality reduction.

# Combining Count Matrices Script
The combine_matrices.R script is used to merge two count matrices into a single combined matrix. This process is essential for aligning datasets before performing joint analyses or comparisons, especially in single-cell RNA sequencing (scRNA-seq) studies.
Prerequisites

    R: Ensure R is installed on your system.
    Input Files: Two tab-separated count matrices of the Plasmodium falciparum and Plasmodium berghei

To perform Combining Count Matrices in this project, we use this bash script provided in the script directory:

Rscript combine_matrices.R --file1 counts1.txt --file2 counts2.txt --output combined_counts.txt

Parameters

    --file1 <counts1.txt>: The path to the  Plasmodium Falciparium count matrix file. This file contains count data for the Plasmodium Falciparium.
    --file2 <counts2.txt>: The path to the Plasmodium berghei count matrix file. This file contains count data for the Plasmodium berghei.
    --output <combined_counts.txt>: The path to the output file where the combined count matrix will be saved.

# Normalizing Count Data Script
The normalize_counts.R script is designed to process a count matrix by ensuring all data columns are numeric, normalizing the counts using DESeq2, and saving the results to a CSV file. This is crucial for preparing count data for downstream analysis in single-cell RNA sequencing (scRNA-seq) studies.
Prerequisites

    R: Ensure R and the DESeq2 package are installed on your system.
    Input File: the combined count matrix file from the above step

To perform GNormalizing Count Data in this project, we use this bash script provided in the script directory:

./normalize_combine_counts.R combined_counts.txt

Parameters

    <combined_counts.txt>: The path to the combined count matrix file. This file should contain the count data for all the Plasmodium falciparium and berghei samples and genes.
Expected Output

After running the script, you will obtain:

    combined_counts_normalized.csv: A CSV file containing the normalized count data. The normalization is performed using the DESeq2 package, preparing the data for further analysis, such as dimensionality reduction and clustering.
    
# Principal Component Analysis (PCA) with Seurat
This section of the pipeline describes how to perform Principal Component Analysis (PCA) on single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. PCA is a crucial step in the dimensionality reduction process, helping to identify the most significant sources of variation in your data. To perform Principal Component Analysis (PCA) in this project, we use this bash script provided in the script directory:

Usage: Rscript pca_analysis.R <path_to_count_data.csv>

Parameters

    <path_to_count_data.csv>: The path to the CSV file containing the normalized count matrix. This file should be the output from the normalization step and will be used as input for the PCA analysis.

Output

    PCA Visualization: A PCA plot that visually represents the main components of variation within your scRNA-seq data.
    Scree Plot: An elbow plot that helps you decide how many PCs to retain.
Plasmodium Falciparium PCA:
<img width="197" alt="image" src="https://github.com/user-attachments/assets/55a00492-6c89-46f3-bc41-21c50f99fade">


Plasmodium berghei PCA:
<img width="197" alt="image" src="https://github.com/user-attachments/assets/e7831b32-7d41-492a-85ce-23a137c98e9e">



This script is a key component of the preprocessing and analysis pipeline for scRNA-seq data, setting the stage for further clustering, visualization, and interpretation of single-cell data.






