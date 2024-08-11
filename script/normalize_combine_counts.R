#!/usr/bin/env Rscript

# Install required packages if not already installed
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("DESeq2")
}

# Load necessary libraries
library(optparse)
library(DESeq2)

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path to the combined count matrix file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "normalized_counts.csv", help = "Output file for normalized counts", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input file path is provided
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file path must be supplied.", call. = FALSE)
}

# Read the count matrix
countData <- read.table(opt$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Ensure the count matrix is numeric and convert to integer
countData[] <- lapply(countData, function(x) {
  if (is.factor(x)) {
    as.integer(as.character(x))
  } else if (is.numeric(x)) {
    as.integer(x)
  } else {
    as.integer(as.character(replace(x, !is.numeric(x), 0)))
  }
})

if (!all(sapply(countData, is.numeric))) {
  stop("All columns in the count matrix must be numeric.", call. = FALSE)
}

# Define sample information and ensure unique row names
sample_names <- make.unique(colnames(countData))
colnames(countData) <- sample_names
colData <- data.frame(
  condition = factor(rep("Sample", length(sample_names))),
  row.names = sample_names
)

# Create DESeq2 dataset with design ~ 1
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 1)

# Run normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts to file
write.csv(as.data.frame(normalized_counts), file = opt$output, row.names = TRUE, quote = FALSE)

cat("Normalization completed. Normalized counts are in", opt$output, "\n")
