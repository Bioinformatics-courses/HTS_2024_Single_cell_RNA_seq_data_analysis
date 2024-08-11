#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-f", "--file1"), type = "character", default = NULL, help = "Path to first count matrix file", metavar = "character"),
  make_option(c("-s", "--file2"), type = "character", default = NULL, help = "Path to second count matrix file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "combined_counts.txt", help = "Output file for combined counts", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if both file paths are provided
if (is.null(opt$file1) | is.null(opt$file2)) {
  print_help(opt_parser)
  stop("Both file paths must be supplied.", call. = FALSE)
}

# Read the count matrices
countData1 <- read.table(opt$file1, header = TRUE, row.names = 1, sep = "\t")
countData2 <- read.table(opt$file2, header = TRUE, row.names = 1, sep = "\t")

# Ensure both datasets have the same features (genes)
all_genes <- union(rownames(countData1), rownames(countData2))

# Add missing features (genes) with zeros
countData1 <- countData1[all_genes, ]
countData2 <- countData2[all_genes, ]

countData1[is.na(countData1)] <- 0
countData2[is.na(countData2)] <- 0

# Combine the count matrices
combined_counts <- cbind(countData1, countData2)

# Write the combined counts to a file
write.table(combined_counts, file = opt$output, sep = "\t", quote = FALSE, col.names = NA)

cat("Combined count matrix saved to", opt$output, "\n")
