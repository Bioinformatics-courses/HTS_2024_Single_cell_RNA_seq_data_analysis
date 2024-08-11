#!/usr/bin/env Rscript

# Load necessary libraries
library(Seurat)

# Function to display usage message
usage <- function() {
    cat("Usage: Rscript pca_analysis.R <path_to_count_data.csv>\n")
    quit(status = 1)
}

# Ensure arguments are passed
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    usage()
}

# Path to the count data file
count_data_path <- args[1]

# Load count data
count_data <- read.csv(count_data_path, row.names = 1)

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = count_data)

# Perform normalization
seurat_object <- NormalizeData(seurat_object)

# Find variable features (genes) for PCA
seurat_object <- FindVariableFeatures(seurat_object)

# Scale data (this step is necessary for PCA)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(seurat_object))

# Determine the maximum number of principal components to compute
n_cells <- ncol(seurat_object)
max_pcs <- min(10, n_cells - 1)  # Compute up to 10 PCs or n_cells-1, whichever is smaller

# Perform PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object), npcs = max_pcs)

# Visualize PCA results
DimPlot(seurat_object, reduction = "pca")

# View PCA results
print(seurat_object[["pca"]])

# Plot a scree plot to determine the number of components to keep
ElbowPlot(seurat_object, reduction = "pca")
