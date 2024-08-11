# Load necessary libraries
library(Seurat)

# Define the path to the count matrix
count_matrix_path <- "/home/mikirium/scRNA-seq/output_counts/count_matrix.txt"

# Read the count matrix
count_matrix <- read.table(count_matrix_path, header = TRUE, row.names = 1, sep = "\t")

# Check for non-numeric columns and remove them
numeric_cols <- sapply(count_matrix, is.numeric)
count_matrix <- count_matrix[, numeric_cols]

# Extract numeric counts columns and ensure they are in matrix format
counts <- as.matrix(count_matrix)

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = counts)

# Filter cells to include only those with non-zero counts
seurat_object <- subset(seurat_object, subset = nCount_RNA > 0)

# Normalize the data
seurat_object <- NormalizeData(seurat_object)

# Find variable features
seurat_object <- FindVariableFeatures(seurat_object)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Save the processed Seurat object
saveRDS(seurat_object, file = "/home/mikirium/HTS_2024_Single_cell_RNA_seq_data_analysis/script/seurat_object_processed.rds")

cat("Processed Seurat object saved.\n")
