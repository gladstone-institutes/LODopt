## Generate mock_seurat_data for the Seurat workflow vignette
##
## This script creates a small simulated dataset that mimics what a user would
## have after a typical scRNA-seq analysis: per-cell metadata (cell type, sample,
## group) and a gene expression matrix.

library(LODopt)
library(Matrix)

set.seed(42)

# --- Simulation parameters ---
K <- 6 # cell types
nsamp <- 12 # samples (6 per group)

# Proportions and variability
alpha <- c(3, 5, 8, 2, 4, 6)
props <- alpha / sum(alpha)
size <- rep(10, K)

# Clusters 1 and 4 change between groups
change_mean <- rep(1, K)
change_mean[c(1, 4)] <- c(0.3, 2.5)

# Simulate cluster x sample counts (tissue-level then downsampled)
sim <- simulate_cellCounts_fromTissue(
  props = props,
  nsamp = nsamp,
  size = size,
  depth = 1e8, # smaller depth for a compact dataset
  change_mean = change_mean
)

counts_mat <- sim$counts # K x nsamp integer matrix
cell_type_names <- c("B_cells", "T_cells_CD4", "T_cells_CD8",
                     "Monocytes", "NK_cells", "Dendritic_cells")
rownames(counts_mat) <- cell_type_names
colnames(counts_mat) <- paste0("sample_", seq_len(nsamp))

# --- Expand to per-cell data ---
cell_types <- character()
sample_ids <- character()
for (i in seq_len(K)) {
  for (j in seq_len(nsamp)) {
    n <- counts_mat[i, j]
    if (n > 0) {
      cell_types <- c(cell_types, rep(cell_type_names[i], n))
      sample_ids <- c(sample_ids, rep(colnames(counts_mat)[j], n))
    }
  }
}

n_cells <- length(cell_types)
barcodes <- paste0("cell_", seq_len(n_cells))

# Sample-level group assignment
group_map <- setNames(
  rep(c("control", "treated"), each = nsamp / 2),
  colnames(counts_mat)
)

cell_metadata <- data.frame(
  cell_type = cell_types,
  sample_id = sample_ids,
  group = group_map[sample_ids],
  row.names = barcodes,
  stringsAsFactors = FALSE
)

# --- Small random gene expression matrix (sparse) ---
n_genes <- 50
gene_names <- paste0("Gene", seq_len(n_genes))

# Sparse matrix: ~10% non-zero entries
counts_matrix <- Matrix::rsparsematrix(
  nrow = n_genes, ncol = n_cells,
  density = 0.10
)
# Make all values positive integers (typical for UMI counts)
counts_matrix@x <- abs(round(counts_matrix@x * 5)) + 1
dimnames(counts_matrix) <- list(gene_names, barcodes)

# --- Package as a list ---
mock_seurat_data <- list(
  counts_matrix = counts_matrix,
  cell_metadata = cell_metadata
)

usethis::use_data(mock_seurat_data, overwrite = TRUE)
