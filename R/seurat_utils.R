#' Prepare LODopt input from a Seurat object
#'
#' Converts a Seurat object into a SummarizedExperiment ready for
#' \code{\link{logodds_optimized_normFactors}}.
#'
#' @param seurat_obj A Seurat object
#' @param cluster_col Column in meta.data with cluster/cell type labels
#' @param sample_col Column in meta.data with sample IDs
#' @param pheno_cols Character vector of columns in meta.data for phenotype data (must be sample-level, not cell-level)
#' @param model_formula Model formula string (RHS only, e.g., "groupid")
#' @param coef_of_interest Index or name of coefficient of interest
#' @param reference_levels Named list: list(variable = reference_level)
#' @param random_seed Integer random seed
#' @param unchanged_clusters NULL, or character/integer vector of unchanged clusters
#' @param sanitize_names Logical, replace spaces/special characters in cluster names (default TRUE)
#'
#' @return A SummarizedExperiment ready for \code{\link{logodds_optimized_normFactors}}
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
prepare_from_seurat <- function(seurat_obj, cluster_col, sample_col, pheno_cols,
                                model_formula, coef_of_interest,
                                reference_levels = NULL, random_seed = 42,
                                unchanged_clusters = NULL,
                                sanitize_names = TRUE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package required. Install with: install.packages('SeuratObject')")
  }

  # Extract metadata
  meta <- seurat_obj@meta.data

  # Cross-tabulate: samples (rows) x clusters (cols)
  count_table <- table(meta[[sample_col]], meta[[cluster_col]])
  counts <- t(as.data.frame.matrix(count_table)) # clusters x samples

  # Sanitize cluster names
  if (sanitize_names) {
    rownames(counts) <- gsub("[^A-Za-z0-9_.]", "_", rownames(counts))
  }

  # Build per-sample phenotype data
  pheno_data <- meta[, c(sample_col, pheno_cols), drop = FALSE]
  pheno_data <- pheno_data[!duplicated(pheno_data[[sample_col]]), , drop = FALSE]
  rownames(pheno_data) <- pheno_data[[sample_col]]
  pheno_data[[sample_col]] <- NULL
  pheno_data <- pheno_data[colnames(counts), , drop = FALSE]

  # Build reference_levels_of_variables in LODopt format
  ref_levels <- NULL
  if (!is.null(reference_levels)) {
    ref_levels <- lapply(names(reference_levels), function(v) c(v, reference_levels[[v]]))
  }

  # Construct SummarizedExperiment
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = pheno_data,
    metadata = list(
      modelFormula = model_formula,
      coef_of_interest_index = coef_of_interest,
      reference_levels_of_variables = ref_levels,
      random_seed = random_seed,
      unchanged_cluster_indices = unchanged_clusters
    )
  )
}
