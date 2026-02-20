library(testthat)
library(LODopt)

# Helper: create a standard test SummarizedExperiment
make_test_se <- function(nsamp = 30, K = 25, seed = 123) {
  set.seed(seed)
  alpha <- 10^runif(K, min = log10(0.5), max = log10(10))
  p <- dirmult::rdirichlet(alpha = alpha) |> sort()
  p <- p[p > 0.001]
  p <- p / sum(p)
  size <- rep(10, length(p))
  change_mean <- rep(1, length(p))
  change_mean[c(1, 3, 8, 15)] <- c(0.2, 2, 0.2, 2)
  depth <- 1e9

  counts_res <- simulate_cellCounts_fromTissue(
    props = p, nsamp = nsamp, depth = depth,
    size = size, change_mean = change_mean
  )
  counts <- counts_res$counts

  pheno_data <- data.frame(
    sampleID = paste0("S", 1:nsamp),
    groupid = c(rep("group0", nsamp / 2), rep("group1", nsamp / 2))
  )
  pheno_data <- tibble::column_to_rownames(pheno_data, "sampleID")

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = pheno_data,
    metadata = list(
      modelFormula = "groupid",
      coef_of_interest_index = 2,
      reference_levels_of_variables = list(c("groupid", "group0")),
      random_seed = 123456,
      unchanged_cluster_indices = NULL
    )
  )
}

# =============================================================================
# Simulation tests
# =============================================================================

test_that("simulate_cellCounts_fromTissue works with basic input", {
  set.seed(123)
  K <- 25
  nsamp <- 30
  alpha <- 10^runif(K, min = log10(0.5), max = log10(10))
  p <- dirmult::rdirichlet(alpha = alpha) |> sort()
  p <- p[p > 0.001]
  p <- p / sum(p)
  size <- rep(10, length(p))
  change_mean <- rep(1, length(p))
  change_mean[c(1, 3, 8, 15)] <- c(0.2, 2, 0.2, 2)
  depth <- 1e9

  counts_res <- simulate_cellCounts_fromTissue(
    props = p, nsamp = nsamp, depth = depth,
    size = size, change_mean = change_mean
  )
  expect_type(counts_res, "list")
  expect_equal(ncol(counts_res$counts), nsamp)
  expect_true(all(counts_res$counts >= 0))
  expect_equal(length(counts_res$theoretical_log_odds_ratios), length(p))
})

# =============================================================================
# Full pipeline test
# =============================================================================

test_that("logodds_optimized_normFactors works with basic input", {
  cellcomp_se <- make_test_se()
  results <- logodds_optimized_normFactors(cellcomp_se, verbose = FALSE)

  expect_type(results, "list")
  expect_s3_class(results$results, "data.frame")
  expect_equal(nrow(results$results), nrow(cellcomp_se))
  expect_true(all(c(
    "cluster_id", "comparison", "estimates",
    "estimates_significance", "adjusted_pvalue"
  ) %in%
    colnames(results$results)))
  # adjusted p-values >= raw p-values
  non_na <- !is.na(results$results$estimates_significance)
  expect_true(all(
    results$results$adjusted_pvalue[non_na] >=
      results$results$estimates_significance[non_na] - 1e-10
  ))
  # optim_factor is a data frame
  expect_s3_class(results$optim_factor, "data.frame")
  expect_true("sampleid" %in% colnames(results$optim_factor))
  expect_true("optim.norm.factor" %in% colnames(results$optim_factor))
})

# =============================================================================
# Input validation tests
# =============================================================================

test_that("logodds_optimized_normFactors rejects non-SE input", {
  expect_error(
    logodds_optimized_normFactors(data.frame(a = 1)),
    "must be a SummarizedExperiment"
  )
})

test_that("logodds_optimized_normFactors rejects missing counts assay", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(data = matrix(1:4, 2, 2)),
    metadata = list(
      modelFormula = "x", coef_of_interest_index = 2,
      reference_levels_of_variables = NULL,
      random_seed = 1, unchanged_cluster_indices = NULL
    )
  )
  expect_error(
    logodds_optimized_normFactors(se),
    "assay named 'counts'"
  )
})

test_that("logodds_optimized_normFactors rejects missing metadata fields", {
  counts <- matrix(1:4, 2, 2)
  colnames(counts) <- c("S1", "S2")
  rownames(counts) <- c("c1", "c2")
  pheno <- data.frame(row.names = c("S1", "S2"), group = c("a", "b"))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = pheno,
    metadata = list(modelFormula = "group")
  )
  expect_error(
    logodds_optimized_normFactors(se),
    "Missing required metadata fields"
  )
})

test_that("logodds_optimized_normFactors rejects negative counts", {
  counts <- matrix(c(-1, 2, 3, 4), 2, 2)
  colnames(counts) <- c("S1", "S2")
  rownames(counts) <- c("c1", "c2")
  pheno <- data.frame(row.names = c("S1", "S2"), group = c("a", "b"))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = pheno,
    metadata = list(
      modelFormula = "group", coef_of_interest_index = 2,
      reference_levels_of_variables = NULL,
      random_seed = 1, unchanged_cluster_indices = NULL
    )
  )
  expect_error(
    logodds_optimized_normFactors(se),
    "non-negative"
  )
})

# =============================================================================
# Name-based argument tests
# =============================================================================

test_that("unchanged_cluster_indices accepts character names", {
  cellcomp_se <- make_test_se()
  cluster_names <- rownames(cellcomp_se)
  # Use first 5 cluster names
  cellcomp_se@metadata$unchanged_cluster_indices <- cluster_names[1:5]

  results <- logodds_optimized_normFactors(cellcomp_se, verbose = FALSE)
  expect_type(results, "list")
  expect_s3_class(results$results, "data.frame")
})

test_that("unchanged_cluster_indices rejects invalid names", {
  cellcomp_se <- make_test_se()
  cellcomp_se@metadata$unchanged_cluster_indices <- c("nonexistent_cluster")

  expect_error(
    logodds_optimized_normFactors(cellcomp_se, verbose = FALSE),
    "not found"
  )
})

test_that("coef_of_interest_index accepts character name", {
  cellcomp_se <- make_test_se()
  cellcomp_se@metadata$coef_of_interest_index <- "groupidgroup1"

  results <- logodds_optimized_normFactors(cellcomp_se, verbose = FALSE)
  expect_type(results, "list")
  expect_s3_class(results$results, "data.frame")
})

# =============================================================================
# Seurat conversion tests
# =============================================================================

test_that("prepare_from_seurat works with mock Seurat object", {
  skip_if_not_installed("SeuratObject")

  # Create a minimal Seurat-like object
  n_cells <- 200
  meta <- data.frame(
    cell_type = sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE),
    sample_id = sample(paste0("S", 1:6), n_cells, replace = TRUE),
    group = NA_character_,
    stringsAsFactors = FALSE
  )
  # Assign group based on sample
  meta$group <- ifelse(meta$sample_id %in% c("S1", "S2", "S3"), "ctrl", "treat")
  rownames(meta) <- paste0("cell_", seq_len(n_cells))

  # Create a minimal count matrix
  counts_mat <- matrix(
    rpois(10 * n_cells, lambda = 5),
    nrow = 10, ncol = n_cells,
    dimnames = list(paste0("Gene", 1:10), rownames(meta))
  )

  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts_mat, meta.data = meta)

  se <- prepare_from_seurat(
    seurat_obj,
    cluster_col = "cell_type",
    sample_col = "sample_id",
    pheno_cols = "group",
    model_formula = "group",
    coef_of_interest = 2,
    reference_levels = list(group = "ctrl"),
    random_seed = 42
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_true("counts" %in% SummarizedExperiment::assayNames(se))
  expect_equal(ncol(se), 6) # 6 samples
  expect_equal(nrow(se), 3) # 3 cell types
  expect_equal(se@metadata$modelFormula, "group")
})

# =============================================================================
# Error handling test
# =============================================================================

test_that("GLMM failure for degenerate cluster returns NA, not error", {
  cellcomp_se <- make_test_se()
  # Set one cluster to all zeros - should cause GLMM issues
  counts <- SummarizedExperiment::assays(cellcomp_se)$counts
  counts[1, ] <- 0
  SummarizedExperiment::assays(cellcomp_se)$counts <- counts

  # Should not error, but the zero cluster may get NA
  results <- logodds_optimized_normFactors(cellcomp_se, verbose = FALSE)
  expect_type(results, "list")
  expect_s3_class(results$results, "data.frame")
})
