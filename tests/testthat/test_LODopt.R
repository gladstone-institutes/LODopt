library(testthat)
library(LODopt)

test_that("simulate_cellCounts_fromTissue works with basic input", {
  # Create test data
  set.seed(123)
   # No of clusters in single cell dataset
   K=25
   nsamp = 30
   alpha = 10^runif(K, min=log10(0.5), max = log10(10))
   p <- dirmult::rdirichlet(alpha = alpha) |> sort()
   p <- p[p > 0.001]
   p <- p/sum(p)
   size <- rep(10, length(p))
   change_mean = rep(1, length(p))
   change_mean[c(1,3,8,15)] = c(0.2, 2, 0.2, 2)
   depth = 1e9
   # Simulate counts
   counts_res <- simulate_cellCounts_fromTissue(props=p,nsamp=nsamp,depth=depth, size = size, change_mean = change_mean)
   expect_is(counts_res, "list")
   expect_equal(ncol(counts_res$counts), nsamp)
  # # Assertions
  # expect_is(results, "data.frame")
  # expect_equal(nrow(results), 5)
  # expect_true(all(c("celltype", "coefficient", "p_value", "adj_p_value") %in% colnames(results)))
  # expect_true(all(!is.na(results$celltype)))
})


test_that("logodds_optimized_normFactors works with basic input", {
  # Create test data
  set.seed(123)
  # No of clusters in single cell dataset
  K=25
  nsamp = 30
  alpha = 10^runif(K, min=log10(0.5), max = log10(10))
  p <- dirmult::rdirichlet(alpha = alpha) |> sort()
  p <- p[p > 0.001]
  p <- p/sum(p)
  size <- rep(10, length(p))
  change_mean = rep(1, length(p))
  change_mean[c(1,3,8,15)] = c(0.2, 2, 0.2, 2)
  depth = 1e9
  # Simulate counts
  counts_res <- simulate_cellCounts_fromTissue(props=p,nsamp=nsamp,depth=depth, size = size, change_mean = change_mean)
  counts <- counts_res$counts

  pheno_data <- data.frame(sampleID = paste0("S", 1:30),
                           groupid = c(rep("group0", 15), rep("group1", 15)))
  pheno_data %<>% tibble::column_to_rownames("sampleID")



  require(SummarizedExperiment)
  model_formula <- "groupid"
  cellcomp_se <- SummarizedExperiment(assays = list(counts=counts),
                                      colData = pheno_data,
                                      metadata = list(modelFormula = model_formula,
                                                      coef_of_interest_index = 2,
                                                      reference_levels_of_variables = list(c("groupid", "group0")),
                                                      random_seed = 123456,
                                                      unchanged_cluster_indices = NULL))
  results <- logodds_optimized_normFactors(cellcomp_se)
  expect_is(results, "list")
  # # Assertions
  expect_is(results$results, "data.frame")
  expect_equal(nrow(results$results), nrow(counts))
  expect_true(all(c("cluster_id", "comparison", "estimates", "estimates_significance") %in% colnames(results$results)))
  # expect_true(all(!is.na(results$celltype)))
})

# test_that("celltype_association handles invalid input", {
#   # Test with mismatched dimensions
#   count_mat <- matrix(1:20, nrow = 4, ncol = 5)
#   pheno <- data.frame(condition = factor(rep(c("A", "B"), each = 2)))
#   design_mat <- model.matrix(~ condition, data = pheno)
#
#   expect_error(celltype_association(count_mat, pheno, design_mat))
# })
#
# test_that("create_design_matrix works correctly", {
#   pheno <- data.frame(
#     condition = factor(c("A", "B", "A", "B")),
#     age = c(25, 30, 35, 40)
#   )
#
#   design_mat <- create_design_matrix(pheno, ~ condition + age)
#
#   expect_is(design_mat, "matrix")
#   expect_equal(nrow(design_mat), nrow(pheno))
#   expect_true(ncol(design_mat) >= 2)
# })
