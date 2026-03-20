# Prepare LODopt input from a Seurat object

Converts a Seurat object into a SummarizedExperiment ready for
[`logodds_optimized_normFactors`](https://gladstone-institutes.github.io/LODopt/reference/logodds_optimized_normFactors.md).

## Usage

``` r
prepare_from_seurat(
  seurat_obj,
  cluster_col,
  sample_col,
  pheno_cols,
  model_formula,
  coef_of_interest,
  reference_levels = NULL,
  random_seed = 42,
  unchanged_clusters = NULL,
  sanitize_names = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cluster_col:

  Column in meta.data with cluster/cell type labels

- sample_col:

  Column in meta.data with sample IDs

- pheno_cols:

  Character vector of columns in meta.data for phenotype data (must be
  sample-level, not cell-level)

- model_formula:

  Model formula string (RHS only, e.g., "groupid")

- coef_of_interest:

  Index or name of coefficient of interest

- reference_levels:

  Named list: list(variable = reference_level)

- random_seed:

  Integer random seed

- unchanged_clusters:

  NULL, or character/integer vector of unchanged clusters

- sanitize_names:

  Logical, replace spaces/special characters in cluster names (default
  TRUE)

## Value

A SummarizedExperiment ready for
[`logodds_optimized_normFactors`](https://gladstone-institutes.github.io/LODopt/reference/logodds_optimized_normFactors.md)
