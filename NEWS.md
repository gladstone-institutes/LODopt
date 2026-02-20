# LODopt 1.1.0

## New features

* Added `prepare_from_seurat()` for converting Seurat objects directly to `SummarizedExperiment` input (#3).
* `coef_of_interest_index` and `unchanged_cluster_indices` now accept character names in addition to integer indices (#5).
* Added BH-adjusted p-values (`adjusted_pvalue` column) to results (#7).
* Optimization automatically retries with new starting values on convergence failure (#8).
* Added `verbose` parameter to `logodds_optimized_normFactors()` to control progress messages.
* Added input validation for `SummarizedExperiment` structure and metadata fields (#4).
* GLMM fitting is wrapped in `tryCatch` so a single cluster failure no longer aborts the full analysis.
* Added built-in `mock_seurat_data` example dataset and a Seurat workflow vignette.

## Bug fixes

* Fixed operator precedence bug in `simulate_cellCounts_fromTissue()` that could produce incorrect tissue-level counts.
* Fixed return value and documentation mismatches (#6).

## Performance

* Precompute centering matrix in `estimate_variance_w_args()` to avoid redundant allocation during optimization.

## Infrastructure

* Added CI/CD workflow (GitHub Actions R CMD check on macOS, Windows, Ubuntu).
* Added testthat test suite with broad coverage.
* Cleaned up DESCRIPTION dependencies.

# LODopt 1.0.0

* Initial release.
