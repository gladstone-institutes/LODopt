# Introduction to LODopt

## Introduction

The `LODopt` package provides tools for analyzing associations between
cell type abundances and experimental conditions in single-cell RNA-seq
data.

First, let’s load the necessary libraries:

``` r
library(LODopt)
library(dplyr)
library(magrittr)
library(SummarizedExperiment)
```

## Quick Start

### Example Data

We simulate data for 30 samples with 25 clusters, where 4 clusters (c0,
c2, c7 and c14) have changed their abundance.

``` r
# Create example count matrix
# Create test data
set.seed(123)
# No of clusters in single cell dataset
K <- 25
# No of cells in tissue
depth <- 1e9
# Total number of samples
nsamp <- 30
alpha <- 10^runif(K, min = log10(0.5), max = log10(10))
p <- dirmult::rdirichlet(alpha = alpha) |> sort()
p <- p[p > 0.001]
p <- p / sum(p)
size <- rep(10, length(p))
change_mean <- rep(1, length(p))
## changed cluster indices
change_mean[c(1, 3, 8, 15)] <- c(0.2, 2, 0.2, 2)
# Simulate counts
counts_res <- simulate_cellCounts_fromTissue(
  props = p,
  nsamp = nsamp,
  depth = depth,
  size = size,
  change_mean = change_mean
)
```

### Set up the required SummarizedExperiment object

``` r
counts <- counts_res$counts

pheno_data <- data.frame(
  sampleID = paste0("S", 1:30),
  groupid = c(rep("group0", 15), rep("group1", 15))
) %>%
  tibble::column_to_rownames("sampleID")

model_formula <- "groupid"
cellcomp_se <- SummarizedExperiment(
  assays = list(counts = counts),
  colData = pheno_data,
  metadata = list(
    modelFormula = model_formula,
    coef_of_interest_index = 2,
    reference_levels_of_variables = list(c("groupid", "group0")),
    random_seed = 123456,
    unchanged_cluster_indices = NULL
  )
)
```

### Running Association Analysis

``` r
cellcomp_res <- logodds_optimized_normFactors(cellcomp_se)
#> Running iteration no. 1
#> Optimization did not converge, retrying (1/5)
#> boundary (singular) fit: see help('isSingular')
#> Running iteration no. 2
#> Found stable solution
```

### Summarizing Results

The estimate column has the log odds ratio estimates for the
cluster-specific association with groupid while the
estimates_significance has the corresponding raw p-values for the
corresponding null hypothesis that the log odds ratio is equal to 0.

``` r
print(cellcomp_res$results)
#>    cluster_id    comparison   estimates estimates_significance adjusted_pvalue
#> 1          c0 groupidgroup1 -1.52150314           1.195556e-17    1.494446e-16
#> 2          c1 groupidgroup1 -0.17537620           2.076247e-01    6.663173e-01
#> 3          c2 groupidgroup1  0.79158810           1.723059e-09    1.076912e-08
#> 4          c3 groupidgroup1  0.04338405           7.257719e-01    8.072362e-01
#> 5          c4 groupidgroup1 -0.04461750           6.811076e-01    8.072362e-01
#> 6          c5 groupidgroup1 -0.13531225           2.398742e-01    6.663173e-01
#> 7          c6 groupidgroup1 -0.03062671           7.410105e-01    8.072362e-01
#> 8          c7 groupidgroup1 -1.75974180           5.271561e-50    1.317890e-48
#> 9          c8 groupidgroup1  0.12793866           1.882550e-01    6.663173e-01
#> 10         c9 groupidgroup1  0.14299432           2.226219e-01    6.663173e-01
#> 11        c10 groupidgroup1 -0.09383282           5.157426e-01    7.584450e-01
#> 12        c11 groupidgroup1  0.09462389           3.953366e-01    7.584450e-01
#> 13        c12 groupidgroup1  0.02717145           8.180407e-01    8.521257e-01
#> 14        c13 groupidgroup1 -0.08186446           4.474527e-01    7.584450e-01
#> 15        c14 groupidgroup1  0.68476298           7.443181e-12    6.202651e-11
#> 16        c15 groupidgroup1 -0.04887438           7.426573e-01    8.072362e-01
#> 17        c16 groupidgroup1  0.07879893           4.974695e-01    7.584450e-01
#> 18        c17 groupidgroup1 -0.12203460           3.191846e-01    7.584450e-01
#> 19        c18 groupidgroup1 -0.14284924           2.136904e-01    6.663173e-01
#> 20        c19 groupidgroup1  0.01474030           8.988846e-01    8.988846e-01
#> 21        c20 groupidgroup1  0.08825082           4.824833e-01    7.584450e-01
#> 22        c21 groupidgroup1  0.05926056           4.994481e-01    7.584450e-01
#> 23        c22 groupidgroup1  0.09703444           4.711025e-01    7.584450e-01
#> 24        c23 groupidgroup1 -0.06478133           5.823189e-01    7.662091e-01
#> 25        c24 groupidgroup1  0.08662988           5.729970e-01    7.662091e-01
```

## Conclusion

The `LODopt` package provides a simple interface for testing
associations between cell type abundances and experimental conditions,
with support for multiple statistical methods and comprehensive result
visualization.
