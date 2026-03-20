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
#> boundary (singular) fit: see help('isSingular')
#> Running iteration no. 2
#> Running iteration no. 3
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
#> 1          c0 groupidgroup1 -1.52085782           1.084964e-17    1.356205e-16
#> 2          c1 groupidgroup1 -0.17332361           2.126613e-01    6.943158e-01
#> 3          c2 groupidgroup1  0.79403652           1.988380e-09    1.242738e-08
#> 4          c3 groupidgroup1  0.04600408           7.111394e-01    8.281567e-01
#> 5          c4 groupidgroup1 -0.04198163           6.982457e-01    8.281567e-01
#> 6          c5 groupidgroup1 -0.13280991           2.520703e-01    7.001954e-01
#> 7          c6 groupidgroup1 -0.02810406           7.619042e-01    8.281567e-01
#> 8          c7 groupidgroup1 -1.75694672           2.737954e-50    6.844884e-49
#> 9          c8 groupidgroup1  0.13052353           1.826849e-01    6.943158e-01
#> 10         c9 groupidgroup1  0.14558659           2.127333e-01    6.943158e-01
#> 11        c10 groupidgroup1 -0.09111349           5.277537e-01    7.761083e-01
#> 12        c11 groupidgroup1  0.09738981           3.801385e-01    7.535890e-01
#> 13        c12 groupidgroup1  0.02982580           8.004104e-01    8.337608e-01
#> 14        c13 groupidgroup1 -0.07925789           4.655635e-01    7.535890e-01
#> 15        c14 groupidgroup1  0.68747633           9.859540e-12    8.216283e-11
#> 16        c15 groupidgroup1 -0.04620725           7.566816e-01    8.281567e-01
#> 17        c16 groupidgroup1  0.08169586           4.822969e-01    7.535890e-01
#> 18        c17 groupidgroup1 -0.11911308           3.284509e-01    7.535890e-01
#> 19        c18 groupidgroup1 -0.13997526           2.221811e-01    6.943158e-01
#> 20        c19 groupidgroup1  0.01767191           8.789407e-01    8.789407e-01
#> 21        c20 groupidgroup1  0.09132718           4.686408e-01    7.535890e-01
#> 22        c21 groupidgroup1  0.06219109           4.732152e-01    7.535890e-01
#> 23        c22 groupidgroup1  0.10016432           4.565396e-01    7.535890e-01
#> 24        c23 groupidgroup1 -0.06190911           5.998054e-01    7.892177e-01
#> 25        c24 groupidgroup1  0.08963392           5.594980e-01    7.770806e-01
```

## Conclusion

The `LODopt` package provides a simple interface for testing
associations between cell type abundances and experimental conditions,
with support for multiple statistical methods and comprehensive result
visualization.
