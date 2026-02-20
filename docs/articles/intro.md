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
#> Optimization did not converge, retrying (1/5)
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
#> 1          c0 groupidgroup1 -1.51858161           1.545083e-17    1.931353e-16
#> 2          c1 groupidgroup1 -0.17351938           2.118907e-01    6.868771e-01
#> 3          c2 groupidgroup1  0.79359832           1.520303e-09    9.501892e-09
#> 4          c3 groupidgroup1  0.04533991           7.134683e-01    8.230478e-01
#> 5          c4 groupidgroup1 -0.04270437           6.942090e-01    8.230478e-01
#> 6          c5 groupidgroup1 -0.13340942           2.472757e-01    6.868771e-01
#> 7          c6 groupidgroup1 -0.02868183           7.572040e-01    8.230478e-01
#> 8          c7 groupidgroup1 -1.75753960           2.666082e-50    6.665206e-49
#> 9          c8 groupidgroup1  0.12982573           1.820380e-01    6.868771e-01
#> 10         c9 groupidgroup1  0.14486437           2.167701e-01    6.868771e-01
#> 11        c10 groupidgroup1 -0.09198216           5.240377e-01    7.706437e-01
#> 12        c11 groupidgroup1  0.09652848           3.885101e-01    7.626061e-01
#> 13        c12 groupidgroup1  0.02899070           8.057244e-01    8.392963e-01
#> 14        c13 groupidgroup1 -0.07997610           4.571285e-01    7.626061e-01
#> 15        c14 groupidgroup1  0.68668034           7.456407e-12    6.213672e-11
#> 16        c15 groupidgroup1 -0.04694615           7.519014e-01    8.230478e-01
#> 17        c16 groupidgroup1  0.08071082           4.880679e-01    7.626061e-01
#> 18        c17 groupidgroup1 -0.12016532           3.272104e-01    7.626061e-01
#> 19        c18 groupidgroup1 -0.14102995           2.199268e-01    6.868771e-01
#> 20        c19 groupidgroup1  0.01669206           8.850909e-01    8.850909e-01
#> 21        c20 groupidgroup1  0.09026214           4.735210e-01    7.626061e-01
#> 22        c21 groupidgroup1  0.06121635           4.860588e-01    7.626061e-01
#> 23        c22 groupidgroup1  0.09901983           4.624356e-01    7.626061e-01
#> 24        c23 groupidgroup1 -0.06289889           5.943226e-01    7.820034e-01
#> 25        c24 groupidgroup1  0.08879573           5.626717e-01    7.814885e-01
```

## Conclusion

The `LODopt` package provides a simple interface for testing
associations between cell type abundances and experimental conditions,
with support for multiple statistical methods and comprehensive result
visualization.
