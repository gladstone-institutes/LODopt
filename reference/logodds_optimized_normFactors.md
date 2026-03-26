# logodds_optimized_normFactors

The primary LODopt method

## Usage

``` r
logodds_optimized_normFactors(cellcomp_se, verbose = TRUE)
```

## Arguments

- cellcomp_se:

  a SummarizedExperiment object that includes

  - counts: an assay named counts with the number of cells from each
    sample in each of the clusters. The counts is a data frame with
    rownames set to the cluster names (without spaces or special
    characters) and column names representing the sample names

  - data for the association analyses assigned to colData as a data
    frame with row names the same as the column names of the counts data
    frame

  - Five parameters assigned to the metadata slot

    - modelFormula representing the model to be tested using the
      variable in the colData

    - coef_of_interest_index representing the index of the coefficient
      of interest to be tested in the model

    - reference_levels_of_variables is a list of vectors.Each vector
      will have two elements - the name of a categorical variable (a
      column name in colData) and the reference level to set this
      variable to

    - random_seed integer random seed

    - unchanged_cluster_indices cluster indices of the clusters known
      not to change, set to NULL if not known

- verbose:

  Logical; if TRUE (default), prints progress messages

## Value

A list with two elements

- results: data frame of results of differential analyses with 4
  columns - cluster_id, comparison, estimates (log odds ratio) and
  estimate_significance (pvalue)

- optim_factor: a data frame with columns sampleid and optim.norm.factor

## Examples

``` r
# Create example data
set.seed(123)
# No of clusters in single cell dataset
K <- 25
nsamp <- 30
alpha <- 10^runif(K, min = log10(0.5), max = log10(10))
p <- dirmult::rdirichlet(alpha = alpha) |> sort()
p <- p[p > 0.001]
p <- p / sum(p)
size <- rep(10, length(p))
change_mean <- rep(1, length(p))
## clusters 1, 3, 8 and 15 are changed. Clusters 1 and 8 reduce in abundance while clusters 3 and 15 increase
change_mean[c(1, 3, 8, 15)] <- c(0.2, 2, 0.2, 2)
depth <- 1e9
# Simulate counts
counts_res <- simulate_cellCounts_fromTissue(
  props = p,
  nsamp = nsamp,
  depth = depth,
  size = size,
  change_mean = change_mean
)
counts <- counts_res$counts

pheno_data <- data.frame(
  sampleID = paste0("S", 1:30),
  groupid = c(rep("group0", 15), rep("group1", 15))
)
require(magrittr)
#> Loading required package: magrittr
pheno_data %<>% tibble::column_to_rownames("sampleID")
require(SummarizedExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> 
#> Attaching package: ‘GenomicRanges’
#> The following object is masked from ‘package:magrittr’:
#> 
#>     subtract
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
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
cellcomp_res <- logodds_optimized_normFactors(cellcomp_se)
#> Running iteration no. 1
#> boundary (singular) fit: see help('isSingular')
#> Running iteration no. 2
#> Running iteration no. 3
#> Found stable solution
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
