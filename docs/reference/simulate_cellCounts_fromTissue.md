# simulate_cellCounts_fromTissue

Simulate single cell assay counts of cell-types derived from a given
tissue for each of nsamp samples starting by simulating tissue level
counts and then downsampling to obtain single cell assay counts

## Usage

``` r
simulate_cellCounts_fromTissue(
  props,
  nsamp,
  size = NULL,
  depth = 1e+09,
  change_mean
)
```

## Arguments

- props:

  A vector of proportions of cell-types in the tissue. Values should be
  positive and add up to 1.

- nsamp:

  An integer representing the number of samples

- size:

  A numeric negative binomial distribution parameter representing the
  biological coefficient of variation for cell counts in the tissue

- depth:

  Expected number of cells in tissue (default = 1e9)

- change_mean:

  a positive numeric vector of the same size as proportion representing
  the change in the proportion of individual cell-types. A value of 1
  implies no change, a value less than 1 implies depletion while values
  above 1 represent increase in the number of the associated cell-type

## Value

A list with the following elements:

- counts: count matrix with rows corresponding to the clusters and nsamp
  columns

- observed_log_odds_ratios: observed log odds ratios based on simulated
  tissue level counts

- theoretical_log_odds_ratios: theoretical log odds ratios based on
  provided parameters

- theoretical_log_absolute_abundance_ratios: theoretical log of absolute
  abundance ratios based on provided parameters

- observed_log_absolute_abundance_ratios: observed log of absolute
  abundance ratios based on simulated tissue level counts

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
change_mean[c(1, 3, 8, 15)] <- c(0.2, 2, 0.2, 2)
depth <- 1e9
# Simulate counts
counts_res <- simulate_cellCounts_fromTissue(props = p, nsamp = nsamp, depth = depth, size = size, change_mean = change_mean)
print(counts_res$counts[1:3, 1:3])
#>    S1 S2 S3
#> c0  7 10  6
#> c1 37 41 61
#> c2 28 29 48
```
