# estimate_variance_w_args

Total (across a specified set of clusters) of the sample variance of the
log odds of clusters membership across all samples within each of those
clusters

## Usage

``` r
estimate_variance_w_args(
  pars,
  counts,
  total_cells,
  optim_clusters,
  centering_matrix = NULL
)
```

## Arguments

- pars:

  A numeric vector of length equal to no. of samples, representing their
  corresponding normalization factors

- counts:

  A cell count matrix with number of columns representing the number of
  samples and number of rows representing the number of cell-types or
  clusters

- total_cells:

  An integer vector representing the total number of cells per sample

- optim_clusters:

  An integer vector representing the indices of the clusters over which
  the total variance is to be calculated

- centering_matrix:

  Optional precomputed centering matrix (nsamp x nsamp). If NULL,
  computed internally.

## Value

A numeric value representing the sum of variances of log odds of cluster
membership across all clusters specified in optim_clusters
