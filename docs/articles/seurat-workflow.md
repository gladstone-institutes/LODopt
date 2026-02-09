# Seurat Workflow for LODopt

## Introduction

If you have a **Seurat object** from your scRNA-seq analysis pipeline,
LODopt provides
[`prepare_from_seurat()`](https://gladstone-institutes.github.io/LODopt/reference/prepare_from_seurat.md)
to convert it directly into the `SummarizedExperiment` format that
[`logodds_optimized_normFactors()`](https://gladstone-institutes.github.io/LODopt/reference/logodds_optimized_normFactors.md)
expects. This vignette walks through the full workflow.

## Load packages

``` r
library(LODopt)
library(SeuratObject)
```

## Load built-in example data

LODopt ships with `mock_seurat_data`, a small simulated dataset with 6
cell types across 12 samples (6 control, 6 treated). Two cell types —
B_cells and Monocytes — have differential abundance between groups.

``` r
data("mock_seurat_data")
str(mock_seurat_data, max.level = 1)
#> List of 2
#>  $ counts_matrix:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ cell_metadata:'data.frame':   13279 obs. of  3 variables:
```

The list contains:

- **`counts_matrix`** — a sparse gene-by-cell expression matrix
- **`cell_metadata`** — a data frame with `cell_type`, `sample_id`, and
  `group` columns

## Create a Seurat object

This is the step where you would normally already have a Seurat object
from your own analysis. Here we construct one from the example data:

``` r
seurat_obj <- CreateSeuratObject(
  counts = mock_seurat_data$counts_matrix,
  meta.data = mock_seurat_data$cell_metadata
)
seurat_obj
#> An object of class Seurat 
#> 50 features across 13279 samples within 1 assay 
#> Active assay: RNA (50 features, 0 variable features)
#>  1 layer present: counts
```

Let’s check the metadata:

``` r
head(seurat_obj@meta.data)
#>        orig.ident nCount_RNA nFeature_RNA cell_type sample_id   group
#> cell_1       cell         32            8   B_cells  sample_1 control
#> cell_2       cell          3            2   B_cells  sample_1 control
#> cell_3       cell         24            5   B_cells  sample_1 control
#> cell_4       cell         34            6   B_cells  sample_1 control
#> cell_5       cell         32            7   B_cells  sample_1 control
#> cell_6       cell         21            4   B_cells  sample_1 control
table(seurat_obj$sample_id, seurat_obj$group)
#>            
#>             control treated
#>   sample_1      835       0
#>   sample_10       0    1004
#>   sample_11       0    1424
#>   sample_12       0     700
#>   sample_2      894       0
#>   sample_3      870       0
#>   sample_4     1237       0
#>   sample_5      570       0
#>   sample_6     1354       0
#>   sample_7        0    1470
#>   sample_8        0    1896
#>   sample_9        0    1025
```

## Convert to SummarizedExperiment

[`prepare_from_seurat()`](https://gladstone-institutes.github.io/LODopt/reference/prepare_from_seurat.md)
handles the conversion. You need to specify:

- **`cluster_col`** — the metadata column with cell type labels
- **`sample_col`** — the metadata column with sample identifiers
- **`pheno_cols`** — sample-level phenotype columns to include
- **`model_formula`** — the right-hand side of the model formula
- **`coef_of_interest`** — which coefficient to test (index or name)
- **`reference_levels`** — reference level for categorical variables

``` r
cellcomp_se <- prepare_from_seurat(
  seurat_obj,
  cluster_col = "cell_type",
  sample_col = "sample_id",
  pheno_cols = "group",
  model_formula = "group",
  coef_of_interest = 2,
  reference_levels = list(group = "control"),
  random_seed = 42,
  unchanged_clusters = NULL
)
cellcomp_se
#> class: SummarizedExperiment 
#> dim: 6 12 
#> metadata(5): modelFormula coef_of_interest_index
#>   reference_levels_of_variables random_seed unchanged_cluster_indices
#> assays(1): counts
#> rownames(6): B_cells Dendritic_cells ... T_cells_CD4 T_cells_CD8
#> rowData names(0):
#> colnames(12): sample_1 sample_10 ... sample_8 sample_9
#> colData names(1): group
```

## Run differential abundance analysis

Now pass the `SummarizedExperiment` to
[`logodds_optimized_normFactors()`](https://gladstone-institutes.github.io/LODopt/reference/logodds_optimized_normFactors.md):

``` r
results <- logodds_optimized_normFactors(cellcomp_se)
#> Running iteration no. 1
#> Running iteration no. 2
#> Found stable solution
```

## Examine results

The output is a list with `results` (per-cluster statistics) and
`optim_factor` (estimated normalization factors).

``` r
results$results
#>        cluster_id   comparison   estimates estimates_significance
#> 1         B_cells grouptreated -1.40405018           1.222965e-22
#> 2 Dendritic_cells grouptreated -0.23437437           9.607062e-02
#> 3       Monocytes grouptreated  0.87963908           2.887086e-04
#> 4        NK_cells grouptreated -0.04957772           7.023453e-01
#> 5     T_cells_CD4 grouptreated  0.17724844           1.535299e-01
#> 6     T_cells_CD8 grouptreated  0.11262103           3.598434e-01
#>   adjusted_pvalue
#> 1    7.337790e-22
#> 2    1.921412e-01
#> 3    8.661258e-04
#> 4    7.023453e-01
#> 5    2.302948e-01
#> 6    4.318121e-01
```

The key columns are:

- **`estimates`** — log-odds ratio for each cell type
- **`estimates_significance`** — raw p-value
- **`adjusted_pvalue`** — BH-adjusted p-value

``` r
results$optim_factor
#>     sampleid optim.norm.factor
#> 1   sample_1         0.9646303
#> 2  sample_10         0.9138434
#> 3  sample_11         1.0161021
#> 4  sample_12         1.1058483
#> 5   sample_2         1.0370502
#> 6   sample_3         1.0025196
#> 7   sample_4         0.9981553
#> 8   sample_5         0.9726147
#> 9   sample_6         1.0107198
#> 10  sample_7         0.9769640
#> 11  sample_8         1.0138792
#> 12  sample_9         0.9990960
```

## Conclusion

The Seurat workflow is straightforward:

1.  Start with a Seurat object containing cell type and sample
    annotations
2.  Use
    [`prepare_from_seurat()`](https://gladstone-institutes.github.io/LODopt/reference/prepare_from_seurat.md)
    to convert to a `SummarizedExperiment`
3.  Run
    [`logodds_optimized_normFactors()`](https://gladstone-institutes.github.io/LODopt/reference/logodds_optimized_normFactors.md)
    to test for differential abundance

For more details on the underlying method, see
[`vignette("intro", package = "LODopt")`](https://gladstone-institutes.github.io/LODopt/articles/intro.md).
