# Mock single-cell data for Seurat workflow examples

A small simulated dataset representing post-analysis scRNA-seq data,
suitable for demonstrating the
[`prepare_from_seurat`](https://gladstone-institutes.github.io/LODopt/reference/prepare_from_seurat.md)
workflow. Generated using
[`simulate_cellCounts_fromTissue`](https://gladstone-institutes.github.io/LODopt/reference/simulate_cellCounts_fromTissue.md)
with 6 cell types, 12 samples (6 control, 6 treated), and ~1500 cells
total. Two cell types (B_cells and Monocytes) have differential
abundance between groups.

## Usage

``` r
mock_seurat_data
```

## Format

A list with two elements:

- counts_matrix:

  A sparse gene-by-cell matrix (50 genes x ~1500 cells) of class
  `dgCMatrix`, suitable for
  [`SeuratObject::CreateSeuratObject()`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html).

- cell_metadata:

  A data frame with one row per cell and three columns:

  cell_type

  :   Character. Cell type label (one of B_cells, T_cells_CD4,
      T_cells_CD8, Monocytes, NK_cells, Dendritic_cells).

  sample_id

  :   Character. Sample identifier (sample_1 through sample_12).

  group

  :   Character. Experimental group ("control" or "treated").

  Row names are cell barcodes (cell_1, cell_2, ...).

## Source

Simulated using `data-raw/mock_seurat_data.R`.
