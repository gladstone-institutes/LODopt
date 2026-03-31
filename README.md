[![DOI](https://zenodo.org/badge/1003934501.svg)](https://doi.org/10.5281/zenodo.19239198) ![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg) [![GitHub issues](https://img.shields.io/github/issues/gladstone-institutes/LODopt)](https://github.com/gladstone-institutes/LODopt/issues) ![R >= 4.4.0](https://img.shields.io/badge/R-4.4.0-blue) ![RStudio](https://img.shields.io/badge/RStudio-2024.12.0%2B467-blue)


# LODopt
*Authors*: Reuben Thomas, Ayushi Agrawal, Natalie Gill and Michela Traglia

## Background
LODopt is a statistical method to estimate changes in the absolute cell composition of different cell-types in a tissue using single cell assay-based counts of cells assigned to a given clusters/cell types per biological unit/sample. It is based on the optimal identification of sample-specific normalization factors to account for input tissue size differences followed by fitting a Generalized Linear Mixed Effects Model (GLMM) using the observed counts per sample per cluster. The range of scenarios over which the inferences derived from LODopt represent changes in (unobserved) absolute counts of different cell-types is much broader than those for other available approaches. Based on tests using simulated data, in comparison to exisiting approaches, LODOpt maintains the correct Type I error, has minimal bias in the estimated parameters of interest and comparable statistical power.


## Installation   
Currently the only way to install is by using the package`devtools`:    
```r
devtools::install_github("gladstone-institutes/LODopt")
```
If you get an error message and everything is spelled correctly, follow these steps before trying again:
```r
#set config
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")

#Go to github page to generate token
usethis::create_github_token() 

#paste your PAT into pop-up that follows...
credentials::set_github_pat()
```
## Tutorials
- [Simulated Data Workflow](https://gladstone-institutes.github.io/LODopt/articles/intro.html) — introduction to using LODopt with simulated single-cell data and a SummarizedExperiment object
- [Seurat Workflow](https://gladstone-institutes.github.io/LODopt/articles/seurat-workflow.html) — converting a Seurat object and running differential abundance analysis

## Citation
Thomas, R., Agrawal, A., Gill, N., & Traglia, M. (2026). gladstone-institutes/LODopt: v1.1.0. Zenodo. https://doi.org/10.5281/zenodo.19239199

## Suggested Methods Description
The LODopt [Thomas, R., Agrawal, A., Gill, N., & Traglia, M. (2026). gladstone-institutes/LODopt: v1.1.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19239199] method estimates the differences in log odds of cluster/cell-type membership between genotypes using a binomial generalized linear mixed effects model of the proportion of cluster-specific cells per mouse. The total number of cells per mouse are adjusted for between-sample differences by multiplicative normalization factors. The estimated proportions of cells per mouse or log odds of cluster membership are consequently also functions of these normalization factors. The optimal choice of these factors are determined by minimizing the total (over all clusters) variance of the log odds of cluster membership across mice. A variation of this method with unit normalization factors was previously used in Koutsodendris et al. [Koutsodendris, Nicole, Jessica Blumenfeld, Ayushi Agrawal, Michela Traglia, Brian Grone, Misha Zilberter, Oscar Yip et al. "Neuronal APOE4 removal protects against tau-mediated gliosis, neurodegeneration and myelin deficits." Nature aging 3, no. 3 (2023): 275-296.]  

## AI Disclosure Statement

Generative AI tools (Claude Code, Anthropic) were used as coding assistants during the development of this package. The authors maintain full responsibility for the accuracy, reproducibility, and scientific validity of all code. AI-assisted outputs were reviewed and validated against expected behavior before integration. The research questions, analytical approaches, parameter selections, and scientific interpretations were determined independently by the authors without AI input.
