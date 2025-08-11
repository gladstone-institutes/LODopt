![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg) [![GitHub issues](https://img.shields.io/github/issues/gladstone-institutes/LODopt)](https://github.com/gladstone-institutes/LODopt/issues) ![R >= 4.4.0](https://img.shields.io/badge/R-4.4.0-blue) ![RStudio](https://img.shields.io/badge/RStudio-2024.12.0%2B467-blue)


# LODopt
*Authors*: Reuben Thomas, Ayushi Agrawal, Natalie Elphick and Michela Traglia

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
## Tutorial
This [tutorial](https://gladstone-institutes.github.io/LODopt/articles/intro.html) introduces how to use the LODopt R package.

