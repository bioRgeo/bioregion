README
================

# bioRgeo <img src="figures/France.PNG" align="right" alt="" width="200" />

# An R package to compare several methods of bioregionalization and estimate the contributions of species

[![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
[![Travis build status](https://travis-ci.org/bioRgeo/bioRgeo)](https://travis-ci.org/bioRgeo/bioRgeo)

`bioRgeo` is an R package designed to compare several methods of
clustering based on multivariate or network approaches.

## 1 Installation

``` r
devtools::install_github("bioRgeo/bioRgeo", build_vignettes = TRUE)
```

## 2 Tutorial

A tutorial vignette showing the main steps of the processing can be
visualised with the following command line:

``` r
vignette("tutorial", package = "bioRgeo")
```

    ## Warning: vignette 'tutorial' not found

html of the vignette can also be accessed **[here](https://biorgeo.github.io/bioRgeo/articles/tutorial.html)**

## 3 Dependencies

`bioRgeo` depends on `dplyr`, `ecodist`, `reshape2`, `fastcluster`,
`DescTools`, `ade4`, `cluster`, `sf`, `ggplot2`, `Rcpp`, `SMUT` and
`igraph`.
