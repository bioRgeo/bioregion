
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bioRgeo <img src="man/figures/logo.png" align="right" alt="" width="200" />

<!-- badges: start -->

[![R build
status](https://github.com/bioRgeo/bioRgeo/workflows/R-CMD-check/badge.svg)](https://github.com/bioRgeo/bioRgeo/actions)
[![version](https://img.shields.io/github/v/release/bioRgeo/bioRgeo?label=version&style=flat&logo=github)](https://github.com/bioRgeo/bioRgeo)
[![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
<!-- badges: end -->

This **R package** gathers a comprehensive set of algorithms to perform
bioregionalization annalyses. <br> The different algorithms come from
the network literature or can be (non)-hierarchical.

# :arrow_double_down: Installation

The package is not on CRAN yet and is still under active development.
You can install the development version from the GitHub repository with
the following command:

``` r
# install.packages("devtools")
devtools::install_github("bioRgeo/bioRgeo")
```

# :scroll: Vignettes

We wrote several vignettes that will help you using the **bioRgeo R
package**. Vignettes available are the following ones: <br>

- **[1. Installation of the executable binary
  files](https://biorgeo.github.io/bioRgeo/articles/a1_install_executable_binary_files.html)**  
- **[2. Matrix and network
  formats](https://biorgeo.github.io/bioRgeo/articles/a2_matrix_and_network_formats.html)**
- **[3. Pairwise similarity/dissimilarity
  metrics](https://biorgeo.github.io/bioRgeo/articles/a3_pairwise_metrics.html)**
- **[4.1 Hierarchical
  clustering](https://biorgeo.github.io/bioRgeo/articles/a4_1_hierarchical_clustering.html)**
- **[4.2 Non-hierarchical
  clustering](https://biorgeo.github.io/bioRgeo/articles/a4_2_non_hierarchical_clustering.html)**
- **[4.3 Network
  clustering](https://biorgeo.github.io/bioRgeo/articles/a4_3_network_clustering.html)**
- **[4.4
  Microbenchmark](https://biorgeo.github.io/bioRgeo/articles/a4_4_microbenchmark.html)**

Alternatively, if you prefer to view the vignettes in R, you can install
the package with `build_vignettes = TRUE`. But be aware that some
vignettes can be slow to generate.

``` r
remotes::install_github("bioRgeo/bioRgeo",
                        dependencies = TRUE, upgrade = "ask", 
                        build_vignettes = TRUE)

vignette("bioRgeo")
```

# :desktop_computer: Functions

An overview of all functions and data is given
**[here](https://biorgeo.github.io/bioRgeo/reference/index.html)**.

# :bug: Find a bug?

Thank you for finding it. Head over to the GitHub Issues tab and let us
know about it. Alternatively, you can also send us an e-mail. We will
try to get to it as soon as we can!

# References and dependencies

`bioRgeo` depends on `ape`, `bipartite`, `cluster`, `data.table`,
`dbscan`, `dendextend`, `dynamicTreeCut`, `earth`, `fastcluster`,
`ggplot2`, `grDevices`, `igraph`, `mathjaxr`, `Matrix`, `Rcpp`,
`Rdpack`, `rlang`, `rmarkdown`, `sf`, `stats`, `tidyr` and `utils`.
