# bioregion 1.2.0.9000

This is a list of changes made in the development/GitHub version of the package 
between bioregion 1.2.0 (CRAN release 2025-01-31) and the next CRAN release.

* Added `as_bioregion_pairwise()` to replace and improve upon
`betapart_to_bioregion()`, which is now deprecated.

* Added a comparison with other R packages for the computation of dissimilarity
  metrics in tutorial 3 (Pairwise similarity/dissimilarity metrics).

# bioregion 1.2.0

This is a list of changes made between **bioregion 1.1.1** 
(CRAN release 2024-04-19) and **bioregion 1.2.0** (CRAN release 2025-01-31).

* Added affinity propagation algorithm (`nhclu_affprop()`).

* Added a new method in `hclu_hierarclust()` to construct a consensus tree called
Iterative Hierarchical Consensus Tree (IHCT). This resolves issues related to 
the order of sites in the distance matrix and builds a consensus hierarchical 
tree with meaningful topology.

* Made many changes to functions related to `hclu_hierarclust()` due to 
this major update.

* Updated generic functions to provide `plot` and `print` methods for 
`hclu_diana()`.

* Added `site_species_metrics()` to the package and workflow.

* Added `bioregion_metrics()` to the package and workflow.

* Renamed `subset_node()` to `site_species_subset()`.

* Added indices `Cz` to `site_species_metrics()`.

* Updated `install_binaries()`:
  - Archive `bin.zip` now stored on GitHub and backed up on NextCloud.
  - Added Infomap version 2.8.0.
  - Added argument `download_only` to execute only the download step.  
&nbsp;

* Added `check_install` argument to `netclu_infomap()`, `netclu_louvain()`, 
and `netclu_oslom()`.

* Added `betapart_to_bioregion()` to the package.

* Added `compare_bioregionalizations()` to the package.

* Added `bioregionalization_metrics()` to the package.

* Updated documentation, vignettes, and tests.

* Modified the way seeds are generated for `nhclu_clara()` and 
`nhclu_clarans()`.
   
# bioregion 1.1.1

This is a list of changes made between **bioregion 1.1.0** 
(CRAN release 2024-03-19) and **bioregion 1.1.1** (CRAN release 2024-04-19).

* Added hierarchy support for Louvain (C++).

* Added `seed` argument to stochastic algorithms (except Louvain C++).

* Added `cut_weight` argument to `netclu_*` functions.

* Changed value for sites without clusters from `0` to `NA`.

* Updated automated tests (code coverage > 60%).

* Standardized controls, inputs, and outputs.

* Fixed a bug in `find_optimal_n()` for cases where partition metrics 
did not vary.

# bioregion 1.1.0

This is a list of changes made between **bioregion 1.0.0** 
(CRAN release 2023-04-15) and **bioregion 1.1.0** (CRAN release 2024-03-19).

* Added the `resolution` parameter to the igraph Louvain implementation.

* Added options to `mat_to_net()` to exclude diagonal and lower triangular 
matrices using `include_diag` and `include_lower`.

* Added a function to extract a subset of nodes (sites or species) from 
`bioregion.clusters` objects containing both types.

* Added a generic function to maintain attributes of `bioregion.pairwise.metric`
objects and track the number of sites and species.

* Added new functions: `nhclu_clara()` and `nhclu_clarans()`.

* Edited vignettes to document new functions.

* Modified controls for `bioregion.pairwise.metric` objects.

* Added the `include_formula` argument to 
`similarity_dissimilarity_conversion()` to (not) select formula metrics.

* Allowed negative values in `similarity()` with the Euclidean metric.

# bioregion 1.0.0 

First release on CRAN

