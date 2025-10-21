# bioregion 1.2.0.9000

This is a list of changes made in the development/GitHub version of the package 
between bioregion 1.2.0 (CRAN release 2025-01-31) and the next CRAN release.


***Function changes***

* Added the `inputs$data_type` field to all clustering outputs to explicitly
  track whether original co-occurrence data were occurrence-based or
  abundance-based. This field is automatically determined based on the algorithm
  type and the similarity/dissimilarity metric used.

* `site_species_metrics()` is now vectorized for much faster performance. It
  also allows providing only `comat`, `net`, or both. When only one is provided,
  the other is automatically created using appropriate weight handling.

* Renamed the class `bioregion.pairwise.metric` to `bioregion.pairwise`.

* Added a `verbose` argument to all talkative functions allowing users to
  control the display of progress messages.

***New features***

* Added export of the function `exportGDF()` with documentation and tests.

* Added `bioregion_colors()` to provide consistent bioregion color palettes for
  use across multiple visualizations (maps, networks, graphs, etc.).

* Updated `map_bioregions()` to handle bioregion colors.

* Updated `site_species_metrics()` to handle multiple bioregionalizations
  simultaneously.

* Added a generic function `summary()` for a clearer display of results.

* Added `bind_pairwise()` to combine pairwise (dis)similarity objects.

* Added `as_bioregion_pairwise()` to replace and improve upon
  `betapart_to_bioregion()`, which is now deprecated.

* Added a comparison with other R packages for computing dissimilarity metrics
  in tutorial 3 (*Pairwise similarity/dissimilarity metrics*).

***Bug fixes***

* Fixed incorrect weight detection in `site_species_metrics()` for unipartite
  networks.

* Modified the `keep_trials` argument in `hclu_hierarclust()` and fixed a
  potential issue with randomized matrix storage.

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

* Added a generic function to maintain attributes of `bioregion.pairwise`
objects and track the number of sites and species.

* Added new functions: `nhclu_clara()` and `nhclu_clarans()`.

* Edited vignettes to document new functions.

* Modified controls for `bioregion.pairwise` objects.

* Added the `include_formula` argument to 
`similarity_dissimilarity_conversion()` to (not) select formula metrics.

* Allowed negative values in `similarity()` with the Euclidean metric.

# bioregion 1.0.0 

First release on CRAN

