# Changelog

## bioregion 1.2.0.9000

This is a list of changes made in the development/GitHub version of the
package between bioregion 1.2.0 (CRAN release 2025-01-31) and the next
CRAN release.

***Function changes***

- Added the `inputs$data_type` field to all clustering outputs to
  explicitly track whether original co-occurrence data were
  occurrence-based or abundance-based. This field is automatically
  determined based on the algorithm type and the
  similarity/dissimilarity metric used.

- Added the `inputs$node_type` field to all clustering outputs to
  explicitly indicate whether the clustering includes only sites or both
  sites and species. This changes include the hidden `node_type`
  attributes.

- [`site_species_metrics()`](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md)
  has been thoroughly reformatted and now provides: species-to-bioregion
  indices, species-to-bioregionalization indices, site-to-species
  cluster indices, site-to-species clustering indices, site-to-bioregion
  indices, and site-to-bioregionalization indices, all of which are
  rigorously defined in the corresponding vignette.

- [`site_species_subset()`](https://bioRgeo.github.io/bioregion/reference/site_species_subset.md)
  has been simplified taking advantage on `inputs$node_type`.

- Renamed the class `bioregion.pairwise.metric` to `bioregion.pairwise`.

- Added a `verbose` argument to all talkative functions allowing users
  to control the display of progress messages.

***New features***

- Added export of the function
  [`exportGDF()`](https://bioRgeo.github.io/bioregion/reference/exportGDF.md)
  with documentation and tests.

- Added
  [`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md)
  to provide consistent bioregion color palettes for use across multiple
  visualizations (maps, networks, graphs, etc.).

- Updated
  [`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
  to handle bioregion colors.

- Added a generic function
  [`summary()`](https://rdrr.io/r/base/summary.html) for a clearer
  display of results.

- Added
  [`bind_pairwise()`](https://bioRgeo.github.io/bioregion/reference/bind_pairwise.md)
  to combine pairwise (dis)similarity objects.

- Added
  [`as_bioregion_pairwise()`](https://bioRgeo.github.io/bioregion/reference/as_bioregion_pairwise.md)
  to replace and improve upon
  [`betapart_to_bioregion()`](https://bioRgeo.github.io/bioregion/reference/betapart_to_bioregion.md),
  which is now deprecated.

- Added a comparison with other R packages for computing dissimilarity
  metrics in tutorial 3 (*Pairwise similarity/dissimilarity metrics*).

***Bug fixes***

- Fixed a problem in `inputs$pairwise_metric` when numeric `index` in
  all clustering outputs.

- Modified the `keep_trials` argument in
  [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)
  and fixed a potential issue with randomized matrix storage.

## bioregion 1.2.0

This is a list of changes made between **bioregion 1.1.1** (CRAN release
2024-04-19) and **bioregion 1.2.0** (CRAN release 2025-01-31).

- Added affinity propagation algorithm
  ([`nhclu_affprop()`](https://bioRgeo.github.io/bioregion/reference/nhclu_affprop.md)).

- Added a new method in
  [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)
  to construct a consensus tree called Iterative Hierarchical Consensus
  Tree (IHCT). This resolves issues related to the order of sites in the
  distance matrix and builds a consensus hierarchical tree with
  meaningful topology.

- Made many changes to functions related to
  [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)
  due to this major update.

- Updated generic functions to provide `plot` and `print` methods for
  [`hclu_diana()`](https://bioRgeo.github.io/bioregion/reference/hclu_diana.md).

- Added
  [`site_species_metrics()`](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md)
  to the package and workflow.

- Added
  [`bioregion_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregion_metrics.md)
  to the package and workflow.

- Renamed `subset_node()` to
  [`site_species_subset()`](https://bioRgeo.github.io/bioregion/reference/site_species_subset.md).

- Added indices `Cz` to
  [`site_species_metrics()`](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md).

- Updated
  [`install_binaries()`](https://bioRgeo.github.io/bioregion/reference/install_binaries.md):

  - Archive `bin.zip` now stored on GitHub and backed up on NextCloud.
  - Added Infomap version 2.8.0.
  - Added argument `download_only` to execute only the download step.  
     

- Added `check_install` argument to
  [`netclu_infomap()`](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md),
  [`netclu_louvain()`](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.md),
  and
  [`netclu_oslom()`](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.md).

- Added
  [`betapart_to_bioregion()`](https://bioRgeo.github.io/bioregion/reference/betapart_to_bioregion.md)
  to the package.

- Added
  [`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md)
  to the package.

- Added
  [`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)
  to the package.

- Updated documentation, vignettes, and tests.

- Modified the way seeds are generated for
  [`nhclu_clara()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
  and
  [`nhclu_clarans()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md).

## bioregion 1.1.1

This is a list of changes made between **bioregion 1.1.0** (CRAN release
2024-03-19) and **bioregion 1.1.1** (CRAN release 2024-04-19).

- Added hierarchy support for Louvain (C++).

- Added `seed` argument to stochastic algorithms (except Louvain C++).

- Added `cut_weight` argument to `netclu_*` functions.

- Changed value for sites without clusters from `0` to `NA`.

- Updated automated tests (code coverage \> 60%).

- Standardized controls, inputs, and outputs.

- Fixed a bug in
  [`find_optimal_n()`](https://bioRgeo.github.io/bioregion/reference/find_optimal_n.md)
  for cases where partition metrics did not vary.

## bioregion 1.1.0

This is a list of changes made between **bioregion 1.0.0** (CRAN release
2023-04-15) and **bioregion 1.1.0** (CRAN release 2024-03-19).

- Added the `resolution` parameter to the igraph Louvain implementation.

- Added options to
  [`mat_to_net()`](https://bioRgeo.github.io/bioregion/reference/mat_to_net.md)
  to exclude diagonal and lower triangular matrices using `include_diag`
  and `include_lower`.

- Added a function to extract a subset of nodes (sites or species) from
  `bioregion.clusters` objects containing both types.

- Added a generic function to maintain attributes of
  `bioregion.pairwise` objects and track the number of sites and
  species.

- Added new functions:
  [`nhclu_clara()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
  and
  [`nhclu_clarans()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md).

- Edited vignettes to document new functions.

- Modified controls for `bioregion.pairwise` objects.

- Added the `include_formula` argument to
  `similarity_dissimilarity_conversion()` to (not) select formula
  metrics.

- Allowed negative values in
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  with the Euclidean metric.

## bioregion 1.0.0

First release on CRAN
