# Calculate contribution metrics of sites and species to each bioregion

This function computes metrics to assess the contribution of a given
species or site to each bioregion.

## Usage

``` r
site_species_metrics(
  bioregionalization,
  bioregion_metrics = c("Specificity", "NSpecificity", "Fidelity", "IndVal", "NIndVal",
    "Rho"),
  bioregionalization_metrics = "P",
  data_type = "auto",
  cluster_on = "site",
  comat,
  similarity = NULL,
  include_cluster = FALSE,
  index = names(similarity)[3],
  verbose = TRUE
)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- bioregion_metrics:

  A `character` vector or a single `character` string specifying the
  metrics to compute. Several metrics belonging to different categories
  are available: species to bioregions (`"Specificity"`,
  `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, `"IndVal"`,
  `"NIndVal"`, `"Rho"` and `"CoreTerms"`), site to bioregions
  (`"Richness"`, `"Rich_Endemics"`, `"Prop_Endemics"`, `"MeanSim"` and
  `"SdSim"`), or site to chorotypes (`"Specificity"`, `"NSpecificity"`,
  `"Fidelity"`, `"NFidelity"`, `"IndVal"`, `"NIndVal"`, `"Rho"` and
  `"CoreTerms"`) when a cluster has been assigned to the species in
  `bioregionalization` (see Details). If `"all"` is specified, all
  metrics will be calculated.

- bioregionalization_metrics:

  A `character` vector or a single `character` string specifying the
  bioregionalization metrics to compute. Some aggregated metrics (such
  as the participation coefficient `"P"` or Silhouette index
  `"Silhouette"`) can be derived from the `bioregion_metrics` (see
  Details). If `"all"` is specified , all bioregionalization metrics
  will be calculated.

- data_type:

  A `character` string specifying which type of data should be
  considered to compute the related metrics (`"Specificity"`,
  ``` "NSpecificity"``"Fidelity" ```, `"NFidelity"`, `"IndVal"`,
  `"NIndVal"`, `"Rho"`, `"CoreTerms"`): occurrences or abundances. By
  default (`"auto"`), the type of data is inferred from
  `bioregionalization` and/or the provided co-occurrence matrix
  (argument `comat`). Other available options are `"occurence"`,
  `"abundance"` or `"both"` (see Details).

- cluster_on:

  A `character` string specifying whether the related metrics
  (`"Specificity"`, `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`,
  `"IndVal"`, `"NIndVal"`, `"Rho"` and `"CoreTerms"`) should be computed
  as contributions from species to bioregions (clusters based on
  `"site"` by default), from site to chorotypes (clusters based on
  `"species"`) or for `"both"`. The latter type of contribution is only
  available if a cluster has been assigned to the species in
  `bioregionalization` (see Details).

- comat:

  A co-occurrence `matrix` with sites as rows and species as columns.

- similarity:

  The output object from
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md).

- include_cluster:

  A `boolean` indicating whether to include an additional column
  `Assigned` in the site-to-bioregions output, indicating for each site
  whether the considered bioregion is its bioregion (`FALSE` by
  default).

- index:

  The name or number of the column to use as similarity. By default, the
  third column name of `similarity` is used.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` containing between one and six elements (listed below),
depending on the selected metrics (`bioregion_metrics` and/or
`bioregionalization_metrics`) and the type of cluster (bioregion when
`cluster_on = "site"` and/or chorotypes when `cluster_on = "species"`).

- **species_bioregions**: A `data.frame` containing the
  species-to-bioregions metric(s) based on `comat`.

- **species_bioregionalization**: A `data.frame` containing the
  species-to-bioregionalization metric(s) based on `comat`.

- **site_bioregions**: A `data.frame` containing the site-to-bioregions
  metric(s) based on `comat` or `similarity`.

- **site_bioregionalization**: A `data.frame` containing the
  site-to-bioregionalization metric(s) based on `comat` or `similarity`.

- **site_chorotypes**: A `data.frame` containing the site-to-chorotypes
  metric(s) based on `comat`.

- **site_chorological**: A `data.frame` containing the
  site-to-chorological classification metric(s) based on `comat`.

Note that if `bioregionalization` contains more than one partition
(i.e., if `dim(bioregionalization$clusters) > 2`), a list of lists will
be returned, with one sublist per partition.

## Details

The **first type** of metrics provided by this function is based on the
contribution of each species to a given bioregion (`cluster_on = "site"`
or `cluster_on = "both"`). This is calculated from the bioregion
assigned to each site (as defined in `bioregionalization`) and an
associated site–species co-occurrence matrix `comat` (a strict match
between site IDs is required).

Occurrence-based metrics (`data_type = "occurrence"` or
`data_type = "both"`) are derived from three core terms:

- **n_sb**: Number of sites belonging to bioregion **b** in which
  species **s** is present.

- **n_s**: Total number of sites in which species **s** is present.

- **n_b**: Total number of sites belonging to bioregion **b**.

These **species-to-bioregions** metrics include:

- [**Specificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#specificity-occurrence),
  as described in De Cáceres M & Legendre P (2009).

- [**NSpecificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nspecificity-occurrence)
  is the normalized version of **Specificity** accounting for bioregion
  size, as described in De Cáceres M & Legendre P (2009).

- [**Fidelity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#fidelity-occurrence)
  as described in De Cáceres M & Legendre P (2009).

- [**IndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#indval-occurrence),
  as described in De Cáceres M & Legendre P (2009).

- [**NIndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nindval-occurrence)
  is the normalized version of **IndVal** accounting for bioregion size,
  as described in De Cáceres M & Legendre P (2009).

- [**Rho**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#rho-occurrence),
  as described in Lenormand *et al.* (2019).

**Species-to-bioregionalization** metrics can also be computed, such as:

- [**Participation**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#participation-occurrence),
  as described in Denelle *et al.* (2020).

Abundance-weighted versions of these metrics (`data_type = "abundance"`
or `data_type = "both"`) can also be derived using the following
analogous core terms:

- **w_sb**: Sum of abundances of species **s** in sites of bioregion
  **b**.

- **w_s**: Total abundance of species **s**.

- **w_b**: Total abundance of all species present in sites of bioregion
  **b**.

These abundance-weighted terms correspond directly to their
occurrence-based counterparts and allow computing abundance versions of
the same metrics. Detailed formulas and examples are provided in the
package vignette.

The **second type** of metrics provided by this function is based on the
contribution of each site to a given bioregion (`cluster_on = "site"` or
`cluster_on = "both"`). This is calculated from the bioregion assigned
to each site (as defined in `bioregionalization`), an associated
site–species co-occurrence matrix `comat` (for `"Richness"`,
`"Rich_Endemics"` and `"Prop_Endemics"`, a strict match between site IDs
is required) and/or a site–site similarity metric `similarity` (for
`"MeanSim"` and `"SdSim"`, a strict match between site IDs is required).

These **site-to-bioregions** metrics include:

- [**Richness**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#richness):
  The number of species present in a site.

- [**Rich_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#richendemics):
  The number of species present in a site that are endemic to a
  bioregion.

- [**Prop_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#propendemics):
  The proportion of species present in a site that are endemic to a
  bioregion.

- [**MeanSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#meansim):
  The mean similarity of each site to the sites of every bioregion.

- [**SdSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#sdsim):
  The corresponding standard deviation of similarity values.

**Site-to-bioregionalization** metrics can also be computed, such as:

- [**Silhouette**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#silhouette),
  as described in Rousseeuw (1987) (similarity-based version).

The **third type** of metrics provided by this function is based on the
contribution of each site to a given *chorotypes*
(`cluster_on = "species"` or `cluster = "both"`, only when a cluster has
been assigned to species in `bioregionalization`). This is calculated
from the chorotype assigned to each species (as defined in
`bioregionalization`) and an associated site–species co-occurrence
matrix (a strict match between species IDs is required).

In this case, occurrence-based metrics (`data_type = "occurrence"` or
`data_type = "both"`) are derived from three core terms:

- **n_gc**: Number of species belonging to chorotype **c** that are
  present in site **g**.

- **n_g**: Total number of species present in site **g**.

- **n_c**: Total number of species belonging to chorotype **c**.

As for the first type, all metrics (including their abundance-weighted
versions when `data_type = "abundance"` or `data_type = "both"`) can be
derived from these **site-to-chorotypes** relationships. Further
details, mathematical definitions, and examples are provided in the
package vignette.

## Note

If `data_type = "auto"`, the choice between occurrence- or abundance-
based metrics will be determined automatically from the input data, and
a message will explain the choice made.

## References

De Cáceres M & Legendre P (2009) Associations between species and groups
of sites: indices and statistical inference. *Ecology* 90, 3566–3574.

Denelle P, Violle C & Munoz F (2020) Generalist plants are more
competitive and more functionally similar to each other than specialist
plants: insights from network analyses. *Journal of Biogeography* 47,
1922–-1933.

Lenormand M, Papuga G, Argagnon O, Soubeyrand M, Alleaume S & Luque S
(2019) Biogeographical network analysis of plant species distribution in
the Mediterranean region. *Ecology and Evolution* 9, 237–250.

Rousseeuw PJ (1987) Silhouettes: A graphical aid to the interpretation
and validation of cluster analysis. *Journal of Computational and
Applied Mathematics* 20, 53–65.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html>.

Associated functions:
[bioregion_metrics](https://bioRgeo.github.io/bioregion/reference/bioregion_metrics.md)
[bioregionalization_metrics](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
data(fishmat)

fishsim <- similarity(fishmat, metric = "Jaccard")

bioregionalization <- hclu_hierarclust(similarity_to_dissimilarity(fishsim),
                                       index = "Jaccard",
                                       method = "average",
                                       randomize = TRUE,
                                       optimal_tree_method = "best",
                                       n_clust = c(1,2,3),
                                       verbose = FALSE)
                                     
ind <- site_species_metrics(bioregionalization = bioregionalization,
                             bioregion_metrics = "all",
                             bioregionalization_metrics = "all",
                             data_type = "auto",
                             cluster_on = "site",
                             comat = fishmat,
                             similarity = fishsim,
                             include_cluster = TRUE,
                             index = 3,
                             verbose = TRUE)
#> The bioregionalization is based on occurence data and comat is based on occurence data so occurrence-based metrics will be computed.
                                     
```
