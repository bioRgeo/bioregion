# Calculate metrics for sites and species relative to bioregions and chorotypes

This function computes metrics that quantify how species and sites
relate to clusters (bioregions or chorotypes). Depending on the type of
clustering, metrics can measure how species are distributed across
bioregions (site clusters), how sites relate to chorotypes (species
clusters), or both.

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
  metrics to compute for each cluster. Available metrics depend on the
  type of clustering (see arg `cluster_on`):

  - **When sites are clustered into bioregions** (default case):
    species-level metrics include `"Specificity"`, `"NSpecificity"`,
    `"Fidelity"`, `"IndVal"`, `"NIndVal"`, `"Rho"`, and `"CoreTerms"`.
    Site-level metrics include `"Richness"`, `"Rich_Endemics"`,
    `"Prop_Endemics"`, `"MeanSim"`, and `"SdSim"`.

  - **When species are clustered into chorotypes** (e.g., bipartite
    network clustering): site-level metrics include `"Specificity"`,
    `"NSpecificity"`, `"Fidelity"`, `"IndVal"`, `"NIndVal"`, `"Rho"`,
    and `"CoreTerms"`.

  Use `"all"` to compute all available metrics. See Details for metric
  descriptions.

- bioregionalization_metrics:

  A `character` vector or a single `character` string specifying summary
  metrics computed across all clusters. These metrics assess how an
  entity (species or site) is distributed across the entire
  bioregionalization, rather than relative to each individual cluster:

  - `"P"`: Participation coefficient measuring how evenly a species or
    site is distributed across clusters (0 = restricted to one cluster,
    1 = evenly spread).

  - `"Silhouette"`: How well a site fits its assigned bioregion compared
    to the nearest alternative bioregion (requires similarity data).

  Use `"all"` to compute all available metrics.

- data_type:

  A `character` string specifying whether metrics should be computed
  based on presence/absence (`"occurrence"`) or abundance values
  (`"abundance"`). This affects how Specificity, Fidelity, IndVal, Rho
  and CoreTerms are calculated:

  - `"auto"` (default): Automatically detected from input data
    (`bioregionalization` and/or `comat`).

  - `"occurrence"`: Metrics based on presence/absence only.

  - `"abundance"`: Metrics weighted by abundance values.

  - `"both"`: Compute both versions of the metrics.

- cluster_on:

  A `character` string specifying what was clustered in the
  bioregionalization, which determines what types of metrics can be
  computed:

  - `"site"` (default): Sites were clustered into bioregions. Metrics
    describe how each **species** is distributed across bioregions.

  - `"species"`: Species were clustered into chorotypes. Metrics
    describe how each **site** relates to chorotypes. Only available
    when species have been assigned to clusters (e.g., bipartite network
    clustering).

  - `"both"`: Compute metrics for both perspectives. Only available when
    both sites and species have cluster assignments.

- comat:

  A site-species `matrix` with sites as rows and species as columns.
  Values can be occurrence (1/0) or abundance. Required for most
  metrics.

- similarity:

  A site-by-site similarity object from
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md).
  Required only for similarity-based metrics (`"MeanSim"`, `"SdSim"`,
  `"Silhouette"`).

- include_cluster:

  A `boolean` indicating whether to add an `Assigned` column in the
  output, marking `TRUE` for rows where the site belongs to the
  bioregion being evaluated. Useful for quickly identifying a site's own
  bioregion. Default is `FALSE`.

- index:

  The name or number of the column to use as similarity. By default, the
  third column name of `similarity` is used.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` containing one or more `data.frame` elements, depending on the
selected metrics and clustering type:

**When sites are clustered (`cluster_on = "site"`):**

- **species_bioregions**: Metrics for each species x bioregion
  combination (e.g., Specificity, IndVal). One row per species x
  bioregion pair.

- **species_bioregionalization**: Summary metrics for each species
  across all bioregions (e.g., Participation coefficient). One row per
  species.

- **site_bioregions**: Metrics for each site x bioregion combination
  (e.g., MeanSim, Richness). One row per site x bioregion pair.

- **site_bioregionalization**: Summary metrics for each site (e.g.,
  Silhouette). One row per site.

**When species are clustered (`cluster_on = "species"`):**

- **site_chorotypes**: Metrics for each site x chorotype combination
  (e.g., Specificity, IndVal). One row per site x chorotype pair.

- **site_chorological**: Summary metrics for each site across all
  chorotypes (e.g., Participation coefficient). One row per site.

Note that if `bioregionalization` contains multiple partitions (i.e., if
`dim(bioregionalization$clusters) > 2`), a nested list will be returned,
with one sublist per partition.

## Details

This function computes metrics that characterize the relationship
between species, sites, and clusters. The available metrics depend on
whether you clustered sites (into bioregions) or species (into
chorotypes).

### — 1. Understanding the two perspectives —

- **Bioregions** are clusters of sites with similar species composition.

- **Chorotypes** are clusters of species with similar distributions.

In general, the package is designed to cluster sites into bioregions.
However, it is possible to group species into clusters. We call these
species clusters 'chorotypes', following conceptual definitions in the
biogeographical literature, to avoid any confusion in the calculation of
metrics.

In some cases, such as bipartite network clustering, both species and
sites receive the same clusters. We maintain the name distinction in the
calculation of metrics - but remember that in this case BIOREGION IDs =
CHOROTYPE IDs. The `cluster_on` argument determines which perspective to
use.

### — 2. Metrics when sites are clustered (`cluster_on = "site"` or `cluster_on = "both"`) —

**Species-per-bioregion metrics** quantify how each species is
distributed across bioregions.

These metrics are derived from three core terms ([see the online
vignette for a visual
diagram](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#metric-components)):

- **n_sb**: Number of sites in bioregion **b** where species **s** is
  present

- **n_s**: Total number of sites in which species **s** is present.

- **n_b**: Total number of sites in bioregion **b**.

Abundance version of these core terms can also be calculated when
`data_type = "abundance"` (or `data_type = "auto"` and
`bioregionalization was based on abundance`):

- **w_sb**: Sum of abundances of species **s** in sites of bioregion
  **b**.

- **w_s**: Total abundance of species **s**.

- **w_b**: Total abundance of all species present in sites of bioregion
  **b**.

The species-per-bioregion metrics are (click on metric names to access
formulas):

- [**Specificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#specificity-occurrence):
  Fraction of a species' occurrences found in a given bioregion (De
  Cáceres & Legendre 2009). A value of 1 means the species occurs only
  in that bioregion.

- [**NSpecificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nspecificity-occurrence):
  Normalized specificity that accounts for differences in bioregion size
  (De Cáceres & Legendre 2009).

- [**Fidelity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#fidelity-occurrence):
  Fraction of sites in a bioregion where the species occurs (De Cáceres
  & Legendre 2009). A value of 1 means the species is present in all
  sites of that bioregion.

- [**IndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#indval-occurrence):
  Indicator Value = Specificity × Fidelity (De Cáceres & Legendre 2009).
  High values identify species that are both restricted to and frequent
  within a bioregion.

- [**NIndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nindval-occurrence):
  Normalized IndVal accounting for bioregion size (De Cáceres & Legendre
  2009).

- [**Rho**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#rho-occurrence):
  Standardized contribution index comparing observed vs. expected
  co-occurrence under random association (Lenormand 2019).

- **CoreTerms**: Raw counts (n, n_b, n_s, n_sb) for custom calculations.

These metrics can be found in the output slot `species_bioregions`.

**Site-per-bioregion metrics** characterize sites relative to
bioregions:

- [**Richness**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics):
  Number of species in the site.

- [**Rich_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics):
  Number of species in the site that are endemic to one bioregion.

- [**Prop_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics):
  Proportion of endemic species in the site.

- [**MeanSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#meansim):
  Mean similarity of a site to all sites in each bioregion.

- [**SdSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#sdsim):
  Standard deviation of similarity values.

These metrics can be found in the output slot `site_bioregions`.

**Summary metrics across the whole bioregionalization:**

These metrics summarize how an entity (species or site) is distributed
across all clusters, rather than in relation to each individual cluster.

*Species-level summary metric:*

- [**P**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#p-occurrence-1)
  (Participation): Evenness of species distribution across bioregions
  (Denelle et al. 2020). Found in output slot
  `species_bioregionalization`.

*Site-level summary metric:*

- [**Silhouette**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#silhouette):
  How well a site fits its assigned bioregion vs. the nearest
  alternative (Rousseeuw 1987). Found in output slot
  `site_bioregionalization`.

### — 3. Metrics when species are clustered (`cluster_on = "species"` or `cluster_on = "both"`) —

**Site-per-chorotype metrics** quantify how each site relates to species
clusters (chorotypes).

The same metrics as above (Specificity, Fidelity, IndVal, etc.) can be
computed, but their interpretation is inverted. These metrics are based
on the following core terms:

- **n_gc**: Number of species belonging to chorotype **c** that are
  present in site **g**.

- **n_g**: Total number of species present in site **g**.

- **n_c**: Total number of species belonging to chorotype **c**.

Abundance version of these core terms can also be calculated when
`data_type = "abundance"` (or `data_type = "auto"` and
`bioregionalization was based on abundance`).

Their interpretation changes, for example:

- **Specificity**: Fraction of a site's species belonging to a
  chorotype.

- **Fidelity**: Fraction of a chorotype's species present in the site.

- **IndVal**: Indicator value for site-chorotype associations.

- **P**: Evenness of sites across chorotypes

## Note

If `data_type = "auto"`, the choice between occurrence- or abundance-
based metrics will be determined automatically from the input data, and
a message will explain the choice made.

Strict matching between entity IDs (site and species IDs) in
`bioregionalization` and in `comat` / `similarity` is required.

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
