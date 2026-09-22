# Calculate metrics for a bioregionalization

This function calculates metrics at the bioregionalization level. These
evaluation metrics can be used to assess and select a partition based on
the optimal number of clusters.

## Usage

``` r
bioregionalization_metrics(
  bioregionalization,
  eval_metrics = c("prop_between_dissim", "anosim"),
  dissimilarity,
  dissimilarity_index = names(dissimilarity)[3],
  comat = NULL,
  anosim_permutations = 1,
  eval_metric = NULL,
  net = NULL,
  site_col = NULL,
  species_col = NULL
)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- eval_metrics:

  A `character` vector or a single `character` string indicating the
  metric(s) to be calculated. Available options are
  `"prop_between_dissim"`, `"anosim"`, `"mean_endemics"` or
  `"tot_endemics"`. Use `"all"` to compute all available metrics. See
  Details for metric descriptions.

- dissimilarity:

  A site-by-site dissimilarity object from
  [`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md).
  Required only for `"prop_between_dissim"` and `"anosim"`.

- dissimilarity_index:

  The name or number of the column to use as dissimilarity. By default,
  the third column name of `dissimilarity` is used.

- comat:

  A site-species `matrix` with sites as rows and species as columns.
  Should be provided if `eval_metrics` includes `"avg_endemism"` or
  `"tot_endemism"`.

- anosim_permutations:

  The number of permutations used to compute the p-value associated with
  the ANOSIM statistic. Defaults to 1, in which case no p-value is
  returned.

- eval_metric:

  Deprecated.

- net:

  Deprecated.

- site_col:

  Deprecated.

- species_col:

  Deprecated.

## Value

A `data.frame` containing the `eval_metrics` values for a
bioregionalization and its partition(s).

## Details

**Evaluation metrics:**

- [**prop_between_dissim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#prop_between_dissim):
  The proportion of total dissimilarity occurring between bioregions,
  following Holt et al. (2013), calculated as the sum of
  between-bioregion dissimilarities divided by the total sum of
  dissimilarities.

- [**anosim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#anosim):
  The Analysis of Similarities (ANOSIM) statistic, based on
  [vegan::anosim()](https://vegandevs.github.io/vegan/reference/anosim.html).
  It measures the separation between within- and between-bioregion
  dissimilarities. When `permutation > 1`, a p-value is calculated using
  permutation testing.

- [**mean_endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#mean_endemics):
  The mean proportion of endemic species across bioregions, following
  Kreft & Jetz (2010). For each bioregion, the proportion of endemic
  species is calculated and then averaged across bioregions.

- [**tot_endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#tot_endemics):
  The proportion of endemic species across all bioregion. It is
  calculated as the total number of endemic species divided by the total
  number of species. Endemic species are those occurring in only one
  bioregion.

## References

Holt BG, Lessard J, Borregaard MK, Fritz SA, Araújo MB, Dimitrov D,
Fabre P, Graham CH, Graves GR, Jønsson Ka, Nogués-Bravo D, Wang Z,
Whittaker RJ, Fjeldså J & Rahbek C (2013) An update of Wallace's
zoogeographic regions of the world. *Science* 339, 74-78.

Kreft H & Jetz W (2010) A framework for delineating biogeographical
regions based on species distributions. *Journal of Biogeography* 37,
2029-2053.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregionalization>.

Associated functions:
[compare_bioregionalizations](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md)
[find_optimal_n](https://bioRgeo.github.io/bioregion/reference/find_optimal_n.md)

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "all")

# User-defined number of clusters
bioreg <- hclu_hierarclust(dissim, 
                           n_clust = 10:15, 
                           index = "Simpson",
                           verbose = FALSE)
#> Warning: The requested number of cluster could not be found for k = 15. Closest number found: 14

met <- bioregionalization_metrics(bioreg,
                                  eval_metrics = "all",
                                  dissimilarity = dissim,
                                  comat = comat)
met
#>   partition n_bioregions prop_between_dissim    anosim mean_endemics
#> 1      K_10           10           0.9481112 0.6405984             0
#> 2      K_11           11           0.9573228 0.6034903             0
#> 3      K_12           12           0.9679755 0.6109700             0
#> 4      K_13           13           0.9768659 0.6856967             0
#> 5    K_14_1           14           0.9821236 0.6588603             0
#> 6    K_14_2           14           0.9821236 0.6588603             0
#>   tot_endemics
#> 1            0
#> 2            0
#> 3            0
#> 4            0
#> 5            0
#> 6            0
```
