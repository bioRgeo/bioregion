# Calculate metrics for one or several bioregionalizations

This function calculates metrics for one or several bioregionalizations,
typically based on outputs from `netclu_`, `hclu_`, or `nhclu_`
functions. Some metrics may require users to provide either a similarity
or dissimilarity matrix, or the initial species-site table.

## Usage

``` r
bioregionalization_metrics(
  bioregionalization,
  dissimilarity = NULL,
  dissimilarity_index = NULL,
  net = NULL,
  site_col = 1,
  species_col = 2,
  eval_metric = "all"
)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- dissimilarity:

  A `dist` object or a `bioregion.pairwise` object (output from
  [`similarity_to_dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md)).
  Required if `eval_metric` includes `"pc_distance"` and `tree` is not a
  `bioregion.hierar.tree` object.

- dissimilarity_index:

  A `character` string indicating the dissimilarity (beta-diversity)
  index to use if dissimilarity is a `data.frame` with multiple
  dissimilarity indices.

- net:

  The site-species network (i.e., bipartite network). Should be provided
  as a `data.frame` if `eval_metric` includes `"avg_endemism"` or
  `"tot_endemism"`.

- site_col:

  The name or index of the column representing site nodes (i.e., primary
  nodes). Should be provided if `eval_metric` includes `"avg_endemism"`
  or `"tot_endemism"`.

- species_col:

  The name or index of the column representing species nodes (i.e.,
  feature nodes). Should be provided if `eval_metric` includes
  `"avg_endemism"` or `"tot_endemism"`.

- eval_metric:

  A `character` vector or a single `character` string indicating the
  metric(s) to be calculated to assess the effect of different numbers
  of clusters. Available options are `"pc_distance"`, `"anosim"`,
  `"avg_endemism"`, or `"tot_endemism"`. If `"all"` is specified, all
  metrics will be calculated.

## Value

A `list` of class `bioregion.bioregionalization.metrics` with two to
three elements:

- `args`: Input arguments.

- `evaluation_df`: A `data.frame` containing the `eval_metric` values
  for all explored numbers of clusters.

- `endemism_results`: If endemism calculations are requested, a list
  with the endemism results for each bioregionalization.

## Details

**Evaluation metrics:**

- `pc_distance`: This metric, as used by Holt et al. (2013), is the
  ratio of the between-cluster sum of dissimilarities (beta-diversity)
  to the total sum of dissimilarities for the full dissimilarity matrix.
  It is calculated in two steps:

  - Compute the total sum of dissimilarities by summing all elements of
    the dissimilarity matrix.

  - Compute the between-cluster sum of dissimilarities by setting
    within-cluster dissimilarities to zero and summing the matrix. The
    `pc_distance` ratio is obtained by dividing the between-cluster sum
    of dissimilarities by the total sum of dissimilarities.

- `anosim`: This metric is the statistic used in the Analysis of
  Similarities, as described in Castro-Insua et al. (2018). It compares
  between-cluster and within-cluster dissimilarities. The statistic is
  computed as: R = (r_B - r_W) / (N (N-1) / 4), where r_B and r_W are
  the average ranks of between-cluster and within-cluster
  dissimilarities, respectively, and N is the total number of sites.
  Note: This function does not estimate significance; for significance
  testing, use
  [vegan::anosim()](https://vegandevs.github.io/vegan/reference/anosim.html).

- `avg_endemism`: This metric is the average percentage of endemism in
  clusters, as recommended by Kreft & Jetz (2010). It is calculated as:
  End_mean = sum_i (E_i / S_i) / K, where E_i is the number of endemic
  species in cluster i, S_i is the number of species in cluster i, and K
  is the total number of clusters.

- `tot_endemism`: This metric is the total endemism across all clusters,
  as recommended by Kreft & Jetz (2010). It is calculated as: End_tot =
  E / C, where E is the total number of endemic species (i.e., species
  found in only one cluster) and C is the number of non-endemic species.

## References

Castro-Insua A, Gómez-Rodríguez C & Baselga A (2018) Dissimilarity
measures affected by richness differences yield biased delimitations of
biogeographic realms. *Nature Communications* 9, 9-11.

Holt BG, Lessard J, Borregaard MK, Fritz SA, Araújo MB, Dimitrov D,
Fabre P, Graham CH, Graves GR, Jønsson Ka, Nogués-Bravo D, Wang Z,
Whittaker RJ, Fjeldså J & Rahbek C (2013) An update of Wallace's
zoogeographic regions of the world. *Science* 339, 74-78.

Kreft H & Jetz W (2010) A framework for delineating biogeographical
regions based on species distributions. *Journal of Biogeography* 37,
2029-2053.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html#optimaln>.

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

comnet <- mat_to_net(comat)

dissim <- dissimilarity(comat, metric = "all")

# User-defined number of clusters
tree1 <- hclu_hierarclust(dissim, 
                          n_clust = 10:15, 
                          index = "Simpson")
#> Building the iterative hierarchical consensus tree... Note that this process can take time especially if you have a lot of sites.
#> Error in if (nrow(dist_matrix) > 2) {    subclusters <- stable_binary_split(dist_matrix, method = method,         n_runs = n_runs, binsplit = "tree")} else {    subclusters <- list(attr(dist_matrix, "Labels")[1], attr(dist_matrix,         "Labels")[2])}: argument is of length zero
tree1
#> Error: object 'tree1' not found

a <- bioregionalization_metrics(tree1, 
                                dissimilarity = dissim, 
                                net = comnet,
                                site_col = "Node1", 
                                species_col = "Node2",
                                eval_metric = c("tot_endemism", 
                                                "avg_endemism",
                                                "pc_distance", 
                                                "anosim"))
#> Error: object 'tree1' not found
a
#> Error: object 'a' not found
```
