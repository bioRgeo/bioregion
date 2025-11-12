# Non-hierarchical clustering: Partitioning Around Medoids

This function performs non-hierarchical clustering based on
dissimilarity using partitioning around medoids (PAM).

## Usage

``` r
nhclu_pam(
  dissimilarity,
  index = names(dissimilarity)[3],
  seed = NULL,
  n_clust = c(1, 2, 3),
  variant = "faster",
  nstart = 1,
  cluster_only = FALSE,
  algorithm_in_output = TRUE,
  ...
)
```

## Arguments

- dissimilarity:

  The output object from
  [`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
  or
  [`similarity_to_dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md),
  or a `dist` object. If a `data.frame` is used, the first two columns
  should represent pairs of sites (or any pair of nodes), and the
  subsequent column(s) should contain the dissimilarity indices.

- index:

  The name or number of the dissimilarity column to use. By default, the
  third column name of `dissimilarity` is used.

- seed:

  A value for the random number generator (`NULL` for random by
  default).

- n_clust:

  An `integer` vector or a single `integer` value specifying the
  requested number(s) of clusters.

- variant:

  A `character` string specifying the PAM variant to use. Defaults to
  `faster`. Available options are `original`, `o_1`, `o_2`, `f_3`,
  `f_4`, `f_5`, or `faster`. See
  [pam](https://rdrr.io/pkg/cluster/man/pam.html) for more details.

- nstart:

  An `integer` specifying the number of random starts for the PAM
  algorithm. Defaults to 1 (for the `faster` variant).

- cluster_only:

  A `boolean` specifying whether only the clustering results should be
  returned from the [pam](https://rdrr.io/pkg/cluster/man/pam.html)
  function. Setting this to `TRUE` makes the function more efficient.

- algorithm_in_output:

  A `boolean` indicating whether the original output of
  [pam](https://rdrr.io/pkg/cluster/man/pam.html) should be included in
  the result. Defaults to `TRUE` (see Value).

- ...:

  Additional arguments to pass to `pam()` (see
  [pam](https://rdrr.io/pkg/cluster/man/pam.html)).

## Value

A `list` of class `bioregion.clusters` with five components:

1.  **name**: A `character` string containing the name of the algorithm.

2.  **args**: A `list` of input arguments as provided by the user.

3.  **inputs**: A `list` of characteristics of the clustering process.

4.  **algorithm**: A `list` of all objects associated with the
    clustering procedure, such as original cluster objects (only if
    `algorithm_in_output = TRUE`).

5.  **clusters**: A `data.frame` containing the clustering results.

If `algorithm_in_output = TRUE`, the `algorithm` slot includes the
output of [pam](https://rdrr.io/pkg/cluster/man/pam.html).

## Details

This method partitions the data into the chosen number of clusters based
on the input dissimilarity matrix. It is more robust than k-means
because it minimizes the sum of dissimilarities between cluster centers
(medoids) and points assigned to the cluster. In contrast, k-means
minimizes the sum of squared Euclidean distances, which makes it
unsuitable for dissimilarity matrices that are not based on Euclidean
distances.

## References

Kaufman L & Rousseeuw PJ (2009) Finding groups in data: An introduction
to cluster analysis. In & Sons. JW (ed.), Finding groups in data: An
introduction to cluster analysis.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html>.

Associated functions:
[nhclu_clara](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
[nhclu_clarans](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md)
[nhclu_dbscan](https://bioRgeo.github.io/bioregion/reference/nhclu_dbscan.md)
[nhclu_kmeans](https://bioRgeo.github.io/bioregion/reference/nhclu_kmeans.md)
[nhclu_affprop](https://bioRgeo.github.io/bioregion/reference/nhclu_affprop.md)

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

comnet <- mat_to_net(comat)
dissim <- dissimilarity(comat, metric = "all")

clust <- nhclu_pam(dissim, n_clust = 2:15, index = "Simpson")
   
```
