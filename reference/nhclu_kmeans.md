# Non-hierarchical clustering: K-means analysis

This function performs non-hierarchical clustering based on
dissimilarity using a k-means analysis.

## Usage

``` r
nhclu_kmeans(
  dissimilarity,
  index = names(dissimilarity)[3],
  seed = NULL,
  n_clust = c(1, 2, 3),
  iter_max = 10,
  nstart = 10,
  algorithm = "Hartigan-Wong",
  algorithm_in_output = TRUE
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

- iter_max:

  An `integer` specifying the maximum number of iterations for the
  k-means method (see [kmeans](https://rdrr.io/r/stats/kmeans.html)).

- nstart:

  An `integer` specifying how many random sets of `n_clust` should be
  selected as starting points for the k-means analysis (see
  [kmeans](https://rdrr.io/r/stats/kmeans.html)).

- algorithm:

  A `character` specifying the algorithm to use for k-means (see
  [kmeans](https://rdrr.io/r/stats/kmeans.html)). Available options are
  Hartigan-Wong, Lloyd, Forgy, and MacQueen.

- algorithm_in_output:

  A `boolean` indicating whether the original output of
  [kmeans](https://rdrr.io/r/stats/kmeans.html) should be included in
  the output. Defaults to `TRUE` (see Value).

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
output of [kmeans](https://rdrr.io/r/stats/kmeans.html).

## Details

This method partitions data into k groups such that the sum of squares
of Euclidean distances from points to the assigned cluster centers is
minimized. K-means cannot be applied directly to dissimilarity or
beta-diversity metrics because these distances are not Euclidean.
Therefore, it first requires transforming the dissimilarity matrix using
Principal Coordinate Analysis (PCoA) with
[pcoa](https://rdrr.io/pkg/ape/man/pcoa.html), and then applying k-means
to the coordinates of points in the PCoA.

Because this additional transformation alters the initial dissimilarity
matrix, the partitioning around medoids method
([nhclu_pam](https://bioRgeo.github.io/bioregion/reference/nhclu_pam.md))
is preferred.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html>.

Associated functions:
[nhclu_clara](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
[nhclu_clarans](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md)
[nhclu_dbscan](https://bioRgeo.github.io/bioregion/reference/nhclu_dbscan.md)
[nhclu_pam](https://bioRgeo.github.io/bioregion/reference/nhclu_pam.md)
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

clust <- nhclu_kmeans(dissim, n_clust = 2:10, index = "Simpson")
```
