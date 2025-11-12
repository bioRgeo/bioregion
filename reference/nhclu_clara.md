# Non-hierarchical clustering: CLARA

This function performs non-hierarchical clustering based on
dissimilarity using partitioning around medoids, implemented via the
Clustering Large Applications (CLARA) algorithm.

## Usage

``` r
nhclu_clara(
  dissimilarity,
  index = names(dissimilarity)[3],
  seed = NULL,
  n_clust = c(1, 2, 3),
  maxiter = 0,
  initializer = "LAB",
  fasttol = 1,
  numsamples = 5,
  sampling = 0.25,
  independent = FALSE,
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

  A value for the random number generator (set to `NULL` for random
  initialization by default).

- n_clust:

  An `integer` vector or a single `integer` specifying the desired
  number(s) of clusters.

- maxiter:

  An `integer` defining the maximum number of iterations.

- initializer:

  A `character` string, either `"BUILD"` (used in the classic PAM
  algorithm) or `"LAB"` (Linear Approximate BUILD).

- fasttol:

  A positive `numeric` value defining the tolerance for fast swapping
  behavior. Defaults to 1.

- numsamples:

  A positive `integer` specifying the number of samples to draw.

- sampling:

  A positive `numeric` value defining the sampling rate.

- independent:

  A `boolean` indicating whether the previous medoids are excluded in
  the next sample. Defaults to `FALSE`.

- algorithm_in_output:

  A `boolean` indicating whether the original output of
  [fastclara](https://rdrr.io/pkg/fastkmedoids/man/fastclara.html)
  should be included in the output. Defaults to `TRUE` (see Value).

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
output of
[fastclara](https://rdrr.io/pkg/fastkmedoids/man/fastclara.html).

## Details

Based on [fastkmedoids](https://cran.r-project.org/package=fastkmedoids)
package
([fastclara](https://rdrr.io/pkg/fastkmedoids/man/fastclara.html)).

## References

Schubert E & Rousseeuw PJ (2019) Faster k-Medoids Clustering: Improving
the PAM, CLARA, and CLARANS Algorithms. *Similarity Search and
Applications* 11807, 171-187.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html>.

Associated functions:
[nhclu_clarans](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md)
[nhclu_dbscan](https://bioRgeo.github.io/bioregion/reference/nhclu_dbscan.md)
[nhclu_kmeans](https://bioRgeo.github.io/bioregion/reference/nhclu_kmeans.md)
[nhclu_pam](https://bioRgeo.github.io/bioregion/reference/nhclu_pam.md)
[nhclu_affprop](https://bioRgeo.github.io/bioregion/reference/nhclu_affprop.md)

## Author

Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "all")

#clust <- nhclu_clara(dissim, index = "Simpson", n_clust = 5)
   
```
