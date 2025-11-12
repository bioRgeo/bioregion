# Non-hierarchical clustering: Affinity Propagation

This function performs non-hierarchical clustering using the Affinity
Propagation algorithm.

## Usage

``` r
nhclu_affprop(
  similarity,
  index = names(similarity)[3],
  seed = NULL,
  p = NA,
  q = NA,
  maxits = 1000,
  convits = 100,
  lam = 0.9,
  details = FALSE,
  nonoise = FALSE,
  K = NULL,
  prc = NULL,
  bimaxit = NULL,
  exact = NULL,
  algorithm_in_output = TRUE,
  verbose = TRUE
)
```

## Arguments

- similarity:

  The output object from
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md),
  or a `dist` object. If a `data.frame` is used, the first two columns
  should represent pairs of sites (or any pair of nodes), and the
  subsequent column(s) should contain the similarity indices.

- index:

  The name or number of the similarity column to use. By default, the
  third column name of `similarity` is used.

- seed:

  The seed for the random number generator used when `nonoise = FALSE`.

- p:

  Input preference, which can be a vector specifying individual
  preferences for each data point. If scalar, the same value is used for
  all data points. If `NA`, exemplar preferences are initialized based
  on the distribution of non-Inf values in the similarity matrix,
  controlled by `q`.

- q:

  If `p = NA`, exemplar preferences are initialized according to the
  distribution of non-Inf values in the similarity matrix. By default,
  the median is used. A value between 0 and 1 specifies the sample
  quantile, where `q = 0.5` results in the median.

- maxits:

  The maximum number of iterations to execute.

- convits:

  The algorithm terminates if the exemplars do not change for `convits`
  iterations.

- lam:

  The damping factor, a value in the range \[0.5, 1). Higher values
  correspond to heavier damping, which may help prevent oscillations.

- details:

  If `TRUE`, detailed information about the algorithm's progress is
  stored in the output object.

- nonoise:

  If `TRUE`, disables the addition of a small amount of noise to the
  similarity object, which prevents degenerate cases.

- K:

  The desired number of clusters. If not `NULL`, the function
  [apclusterK](https://rdrr.io/pkg/apcluster/man/apclusterK-methods.html)
  is called.

- prc:

  A parameter needed when `K` is not `NULL`. The algorithm stops if the
  number of clusters deviates by less than `prc` percent from the
  desired value `K`. Set to 0 to enforce exactly `K` clusters.

- bimaxit:

  A parameter needed when `K` is not `NULL`. Specifies the maximum
  number of bisection steps to perform. No warning is issued if the
  number of clusters remains outside the desired range.

- exact:

  A flag indicating whether to compute the initial preference range
  exactly.

- algorithm_in_output:

  A `boolean` indicating whether to include the original output of
  [apcluster](https://rdrr.io/pkg/apcluster/man/apcluster-methods.html)
  in the result. Defaults to `TRUE`.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` of class `bioregion.clusters` with five slots:

1.  **name**: A `character` string containing the name of the algorithm.

2.  **args**: A `list` of input arguments as provided by the user.

3.  **inputs**: A `list` describing the characteristics of the
    clustering process.

4.  **algorithm**: A `list` of objects associated with the clustering
    procedure, such as original cluster objects (if
    `algorithm_in_output = TRUE`).

5.  **clusters**: A `data.frame` containing the clustering results.

If `algorithm_in_output = TRUE`, the `algorithm` slot includes the
output of
[apcluster](https://rdrr.io/pkg/apcluster/man/apcluster-methods.html).

## Details

This function is based on the
[apcluster](https://cran.r-project.org/package=apcluster) package
([apcluster](https://rdrr.io/pkg/apcluster/man/apcluster-methods.html)).

## References

Frey B & Dueck D (2007) Clustering by Passing Messages Between Data
Points. *Science* 315, 972-976.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html>.

Associated functions:
[nhclu_clara](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
[nhclu_clarans](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md)
[nhclu_dbscan](https://bioRgeo.github.io/bioregion/reference/nhclu_dbscan.md)
[nhclu_kmeans](https://bioRgeo.github.io/bioregion/reference/nhclu_kmeans.md)
nhclu_affprop

## Author

Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)

## Examples

``` r
comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
prob = 1/1:1001), 10, 12)
rownames(comat_1) <- paste0("Site", 1:10)
colnames(comat_1) <- paste0("Species", 1:12)
comat_1 <- cbind(comat_1,
                 matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 1:10),
                                        paste0("Species", 13:20))))
                                        
comat_2 <- matrix(sample(0:1000, 
                         size = 10*12, 
                         replace = TRUE, 
                         prob = 1/1:1001), 
                  10, 12)
rownames(comat_2) <- paste0("Site", 11:20)
colnames(comat_2) <- paste0("Species", 9:20)
comat_2 <- cbind(matrix(0, 10, 8, 
                        dimnames = list(paste0("Site", 11:20),
                                        paste0("Species", 1:8))),
                 comat_2)
                 
comat <- rbind(comat_1, comat_2)

dissim <- dissimilarity(comat, metric = "Simpson")
sim <- dissimilarity_to_similarity(dissim)

clust1 <- nhclu_affprop(sim)

clust2 <- nhclu_affprop(sim, q = 1)

# Fixed number of clusters 
clust3 <- nhclu_affprop(sim, K = 2, prc = 10, bimaxit = 20, exact = FALSE)
#> Trying p = 0.9930872 
#>    Number of clusters: 6 
#> Trying p = 0.9308716 
#>    Number of clusters: 6 
#> Trying p = 0.3087157 
#>    Number of clusters: 1 
#> Trying p = 0.6543579 (bisection step no. 1 )
#>    Number of clusters: 1 
#> Trying p = 0.8271789 (bisection step no. 2 )
#>    Number of clusters: 3 
#> Trying p = 0.7407684 (bisection step no. 3 )
#>    Number of clusters: 1 
#> Trying p = 0.7839737 (bisection step no. 4 )
#>    Number of clusters: 1 
#> Trying p = 0.8055763 (bisection step no. 5 )
#>    Number of clusters: 3 
#> Trying p = 0.794775 (bisection step no. 6 )
#>    Number of clusters: 2 
#> 
#> Number of clusters: 2 for p = 0.794775 
```
