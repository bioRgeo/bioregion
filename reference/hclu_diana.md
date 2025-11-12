# Divisive hierarchical clustering based on dissimilarity or beta-diversity

This function computes a divisive hierarchical clustering from a
dissimilarity (beta-diversity) `data.frame`, calculates the cophenetic
correlation coefficient, and can generate clusters from the tree if
requested by the user. The function implements randomization of the
dissimilarity matrix to generate the tree, with a selection method based
on the optimal cophenetic correlation coefficient. Typically, the
dissimilarity `data.frame` is a `bioregion.pairwise` object obtained by
running `similarity` or `similarity` followed by
`similarity_to_dissimilarity`.

## Usage

``` r
hclu_diana(
  dissimilarity,
  index = names(dissimilarity)[3],
  n_clust = NULL,
  cut_height = NULL,
  find_h = TRUE,
  h_max = 1,
  h_min = 0,
  verbose = TRUE
)
```

## Arguments

- dissimilarity:

  The output object from
  [`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
  or
  [`similarity_to_dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md),
  or a `dist` object. If a `data.frame` is used, the first two columns
  represent pairs of sites (or any pair of nodes), and the remaining
  column(s) contain the dissimilarity indices.

- index:

  The name or number of the dissimilarity column to use. By default, the
  third column name of `dissimilarity` is used.

- n_clust:

  An `integer` vector or a single `integer` indicating the number of
  clusters to be obtained from the hierarchical tree, or the output from
  [bioregionalization_metrics](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md).
  Should not be used concurrently with `cut_height`.

- cut_height:

  A `numeric` vector indicating the height(s) at which the tree should
  be cut. Should not be used concurrently with `n_clust`.

- find_h:

  A `boolean` indicating whether the cutting height should be determined
  for the requested `n_clust`.

- h_max:

  A `numeric` value indicating the maximum possible tree height for the
  chosen `index`.

- h_min:

  A `numeric` value indicating the minimum possible height in the tree
  for the chosen `index`.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` of class `bioregion.clusters` with five slots:

1.  **name**: A `character` string containing the name of the algorithm.

2.  **args**: A `list` of input arguments as provided by the user.

3.  **inputs**: A `list` describing the characteristics of the
    clustering process.

4.  **algorithm**: A `list` containing all objects associated with the
    clustering procedure, such as the original cluster objects.

5.  **clusters**: A `data.frame` containing the clustering results.

## Details

The function is based on
[diana](https://rdrr.io/pkg/cluster/man/diana.html). Chapter 6 of
Kaufman & Rousseeuw (1990) fully details the functioning of the diana
algorithm.

To find an optimal number of clusters, see
[`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)

## References

Kaufman L & Rousseeuw PJ (2009) Finding groups in data: An introduction
to cluster analysis. In & Sons. JW (ed.), *Finding groups in data: An
introduction to cluster analysis*.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html>.

Associated functions:
[cut_tree](https://bioRgeo.github.io/bioregion/reference/cut_tree.md)

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

data("fishmat")
fishdissim <- dissimilarity(fishmat)
fish_diana <- hclu_diana(fishdissim, index = "Simpson")
#> Output tree has a 0.51 cophenetic correlation coefficient with the initial dissimilarity matrix

```
