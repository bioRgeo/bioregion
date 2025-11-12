# Cut a hierarchical tree

This function is designed to work on a hierarchical tree and cut it at
user-selected heights. It works with outputs from either
`hclu_hierarclust` or `hclust` objects. The function allows for cutting
the tree based on the chosen number(s) of clusters or specified
height(s). Additionally, it includes a procedure to automatically
determine the cutting height for the requested number(s) of clusters.

## Usage

``` r
cut_tree(
  tree,
  n_clust = NULL,
  cut_height = NULL,
  find_h = TRUE,
  h_max = 1,
  h_min = 0,
  dynamic_tree_cut = FALSE,
  dynamic_method = "tree",
  dynamic_minClusterSize = 5,
  dissimilarity = NULL,
  show_hierarchy = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- tree:

  A `bioregion.hierar.tree` or an `hclust` object.

- n_clust:

  An `integer` vector or a single `integer` indicating the number of
  clusters to be obtained from the hierarchical tree, or the output from
  [`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md).
  This should not be used concurrently with `cut_height`.

- cut_height:

  A `numeric` vector specifying the height(s) at which the tree should
  be cut. This should not be used concurrently with `n_clust` or
  `optim_method`.

- find_h:

  A `boolean` indicating whether the cutting height should be determined
  for the requested `n_clust`.

- h_max:

  A `numeric` value indicating the maximum possible tree height for
  determining the cutting height when `find_h = TRUE`.

- h_min:

  A `numeric` value specifying the minimum possible height in the tree
  for determining the cutting height when `find_h = TRUE`.

- dynamic_tree_cut:

  A `boolean` indicating whether the dynamic tree cut method should be
  used. If `TRUE`, `n_clust` and `cut_height` are ignored.

- dynamic_method:

  A `character` string specifying the method to be used for dynamically
  cutting the tree: either `"tree"` (clusters searched only within the
  tree) or `"hybrid"` (clusters searched in both the tree and the
  dissimilarity matrix).

- dynamic_minClusterSize:

  An `integer` indicating the minimum cluster size for the dynamic tree
  cut method (see
  [dynamicTreeCut::cutreeDynamic()](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)).

- dissimilarity:

  Relevant only if `dynamic_method = "hybrid"`. Provide the
  dissimilarity `data.frame` used to build the `tree`.

- show_hierarchy:

  A `boolean` specifying if the hierarchy of clusters should be
  identifiable in the outputs (`FALSE` by default).

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

- ...:

  Additional arguments passed to
  [dynamicTreeCut::cutreeDynamic()](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)
  to customize the dynamic tree cut method.

## Value

If `tree` is an output from
[`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md),
the same object is returned with updated content (i.e., `args` and
`clusters`). If `tree` is an `hclust` object, a `data.frame` containing
the clusters is returned.

## Details

The function supports two main methods for cutting the tree. First, the
tree can be cut at a uniform height (specified by `cut_height` or
determined automatically for the requested `n_clust`). Second, the
dynamic tree cut method (Langfelder et al., 2008) can be applied, which
adapts to the shape of branches in the tree, cutting at varying heights
based on cluster positions.

The dynamic tree cut method has two variants:

- The tree-based variant (`dynamic_method = "tree"`) uses a top-down
  approach, relying solely on the tree and the order of clustered
  objects.

- The hybrid variant (`dynamic_method = "hybrid"`) employs a bottom-up
  approach, leveraging both the tree and the dissimilarity matrix to
  identify clusters based on dissimilarity among sites. This approach is
  useful for detecting outliers within clusters.

## Note

The `find_h` argument is ignored if `dynamic_tree_cut = TRUE`, as
cutting heights cannot be determined in this case.

## References

Langfelder P, Zhang B & Horvath S (2008) Defining clusters from a
hierarchical cluster tree: the Dynamic Tree Cut package for R.
*BIOINFORMATICS* 24, 719-720.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html>.

Associated functions:
[hclu_hierarclust](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)

## Author

Pierre Denelle (<pierre.denelle@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

simil <- similarity(comat, metric = "all")
dissimilarity <- similarity_to_dissimilarity(simil)

# User-defined number of clusters
tree1 <- hclu_hierarclust(dissimilarity,
                          n_clust = 5)
#> Building the iterative hierarchical consensus tree... Note that this process can take time especially if you have a lot of sites.
#> Error in if (nrow(dist_matrix) > 2) {    subclusters <- stable_binary_split(dist_matrix, method = method,         n_runs = n_runs, binsplit = "tree")} else {    subclusters <- list(attr(dist_matrix, "Labels")[1], attr(dist_matrix,         "Labels")[2])}: argument is of length zero
tree2 <- cut_tree(tree1, cut_height = .05)
#> Error: object 'tree1' not found
tree3 <- cut_tree(tree1, n_clust = c(3, 5, 10))
#> Error: object 'tree1' not found
tree4 <- cut_tree(tree1, cut_height = c(.05, .1, .15, .2, .25))
#> Error: object 'tree1' not found
tree5 <- cut_tree(tree1, n_clust = c(3, 5, 10), find_h = FALSE)
#> Error: object 'tree1' not found

hclust_tree <- tree2$algorithm$final.tree
#> Error: object 'tree2' not found
clusters_2 <- cut_tree(hclust_tree, n_clust = 10)
#> Error: object 'hclust_tree' not found

cluster_dynamic <- cut_tree(tree1, dynamic_tree_cut = TRUE,
                            dissimilarity = dissimilarity)
#> Error: object 'tree1' not found
```
