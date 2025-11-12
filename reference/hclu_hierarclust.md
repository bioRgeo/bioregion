# Hierarchical clustering based on dissimilarity or beta-diversity

This function generates a hierarchical tree from a dissimilarity
(beta-diversity) `data.frame`, calculates the cophenetic correlation
coefficient, and optionally retrieves clusters from the tree upon user
request. The function includes a randomization process for the
dissimilarity matrix to generate the tree, with two methods available
for constructing the final tree. Typically, the dissimilarity
`data.frame` is a `bioregion.pairwise` object obtained by running
`similarity`, or by running `similarity` followed by
`similarity_to_dissimilarity`.

## Usage

``` r
hclu_hierarclust(
  dissimilarity,
  index = names(dissimilarity)[3],
  method = "average",
  randomize = TRUE,
  n_runs = 100,
  keep_trials = "no",
  optimal_tree_method = "iterative_consensus_tree",
  n_clust = NULL,
  cut_height = NULL,
  find_h = TRUE,
  h_max = 1,
  h_min = 0,
  consensus_p = 0.5,
  show_hierarchy = FALSE,
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
  represent pairs of sites (or any pair of nodes), and the subsequent
  column(s) contain the dissimilarity indices.

- index:

  The name or number of the dissimilarity column to use. By default, the
  third column name of `dissimilarity` is used.

- method:

  The name of the hierarchical classification method, as in
  [hclust](https://rdrr.io/pkg/fastcluster/man/hclust.html). Should be
  one of `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`
  (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC), or
  `"centroid"` (= UPGMC).

- randomize:

  A `boolean` indicating whether the dissimilarity matrix should be
  randomized to account for the order of sites in the dissimilarity
  matrix.

- n_runs:

  The number of trials for randomizing the dissimilarity matrix.

- keep_trials:

  A `character` string indicating whether random trial results
  (including the randomized matrix, the associated tree and metrics for
  that tree) should be stored in the output object. Possible values are
  `"no"` (default), `"all"` or `"metrics"`. Note that this parameter is
  automatically set to `"no"` if
  `optimal_tree_method = "iterative_consensus_tree"`.

- optimal_tree_method:

  A `character` string indicating how the final tree should be obtained
  from all trials. Possible values are `"iterative_consensus_tree"`
  (default), `"best"` or `"consensus"`. **We recommend
  `"iterative_consensus_tree"`. See Details.**

- n_clust:

  An `integer` vector or a single `integer` indicating the number of
  clusters to be obtained from the hierarchical tree, or the output from
  [bioregionalization_metrics](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md).
  This parameter should not be used simultaneously with `cut_height`.

- cut_height:

  A `numeric` vector indicating the height(s) at which the tree should
  be cut. This parameter should not be used simultaneously with
  `n_clust`.

- find_h:

  A `boolean` indicating whether the height of the cut should be found
  for the requested `n_clust`.

- h_max:

  A `numeric` value indicating the maximum possible tree height for the
  chosen `index`.

- h_min:

  A `numeric` value indicating the minimum possible height in the tree
  for the chosen `index`.

- consensus_p:

  A `numeric` value (applicable only if
  `optimal_tree_method = "consensus"`) indicating the threshold
  proportion of trees that must support a region/cluster for it to be
  included in the final consensus tree.

- show_hierarchy:

  A `boolean` specifying if the hierarchy of clusters should be
  identifiable in the outputs (`FALSE` by default). This argument is
  only used if the tree is cut (i.e., `n_clust` or `cut_height` is
  provided).

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

In the `algorithm` slot, users can find the following elements:

- `trials`: A list containing all randomization trials. Each trial
  includes the dissimilarity matrix with randomized site order, the
  associated tree, and the cophenetic correlation coefficient for that
  tree.

- `final.tree`: An `hclust` object representing the final hierarchical
  tree to be used.

- `final.tree.coph.cor`: The cophenetic correlation coefficient between
  the initial dissimilarity matrix and the `final.tree`.

## Details

The function is based on
[hclust](https://rdrr.io/pkg/fastcluster/man/hclust.html). The default
method for the hierarchical tree is `average`, i.e. UPGMA as it has been
recommended as the best method to generate a tree from beta diversity
dissimilarity (Kreft & Jetz, 2010).

Clusters can be obtained by two methods:

- Specifying a desired number of clusters in `n_clust`

- Specifying one or several heights of cut in `cut_height`

To find an optimal number of clusters, see
[`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)

It is important to pay attention to the fact that the order of rows in
the input distance matrix influences the tree topology as explained in
Dapporto (2013). To address this, the function generates multiple trees
by randomizing the distance matrix.

Two methods are available to obtain the final tree:

- `optimal_tree_method = "iterative_consensus_tree"`: The Iterative
  Hierarchical Consensus Tree (IHCT) method reconstructs a consensus
  tree by iteratively splitting the dataset into two subclusters based
  on the pairwise dissimilarity of sites across `n_runs` trees based on
  `n_runs` randomizations of the distance matrix. At each iteration, it
  identifies the majority membership of sites into two stable groups
  across all trees, calculates the height based on the selected linkage
  method (`method`), and enforces monotonic constraints on node heights
  to produce a coherent tree structure. This approach provides a robust,
  hierarchical representation of site relationships, balancing cluster
  stability and hierarchical constraints.

- `optimal_tree_method = "best"`: This method selects one tree among
  with the highest cophenetic correlation coefficient, representing the
  best fit between the hierarchical structure and the original distance
  matrix.

- `optimal_tree_method = "consensus"`: This method constructs a
  consensus tree using phylogenetic methods with the function
  [consensus](https://rdrr.io/pkg/ape/man/consensus.html). When using
  this option, you must set the `consensus_p` parameter, which indicates
  the proportion of trees that must contain a region/cluster for it to
  be included in the final consensus tree. Consensus trees lack an
  inherent height because they represent a majority structure rather
  than an actual hierarchical clustering. To assign heights, we use a
  non-negative least squares method
  ([nnls.tree](https://klausvigo.github.io/phangorn/reference/designTree.html))
  based on the initial distance matrix, ensuring that the consensus tree
  preserves approximate distances among clusters.

We recommend using the `"iterative_consensus_tree"` as all the branches
of this tree will always reflect the majority decision among many
randomized versions of the distance matrix. This method is inspired by
Dapporto et al. (2015), which also used the majority decision among many
randomized versions of the distance matrix, but it expands it to
reconstruct the entire topology of the tree iteratively.

We do not recommend using the basic `consensus` method because in many
contexts it provides inconsistent results, with a meaningless tree
topology and a very low cophenetic correlation coefficient.

For a fast exploration of the tree, we recommend using the `best` method
which will only select the tree with the highest cophenetic correlation
coefficient among all randomized versions of the distance matrix.

## References

Kreft H & Jetz W (2010) A framework for delineating biogeographical
regions based on species distributions. *Journal of Biogeography* 37,
2029-2053.

Dapporto L, Ramazzotti M, Fattorini S, Talavera G, Vila R & Dennis, RLH
(2013) Recluster: an unbiased clustering procedure for beta-diversity
turnover. *Ecography* 36, 1070–1075.

Dapporto L, Ciolli G, Dennis RLH, Fox R & Shreeve TG (2015) A new
procedure for extrapolating turnover regionalization at mid-small
spatial scales, tested on British butterflies. *Methods in Ecology and
Evolution* 6 , 1287–1297.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html>.

Associated functions:
[cut_tree](https://bioRgeo.github.io/bioregion/reference/cut_tree.md)

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

dissim <- dissimilarity(comat, metric = "Simpson")

# User-defined number of clusters
tree1 <- hclu_hierarclust(dissim, 
                          n_clust = 5)
#> Building the iterative hierarchical consensus tree... Note that this process can take time especially if you have a lot of sites.
#> Error in if (nrow(dist_matrix) > 2) {    subclusters <- stable_binary_split(dist_matrix, method = method,         n_runs = n_runs, binsplit = "tree")} else {    subclusters <- list(attr(dist_matrix, "Labels")[1], attr(dist_matrix,         "Labels")[2])}: argument is of length zero
tree1
#> Error: object 'tree1' not found
plot(tree1)
#> Error: object 'tree1' not found
str(tree1)
#> Error: object 'tree1' not found
tree1$clusters
#> Error: object 'tree1' not found

# User-defined height cut
# Only one height
tree2 <- hclu_hierarclust(dissim, 
                          cut_height = .05)
#> Building the iterative hierarchical consensus tree... Note that this process can take time especially if you have a lot of sites.
#> Error in if (nrow(dist_matrix) > 2) {    subclusters <- stable_binary_split(dist_matrix, method = method,         n_runs = n_runs, binsplit = "tree")} else {    subclusters <- list(attr(dist_matrix, "Labels")[1], attr(dist_matrix,         "Labels")[2])}: argument is of length zero
tree2
#> Error: object 'tree2' not found
tree2$clusters
#> Error: object 'tree2' not found

# Multiple heights
tree3 <- hclu_hierarclust(dissim, 
                          cut_height = c(.05, .15, .25))
#> Building the iterative hierarchical consensus tree... Note that this process can take time especially if you have a lot of sites.
#> Error in if (nrow(dist_matrix) > 2) {    subclusters <- stable_binary_split(dist_matrix, method = method,         n_runs = n_runs, binsplit = "tree")} else {    subclusters <- list(attr(dist_matrix, "Labels")[1], attr(dist_matrix,         "Labels")[2])}: argument is of length zero

tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#> Error: object 'tree3' not found
# Info on each partition can be found in table cluster_info
tree3$cluster_info
#> Error: object 'tree3' not found
plot(tree3)
#> Error: object 'tree3' not found
```
