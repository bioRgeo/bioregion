# Finding communities using the Leiden algorithm

This function finds communities in a (un)weighted undirected network
based on the Leiden algorithm of Traag, van Eck & Waltman.

## Usage

``` r
netclu_leiden(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  seed = NULL,
  objective_function = "CPM",
  resolution_parameter = 1,
  beta = 0.01,
  n_iterations = 2,
  vertex_weights = NULL,
  bipartite = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  algorithm_in_output = TRUE
)
```

## Arguments

- net:

  The output object from
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md).
  If a `data.frame` is used, the first two columns represent pairs of
  sites (or any pair of nodes), and the next column(s) are the
  similarity indices.

- weight:

  A `boolean` indicating if the weights should be considered if there
  are more than two columns.

- cut_weight:

  A minimal weight value. If `weight` is TRUE, the links between sites
  with a weight strictly lower than this value will not be considered
  (`0` by default).

- index:

  The name or number of the column to use as weight. By default, the
  third column name of `net` is used.

- seed:

  The random number generator seed (NULL for random by default).

- objective_function:

  A string indicating the objective function to use, either the Constant
  Potts Model ("CPM") or "modularity" ("CPM" by default).

- resolution_parameter:

  The resolution parameter to use. Higher resolutions lead to smaller
  communities, while lower resolutions lead to larger communities.

- beta:

  A parameter affecting the randomness in the Leiden algorithm. This
  affects only the refinement step of the algorithm.

- n_iterations:

  The number of iterations for the Leiden algorithm. Each iteration may
  further improve the partition.

- vertex_weights:

  The vertex weights used in the Leiden algorithm. If not provided, they
  will be automatically determined based on the objective_function.
  Please see the details of this function to understand how to interpret
  the vertex weights.

- bipartite:

  A `boolean` indicating if the network is bipartite (see Details).

- site_col:

  The name or number for the column of site nodes (i.e., primary nodes).

- species_col:

  The name or number for the column of species nodes (i.e., feature
  nodes).

- return_node_type:

  A `character` indicating what types of nodes ("site", "species", or
  "both") should be returned in the output (`"both"` by default).

- algorithm_in_output:

  A `boolean` indicating if the original output of
  [cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)
  should be returned in the output (`TRUE` by default, see Value).

## Value

A `list` of class `bioregion.clusters` with five slots:

1.  **name**: A `character` containing the name of the algorithm.

2.  **args**: A `list` of input arguments as provided by the user.

3.  **inputs**: A `list` of characteristics of the clustering process.

4.  **algorithm**: A `list` of all objects associated with the
    clustering procedure, such as original cluster objects (only if
    `algorithm_in_output = TRUE`).

5.  **clusters**: A `data.frame` containing the clustering results.

In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can find
the output of
[cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html).

## Details

This function is based on the Leiden algorithm (Traag et al., 2019) as
implemented in the [igraph](https://cran.r-project.org/package=igraph)
package
([cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)).

## Note

Although this algorithm was not primarily designed to deal with
bipartite networks, it is possible to consider the bipartite network as
a unipartite network (`bipartite = TRUE`).

Do not forget to indicate which of the first two columns is dedicated to
the site nodes (i.e., primary nodes) and species nodes (i.e. feature
nodes) using the arguments `site_col` and `species_col`. The type of
nodes returned in the output can be chosen with the argument
`return_node_type` equal to `"both"` to keep both types of nodes,
`"site"` to preserve only the site nodes, and `"species"` to preserve
only the species nodes.

## References

Traag VA, Waltman L & Van Eck NJ (2019) From Louvain to Leiden:
guaranteeing well-connected communities. *Scientific reports* 9, 5233.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html>.

Associated functions:
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.md)
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

net <- similarity(comat, metric = "Simpson")
com <- netclu_leiden(net)

net_bip <- mat_to_net(comat, weight = TRUE)
clust2 <- netclu_leiden(net_bip, bipartite = TRUE)
```
