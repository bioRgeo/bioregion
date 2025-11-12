# Community structure detection in weighted bipartite networks via modularity optimization

This function takes a bipartite weighted graph and computes modules by
applying Newman’s modularity measure in a bipartite weighted version.

## Usage

``` r
netclu_beckett(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  seed = NULL,
  forceLPA = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  algorithm_in_output = TRUE
)
```

## Arguments

- net:

  A `data.frame` representing a bipartite network with the first two
  columns representing undirected links between pairs of nodes, and the
  next column(s) representing the weights of the links.

- weight:

  A `boolean` indicating whether weights should be considered if there
  are more than two columns (see Note).

- cut_weight:

  A minimal weight value. If `weight` is TRUE, links with weights
  strictly lower than this value will not be considered (`0` by
  default).

- index:

  The name or number of the column to use as weight. By default, the
  third column name of `net` is used.

- seed:

  The seed for the random number generator (`NULL` for random by
  default).

- forceLPA:

  A `boolean` indicating whether the even faster pure LPA-algorithm of
  Beckett should be used. DIRT-LPA (the default) is less likely to get
  trapped in a local minimum but is slightly slower. Defaults to
  `FALSE`.

- site_col:

  The name or number of the column for site nodes (i.e., primary nodes).

- species_col:

  The name or number of the column for species nodes (i.e., feature
  nodes).

- return_node_type:

  A `character` indicating which types of nodes (`"site"`, `"species"`,
  or `"both"`) should be returned in the output (`"both"` by default).

- algorithm_in_output:

  A `boolean` indicating whether the original output of
  [computeModules](https://rdrr.io/pkg/bipartite/man/computeModules.html)
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

If `algorithm_in_output = TRUE`, users can find the output of
[computeModules](https://rdrr.io/pkg/bipartite/man/computeModules.html)
in the `algorithm` slot.

## Details

This function is based on the modularity optimization algorithm provided
by Stephen Beckett (Beckett, 2016) as implemented in the
[bipartite](https://cran.r-project.org/package=bipartite) package
([computeModules](https://rdrr.io/pkg/bipartite/man/computeModules.html)).

## Note

Beckett's algorithm is designed to handle weighted bipartite networks.
If `weight = FALSE`, a weight of 1 will be assigned to each pair of
nodes. Ensure that the `site_col` and `species_col` arguments correctly
identify the respective columns for site nodes (primary nodes) and
species nodes (feature nodes). The type of nodes returned in the output
can be selected using the `return_node_type` argument: `"both"` to
include both node types, `"site"` to return only site nodes, or
`"species"` to return only species nodes.

## References

Beckett SJ (2016) Improved community detection in weighted bipartite
networks. *Royal Society Open Science* 3, 140536.

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
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

com <- netclu_beckett(net)
```
