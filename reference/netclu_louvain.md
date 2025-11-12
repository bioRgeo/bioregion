# Louvain community finding

This function finds communities in a (un)weighted undirected network
based on the Louvain algorithm.

## Usage

``` r
netclu_louvain(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  lang = "igraph",
  resolution = 1,
  seed = NULL,
  q = 0,
  c = 0.5,
  k = 1,
  bipartite = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  binpath = "tempdir",
  check_install = TRUE,
  path_temp = "louvain_temp",
  delete_temp = TRUE,
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

- lang:

  A string indicating which version of Louvain should be used
  (`"igraph"` or `"cpp"`, see Details).

- resolution:

  A resolution parameter to adjust the modularity (1 is chosen by
  default, see Details).

- seed:

  The random number generator seed (only when `lang = "igraph"`, NULL
  for random by default).

- q:

  The quality function used to compute the partition of the graph
  (modularity is chosen by default, see Details).

- c:

  The parameter for the Owsinski-Zadrozny quality function (between 0
  and 1, 0.5 is chosen by default).

- k:

  The kappa_min value for the Shi-Malik quality function (it must be \>
  0, 1 is chosen by default).

- bipartite:

  A `boolean` indicating if the network is bipartite (see Details).

- site_col:

  The name or number for the column of site nodes (i.e., primary nodes).

- species_col:

  The name or number for the column of species nodes (i.e., feature
  nodes).

- return_node_type:

  A `character` indicating what types of nodes (`"site"`, `"species"`,
  or `"both"`) should be returned in the output (`"both"` by default).

- binpath:

  A `character` indicating the path to the bin folder (see
  [install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  and Details).

- check_install:

  A `boolean` indicating if the function should check that Louvain has
  been properly installed (see
  [install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  and Details).

- path_temp:

  A `character` indicating the path to the temporary folder (see
  Details).

- delete_temp:

  A `boolean` indicating if the temporary folder should be removed (see
  Details).

- algorithm_in_output:

  A `boolean` indicating if the original output of
  [cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html)
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
[cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html)
if `lang = "igraph"` and the following element if `lang = "cpp"`:

- `cmd`: The command line used to run Louvain.

- `version`: The Louvain version.

- `web`: The Louvain's website.

## Details

Louvain is a network community detection algorithm proposed in (Blondel
et al., 2008). This function offers two implementations of the Louvain
algorithm (controlled by the `lang` parameter): the
[igraph](https://cran.r-project.org/package=igraph) implementation
([cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html))
and the C++ implementation (<https://sourceforge.net/projects/louvain/>,
version 0.3).

The [igraph](https://cran.r-project.org/package=igraph) implementation
allows adjustment of the resolution parameter of the modularity function
(`resolution` argument) used internally by the algorithm. Lower values
typically yield fewer, larger clusters. The original definition of
modularity is recovered when the resolution parameter is set to 1 (by
default).

The C++ implementation provides several quality functions: `q = 0` for
the classical Newman-Girvan criterion (Modularity), `q = 1` for the
Zahn-Condorcet criterion, `q = 2` for the Owsinski-Zadrozny criterion
(parameterized by `c`), `q = 3` for the Goldberg Density criterion,
`q = 4` for the A-weighted Condorcet criterion, `q = 5` for the
Deviation to Indetermination criterion, `q = 6` for the Deviation to
Uniformity criterion, `q = 7` for the Profile Difference criterion,
`q = 8` for the Shi-Malik criterion (parameterized by `k`), and `q = 9`
for the Balanced Modularity criterion.

The C++ version is based on version 0.3
(<https://sourceforge.net/projects/louvain/>). Binary files are required
to run it, and can be installed with
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md).

**If you changed the default path to the `bin` folder while running
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md),
PLEASE MAKE SURE to set `binpath` accordingly.**

**If you did not use
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
to change the permissions or test the binary files, PLEASE MAKE SURE to
set `check_install` accordingly.**

The C++ version generates temporary folders and/or files in the
`path_temp` folder ("louvain_temp" with a unique timestamp located in
the bin folder in `binpath` by default). This temporary folder is
removed by default (`delete_temp = TRUE`).

## Note

Although this algorithm was not primarily designed to deal with
bipartite networks, it is possible to consider the bipartite network as
a unipartite network (`bipartite = TRUE`).

Do not forget to indicate which of the first two columns is dedicated to
the site nodes (i.e., primary nodes) and species nodes (i.e., feature
nodes) using the arguments `site_col` and `species_col`. The type of
nodes returned in the output can be chosen with the argument
`return_node_type` equal to `"both"` to keep both types of nodes,
`"site"` to preserve only the site nodes, and `"species"` to preserve
only the species nodes.

## References

Blondel VD, Guillaume JL, Lambiotte R & Mech ELJS (2008) Fast unfolding
of communities in large networks. *J. Stat. Mech.* 10, P10008.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html>.

Associated functions:
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)
[netclu_greedy](https://bioRgeo.github.io/bioregion/reference/netclu_greedy.md)
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
com <- netclu_louvain(net, lang = "igraph")
```
