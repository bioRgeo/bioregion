# Infomap community finding

This function finds communities in a (un)weighted (un)directed network
based on the Infomap algorithm
(<https://github.com/mapequation/infomap>).

## Usage

``` r
netclu_infomap(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  seed = NULL,
  nbmod = 0,
  markovtime = 1,
  numtrials = 1,
  twolevel = FALSE,
  show_hierarchy = FALSE,
  directed = FALSE,
  bipartite_version = FALSE,
  bipartite = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  version = "2.8.0",
  binpath = "tempdir",
  check_install = TRUE,
  path_temp = "infomap_temp",
  delete_temp = TRUE
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

  The seed for the random number generator (`NULL` for random by
  default).

- nbmod:

  Penalize solutions the more they differ from this number (`0` by
  default for no preferred number of modules).

- markovtime:

  Scales link flow to change the cost of moving between modules, higher
  values result in fewer modules (`1` by default).

- numtrials:

  For the number of trials before picking up the best solution.

- twolevel:

  A `boolean` indicating if the algorithm should optimize a two-level
  partition of the network (`FALSE` by default for multi-level).

- show_hierarchy:

  A `boolean` specifying if the hierarchy of community should be
  identifiable in the outputs (`FALSE` by default).

- directed:

  A `boolean` indicating if the network is directed (from column 1 to
  column 2).

- bipartite_version:

  A `boolean` indicating if the bipartite version of Infomap should be
  used (see Note).

- bipartite:

  A `boolean` indicating if the network is bipartite (see Note).

- site_col:

  The name or number for the column of site nodes (i.e. primary nodes).

- species_col:

  The name or number for the column of species nodes (i.e. feature
  nodes).

- return_node_type:

  A `character` indicating what types of nodes (`"site"`, `"species"`,
  or `"both"`) should be returned in the output (`"both"` by default).

- version:

  A `character` indicating the Infomap version to use.

- binpath:

  A `character` indicating the path to the bin folder (see
  [install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  and Details).

- check_install:

  A `boolean` indicating if the function should check that the Infomap
  has been properly installed (see
  [install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  and Details).

- path_temp:

  A `character` indicating the path to the temporary folder (see
  Details).

- delete_temp:

  A `boolean` indicating if the temporary folder should be removed (see
  Details).

## Value

A `list` of class `bioregion.clusters` with five slots:

1.  **name**: A `character` containing the name of the algorithm.

2.  **args**: A `list` of input arguments as provided by the user.

3.  **inputs**: A `list` of characteristics of the clustering process.

4.  **algorithm**: A `list` of all objects associated with the
    clustering procedure, such as original cluster objects.

5.  **clusters**: A `data.frame` containing the clustering results.

In the `algorithm` slot, users can find the following elements:

- `cmd`: The command line used to run Infomap.

- `version`: The Infomap version.

- `web`: Infomap's GitHub repository.

## Details

Infomap is a network clustering algorithm based on the Map equation
proposed in Rosvall & Bergstrom (2008) that finds communities in
(un)weighted and (un)directed networks.

This function is based on the C++ version of Infomap
(<https://github.com/mapequation/infomap/releases>). This function needs
binary files to run. They can be installed with
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md).

**If you changed the default path to the `bin` folder while running
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
PLEASE MAKE SURE to set `binpath` accordingly.**

**If you did not use
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
to change the permissions and test the binary files PLEASE MAKE SURE to
set `check_install` accordingly.**

The C++ version of Infomap generates temporary folders and/or files that
are stored in the `path_temp` folder ("infomap_temp" with a unique
timestamp located in the bin folder in `binpath` by default). This
temporary folder is removed by default (`delete_temp = TRUE`).

Several versions of Infomap are available in the package. See
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
for more details.

## Note

Infomap has been designed to deal with bipartite networks. To use this
functionality, set the `bipartite_version` argument to TRUE in order to
approximate a two-step random walker (see
<https://www.mapequation.org/infomap/> for more information). Note that
a bipartite network can also be considered as a unipartite network
(`bipartite = TRUE`).

In both cases, do not forget to indicate which of the first two columns
is dedicated to the site nodes (i.e., primary nodes) and species nodes
(i.e. feature nodes) using the arguments `site_col` and `species_col`.
The type of nodes returned in the output can be chosen with the argument
`return_node_type` equal to `"both"` to keep both types of nodes,
`"site"` to preserve only the site nodes, and `"species"` to preserve
only the species nodes.

## References

Rosvall M & Bergstrom CT (2008) Maps of random walks on complex networks
reveal community structure. *Proceedings of the National Academy of
Sciences* 105, 1118-1123.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html>.

Associated functions:
[netclu_greedy](https://bioRgeo.github.io/bioregion/reference/netclu_greedy.md)
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
com <- netclu_infomap(net)
#> Infomap 2.8.0 is not installed... Please have a look at https//bioRgeo.github.io/bioregion/articles/a1_install_binary_files.html for more details.
#> It should be located in /tmp/Rtmpxc9hSF/bin/INFOMAP/2.8.0/
```
