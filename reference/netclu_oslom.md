# OSLOM community finding

This function finds communities in a (un)weighted (un)directed network
based on the OSLOM algorithm (<http://oslom.org/>, version 2.4).

## Usage

``` r
netclu_oslom(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  seed = NULL,
  reassign = "no",
  r = 10,
  hr = 50,
  t = 0.1,
  cp = 0.5,
  directed = FALSE,
  bipartite = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  binpath = "tempdir",
  check_install = TRUE,
  path_temp = "oslom_temp",
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
  with a weight strictly lower than this value will not be considered (0
  by default).

- index:

  Name or number of the column to use as weight. By default, the third
  column name of `net` is used.

- seed:

  For the random number generator (NULL for random by default).

- reassign:

  A `character` indicating if the nodes belonging to several community
  should be reassigned and what method should be used (see Note).

- r:

  The number of runs for the first hierarchical level (10 by default).

- hr:

  The number of runs for the higher hierarchical level (50 by default, 0
  if you are not interested in hierarchies).

- t:

  The p-value, the default value is 0.10. Increase this value if you
  want more modules.

- cp:

  Kind of resolution parameter used to decide between taking some
  modules or their union (default value is 0.5; a bigger value leads to
  bigger clusters).

- directed:

  A `boolean` indicating if the network is directed (from column 1 to
  column 2).

- bipartite:

  A `boolean` indicating if the network is bipartite (see Details).

- site_col:

  Name or number for the column of site nodes (i.e. primary nodes).

- species_col:

  Name or number for the column of species nodes (i.e. feature nodes).

- return_node_type:

  A `character` indicating what types of nodes (`site`, `species`, or
  `both`) should be returned in the output (`return_node_type = "both"`
  by default).

- binpath:

  A `character` indicating the path to the bin folder (see
  [install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  and Details).

- check_install:

  A `boolean` indicating if the function should check that the OSLOM has
  been properly installed (see
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
    clustering procedure, such as original cluster objects (only if
    `algorithm_in_output = TRUE`).

5.  **clusters**: A `data.frame` containing the clustering results.

In the `algorithm` slot, users can find the following elements:

- `cmd`: The command line used to run OSLOM.

- `version`: The OSLOM version.

- `web`: The OSLOM's web site.

## Details

OSLOM is a network community detection algorithm proposed in
Lancichinetti et al. (2011) that finds statistically significant
(overlapping) communities in (un)weighted and (un)directed networks.

This function is based on the 2.4 C++ version of OSLOM
(<http://www.oslom.org/software.htm>). This function needs files to run.
They can be installed with
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md).

**If you changed the default path to the `bin` folder while running
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md),
PLEASE MAKE SURE to set `binpath` accordingly.**

**If you did not use
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
to change the permissions and test the binary files, PLEASE MAKE SURE to
set `check_install` accordingly.**

The C++ version of OSLOM generates temporary folders and/or files that
are stored in the `path_temp` folder (folder "oslom_temp" with a unique
timestamp located in the bin folder in `binpath` by default). This
temporary folder is removed by default (`delete_temp = TRUE`).

## Note

Although this algorithm was not primarily designed to deal with
bipartite networks, it is possible to consider the bipartite network as
unipartite network (`bipartite = TRUE`). Do not forget to indicate which
of the first two columns is dedicated to the site nodes (i.e. primary
nodes) and species nodes (i.e. feature nodes) using the arguments
`site_col` and `species_col`. The type of nodes returned in the output
can be chosen with the argument `return_node_type` equal to `both` to
keep both types of nodes, `sites` to preserve only the sites nodes, and
`species` to preserve only the species nodes.

Since OSLOM potentially returns overlapping communities, we propose two
methods to reassign the 'overlapping' nodes: randomly
(`reassign = "random"`) or based on the closest candidate community
(`reassign = "simil"`) (only for weighted networks, in this case the
closest candidate community is determined with the average similarity).
By default, `reassign = "no"` and all the information will be provided.
The number of partitions will depend on the number of overlapping
modules (up to three). The suffix `_semel`, `_bis`, and `_ter` are added
to the column names. The first partition (`_semel`) assigns a module to
each node. A value of `NA` in the second (`_bis`) and third (`_ter`)
columns indicates that no overlapping module was found for this node
(i.e. non-overlapping nodes).

## References

Lancichinetti A, Radicchi F, Ramasco JJ & Fortunato S (2011) Finding
statistically significant communities in networks. *PLOS ONE* 6, e18961.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html>.

Associated functions:
[netclu_greedy](https://bioRgeo.github.io/bioregion/reference/netclu_greedy.md)
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.md)

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
com <- netclu_oslom(net)
#> OSLOM is not installed... Please have a look at https://bioRgeo.github.io/bioregion/articles/a1_install_binary_files.html for more details.
#> It should be located in /tmp/Rtmp90CR88/bin/OSLOM/
```
