# Extract a subset of sites or species from a `bioregion.clusters` object

This function extracts a subset of nodes based on their type (`"site"`
or `"species"`) from a `bioregion.clusters` object, which contains both
types of nodes (sites and species).

## Usage

``` r
site_species_subset(clusters, node_type = "site")
```

## Arguments

- clusters:

  An object of class `bioregion.clusters`.

- node_type:

  A `character` string indicating the type of nodes to extract. Possible
  values are `"site"` or `"species"`. The default is `"site"`.

## Value

An object of class `bioregion.clusters` containing only the specified
node type (sites or species).

## Note

Some `bioregion.clusters` objects may contain both types of nodes (sites
and species). This information is available in the `$inputs$node_type`
slot.

This function allows you to extract a specific type of node (either
sites or species) from any `bioregion.clusters` object that includes
both.

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

clusters <- netclu_louvain(net, lang = "igraph", bipartite = TRUE)

clusters_sites <- site_species_subset(clusters, node_type = "site")
```
