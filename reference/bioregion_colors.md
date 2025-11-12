# Add color palettes to bioregion cluster objects

This function assigns colors to clusters in a `bioregion.clusters`
object using color palettes from the `rcartocolor` package. It handles
large numbers of clusters by assigning vivid colors to the most
important clusters (based on size), grey shades to less important
clusters, and optionally black to insignificant clusters.

## Usage

``` r
bioregion_colors(
  clusters,
  palette = "Vivid",
  cluster_ordering = "n_sites",
  cutoff_insignificant = NULL
)
```

## Arguments

- clusters:

  An object of class `bioregion.clusters`, typically output from
  clustering functions such as
  [`netclu_infomap()`](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md),
  [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md),
  or
  [`nhclu_pam()`](https://bioRgeo.github.io/bioregion/reference/nhclu_pam.md).

- palette:

  A `character` string indicating which color palette from `rcartocolor`
  to use. Default is `"Vivid"`. Other qualitative palettes include
  `"Bold"`, `"Prism"`, `"Safe"`, `"Antique"`, and `"Pastel"`.

- cluster_ordering:

  A `character` string indicating the criterion for ranking clusters to
  determine color assignment priority. Options are:

  - `"n_sites"` (default): Rank by number of sites in each cluster

  - `"n_species"`: Rank by number of species (bipartite networks only)

  - `"n_both"`: Rank by combined sites + species (bipartite networks
    only)

  Larger clusters (by the chosen criterion) receive vivid colors first.

- cutoff_insignificant:

  A `numeric` value or `NULL` (default). When specified, clusters with
  values at or below this threshold (based on the `cluster_ordering`
  criterion) are considered insignificant and colored black, reducing
  visual clutter on maps. If `NULL`, all clusters receive distinct
  colors.

## Value

A modified `bioregion.clusters` object with two additional elements:

- `colors`: A `list` where each element corresponds to a partition
  (bioregionalization). Each list element is a `data.frame` with two
  columns:

  - `cluster` (`character`): Cluster identifier for that partition

  - `color` (`character`): Hex color code (e.g., "#FF5733")

- `clusters_colors`: A `data.frame` with the same structure as the
  `clusters` element, but with cluster IDs replaced by their
  corresponding hex color codes for direct use in plotting functions.

## Details

The function uses a two-step algorithm to assign colors:

**Step 1: Identify insignificant clusters** (if `cutoff_insignificant`
is specified)

Insignificant clusters are those with a marginal size compared to
others. This is a subjective threshold set by the user. All such
clusters are assigned the color black (#000000) to minimize their visual
impact. Clusters with values at or below the threshold are assigned
black (#000000).

**Step 2: Assign colors to significant clusters**

Remaining clusters are ranked by the `cluster_ordering` criterion:

- **Top clusters** (up to 12): Receive distinct colors from the chosen
  palette. This limit is because above 12 the human eye struggles to
  distinguish between colors.

- **Remaining clusters** (beyond top 12): Receive shades of grey from
  light (#CCCCCC) to dark (#404040), maintaining visual distinction but
  with less prominence.

**Multiple partitions**: If the cluster object contains multiple
partitions (e.g., from hierarchical clustering with different k values),
colors are assigned independently for each partition. Each partition
gets its own color scale optimized for the number of clusters in that
partition.

## Note

The colored cluster object can be directly used with
[`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md),
which will automatically detect and apply the color scheme when present.

## References

Color palettes from the `rcartocolor` package: Nowosad J (2018).
"CARTOColors: color palettes inspired by CARTO."
<https://github.com/Nowosad/rcartocolor>

## See also

[`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
for visualizing colored clusters on maps

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
data(fishmat)
data(fishsf)

# Basic example with few clusters
sim <- similarity(fishmat, metric = "Simpson")
clust <- netclu_greedy(sim)
clust_colored <- bioregion_colors(clust)
print(clust_colored)
#> Clustering results for algorithm : netclu_greedy 
#>  - Number of sites:  338 
#> Clustering results:
#>  - Number of partitions:  1 
#>  - Number of clusters:  4 
#>  - Color palette assigned:
#>    *  K_4 : 4 vivid colors 

if (FALSE) { # \dontrun{
# Map with automatic colors
map_bioregions(clust_colored, fishsf)

# Example with many clusters and cutoff
dissim <- similarity_to_dissimilarity(sim)
clust <- hclu_hierarclust(dissim,
                          optimal_tree_method = "best",
                          n_clust = 15)
clust_colored2 <- bioregion_colors(clust, 
                                   cluster_ordering = "n_sites",
                                   cutoff_insignificant = 1)
map_bioregions(clust_colored2, fishsf)

# Example with different palette
clust_colored3 <- bioregion_colors(clust, palette = "Bold")
map_bioregions(clust_colored3, fishsf)


# Example with bipartite network
clust_bip <- netclu_greedy(fishdf, bipartite = TRUE)
clust_bip_colored <- bioregion_colors(clust_bip, 
                                      cluster_ordering = "n_both")
map_bioregions(clust_bip_colored, fishsf)
                                      
} # }
```
