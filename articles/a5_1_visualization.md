# 5.1 Visualization

In this vignette, we aim at illustrating how to plot the bioregions
identified with the diferent algorithms available in `bioregion`.  
Using datasets coming along with `bioregion`, we show different methods
to plot your results.

## 1. Data

For this vignette, we mostly rely on the dataset describing the
distribution of vascular plants in the Mediterranean part of France. We
first load the matrix format of this dataset, computes the dissimilarity
matrix out of it and also load the data.frame format of the data.

``` r

data(vegedf)
data(vegemat)

vegedissim <- dissimilarity(vegemat, metric = "Simpson")
```

Since we aim at plotting the result, we also need the object `vegesf`
linking each site of the dataset to a geometry.

``` r

data(vegesf)
```

We also import the world coastlines, available from the `rnaturalearth`
R package.

``` r

world <- rnaturalearth::ne_coastline(returnclass = "sf", scale = "medium")
# Align the CRS of both objects
vegesf <- st_transform(vegesf, crs = st_crs(world))
```

## 2. Summary of visualization methods in the bioregion R package

The package provides functions for three visualization methods:

1.  Maps (directly in R):
    [`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md):
    [section 3](#maps)
2.  Bar plots and trees of clusters (directly in R)
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html): [section
    4](#clusterplots)
3.  Export of networks (for visualization in external software Gephi)
    `exportGDF`: [section 5](#networkplots)

To facilitate graphical analysis, the package also provides a function
to make sure colors are consistent across all visualizations:
`bioregion_colors`. This function assigns fixed colors to bioregions and
ensures that these colors remain consistent across different types of
plots. The same bioregion should have the same color in all your figures
(maps, summary plots, networks, etc.).

#### Example visualization with a colorblind-safe palette

``` r

vege_nhclu_kmeans <- nhclu_kmeans(vegedissim,
                                  n_clust = 3)


# Apply the Safe palette
vege_kmeans_safe <- bioregion_colors(vege_nhclu_kmeans, 
                                     palette = "Safe")

# View the color assignments
vege_kmeans_safe$colors
```

    ## $K_3
    ##   cluster   color
    ## 1       1 #88CCEE
    ## 2       2 #CC6677
    ## 3       3 #DDCC77

``` r

# 1. Map with Safe colors
map_bioregions(vege_kmeans_safe, 
               map = vegesf)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-4-1.png)

``` r

# 2. Bar plot of clusters
plot(vege_kmeans_safe)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-4-2.png)

Note that `vege_kmeans_safe` is identical to `vege_kmeans` with two new
slots: \* `colors`: a small data.frame summarising the color assigned to
each bioregion \* `clusters_colors`: a data.frame identical to the
clusters slot, but with content modified to include colors

The same color scheme can also be applied when exporting networks to
Gephi. The
[`exportGDF()`](https://bioRgeo.github.io/bioregion/reference/exportGDF.md)
function automatically uses the colors defined by
[`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md),
ensuring consistency between your maps and network visualizations.

## 3. Maps

### 3.1 `map_bioregions()`

The simplest method is to use the function
[`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
from the package. This function can directly provide a plot of each site
colored according to the cluster they belong to.

Let’s take an example with a K-means clustering, with a number of
clusters set to 3.

``` r

vege_nhclu_kmeans <- nhclu_kmeans(vegedissim, 
                                  n_clust = 3, 
                                  index = "Simpson")
```

[`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
function can now simply takes the object `vege_nhclu_kmeans`, which is
of `bioregion.clusters` class, and the spatial distribution of sites,
stored in `vegesf`.

``` r

map_bioregions(vege_nhclu_kmeans, 
               map = vegesf)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-6-1.png)

If we define the color set beforehand, the maps will follow the chosen
colors:

``` r

# Notice how we overwrite the object when using bioregion_colors.
# In reality, it is just appending colors to the object
vege_nhclu_kmeans <- bioregion_colors(vege_nhclu_kmeans,
                                      palette = "Antique")
map_bioregions(vege_nhclu_kmeans, 
               map = vegesf)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-7-1.png)

### 3.2 Custom plot

If you want to customize yourself the plot and not simply rely on the
default option,
[`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
gives you the possibility to extract each site as well as its geometry
and cluster number.  
For this purpose, you can set the arguments like in the chunk below:

``` r

custom <- map_bioregions(vege_nhclu_kmeans, 
                         map = vegesf,
                         map_as_output = TRUE, 
                         plot = FALSE)
# Notice the partition is in output, as well as the colors
custom
```

    ## Simple feature collection with 715 features and 3 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 1.686171 ymin: 42.29604 xmax: 7.798711 ymax: 45.13742
    ## Geodetic CRS:  WGS 84
    ## First 10 features:
    ##    Site K_3 K_3_color                       geometry
    ## 1    35   1   #AF6458 POLYGON ((6.098099 45.13742...
    ## 2    36   1   #AF6458 POLYGON ((6.22521 45.13381,...
    ## 3    37   1   #AF6458 POLYGON ((6.352304 45.13007...
    ## 4    38   1   #AF6458 POLYGON ((6.47938 45.12617,...
    ## 5    39   1   #AF6458 POLYGON ((6.606438 45.12213...
    ## 6    84   1   #AF6458 POLYGON ((6.093117 45.04745...
    ## 7    85   1   #AF6458 POLYGON ((6.220024 45.04385...
    ## 8    86   1   #AF6458 POLYGON ((6.346914 45.04011...
    ## 9    87   1   #AF6458 POLYGON ((6.473786 45.03622...
    ## 10   88   1   #AF6458 POLYGON ((6.600641 45.03219...

``` r

# Crop world coastlines to the extent of the sf object of interest
europe <- sf::st_crop(world, sf::st_bbox(custom))


# Plot
# To get identical colors identical to other plots, we
# fill with the color and use scale_fill_identity()
ggplot(custom) +
  geom_sf(aes(fill = K_3_color), show.legend = FALSE) +
  geom_sf(data = europe) +
  scale_fill_identity() +
  labs(title = "Kmeans with 3 clusters") +
  theme_void()
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-8-1.png)

Another option which adds a legend consists in filling the sites by
bioregion (`K_3`) instead of by color, and give `scale_fill_manual()`
the color of each bioregion:

``` r

# Get the bioregion-color table for the partition we want to plot
# (here, K_3)
bioregion_cols <- vege_nhclu_kmeans$colors$K_3


# Named vector: bioregion -> color
legend_colors <- setNames(bioregion_cols$color, bioregion_cols$cluster)

ggplot(custom) +
  geom_sf(aes(fill = K_3)) +
  geom_sf(data = europe) +
  scale_fill_manual(values = legend_colors, name = "Bioregion") +
  labs(title = "Kmeans with 3 clusters") +
  theme_void()
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-9-1.png)

### 3.3 Plot with facets

Two methods are available:

- the simplest consists in using
  [`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
  directly:

``` r

vege_tree <- hclu_hierarclust(vegedissim, n_clust = c(3, 5, 8))
vege_tree <- bioregion_colors(vege_tree)
map_bioregions(vege_tree, map = vegesf)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-10-1.png)

- a more complex method consists in exporting the object and using it
  with ggplot. For this, you will need to transform from wide to long
  format, we provide an example code below:

``` r

vege_map <- map_bioregions(vege_tree, map = vegesf,
                           map_as_output = TRUE, plot = FALSE)
partitions <- colnames(vege_tree$clusters)[-1]

# One block of rows per partition, stacked
vege_long <- do.call(rbind, 
lapply(partitions, function(p) {
  data.frame(partition = p,
             bioregion = vege_map[[p]],
             color = vege_map[[paste0(p, "_color")]],
             geometry = sf::st_geometry(vege_map))
})) |> sf::st_as_sf()

ggplot(vege_long) +
  geom_sf(aes(fill = color), color = NA) +
  geom_sf(data = europe) +
  scale_fill_identity() +
  facet_wrap(~ partition) +
  theme_void()
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-11-1.png)

## 4. Plot cluster objects

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws three
types of plots for clustering results: trees (for hierarchical
clustering and nested partitions), bars showing the size of each cluster
(for all algorithms), and trees of the clusters, with the bar of each
cluster at the tip of its branch. The type is chosen automatically, or
with the argument `plot_type`. All plots use the colors of
[`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md),
so that they match the maps.

Three types of plots are available:

1.  Bar plots: plots showing the size of clusters, filled with one
    square per site (and one circle per species for bipartite clusters).
    This method works for all clustering methods.

![](a5_1_visualization_files/figure-html/unnamed-chunk-13-1.png)

2.  Tree plots: plots showing hierarchical trees, only working for
    hierarchical clustering methods
    ([`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md),
    [`hclu_diana()`](https://bioRgeo.github.io/bioregion/reference/hclu_diana.md),
    [`netclu_infomap()`](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)).

![](a5_1_visualization_files/figure-html/unnamed-chunk-14-1.png)

3.  Tree-cluster plots (or tree-bars): plots combining trees with bar
    plots, only for hierarchical clustering methods.

![](a5_1_visualization_files/figure-html/unnamed-chunk-15-1.png)

Plot types can be changed with argument `plot_type`. When
`plot_type = "auto"` (default), the plot is automatically chosen
depending on the clustering algorithm: bars for non-hierarchical and
trees for hierarchical clusterings.

### 4.1 Bar plots

#### 4.1.1 Simple bar plots (all cluster types)

For all algorithms, clusters can be shown as bars, one per cluster, with
one square per site. Clusters are sorted by decreasing number of sites:

``` r

vege_kmeans <- nhclu_kmeans(vegedissim, n_clust = 3)

plot(vege_kmeans)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-16-1.png)

It works when there are multiple partitions with multiple cuts:

``` r

vege_kmeans <- nhclu_kmeans(vegedissim, n_clust = c(3, 7, 14))

plot(vege_kmeans)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-17-1.png)

With many sites, each square can stand for several sites. This is chosen
automatically, or you can define it with `items_per_unit`.

``` r

plot(vege_kmeans, items_per_unit = 5, horizontal = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-18-1.png)

With many clusters, the smallest clusters beyond `max_clusters` (25 by
default) are grouped in a last, wide bar, `"other"`:

``` r

vege_kmeans_15 <- nhclu_kmeans(vegedissim, n_clust = 15)

plot(vege_kmeans_15, max_clusters = 8, horizontal = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-19-1.png)

With several partitions, each partition is drawn in its own panel. By
default, all panels share the same axis of the numbers of sites, so that
they can be compared; with `scales = "free"`, each panel has its own
axis:

``` r

plot(vege_kmeans)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-20-1.png)

``` r

plot(vege_kmeans, scales = "free")
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-20-2.png)

Titles and axes are set as in other plots of R, with `main`, `xlab`,
`ylab` and `axes`:

``` r

plot(vege_nhclu_kmeans, main = "k-means, 3 clusters", ylab = "Sites")
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-21-1.png)

#### 4.1.2 Bipartite bar plots (bipartite clustering)

In bipartite clustering, i.e. when both sites AND species are clustered,
we can distinguish on the bar plot the contents in terms of sites and
species.

Species are drawn as circles, below the sites by default
(`bipartite_layout = "mirror"`), above them in the same bar (`"stack"`),
or in a second bar (`"side"`):

``` r

data(fishdf)
fish_greedy <- netclu_greedy(fishdf, bipartite = TRUE)
plot(fish_greedy)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-22-1.png)

``` r

plot(fish_greedy, bipartite_layout = "stack")
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-22-2.png)

``` r

plot(fish_greedy, bipartite_layout = "side")
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-22-3.png)

#### 4.1.3 Plotting one partition with colors of another partition (try it if you do hierarchical clustering!!)

You can plot a partition with the colors of another partition,
e.g. showing how child partitions are grouped in their parent partition:

``` r

# We cut the tree twice, to have a parent partition and a child partition
tree <- cut_tree(tree, cut_height = c(.5, .6))
plot(tree, plot_type = "bars",
     partitions = 2, # Bars to be drawn
     color_partition = 1) # Colors to apply
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-23-1.png) Here we
see which clusters belong to the orange cluster in the parent partition,
and which clusters belong to the purple cluster.

Now, we can do it the other way around: how the parent clusters are
split into children clusters.

``` r

plot(tree, plot_type = "bars", partitions = 1, color_partition = 2)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-24-1.png)

This allows you to represent more layers of information in a single
graph.

### 4.2 Hierarchical / tree-like clusters

Tree-like clustering results can be visualized with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). If the tree
was cut, clustering results will be displayed as colored bars below the
tree.

``` r

tree <- hclu_hierarclust(vegedissim)
tree <- cut_tree(tree,
                 cut_height = .6)

plot(tree)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-25-1.png)

On large trees, labels are difficult to read, so we can turn them off:

``` r

plot(tree, labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-26-1.png)

If there are multiple cut heights, multiple bars display from coarsest
to finest:

``` r

tree <- cut_tree(tree,
                 cut_height = c(.4, .5, .6))

plot(tree, labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-27-1.png)

The tree can be flipped (i.e., mirrored) if need be:

``` r

plot(tree, flip = TRUE,
     labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-28-1.png)

To read complex labels, it is possible to plot the tree horizontally in
long page format (note that you can also flip the horizontal tree if you
want):

``` r

plot(tree, horizontal = TRUE,
     cex = .4) # Reduce label font size
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-29-1.png)

The clusters of finer partitions can be collapsed into triangles, as
wide as their number of sites and as high as the height at which their
sites are joined. The argument `partitions` chooses the partitions to
show, by their number (1 for the first partition, 2 for the second,
etc.) or their name; the clusters of the last one are collapsed:

``` r

# No partitions: collapse only the finest one
plot(tree, collapse = TRUE, labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-30-1.png)

``` r

# Specify the partitions to display with arg partitions
# e.g. here, only partitions 1 and 2
plot(tree, partitions = 1:2, collapse = TRUE, labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-30-2.png)

The tree can also be drawn down to the clusters only, with the size of
each cluster at the tip of its branch (one square per site). The height
axis then starts just below the lowest cut, to spread the merges between
clusters:

``` r

plot(tree, plot_type = "cluster_tree", horizontal = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-31-1.png)

Again, you can choose which partition to display at the finest level
with arg `partitions`

``` r

plot(tree, partitions = 2, plot_type = "cluster_tree", horizontal = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-32-1.png)

When the merges between clusters are crowded at the top of the tree,
they can be evenly spaced, keeping their order (`heights = "rank"`), or
drawn with one step per partition (`heights = "levels"`):

``` r

plot(tree, plot_type = "cluster_tree", heights = "rank", horizontal = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-33-1.png)

*Notice that when heights correspond to rank, the axis disappears,
because it loses the ‘dissimilarity’ meaning*.

The tree can displayed many layers of cuts, which becomes fascinating as
it shows how the hierarchy evolves as new cuts are drawn:

``` r

tree <- cut_tree(tree,
                 cut_height = seq(0, .6, length = 30))

plot(tree, labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-34-1.png)

Note, however, that only 12 colors can be displayed simultaneously, and
so as the number of clusters increases, there are more and more clusters
displayed in shades of grey, which makes them more difficult to read
past a certain point.

You can predefine your colors before plotting the tree with
[`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md):

``` r

tree <- cut_tree(tree,
                 cut_height = .5)
# Choose a colorblind-friendly palette
tree <- bioregion_colors(tree, palette = "Safe")
plot(tree, horizontal = TRUE,
     labels = FALSE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-35-1.png)

#### How to plot Infomap (Map Equation) and other non-standard hierarchical trees

Some algorithms produce nested partitions without building a dendrogram.
For example,
[`netclu_infomap()`](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)
with `show_hierarchy = TRUE` returns several levels of clusters.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws these
levels as a tree, with one step per level. The main difference with
previous trees is that infomap is not a “height-based” tree. Branch
lengths have no meaning in infomap trees. They rather represent nested
partitions, called *levels*:

- the first level divides the sites into a few large bioregions
- the second level divides each of these bioregions into smaller
  sub-bioregions
- and so on, until the algorithm finds no further meaningful subdivision

The tree simply shows how the clusters of one level fit inside the
clusters of the level above, like countries divided into regions, and
regions into departments. Consequently, it is NOT a tree you can cut,
and, contrary to standard hierarchical trees, a cluster can split into
more than two sub-clusters. Furthermore, some clusters are not
subdivided further: their branches stop earlier than others. In the
output, these clusters keep the same name at all deeper levels.

``` r

install_binaries(verbose = FALSE)
data(fishdf)

fish_infomap <- netclu_infomap(fishdf, bipartite_version = TRUE,
                               markovtime = .7, # we specify a lower markov time
                               # than the default here, to get at least 2 levels
                               show_hierarchy = TRUE)

# The default plot is a cluster tree
plot(fish_infomap, horizontal = TRUE, cex = 1)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-36-1.png)

``` r

# You can ask for a simple tree with plot_type = "tree"
plot(fish_infomap, plot_type = "tree", 
     horizontal = TRUE, cex = .3)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-36-2.png)

``` r

# You can also collapse the tree leaves to make it more readable
plot(fish_infomap, plot_type = "tree", 
     horizontal = TRUE, cex = .3,
     collapse = TRUE)
```

![](a5_1_visualization_files/figure-html/unnamed-chunk-37-1.png)

And that’s all for now for cluster plots, but keep checking this file
regularly, as we might come up with more plot ideas!

## 5. Exporting networks for visualization in Gephi

In addition to visualizing bioregions on maps, you may want to visualize
the network structure of your data using specialized network
visualization software like [Gephi](https://gephi.org/). The `exportGDF`
function allows you to export your data in GDF (Graph Data Format),
which can be directly imported into Gephi.

For these examples, we will use the fish dataset which provides a
clearer network structure that is easier to visualize:

``` r

# Load the fish dataset
data(fishdf)
data(fishmat)
```

### 5.1 Basic network export

Let’s start by creating a simple bipartite network from the fish data
and exporting it without weights or colors:

``` r

# Export the entire bipartite (site-species) fish network to GDF format 
exportGDF(fishdf,
          col1 = "Site",
          col2 = "Species",
          file = "fish_network_basic.gdf")
```

This creates a GDF file with nodes (sites and species) and edges
(interactions) that can be opened in Gephi. The fish dataset contains
2703 site-species interactions across 338 sites and 195 species.

In gephi, this network can be visualised as follows, after applying a
layout algorithm:

![Network fish basic](figures/05_network_fish_basic.png)

### 5.2 Unipartite network with edge weights

In addition to bipartite networks, you can create unipartite networks
where sites are connected based on shared species composition. Edge
weights represent jaccard similarity between sites:

``` r

# Create a site-to-site network based on shared species
# First, convert the matrix to a dissimilarity object
fish_dissim_jac <- dissimilarity(fishmat, metric = "Jaccard")

# Convert dissimilarity to similarity (shared species proportion)
# Jaccard similarity = 1 - Jaccard dissimilarity
# This is because links in a network represent similarity, not dissimilarity
fish_sim <- dissimilarity_to_similarity(fish_dissim_jac)

# Remove all links with weight 0 BEFORE exporting to GDF
fish_sim <- fish_sim[which(fish_sim$Jaccard > 0), ]  # Remove zero similarity links

# Export the unipartite network with shared species as weights
exportGDF(fish_sim,
          col1 = "Site1",
          col2 = "Site2",
          weight = "Jaccard",
          file = "fish_network_unipartite_weighted.gdf")
```

![Network fish unipartite
weighted](figures/05_network_fish_unipartite_weighted.png)

This creates a unipartite network where nodes are sites and edges
represent the Jaccard similarity index. The edge weight indicates how
similar is each each pair of sites. In Gephi, these weights can be used
to size edges or influence layout algorithms, helping to visualize which
sites have similar fish communities.

### 5.3 Network with node attributes and colors

To create more sophisticated visualizations, you can add node attributes
and colors. For example, we can color sites according to their bioregion
assignment:

``` r

# Install binaries if not already installed
install_binaries(binpath = "tempdir")

# Run a bipartite network clustering analysis
fish_netclu <- netclu_infomap(fishdf, # Run Infomap on the bipartite network
                             weight = FALSE, # Do not use weights (no abundance here) for clustering
                             binpath = "tempdir", # Path to the Infomap binary
                             numtrials = 100, # Number of trials for Infomap
                             show_hierarchy = TRUE, # Show hierarchical levels
                             bipartite = TRUE, # This is a bipartite network
                             bipartite_version = FALSE, # Use standard Infomap, not bipartite version
                             # Because in our experience the standard version produces better results even
                             # for bipartite biogeographical networks
                             site_col = "Site", # Column name for sites
                             species_col = "Species", # Column name for species
                             seed = 123) # Set seed for reproducibility

fish_netclu <- bioregion_colors(fish_netclu, # Assign colors to bioregions
                                   palette = "Vivid", # Color palette from rcartocolor
                                   cluster_ordering = "n_both") # Order clusters by total size (sites + species)

# Let's make the map for the first bioregionalization
map_bioregions(fish_netclu, map = fishsf,
               partition_index = 1)

# Export the first bioregionalization to GDF with colors and attributes
exportGDF(fishdf, 
          col1 = "Site", 
          col2 = "Species",
          bioregions = fish_netclu,
          file = "fish_network_colored.gdf")
# Note: the first bioregionalization is exported by default, but you 
# can specify another one with the 'cluster_column' argument
```

![Network fish bipartite
infomap](figures/05_map_fish_bipartite_infomap.png)

![Network fish bipartite
colored](figures/05_network_fish_bipartite_colored.png)

### 5.4 Visualizing in Gephi

Once you have exported your network to GDF format, you can visualize it
in Gephi:

1.  Open Gephi and go to **File \> Open** and select your GDF file
2.  The network will be imported with:
    - All nodes (sites and species) and edges (occurrences here)
    - Edge weights representing abundance (if specified)
    - Node colors indicating bioregions for sites (if specified)
    - Node attributes such as bioregion and node type (if specified)
3.  **Layout**: Apply a ForceAtlas 2 algorithm to arrange nodes
4.  **Appearance**:
    - Node colors are already set based on bioregions
    - Adjust node sizes based on degree (= occurrence/richness in our
      example here) or other attributes
5.  Now, spend time on Gephi to make it look gorgeous ;)
