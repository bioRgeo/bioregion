# Create a map of bioregions

This plot function can be used to visualize bioregions based on a
`bioregion.clusters` object combined with a spatial object (`sf` or
`terra`).

## Usage

``` r
map_bioregions(
  bioregionalization,
  map,
  partition_index = NULL,
  map_as_output = FALSE,
  plot = TRUE,
  clusters = NULL,
  geometry = NULL,
  write_clusters = NULL,
  ...
)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- map:

  A spatial object that can be handled by `sf` or `terra`. The first
  attribute or layer should correspond to the sites' ID (see Details).

- partition_index:

  An `integer`, `character`, or `NULL` specifying which
  `bioregionalization`'s partition(s) to plot. By default (`NULL`), all
  partitions are plotted. If an `integer` or vector of `integers` is
  provided, partition(s) are selected by column number(s) in the
  `bioregionalization` data.frame (starting from 1 after the ID column).
  If a `character` or vector of `characters`, partition(s) are selected
  by name(s) matching column names in `bioregionalization`.

- map_as_output:

  A `boolean` indicating if the `sf` `data.frame` object used for the
  plot should be returned.

- plot:

  A `boolean` indicating if the plot should be drawn.

- clusters:

  Deprecated. Use `bioregionalization` instead. The former
  `bioregionalization` has been replaced by `partition_index`.

- geometry:

  Deprecated. Use `map` instead.

- write_clusters:

  Deprecated. Use `map_as_output` instead.

- ...:

  Further arguments to be passed to
  [`sf::plot()`](https://r-spatial.github.io/sf/reference/plot.html).

## Value

One or several maps of bioregions if `plot = TRUE` and the `sf`
`data.frame` object used for the plot if `map_as_output = TRUE`.

## Details

The site IDs in `bioregionalization` and `map` should correspond. They
must have the same type (i.e., `character` if `bioregionalization` is a
`bioregion.clusters` object), and the sites in `bioregionalization`
should be included among the sites in `map`. If `map` is an `sf` or a
`SpatVector` (`terra`) object, it should contain an attribute table with
the IDs in the first column. If `map` is a `SpatRaster` (`terra`)
object, it should contain the IDs in the first layer.

If the `bioregionalization` object contains both types of nodes (sites
and species), only site will be mapped. The function automatically
filters to site nodes using the `node_type` attribute.

**Colors**: If the `bioregionalization` object contains colors (added
via
[`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md)),
these colors will be automatically used for plotting. Otherwise, the
default `sf` color scheme will be applied.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a5_1_visualization.html>.

Associated functions:
[bioregion_colors](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
data(fishmat)
data(fishsf)

net <- similarity(fishmat, metric = "Simpson")
clu <- netclu_greedy(net)
mapclu <- map_bioregions(clu, 
                         map = fishsf, 
                         map_as_output = TRUE, 
                         plot = FALSE)

# With colors
clu_colored <- bioregion_colors(clu)
mapclu <- map_bioregions(clu_colored, 
                         map = fishsf, 
                         map_as_output = TRUE, 
                         plot = FALSE)
           
```
