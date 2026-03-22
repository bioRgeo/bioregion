# Calculate metrics for bioregions

This function calculates the number of sites, the number of species, the
number of endemic species and the proportion of endemism per bioregion.
The spatial coherence can be optionally computed if a spatial object is
provided.

## Usage

``` r
bioregion_metrics(bioregionalization, comat, map = NULL, col_bioregion = NULL)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- comat:

  A site-species `matrix` with sites as rows and species as columns.

- map:

  A spatial object that can be handled by `sf` or `terra`. The first
  attribute or layer should correspond to the sites' ID (see Details).
  Needed only for the spatial coherence (`NULL` by default).

- col_bioregion:

  Deprecated.

## Value

A `data.frame` with 5 columns (Bioregion ID and metrics described below)
or 7 if spatial coherence is computed.

- **NbSites**: Number of sites per bioregion

- **Richness**: Number of distinct species per bioregion.

- **Rich_Endemics**: Number of species found only in the bioregion.

- **Prop_Endemics**: Fraction of endemics species.

- [**SC_size**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregion-metrics-spatial-coherence):
  Spatial coherence based on size, fraction of the number of site
  contained in the bioregion's largest contiguous patch.

- [**SC_area**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregion-metrics-spatial-coherence):
  Spatial coherence based on area, fraction of the bioregion area
  contained in its largest contiguous patch.

Note that if `bioregionalization` contains multiple partitions (i.e., if
`dim(bioregionalization$clusters) > 2`), a list will be returned.

## Details

`map` should be the output of
`map_bioregions(bioregionalization, geometry, write_clusters = TRUE)`

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html>.

Associated functions:
[site_species_metrics](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md)
[bioregionalization_metrics](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)

## Author

Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)

## Examples

``` r
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

net <- similarity(comat, metric = "Simpson")
clust <- netclu_louvain(net)

bioregion_metrics(bioregionalization = clust, 
                  comat = comat) 
#>   Bioregion NbSites Richness Rich_Endemics Prop_Endemics
#> 1         1       5       10            10             1
```
