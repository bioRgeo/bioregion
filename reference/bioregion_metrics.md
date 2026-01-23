# Calculate contribution metrics for bioregions

This function calculates the number of sites per bioregion, as well as
the number of species these sites have, the number of endemic species,
and the proportion of endemism.

## Usage

``` r
bioregion_metrics(bioregionalization, comat, map = NULL, col_bioregion = NULL)
```

## Arguments

- bioregionalization:

  A `bioregion.clusters` object.

- comat:

  A co-occurrence `matrix` with sites as rows and species as columns.

- map:

  A spatial `sf data.frame` with sites and bioregions. It is the output
  of the function `map_bioregions`. `NULL` by default.

- col_bioregion:

  An `integer` specifying the column position of the bioregion.

## Value

A `data.frame` with 5 columns, or 6 if spatial coherence is computed.

## Details

Endemic species are species found only in the sites belonging to one
bioregion.

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
#>   Bioregion Site_number Species_number Endemics Percentage_Endemic
#> 1         1           5             10       10                100
```
