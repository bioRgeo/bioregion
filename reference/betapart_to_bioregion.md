# Convert betapart dissimilarity to bioregion dissimilarity (DEPRECATED)

This function converts dissimilarity results produced by the betapart
package (and packages using betapart, such as phyloregion) into a
dissimilarity object compatible with the bioregion package. This
function only converts object types to make them compatible with
bioregion; it does not modify the beta-diversity values. This function
allows the inclusion of phylogenetic beta diversity to compute
bioregions with bioregion.

## Usage

``` r
betapart_to_bioregion(betapart_result)
```

## Arguments

- betapart_result:

  An object produced by the betapart package (e.g., using the
  `beta.pair` function).

## Value

A dissimilarity object of class `bioregion.pairwise`, compatible with
the bioregion package.

## See also

This function is deprecated, use
[as_bioregion_pairwise](https://bioRgeo.github.io/bioregion/reference/as_bioregion_pairwise.md)
instead.

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

if (FALSE) { # \dontrun{
beta_div <- betapart::beta.pair.abund(comat)
betapart_to_bioregion(beta_div)
} # }
```
