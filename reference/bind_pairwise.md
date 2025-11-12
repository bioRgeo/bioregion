# Combine and enrich bioregion (dis)similarity object(s)

Combine two `bioregion.pairwise` objects and/or compute new pairwise
metrics based on the columns of the object(s).

## Usage

``` r
bind_pairwise(primary_metrics, secondary_metrics, new_metrics = NULL)
```

## Arguments

- primary_metrics:

  A `bioregion.pairwise` object. This is the main set of pairwise
  metrics that will be used as a base for the combination.

- secondary_metrics:

  A second `bioregion.pairwise` object to be combined with
  `primary_metrics`. It must have the same sites identifiers and
  pairwise structure. Can be set to `NULL` if `new_metrics` is
  specified.

- new_metrics:

  A `character` vector or a single `character` string specifying custom
  formula(s) based on the column names of `primary_metrics` and
  `secondary_metrics` (see Details).

## Value

A new `bioregion.pairwise` object containing the combined and/or
enriched data. It includes all original metrics from the inputs, as well
as any newly computed metrics.

## Details

When both `primary_metrics` and `secondary_metrics` are provided and if
the pairwise structure is identical the function combine the two
objects. If `new_metrics` is provided, each formula is evaluated based
on the column names of `primary_metrics` (and `secondary_metrics` if
provided).

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html>.

Associated functions:
[dissimilarity](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
[similarity](https://bioRgeo.github.io/bioregion/reference/similarity.md)
[as_bioregion_pairwise](https://bioRgeo.github.io/bioregion/reference/as_bioregion_pairwise.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("s", 1:5)
colnames(comat) <- paste0("sp", 1:10)

sim <- bind_pairwise(primary_metrics = similarity(comat, 
                                                               metric = "abc"),
                                  secondary_metrics = similarity(comat, 
                                                                 metric = "Simpson"),
                                  new_metrics = "1 - (b + c) / (a + b + c)")
```
