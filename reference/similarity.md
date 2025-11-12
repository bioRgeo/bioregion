# Compute similarity metrics between sites based on species composition

This function generates a `data.frame` where each row provides one or
several similarity metrics between pairs of sites, based on a
co-occurrence `matrix` with sites as rows and species as columns.

## Usage

``` r
similarity(comat, metric = "Simpson", formula = NULL, method = "prodmat")
```

## Arguments

- comat:

  A co-occurrence `matrix` with sites as rows and species as columns.

- metric:

  A `character` vector or a single `character` string specifying the
  metrics to compute (see Details). Available options are `"abc"`,
  `"ABC"`, `"Jaccard"`, `"Jaccardturn"`, `"Sorensen"`, `"Simpson"`,
  `"Bray"`, `"Brayturn"`, and `"Euclidean"`. If `"all"` is specified,
  all metrics will be calculated. Can be set to `NULL` if `formula` is
  used.

- formula:

  A `character` vector or a single `character` string specifying custom
  formula(s) based on the `a`, `b`, `c`, `A`, `B`, and `C` quantities
  (see Details). The default is `NULL`.

- method:

  A `character` string specifying the method to compute `abc` (see
  Details). The default is `"prodmat"`, which is more efficient but
  memory-intensive. Alternatively, `"loops"` is less memory-intensive
  but slower.

## Value

A `data.frame` with the additional class `bioregion.pairwise`,
containing one or several similarity metrics between pairs of sites. The
first two columns represent the pairs of sites. There is one column per
similarity metric provided in `metric` and `formula`, except for the
`abc` and `ABC` metrics, which are stored in three separate columns (one
for each letter).

## Details

With `a` the number of species shared by a pair of sites, `b` species
only present in the first site and `c` species only present in the
second site.

Jaccard = 1 - (b + c) / (a + b + c)

Jaccardturn = 1 - 2min(b, c) / (a + 2min(b, c)) (Baselga, 2012)

Sorensen = 1 - (b + c) / (2a + b + c)

Simpson = 1 - min(b, c) / (a + min(b, c))

If abundances data are available, Bray-Curtis and its turnover component
can also be computed with the following equation:

Bray = 1 - (B + C) / (2A + B + C)

Brayturn = 1 - min(B, C) / (A + min(B, C)) (Baselga, 2013)

with `A` the sum of the lesser values for common species shared by a
pair of sites. `B` and `C` are the total number of specimens counted at
both sites minus `A`.

`formula` can be used to compute customized metrics with the terms `a`,
`b`, `c`, `A`, `B`, and `C`. For example
`formula = c("1 - pmin(b,c) / (a + pmin(b,c))", "1 - (B + C) / (2*A + B + C)")`
will compute the Simpson and Bray-Curtis similarity metrics,
respectively. Note that `pmin` is used in the Simpson formula because
`a`, `b`, `c`, `A`, `B` and `C` are `numeric` vectors.

Euclidean computes the Euclidean similarity between each pair of sites
following this equation:

Euclidean = 1 / (1 + d_ij)

Where d_ij is the Euclidean distance between site i and site j in terms
of species composition.

## References

Baselga A (2012) The Relationship between Species Replacement,
Dissimilarity Derived from Nestedness, and Nestedness. *Global Ecology
and Biogeography* 21, 1223–1232.

Baselga A (2013) Separating the two components of abundance-based
dissimilarity: balanced changes in abundance vs. abundance gradients.
*Methods in Ecology and Evolution* 4, 552–557.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html>.

Associated functions:
[dissimilarity](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
[similarity_to_dissimilarity](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("s", 1:5)
colnames(comat) <- paste0("sp", 1:10)

sim <- similarity(comat, metric = c("abc", "ABC", "Simpson", "Brayturn"))

sim <- similarity(comat, metric = "all",
formula = "1 - (b + c) / (a + b + c)")
```
