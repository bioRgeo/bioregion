# Convert a matrix or list of matrices to a bioregion (dis)similarity object

Converts a (dis)similarity `matrix` or a `list` of such matrices into a
`bioregion.pairwise` object compatible with the `bioregion` package. The
input can come from base R, `dist` objects, or outputs from other
packages.

## Usage

``` r
as_bioregion_pairwise(
  mat,
  metric_name = NULL,
  pkg = NULL,
  is_similarity = FALSE
)
```

## Arguments

- mat:

  A `matrix`, a `dist` object, or a `list` of these representing
  pairwise similarity or dissimilarity values to convert into a
  `bioregion.pairwise` object. This function can also directly handle
  outputs from other R packages (see the `pkg` argument).

- metric_name:

  Optional `character` vector or single `character` string specifying
  the name of the (dis)similarity metric(s), which will appear as column
  names in the output (see Note).

- pkg:

  An optional `character` string indicating the name of the package from
  which `mat` was generated (`NULL` by default, see Details). Available
  options are `"adespatial"`, `"betapart"`, `"ecodist"`, or `"vegan"`.

- is_similarity:

  A `logical` value indicating whether the input data represents
  similarity (`TRUE`) or dissimilarity (`FALSE`).

## Value

A dissimilarity or similarity object of class `bioregion.pairwise`,
compatible with the `bioregion` package.

## Details

This function can directly handle outputs from ten functions across four
packages:

- **adespatial**:
  [beta.div](http://adeverse.github.io/adespatial/reference/beta.div.md),
  [beta.div.comp](http://adeverse.github.io/adespatial/reference/beta.div.comp.md)

- **betapart**:
  [beta.pair](https://rdrr.io/pkg/betapart/man/beta.pair.html),
  [beta.pair.abund](https://rdrr.io/pkg/betapart/man/beta.pair.abund.html),
  [betapart.core](https://rdrr.io/pkg/betapart/man/betapart.core.html),
  [betapart.core.abund](https://rdrr.io/pkg/betapart/man/betapart.core.abund.html)

- **ecodist**:
  [distance](https://rdrr.io/pkg/ecodist/man/distance.html),
  [bcdist](https://rdrr.io/pkg/ecodist/man/bcdist.html)

- **vegan**:
  [vegdist](https://vegandevs.github.io/vegan/reference/vegdist.html),
  [designdist](https://vegandevs.github.io/vegan/reference/designdist.html)

See the documentation of these packages for more information:

- https://cran.r-project.org/package=adespatial

- https://cran.r-project.org/package=betapart

- https://cran.r-project.org/package=ecodist

- https://cran.r-project.org/package=vegan

## Note

If no specific package is specified (i.e., `pkg = NULL`), site names
will be based on the row names of the first matrix. If row names are
`NULL`, they will be generated automatically. If `mat` is a named list,
those names will be used as column names only if `metric_name = NULL`.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html>.

Associated functions:
[dissimilarity](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
[similarity](https://bioRgeo.github.io/bioregion/reference/similarity.md)
[bind_pairwise](https://bioRgeo.github.io/bioregion/reference/bind_pairwise.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
mat <- matrix(runif(100), 10, 10)
rownames(mat) <- paste0("s",1:10)

pair <- as_bioregion_pairwise(list(mat,mat,mat), 
                              metric_name = NULL,
                              pkg = NULL,
                              is_similarity = FALSE)
                              
pair
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  10 
#>  - Total number of species:  NA 
#>  - Number of rows:  45 
#>  - Number of dissimilarity metrics:  3 
#> 
#> 
#>    Site1 Site2    Metric1    Metric2    Metric3
#> 2     s1    s2 0.03424133 0.03424133 0.03424133
#> 3     s1    s3 0.73531960 0.73531960 0.73531960
#> 4     s1    s4 0.30083081 0.30083081 0.30083081
#> 5     s1    s5 0.64167935 0.64167935 0.64167935
#> 6     s1    s6 0.17467589 0.17467589 0.17467589
#> 7     s1    s7 0.57004495 0.57004495 0.57004495
#> 8     s1    s8 0.02806097 0.02806097 0.02806097
#> 9     s1    s9 0.26137136 0.26137136 0.26137136
#> 10    s1   s10 0.51855664 0.51855664 0.51855664
#> 13    s2    s3 0.19595673 0.19595673 0.19595673
#> 14    s2    s4 0.63646561 0.63646561 0.63646561
#> 15    s2    s5 0.66028435 0.66028435 0.66028435
#> 16    s2    s6 0.53157354 0.53157354 0.53157354
#> 17    s2    s7 0.33571908 0.33571908 0.33571908
#> 18    s2    s8 0.46598719 0.46598719 0.46598719
#> 19    s2    s9 0.29005016 0.29005016 0.29005016
#> 20    s2   s10 0.84612005 0.84612005 0.84612005
#> 24    s3    s4 0.47902455 0.47902455 0.47902455
#> 25    s3    s5 0.09602416 0.09602416 0.09602416
#> 26    s3    s6 0.49363702 0.49363702 0.49363702
#> 27    s3    s7 0.59626279 0.59626279 0.59626279
#> 28    s3    s8 0.39003139 0.39003139 0.39003139
#> 29    s3    s9 0.48007517 0.48007517 0.48007517
#> 30    s3   s10 0.71826972 0.71826972 0.71826972
#> 35    s4    s5 0.76560016 0.76560016 0.76560016
#> 36    s4    s6 0.77930863 0.77930863 0.77930863
#> 37    s4    s7 0.19151803 0.19151803 0.19151803
#> 38    s4    s8 0.02006522 0.02006522 0.02006522
#> 39    s4    s9 0.92000555 0.92000555 0.92000555
#> 40    s4   s10 0.24131402 0.24131402 0.24131402
#> 46    s5    s6 0.20417834 0.20417834 0.20417834
#> 47    s5    s7 0.94776394 0.94776394 0.94776394
#> 48    s5    s8 0.37697093 0.37697093 0.37697093
#> 49    s5    s9 0.40072018 0.40072018 0.40072018
#> 50    s5   s10 0.54704337 0.54704337 0.54704337
#> 57    s6    s7 0.54248041 0.54248041 0.54248041
#> 58    s6    s8 0.55991284 0.55991284 0.55991284
#> 59    s6    s9 0.21317271 0.21317271 0.21317271
#> 60    s6   s10 0.83480182 0.83480182 0.83480182
#> 68    s7    s8 0.85708359 0.85708359 0.85708359
#> 69    s7    s9 0.67176682 0.67176682 0.67176682
#> 70    s7   s10 0.02795603 0.02795603 0.02795603
#> 79    s8    s9 0.05861411 0.05861411 0.05861411
#> 80    s8   s10 0.46938430 0.46938430 0.46938430
#> 90    s9   s10 0.80568003 0.80568003 0.80568003
```
