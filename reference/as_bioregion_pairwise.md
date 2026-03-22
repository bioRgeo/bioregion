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
                              
```
