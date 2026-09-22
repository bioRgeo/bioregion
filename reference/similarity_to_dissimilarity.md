# Convert similarity metrics to dissimilarity metrics

This function converts a `data.frame` of similarity metrics between
sites into dissimilarity metrics (beta diversity).

## Usage

``` r
similarity_to_dissimilarity(similarity, include_formula = TRUE)
```

## Arguments

- similarity:

  The output object from
  [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  or
  [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md).

- include_formula:

  A `boolean` indicating whether metrics based on custom formula(s)
  should also be converted (see Details). The default is `TRUE`.

## Value

A `data.frame` with additional class `bioregion.pairwise`, providing
dissimilarity metric(s) between each pair of sites based on a similarity
object.

## Note

The behavior of this function changes depending on column names. Columns
`Site1` and `Site2` are copied identically. If there are columns called
`a`, `b`, `c`, `A`, `B`, `C` they will also be copied identically. If
there are columns based on your own formula (argument `formula` in
[`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md))
or not in the original list of similarity metrics (argument `metrics` in
[`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md))
and if the argument `include_formula` is set to `FALSE`, they will also
be copied identically. Otherwise there are going to be converted like
they other columns (default behavior).

If a column is called `Euclidean`, its distance will be calculated based
on the following formula:

Euclidean distance = (1 - Euclidean similarity) / Euclidean similarity

Otherwise, all other columns will be transformed into dissimilarity with
the following formula:

dissimilarity = 1 - similarity

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html>.

Associated functions:
[dissimilarity](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
similarity_to_dissimilarity

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

simil <- similarity(comat, metric = "all")
simil
#> Data.frame of similarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of similarity metrics:  7 
#> 
#> 
#>    Site1 Site2 Jaccard Jaccardturn  Sorensen Simpson       Bray   Brayturn
#> 2     s1    s2     0.7   0.7777778 0.8235294   0.875 0.29291154 0.38580247
#> 3     s1    s3     0.8   1.0000000 0.8888889   1.000 0.32753623 0.69753086
#> 4     s1    s4     0.6   0.6000000 0.7500000   0.750 0.08715084 0.12037037
#> 5     s1    s5     0.7   0.7777778 0.8235294   0.875 0.19815668 0.46450617
#> 8     s2    s3     0.9   1.0000000 0.9473684   1.000 0.40113529 0.60056657
#> 9     s2    s4     0.7   0.7777778 0.8235294   0.875 0.08995911 0.09348442
#> 10    s2    s5     1.0   1.0000000 1.0000000   1.000 0.38503914 0.62700661
#> 14    s3    s4     0.8   1.0000000 0.8888889   1.000 0.34357714 0.48949212
#> 15    s3    s5     0.9   1.0000000 0.9473684   1.000 0.14615726 0.15577652
#> 20    s4    s5     0.7   0.7777778 0.8235294   0.875 0.20328426 0.31436077
#>       Euclidean a b c   A    B    C
#> 2  0.0017127162 7 1 2 250  398  809
#> 3  0.0009888118 8 0 2 452  196 1660
#> 4  0.0014871408 6 2 2  78  570 1064
#> 5  0.0009165988 7 1 2 301  347 2089
#> 8  0.0011149902 9 0 1 636  423 1476
#> 9  0.0012273278 7 2 1  99  960 1043
#> 10 0.0010604161 9 0 0 664  395 1726
#> 14 0.0010485520 8 2 0 559 1553  583
#> 15 0.0006832077 9 1 0 329 1783 2061
#> 20 0.0008718495 7 1 2 359  783 2031

dissimilarity <- similarity_to_dissimilarity(simil)
dissimilarity
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2 Jaccard Jaccardturn   Sorensen Simpson      Bray  Brayturn
#> 2     s1    s2     0.3   0.2222222 0.17647059   0.125 0.7070885 0.6141975
#> 3     s1    s3     0.2   0.0000000 0.11111111   0.000 0.6724638 0.3024691
#> 4     s1    s4     0.4   0.4000000 0.25000000   0.250 0.9128492 0.8796296
#> 5     s1    s5     0.3   0.2222222 0.17647059   0.125 0.8018433 0.5354938
#> 8     s2    s3     0.1   0.0000000 0.05263158   0.000 0.5988647 0.3994334
#> 9     s2    s4     0.3   0.2222222 0.17647059   0.125 0.9100409 0.9065156
#> 10    s2    s5     0.0   0.0000000 0.00000000   0.000 0.6149609 0.3729934
#> 14    s3    s4     0.2   0.0000000 0.11111111   0.000 0.6564229 0.5105079
#> 15    s3    s5     0.1   0.0000000 0.05263158   0.000 0.8538427 0.8442235
#> 20    s4    s5     0.3   0.2222222 0.17647059   0.125 0.7967157 0.6856392
#>    Euclidean a b c   A    B    C
#> 2   582.8679 7 1 2 250  398  809
#> 3  1010.3148 8 0 2 452  196 1660
#> 4   671.4313 6 2 2  78  570 1064
#> 5  1089.9899 7 1 2 301  347 2089
#> 8   895.8689 9 0 1 636  423 1476
#> 9   813.7782 7 2 1  99  960 1043
#> 10  942.0260 9 0 0 664  395 1726
#> 14  952.6962 8 2 0 559 1553  583
#> 15 1462.6838 9 1 0 329 1783 2061
#> 20 1145.9869 7 1 2 359  783 2031
```
