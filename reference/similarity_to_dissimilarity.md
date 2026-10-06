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
#>    Site1 Site2 Jaccard Jaccardturn  Sorensen   Simpson      Bray   Brayturn
#> 2     s1    s2     1.0         1.0 1.0000000 1.0000000 0.1204139 0.12090680
#> 3     s1    s3     1.0         1.0 1.0000000 1.0000000 0.1818182 0.19080302
#> 4     s1    s4     0.9         1.0 0.9473684 1.0000000 0.1286307 0.19159456
#> 5     s1    s5     0.9         1.0 0.9473684 1.0000000 0.5017964 0.52342286
#> 8     s2    s3     1.0         1.0 1.0000000 1.0000000 0.1372742 0.14344544
#> 9     s2    s4     0.9         1.0 0.9473684 1.0000000 0.1268252 0.18788628
#> 10    s2    s5     0.9         1.0 0.9473684 1.0000000 0.0769462 0.08060453
#> 14    s3    s4     0.9         1.0 0.9473684 1.0000000 0.3671668 0.51421508
#> 15    s3    s5     0.9         1.0 0.9473684 1.0000000 0.2521902 0.27659574
#> 20    s4    s5     0.8         0.8 0.8888889 0.8888889 0.1357928 0.21384425
#>       Euclidean  a b c   A    B    C
#> 2  0.0007855370 10 0 0 192 1409 1396
#> 3  0.0007790752 10 0 0 278 1323 1179
#> 4  0.0008803537  9 1 0 155 1446  654
#> 5  0.0012199128  9 1 0 838  763  901
#> 8  0.0008423520 10 0 0 209 1379 1248
#> 9  0.0010275651  9 1 0 152 1436  657
#> 10 0.0008753303  9 1 0 128 1460 1611
#> 14 0.0011580265  9 1 0 416 1041  393
#> 15 0.0008849734  9 1 0 403 1054 1336
#> 20 0.0010725430  8 1 1 173  636 1566

dissimilarity <- similarity_to_dissimilarity(simil)
dissimilarity
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2 Jaccard Jaccardturn   Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2     0.0         0.0 0.00000000 0.0000000 0.8795861 0.8790932
#> 3     s1    s3     0.0         0.0 0.00000000 0.0000000 0.8181818 0.8091970
#> 4     s1    s4     0.1         0.0 0.05263158 0.0000000 0.8713693 0.8084054
#> 5     s1    s5     0.1         0.0 0.05263158 0.0000000 0.4982036 0.4765771
#> 8     s2    s3     0.0         0.0 0.00000000 0.0000000 0.8627258 0.8565546
#> 9     s2    s4     0.1         0.0 0.05263158 0.0000000 0.8731748 0.8121137
#> 10    s2    s5     0.1         0.0 0.05263158 0.0000000 0.9230538 0.9193955
#> 14    s3    s4     0.1         0.0 0.05263158 0.0000000 0.6328332 0.4857849
#> 15    s3    s5     0.1         0.0 0.05263158 0.0000000 0.7478098 0.7234043
#> 20    s4    s5     0.2         0.2 0.11111111 0.1111111 0.8642072 0.7861557
#>    Euclidean  a b c   A    B    C
#> 2  1272.0145 10 0 0 192 1409 1396
#> 3  1282.5732 10 0 0 278 1323 1179
#> 4  1134.9070  9 1 0 155 1446  654
#> 5   818.7307  9 1 0 838  763  901
#> 8  1186.1522 10 0 0 209 1379 1248
#> 9   972.1744  9 1 0 152 1436  657
#> 10 1141.4259  9 1 0 128 1460 1611
#> 14  862.5381  9 1 0 416 1041  393
#> 15 1128.9774  9 1 0 403 1054 1336
#> 20  931.3635  8 1 1 173  636 1566
```
