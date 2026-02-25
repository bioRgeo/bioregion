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
#>    Site1 Site2 Jaccard Jaccardturn  Sorensen   Simpson       Bray  Brayturn
#> 2     s1    s2     0.9        1.00 0.9473684 1.0000000 0.25078370 0.5405405
#> 3     s1    s3     0.7        1.00 0.8235294 1.0000000 0.27382550 0.6244898
#> 4     s1    s4     1.0        1.00 1.0000000 1.0000000 0.47066246 0.7612245
#> 5     s1    s5     0.9        1.00 0.9473684 1.0000000 0.14204004 0.3040816
#> 8     s2    s3     0.6        0.75 0.7500000 0.8571429 0.04437401 0.2837838
#> 9     s2    s4     0.9        1.00 0.9473684 1.0000000 0.11906677 0.5000000
#> 10    s2    s5     0.8        0.80 0.8888889 0.8888889 0.08769932 0.5202703
#> 14    s3    s4     0.7        1.00 0.8235294 1.0000000 0.56830986 0.7369863
#> 15    s3    s5     0.6        0.75 0.7500000 0.8571429 0.08827915 0.0920398
#> 20    s4    s5     0.9        1.00 0.9473684 1.0000000 0.27376989 0.3378995
#>       Euclidean  a b c   A    B    C
#> 2  0.0037812677  9 1 0  80  410   68
#> 3  0.0011706950  7 3 0 306  184 1439
#> 4  0.0025378580 10 0 0 373  117  722
#> 5  0.0011303568  9 1 0 149  341 1459
#> 8  0.0010360406  6 3 1  42  106 1703
#> 9  0.0017050975  9 0 1  74   74 1021
#> 10 0.0011410368  8 1 1  77   71 1531
#> 14 0.0015001050  7 0 3 807  938  288
#> 15 0.0007942764  6 1 3 148 1597 1460
#> 20 0.0011221243  9 1 0 370  725 1238

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
#> 2     s1    s2     0.1        0.00 0.05263158 0.0000000 0.7492163 0.4594595
#> 3     s1    s3     0.3        0.00 0.17647059 0.0000000 0.7261745 0.3755102
#> 4     s1    s4     0.0        0.00 0.00000000 0.0000000 0.5293375 0.2387755
#> 5     s1    s5     0.1        0.00 0.05263158 0.0000000 0.8579600 0.6959184
#> 8     s2    s3     0.4        0.25 0.25000000 0.1428571 0.9556260 0.7162162
#> 9     s2    s4     0.1        0.00 0.05263158 0.0000000 0.8809332 0.5000000
#> 10    s2    s5     0.2        0.20 0.11111111 0.1111111 0.9123007 0.4797297
#> 14    s3    s4     0.3        0.00 0.17647059 0.0000000 0.4316901 0.2630137
#> 15    s3    s5     0.4        0.25 0.25000000 0.1428571 0.9117208 0.9079602
#> 20    s4    s5     0.1        0.00 0.05263158 0.0000000 0.7262301 0.6621005
#>    Euclidean  a b c   A    B    C
#> 2   263.4616  9 1 0  80  410   68
#> 3   853.1934  7 3 0 306  184 1439
#> 4   393.0331 10 0 0 373  117  722
#> 5   883.6764  9 1 0 149  341 1459
#> 8   964.2132  6 3 1  42  106 1703
#> 9   585.4767  9 0 1  74   74 1021
#> 10  875.3959  8 1 1  77   71 1531
#> 14  665.6200  7 0 3 807  938  288
#> 15 1258.0076  6 1 3 148 1597 1460
#> 20  890.1668  9 1 0 370  725 1238
```
