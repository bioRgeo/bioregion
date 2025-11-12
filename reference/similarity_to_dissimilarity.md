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
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen Simpson       Bray  Brayturn
#> 2     s1    s2 0.7777778   0.7777778 0.8750000   0.875 0.28300259 0.5754386
#> 3     s1    s3 0.7000000   0.7777778 0.8235294   0.875 0.26832845 0.6421053
#> 4     s1    s4 0.7777778   0.7777778 0.8750000   0.875 0.26045884 0.6771930
#> 5     s1    s5 0.8000000   1.0000000 0.8888889   1.000 0.09679666 0.4877193
#> 8     s2    s3 0.7000000   0.7777778 0.8235294   0.875 0.32360471 0.3615561
#> 9     s2    s4 0.7777778   0.7777778 0.8750000   0.875 0.43650410 0.5171625
#> 10    s2    s5 0.8000000   1.0000000 0.8888889   1.000 0.22016758 0.4359268
#> 14    s3    s4 0.7000000   0.7777778 0.8235294   0.875 0.26274165 0.2771084
#> 15    s3    s5 0.9000000   1.0000000 0.9473684   1.000 0.07801418 0.1325301
#> 20    s4    s5 0.8000000   1.0000000 0.8888889   1.000 0.09302326 0.1470343
#>       Euclidean a b c   A    B    C
#> 2  0.0023121341 7 1 1 164  121  710
#> 3  0.0017037657 7 1 2 183  102  896
#> 4  0.0013878424 7 1 1 193   92 1004
#> 5  0.0008442520 8 0 2 139  146 2448
#> 8  0.0014929371 7 1 2 316  558  763
#> 9  0.0015926682 7 1 1 452  422  745
#> 10 0.0008229557 8 0 2 381  493 2206
#> 14 0.0011348981 7 2 1 299  780  898
#> 15 0.0007504370 9 0 1 143  936 2444
#> 20 0.0006782972 8 0 2 176 1021 2411

dissimilarity <- similarity_to_dissimilarity(simil)
dissimilarity
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn   Sorensen Simpson      Bray  Brayturn
#> 2     s1    s2 0.2222222   0.2222222 0.12500000   0.125 0.7169974 0.4245614
#> 3     s1    s3 0.3000000   0.2222222 0.17647059   0.125 0.7316716 0.3578947
#> 4     s1    s4 0.2222222   0.2222222 0.12500000   0.125 0.7395412 0.3228070
#> 5     s1    s5 0.2000000   0.0000000 0.11111111   0.000 0.9032033 0.5122807
#> 8     s2    s3 0.3000000   0.2222222 0.17647059   0.125 0.6763953 0.6384439
#> 9     s2    s4 0.2222222   0.2222222 0.12500000   0.125 0.5634959 0.4828375
#> 10    s2    s5 0.2000000   0.0000000 0.11111111   0.000 0.7798324 0.5640732
#> 14    s3    s4 0.3000000   0.2222222 0.17647059   0.125 0.7372583 0.7228916
#> 15    s3    s5 0.1000000   0.0000000 0.05263158   0.000 0.9219858 0.8674699
#> 20    s4    s5 0.2000000   0.0000000 0.11111111   0.000 0.9069767 0.8529657
#>    Euclidean a b c   A    B    C
#> 2   431.5009 7 1 1 164  121  710
#> 3   585.9351 7 1 2 183  102  896
#> 4   719.5429 7 1 1 193   92 1004
#> 5  1183.4805 8 0 2 139  146 2448
#> 8   668.8206 7 1 2 316  558  763
#> 9   626.8772 7 1 1 452  422  745
#> 10 1214.1322 8 0 2 381  493 2206
#> 14  880.1364 7 2 1 299  780  898
#> 15 1331.5570 9 0 1 143  936 2444
#> 20 1473.2800 8 0 2 176 1021 2411
```
