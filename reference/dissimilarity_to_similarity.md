# Convert dissimilarity metrics to similarity metrics

This function converts a `data.frame` of dissimilarity metrics (beta
diversity) between sites into similarity metrics.

## Usage

``` r
dissimilarity_to_similarity(dissimilarity, include_formula = TRUE)
```

## Arguments

- dissimilarity:

  the output object from
  [`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
  or
  [`similarity_to_dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md).

- include_formula:

  a `boolean` indicating whether metrics based on custom formula(s)
  should also be converted (see Details). The default is `TRUE`.

## Value

A `data.frame` with the additional class `bioregion.pairwise`, providing
similarity metrics for each pair of sites based on a dissimilarity
object.

## Note

The behavior of this function changes depending on column names. Columns
`Site1` and `Site2` are copied identically. If there are columns called
`a`, `b`, `c`, `A`, `B`, `C` they will also be copied identically. If
there are columns based on your own formula (argument `formula` in
[`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md))
or not in the original list of dissimilarity metrics (argument `metrics`
in
[`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md))
and if the argument `include_formula` is set to `FALSE`, they will also
be copied identically. Otherwise there are going to be converted like
they other columns (default behavior).

If a column is called `Euclidean`, the similarity will be calculated
based on the following formula:

Euclidean similarity = 1 / (1 - Euclidean distance)

Otherwise, all other columns will be transformed into dissimilarity with
the following formula:

similarity = 1 - dissimilarity

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html>.

Associated functions:
[similarity](https://bioRgeo.github.io/bioregion/reference/similarity.md)
dissimilarity_to_similarity

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

dissimil <- dissimilarity(comat, metric = "all")
dissimil
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2 0.2000000   0.0000000 0.1111111 0.0000000 0.5068871 0.1516588
#> 3     s1    s3 0.2000000   0.0000000 0.1111111 0.0000000 0.8447894 0.8341232
#> 4     s1    s4 0.2000000   0.0000000 0.1111111 0.0000000 0.6442396 0.5426540
#> 5     s1    s5 0.3000000   0.0000000 0.1764706 0.0000000 0.8426362 0.7227488
#> 8     s2    s3 0.4000000   0.4000000 0.2500000 0.2500000 0.9099338 0.8583333
#> 9     s2    s4 0.2222222   0.2222222 0.1250000 0.1250000 0.5735381 0.4555053
#> 10    s2    s5 0.3333333   0.2500000 0.2000000 0.1428571 0.7011933 0.6961165
#> 14    s3    s4 0.2222222   0.2222222 0.1250000 0.1250000 0.8285214 0.7958333
#> 15    s3    s5 0.3333333   0.2500000 0.2000000 0.1428571 0.4446602 0.1062500
#> 20    s4    s5 0.5000000   0.4444444 0.3333333 0.2857143 0.7546296 0.6802413
#>    Euclidean a b c   A   B   C
#> 2   333.8053 8 2 0 358  64 672
#> 3   399.5973 8 2 0  70 352 410
#> 4   310.5978 8 2 0 193 229 470
#> 5   609.6925 7 3 0 117 305 948
#> 8   611.8660 6 2 2  68 962 412
#> 9   485.2061 7 1 1 361 669 302
#> 10  657.9126 6 2 1 313 717 752
#> 14  456.2510 7 1 1  98 382 565
#> 15  395.4605 6 2 1 429  51 636
#> 20  590.0085 5 3 2 212 451 853

similarity <- dissimilarity_to_similarity(dissimil)
similarity
#> Data.frame of similarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of similarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen   Simpson       Bray  Brayturn
#> 2     s1    s2 0.8000000   1.0000000 0.8888889 1.0000000 0.49311295 0.8483412
#> 3     s1    s3 0.8000000   1.0000000 0.8888889 1.0000000 0.15521064 0.1658768
#> 4     s1    s4 0.8000000   1.0000000 0.8888889 1.0000000 0.35576037 0.4573460
#> 5     s1    s5 0.7000000   1.0000000 0.8235294 1.0000000 0.15736382 0.2772512
#> 8     s2    s3 0.6000000   0.6000000 0.7500000 0.7500000 0.09006623 0.1416667
#> 9     s2    s4 0.7777778   0.7777778 0.8750000 0.8750000 0.42646190 0.5444947
#> 10    s2    s5 0.6666667   0.7500000 0.8000000 0.8571429 0.29880668 0.3038835
#> 14    s3    s4 0.7777778   0.7777778 0.8750000 0.8750000 0.17147857 0.2041667
#> 15    s3    s5 0.6666667   0.7500000 0.8000000 0.8571429 0.55533981 0.8937500
#> 20    s4    s5 0.5000000   0.5555556 0.6666667 0.7142857 0.24537037 0.3197587
#>      Euclidean a b c   A   B   C
#> 2  0.002986810 8 2 0 358  64 672
#> 3  0.002496272 8 2 0  70 352 410
#> 4  0.003209265 8 2 0 193 229 470
#> 5  0.001637485 7 3 0 117 305 948
#> 8  0.001631678 6 2 2  68 962 412
#> 9  0.002056741 7 1 1 361 669 302
#> 10 0.001517652 6 2 1 313 717 752
#> 14 0.002186983 7 1 1  98 382 565
#> 15 0.002522319 6 2 1 429  51 636
#> 20 0.001692023 5 3 2 212 451 853
```
