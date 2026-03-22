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
#>    Site1 Site2   Jaccard Jaccardturn   Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2 0.1111111   0.0000000 0.05882353 0.0000000 0.8867588 0.8813084
#> 3     s1    s3 0.2000000   0.2000000 0.11111111 0.1111111 0.9180088 0.8806479
#> 4     s1    s4 0.1000000   0.0000000 0.05263158 0.0000000 0.9252802 0.9232737
#> 5     s1    s5 0.3000000   0.2222222 0.17647059 0.1250000 0.8639281 0.8634021
#> 8     s2    s3 0.3000000   0.2222222 0.17647059 0.1250000 0.5628019 0.3233645
#> 9     s2    s4 0.2000000   0.0000000 0.11111111 0.0000000 0.7614918 0.7429907
#> 10    s2    s5 0.4000000   0.4000000 0.25000000 0.2500000 0.8272158 0.8196262
#> 14    s3    s4 0.1000000   0.0000000 0.05263158 0.0000000 0.8217366 0.7491909
#> 15    s3    s5 0.1111111   0.0000000 0.05882353 0.0000000 0.8068115 0.7173540
#> 20    s4    s5 0.2000000   0.0000000 0.11111111 0.0000000 0.9666667 0.9656357
#>    Euclidean a b c   A    B    C
#> 2  1054.2822 8 1 0 127 1046  943
#> 3  1368.0040 8 1 1 140 1033 2102
#> 4  1120.2370 9 0 1  90 1083 1146
#> 5  1164.6669 7 2 1 159 1014 1005
#> 8   944.2245 7 1 2 724  346 1518
#> 9   802.8474 8 0 2 275  795  961
#> 10 1047.7843 6 2 2 193  877  971
#> 14 1217.2280 9 0 1 310 1932  926
#> 15 1277.3739 8 1 0 329 1913  835
#> 20 1082.3881 8 2 0  40 1196 1124

similarity <- dissimilarity_to_similarity(dissimil)
similarity
#> Data.frame of similarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of similarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen   Simpson       Bray   Brayturn
#> 2     s1    s2 0.8888889   1.0000000 0.9411765 1.0000000 0.11324119 0.11869159
#> 3     s1    s3 0.8000000   0.8000000 0.8888889 0.8888889 0.08199122 0.11935209
#> 4     s1    s4 0.9000000   1.0000000 0.9473684 1.0000000 0.07471980 0.07672634
#> 5     s1    s5 0.7000000   0.7777778 0.8235294 0.8750000 0.13607189 0.13659794
#> 8     s2    s3 0.7000000   0.7777778 0.8235294 0.8750000 0.43719807 0.67663551
#> 9     s2    s4 0.8000000   1.0000000 0.8888889 1.0000000 0.23850824 0.25700935
#> 10    s2    s5 0.6000000   0.6000000 0.7500000 0.7500000 0.17278424 0.18037383
#> 14    s3    s4 0.9000000   1.0000000 0.9473684 1.0000000 0.17826337 0.25080906
#> 15    s3    s5 0.8888889   1.0000000 0.9411765 1.0000000 0.19318849 0.28264605
#> 20    s4    s5 0.8000000   1.0000000 0.8888889 1.0000000 0.03333333 0.03436426
#>       Euclidean a b c   A    B    C
#> 2  0.0009476138 8 1 0 127 1046  943
#> 3  0.0007304580 8 1 1 140 1033 2102
#> 4  0.0008918721 9 0 1  90 1083 1146
#> 5  0.0008578780 7 2 1 159 1014 1005
#> 8  0.0010579497 7 1 2 724  346 1518
#> 9  0.0012440172 8 0 2 275  795  961
#> 10 0.0009534849 6 2 2 193  877  971
#> 14 0.0008208644 9 0 1 310 1932  926
#> 15 0.0007822438 8 1 0 329 1913  835
#> 20 0.0009230303 8 2 0  40 1196 1124
```
