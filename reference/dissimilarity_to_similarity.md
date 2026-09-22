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
#>    Site1 Site2 Jaccard Jaccardturn   Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2     0.1   0.0000000 0.05263158 0.0000000 0.4074559 0.1996466
#> 3     s1    s3     0.3   0.2222222 0.17647059 0.1250000 0.7809187 0.2662722
#> 4     s1    s4     0.2   0.2000000 0.11111111 0.1111111 0.8180058 0.8047767
#> 5     s1    s5     0.1   0.0000000 0.05263158 0.0000000 0.7957046 0.7938689
#> 8     s2    s3     0.2   0.0000000 0.11111111 0.0000000 0.6136054 0.1597633
#> 9     s2    s4     0.1   0.0000000 0.05263158 0.0000000 0.7016177 0.5600707
#> 10    s2    s5     0.0   0.0000000 0.00000000 0.0000000 0.8240741 0.7650177
#> 14    s3    s4     0.3   0.2222222 0.17647059 0.1250000 0.8946541 0.6035503
#> 15    s3    s5     0.2   0.0000000 0.11111111 0.0000000 0.8923767 0.6449704
#> 20    s4    s5     0.1   0.0000000 0.05263158 0.0000000 0.6642265 0.6363636
#>    Euclidean  a b c   A   B    C
#> 2   308.6341  9 0 1 453 510  113
#> 3   475.1147  7 2 1 124 839   45
#> 4   758.7345  8 1 1 188 775  915
#> 5   684.5444  9 0 1 195 768  751
#> 8   222.7442  8 2 0 142 424   27
#> 9   579.9302  9 1 0 249 317  854
#> 10  486.5450 10 0 0 133 433  813
#> 14  613.3857  7 1 2  67 102 1036
#> 15  425.3504  8 0 2  60 109  886
#> 20  642.2904  9 0 1 344 759  602

similarity <- dissimilarity_to_similarity(dissimil)
similarity
#> Data.frame of similarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of similarity metrics:  7 
#> 
#> 
#>    Site1 Site2 Jaccard Jaccardturn  Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2     0.9   1.0000000 0.9473684 1.0000000 0.5925441 0.8003534
#> 3     s1    s3     0.7   0.7777778 0.8235294 0.8750000 0.2190813 0.7337278
#> 4     s1    s4     0.8   0.8000000 0.8888889 0.8888889 0.1819942 0.1952233
#> 5     s1    s5     0.9   1.0000000 0.9473684 1.0000000 0.2042954 0.2061311
#> 8     s2    s3     0.8   1.0000000 0.8888889 1.0000000 0.3863946 0.8402367
#> 9     s2    s4     0.9   1.0000000 0.9473684 1.0000000 0.2983823 0.4399293
#> 10    s2    s5     1.0   1.0000000 1.0000000 1.0000000 0.1759259 0.2349823
#> 14    s3    s4     0.7   0.7777778 0.8235294 0.8750000 0.1053459 0.3964497
#> 15    s3    s5     0.8   1.0000000 0.8888889 1.0000000 0.1076233 0.3550296
#> 20    s4    s5     0.9   1.0000000 0.9473684 1.0000000 0.3357735 0.3636364
#>      Euclidean  a b c   A   B    C
#> 2  0.003229619  9 0 1 453 510  113
#> 3  0.002100334  7 2 1 124 839   45
#> 4  0.001316249  8 1 1 188 775  915
#> 5  0.001458695  9 0 1 195 768  751
#> 8  0.004469389  8 2 0 142 424   27
#> 9  0.001721377  9 1 0 249 317  854
#> 10 0.002051093 10 0 0 133 433  813
#> 14 0.001627642  7 1 2  67 102 1036
#> 15 0.002345488  8 0 2  60 109  886
#> 20 0.001554508  9 0 1 344 759  602
```
