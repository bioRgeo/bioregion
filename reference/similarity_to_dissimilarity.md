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
#>    Site1 Site2 Jaccard Jaccardturn  Sorensen   Simpson       Bray   Brayturn
#> 2     s1    s2    0.50   0.5555556 0.6666667 0.7142857 0.06253006 0.08322663
#> 3     s1    s3    0.50   0.5555556 0.6666667 0.7142857 0.31105901 0.31562252
#> 4     s1    s4    0.80   1.0000000 0.8888889 1.0000000 0.07371484 0.09947644
#> 5     s1    s5    0.70   0.7777778 0.8235294 0.8750000 0.41261667 0.49383667
#> 8     s2    s3    0.75   0.7500000 0.8571429 0.8571429 0.05386876 0.07042254
#> 9     s2    s4    0.70   1.0000000 0.8235294 1.0000000 0.67831715 0.68586387
#> 10    s2    s5    0.60   0.7500000 0.7500000 0.8571429 0.09266409 0.15364917
#> 14    s3    s4    0.70   1.0000000 0.8235294 1.0000000 0.06716049 0.08900524
#> 15    s3    s5    0.60   0.7500000 0.7500000 0.8571429 0.12573290 0.15305313
#> 20    s4    s5    0.90   1.0000000 0.9473684 1.0000000 0.20287602 0.34162304
#>       Euclidean a b c   A    B    C
#> 2  0.0010904639 5 3 2  65 1233  716
#> 3  0.0010783034 5 3 2 398  900  863
#> 4  0.0011471735 8 0 2  76 1222  688
#> 5  0.0013534865 7 1 2 641  657 1168
#> 8  0.0008731678 6 1 1  55  726 1206
#> 9  0.0042734259 7 0 3 524  257  240
#> 10 0.0008903278 6 1 3 120  661 1689
#> 14 0.0009231102 7 0 3  68 1193  696
#> 15 0.0007119358 6 1 3 193 1068 1616
#> 20 0.0009615634 9 1 0 261  503 1548

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
#> 2     s1    s2    0.50   0.4444444 0.33333333 0.2857143 0.9374699 0.9167734
#> 3     s1    s3    0.50   0.4444444 0.33333333 0.2857143 0.6889410 0.6843775
#> 4     s1    s4    0.20   0.0000000 0.11111111 0.0000000 0.9262852 0.9005236
#> 5     s1    s5    0.30   0.2222222 0.17647059 0.1250000 0.5873833 0.5061633
#> 8     s2    s3    0.25   0.2500000 0.14285714 0.1428571 0.9461312 0.9295775
#> 9     s2    s4    0.30   0.0000000 0.17647059 0.0000000 0.3216828 0.3141361
#> 10    s2    s5    0.40   0.2500000 0.25000000 0.1428571 0.9073359 0.8463508
#> 14    s3    s4    0.30   0.0000000 0.17647059 0.0000000 0.9328395 0.9109948
#> 15    s3    s5    0.40   0.2500000 0.25000000 0.1428571 0.8742671 0.8469469
#> 20    s4    s5    0.10   0.0000000 0.05263158 0.0000000 0.7971240 0.6583770
#>    Euclidean a b c   A    B    C
#> 2   916.0409 5 3 2  65 1233  716
#> 3   926.3828 5 3 2 398  900  863
#> 4   870.7078 8 0 2  76 1222  688
#> 5   737.8326 7 1 2 641  657 1168
#> 8  1144.2552 6 1 1  55  726 1206
#> 9   233.0043 7 0 3 524  257  240
#> 10 1122.1818 6 1 3 120  661 1689
#> 14 1082.2943 7 0 3  68 1193  696
#> 15 1403.6210 6 1 3 193 1068 1616
#> 20 1038.9731 9 1 0 261  503 1548
```
