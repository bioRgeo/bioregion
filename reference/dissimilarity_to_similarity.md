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
#>    Site1 Site2   Jaccard Jaccardturn   Sorensen   Simpson      Bray    Brayturn
#> 2     s1    s2 0.3333333   0.0000000 0.20000000 0.0000000 0.8376812 0.770491803
#> 3     s1    s3 0.4000000   0.0000000 0.25000000 0.0000000 0.8650576 0.779372197
#> 4     s1    s4 0.6000000   0.5000000 0.42857143 0.3333333 0.8733866 0.798828125
#> 5     s1    s5 0.4000000   0.0000000 0.25000000 0.0000000 0.9630390 0.921965318
#> 8     s2    s3 0.1000000   0.0000000 0.05263158 0.0000000 0.7408469 0.332786885
#> 9     s2    s4 0.3000000   0.2222222 0.17647059 0.1250000 0.5525847 0.509765625
#> 10    s2    s5 0.1000000   0.0000000 0.05263158 0.0000000 0.8933054 0.852601156
#> 14    s3    s4 0.2000000   0.0000000 0.11111111 0.0000000 0.7186987 0.164062500
#> 15    s3    s5 0.0000000   0.0000000 0.00000000 0.0000000 0.7608620 0.005780347
#> 20    s4    s5 0.2000000   0.0000000 0.11111111 0.0000000 0.4848485 0.361271676
#>    Euclidean  a b c   A    B    C
#> 2   921.9192  6 0 3 140  975  470
#> 3  1419.8761  6 0 4 246  869 2285
#> 4   981.7734  4 2 4 103 1012  409
#> 5  1024.6809  6 0 4  27 1088  319
#> 8  1141.2909  9 0 1 407  203 2124
#> 9   271.4295  7 2 1 251  359  261
#> 10  427.1112  9 0 1  51  559  295
#> 14 1100.5658  8 2 0 428 2103   84
#> 15 1116.5809 10 0 0 344 2187    2
#> 20  193.4322  8 0 2 221  291  125

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
#> 2     s1    s2 0.6666667   1.0000000 0.8000000 1.0000000 0.16231884 0.22950820
#> 3     s1    s3 0.6000000   1.0000000 0.7500000 1.0000000 0.13494240 0.22062780
#> 4     s1    s4 0.4000000   0.5000000 0.5714286 0.6666667 0.12661340 0.20117188
#> 5     s1    s5 0.6000000   1.0000000 0.7500000 1.0000000 0.03696099 0.07803468
#> 8     s2    s3 0.9000000   1.0000000 0.9473684 1.0000000 0.25915314 0.66721311
#> 9     s2    s4 0.7000000   0.7777778 0.8235294 0.8750000 0.44741533 0.49023438
#> 10    s2    s5 0.9000000   1.0000000 0.9473684 1.0000000 0.10669456 0.14739884
#> 14    s3    s4 0.8000000   1.0000000 0.8888889 1.0000000 0.28130135 0.83593750
#> 15    s3    s5 1.0000000   1.0000000 1.0000000 1.0000000 0.23913799 0.99421965
#> 20    s4    s5 0.8000000   1.0000000 0.8888889 1.0000000 0.51515152 0.63872832
#>       Euclidean  a b c   A    B    C
#> 2  0.0010835185  6 0 3 140  975  470
#> 3  0.0007037912  6 0 4 246  869 2285
#> 4  0.0010175286  4 2 4 103 1012  409
#> 5  0.0009749621  6 0 4  27 1088  319
#> 8  0.0008754337  9 0 1 407  203 2124
#> 9  0.0036706738  7 2 1 251  359  261
#> 10 0.0023358416  9 0 1  51  559  295
#> 14 0.0009077987  8 2 0 428 2103   84
#> 15 0.0008947898 10 0 0 344 2187    2
#> 20 0.0051431821  8 0 2 221  291  125
```
