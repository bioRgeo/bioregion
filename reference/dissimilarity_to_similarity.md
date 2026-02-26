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
#> 2     s1    s2     0.2   0.2000000 0.11111111 0.1111111 0.6091603 0.4498567
#> 3     s1    s3     0.3   0.2222222 0.17647059 0.1250000 0.5096744 0.3902582
#> 4     s1    s4     0.0   0.0000000 0.00000000 0.0000000 0.6015850 0.5635359
#> 5     s1    s5     0.1   0.0000000 0.05263158 0.0000000 0.5075301 0.3234483
#> 8     s2    s3     0.3   0.2222222 0.17647059 0.1250000 0.6922581 0.6583095
#> 9     s2    s4     0.2   0.2000000 0.11111111 0.1111111 0.4150430 0.0752149
#> 10    s2    s5     0.1   0.0000000 0.05263158 0.0000000 0.7224174 0.7170487
#> 14    s3    s4     0.3   0.2222222 0.17647059 0.1250000 0.7174926 0.6085681
#> 15    s3    s5     0.2   0.0000000 0.11111111 0.0000000 0.3551046 0.2986207
#> 20    s4    s5     0.1   0.0000000 0.05263158 0.0000000 0.6481647 0.4579310
#>    Euclidean a b c    A    B    C
#> 2  1091.2846 8 1 1  768 1766  628
#> 3   859.1438 7 2 1 1039 1495  665
#> 4  1425.3350 9 0 0 1106 1428 1912
#> 5   891.0376 9 0 1  981 1553  469
#> 8  1113.2457 7 2 1  477  919 1227
#> 9  1038.9976 8 1 1 1291  105 1727
#> 10 1083.7629 9 0 1  395 1001 1055
#> 14 1513.4973 7 1 2  667 1037 2351
#> 15  512.4080 8 0 2 1017  687  433
#> 20 1421.0961 9 0 1  786 2232  664

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
#> 2     s1    s2     0.8   0.8000000 0.8888889 0.8888889 0.3908397 0.5501433
#> 3     s1    s3     0.7   0.7777778 0.8235294 0.8750000 0.4903256 0.6097418
#> 4     s1    s4     1.0   1.0000000 1.0000000 1.0000000 0.3984150 0.4364641
#> 5     s1    s5     0.9   1.0000000 0.9473684 1.0000000 0.4924699 0.6765517
#> 8     s2    s3     0.7   0.7777778 0.8235294 0.8750000 0.3077419 0.3416905
#> 9     s2    s4     0.8   0.8000000 0.8888889 0.8888889 0.5849570 0.9247851
#> 10    s2    s5     0.9   1.0000000 0.9473684 1.0000000 0.2775826 0.2829513
#> 14    s3    s4     0.7   0.7777778 0.8235294 0.8750000 0.2825074 0.3914319
#> 15    s3    s5     0.8   1.0000000 0.8888889 1.0000000 0.6448954 0.7013793
#> 20    s4    s5     0.9   1.0000000 0.9473684 1.0000000 0.3518353 0.5420690
#>       Euclidean a b c    A    B    C
#> 2  0.0009155123 8 1 1  768 1766  628
#> 3  0.0011625964 7 2 1 1039 1495  665
#> 4  0.0007010975 9 0 0 1106 1428 1912
#> 5  0.0011210290 9 0 1  981 1553  469
#> 8  0.0008974681 7 2 1  477  919 1227
#> 9  0.0009615407 8 1 1 1291  105 1727
#> 10 0.0009218604 9 0 1  395 1001 1055
#> 14 0.0006602851 7 1 2  667 1037 2351
#> 15 0.0019477685 8 0 2 1017  687  433
#> 20 0.0007031874 9 0 1  786 2232  664
```
