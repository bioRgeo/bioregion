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
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen Simpson       Bray   Brayturn
#> 2     s1    s2 0.9000000   1.0000000 0.9473684   1.000 0.23777506 0.71639042
#> 3     s1    s3 0.8000000   1.0000000 0.8888889   1.000 0.08013449 0.17023810
#> 4     s1    s4 0.8000000   1.0000000 0.8888889   1.000 0.22286477 0.28353141
#> 5     s1    s5 0.8000000   1.0000000 0.8888889   1.000 0.13740960 0.66773163
#> 8     s2    s3 0.7000000   0.7777778 0.8235294   0.875 0.33405640 0.42541436
#> 9     s2    s4 0.7000000   0.7777778 0.8235294   0.875 0.07705628 0.16390424
#> 10    s2    s5 0.8888889   1.0000000 0.9411765   1.000 0.51635514 0.70607029
#> 14    s3    s4 0.6000000   0.6000000 0.7500000   0.750 0.01380898 0.02142857
#> 15    s3    s5 0.6000000   0.6000000 0.7500000   0.750 0.20468343 0.37699681
#> 20    s4    s5 0.6000000   0.6000000 0.7500000   0.750 0.10000000 0.33226837
#>       Euclidean a b c   A    B    C
#> 2  0.0008120233 9 1 0 389 2340  154
#> 3  0.0006859250 8 2 0 143 2586  697
#> 4  0.0006581515 8 2 0 501 2228 1266
#> 5  0.0007696664 8 2 0 209 2520  104
#> 8  0.0018542552 7 2 1 231  312  609
#> 9  0.0009004376 7 2 1  89  454 1678
#> 10 0.0053594716 8 1 0 221  322   92
#> 14 0.0007889834 6 2 2  18  822 1749
#> 15 0.0016270313 6 2 2 118  722  195
#> 20 0.0009005633 6 2 2 104 1663  209

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
#> 2     s1    s2 0.1000000   0.0000000 0.05263158   0.000 0.7622249 0.2836096
#> 3     s1    s3 0.2000000   0.0000000 0.11111111   0.000 0.9198655 0.8297619
#> 4     s1    s4 0.2000000   0.0000000 0.11111111   0.000 0.7771352 0.7164686
#> 5     s1    s5 0.2000000   0.0000000 0.11111111   0.000 0.8625904 0.3322684
#> 8     s2    s3 0.3000000   0.2222222 0.17647059   0.125 0.6659436 0.5745856
#> 9     s2    s4 0.3000000   0.2222222 0.17647059   0.125 0.9229437 0.8360958
#> 10    s2    s5 0.1111111   0.0000000 0.05882353   0.000 0.4836449 0.2939297
#> 14    s3    s4 0.4000000   0.4000000 0.25000000   0.250 0.9861910 0.9785714
#> 15    s3    s5 0.4000000   0.4000000 0.25000000   0.250 0.7953166 0.6230032
#> 20    s4    s5 0.4000000   0.4000000 0.25000000   0.250 0.9000000 0.6677316
#>    Euclidean a b c   A    B    C
#> 2  1230.4918 9 1 0 389 2340  154
#> 3  1456.8854 8 2 0 143 2586  697
#> 4  1518.4071 8 2 0 501 2228 1266
#> 5  1298.2642 8 2 0 209 2520  104
#> 8   538.3001 7 2 1 231  312  609
#> 9  1109.5711 7 2 1  89  454 1678
#> 10  185.5856 8 1 0 221  322   92
#> 14 1266.4537 6 2 2  18  822 1749
#> 15  613.6163 6 2 2 118  722  195
#> 20 1109.4161 6 2 2 104 1663  209
```
