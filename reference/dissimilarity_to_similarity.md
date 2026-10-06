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
#> 2     s1    s2 0.3333333   0.2500000 0.20000000 0.1428571 0.9659533 0.938704028
#> 3     s1    s3 0.1111111   0.0000000 0.05882353 0.0000000 0.5629782 0.000000000
#> 4     s1    s4 0.4000000   0.4000000 0.25000000 0.2500000 0.9525723 0.941176471
#> 5     s1    s5 0.2222222   0.2222222 0.12500000 0.1250000 0.6687756 0.629783694
#> 8     s2    s3 0.4000000   0.2500000 0.25000000 0.1428571 0.8072084 0.007005254
#> 9     s2    s4 0.3333333   0.2500000 0.20000000 0.1428571 0.3087675 0.047285464
#> 10    s2    s5 0.5000000   0.4444444 0.33333333 0.2857143 0.9774394 0.964973730
#> 14    s3    s4 0.3000000   0.2222222 0.17647059 0.1250000 0.7450111 0.197407777
#> 15    s3    s5 0.1111111   0.0000000 0.05882353 0.0000000 0.6511592 0.054908486
#> 20    s4    s5 0.2222222   0.2222222 0.12500000 0.1250000 0.9074830 0.898305085
#>    Euclidean a b c    A    B    C
#> 2  1035.4641 6 2 1   35 1450  536
#> 3  1524.1890 8 0 1 1485    0 3826
#> 4  1195.4020 6 2 2   59 1426  944
#> 5   960.5910 7 1 1  445 1040  757
#> 8  1922.9930 6 1 3  567    4 4744
#> 9   293.0051 6 1 2  544   27  459
#> 10  954.0865 5 2 3   20  551 1182
#> 14 1913.5862 7 2 1  805 4506  198
#> 15 1787.7721 8 1 0 1136 4175   66
#> 20 1116.1868 7 1 1  102  901 1100

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
#> 2     s1    s2 0.6666667   0.7500000 0.8000000 0.8571429 0.03404669 0.06129597
#> 3     s1    s3 0.8888889   1.0000000 0.9411765 1.0000000 0.43702178 1.00000000
#> 4     s1    s4 0.6000000   0.6000000 0.7500000 0.7500000 0.04742765 0.05882353
#> 5     s1    s5 0.7777778   0.7777778 0.8750000 0.8750000 0.33122441 0.37021631
#> 8     s2    s3 0.6000000   0.7500000 0.7500000 0.8571429 0.19279157 0.99299475
#> 9     s2    s4 0.6666667   0.7500000 0.8000000 0.8571429 0.69123253 0.95271454
#> 10    s2    s5 0.5000000   0.5555556 0.6666667 0.7142857 0.02256063 0.03502627
#> 14    s3    s4 0.7000000   0.7777778 0.8235294 0.8750000 0.25498891 0.80259222
#> 15    s3    s5 0.8888889   1.0000000 0.9411765 1.0000000 0.34884078 0.94509151
#> 20    s4    s5 0.7777778   0.7777778 0.8750000 0.8750000 0.09251701 0.10169492
#>       Euclidean a b c    A    B    C
#> 2  0.0009648187 6 2 1   35 1450  536
#> 3  0.0006556565 8 0 1 1485    0 3826
#> 4  0.0008358394 6 2 2   59 1426  944
#> 5  0.0010399432 7 1 1  445 1040  757
#> 8  0.0005197524 6 1 3  567    4 4744
#> 9  0.0034013013 6 1 2  544   27  459
#> 10 0.0010470256 5 2 3   20  551 1182
#> 14 0.0005223061 7 2 1  805 4506  198
#> 15 0.0005590427 8 1 0 1136 4175   66
#> 20 0.0008951054 7 1 1  102  901 1100
```
