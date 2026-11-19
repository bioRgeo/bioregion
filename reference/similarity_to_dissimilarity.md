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
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2 0.7000000   0.7777778 0.8235294 0.8750000 0.3521514 0.6461366
#> 3     s1    s3 0.8000000   0.8000000 0.8888889 0.8888889 0.2485844 0.3824042
#> 4     s1    s4 0.6000000   0.7500000 0.7500000 0.8571429 0.4523741 0.7204400
#> 5     s1    s5 0.8888889   1.0000000 0.9411765 1.0000000 0.1906336 0.2776886
#> 8     s2    s3 0.7000000   0.7777778 0.8235294 0.8750000 0.4115630 0.4703247
#> 9     s2    s4 0.6666667   0.7500000 0.8000000 0.8571429 0.4264113 0.4736842
#> 10    s2    s5 0.6000000   0.6000000 0.7500000 0.7500000 0.2917251 0.3493841
#> 14    s3    s4 0.6000000   0.7500000 0.7500000 0.8571429 0.3126396 0.3208066
#> 15    s3    s5 0.7000000   0.7777778 0.8235294 0.8750000 0.1737678 0.1811847
#> 20    s4    s5 0.6666667   0.7500000 0.8000000 0.8571429 0.3252033 0.3483043
#>       Euclidean a b c   A    B    C
#> 2  0.0009783967 7 2 1 577 1807  316
#> 3  0.0008281358 8 1 1 439 1945  709
#> 4  0.0009848376 6 3 1 786 1598  305
#> 5  0.0007202679 8 1 0 346 2038  900
#> 8  0.0022786807 7 1 2 420  473  728
#> 9  0.0019728019 6 2 1 423  470  668
#> 10 0.0011094518 6 2 2 312  581  934
#> 14 0.0015772989 6 3 1 350  798  741
#> 15 0.0010265430 7 2 1 208  940 1038
#> 20 0.0011061405 6 1 2 380  711  866

dissimilarity <- similarity_to_dissimilarity(simil)
dissimilarity
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn   Sorensen   Simpson      Bray  Brayturn
#> 2     s1    s2 0.3000000   0.2222222 0.17647059 0.1250000 0.6478486 0.3538634
#> 3     s1    s3 0.2000000   0.2000000 0.11111111 0.1111111 0.7514156 0.6175958
#> 4     s1    s4 0.4000000   0.2500000 0.25000000 0.1428571 0.5476259 0.2795600
#> 5     s1    s5 0.1111111   0.0000000 0.05882353 0.0000000 0.8093664 0.7223114
#> 8     s2    s3 0.3000000   0.2222222 0.17647059 0.1250000 0.5884370 0.5296753
#> 9     s2    s4 0.3333333   0.2500000 0.20000000 0.1428571 0.5735887 0.5263158
#> 10    s2    s5 0.4000000   0.4000000 0.25000000 0.2500000 0.7082749 0.6506159
#> 14    s3    s4 0.4000000   0.2500000 0.25000000 0.1428571 0.6873604 0.6791934
#> 15    s3    s5 0.3000000   0.2222222 0.17647059 0.1250000 0.8262322 0.8188153
#> 20    s4    s5 0.3333333   0.2500000 0.20000000 0.1428571 0.6747967 0.6516957
#>    Euclidean a b c   A    B    C
#> 2  1021.0803 7 2 1 577 1807  316
#> 3  1206.5314 8 1 1 439 1945  709
#> 4  1014.3959 6 3 1 786 1598  305
#> 5  1387.3723 8 1 0 346 2038  900
#> 8   437.8504 7 1 2 420  473  728
#> 9   505.8933 6 2 1 423  470  668
#> 10  900.3460 6 2 2 312  581  934
#> 14  632.9953 6 3 1 350  798  741
#> 15  973.1434 7 2 1 208  940 1038
#> 20  903.0443 6 1 2 380  711  866
```
