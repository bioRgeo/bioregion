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
#> 2     s1    s2 0.7777778   0.7777778 0.8750000   0.875 0.03769740 0.03936170
#> 3     s1    s3 0.7777778   0.7777778 0.8750000   0.875 0.26734843 0.52890173
#> 4     s1    s4 0.7000000   0.7777778 0.8235294   0.875 0.13344739 0.15249267
#> 5     s1    s5 0.6000000   0.6000000 0.7500000   0.750 0.04441584 0.06744868
#> 8     s2    s3 0.6000000   0.6000000 0.7500000   0.750 0.09331260 0.17341040
#> 9     s2    s4 0.7000000   0.7777778 0.8235294   0.875 0.10110865 0.12127660
#> 10    s2    s5 0.6000000   0.6000000 0.7500000   0.750 0.04166667 0.06702128
#> 14    s3    s4 0.7000000   0.7777778 0.8235294   0.875 0.32269717 0.77456647
#> 15    s3    s5 0.6000000   0.6000000 0.7500000   0.750 0.15637860 0.54913295
#> 20    s4    s5 0.7000000   0.7777778 0.8235294   0.875 0.51603413 0.66692015
#>       Euclidean a b c   A   B    C
#> 2  0.0010230218 7 1 1  37 986  903
#> 3  0.0017174666 7 1 1 183 840  163
#> 4  0.0011814787 7 1 2 156 867 1159
#> 5  0.0007334472 6 2 2  69 954 2015
#> 8  0.0013779355 6 2 2  60 880  286
#> 9  0.0011237262 7 1 2 114 826 1201
#> 10 0.0007315406 6 2 2  63 877 2021
#> 14 0.0018739105 7 1 2 268  78 1047
#> 15 0.0008794258 6 2 2 190 156 1894
#> 20 0.0012955458 7 2 1 877 438 1207

dissimilarity <- similarity_to_dissimilarity(simil)
dissimilarity
#> Data.frame of dissimilarity between sites
#>  - Total number of sites:  5 
#>  - Total number of species:  10 
#>  - Number of rows:  10 
#>  - Number of dissimilarity metrics:  7 
#> 
#> 
#>    Site1 Site2   Jaccard Jaccardturn  Sorensen Simpson      Bray  Brayturn
#> 2     s1    s2 0.2222222   0.2222222 0.1250000   0.125 0.9623026 0.9606383
#> 3     s1    s3 0.2222222   0.2222222 0.1250000   0.125 0.7326516 0.4710983
#> 4     s1    s4 0.3000000   0.2222222 0.1764706   0.125 0.8665526 0.8475073
#> 5     s1    s5 0.4000000   0.4000000 0.2500000   0.250 0.9555842 0.9325513
#> 8     s2    s3 0.4000000   0.4000000 0.2500000   0.250 0.9066874 0.8265896
#> 9     s2    s4 0.3000000   0.2222222 0.1764706   0.125 0.8988914 0.8787234
#> 10    s2    s5 0.4000000   0.4000000 0.2500000   0.250 0.9583333 0.9329787
#> 14    s3    s4 0.3000000   0.2222222 0.1764706   0.125 0.6773028 0.2254335
#> 15    s3    s5 0.4000000   0.4000000 0.2500000   0.250 0.8436214 0.4508671
#> 20    s4    s5 0.3000000   0.2222222 0.1764706   0.125 0.4839659 0.3330798
#>    Euclidean a b c   A   B    C
#> 2   976.4963 7 1 1  37 986  903
#> 3   581.2530 7 1 1 183 840  163
#> 4   845.3969 7 1 2 156 867 1159
#> 5  1362.4247 6 2 2  69 954 2015
#> 8   724.7234 6 2 2  60 880  286
#> 9   888.8965 7 1 2 114 826 1201
#> 10 1365.9780 6 2 2  63 877 2021
#> 14  532.6434 7 1 2 268  78 1047
#> 15 1136.1056 6 2 2 190 156 1894
#> 20  770.8755 7 2 1 877 438 1207
```
