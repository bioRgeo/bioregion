# 2. Matrix and network formats

## 1. Load data

The `bioregion`’s package contains as example dataset the spatial
distribution of Mediterranean vegetation. This dataset has been analyzed
in [this article](https://onlinelibrary.wiley.com/doi/10.1002/ece3.4718)
and contains the abundance of 3,697 species in 715 sites. This dataset
is composed of three files,
[vegedf](https://bioRgeo.github.io/bioregion/reference/vegedf.html) a
`data.frame` with 460,878 rows and 3 columns (Site, Species and
Abundance),

``` r
data(vegedf)
head(vegedf)
```

    ##   Site Species Abundance
    ## 1   35   10017         1
    ## 2   35   10024        18
    ## 3   35   10034         1
    ## 4   35   10035         1
    ## 5   35   10056         2
    ## 6   35   10080         3

``` r
dim(vegedf)
```

    ## [1] 460878      3

``` r
sum(!duplicated(vegedf[,1]))
```

    ## [1] 715

``` r
sum(!duplicated(vegedf[,2]))
```

    ## [1] 3697

[vegemat](https://bioRgeo.github.io/bioregion/reference/vegemat.html) a
co-occurrence `matrix` containing the same information gathered in a
`matrix` with 715 rows and 3,697 columns,

``` r
data(vegemat)
vegemat[1:10,1:10]
```

    ##     Species
    ## Site 10001 10002 10003 10004 10005 10006 10007 10008 10009 10010
    ##   35     0     0     0     0     0     0     0     0     0     0
    ##   36     2     0     0     0     0     0     1    12     0     0
    ##   37     0     0     0     0     0     0     0     0     0     0
    ##   38     0     0     0     0     0     0     0     0     0     0
    ##   39     5     0     0     0     0     0     0     2     0     0
    ##   84     0     0     0     0     0     0     0     0     0     0
    ##   85     3     0     0     0     0     0     1     7     0     0
    ##   86     0     0     0     2     0     0     2    22     0     0
    ##   87    16     0     0     0     0     0     2    54     0     0
    ##   88   228     0     0     0     0     0     0     5     0     0

``` r
dim(vegemat)
```

    ## [1]  715 3697

and [vegesf](https://bioRgeo.github.io/bioregion/reference/vegesf.html)
a spatial object containing the geometry of the 715 sites.

## 2. From matrix to network

The function
[mat_to_net](https://bioRgeo.github.io/bioregion/reference/mat_to_net.html)
transforms a co-occurrence `matrix` such as `vegemat` into a network
represented by a `data.frame` (such as
[vegedf](https://bioRgeo.github.io/bioregion/reference/vegedf.html) in
this case). If `weight = TRUE` a third column is added with the values
contained in the `matrix`.

``` r
net <- mat_to_net(vegemat, weight = TRUE, remove_zeroes = FALSE)
```

In line with the network format, the two first columns are named `Node1`
and `Node2` by default.

``` r
head(net)
```

    ##   Node1 Node2 Weight
    ## 1    35 10001      0
    ## 2    35 10002      0
    ## 3    35 10003      0
    ## 4    35 10004      0
    ## 5    35 10005      0
    ## 6    35 10006      0

``` r
dim(net)
```

    ## [1] 2643355       3

If `remove_zeroes = TRUE` the pairs of nodes with a weight equal to 0
will be removed from the output.

``` r
net <- mat_to_net(vegemat, weight = TRUE, remove_zeroes = TRUE)
```

``` r
head(net)
```

    ##    Node1 Node2 Weight
    ## 17    35 10017      1
    ## 24    35 10024     18
    ## 34    35 10034      1
    ## 35    35 10035      1
    ## 56    35 10056      2
    ## 80    35 10080      3

``` r
dim(net)
```

    ## [1] 460878      3

## 3. From network to matrix

The function
[net_to_mat](https://bioRgeo.github.io/bioregion/reference/net_to_mat.html)
does the opposite. It transforms a network represented by a two- or a
three-columns `data.frame` (such as
[vegedf](https://bioRgeo.github.io/bioregion/reference/vegedf.html))
into a co-occurrence `matrix` (such as
[vegemat](https://bioRgeo.github.io/bioregion/reference/vegemat.html) in
this case).

``` r
mat <- net_to_mat(vegedf, weight = TRUE, squared = FALSE, symmetrical = FALSE, missing_value = 0)
```

``` r
mat[1:5,1:5]
```

    ##    10017 10024 10034 10035 10056
    ## 35     1    18     1     1     2
    ## 36   252    57    72    19    75
    ## 37    66     1    13    23    43
    ## 38    17     1     5    89    27
    ## 39    17    17    34     3     8

``` r
dim(mat)
```

    ## [1]  715 3697

If `squared = TRUE` a squared matrix will be generated, the rownames and
colnames will correspond to the concatenation without duplicates of the
two first columns of the `data.frame`.

``` r
mat <- net_to_mat(vegedf, weight = TRUE, squared = TRUE, symmetrical = FALSE, missing_value = 0)
```

``` r
mat[1:5,1:5]
```

    ##    35 36 37 38 39
    ## 35  0  0  0  0  0
    ## 36  0  0  0  0  0
    ## 37  0  0  0  0  0
    ## 38  0  0  0  0  0
    ## 39  0  0  0  0  0

``` r
dim(mat)
```

    ## [1] 4412 4412

The argument `missing_value` defines the value to assign to the pairs of
nodes not present in the input network. The default value is 0 but any
other numeric value can be used.

``` r
temp <- data.frame(Site=c("35","36","36","38","39"), Species=c("36","35","37","37","39"), Abundance=c(1,2,3,4,0))
net <- rbind(temp,vegedf)
mat <- net_to_mat(net, weight = TRUE, squared = TRUE, symmetrical = FALSE, missing_value = -1)
```

``` r
mat[1:5,1:5]
```

    ##    35 36 38 39 37
    ## 35 -1  1 -1 -1 -1
    ## 36  2 -1 -1 -1  3
    ## 38 -1 -1 -1 -1  4
    ## 39 -1 -1 -1  0 -1
    ## 37 -1 -1 -1 -1 -1

Finally, if `squared = TRUE` it is possible to get a symmetrical matrix
as output (`symmetrical = TRUE`). In this case the resulting squared
matrix will be symmetrical, except for the symmetrical pairs of nodes
already present in the input network (35 \<-\> 36) in the example below.

``` r
mat <- net_to_mat(net, weight = TRUE, squared = TRUE, symmetrical = TRUE, missing_value = 0)
```

``` r
mat[1:5,1:5]
```

    ##    35 36 38 39 37
    ## 35  0  1  0  0  0
    ## 36  2  0  0  0  3
    ## 38  0  0  0  0  4
    ## 39  0  0  0  0  0
    ## 37  0  3  4  0  0
