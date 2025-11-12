# 5.3 Compare bioregionalizations

In this vignette, we aim at comparing the assignment of sites into
different bioregions across multiple bioregionalizations, using the
function
[`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md).

## 1. Data

We use the vegetation dataset that comes with `bioregion`.

``` r
data("vegedf")
data("vegemat")

# Calculation of (dis)similarity matrices
vegedissim <- dissimilarity(vegemat, metric = c("Simpson"))
vegesim <- dissimilarity_to_similarity(vegedissim)
```

## 2. Bioregionalization

We use the same three bioregionalization algorithms as in the
[visualization
vignette](https://biorgeo.github.io/bioregion/articles/a5_visualization.html),
i.e. a non-hierarchical, hierarchical and network bioregionalizations.  
We chose 3 bioregions for the non-hierarchical and hierarchical
bioregionalizations.  

``` r
# Non hierarchical bioregionalization
vege_nhclu_kmeans <- nhclu_kmeans(vegedissim, n_clust = 3, index = "Simpson")
vege_nhclu_kmeans$cluster_info # 3
```

    ##     partition_name n_clust
    ## K_3            K_3       3

``` r
# Hierarchical bioregionalization
set.seed(1)
vege_hclu_hierarclust <- hclu_hierarclust(dissimilarity = vegedissim,
                                          method = "mcquitty", n_clust = 3,
                                          optimal_tree_method = "best")
vege_hclu_hierarclust$cluster_info # 3
```

    ##   partition_name n_clust requested_n_clust output_cut_height
    ## 1            K_3       3                 3             0.625

``` r
# Network bioregionalization
set.seed(1)
vege_netclu_walktrap <- netclu_walktrap(vegesim,
                                        index = names(vegesim)[3])
vege_netclu_walktrap$cluster_info # 3
```

    ##     partition_name n_clust
    ## K_3            K_3       3

## 3. Compare the bioregionalizations

Before comparing the bioregionalizations, we build a common `data.frame`
containing the three distinct bioregionalizations.  

``` r
comp <- dplyr::left_join(vege_hclu_hierarclust$clusters,
                         vege_netclu_walktrap$clusters,
                         by = "ID")
colnames(comp) <- c("ID", "K_3_hclu", "K_3_netclu")
comp <- dplyr::left_join(comp,
                         vege_nhclu_kmeans$clusters,
                         by = "ID")
colnames(comp) <- c("ID", "K_3_hclu", "K_3_netclu", "K_3_nhclu")

head(comp)
```

    ##     ID K_3_hclu K_3_netclu K_3_nhclu
    ## 1  512        1          1         1
    ## 2  799        1          3         1
    ## 3  375        2          1         2
    ## 4  476        2          1         2
    ## 5  971        1          3         1
    ## 6 1282        1          3         1

We can now run the function
[`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md).

``` r
hclu_vs_netclu <- compare_bioregionalizations(
  bioregionalizations = comp[, c("K_3_hclu", "K_3_netclu", "K_3_nhclu")],
  store_pairwise_membership = TRUE,
  cor_frequency = TRUE,
  store_confusion_matrix = TRUE)

str(hclu_vs_netclu)
```

    ## List of 7
    ##  $ args                         :List of 4
    ##   ..$ indices                  : chr [1:2] "rand" "jaccard"
    ##   ..$ cor_frequency            : logi TRUE
    ##   ..$ store_pairwise_membership: logi TRUE
    ##   ..$ store_confusion_matrix   : logi TRUE
    ##  $ inputs                       : Named int [1:2] 715 3
    ##   ..- attr(*, "names")= chr [1:2] "number_items" "number_bioregionalizations"
    ##  $ pairwise_membership          : logi [1:255255, 1:3] TRUE FALSE FALSE TRUE TRUE TRUE ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:255255] "1_2" "1_3" "1_4" "1_5" ...
    ##   .. ..$ : chr [1:3] "K_3_hclu" "K_3_netclu" "K_3_nhclu"
    ##  $ freq_item_pw_membership      : Named num [1:255255] 2 1 1 2 2 2 2 2 2 2 ...
    ##   ..- attr(*, "names")= chr [1:255255] "1_2" "1_3" "1_4" "1_5" ...
    ##  $ bioregionalization_freq_cor  : Named num [1:3] 0.851 0.737 0.893
    ##   ..- attr(*, "names")= chr [1:3] "K_3_hclu" "K_3_netclu" "K_3_nhclu"
    ##  $ confusion_matrix             :List of 3
    ##   ..$ K_3_hclu%K_3_netclu : Named int [1:4] 63316 41892 36800 113247
    ##   .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
    ##   ..$ K_3_hclu%K_3_nhclu  : Named int [1:4] 85543 19665 11248 138799
    ##   .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
    ##   ..$ K_3_netclu%K_3_nhclu: Named int [1:4] 66314 33802 30477 124662
    ##   .. ..- attr(*, "names")= chr [1:4] "a" "b" "c" "d"
    ##  $ bioregionalization_comparison:'data.frame':   3 obs. of  3 variables:
    ##   ..$ bioregionalization_comparison: chr [1:3] "K_3_hclu%K_3_netclu" "K_3_hclu%K_3_nhclu" "K_3_netclu%K_3_nhclu"
    ##   ..$ rand                         : num [1:3] 0.692 0.879 0.748
    ##   ..$ jaccard                      : num [1:3] 0.446 0.735 0.508
    ##  - attr(*, "class")= chr [1:2] "bioregion.bioregionalization.comparison" "list"

[`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md)
produces several outputs which: - look within each bioregionalization
how sites are assigned to bioregions - compare different
bioregionalizations by analysing whether they produce similar pairwise
memberships

Let’s first look at pairwise membership within bioregionalization.

### 3.1 Pairwise membership

The number of pairwise combinations for \\n\\ sites equals \\n(n-1)/2\\.
So in our case, where we have 715 sites, we do end up with 2.55255^{5}
pairwise combinations.

``` r
nrow(hclu_vs_netclu$pairwise_membership) == nrow(comp)*(nrow(comp)-1)/2
```

    ## [1] TRUE

Pairwise memberships look for each pairs of site whether they are
assigned to the same or to a different bioregion. Let’s look at the
sites 1 and 9 across the different bioregionalization:

``` r
comp[c(1, 9), ]
```

    ##    ID K_3_hclu K_3_netclu K_3_nhclu
    ## 1 512        1          1         1
    ## 9 948        1          3         1

We can see that the sites 1 and 9 are classified in the same bioregion
in the first two bioregionalizations, but not in the third one.  
The `$pairwise_membership` output of
[`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md)
shows this as a `TRUE/FALSE` statement.

``` r
hclu_vs_netclu$pairwise_membership[8:10, ]
```

    ##      K_3_hclu K_3_netclu K_3_nhclu
    ## 1_9      TRUE      FALSE      TRUE
    ## 1_10     TRUE      FALSE      TRUE
    ## 1_11     TRUE      FALSE      TRUE

The number of times each pair of sites are clustered together (i.e. the
sum of rows of the table in `$pairwise_membership`) is available in the
`$freq_item_pw_membership` output:

``` r
hclu_vs_netclu$freq_item_pw_membership[c(1, 8)]
```

    ## 1_2 1_9 
    ##   2   2

The sites 1 and 2 were never classified in the same bioregion across the
three bioregionalizations. Sites 1 and 9 were classified in the same
bioregion in two bioregionalizations. If we look at the total
frequencies:

``` r
table(hclu_vs_netclu$freq_item_pw_membership)
```

    ## 
    ##      0      1      2      3 
    ## 111723  41539  45403  56590

we see that the most dominant situation is when sites are never assigned
to the same bioregion.

### 3.2 Confusion matrix

The confusion matrix allows to compare different bioregionalizations by
looking at the similarity of their pairwise memberships. To do so, the
function computes a confusion matrix with four elements:

. \\a\\ number of pairs of sites grouped in bioregionalization 1 and in
bioregionalization 2 . \\b\\ number of pairs of sites grouped in
bioregionalization 1 but not in bioregionalization 2 . \\c\\ number of
pairs of sites not grouped in bioregionalization 1 but grouped in
bioregionalization 2 . \\d\\ number of pairs of sites not grouped in
both bioregionalization 1 & 2

``` r
hclu_vs_netclu$confusion_matrix
```

    ## $`K_3_hclu%K_3_netclu`
    ##      a      b      c      d 
    ##  63316  41892  36800 113247 
    ## 
    ## $`K_3_hclu%K_3_nhclu`
    ##      a      b      c      d 
    ##  85543  19665  11248 138799 
    ## 
    ## $`K_3_netclu%K_3_nhclu`
    ##      a      b      c      d 
    ##  66314  33802  30477 124662

Based on the confusion matrices, we can compute a range of indices to
indicate the agreement among bioregionalizations. As of now, we have
implemented:  
  
*Rand index* \\(a+d)/(a+b+c+d)\\ The Rand index measures agreement among
bioregionalizations by accounting for both the pairs of sites that are
grouped, but also the pairs of sites that are not grouped.  
  
*Jaccard index* \\a/(a+b+c)\\ The Jaccard index measures agreement among
bioregionalizations by only accounting for pairs of sites that are
grouped.

These two metrics are complementary, because the Jaccard index will tell
if bioregionalizations are similar in their clustering structure,
whereas the Rand index will tell if bioregionalizations are similar not
only in the pairs of items clustered together, but also in terms of the
pairs of sites that are not clustered together. For example, take two
bioregionalizations which never group together the same pairs of sites.
Their Jaccard index will be 0, whereas the Rand index can be \> 0 due to
the sites that are not grouped together.

Additional indices can be manually computed by the users on the basis of
the list of confusion matrices.

In some cases, users may be interested in finding which of the
bioregionalizations is most representative of all bioregionalizations To
find it out, we can compare the pairwise membership of each
bioregionalization with the total frequency of pairwise membership
across all bioregionalizations. This correlation can be requested with
`cor_frequency = TRUE`.

``` r
hclu_vs_netclu$bioregionalization_freq_cor
```

    ##   K_3_hclu K_3_netclu  K_3_nhclu 
    ##  0.8507872  0.7365956  0.8934004

Here the third bioregionalization is the most representative of all
bioregionalizations.
