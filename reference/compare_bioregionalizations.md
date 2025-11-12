# Compare cluster memberships among multiple bioregionalizations

This function computes pairwise comparisons for several
bioregionalizations, usually outputs from `netclu_`, `hclu_`, or
`nhclu_` functions. It also provides the confusion matrix from pairwise
comparisons, enabling the user to compute additional comparison metrics.

## Usage

``` r
compare_bioregionalizations(
  bioregionalizations,
  indices = c("rand", "jaccard"),
  cor_frequency = FALSE,
  store_pairwise_membership = TRUE,
  store_confusion_matrix = TRUE,
  verbose = TRUE
)
```

## Arguments

- bioregionalizations:

  A `data.frame` object where each row corresponds to a site, and each
  column to a bioregionalization.

- indices:

  `NULL` or `character`. Indices to compute for the pairwise comparison
  of bioregionalizations. Currently available metrics are `"rand"` and
  `"jaccard"`.

- cor_frequency:

  A `boolean`. If `TRUE`, computes the correlation between each
  bioregionalization and the total frequency of co-membership of items
  across all bioregionalizations. This is useful for identifying which
  bioregionalization(s) is(are) most representative of all computed
  bioregionalizations.

- store_pairwise_membership:

  A `boolean`. If `TRUE`, stores the pairwise membership of items in the
  output object.

- store_confusion_matrix:

  A `boolean`. If `TRUE`, stores the confusion matrices of pairwise
  bioregionalization comparisons in the output object.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` containing 4 to 7 elements:

1.  **args**: A `list` of user-provided arguments.

2.  **inputs**: A `list` containing information on the input
    bioregionalizations, such as the number of items clustered.

3.  **pairwise_membership** (optional): If
    `store_pairwise_membership = TRUE`, a `boolean matrix` where `TRUE`
    indicates two items are in the same cluster, and `FALSE` indicates
    they are not.

4.  **freq_item_pw_membership**: A `numeric vector` containing the
    number of times each item pair is clustered together, corresponding
    to the sum of rows in `pairwise_membership`.

5.  **bioregionalization_freq_cor** (optional): If
    `cor_frequency = TRUE`, a `numeric vector` of correlations between
    individual bioregionalizations and the total frequency of pairwise
    membership.

6.  **confusion_matrix** (optional): If `store_confusion_matrix = TRUE`,
    a `list` of confusion matrices for each pair of bioregionalizations.

7.  **bioregionalization_comparison**: A `data.frame` containing
    comparison results, where the first column indicates the
    bioregionalizations compared, and the remaining columns contain the
    requested `indices`.

## Details

This function operates in two main steps:

1.  Within each bioregionalization, the function compares all pairs of
    items and documents whether they are clustered together (`TRUE`) or
    separately (`FALSE`). For example, if site 1 and site 2 are
    clustered in the same cluster in bioregionalization 1, their
    pairwise membership `site1_site2` will be `TRUE`. This output is
    stored in the `pairwise_membership` slot if
    `store_pairwise_membership = TRUE`.

2.  Across all bioregionalizations, the function compares their pairwise
    memberships to determine similarity. For each pair of
    bioregionalizations, it computes a confusion matrix with the
    following elements:

- `a`: Number of item pairs grouped in both bioregionalizations.

- `b`: Number of item pairs grouped in the first but not in the second
  bioregionalization.

- `c`: Number of item pairs grouped in the second but not in the first
  bioregionalization.

- `d`: Number of item pairs not grouped in either bioregionalization.

The confusion matrix is stored in `confusion_matrix` if
`store_confusion_matrix = TRUE`.

Based on these confusion matrices, various indices can be computed to
measure agreement among bioregionalizations. The currently implemented
indices are:

- **Rand index**: `(a + d) / (a + b + c + d)` Measures agreement by
  considering both grouped and ungrouped item pairs.

- **Jaccard index**: `a / (a + b + c)` Measures agreement based only on
  grouped item pairs.

These indices are complementary: the Jaccard index evaluates clustering
similarity, while the Rand index considers both clustering and
separation. For example, if two bioregionalizations never group the same
pairs, their Jaccard index will be 0, but their Rand index may be \> 0
due to ungrouped pairs.

Users can compute additional indices manually using the list of
confusion matrices.

To identify which bioregionalization is most representative of the
others, the function can compute the correlation between the pairwise
membership of each bioregionalization and the total frequency of
pairwise membership across all bioregionalizations. This is enabled by
setting `cor_frequency = TRUE`.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a5_3_compare_bioregionalizations.html>.

Associated functions:
[bioregionalization_metrics](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
# We here compare three different bioregionalizations
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "Simpson")
bioregion1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

net <- similarity(comat, metric = "Simpson")
bioregion2 <- netclu_greedy(net)
bioregion3 <- netclu_walktrap(net)

# Make one single data.frame with the bioregionalizations to compare
compare_df <- merge(bioregion1$clusters, bioregion2$clusters, by = "ID")
compare_df <- merge(compare_df, bioregion3$clusters, by = "ID")
colnames(compare_df) <- c("Site", "Hclu", "Greedy", "Walktrap")
rownames(compare_df) <- compare_df$Site
compare_df <- compare_df[, c("Hclu", "Greedy", "Walktrap")]

# Running the function
compare_bioregionalizations(compare_df)
#> 2025-11-12 14:36:35.4352 - Computing pairwise membership comparisons for each bioregionalization...
#> 2025-11-12 14:36:35.437203 - Comparing memberships among bioregionalizations...
#> 2025-11-12 14:36:35.437917 - Computing Rand index...
#> 2025-11-12 14:36:35.438319 - Computing Jaccard index...
#> $args
#> $args$indices
#> [1] "rand"    "jaccard"
#> 
#> $args$cor_frequency
#> [1] FALSE
#> 
#> $args$store_pairwise_membership
#> [1] TRUE
#> 
#> $args$store_confusion_matrix
#> [1] TRUE
#> 
#> 
#> $inputs
#>               number_items number_bioregionalizations 
#>                         20                          3 
#> 
#> $pairwise_membership
#>        Hclu Greedy Walktrap
#> 1_2   FALSE   TRUE     TRUE
#> 1_3   FALSE   TRUE     TRUE
#> 1_4   FALSE   TRUE     TRUE
#> 1_5    TRUE   TRUE     TRUE
#> 1_6   FALSE   TRUE     TRUE
#> 1_7   FALSE   TRUE     TRUE
#> 1_8   FALSE   TRUE     TRUE
#> 1_9    TRUE   TRUE     TRUE
#> 1_10  FALSE   TRUE     TRUE
#> 1_11   TRUE   TRUE     TRUE
#> 1_12  FALSE   TRUE     TRUE
#> 1_13  FALSE   TRUE     TRUE
#> 1_14   TRUE   TRUE     TRUE
#> 1_15   TRUE   TRUE     TRUE
#> 1_16  FALSE   TRUE     TRUE
#> 1_17   TRUE   TRUE     TRUE
#> 1_18  FALSE   TRUE     TRUE
#> 1_19  FALSE  FALSE     TRUE
#> 1_20  FALSE   TRUE     TRUE
#> 2_3   FALSE   TRUE     TRUE
#> 2_4    TRUE   TRUE     TRUE
#> 2_5   FALSE   TRUE     TRUE
#> 2_6    TRUE   TRUE     TRUE
#> 2_7    TRUE   TRUE     TRUE
#> 2_8    TRUE   TRUE     TRUE
#> 2_9   FALSE   TRUE     TRUE
#> 2_10  FALSE   TRUE     TRUE
#> 2_11  FALSE   TRUE     TRUE
#> 2_12   TRUE   TRUE     TRUE
#> 2_13  FALSE   TRUE     TRUE
#> 2_14  FALSE   TRUE     TRUE
#> 2_15  FALSE   TRUE     TRUE
#> 2_16  FALSE   TRUE     TRUE
#> 2_17  FALSE   TRUE     TRUE
#> 2_18  FALSE   TRUE     TRUE
#> 2_19   TRUE  FALSE     TRUE
#> 2_20  FALSE   TRUE     TRUE
#> 3_4   FALSE   TRUE     TRUE
#> 3_5   FALSE   TRUE     TRUE
#> 3_6   FALSE   TRUE     TRUE
#> 3_7   FALSE   TRUE     TRUE
#> 3_8   FALSE   TRUE     TRUE
#> 3_9   FALSE   TRUE     TRUE
#> 3_10   TRUE   TRUE     TRUE
#> 3_11  FALSE   TRUE     TRUE
#> 3_12  FALSE   TRUE     TRUE
#> 3_13   TRUE   TRUE     TRUE
#> 3_14  FALSE   TRUE     TRUE
#> 3_15  FALSE   TRUE     TRUE
#> 3_16   TRUE   TRUE     TRUE
#> 3_17  FALSE   TRUE     TRUE
#> 3_18   TRUE   TRUE     TRUE
#> 3_19  FALSE  FALSE     TRUE
#> 3_20   TRUE   TRUE     TRUE
#> 4_5   FALSE   TRUE     TRUE
#> 4_6    TRUE   TRUE     TRUE
#> 4_7    TRUE   TRUE     TRUE
#> 4_8    TRUE   TRUE     TRUE
#> 4_9   FALSE   TRUE     TRUE
#> 4_10  FALSE   TRUE     TRUE
#> 4_11  FALSE   TRUE     TRUE
#> 4_12   TRUE   TRUE     TRUE
#> 4_13  FALSE   TRUE     TRUE
#> 4_14  FALSE   TRUE     TRUE
#> 4_15  FALSE   TRUE     TRUE
#> 4_16  FALSE   TRUE     TRUE
#> 4_17  FALSE   TRUE     TRUE
#> 4_18  FALSE   TRUE     TRUE
#> 4_19   TRUE  FALSE     TRUE
#> 4_20  FALSE   TRUE     TRUE
#> 5_6   FALSE   TRUE     TRUE
#> 5_7   FALSE   TRUE     TRUE
#> 5_8   FALSE   TRUE     TRUE
#> 5_9    TRUE   TRUE     TRUE
#> 5_10  FALSE   TRUE     TRUE
#> 5_11   TRUE   TRUE     TRUE
#> 5_12  FALSE   TRUE     TRUE
#> 5_13  FALSE   TRUE     TRUE
#> 5_14   TRUE   TRUE     TRUE
#> 5_15   TRUE   TRUE     TRUE
#> 5_16  FALSE   TRUE     TRUE
#> 5_17   TRUE   TRUE     TRUE
#> 5_18  FALSE   TRUE     TRUE
#> 5_19  FALSE  FALSE     TRUE
#> 5_20  FALSE   TRUE     TRUE
#> 6_7    TRUE   TRUE     TRUE
#> 6_8    TRUE   TRUE     TRUE
#> 6_9   FALSE   TRUE     TRUE
#> 6_10  FALSE   TRUE     TRUE
#> 6_11  FALSE   TRUE     TRUE
#> 6_12   TRUE   TRUE     TRUE
#> 6_13  FALSE   TRUE     TRUE
#> 6_14  FALSE   TRUE     TRUE
#> 6_15  FALSE   TRUE     TRUE
#> 6_16  FALSE   TRUE     TRUE
#> 6_17  FALSE   TRUE     TRUE
#> 6_18  FALSE   TRUE     TRUE
#> 6_19   TRUE  FALSE     TRUE
#> 6_20  FALSE   TRUE     TRUE
#> 7_8    TRUE   TRUE     TRUE
#> 7_9   FALSE   TRUE     TRUE
#> 7_10  FALSE   TRUE     TRUE
#> 7_11  FALSE   TRUE     TRUE
#> 7_12   TRUE   TRUE     TRUE
#> 7_13  FALSE   TRUE     TRUE
#> 7_14  FALSE   TRUE     TRUE
#> 7_15  FALSE   TRUE     TRUE
#> 7_16  FALSE   TRUE     TRUE
#> 7_17  FALSE   TRUE     TRUE
#> 7_18  FALSE   TRUE     TRUE
#> 7_19   TRUE  FALSE     TRUE
#> 7_20  FALSE   TRUE     TRUE
#> 8_9   FALSE   TRUE     TRUE
#> 8_10  FALSE   TRUE     TRUE
#> 8_11  FALSE   TRUE     TRUE
#> 8_12   TRUE   TRUE     TRUE
#> 8_13  FALSE   TRUE     TRUE
#> 8_14  FALSE   TRUE     TRUE
#> 8_15  FALSE   TRUE     TRUE
#> 8_16  FALSE   TRUE     TRUE
#> 8_17  FALSE   TRUE     TRUE
#> 8_18  FALSE   TRUE     TRUE
#> 8_19   TRUE  FALSE     TRUE
#> 8_20  FALSE   TRUE     TRUE
#> 9_10  FALSE   TRUE     TRUE
#> 9_11   TRUE   TRUE     TRUE
#> 9_12  FALSE   TRUE     TRUE
#> 9_13  FALSE   TRUE     TRUE
#> 9_14   TRUE   TRUE     TRUE
#> 9_15   TRUE   TRUE     TRUE
#> 9_16  FALSE   TRUE     TRUE
#> 9_17   TRUE   TRUE     TRUE
#> 9_18  FALSE   TRUE     TRUE
#> 9_19  FALSE  FALSE     TRUE
#> 9_20  FALSE   TRUE     TRUE
#> 10_11 FALSE   TRUE     TRUE
#> 10_12 FALSE   TRUE     TRUE
#> 10_13  TRUE   TRUE     TRUE
#> 10_14 FALSE   TRUE     TRUE
#> 10_15 FALSE   TRUE     TRUE
#> 10_16  TRUE   TRUE     TRUE
#> 10_17 FALSE   TRUE     TRUE
#> 10_18  TRUE   TRUE     TRUE
#> 10_19 FALSE  FALSE     TRUE
#> 10_20  TRUE   TRUE     TRUE
#> 11_12 FALSE   TRUE     TRUE
#> 11_13 FALSE   TRUE     TRUE
#> 11_14  TRUE   TRUE     TRUE
#> 11_15  TRUE   TRUE     TRUE
#> 11_16 FALSE   TRUE     TRUE
#> 11_17  TRUE   TRUE     TRUE
#> 11_18 FALSE   TRUE     TRUE
#> 11_19 FALSE  FALSE     TRUE
#> 11_20 FALSE   TRUE     TRUE
#> 12_13 FALSE   TRUE     TRUE
#> 12_14 FALSE   TRUE     TRUE
#> 12_15 FALSE   TRUE     TRUE
#> 12_16 FALSE   TRUE     TRUE
#> 12_17 FALSE   TRUE     TRUE
#> 12_18 FALSE   TRUE     TRUE
#> 12_19  TRUE  FALSE     TRUE
#> 12_20 FALSE   TRUE     TRUE
#> 13_14 FALSE   TRUE     TRUE
#> 13_15 FALSE   TRUE     TRUE
#> 13_16  TRUE   TRUE     TRUE
#> 13_17 FALSE   TRUE     TRUE
#> 13_18  TRUE   TRUE     TRUE
#> 13_19 FALSE  FALSE     TRUE
#> 13_20  TRUE   TRUE     TRUE
#> 14_15  TRUE   TRUE     TRUE
#> 14_16 FALSE   TRUE     TRUE
#> 14_17  TRUE   TRUE     TRUE
#> 14_18 FALSE   TRUE     TRUE
#> 14_19 FALSE  FALSE     TRUE
#> 14_20 FALSE   TRUE     TRUE
#> 15_16 FALSE   TRUE     TRUE
#> 15_17  TRUE   TRUE     TRUE
#> 15_18 FALSE   TRUE     TRUE
#> 15_19 FALSE  FALSE     TRUE
#> 15_20 FALSE   TRUE     TRUE
#> 16_17 FALSE   TRUE     TRUE
#> 16_18  TRUE   TRUE     TRUE
#> 16_19 FALSE  FALSE     TRUE
#> 16_20  TRUE   TRUE     TRUE
#> 17_18 FALSE   TRUE     TRUE
#> 17_19 FALSE  FALSE     TRUE
#> 17_20 FALSE   TRUE     TRUE
#> 18_19 FALSE  FALSE     TRUE
#> 18_20  TRUE   TRUE     TRUE
#> 19_20 FALSE  FALSE     TRUE
#> 
#> $freq_item_pw_membership
#>   1_2   1_3   1_4   1_5   1_6   1_7   1_8   1_9  1_10  1_11  1_12  1_13  1_14 
#>     2     2     2     3     2     2     2     3     2     3     2     2     3 
#>  1_15  1_16  1_17  1_18  1_19  1_20   2_3   2_4   2_5   2_6   2_7   2_8   2_9 
#>     3     2     3     2     1     2     2     3     2     3     3     3     2 
#>  2_10  2_11  2_12  2_13  2_14  2_15  2_16  2_17  2_18  2_19  2_20   3_4   3_5 
#>     2     2     3     2     2     2     2     2     2     2     2     2     2 
#>   3_6   3_7   3_8   3_9  3_10  3_11  3_12  3_13  3_14  3_15  3_16  3_17  3_18 
#>     2     2     2     2     3     2     2     3     2     2     3     2     3 
#>  3_19  3_20   4_5   4_6   4_7   4_8   4_9  4_10  4_11  4_12  4_13  4_14  4_15 
#>     1     3     2     3     3     3     2     2     2     3     2     2     2 
#>  4_16  4_17  4_18  4_19  4_20   5_6   5_7   5_8   5_9  5_10  5_11  5_12  5_13 
#>     2     2     2     2     2     2     2     2     3     2     3     2     2 
#>  5_14  5_15  5_16  5_17  5_18  5_19  5_20   6_7   6_8   6_9  6_10  6_11  6_12 
#>     3     3     2     3     2     1     2     3     3     2     2     2     3 
#>  6_13  6_14  6_15  6_16  6_17  6_18  6_19  6_20   7_8   7_9  7_10  7_11  7_12 
#>     2     2     2     2     2     2     2     2     3     2     2     2     3 
#>  7_13  7_14  7_15  7_16  7_17  7_18  7_19  7_20   8_9  8_10  8_11  8_12  8_13 
#>     2     2     2     2     2     2     2     2     2     2     2     3     2 
#>  8_14  8_15  8_16  8_17  8_18  8_19  8_20  9_10  9_11  9_12  9_13  9_14  9_15 
#>     2     2     2     2     2     2     2     2     3     2     2     3     3 
#>  9_16  9_17  9_18  9_19  9_20 10_11 10_12 10_13 10_14 10_15 10_16 10_17 10_18 
#>     2     3     2     1     2     2     2     3     2     2     3     2     3 
#> 10_19 10_20 11_12 11_13 11_14 11_15 11_16 11_17 11_18 11_19 11_20 12_13 12_14 
#>     1     3     2     2     3     3     2     3     2     1     2     2     2 
#> 12_15 12_16 12_17 12_18 12_19 12_20 13_14 13_15 13_16 13_17 13_18 13_19 13_20 
#>     2     2     2     2     2     2     2     2     3     2     3     1     3 
#> 14_15 14_16 14_17 14_18 14_19 14_20 15_16 15_17 15_18 15_19 15_20 16_17 16_18 
#>     3     2     3     2     1     2     2     3     2     1     2     2     3 
#> 16_19 16_20 17_18 17_19 17_20 18_19 18_20 19_20 
#>     1     3     2     1     2     1     3     1 
#> 
#> $confusion_matrix
#> $confusion_matrix$`Hclu%Greedy`
#>   a   b   c   d 
#>  51   6 120  13 
#> 
#> $confusion_matrix$`Hclu%Walktrap`
#>   a   b   c   d 
#>  57   0 133   0 
#> 
#> $confusion_matrix$`Greedy%Walktrap`
#>   a   b   c   d 
#> 171   0  19   0 
#> 
#> 
#> $bioregionalization_comparison
#>   bioregionalization_comparison      rand   jaccard
#> 1                   Hclu%Greedy 0.3368421 0.2881356
#> 2                 Hclu%Walktrap 0.3000000 0.3000000
#> 3               Greedy%Walktrap 0.9000000 0.9000000
#> 
#> attr(,"class")
#> [1] "bioregion.bioregionalization.comparison"
#> [2] "list"                                   

# Find out which bioregionalizations are most representative
compare_bioregionalizations(compare_df,
                            cor_frequency = TRUE)
#> 2025-11-12 14:36:35.441756 - Computing pairwise membership comparisons for each bioregionalization...
#> 2025-11-12 14:36:35.443279 - Comparing memberships among bioregionalizations...
#> 2025-11-12 14:36:35.443913 - Computing Rand index...
#> 2025-11-12 14:36:35.44426 - Computing Jaccard index...
#> 2025-11-12 14:36:35.444576 - Computing the correlation between each bioregionalization and the vector of frequency of pairwise membership...
#> $args
#> $args$indices
#> [1] "rand"    "jaccard"
#> 
#> $args$cor_frequency
#> [1] TRUE
#> 
#> $args$store_pairwise_membership
#> [1] TRUE
#> 
#> $args$store_confusion_matrix
#> [1] TRUE
#> 
#> 
#> $inputs
#>               number_items number_bioregionalizations 
#>                         20                          3 
#> 
#> $pairwise_membership
#>        Hclu Greedy Walktrap
#> 1_2   FALSE   TRUE     TRUE
#> 1_3   FALSE   TRUE     TRUE
#> 1_4   FALSE   TRUE     TRUE
#> 1_5    TRUE   TRUE     TRUE
#> 1_6   FALSE   TRUE     TRUE
#> 1_7   FALSE   TRUE     TRUE
#> 1_8   FALSE   TRUE     TRUE
#> 1_9    TRUE   TRUE     TRUE
#> 1_10  FALSE   TRUE     TRUE
#> 1_11   TRUE   TRUE     TRUE
#> 1_12  FALSE   TRUE     TRUE
#> 1_13  FALSE   TRUE     TRUE
#> 1_14   TRUE   TRUE     TRUE
#> 1_15   TRUE   TRUE     TRUE
#> 1_16  FALSE   TRUE     TRUE
#> 1_17   TRUE   TRUE     TRUE
#> 1_18  FALSE   TRUE     TRUE
#> 1_19  FALSE  FALSE     TRUE
#> 1_20  FALSE   TRUE     TRUE
#> 2_3   FALSE   TRUE     TRUE
#> 2_4    TRUE   TRUE     TRUE
#> 2_5   FALSE   TRUE     TRUE
#> 2_6    TRUE   TRUE     TRUE
#> 2_7    TRUE   TRUE     TRUE
#> 2_8    TRUE   TRUE     TRUE
#> 2_9   FALSE   TRUE     TRUE
#> 2_10  FALSE   TRUE     TRUE
#> 2_11  FALSE   TRUE     TRUE
#> 2_12   TRUE   TRUE     TRUE
#> 2_13  FALSE   TRUE     TRUE
#> 2_14  FALSE   TRUE     TRUE
#> 2_15  FALSE   TRUE     TRUE
#> 2_16  FALSE   TRUE     TRUE
#> 2_17  FALSE   TRUE     TRUE
#> 2_18  FALSE   TRUE     TRUE
#> 2_19   TRUE  FALSE     TRUE
#> 2_20  FALSE   TRUE     TRUE
#> 3_4   FALSE   TRUE     TRUE
#> 3_5   FALSE   TRUE     TRUE
#> 3_6   FALSE   TRUE     TRUE
#> 3_7   FALSE   TRUE     TRUE
#> 3_8   FALSE   TRUE     TRUE
#> 3_9   FALSE   TRUE     TRUE
#> 3_10   TRUE   TRUE     TRUE
#> 3_11  FALSE   TRUE     TRUE
#> 3_12  FALSE   TRUE     TRUE
#> 3_13   TRUE   TRUE     TRUE
#> 3_14  FALSE   TRUE     TRUE
#> 3_15  FALSE   TRUE     TRUE
#> 3_16   TRUE   TRUE     TRUE
#> 3_17  FALSE   TRUE     TRUE
#> 3_18   TRUE   TRUE     TRUE
#> 3_19  FALSE  FALSE     TRUE
#> 3_20   TRUE   TRUE     TRUE
#> 4_5   FALSE   TRUE     TRUE
#> 4_6    TRUE   TRUE     TRUE
#> 4_7    TRUE   TRUE     TRUE
#> 4_8    TRUE   TRUE     TRUE
#> 4_9   FALSE   TRUE     TRUE
#> 4_10  FALSE   TRUE     TRUE
#> 4_11  FALSE   TRUE     TRUE
#> 4_12   TRUE   TRUE     TRUE
#> 4_13  FALSE   TRUE     TRUE
#> 4_14  FALSE   TRUE     TRUE
#> 4_15  FALSE   TRUE     TRUE
#> 4_16  FALSE   TRUE     TRUE
#> 4_17  FALSE   TRUE     TRUE
#> 4_18  FALSE   TRUE     TRUE
#> 4_19   TRUE  FALSE     TRUE
#> 4_20  FALSE   TRUE     TRUE
#> 5_6   FALSE   TRUE     TRUE
#> 5_7   FALSE   TRUE     TRUE
#> 5_8   FALSE   TRUE     TRUE
#> 5_9    TRUE   TRUE     TRUE
#> 5_10  FALSE   TRUE     TRUE
#> 5_11   TRUE   TRUE     TRUE
#> 5_12  FALSE   TRUE     TRUE
#> 5_13  FALSE   TRUE     TRUE
#> 5_14   TRUE   TRUE     TRUE
#> 5_15   TRUE   TRUE     TRUE
#> 5_16  FALSE   TRUE     TRUE
#> 5_17   TRUE   TRUE     TRUE
#> 5_18  FALSE   TRUE     TRUE
#> 5_19  FALSE  FALSE     TRUE
#> 5_20  FALSE   TRUE     TRUE
#> 6_7    TRUE   TRUE     TRUE
#> 6_8    TRUE   TRUE     TRUE
#> 6_9   FALSE   TRUE     TRUE
#> 6_10  FALSE   TRUE     TRUE
#> 6_11  FALSE   TRUE     TRUE
#> 6_12   TRUE   TRUE     TRUE
#> 6_13  FALSE   TRUE     TRUE
#> 6_14  FALSE   TRUE     TRUE
#> 6_15  FALSE   TRUE     TRUE
#> 6_16  FALSE   TRUE     TRUE
#> 6_17  FALSE   TRUE     TRUE
#> 6_18  FALSE   TRUE     TRUE
#> 6_19   TRUE  FALSE     TRUE
#> 6_20  FALSE   TRUE     TRUE
#> 7_8    TRUE   TRUE     TRUE
#> 7_9   FALSE   TRUE     TRUE
#> 7_10  FALSE   TRUE     TRUE
#> 7_11  FALSE   TRUE     TRUE
#> 7_12   TRUE   TRUE     TRUE
#> 7_13  FALSE   TRUE     TRUE
#> 7_14  FALSE   TRUE     TRUE
#> 7_15  FALSE   TRUE     TRUE
#> 7_16  FALSE   TRUE     TRUE
#> 7_17  FALSE   TRUE     TRUE
#> 7_18  FALSE   TRUE     TRUE
#> 7_19   TRUE  FALSE     TRUE
#> 7_20  FALSE   TRUE     TRUE
#> 8_9   FALSE   TRUE     TRUE
#> 8_10  FALSE   TRUE     TRUE
#> 8_11  FALSE   TRUE     TRUE
#> 8_12   TRUE   TRUE     TRUE
#> 8_13  FALSE   TRUE     TRUE
#> 8_14  FALSE   TRUE     TRUE
#> 8_15  FALSE   TRUE     TRUE
#> 8_16  FALSE   TRUE     TRUE
#> 8_17  FALSE   TRUE     TRUE
#> 8_18  FALSE   TRUE     TRUE
#> 8_19   TRUE  FALSE     TRUE
#> 8_20  FALSE   TRUE     TRUE
#> 9_10  FALSE   TRUE     TRUE
#> 9_11   TRUE   TRUE     TRUE
#> 9_12  FALSE   TRUE     TRUE
#> 9_13  FALSE   TRUE     TRUE
#> 9_14   TRUE   TRUE     TRUE
#> 9_15   TRUE   TRUE     TRUE
#> 9_16  FALSE   TRUE     TRUE
#> 9_17   TRUE   TRUE     TRUE
#> 9_18  FALSE   TRUE     TRUE
#> 9_19  FALSE  FALSE     TRUE
#> 9_20  FALSE   TRUE     TRUE
#> 10_11 FALSE   TRUE     TRUE
#> 10_12 FALSE   TRUE     TRUE
#> 10_13  TRUE   TRUE     TRUE
#> 10_14 FALSE   TRUE     TRUE
#> 10_15 FALSE   TRUE     TRUE
#> 10_16  TRUE   TRUE     TRUE
#> 10_17 FALSE   TRUE     TRUE
#> 10_18  TRUE   TRUE     TRUE
#> 10_19 FALSE  FALSE     TRUE
#> 10_20  TRUE   TRUE     TRUE
#> 11_12 FALSE   TRUE     TRUE
#> 11_13 FALSE   TRUE     TRUE
#> 11_14  TRUE   TRUE     TRUE
#> 11_15  TRUE   TRUE     TRUE
#> 11_16 FALSE   TRUE     TRUE
#> 11_17  TRUE   TRUE     TRUE
#> 11_18 FALSE   TRUE     TRUE
#> 11_19 FALSE  FALSE     TRUE
#> 11_20 FALSE   TRUE     TRUE
#> 12_13 FALSE   TRUE     TRUE
#> 12_14 FALSE   TRUE     TRUE
#> 12_15 FALSE   TRUE     TRUE
#> 12_16 FALSE   TRUE     TRUE
#> 12_17 FALSE   TRUE     TRUE
#> 12_18 FALSE   TRUE     TRUE
#> 12_19  TRUE  FALSE     TRUE
#> 12_20 FALSE   TRUE     TRUE
#> 13_14 FALSE   TRUE     TRUE
#> 13_15 FALSE   TRUE     TRUE
#> 13_16  TRUE   TRUE     TRUE
#> 13_17 FALSE   TRUE     TRUE
#> 13_18  TRUE   TRUE     TRUE
#> 13_19 FALSE  FALSE     TRUE
#> 13_20  TRUE   TRUE     TRUE
#> 14_15  TRUE   TRUE     TRUE
#> 14_16 FALSE   TRUE     TRUE
#> 14_17  TRUE   TRUE     TRUE
#> 14_18 FALSE   TRUE     TRUE
#> 14_19 FALSE  FALSE     TRUE
#> 14_20 FALSE   TRUE     TRUE
#> 15_16 FALSE   TRUE     TRUE
#> 15_17  TRUE   TRUE     TRUE
#> 15_18 FALSE   TRUE     TRUE
#> 15_19 FALSE  FALSE     TRUE
#> 15_20 FALSE   TRUE     TRUE
#> 16_17 FALSE   TRUE     TRUE
#> 16_18  TRUE   TRUE     TRUE
#> 16_19 FALSE  FALSE     TRUE
#> 16_20  TRUE   TRUE     TRUE
#> 17_18 FALSE   TRUE     TRUE
#> 17_19 FALSE  FALSE     TRUE
#> 17_20 FALSE   TRUE     TRUE
#> 18_19 FALSE  FALSE     TRUE
#> 18_20  TRUE   TRUE     TRUE
#> 19_20 FALSE  FALSE     TRUE
#> 
#> $freq_item_pw_membership
#>   1_2   1_3   1_4   1_5   1_6   1_7   1_8   1_9  1_10  1_11  1_12  1_13  1_14 
#>     2     2     2     3     2     2     2     3     2     3     2     2     3 
#>  1_15  1_16  1_17  1_18  1_19  1_20   2_3   2_4   2_5   2_6   2_7   2_8   2_9 
#>     3     2     3     2     1     2     2     3     2     3     3     3     2 
#>  2_10  2_11  2_12  2_13  2_14  2_15  2_16  2_17  2_18  2_19  2_20   3_4   3_5 
#>     2     2     3     2     2     2     2     2     2     2     2     2     2 
#>   3_6   3_7   3_8   3_9  3_10  3_11  3_12  3_13  3_14  3_15  3_16  3_17  3_18 
#>     2     2     2     2     3     2     2     3     2     2     3     2     3 
#>  3_19  3_20   4_5   4_6   4_7   4_8   4_9  4_10  4_11  4_12  4_13  4_14  4_15 
#>     1     3     2     3     3     3     2     2     2     3     2     2     2 
#>  4_16  4_17  4_18  4_19  4_20   5_6   5_7   5_8   5_9  5_10  5_11  5_12  5_13 
#>     2     2     2     2     2     2     2     2     3     2     3     2     2 
#>  5_14  5_15  5_16  5_17  5_18  5_19  5_20   6_7   6_8   6_9  6_10  6_11  6_12 
#>     3     3     2     3     2     1     2     3     3     2     2     2     3 
#>  6_13  6_14  6_15  6_16  6_17  6_18  6_19  6_20   7_8   7_9  7_10  7_11  7_12 
#>     2     2     2     2     2     2     2     2     3     2     2     2     3 
#>  7_13  7_14  7_15  7_16  7_17  7_18  7_19  7_20   8_9  8_10  8_11  8_12  8_13 
#>     2     2     2     2     2     2     2     2     2     2     2     3     2 
#>  8_14  8_15  8_16  8_17  8_18  8_19  8_20  9_10  9_11  9_12  9_13  9_14  9_15 
#>     2     2     2     2     2     2     2     2     3     2     2     3     3 
#>  9_16  9_17  9_18  9_19  9_20 10_11 10_12 10_13 10_14 10_15 10_16 10_17 10_18 
#>     2     3     2     1     2     2     2     3     2     2     3     2     3 
#> 10_19 10_20 11_12 11_13 11_14 11_15 11_16 11_17 11_18 11_19 11_20 12_13 12_14 
#>     1     3     2     2     3     3     2     3     2     1     2     2     2 
#> 12_15 12_16 12_17 12_18 12_19 12_20 13_14 13_15 13_16 13_17 13_18 13_19 13_20 
#>     2     2     2     2     2     2     2     2     3     2     3     1     3 
#> 14_15 14_16 14_17 14_18 14_19 14_20 15_16 15_17 15_18 15_19 15_20 16_17 16_18 
#>     3     2     3     2     1     2     2     3     2     1     2     2     3 
#> 16_19 16_20 17_18 17_19 17_20 18_19 18_20 19_20 
#>     1     3     2     1     2     1     3     1 
#> 
#> $bioregionalization_freq_cor
#>      Hclu    Greedy  Walktrap 
#> 0.8347745 0.5409681 0.0000000 
#> 
#> $confusion_matrix
#> $confusion_matrix$`Hclu%Greedy`
#>   a   b   c   d 
#>  51   6 120  13 
#> 
#> $confusion_matrix$`Hclu%Walktrap`
#>   a   b   c   d 
#>  57   0 133   0 
#> 
#> $confusion_matrix$`Greedy%Walktrap`
#>   a   b   c   d 
#> 171   0  19   0 
#> 
#> 
#> $bioregionalization_comparison
#>   bioregionalization_comparison      rand   jaccard
#> 1                   Hclu%Greedy 0.3368421 0.2881356
#> 2                 Hclu%Walktrap 0.3000000 0.3000000
#> 3               Greedy%Walktrap 0.9000000 0.9000000
#> 
#> attr(,"class")
#> [1] "bioregion.bioregionalization.comparison"
#> [2] "list"                                   
                                
```
