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
#> 2026-10-24 07:09:02.178469 - Computing pairwise membership comparisons for each bioregionalization...
#> 2026-10-24 07:09:02.179539 - Comparing memberships among bioregionalizations...
#> 2026-10-24 07:09:02.180448 - Computing Rand index...
#> 2026-10-24 07:09:02.180882 - Computing Jaccard index...
#> Bioregionalization comparison:
#>  - 3 bioregionalizations compared
#>  - 20 items in the clustering
#>  - Requested indices:  rand jaccard 
#>  - Metric summary:
#>           rand   jaccard
#> Min  0.3421053 0.3421053
#> Mean 0.5473684 0.5302256
#> Max  0.9000000 0.9000000
#>  - Item pairwise membership  stored in outputs
#>  - Confusion matrices of bioregionalization comparisons  stored in outputs

# Find out which bioregionalizations are most representative
compare_bioregionalizations(compare_df,
                            cor_frequency = TRUE)
#> 2026-10-24 07:09:02.182863 - Computing pairwise membership comparisons for each bioregionalization...
#> 2026-10-24 07:09:02.183555 - Comparing memberships among bioregionalizations...
#> 2026-10-24 07:09:02.184237 - Computing Rand index...
#> 2026-10-24 07:09:02.184581 - Computing Jaccard index...
#> 2026-10-24 07:09:02.184901 - Computing the correlation between each bioregionalization and the vector of frequency of pairwise membership...
#> Bioregionalization comparison:
#>  - 3 bioregionalizations compared
#>  - 20 items in the clustering
#>  - Requested indices:  rand jaccard 
#>  - Metric summary:
#>           rand   jaccard
#> Min  0.3421053 0.3421053
#> Mean 0.5473684 0.5302256
#> Max  0.9000000 0.9000000
#>  - Correlation between each bioregionalization and the total frequency of item pairwise membership computed:
#>    # Range:  0  -  0.859 
#>    # Bioregionalization(s) most representative (i.e., highest correlation): 
#>  Hclu 
#>  Correlation =  0.859 
#>  - Item pairwise membership  stored in outputs
#>  - Confusion matrices of bioregionalization comparisons  stored in outputs
  
```
