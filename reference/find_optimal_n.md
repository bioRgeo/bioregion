# Search for an optimal number of clusters in a list of bioregionalizations

This function aims to optimize one or several criteria on a set of
ordered bioregionalizations. It is typically used to find one or more
optimal cluster counts on hierarchical trees to cut or ranges of
bioregionalizations from k-means or PAM. Users should exercise caution
in other cases (e.g., unordered bioregionalizations or unrelated
bioregionalizations).

## Usage

``` r
find_optimal_n(
  bioregionalizations,
  metrics_to_use = "all",
  criterion = "elbow",
  step_quantile = 0.99,
  step_levels = NULL,
  step_round_above = TRUE,
  metric_cutoffs = c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999),
  n_breakpoints = 1,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- bioregionalizations:

  A `bioregion.bioregionalization.metrics` object (output from
  [`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md))
  or a `data.frame` with the first two columns named `K`
  (bioregionalization name) and `n_clusters` (number of clusters),
  followed by columns with numeric evaluation metrics.

- metrics_to_use:

  A `character` vector or single string specifying metrics in
  `bioregionalizations` for calculating optimal clusters. Defaults to
  `"all"` (uses all metrics).

- criterion:

  A `character` string specifying the criterion to identify optimal
  clusters. Options include `"elbow"`, `"increasing_step"`,
  `"decreasing_step"`, `"cutoff"`, `"breakpoints"`, `"min"`, or `"max"`.
  Defaults to `"elbow"`. See Details.

- step_quantile:

  For `"increasing_step"` or `"decreasing_step"`, specifies the quantile
  of differences between consecutive bioregionalizations as the cutoff
  to identify significant steps in `eval_metric`.

- step_levels:

  For `"increasing_step"` or `"decreasing_step"`, specifies the number
  of largest steps to retain as cutoffs.

- step_round_above:

  A `boolean` indicating whether the optimal clusters are above (`TRUE`)
  or below (`FALSE`) identified steps. Defaults to `TRUE`.

- metric_cutoffs:

  For `criterion = "cutoff"`, specifies the cutoffs of `eval_metric` to
  extract cluster counts.

- n_breakpoints:

  Specifies the number of breakpoints to find in the curve. Defaults to
  1.

- plot:

  A `boolean` indicating if a plot of the first `eval_metric` with
  identified optimal clusters should be drawn.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

A `list` of class `bioregion.optimal.n` with these elements:

- `args`: Input arguments.

- `evaluation_df`: The input evaluation `data.frame`, appended with
  `boolean` columns for optimal cluster counts.

- `optimal_nb_clusters`: A `list` with optimal cluster counts for each
  metric in `"metrics_to_use"`, based on the chosen `criterion`.

- `plot`: The plot (if requested).

## Details

This function explores evaluation metric ~ cluster relationships,
applying criteria to find optimal cluster counts.

**Note on criteria:** Several criteria can return multiple optimal
cluster counts, emphasizing hierarchical or nested bioregionalizations.
This approach aligns with modern recommendations for biological
datasets, as seen in Ficetola et al. (2017)'s reanalysis of Holt et al.
(2013).

**Criteria for optimal clusters:**

- `elbow`: Identifies the "elbow" point in the evaluation metric curve,
  where incremental improvements diminish. Based on a method to find the
  maximum distance from a straight line linking curve endpoints.

- `increasing_step` or `decreasing_step`: Highlights significant
  increases or decreases in metrics by analyzing pairwise differences
  between bioregionalizations. Users specify `step_quantile` or
  `step_levels`.

- `cutoffs`: Derives clusters from specified metric cutoffs, e.g., as in
  Holt et al. (2013). Adjust cutoffs based on spatial scale.

- `breakpoints`: Uses segmented regression to find breakpoints. Requires
  specifying `n_breakpoints`.

- `min` & `max`: Selects clusters at minimum or maximum metric values.

## Note

Please note that finding the optimal number of clusters is a procedure
which normally requires decisions from the users, and as such can hardly
be fully automatized. Users are strongly advised to read the references
indicated below to look for guidance on how to choose their optimal
number(s) of clusters. Consider the "optimal" numbers of clusters
returned by this function as first approximation of the best numbers for
your bioregionalization.

## References

Holt BG, Lessard J, Borregaard MK, Fritz SA, Araújo MB, Dimitrov D,
Fabre P, Graham CH, Graves GR, Jønsson Ka, Nogués-Bravo D, Wang Z,
Whittaker RJ, Fjeldså J & Rahbek C (2013) An update of Wallace's
zoogeographic regions of the world. *Science* 339, 74-78.

Ficetola GF, Mazel F & Thuiller W (2017) Global determinants of
zoogeographical boundaries. *Nature Ecology & Evolution* 1, 0089.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html#optimaln>.

Associated functions:
[hclu_hierarclust](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)

## Author

Boris Leroy (<leroy.boris@gmail.com>)  
Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

## Examples

``` r
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "all")

# User-defined number of clusters
tree <- hclu_hierarclust(dissim,
                          optimal_tree_method = "best",
                          n_clust = 5:10)
#> Randomizing the dissimilarity matrix with 100 trials
#>  -- range of cophenetic correlation coefficients among trials: 0.787 - 0.7896
#> 
#> Final tree has a 0.7896 cophenetic correlation coefficient with the initial dissimilarity matrix
#> Determining the cut height to reach 5 groups...
#> --> 0.25
#> Determining the cut height to reach 6 groups...
#> --> 0.234375
#> Determining the cut height to reach 7 groups...
#> --> 0.2265625
#> Determining the cut height to reach 8 groups...
#> --> 0.21875
#> Determining the cut height to reach 9 groups...
#> --> 0.203125
#> Determining the cut height to reach 10 groups...
#> --> 0.1875
tree
#> Clustering results for algorithm : hclu_hierarclust 
#>  (hierarchical clustering based on a dissimilarity matrix)
#>  - Number of sites:  20 
#>  - Name of dissimilarity metric:  Jaccard 
#>  - Tree construction method:  average 
#>  - Randomization of the dissimilarity matrix:  yes, number of trials 100 
#>  - Method to compute the final tree:  Tree with the best cophenetic correlation coefficient 
#>  - Cophenetic correlation coefficient:  0.79 
#>  - Number of clusters requested by the user:  5 
#> Clustering results:
#>  - Number of partitions:  6 
#>  - Partitions are hierarchical
#>  - Number of clusters:  5 6 7 8 9 10 
#>  - Height of cut of the hierarchical tree: 0.25 0.234 0.227 0.219 0.203 0.188 

a <- bioregionalization_metrics(tree,
                                dissimilarity = dissim,
                                species_col = "Node2",
                                site_col = "Node1",
                                eval_metric = "anosim")
#> Computing similarity-based metrics...
#>   - anosim OK
                                   
find_optimal_n(a, criterion = 'increasing_step', plot = FALSE)
#> Number of bioregionalizations: 6
#> ...Caveat: be cautious with the interpretation of metric analyses with such a low number of bioregionalizations
#> Searching for potential optimal number(s) of clusters based on the increasing_step method
#>  - Step method
#> Search for an optimal number of clusters:
#>  - 6  partition(s) evaluated
#>  - Range of clusters explored: from  5  to  10 
#>  - Evaluated metric(s):  anosim 
#> 
#> Potential optimal partition(s):
#>  - Criterion chosen to optimise the number of clusters:  increasing_step 
#>    (step quantile chosen:  0.99  (i.e., only the top 1 %  increase  in evaluation metrics  are used as break points for the number of clusters)
#>  - Optimal partition(s) of clusters for each metric:
#> anosim - 9
```
