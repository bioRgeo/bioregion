# 4.1 Hierarchical clustering

Hierarchical clustering consists in creating a hierarchical tree from a
matrix of distances (or beta-diversities). From this hierarchical tree,
clusters can be obtained by cutting the tree.

Although these methods are conceptually simple, their implementation can
be complex and requires important user decisions. Here, we provide a
step-by-step guide to performing hierarchical clustering analyses with
`bioregion`, along with comments on the philosophy on how we designed
the functions.

Hierarchical clustering takes place on the right-hand side of the
`bioregion` conceptual diagram:

![Workflow of the bioregion package for hierarchical
clustering.](../reference/figures/workflow_nonetwork.png)

## 1. Compute dissimilarity indices from input data

To initiate the hierarchical clustering procedure, you need to provide
pairwise distances between sites. These pairwise distances between sites
can be obtained by running
[`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
on a species-site matrix, such as a presence-absence or an abundance
matrix.

In the example below, we use the vegetation dataset from the package to
compute distance metrics.

``` r

# Work with the vegetation dataset we include in the package
data(vegemat)

# This is an abundance matrix where sites are in rows and species in columns
vegemat[1:10, 1:10]
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

We are going to compute the \\\beta\_{sim}\\ diversity metric, which is
a presence-absence dissimilarity index. The formula is as follows:
\\\beta\_{sim} = min(b, c) / (a+min(b, c))\\

Where *a* is the number of species shared by both sites; *b* is the
number of species occurring only in the first site; and *c* is the
number of species only occurring only in the second site.

We typically choose this metric for bioregionalization, because it is
the **turnover** component of the Sorensen index (Baselga, 2012) (in a
nutshell, it tells us how sites are different because they have distinct
species), and because it is less dependent on species richness than the
Jaccard turnover (Leprieur & Oikonomou, 2014). Alternatively, given that
we have abundance data here, we could also use the Bray-Curtis turnover
index (Baselga, 2013). The choice of the distance metric is very
important for the outcome of the clustering procedure, so we recommend
that you choose carefully depending on your research question.

``` r

dissim <- dissimilarity(vegemat)

head(dissim)
```

    ## Data.frame of dissimilarity between sites
    ##  - Total number of sites:  715 
    ##  - Total number of species:  3697 
    ##  - Number of rows:  255255 
    ##  - Number of dissimilarity metrics:  1 
    ## 
    ## 
    ##   Site1 Site2    Simpson
    ## 2    35    36 0.02325581
    ## 3    35    37 0.03100775
    ## 4    35    38 0.05426357
    ## 5    35    39 0.05426357
    ## 6    35    84 0.72093023
    ## 7    35    85 0.08527132

By default, only the Simpson index is computed, but other options are
available in the `metric` argument of dissimilarity(). Furthermore,
users can also write down their own formula to compute any index they
want for in the argument `formula`, see ?dissimilarity().

We are now ready to start the hierarchical clustering procedure with the
object `dissim` we have just created. Alternatively, you can also use
other types of objects in
[`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md),
such as a distance matrix object (class `dist`) or a `data.frame` of
your own crafting (make sure to read the required format carefully in
[`?hclu_hierarclust`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)).

## 2. Hierarchical clustering with basic parameters

Hierarchical clustering, and the associated hierarchical tree, can be
constructed in two ways: - *agglomerative*, where all sites are
initially assigned to their own bioregion and they are progressively
grouped together - *divisive*, where all sites initially belong to the
same unique bioregion and are then progressively divided into different
bioregions

Subsections 2.1 to 2.3 detail the functioning of agglomerative
hierarchical clustering, while [sub-section 2.4.](#id_2.4.) illustrates
the divisive method.

### 2.1 Basic usage

The basic use of the function is as follows:

``` r

tree1 <- hclu_hierarclust(dissim)
```

The function gives us some information as it proceeds. Notably, it talks
about a randomization of the dissimilarity matrix - this is a very
important feature because hierarchical clustering is strongly influenced
by the order of the sites in the distance matrix. Therefore, by default,
the function performs a randomization of the order of sites in the
distance matrix with 30 trials ([more information in the randomization
section](#randomization)). It also tells us that it built a consensus
tree based on an iterative hierarchical consensus algorithm, and it
found a tree with a cophenetic correlation coefficient of 0.7.

We can see type the name of the object in the console to see more
information:

``` r

tree1
```

    ## Clustering results for algorithm : hclu_hierarclust 
    ##  (hierarchical clustering based on a dissimilarity matrix)
    ##  - Number of sites:  715 
    ##  - Name of dissimilarity metric:  Simpson 
    ##  - Tree construction method:  average 
    ##  - Randomization of the dissimilarity matrix:  yes, number of trials 100 
    ##  - Method to compute the final tree:  Iterative hierarchical consensus tree 
    ##  - Cophenetic correlation coefficient:  0.7 
    ## Clustering procedure incomplete - no clusters yet

The last line tells us that the the clustering procedure is incomplete:
the tree has been built, but it has not yet been cut - so there are no
clusters in the object yet.

To cut the tree, we can use the
[`cut_tree()`](https://bioRgeo.github.io/bioregion/reference/cut_tree.md)
function:

``` r

# Ask for 3 clusters
tree1 <- cut_tree(tree1,
                  n_clust = 3)
```

Here, we asked for 3 clusters, and the algorithm automatically finds the
height at which 3 clusters are found (h = 0.547).

``` r

tree1
```

    ## Clustering results for algorithm : hclu_hierarclust 
    ##  (hierarchical clustering based on a dissimilarity matrix)
    ##  - Number of sites:  715 
    ##  - Name of dissimilarity metric:  Simpson 
    ##  - Tree construction method:  average 
    ##  - Randomization of the dissimilarity matrix:  yes, number of trials 100 
    ##  - Method to compute the final tree:  Iterative hierarchical consensus tree 
    ##  - Cophenetic correlation coefficient:  0.7 
    ##  - Number of clusters requested by the user:  3 
    ## Clustering results:
    ##  - Number of partitions:  1 
    ##  - Number of clusters:  3 
    ##  - Height of cut of the hierarchical tree: 0.547

When we type again the name of the object in the console, it gives us
the results of the clustering: we have

- **1 partition**: a partition is a clustering result. We only cut the
  tree once, so we only have 1 partition at the moment
- **4 clusters**: this is the number of clusters in the partition. We
  asked for 3, and we obtained 3, which is good. Sometimes, however, we
  cannot get the number of clusters we asked for - in which case the
  outcome will be indicated.
- **a height of cut at 0.547**: this is the height of cut at which we
  can obtain 4 clusters in our tree.

We can make a quick plot of our partitioned tree with

``` r

# We reduced the size of text labels with cex = .2, because there are too many sites
plot(tree1, cex = .2)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-7-1.png)

Let’s see how it looks like on a map:

``` r

#data(vegesf)
#map_bioregions(tree1, map = vegesf)
```

Now, this is a hierarchical tree, and cutting it only once (= only 1
partition) oversimplifies the result of the tree. Why not cut it
multiple times? For example, we could make deep, intermediate, and
shallow cuts to the tree, likewise to Ficetola *et al.* (2017), which
would allow us to see broad- to fine-scale relationships among sites in
our tree.

We can specify, e.g. 4, 10 and 20 clusters:

``` r

# Ask for 4, 10 and 20 clusters
tree1 <- cut_tree(tree1,
                  n_clust = c(2, 3, 12))

plot(tree1, cex = .2)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-9-1.png)

We can also see directly in the console the hierarchical structure of
the tree, with info on each cluster, by using the convenient
[`summary()`](https://rdrr.io/r/base/summary.html) function.

``` r

summary(tree1)
```

    ## 
    ## Summary of clustering results
    ## =============================
    ## 
    ## Algorithm:  hclu_hierarclust 
    ## Number of sites:  715 
    ## Number of bioregionalizations:  3 
    ## Hierarchical clustering: Yes
    ## 
    ## Bioregionalization 1: K_2
    ## ------------------------- 
    ## Total clusters:  2 
    ## Top 2 clusters by size:
    ##   Cluster 1: 538 items
    ##   Cluster 2: 177 items
    ## 
    ## Bioregionalization 2: K_3
    ## ------------------------- 
    ## Total clusters:  3 
    ## Top 3 clusters by size:
    ##   Cluster 1: 538 items
    ##   Cluster 3: 113 items
    ##   Cluster 2: 64 items
    ## 
    ## Bioregionalization 3: K_12
    ## -------------------------- 
    ## Total clusters:  12 
    ## Top 10 clusters by size:
    ##   Cluster 2: 414 items
    ##   Cluster 8: 99 items
    ##   Cluster 1: 63 items
    ##   Cluster 5: 62 items
    ##   Cluster 3: 36 items
    ##   Cluster 7: 14 items
    ##   Cluster 11: 10 items
    ##   Cluster 6: 8 items
    ##   Cluster 4: 5 items
    ##   Cluster 12: 2 items
    ##   ... and 2 more cluster(s)
    ## 
    ## Hierarchical structure
    ## ======================
    ## 
    ## 1 (n=538)
    ## └─1 (n=538)
    ##   ├─1 (n=63)
    ##   ├─11 (n=10)
    ##   ├─12 (n=2)
    ##   ├─2 (n=414)
    ##   ├─3 (n=36)
    ##   ├─4 (n=5)
    ##   └─6 (n=8)
    ## 
    ## 2 (n=177)
    ## ├─2 (n=64)
    ## │ ├─10 (n=1)
    ## │ ├─5 (n=62)
    ## │ └─9 (n=1)
    ## └─3 (n=113)
    ##   ├─7 (n=14)
    ##   └─8 (n=99)

However, it may be more useful to choose the heights of cut, rather than
the number of clusters. We could, for example, cut the tree at heights
0.4 (shallow cut), 0.5 (intermediate cut) and 0.6 (deep cut):

``` r

tree1 <- cut_tree(tree1,
                  cut_height = c(.4, .5, .6))

plot(tree1, cex = .2)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-11-1.png)

The plot is not easy to read because of the large number of sites. We
can rather extract the information directly from the object:

``` r

tree1
```

    ## Clustering results for algorithm : hclu_hierarclust 
    ##  (hierarchical clustering based on a dissimilarity matrix)
    ##  - Number of sites:  715 
    ##  - Name of dissimilarity metric:  Simpson 
    ##  - Tree construction method:  average 
    ##  - Randomization of the dissimilarity matrix:  yes, number of trials 100 
    ##  - Method to compute the final tree:  Iterative hierarchical consensus tree 
    ##  - Cophenetic correlation coefficient:  0.7 
    ##  - Heights of cut requested by the user:  0.4 0.5 0.6 
    ## Clustering results:
    ##  - Number of partitions:  3 
    ##  - Partitions are hierarchical
    ##  - Number of clusters:  2 9 24 
    ##  - Height of cut of the hierarchical tree: 0.6 0.5 0.4

From the result, we can read that for the deep cut partition (h = 0.6)
we have clusters, for the intermediate cut partition (h = 0.5) we have 9
clusters and for the shallow cut partition (h = 0.4) we have 24
clusters.

Let’s look at the hierarchical structure of clusters with
[`summary()`](https://rdrr.io/r/base/summary.html):

``` r

summary(tree1)
```

    ## 
    ## Summary of clustering results
    ## =============================
    ## 
    ## Algorithm:  hclu_hierarclust 
    ## Number of sites:  715 
    ## Number of bioregionalizations:  3 
    ## Hierarchical clustering: Yes
    ## 
    ## Bioregionalization 1: K_2
    ## ------------------------- 
    ## Total clusters:  2 
    ## Top 2 clusters by size:
    ##   Cluster 1: 538 items
    ##   Cluster 2: 177 items
    ## 
    ## Bioregionalization 2: K_9
    ## ------------------------- 
    ## Total clusters:  9 
    ## Top 9 clusters by size:
    ##   Cluster 2: 450 items
    ##   Cluster 6: 113 items
    ##   Cluster 1: 63 items
    ##   Cluster 4: 62 items
    ##   Cluster 5: 18 items
    ##   Cluster 3: 5 items
    ##   Cluster 9: 2 items
    ##   Cluster 7: 1 items
    ##   Cluster 8: 1 items
    ## 
    ## Bioregionalization 3: K_24
    ## -------------------------- 
    ## Total clusters:  24 
    ## Top 10 clusters by size:
    ##   Cluster 3: 296 items
    ##   Cluster 8: 106 items
    ##   Cluster 12: 97 items
    ##   Cluster 1: 59 items
    ##   Cluster 6: 49 items
    ##   Cluster 4: 33 items
    ##   Cluster 9: 14 items
    ##   Cluster 13: 11 items
    ##   Cluster 22: 10 items
    ##   Cluster 10: 6 items
    ##   ... and 14 more cluster(s)
    ## 
    ## Hierarchical structure
    ## ======================
    ## 
    ## 1 (n=538)
    ## ├─1 (n=63)
    ## │ ├─1 (n=59)
    ## │ └─2 (n=4)
    ## ├─2 (n=450)
    ## │ ├─11 (n=3)
    ## │ ├─13 (n=11)
    ## │ ├─24 (n=1)
    ## │ ├─3 (n=296)
    ## │ ├─4 (n=33)
    ## │ └─8 (n=106)
    ## ├─3 (n=5)
    ## │ └─5 (n=5)
    ## ├─5 (n=18)
    ## │ ├─10 (n=6)
    ## │ ├─22 (n=10)
    ## │ └─7 (n=2)
    ## └─9 (n=2)
    ##   └─23 (n=2)
    ## 
    ## 2 (n=177)
    ## ├─4 (n=62)
    ## │ ├─14 (n=4)
    ## │ ├─15 (n=3)
    ## │ ├─18 (n=4)
    ## │ ├─19 (n=1)
    ## │ ├─20 (n=1)
    ## │ └─6 (n=49)
    ## ├─6 (n=113)
    ## │ ├─12 (n=97)
    ## │ ├─17 (n=2)
    ## │ └─9 (n=14)
    ## ├─7 (n=1)
    ## │ └─16 (n=1)
    ## └─8 (n=1)
    ##   └─21 (n=1)

Here is how the maps look like:

``` r

#for(i in 2:ncol(tree1$clusters)){
#  map_bioregions(tree1$clusters[, c(1, i)], vegesf)
#}
```

In the next section we will see what are the default settings and why we
chose them, and then we will see how to find optimal numbers of
clusters.

### 2.2 Exploring the outputs

To explore the object, you can use
[`str()`](https://rdrr.io/r/utils/str.html) to see the object structure:

``` r

str(tree1)
```

    ##  $ name        : chr "hclu_hierarclust"
    ##  $ args        :List of 16
    ##   ..$ index              : chr "Simpson"
    ##   ..$ method             : chr "average"
    ##   ..$ randomize          : logi TRUE
    ##   ..$ seed               : NULL
    ##   ..$ n_runs             : num 100
    ##   ..$ optimal_tree_method: chr "iterative_consensus_tree"
    ##   ..$ keep_trials        : chr "no"
    ##   ..$ n_clust            : NULL
    ##   ..$ cut_height         : num [1:3] 0.4 0.5 0.6
    ##   ..$ find_h             : logi TRUE
    ##   ..$ h_max              : num 1
    ##   ..$ h_min              : num 0
    ##   ..$ consensus_p        : num 0.5
    ##   ..$ show_hierarchy     : logi FALSE
    ##   ..$ verbose            : logi TRUE
    ##   ..$ dynamic_tree_cut   : logi FALSE
    ##  $ inputs      :List of 9
    ##   ..$ bipartite      : logi FALSE
    ##   ..$ weight         : logi TRUE
    ##   ..$ pairwise       : logi TRUE
    ##   ..$ pairwise_metric: chr "Simpson"
    ##   ..$ dissimilarity  : logi TRUE
    ##   ..$ nb_sites       : int 715
    ##   ..$ data_type      : chr "occurrence"
    ##   ..$ node_type      : chr "site"
    ##   ..$ hierarchical   : logi TRUE
    ##  $ algorithm   :List of 6
    ##   ..$ final.tree         :List of 5
    ##   .. ..- attr(*, "class")= chr "hclust"
    ##   ..$ final.tree.coph.cor: num 0.7
    ##   ..$ final.tree.msd     : num 0.0221
    ##   ..$ trials             : chr "Trials not stored in output"
    ##   ..$ output_n_clust     : Named int [1:3] 2 9 24
    ##   .. ..- attr(*, "names")= chr [1:3] "h_0.6" "h_0.5" "h_0.4"
    ##   ..$ output_cut_height  : num [1:3] 0.6 0.5 0.4
    ##  $ clusters    :'data.frame':    715 obs. of  4 variables:
    ##   ..$ ID  : chr [1:715] "1003" "1004" "1005" "1006" ...
    ##   ..$ K_2 : chr [1:715] "1" "1" "1" "1" ...
    ##   ..$ K_9 : chr [1:715] "1" "1" "1" "1" ...
    ##   ..$ K_24: chr [1:715] "1" "1" "1" "1" ...
    ##   ..- attr(*, "node_type")= chr [1:715] "site" "site" "site" "site" ...
    ##  $ cluster_info:'data.frame':    3 obs. of  3 variables:
    ##   ..$ partition_name      : chr [1:3] "K_2" "K_9" "K_24"
    ##   ..$ n_clust             : int [1:3] 2 9 24
    ##   ..$ requested_cut_height: num [1:3] 0.6 0.5 0.4

It show you the different slots in the object, and how you can access
them. For example, if I want to access the `clusters` slot, I have to
type `tree1$clusters`.

- **name**: the name of the method we are using
- **args**: the arguments you have selected for your tree
- **inputs**: this is mostly for internal use in the package, it
  provides some info about the nature of input data and methods
- **algorithm**: this slot contains detailed information about the
  hierarchical clustering. For example, you can have access to the raw
  tree here, in `hclust` format. To access it, I can type
  `tree1$algorithm$final.tree`
- **clusters**: this is a `data.frame` containing your partitions. The
  first column is your sites, and all the other columns are the
  partitions.
- **cluster_info**: this is a small `data.frame` which will help you
  link your requests with the `clusters` `data.frame`. Its content
  varies depending on your choices; for example, in my case, it looks
  like this:

``` r

tree1$cluster_info
```

    ##       partition_name n_clust requested_cut_height
    ## h_0.6            K_2       2                  0.6
    ## h_0.5            K_9       9                  0.5
    ## h_0.4           K_24      24                  0.4

It shows the name of the partition (corresponding to column names in
`tree1$clusters`), the number of clusters in each partition, and the cut
height I initially requested.

### 2.3 Explanation of the default settings and how to change them

#### 2.3.1 Randomization of the distance matrix

The order of sites in the distance matrix influences the outcome of the
hierarchical tree. Let’s see that with an example:

``` r

# Compute the tree without randomizing the distance matrix
tree2 <- hclu_hierarclust(dissim,
                          randomize = FALSE)
plot(tree2, cex = .1)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-17-1.png)

This is how the tree looks like when the matrix is not randomized. Now
let’s randomize it and regenerate the tree:

``` r

# This line randomizes the order of rows in the distance matrix
dissim_random <- dissim[sample(1:nrow(dissim)), ]

# Recompute the tree
tree3 <- hclu_hierarclust(dissim_random,
                          randomize = FALSE)
plot(tree3, cex = .1)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-18-1.png)

See how the tree looks different? This is problematic because it means
that the outcome is heavily influenced by the order of sites in the
distance matrix.

To address this issue, we have developed an iterative algorithm that
will reconstruct the entire tree from top to bottom, by selecting for
each branch a majority decision among multiple randomizations of the
distance matrix (100 times by default, can be increased). This method is
called **Iterative Hierarchical Consensus Tree** (argument
`optimal_tree_method = "iterative_consensus_tree"`, default value) and
it ensures that you obtain a consensus tree that it will find a majority
decision for each branch of the tree. The tree produced with this method
generally have a better topology than any individual tree. We estimate
the performance of the topology with the **cophenetic correlation
coefficient**, which is the *correlation between the initial distance
\\\beta\_{sim}\\ among sites* and *the cophenetic distance*, which is
the distance at which sites are connected in the tree. It tells us how
representative is the tree of the initial distance matrix.

Although this method performs better than any other available method, it
comes with a computing cost: it needs to randomize the distance matrix
multiple times for each branching of the tree. Therefore, we recommend
using it to obtain a robust tree - be patient in that case. Otherwise,
if you only need a very fast look at a *good* tree, you can simply
select the single best tree among multiple randomization trials. It will
never be as good as the **Iterative Hierarchical Consensus Tree**, but
it will be a best choice for a fast exploration. To do that, choose
`optimal_tree_method = "best"`. It will select the tree that best
represents the distance matrix; i.e., the one that has the highest
cophenetic correlation coefficient among all trials.

Let’s see an example of `optimal_tree_method = "best"`. We can also  
ask the function to keep all individual trees for further exploration  
with `keep_trials = "all"` (including, for each trial, the randomized
matrix,  
the associated tree, and metrics for that tree).  
By default, `keep_trials = "no"` and `keep_trials = "metrics"` can be
used to  
keep only the metrics associated with each tree.

``` r

tree_best <- hclu_hierarclust(dissim_random,
                              randomize = TRUE,
                              optimal_tree_method = "best",
                              keep_trials = "all")
```

Another possible approach is to build a simple consensus tree among all
the trials. However, we generally do not recommend constructing a
consensus tree, because the topology of simple consensus trees can be
very problematic if there are a lot of ties in the distance matrix.
Let’s see it in action here:

``` r

tree_consensus <- hclu_hierarclust(dissim_random,
                                   randomize = TRUE,
                                   optimal_tree_method = "consensus",
                                   keep_trials = "no")
```

See how the cophenetic correlation coefficient for the consensus tree is
terrible compared to the IHCT and best tree ?

This consensus tree has almost no correlation with our initial distance
matrix. This is because its topology is terribly wrong, see how the tree
looks like:

``` r

plot(tree_consensus)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-21-1.png)

Booo, this is just a large rake, not a tree!!!

#### 2.3.2 Tree construction algorithm

By default, the function uses the UPGMA method (Unweighted Pair Group
Method with Arithmetic Mean) because it has been recommended in
bioregionalization for its better performance over other approaches
(Kreft & Jetz, 2010). You can change this method by changing the
argument `method`; all methods implemented in
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) are available.

Note that the current height distances for methods `method = "ward.D"`
and `method = "ward.D2"` may differ from the calculations in
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html), due to the
iterative nature of the algorithm. In addition, using the method
`method = "single"` are much slower than the other approaches, and we
have not yet implemented a workaround to make it faster (do not hesitate
to contact us if you need a faster implementation or have an idea of how
to make it run faster).

#### 2.3.3 Cutting the tree

There are three ways of cutting the tree:

1.  **Specify the expected number of clusters**: you can request a
    specific number of clusters (`n_clust = 5` for example). You can
    also request multiple partitions, each with their own number of
    clusters (`n_clust = c(5, 10, 15)`) for example.

**Note:** When you specify the number of clusters, the function will
search for the associated height of cut automatically; you can disable
this parameter with `find_h = FALSE`. It will search for this *h* value
between `h_max` (default 1) and `h_min` (default 0). These arguments can
be adjusted if you are working with indices whose values do not range
between 0 and 1.

2.  **Specify the height of cut**: you can request the height at which
    you want to cut the tree (e.g., `cut_height = 0.5`). You can also
    request multiple partitions, each with their own cut height
    (`cut_height = c(0.4, 0.5, 0.6)`) for example.

3.  **Use a dynamic tree cut method**: Rather than cutting the entire
    tree at once, this alternative approach consists in cutting
    individual branches at different heights. This method can be
    requested by using `dynamic_tree_cut = TRUE`, and is based on the
    dynamicTreeCut R package.

### 2.4. Divisive clustering

While the agglomerative hierarchical clustering in the previous
subsections followed a bottom-up approach, divisive clustering follows a
top-down approach. This means that in the first step of the clustering,
all sites belong to the same bioregion, and then sites are iteratively
divided into different bioregions until all sites belong to a unique
bioregion.

Divisive clustering is following the DIvisive ANAlysis (DIANA)
clustering algorithm described in (Kaufman & Rousseeuw, 2009).

At each step, the algorithm splits the largest cluster by identifying
the most dissimilar observation (i.e. site) and then putting sites that
are closer to this most dissimilar observation than to the ‘old party’
group into a splinter group. The result is that the large cluster is
split into two clusters.

The function `hclu_diana` performs the Diana divisive clustering.

``` r

# Compute the tree with the Diana algorithm
tree_diana <- hclu_diana(dissim)

plot(tree_diana)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-22-1.png)

## 3. How to find an optimal number of clusters?

![How to find an optimal number of
clusters?](../reference/figures/find_optimal_n.png)

1.  Step 1. **Build a tree** with
    [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)

2.  Step 2. **Explore a range of partitions**, from a minimum (e.g.,
    starting at 2 clusters) up to a maximum (e.g. \\n-1\\ clusters where
    \\n\\ is the number of sites).

3.  Step 3. **Calculate one or several metrics for each partition**, to
    be used as the basis for evaluation plots.

4.  Step 4. **Search for one or several optimal number(s) of clusters
    using evaluation plots**. Different criteria can be applied to
    identify the optimal number(s) of clusters.

5.  Step 5. **Export the optimal partitions from your cluster object.**

### 3.1 A practical example

In this example we will compute the evaluation metric used by Holt *et
al.* (2013), which compares the total dissimilarity with the
inter-cluster dissimilarity (sum of distances between clusters). Then we
will choose the optimal number of clusters as the elbow of the
evaluation plot.

``` r

data(vegemat)

# Calculate dissimilarities
dissim <- dissimilarity(vegemat)

# Step 1 & 2. Compute the tree and cut it into many different partitions
tree4 <- hclu_hierarclust(dissim,
                          n_clust = 2:100)

# Step 3. Calculate the same evaluation metric as Holt et al. 2013
eval_tree4 <- bioregionalization_metrics(tree4, 
                                         eval_metrics = "prop_between_dissim",
                                         dissimilarity = dissim)

# Step 4. Find the optimal number of clusters
opti_n_tree4 <- find_optimal_n(eval_tree4)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-24-1.png)

``` r

opti_n_tree4
```

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  elbow 
    ##  - Optimal partition(s) of clusters for each metric:

``` r

# Step 5. Extract the optimal number of clusters
# We get the name of the correct partition in the next line
K_name <- opti_n_tree4$evaluation_df$partition[opti_n_tree4$evaluation_df$optimal_n_prop_between_dissim]
# Look at the site-cluster table
head(tree4$clusters[, c("ID", K_name)])
```

    ##        ID K_14
    ## 1003 1003    1
    ## 1004 1004    1
    ## 1005 1005    1
    ## 1006 1006    1
    ## 1007 1007    2
    ## 1008 1008    2

``` r

# Make a map of the clusters
#data("vegesf")
#library("sf")
#map_bioregions(tree4$clusters[, c("ID", K_name)], vegesf)
```

Or if you are allergic to lines of code, you could also simply recut
your tree at the identified optimal number of cut-offs with
[`cut_tree()`](https://bioRgeo.github.io/bioregion/reference/cut_tree.md).

### 3.2 Evaluation metrics

Currently, there are four evaluation metrics available in the package
described
[here](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregionalization).

**Important note**

To be able to calculate `prop_between_dissim` and `anosim`, you need to
provide your dissimilarity object to the argument `dissimilarity`. In
addition, to be able to calculate `mean_endemics` and `tot_endemics`,
you need to provide your species-site co-occurrence matrix to the
argument `comat`.

Let’s see that in practice. Depending on the size of your dataset,
computing endemism-based metrics can take a while.

``` r

# Calculate prop_between_dissim and anosim
bioregionalization_metrics(tree4, 
                           eval_metrics = c("prop_between_dissim", "anosim"),
                           dissimilarity = dissim)
```

    ##    partition n_bioregions prop_between_dissim    anosim
    ## 1      K_2_1            2           0.5039255 0.7038453
    ## 2      K_2_2            2           0.5039255 0.7038453
    ## 3        K_4            4           0.6115668 0.7055336
    ## 4        K_5            5           0.6476470 0.7255243
    ## 5        K_6            6           0.6490485 0.7262948
    ## 6        K_7            7           0.6494868 0.7264502
    ## 7        K_8            8           0.6497724 0.7265047
    ## 8        K_9            9           0.6521960 0.7267220
    ## 9       K_10           10           0.7797049 0.7412313
    ## 10      K_11           11           0.7810061 0.7419956
    ## 11      K_12           12           0.7820362 0.7425334
    ## 12      K_13           13           0.7866725 0.7452144
    ## 13      K_14           14           0.8304847 0.7763484
    ## 14      K_15           15           0.8311333 0.7767871
    ## 15      K_16           16           0.8321952 0.7774242
    ## 16      K_17           17           0.8339339 0.7784386
    ## 17      K_18           18           0.8341172 0.7785430
    ## 18      K_19           19           0.8422558 0.7829579
    ## 19      K_20           20           0.8462671 0.7855295
    ## 20      K_21           21           0.8586412 0.7937815
    ## 21      K_22           22           0.8601130 0.7949075
    ## 22      K_23           23           0.8601494 0.7949336
    ## 23      K_24           24           0.8601767 0.7949526
    ## 24      K_25           25           0.8665417 0.7995322
    ## 25      K_26           26           0.8665682 0.7995448
    ## 26      K_27           27           0.8665800 0.7995504
    ## 27    K_29_1           29           0.8670966 0.7996385
    ## 28    K_29_2           29           0.8670966 0.7996385
    ## 29      K_30           30           0.8682584 0.8000748
    ## 30      K_31           31           0.8683940 0.8001369
    ## 31      K_32           32           0.8684170 0.8001458
    ## 32      K_33           33           0.8684199 0.8001461
    ## 33      K_34           34           0.8701128 0.8004195
    ## 34      K_35           35           0.8987808 0.8176325
    ## 35      K_36           36           0.8987950 0.8176437
    ## 36      K_37           37           0.8992759 0.8181359
    ## 37      K_38           38           0.8995875 0.8183957
    ## 38      K_39           39           0.8997950 0.8185833
    ## 39      K_40           40           0.8998033 0.8185881
    ## 40      K_41           41           0.8998339 0.8186134
    ## 41      K_42           42           0.9010559 0.8195488
    ## 42      K_43           43           0.9010752 0.8195687
    ## 43      K_44           44           0.9010917 0.8195774
    ## 44      K_45           45           0.9011191 0.8195968
    ## 45      K_46           46           0.9012094 0.8196792
    ## 46      K_47           47           0.9012203 0.8196847
    ## 47      K_48           48           0.9053491 0.8222804
    ## 48      K_49           49           0.9055867 0.8223740
    ## 49      K_50           50           0.9056719 0.8224158
    ## 50      K_51           51           0.9063492 0.8225512
    ## 51      K_52           52           0.9063730 0.8225596
    ## 52      K_53           53           0.9068066 0.8226861
    ## 53      K_54           54           0.9068429 0.8226932
    ## 54      K_55           55           0.9068585 0.8226927
    ## 55    K_57_1           57           0.9407998 0.8318092
    ## 56    K_57_2           57           0.9407998 0.8318092
    ## 57      K_58           58           0.9408908 0.8320144
    ## 58      K_59           59           0.9409134 0.8320387
    ## 59      K_60           60           0.9411474 0.8322919
    ## 60      K_61           61           0.9411549 0.8322979
    ## 61      K_62           62           0.9417502 0.8330224
    ## 62      K_63           63           0.9417724 0.8330419
    ## 63      K_64           64           0.9418020 0.8330797
    ## 64      K_65           65           0.9419761 0.8332177
    ## 65      K_66           66           0.9419908 0.8332309
    ## 66      K_67           67           0.9420202 0.8332548
    ## 67      K_68           68           0.9431480 0.8342166
    ## 68      K_69           69           0.9431504 0.8342171
    ## 69      K_70           70           0.9431697 0.8342267
    ## 70      K_71           71           0.9475710 0.8369940
    ## 71      K_72           72           0.9504739 0.8397503
    ## 72      K_73           73           0.9504787 0.8397522
    ## 73      K_74           74           0.9505097 0.8397891
    ## 74      K_75           75           0.9526047 0.8419845
    ## 75      K_76           76           0.9526546 0.8421129
    ## 76      K_77           77           0.9534397 0.8428665
    ## 77      K_78           78           0.9534469 0.8428795
    ## 78      K_79           79           0.9534516 0.8428818
    ## 79      K_80           80           0.9547543 0.8441475
    ## 80      K_81           81           0.9547590 0.8441498
    ## 81      K_82           82           0.9547637 0.8441521
    ## 82      K_83           83           0.9548903 0.8442458
    ## 83      K_84           84           0.9552727 0.8444641
    ## 84      K_85           85           0.9552958 0.8444718
    ## 85      K_86           86           0.9597738 0.8475812
    ## 86      K_87           87           0.9597784 0.8475819
    ## 87      K_88           88           0.9598124 0.8476134
    ## 88      K_89           89           0.9601166 0.8478113
    ## 89      K_90           90           0.9601437 0.8478264
    ## 90      K_91           91           0.9603125 0.8480385
    ## 91      K_92           92           0.9603147 0.8480383
    ## 92      K_93           93           0.9619450 0.8493302
    ## 93      K_94           94           0.9624150 0.8501088
    ## 94      K_95           95           0.9637310 0.8511409
    ## 95      K_96           96           0.9637467 0.8511462
    ## 96   K_100_1          100           0.9805066 0.8635191
    ## 97   K_100_2          100           0.9805066 0.8635191
    ## 98   K_100_3          100           0.9805066 0.8635191
    ## 99   K_100_4          100           0.9805066 0.8635191

``` r

# Calculate mean_endemics and tot_endemics
bioregionalization_metrics(tree4,
                           eval_metrics = c("mean_endemics", "tot_endemics"),
                           dissimilarity = NULL,
                           comat = vegemat)
```

    ##    partition n_bioregions mean_endemics tot_endemics
    ## 1      K_2_1            2   0.177309950   0.30916960
    ## 2      K_2_2            2   0.177309950   0.30916960
    ## 3        K_4            4   0.046560673   0.15120368
    ## 4        K_5            5   0.032444543   0.12388423
    ## 5        K_6            6   0.026866440   0.12253178
    ## 6        K_7            7   0.023071178   0.12253178
    ## 7        K_8            8   0.020189595   0.12253178
    ## 8        K_9            9   0.017978651   0.12226129
    ## 9       K_10           10   0.015228298   0.10738437
    ## 10      K_11           11   0.013815277   0.10711388
    ## 11      K_12           12   0.012589709   0.10657290
    ## 12      K_13           13   0.012022577   0.10035164
    ## 13      K_14           14   0.010259001   0.08520422
    ## 14      K_15           15   0.009604125   0.08520422
    ## 15      K_16           16   0.009003867   0.08520422
    ## 16      K_17           17   0.008520651   0.08520422
    ## 17      K_18           18   0.007993111   0.08466324
    ## 18      K_19           19   0.007539460   0.08412226
    ## 19      K_20           20   0.007131691   0.08358128
    ## 20      K_21           21   0.006761121   0.08222883
    ## 21      K_22           22   0.006487602   0.08222883
    ## 22      K_23           23   0.006208174   0.08222883
    ## 23      K_24           24   0.005950796   0.08222883
    ## 24      K_25           25   0.005719272   0.08195834
    ## 25      K_26           26   0.005500096   0.08195834
    ## 26      K_27           27   0.005296389   0.08195834
    ## 27    K_29_1           29   0.004932414   0.08195834
    ## 28    K_29_2           29   0.004932414   0.08195834
    ## 29      K_30           30   0.004766607   0.08141737
    ## 30      K_31           31   0.004613189   0.08141737
    ## 31      K_32           32   0.004470415   0.08141737
    ## 32      K_33           33   0.004334948   0.08141737
    ## 33      K_34           34   0.004225648   0.08114688
    ## 34      K_35           35   0.004094471   0.07952394
    ## 35      K_36           36   0.003982291   0.07952394
    ## 36      K_37           37   0.003863717   0.07898296
    ## 37      K_38           38   0.003748818   0.07871247
    ## 38      K_39           39   0.003641791   0.07844198
    ## 39      K_40           40   0.003550746   0.07844198
    ## 40      K_41           41   0.003468117   0.07844198
    ## 41      K_42           42   0.003413486   0.07844198
    ## 42      K_43           43   0.003334384   0.07844198
    ## 43      K_44           44   0.003267420   0.07844198
    ## 44      K_45           45   0.003195549   0.07844198
    ## 45      K_46           46   0.003099440   0.07790100
    ## 46      K_47           47   0.003033494   0.07790100
    ## 47      K_48           48   0.002947683   0.07654855
    ## 48      K_49           49   0.002879376   0.07627806
    ## 49      K_50           50   0.002877231   0.07600757
    ## 50      K_51           51   0.002859146   0.07600757
    ## 51      K_52           52   0.002809455   0.07600757
    ## 52      K_53           53   0.002758729   0.07600757
    ## 53      K_54           54   0.002709092   0.07600757
    ## 54      K_55           55   0.002660079   0.07600757
    ## 55    K_57_1           57   0.002562755   0.07519610
    ## 56    K_57_2           57   0.002562755   0.07519610
    ## 57      K_58           58   0.002566416   0.07519610
    ## 58      K_59           59   0.002524111   0.07519610
    ## 59      K_60           60   0.002513575   0.07465513
    ## 60      K_61           61   0.002472369   0.07465513
    ## 61      K_62           62   0.002442548   0.07465513
    ## 62      K_63           63   0.002407113   0.07465513
    ## 63      K_64           64   0.002369502   0.07465513
    ## 64      K_65           65   0.002326499   0.07438464
    ## 65      K_66           66   0.002275835   0.07411415
    ## 66      K_67           67   0.002245839   0.07411415
    ## 67      K_68           68   0.002218461   0.07411415
    ## 68      K_69           69   0.002189388   0.07411415
    ## 69      K_70           70   0.002159766   0.07411415
    ## 70      K_71           71   0.002129981   0.07411415
    ## 71      K_72           72   0.001772702   0.05788477
    ## 72      K_73           73   0.001760624   0.05788477
    ## 73      K_74           74   0.001722344   0.05734379
    ## 74      K_75           75   0.001704328   0.05734379
    ## 75      K_76           76   0.001685072   0.05734379
    ## 76      K_77           77   0.001663286   0.05734379
    ## 77      K_78           78   0.001641962   0.05734379
    ## 78      K_79           79   0.001621178   0.05734379
    ## 79      K_80           80   0.001605577   0.05707330
    ## 80      K_81           81   0.001585755   0.05707330
    ## 81      K_82           82   0.001566417   0.05707330
    ## 82      K_83           83   0.001547566   0.05707330
    ## 83      K_84           84   0.001529233   0.05707330
    ## 84      K_85           85   0.001511242   0.05707330
    ## 85      K_86           86   0.001497222   0.05680281
    ## 86      K_87           87   0.001480012   0.05680281
    ## 87      K_88           88   0.001463324   0.05680281
    ## 88      K_89           89   0.001436875   0.05626183
    ## 89      K_90           90   0.001429062   0.05626183
    ## 90      K_91           91   0.001418318   0.05626183
    ## 91      K_92           92   0.001402902   0.05626183
    ## 92      K_93           93   0.001389665   0.05626183
    ## 93      K_94           94   0.001367387   0.05572085
    ## 94      K_95           95   0.001353125   0.05572085
    ## 95      K_96           96   0.001339030   0.05572085
    ## 96   K_100_1          100   0.001285513   0.05490939
    ## 97   K_100_2          100   0.001285513   0.05490939
    ## 98   K_100_3          100   0.001285513   0.05490939
    ## 99   K_100_4          100   0.001285513   0.05490939

### 3.2 Criteria to choose an optimal number of clusters

Choosing the optimal number of clusters is a long-standing issue in the
literature, and there is no absolute and objective answer to this
question. A plethora of methods have been proposed over the years, and
the best approach to tackle this issue is probably to compare the
results of multiple approaches to make an informed decision.

In the `bioregion` package, we have implemented several methods that are
specifically suited for bioregionalization analysis.

For example, a standard criterion used for identifying the optimal
number of clusters is the **elbow method**, which is the default
criterion in
[`find_optimal_n()`](https://bioRgeo.github.io/bioregion/reference/find_optimal_n.md).
However, we recommend moving beyond the paradigm of a single optimal
number of clusters, which is likely an oversimplification of the
hierarchy of the tree. We recommend considering multiple cuts of the
tree, and we provide several methods for doing so: **identifying large
steps in the curve** or **using multiple cutoffs**. Additionally, we
implement other approaches, such as using the maximum or minimum value
of the metrics, or by finding break points in the curve with a segmented
model.

Before we look at these different methods, we will compute all the
evaluation metrics and store them in `eval_tree4`:

``` r

eval_tree4 <- bioregionalization_metrics(tree4, 
                                         eval_metrics = "all",
                                         dissimilarity = dissim, 
                                         comat = vegemat)
```

#### 3.2.1 Elbow method

The elbow method consists in find the ‘elbow’ in the form of the
metric-cluster relationship. This method will typically work for metrics
which have an L-shaped form (typically, `prop_between_dissim` and
endemism metrics), but not for other metrics (e.g. the form of `anosim`
does not necessarily follow an L-shape).

*The rationale behind the elbow method is to find a cutoff above which
the metric values stop increasing significantly, such that adding new
clusters does not provide much significant improvement in the metric.*

The elbow method is the default method in
[`find_optimal_n()`](https://bioRgeo.github.io/bioregion/reference/find_optimal_n.md).
There are no parameters to adjust it. If the curve is not elbow-shaped,
it may give spurious results.

``` r

find_optimal_n(eval_tree4)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-27-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  prop_between_dissim anosim mean_endemics tot_endemics 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  elbow 
    ##  - Optimal partition(s) of clusters for each metric:

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-28-1.png)

In our example above, the optimal number of clusters varies depending on
the metric, from a minimum of 10 to a maximum of 35. The final choice
depends on your metric preferences with respect to metrics, and your
objectives with the clustering. Alternatively, two cut-offs could be
used, a deep cut-off based on the endemism metrics e.g. at a value of
10, and a shallow cutoff based on `prop_between_dissim`, at 14.

#### 3.2.2 Step method

The step method consists in identifying the largest “steps” in metrics,
i.e., the largest increases or decreases in the value of the metric.

To do this, the function calculates all successive differences in
metrics between partitions. It will then keep only the largest positive
differences (`increasing_step`) or negative differences
(`decreasing_step`). `increasing_step` is for increasing metrics
(`prop_between_dissim`) and `decreasing_step` is for decreasing metrics
(`Mean_Endemics` and `tot_endemics`). `anosim` values can either
increase or decrease depending on your dataset, so you would have to
explore both ways.

By default, the function selects the top 1% steps:

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = c("anosim", "prop_between_dissim"),
               criterion = "increasing_step")
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-29-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  anosim prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  increasing_step 
    ##    (step quantile chosen:  0.99  (i.e., only the top 1 %  increase  in evaluation metrics  are used as break points for the number of clusters)
    ##  - Optimal partition(s) of clusters for each metric:

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = c("mean_endemics", "tot_endemics"),
               criterion = "decreasing_step")
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-29-2.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  mean_endemics tot_endemics 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  decreasing_step 
    ##    (step quantile chosen:  0.99  (i.e., only the top 1 %  decrease  in evaluation metrics  are used as break points for the number of clusters)
    ##  - Optimal partition(s) of clusters for each metric:

However, you can adjust it in two different ways. First, choose a number
of steps to select, e.g. to select the largest 3 steps, use
`step_levels = 3`:

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = c("anosim", "prop_between_dissim"),
               criterion = "increasing_step",
               step_levels = 3)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-30-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  anosim prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  increasing_step 
    ##    (step quantile chosen:  0.99  (i.e., only the top 1 %  increase  in evaluation metrics  are used as break points for the number of clusters)
    ##  - Optimal partition(s) of clusters for each metric:

Note that these steps generally correspond to large jumps in the tree,
which is why we like this approach as it fits well with the hierarchical
nature of the tree.

Second, you can set a quantile of steps to select, e.g. to select the 5%
largest steps set the quantile to 0.95 (`step_quantile = 0.95`):

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = c("anosim", "prop_between_dissim"),
               criterion = "increasing_step",
               step_quantile = 0.95)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-31-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  anosim prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  increasing_step 
    ##    (step quantile chosen:  0.95  (i.e., only the top 5 %  increase  in evaluation metrics  are used as break points for the number of clusters)
    ##  - Optimal partition(s) of clusters for each metric:

Finally, a question that may arise is which cluster number to select
when a large step occurs. For example, if the largest step occurs
between partitions with 4 and 5 clusters, should we keep the partition
with 4 clusters, or the partition with 5 clusters?

By default, the function keeps the partition with \\N + 1\\ (5 clusters
in our example above). You can change this by setting
`step_round_above = FALSE`:

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = c("anosim", "prop_between_dissim"),
               criterion = "decreasing_step",
               step_round_above = FALSE)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-32-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  anosim prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  decreasing_step 
    ##    (step quantile chosen:  0.99  (i.e., only the top 1 %  decrease  in evaluation metrics  are used as break points for the number of clusters)
    ##  - Optimal partition(s) of clusters for each metric:

#### 3.2.3 Cutting at different cut-off values

The idea of this method is to select specific metric values at which the
number of clusters should be used. For example, in their study, Holt *et
al.* (2013) used different cutoffs for `BetweenDissim` to find the
global biogeographic regions: 0.90, 0.95, 0.99, 0.999. The higher the
value, the more -diversity is explained, but also the more clusters
there are. Therefore, the choice is a trade-off between the total
-diversity explained and the number of clusters.

Eventually, the choice of these values depends on different factors:

1.  The geographic scope of your study. A global scale study can use
    large cutoffs like Holt *et al.* (2013) and end up with a reasonable
    number of clusters, whereas in regional to local scale studies with
    less endemism and more taxa shared among clusters, these values are
    too high, and other cutoffs should be explored, such as 0.5 and
    0.75.

2.  The characteristics of your study which will increase or decrease
    the degree of endemism among clusters: dispersal capacities of your
    taxonomic group, the connectivity/barriers in your study area, etc.
    Use lower cutoffs when you have a large number of widespread
    species, use higher cutoffs when you have high degrees of endemism.

3.  Using abundance or phylogenetic data to compute -diversity metrics
    may allow you to better distinguish clusters, which in turn will
    allow you to use higher cutoffs.

For example, in our case, a regional-scale study on vegetation, we can
use three cutoffs: 0.6 (deep cutoff), 0.8 (intermediate cutoff), and 0.9
(shallow cutoff).

``` r

find_optimal_n(eval_tree4,
               metrics_to_use = "prop_between_dissim",
               criterion = "cutoff",
               metric_cutoffs = c(.6, .8, .9))
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-33-1.png)

    ## Search for an optimal number of clusters:
    ##  - 99  partition(s) evaluated
    ##  - Range of clusters explored: from  2  to  100 
    ##  - Evaluated metric(s):  prop_between_dissim 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  cutoff 
    ##    --> cutoff(s) chosen:  0.6 0.8 0.9 
    ##  - Optimal partition(s) of clusters for each metric:

#### 3.2.4 Cutting at the maximum or minimum metric value

This criterion finds the maximum (`criterion = "max"`) or minimum
(`criterion = "min"`) value of the metric in the list of partitions and
selects the corresponding partition. Such a criterion can be interesting
in the case of `anosim`, but is probably much less useful for the other
metrics implemented in the package.

#### 3.2.5 Finding break points in the curve

This criterion consists in applying a segmented regression model with
the formula evaluation metric ~ number of clusters. The user can define
the number of breaks to be identified on the curve. Note that such a
model is likely to require a minimum number of points to find an
appropriate number of clusters. In our example here, we make 100 cuts in
the tree to have enough values.

``` r

tree5 <- cut_tree(tree4,
                  cut_height = seq(0, max(tree4$algorithm$final.tree$height), 
                                   length = 100)) 

eval_tree5 <- bioregionalization_metrics(tree5, 
                                         eval_metrics = "all",
                                         dissimilarity = dissim, 
                                         comat = vegemat)

find_optimal_n(eval_tree5,
               criterion = "breakpoints")
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-34-1.png)

    ## Search for an optimal number of clusters:
    ##  - 100  partition(s) evaluated
    ##  - Range of clusters explored: from  1  to  701 
    ##  - Evaluated metric(s):  prop_between_dissim anosim mean_endemics tot_endemics 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  breakpoints 
    ##  - Optimal partition(s) of clusters for each metric:

We can ask for a higher number of breaks:

- 2 breaks

``` r

find_optimal_n(eval_tree5,
               criterion = "breakpoints",
               n_breakpoints = 2)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-35-1.png)

    ## Search for an optimal number of clusters:
    ##  - 100  partition(s) evaluated
    ##  - Range of clusters explored: from  1  to  701 
    ##  - Evaluated metric(s):  prop_between_dissim anosim mean_endemics tot_endemics 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  breakpoints 
    ##  - Optimal partition(s) of clusters for each metric:

- 3 breaks

``` r

find_optimal_n(eval_tree5,
               criterion = "breakpoints",
               n_breakpoints = 3)
```

![](a4_1_hierarchical_clustering_files/figure-html/unnamed-chunk-36-1.png)

    ## Search for an optimal number of clusters:
    ##  - 100  partition(s) evaluated
    ##  - Range of clusters explored: from  1  to  701 
    ##  - Evaluated metric(s):  prop_between_dissim anosim mean_endemics tot_endemics 
    ## 
    ## Potential optimal partition(s):
    ##  - Criterion chosen to optimise the number of clusters:  breakpoints 
    ##  - Optimal partition(s) of clusters for each metric:

Increasing the number of breaks can be useful in situations where you
have, for example, non-linear silhouettes of metric ~ n clusters.

## 4. OPTICS as a semi-hierarchical clustering approach

OPTICS (Ordering Points To Identify the Clustering Structure) is a
semi-hierarchical clustering approach that orders the points in the
datasets so that the closest points become neighbors, calculates a
‘reachability’ distance between each point, and then extracts clusters
from this reachability distance in a hierarchical manner. However, the
hierarchical nature of clusters is not directly provided by the
algorithm in a tree-like output. Hence, users should explore the
‘reachability plot’ to understand the hierarchical nature of their
OPTICS clusters, and read the related publication to further understand
this method (Hahsler *et al.*, 2019).

To run the optics algorithm, use the
[`hclu_optics()`](https://bioRgeo.github.io/bioregion/reference/hclu_optics.md)
function:

``` r

data(vegemat)

# Calculate dissimilarities
dissim <- dissimilarity(vegemat)

clust1 <- hclu_optics(dissim)

clust1
```

    ## Clustering results for algorithm : hclu_optics 
    ##  - Number of sites:  715 
    ## Clustering results:
    ##  - Number of partitions:  1 
    ##  - Number of clusters:  9

## 5. References

Baselga A (2012) The relationship between species replacement,
dissimilarity derived from nestedness, and nestedness. *Global Ecology
and Biogeography* 21, 1223–1232.

Baselga A (2013) Separating the two components of abundance-based
dissimilarity: Balanced changes in abundance vs. Abundance gradients.
*Methods in Ecology and Evolution* 4, 552–557.

Ficetola GF, Mazel F & Thuiller W (2017) Global determinants of
zoogeographical boundaries. *Nature Ecology & Evolution* 1, 0089.

Hahsler M, Piekenbrock M & Doran D (2019) Dbscan: Fast density-based
clustering with r. *Journal of Statistical Software* 91, 1–30.

Holt BG, Lessard J-P, Borregaard MK *et al.* (2013) An update of
Wallace’s zoogeographic regions of the world. *Science* 339, 74–78.

Kaufman L & Rousseeuw PJ (2009) Finding groups in data: An introduction
to cluster analysis. In: *Finding groups in data: An introduction to
cluster analysis.* (ed JW & Sons.). John Wiley & Sons.

Kreft H & Jetz W (2010) A framework for delineating biogeographical
regions based on species distributions. *Journal of Biogeography* 37,
2029–2053.

Leprieur F & Oikonomou A (2014) The need for richness-independent
measures of turnover when delineating biogeographical regions. *Journal
of Biogeography* 41, 417–420.
