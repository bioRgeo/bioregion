#' Hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function generates a hierarchical tree from a dissimilarity
#' (beta-diversity) `data.frame`, calculates the cophenetic correlation
#' coefficient, and optionally retrieves clusters from the tree upon user 
#' request. The function includes a randomization process for the dissimilarity 
#' matrix to generate the tree, with two methods available for constructing the 
#' final tree. Typically, the dissimilarity `data.frame` is a
#' `bioregion.pairwise` object obtained by running `similarity`,
#' or by running `similarity` followed by `similarity_to_dissimilarity`.
#'
#' @param dissimilarity The output object from [dissimilarity()] or
#'  [similarity_to_dissimilarity()], or a `dist` object. 
#'  If a `data.frame` is used, the first two columns represent pairs of sites 
#'  (or any pair of nodes), and the subsequent column(s) contain the 
#'  dissimilarity indices.
#' 
#' @param index The name or number of the dissimilarity column to use. By 
#' default, the third column name of `dissimilarity` is used.
#' 
#' @param method The name of the hierarchical classification method, as in
#' [hclust][fastcluster::hclust]. Should be one of `"ward.D"`,
#' `"ward.D2"`, `"single"`, `"complete"`, `"average"`
#' (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC), or
#' `"centroid"` (= UPGMC).
#' 
#' @param randomize A `boolean` indicating whether the dissimilarity matrix 
#' should be randomized to account for the order of sites in the dissimilarity
#'  matrix.
#'  
#' @param seed A value for the random number generator (`NULL` for random by 
#' default).
#' 
#' @param n_runs The number of trials for randomizing the dissimilarity matrix.
#' 
#' @param keep_trials A `character` string indicating whether random trial results 
#' (including the randomized matrix, the associated tree and metrics for that tree) 
#' should be stored in the output object. Possible values are `"no"` (default), 
#' `"all"` or `"metrics"`. Note that this parameter is automatically set to 
#' `"no"` if `optimal_tree_method = "ihct"`.
#' 
#' @param optimal_tree_method A `character` string indicating how the final tree
#' should be obtained from all trials. Possible values are 
#' `"ihct"` (default), `"best"` or `"consensus"`. `"iterative_consensus_tree"`
#' is still accepted as another name for `"ihct"`, so that code written for
#' earlier versions of bioregion keeps working.
#' **We recommend `"ihct"`. See Details.**
#' 
#' @param n_clust An `integer` vector or a single `integer` indicating the 
#' number of clusters to be obtained from the hierarchical tree, or the output 
#' from [bioregionalization_metrics]. This parameter should not be used 
#' simultaneously with `cut_height`.
#' 
#' @param cut_height A `numeric` vector indicating the height(s) at which the
#' tree should be cut. This parameter should not be used simultaneously with 
#' `n_clust`.
#' 
#' @param find_h A `boolean` indicating whether the height of the cut should be 
#' found for the requested `n_clust`.
#' 
#' @param h_max A `numeric` value indicating the maximum possible tree height 
#' for the chosen `index`.
#' 
#' @param h_min A `numeric` value indicating the minimum possible height in the
#'  tree for the chosen `index`.
#' 
#' @param consensus_p A `numeric` value (applicable only if 
#' `optimal_tree_method = "consensus"`) indicating the threshold proportion of 
#' trees that must support a region/cluster for it to be included in the final 
#' consensus tree.
#' 
#' @param show_hierarchy A `boolean` specifying if the hierarchy of clusters
#' should be identifiable in the outputs (`FALSE` by default). This argument is
#' only used if the tree is cut (i.e., `n_clust` or `cut_height` is provided).
#'
#' @param ihct_top_n_trees An `integer` (applicable only if
#' `optimal_tree_method = "ihct"`) indicating how many of
#' the best randomized trees are used to decide each division of the tree
#' (`2` by default). See Details.
#'
#' @param ihct_variation_drop A `numeric` value between 0 and 1 (applicable
#' only if `optimal_tree_method = "ihct"` and
#' `method = "average"`). This is an optimization parameter: it reduces the
#' number of times the dissimilarity matrix is randomized while the tree is
#' built, which makes the tree faster to obtain for a very small cost in its
#' fit to the data. Its rule is based on variation: new randomizations are made
#' once a group of sites has lost this share of the variation that was left to
#' decide when they were last made. The default `0.2` therefore means
#' "randomize again once a fifth of what was left to decide has been decided".
#' Set it to `0` to randomize at every division, or to `1` to switch this rule
#' off and leave `ihct_sites_drop` to decide on its own. See Details.
#'
#' @param ihct_sites_drop A `numeric` value of 0 or more (applicable only if
#' `optimal_tree_method = "ihct"`, `method = "average"` and
#' `ihct_variation_drop > 0`). This is a second optimization parameter, used
#' together with `ihct_variation_drop`, with a rule based on sites rather than
#' on variation: new randomizations are also made once a group of sites has
#' lost this many sites since they were last made (`10` by default). Set it to
#' `0` or `1` to randomize at every division, or to `Inf` to switch this rule
#' off and leave `ihct_variation_drop` to decide on its own. See Details.
#'
#' @param ihct_height_rule A `character` string (applicable only if
#' `optimal_tree_method = "ihct"`) indicating how the
#' heights of the tree are corrected when a division comes out lower than a
#' division it contains. With `"least_squares"` (default) the heights are moved
#' as little as possible, which fits the dissimilarities better; with
#' `"max_child"` each division is raised to the highest division it contains,
#' as in bioregion 1.4.0 and earlier. See Details.
#'
#' @param ihct_n_workers An `integer` of 1 or more (applicable only if
#' `optimal_tree_method = "ihct"`) indicating how many
#' processes of your computer may build the randomized trees at the same time.
#' With `1` (default) they are built one after another, as before. Higher
#' values are worth it on large matrices only, and they do not change the
#' tree: the same `seed` gives the same result whatever this is set to. See
#' Details.
#'
#' @param verbose A `boolean` indicating whether to
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: A `character` string containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` describing the characteristics of the clustering process.}
#' \item{**algorithm**: A `list` containing all objects associated with the
#'  clustering procedure, such as the original cluster objects.}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' In the `algorithm` slot, users can find the following elements:
#'
#' \itemize{
#' \item{`trials`: A list containing all randomization trials. Each trial
#' includes the dissimilarity matrix with randomized site order, the
#' associated tree, and the cophenetic correlation coefficient for
#' that tree.}
#' \item{`final.tree`: An `hclust` object representing the final
#' hierarchical tree to be used.}
#' \item{`final.tree.coph.cor`: The cophenetic correlation coefficient
#' between the initial dissimilarity matrix and the `final.tree`.}
#' }
#'  
#' @details
#' The function is based on [hclust][fastcluster::hclust].
#' The default method for the hierarchical tree is `average`, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity dissimilarity (Kreft & Jetz, 2010).
#'
#' Clusters can be obtained by two methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in `n_clust`}
#' \item{Specifying one or several heights of cut in `cut_height`}}
#'
#' To find an optimal number of clusters, see [bioregionalization_metrics()]
#' 
#' It is important to pay attention to the fact that the order of rows
#' in the input distance matrix influences the tree topology as explained in 
#' Dapporto (2013). To address this, the function generates multiple trees by 
#' randomizing the distance matrix. 
#' 
#' Two methods are available to obtain the final tree:
#' \itemize{
#' 
#' \item{`optimal_tree_method = "ihct"`: The Iterative 
#' Hierarchical Consensus Tree (IHCT) method reconstructs a consensus tree by 
#' iteratively splitting the dataset into two subclusters based on the pairwise 
#' dissimilarity of sites across `n_runs` trees based on `n_runs` randomizations
#' of the distance matrix. At each iteration, it 
#' identifies the majority membership of sites into two stable groups across
#' all trees,
#' calculates the height based on the selected linkage method (`method`),
#' and enforces monotonic constraints on 
#' node heights to produce a coherent tree structure. 
#' This approach provides a robust, hierarchical representation of site 
#' relationships, balancing 
#' cluster stability and hierarchical constraints.}
#'
#
#' 
#' \item{`optimal_tree_method = "best"`: This method selects one tree among with 
#' the highest cophenetic correlation coefficient, representing the best fit 
#' between the hierarchical structure and the original distance matrix. }
#' 
#' \item{`optimal_tree_method = "consensus"`: This method constructs a consensus 
#' tree using phylogenetic methods with the function 
#' [consensus][ape::consensus].
#' When using this option, you must set the `consensus_p` parameter, which 
#' indicates 
#' the proportion of trees that must contain a region/cluster for it to be 
#' included 
#' in the final consensus tree. 
#' Consensus trees lack an inherent height because they represent a majority 
#' structure rather than an actual hierarchical clustering. To assign heights, 
#' we use a non-negative least squares method ([nnls.tree][phangorn::nnls.tree]) 
#' based on the initial distance matrix, ensuring that the consensus 
#' tree preserves 
#' approximate distances among clusters.}
#' }
#' 
#' 
#' We recommend using `"ihct"` as all the branches of
#' this tree will always reflect the majority decision among many randomized 
#' versions of the distance matrix. This method is inspired by 
#' Dapporto et al. (2015), which also used the majority decision
#' among many randomized versions of the distance matrix, but it expands it 
#' to reconstruct the entire topology of the tree iteratively. 
#' 
#' We do not recommend using the basic `consensus` method because in many 
#' contexts it provides inconsistent results, with a meaningless tree topology
#' and a very low cophenetic correlation coefficient. 
#' 
#' For a fast exploration of the tree, we recommend using the `best` method
#' which will only select the tree with the highest cophenetic correlation
#' coefficient among all randomized versions of the distance matrix. 
#'
#' ~ **Iterative Hierarchical Consensus Tree** details ~
#' 
#' The paragraphs below cover the `ihct_*` arguments of this function. The
#' algorithm itself is walked through step by step in [ihct()], which is the
#' function doing the work here, and which can also be called on its own if
#' you want the tree without the rest of this function.
#'
#' --> *Tree quality*
#' 
#' `ihct_top_n_trees` sets how many of the randomized trees, ranked by their 
#' quality (i.e., how well
#' they fit the dissimilarities with cophenetic correlation coefficient CCC),
#'  are used to decide a division: with `1` the
#' division is the top division of the best tree; with more, sites are grouped
#' according to how often they fall on the same side in these trees.
#' We recommend leaving this to default values, or to a low number of trees,
#' because it provided the best results in
#' our tests (highest CCC).
#'
#' --> *Computation time*
#'
#' The algorithm can be long to run because at every division of the tree it has
#' to randomize `n_runs` trees. On large datasets, building new
#' trees at every division makes this method slow. To make the function usable
#' on large datasets, we provide two optimization parameters. These parameters
#'  make the
#' algorithm reuse previously randomized trees at new division, unless a 
#' threshold of change is reached: 
#' \itemize{
#' \item{
#' `ihct_variation_drop` is based on the amount of variation from the 
#' dissimilarity matrix.
#' It triggers a new tree randomization only when the amount of variation  
#' in the group being divided has reached a threshold since last randomization
#' (default: 20% drop in variation). In other words, new randomizations happen
#' only when 
#' tree divisions reach a certain threshold of variation since the last 
#' randomization. For example, when a tree peels off only 1 site at a time, 
#' this argument makes sure no new randomization trigger unless variability
#' reaches the desired threshold.}
#' \item{
#' `ihct_sites_drop`is based on how many sites the group has lost. 
#' It triggers new randomizations only when a certain number of sites have been
#' excluded since the last randomization.}
#' }
#' 
#' These two parameters with their defaults (`ihct_variation_drop = 0.2`, 
#' `ihct_sites_drop = 10`) result in marginal
#' changes in algorithm performance (loss in CCC <0.001) and make the tree 
#' 1.5 faster to build. In our tests, the larger the datasets, the higher
#' the savings with these two optimization parameters. Note, however,
#' that reusing trees only works with UPGMA currently, so it only applies to
#'  `method = "average"`.
#'
#' The two arguments work as a pair, and each of them can be set so that the
#' other no longer has any effect. A new randomization is made as soon as
#' *either* of them asks for one, so whichever of the two asks more often is
#' the one that decides:
#' \itemize{
#' \item{`ihct_sites_drop = 1` (or `0`) means new randomizations every time 
#' a site is treated (so randomizations at every division, and 
#' `ihct_variation_drop` is never used).}
#' \item{`ihct_variation_drop = 0` likewise means new randomizations at 
#' every division, and
#' `ihct_sites_drop` is then never used.}
#' \item{`ihct_variation_drop = 1` never triggers new randomizations, leaving
#' `ihct_sites_drop` to decide on its own, and `ihct_sites_drop = Inf` never
#' triggers randomizations, leaving `ihct_variation_drop` to decide on its own.}
#' \item{Both switched off (`ihct_variation_drop = 1` and
#' `ihct_sites_drop = Inf`) randomizes once, at the first division, and reuses
#' those trees for the whole tree. This is the fastest setting and the one that
#' fits the data least well.}}
#'
#' To reproduce the tree that bioregion 1.4.0 and earlier produced, for the
#' same `seed`, use:
#'
#' \preformatted{
#' hclu_hierarclust(dissimilarity,
#'                  method = "average",
#'                  optimal_tree_method = "ihct",
#'                  n_runs = 100,
#'                  ihct_top_n_trees = 2,
#'                  ihct_variation_drop = 0,
#'                  ihct_height_rule = "max_child")
#' }
#'
#' --> *Height of nodes in the tree*
#' 
#' The height of nodes in a tree must be monotonous, i.e. a child node cannot 
#' be have a higher height than its parents. However, this situation can happen
#' when building the tree, which is why all tree construction algorithms have a
#' monotonicity section where node height is recalculated.
#' 
#' `ihct_height_rule` decides how we do it in IHCT. `"max_child"` raises every
#'  division to the
#' highest division it contains, which is simple but can push a division far
#' above the dissimilarities it summarizes. `"least_squares"` (default) instead
#' moves the heights as little as possible, which with `method = "average"`
#' gives the heights that fit the dissimilarities best on the topology at hand,
#' so the cophenetic correlation is never below the one `"max_child"` gives and
#' is usually above it. 
#' We recommend leaving to default as in our own testing it provides better 
#' performance (highest CCC).
#'
#' --> *Parallelization for quicker computation time*
#'  
#' Most of the waiting is spent building the randomized trees, and the runs of
#' one group do not depend on each other, so `ihct_n_workers` can share them
#' between several processes of your computer. This only pays on large
#' matrices, where a single run is slow enough to be worth sending to another
#' process: groups of fewer than 200 sites are always done in one process, and
#' small datasets should be left at `ihct_n_workers = 1`. We found that
#' 4 workers give good gains (about three times faster on a 5,000-site
#' matrix); beyond that the processes spend their time waiting for memory
#' rather than computing, and on a 10,000-site matrix going from 4 workers
#' to 8 provided only limited gains while doubling the memory needed. 
#' Each worker also needs its own copy of the
#' dissimilarity matrix on Windows, about 200 MB for 5,000 sites and 800 MB
#' for 10,000, so ask for fewer workers than your memory allows copies.
#' Whatever you set, the tree is the same: the random shuffles are always drawn
#' in the same order by the main process, and only the building of the trees is
#' handed out to workers. 
#' 
#' @references
#' Kreft H & Jetz W (2010) A framework for delineating biogeographical regions
#' based on species distributions. \emph{Journal of Biogeography} 37, 2029-2053.
#' 
#' Dapporto L, Ramazzotti M, Fattorini S, Talavera G, Vila R & Dennis, RLH 
#' (2013) Recluster: an unbiased clustering procedure for beta-diversity 
#' turnover. \emph{Ecography} 36, 1070--1075.
#' 
#' Dapporto L, Ciolli G, Dennis RLH, Fox R & Shreeve TG (2015) A new procedure 
#' for extrapolating turnover regionalization at mid-small spatial scales, 
#' tested on British butterflies. \emph{Methods in Ecology and Evolution} 6
#' , 1287--1297. 
#' 
#' @seealso
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [cut_tree] [ihct]
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "Simpson")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim, 
#'                           n_clust = 5)
#' tree1
#' plot(tree1)
#' str(tree1)
#' tree1$clusters
#' 
#' # User-defined height cut
#' # Only one height
#' tree2 <- hclu_hierarclust(dissim, 
#'                           cut_height = .05)
#' tree2
#' tree2$clusters
#' 
#' # Multiple heights
#' tree3 <- hclu_hierarclust(dissim, 
#'                           cut_height = c(.05, .15, .25))
#' 
#' tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#' # Info on each partition can be found in table cluster_info
#' tree3$cluster_info
#' plot(tree3)
#' 
#' @export
hclu_hierarclust <- function(dissimilarity,
                             index = names(dissimilarity)[3],
                             method = "average",
                             randomize = TRUE,
                             seed = NULL,
                             n_runs = 100,
                             keep_trials = "no",
                             optimal_tree_method = "ihct", 
                             n_clust = NULL,
                             cut_height = NULL,
                             find_h = TRUE,
                             h_max = 1,
                             h_min = 0,
                             consensus_p = 0.5,
                             show_hierarchy = FALSE,
                             ihct_top_n_trees = 2,
                             ihct_variation_drop = 0.2,
                             ihct_sites_drop = 10,
                             ihct_height_rule = "least_squares",
                             ihct_n_workers = 1,
                             verbose = TRUE){
  # 1. Controls ---------------------------------------------------------------
  controls(args = NULL, data = dissimilarity, type = "input_nhandhclu")
  if(!inherits(dissimilarity, "dist")){
    controls(args = NULL, data = dissimilarity, type = "input_dissimilarity")
    controls(args = NULL, data = dissimilarity, 
             type = "input_data_frame_nhandhclu")
    controls(args = index, data = dissimilarity, type = "input_net_index")
    net <- dissimilarity
    # Convert tibble into dataframe
    if(inherits(net, "tbl_df")){
      net <- as.data.frame(net)
    }
    colnameindex <- index
    if(is.numeric(colnameindex)){
      colnameindex <- colnames(net)[index]
      if(is.null(colnameindex)){
        colnameindex <- NA
      }
    }
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_value")
    dist_mat <- net_to_mat(net,
                           weight = TRUE, 
                           squared = TRUE, 
                           symmetrical = TRUE)
  } else {
    controls(args = NULL, data = dissimilarity, type = "input_dist")
    dist_mat <- as.matrix(dissimilarity)
    if(is.null(names(dist_mat))){
      rownames(dist_mat) <- paste0(1:dim(dist_mat)[1])
      colnames(dist_mat) <- rownames(dist_mat)
      if(verbose){
        message("No labels detected, they have been assigned automatically.")
      }
    }
    dissimilarity <- mat_to_net(dist_mat, weight = TRUE)
    if(is.null(index)) {
      colnames(dissimilarity)[3] <- "Dissimilarity"
    } else {
      colnames(dissimilarity)[3] <- index
    }
    colnameindex <- NA
  }
  nsites <- dim(dist_mat)[1]

  
  controls(args = method, data = NULL, type = "character")
  if(!(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                          "mcquitty", "median", "centroid" ))){
    stop(paste0("Please choose method from the following:\n",
                "ward.D, ward.D2, single, complete, average, mcquitty, median ", 
                "or centroid"), 
         call. = FALSE)
  }
  controls(args = randomize, data = NULL, type = "boolean")
  if(randomize & !is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  controls(args = n_runs, data = NULL, type = "strict_positive_integer")
  controls(args = ihct_top_n_trees, data = NULL, type = "strict_positive_integer")
  controls(args = ihct_n_workers, data = NULL, type = "strict_positive_integer")
  controls(args = ihct_variation_drop, data = NULL, type = "positive_numeric")
  if(ihct_variation_drop > 1) {
    stop("ihct_variation_drop must be between 0 and 1.",
         call. = FALSE)
  }
  if(!is.numeric(ihct_sites_drop) || length(ihct_sites_drop) != 1 || is.na(ihct_sites_drop) ||
     ihct_sites_drop < 0) {
    stop(paste0("ihct_sites_drop must be a single number of sites, 0 or more ",
                "(0 or 1 to randomize at every division, Inf to leave the ",
                "randomizations to ihct_variation_drop)."),
         call. = FALSE)
  }
  controls(args = ihct_height_rule, data = NULL, type = "character")
  if(!(ihct_height_rule %in% c("least_squares", "max_child"))){
    stop(paste0("Please choose ihct_height_rule from the following:\n",
                "least_squares or max_child"),
         call. = FALSE)
  }
  controls(args = keep_trials, data = NULL, type = "character")
  if(!(keep_trials %in% c("no", "all", "metrics"))){
    stop(paste0("Please choose keep_trials from the following:\n",
                "no, all or metrics"), 
         call. = FALSE)
  }
  controls(args = optimal_tree_method, data = NULL, type = "character")
  if(!(optimal_tree_method %in% c("ihct", "iterative_consensus_tree",
                                  "best", "consensus"))){
    stop(paste0("Please choose optimal_tree_method from the following:\n",
                "ihct, best or consensus"),
    call. = FALSE)
  }
  # "iterative_consensus_tree" is the name used up to bioregion 1.4.0
  if(optimal_tree_method == "iterative_consensus_tree"){
    optimal_tree_method <- "ihct"
  }

  if(!is.null(n_clust)) {
    if(is.numeric(n_clust)) {
        controls(args = n_clust, 
                 data = NULL, 
                 type = "strict_positive_integer_vector")
    } else if(inherits(n_clust, "bioregion.bioregionalization.metrics")){
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop(paste0("n_clust does not have an optimal number of clusters. ",
                    "Did you specify partition_optimisation = TRUE in ",
                    "bioregionalization_metrics()?"), 
             call. = FALSE)
      }
    } else{
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from bioregionalization_metrics()", 
           call. = FALSE)
    }
    if(!is.null(cut_height)){
      stop(paste0("Please provide either n_clust or cut_height, ",
                  "but not both at the same time."), 
           call. = FALSE)
    }
  }
  if(!is.null(cut_height)){
    controls(args = cut_height, data = NULL, type = "positive_numeric_vector")
  }
  controls(args = find_h, data = NULL, type = "boolean")
  if(find_h){
    controls(args = h_min, data = NULL, type = "positive_numeric")
    controls(args = h_max, data = NULL, type = "positive_numeric")
    if(h_min > h_max){
      stop("h_min must be inferior to h_max.",
           call. = FALSE)
    }
  }
  controls(args = consensus_p, data = NULL, type = "positive_numeric")
  if(consensus_p < 0.5 | consensus_p > 1) {
    stop("consensus_p must be between 0.5 and 1.",
         call. = FALSE)
  }
  controls(args = show_hierarchy, data = NULL, type = "boolean")
  controls(args = verbose, data = NULL, type = "boolean")
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "hclu_hierarclust")
  
  # Outputs args
  outputs$args <- list(index = index,
                       method = method,
                       randomize = randomize,
                       seed = seed,
                       n_runs = n_runs,
                       optimal_tree_method = optimal_tree_method,
                       keep_trials = keep_trials,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min,
                       consensus_p = consensus_p,
                       show_hierarchy = show_hierarchy,
                       ihct_top_n_trees = ihct_top_n_trees,
                       ihct_variation_drop = ihct_variation_drop,
                       ihct_sites_drop = ihct_sites_drop,
                       ihct_height_rule = ihct_height_rule,
                       verbose = verbose)
  
  # Determine pairwise_metric and data_type
  pairwise_metric <- ifelse(!inherits(dissimilarity, "dist"), 
                            colnameindex, 
                            NA)
  data_type <- detect_data_type_from_metric(pairwise_metric)
  
  # Outputs inputs
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = pairwise_metric,
                         dissimilarity = TRUE,
                         nb_sites = nsites,
                         data_type = data_type,
                         node_type = "site")
  
  if(randomize) {
    if(optimal_tree_method == "ihct") {
      if(verbose){
        message(paste0("Building the iterative hierarchical consensus tree...",
                       " Note that this",
                       " process can take time especially if you have a lot of",
                       " sites."))
      }

      if(method == "mcquitty") {
        warning("mcquitty (WPGMA) method may not be properly implemented",
                " in Iterative Hierarchical Tree Construction (IHCT), because of the ",
                "hybrid divise-agglomerative nature of IHCT. ",
                "In WPGMA, heights are updated iteratively from bottom to top as ",
                "the tree is constructed. ",
                "In IHCT, divisions are created from top to bottom, based on a ",
                "majority decision among many trees. ", 
                "Hence, it is not possible to exactly compute WPGMA calculations ",
                "with IHCT - final height calculations are approximated into UPGMA.")
      }

      # On a large matrix, most of the computing time goes into the fresh
      # randomizations ihct_sites_drop asks for: the tree peels sites off a few
      # at a time, so that rule fires on groups that still hold nearly every
      # site, over and over. What it adds to the fit of the tree is too small
      # for us to measure at those sizes (tests on 2,000 to 10,000
      # simulated sites, always under 0.001 of cophenetic correlation), so we
      # recommend to switch it off above 2,000 sites
      if(verbose && method == "average" && ihct_variation_drop > 0 &&
         is.finite(ihct_sites_drop) && nrow(dist_mat) > 2000) {
        message(paste0("This dissimilarity matrix has ", nrow(dist_mat),
                       " sites. On matrices this large, setting ",
                       "ihct_sites_drop = Inf builds the tree about two to ",
                       "three times faster, and in our tests the tree fitted ",
                       "the dissimilarities just as well (cophenetic ",
                       "correlation differed by less than 0.001). See the ",
                       "Details section of ?hclu_hierarclust."))
      }

      if (!is.null(seed)) set.seed(seed) # generate seed
      
      consensus_tree <- ihct(dist_mat,
                             method = method,
                             n_runs = n_runs,
                             top_n_trees = ihct_top_n_trees,
                             variation_drop = ihct_variation_drop,
                             sites_drop = ihct_sites_drop,
                             height_rule = ihct_height_rule,
                             n_workers = ihct_n_workers,
                             verbose = verbose)
      
      if (!is.null(seed)) rm(.Random.seed, envir = globalenv()) # remove seed

      # Compute hierarchical tree
      outputs$algorithm$final.tree <- consensus_tree
      
      evals <- tree_eval(consensus_tree,
                         dist_mat)
      
      outputs$algorithm$final.tree.coph.cor <- evals$cophcor
      outputs$algorithm$final.tree.msd <- evals$msd
      
    } else {
      if(verbose){
        message(paste0("Randomizing the dissimilarity matrix with ", 
                       n_runs,
                       " trials"))
      }

      if (!is.null(seed)) set.seed(seed)
      results <- vector("list", n_runs)
      for (run in 1:n_runs) {
        trial <- list()

        trial$dist_mat <- randomize_dist(dist_mat)

        trial$hierartree <- fastcluster::hclust(stats::as.dist(trial$dist_mat), 
                                                method = method)

        evals <- tree_eval(trial$hierartree, trial$dist_mat)
        trial$cophcor <- evals$cophcor
        trial$msd <- evals$msd
        
        if(keep_trials != "all"){
          trial$dist_mat <- 0
        }

        results[[run]] <- trial
      }
      if (!is.null(seed)) rm(.Random.seed, envir = globalenv())
      
      if (optimal_tree_method == "best") {
        coph.coeffs <- sapply(results, function(x) x$cophcor)
        
        if(verbose){
          message(paste0(" -- range of cophenetic correlation coefficients ",
                         "among trials: ",
                         round(min(coph.coeffs), 4), 
                         " - ", 
                         round(max(coph.coeffs), 4)))
        }


        best.run <- which.max(coph.coeffs)
        final.tree <- results[[best.run]]$hierartree
        final.tree.metrics <- results[[best.run]]
        
      } else if (optimal_tree_method == "consensus") {
        if (n_runs < 2) {
          stop("At least two trees are required to calculate a consensus.",
               call. = FALSE)
        }

        trees <- lapply(results, function(trial) ape::as.phylo(trial$hierartree))

        consensus_tree <- ape::consensus(trees, p = consensus_p)
        consensus_tree <- phangorn::nnls.tree(stats::as.dist(dist_mat), 
                                              consensus_tree, 
                                              method = "ultrametric", 
                                              trace = 0)
        # nnls.tree() results in some branches with NEGATIVE branch lengths
        # at about -1e-16 instead of 0, which
        # can put a node above its parent. 
        # We correct the negative branch lengths to 0 before converting to hclust
        consensus_tree$edge.length[consensus_tree$edge.length < 0] <- 0
        consensus_tree <- ape::multi2di(consensus_tree)
        consensus_tree <- ape::as.hclust.phylo(consensus_tree)

        final.tree <- consensus_tree
        evals <- tree_eval(consensus_tree, dist_mat)
        final.tree.metrics <- list(cophcor = evals$cophcor, 
                                   msd = evals$msd)
      }
      
      # keep_trials
      if(keep_trials == "metrics"){
        for (run in 1:n_runs) {
          results[[run]] <- results[[run]][-c(1,2)]
        }
      }

      outputs$algorithm$final.tree <- final.tree
      outputs$algorithm$final.tree.coph.cor <- final.tree.metrics$cophcor
      outputs$algorithm$final.tree.msd <- final.tree.metrics$msd

    }
    if(verbose){
      message(paste0("\nFinal tree has a ",
                     round(outputs$algorithm$final.tree.coph.cor, 4),
                     " cophenetic correlation coefficient with the initial ",
                     "dissimilarity matrix\n"))
    }

  } else  {
    outputs$algorithm$final.tree <- fastcluster::hclust(stats::as.dist(dist_mat),
                                                        method = method)
    
    
    evals <- tree_eval(outputs$algorithm$final.tree,
                       dist_mat)
    
    outputs$algorithm$final.tree.coph.cor <- evals$cophcor
    outputs$algorithm$final.tree.msd <- evals$msd
    
    if(verbose){
      message(paste0("Output tree has a ",
                     round(outputs$algorithm$final.tree.coph.cor, 4),
                     " cophenetic correlation coefficient with the initial
                   dissimilarity matrix\n"))
    }

  }
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  if(any(!is.null(n_clust) | !is.null(cut_height))){
    outputs <- cut_tree(outputs,
                        n_clust = n_clust,
                        cut_height = cut_height,
                        find_h = find_h,
                        h_max = h_max,
                        h_min = h_min,
                        show_hierarchy = show_hierarchy,
                        verbose = verbose)
    
    # Add hierarchical attribute
    outputs$inputs$hierarchical <- ifelse(ncol(outputs$clusters) > 2,
                                          TRUE,
                                          FALSE)
    
    # Add node_type attribute
    attr(outputs$clusters, "node_type") <- rep("site", dim(outputs$clusters)[1])
    
  } else {
    outputs$clusters <- NA
    outputs$cluster_info <- NA
    outputs$inputs$hierarchical <- FALSE
  }
  
  # Keep trials 
  if(randomize){
    if(keep_trials == "no"){
      outputs$algorithm$trials <- "Trials not stored in output"
    }else{
      if(optimal_tree_method == "ihct"){
        outputs$algorithm$trials <- "Trials not stored in output"
      }else{
        outputs$algorithm$trials <- results
      }
    }
  }
  
  return(outputs)
}

