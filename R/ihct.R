# Iterative Hierarchical Consensus Tree (IHCT)
#
# What the algorithm does, what every argument is for, and the measurements
# behind the defaults are written for users in the roxygen block below
# (?ihct). This header only records what a reader of this file needs and the
# help page does not say.
#
# Implementation notes (2026-09):
# - Sites are handled as integer positions in dist_mat; names are only used to
#   sort sites (the co-assignment step works on sites sorted by name, as we
#   originally implementated) and to label the final tree.
# - The groups still to be divided are kept in a work list instead of using R
#   recursion, so that deep chain-like trees do not exhaust R's call stack.
# - A candidate tree is a list(merge, height, pairs, leaf_site, score): the
#   merge matrix and the heights of an hclust object, the number of site pairs
#   joined at each node, the site of dist_mat every leaf stands for, and the
#   score it was ranked on. Fresh and inherited trees have the same form, so
#   the rest of the code does not have to know where a tree comes from.
# - For UPGMA, the cophenetic correlation of a tree is obtained from the merge
#   sizes and heights only (see tree_fit_score), which avoids computing the
#   full cophenetic matrix of every randomized tree.
# - Pruning a tree gives every node that lost sites the mean dissimilarity
#   between the two groups it still joins, which is the height UPGMA gives it;
#   the other linkage methods define their heights otherwise, which is why
#   trees are only ever passed down with method = "average".
# - Either reuse criterion at its "rebuild always" value makes the other idle
#   and leaves nothing to inherit, so pass_trees_down is FALSE for either of
#   them and the bookkeeping the reuse needs is skipped altogether.
# - The randomized runs of a large group can be shared between worker processes
#   (n_workers). The permutations are always drawn on this process, so the tree
#   is the same however many workers there are (see start_workers).

#' Iterative hierarchical consensus tree from a dissimilarity matrix
#'
#' @description
#' This function builds a hierarchical tree of sites from a matrix of pairwise
#' dissimilarities, using the Iterative Hierarchical Consensus Tree (IHCT)
#' algorithm. Unlike an ordinary hierarchical clustering, the tree does not
#' depend on the order in which the sites happen to be stored in the matrix.
#'
#' Most users should use [hclu_hierarclust()] and not this function. 
#' `ihct()` is an algorithm used in [hclu_hierarclust()], exported for
#' users who want the algorithm itself for their purposes.
#'
#' @param dist_mat A square `matrix` of dissimilarities between sites, with the
#' site names as row and column names. A `dist` object is also accepted and
#' converted.
#'
#' @param method The name of the linkage method, as in
#' [hclust][fastcluster::hclust]. It is used both to build the randomized trees
#' and to give every division of the final tree its height. Should be one of
#' `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`
#' (= UPGMA, the default), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC) or
#' `"centroid"` (= UPGMC). Two of the optimization parameters described in
#' details are only
#' available with some of these.
#'
#' @param n_runs The number of randomized trees built to decide one division
#' (`100` by default). More trees make each division more stable and the whole
#' tree slower to obtain.
#'
#' @param top_n_trees An `integer` indicating how many of the randomized trees,
#' ranked by how well they fit the dissimilarities, decide a division (`2` by
#' default). With `1` the division is read straight off the best tree; with
#' more, sites are grouped according to how often they end up on the same side
#' across those trees. See details.
#'
#' @param variation_drop A `numeric` value between 0 and 1 (`0.2` by default).
#' This is a speed optimization parameter: it lets a group of sites reuse the 
#' randomized trees of
#' the group it was split from instead of building its own, which is faster and
#' costs a very small amount of fit. Its rule is based on variation: fresh
#' trees are built once a group has lost this share of the variation it still
#' held when its current trees were built. Set it to `0` to build fresh trees
#' at every division, or to `1` to switch this rule off and leave `sites_drop`
#' in charge. Only used with `method = "average"`. See Details.
#'
#' @param sites_drop A `numeric` value of 0 or more (`10` by default). This is
#' the second speed optimization parameter, working alongside `variation_drop`,
#'  with a rule based
#' on sites rather than on variation: fresh trees are also built once a group
#' has lost this many sites since its current trees were built. Set it to `0`
#' or `1` to build fresh trees at every division, or to `Inf` to switch this
#' rule off and leave `variation_drop` in charge. Only used with
#' `method = "average"` and `variation_drop > 0`. See Details.
#'
#' @param height_rule A `character` string indicating how the heights of the
#' tree are corrected when a division comes out lower than a division it
#' contains. With `"least_squares"` (default) the heights are moved as little
#' as possible, which fits the dissimilarities better; with `"max_child"` every
#' division is raised to the highest division it contains, as in bioregion
#' 1.4.0 and earlier. See Details.
#'
#' @param tie_block_resolution A `boolean` (`TRUE` by default) deciding what
#' happens to a group of sites whose dissimilarities are all equal,
#' where the randomization has nothing left to choose between. With
#' `TRUE` such a group is resolved directly at that value, which saves its
#' `n_runs` trees and costs nothing in fit; with `FALSE` it is divided from
#' randomized trees like any other group. Only used with `method` set to
#' `"single"`, `"complete"`, `"average"` or `"mcquitty"`. See Details.
#'
#' @param n_workers An `integer` of 1 or more indicating how many processes of
#' your computer build the randomized trees at the same time. With `1`
#' (default) they are built one after another. Higher values are worth it on
#' large matrices only, and they do not change the tree: the same seed gives
#' the same result whatever this is set to. See Details.
#'
#' @param size_parallel An `integer` of 2 or more indicating the smallest group
#' of sites whose randomized trees are worth handing out to the worker
#' processes (`200` by default). Groups smaller than this are always done in
#' one process, because sending them out costs more than it saves. Only used
#' when `n_workers > 1`.
#'
#' @param verbose A `boolean` indicating whether to display progress and
#' information messages. Set to `FALSE` to suppress them.
#'
#' @return
#' An `hclust` object, whose `method` element is set to
#' `"Iterative Hierarchical Consensus Tree"`. Its heights are monotone (no
#' division is lower than a division it contains) and its labels are the site
#' names of `dist_mat`, in alphabetical order. It can be plotted, cut with
#' [cut_tree()] or [cutree][stats::cutree], and used anywhere an `hclust`
#' object is expected.
#'
#' @details
#' ## --- 1. The problem: the order of the sites changes the tree ---
#'
#' A dissimilarity matrix computed from species composition contains a great
#' many identical values, because many pairs of sites share exactly the same
#' number of species. Hierarchical clustering has to break those ties somehow,
#' and [hclust][fastcluster::hclust] breaks them by taking whichever pair comes
#' first in the matrix. Shuffle the rows of the matrix and you get a different
#' tree out of the very same data (Dapporto et al., 2013). IHCT removes that
#' arbitrariness by generating a true consensus tree from many randomizations.
#'
#' ## --- 2. The tree is built from the top down ---
#' 
#' IHCT does not join sites into ever larger groups, the way ordinary
#' hierarchical clustering does. It starts with a single group containing
#' every site, cuts
#' that group in two, then cuts each of the two in two, and carries on until
#' every group holds a single site. Each of those cuts is called a division,
#' and the whole algorithm is about deciding one division well. Divisions
#' are always strictly binary (i.e., one branch is always cut into two 
#' branches).
#'
#' ## --- 3. How one division is decided ---
#'
#' Take a group of sites that has to be divided. Four things happen to it:
#'
#' \enumerate{
#' \item Its sites are put in a random order and a tree is built on them, with
#' [hclust][fastcluster::hclust] and the linkage of `method`. This is repeated
#' `n_runs` times, giving `n_runs` trees that differ only in how the ties were
#' broken.
#' \item Each of those trees is scored by its cophenetic correlation 
#' coefficient, that is,
#' by how closely the distances read off the tree reproduce the dissimilarities
#' the tree was built from. The `top_n_trees` best-scoring trees are kept and
#' the rest are discarded.
#' \item From the kept trees, we take the top division to decide how to split
#' the group of sites into two. With
#' `top_n_trees = 1` we take the top division of the best tree. With more, 
#' we take a consensus division based on how often sites are on the same side
#' across the kept trees. The rule we used is that the consensus is based on 
#' the majority decision, sites are grouped together if theyr are together
#' in > 50% of trees. 
#' \item The division is given a height, computed from the dissimilarities
#' between the two groups it separates in the way `method` prescribes: their
#' mean for `"average"` (UPGMA), their smallest value for `"single"`, their
#' largest for `"complete"`, and so on.
#' }
#'
#' For each of the two groups that come out of the division, the algorithm
#' run them through the same four
#' steps, and so on for their own halves, until we reach single sites.
#'
#' This why the algorithm is called "iterative": it randomizes trees again
#' at every division from top to bottom. It is why this method provides better
#' performance compareds to other approaches: a tie broken one way at the top
#'  of the tree does not force the same
#' choice at the divisions underneath it.
#'
#' Since the randomization uses R's random number generator, call
#' [set.seed()] before `ihct()` if you want the same tree twice.
#'
#' ## --- 4. Computaion time and how to shorten it ---
#'
#' Step 1 above is where nearly all of the computing time goes: `n_runs` trees
#' are built at every single division. However, the trees usually peel a handful of
#' sites off the rest at each division rather than splitting the data in half.
#' This means randomizing everything at every division can be inefficient, resulting
#' in large computation times, making the function too long to use on large datasets.
#'
#' To make the function usable
#' on large datasets, we provide two optimization parameters. These parameters
#'  make the
#' algorithm reuse previously randomized trees at new division, unless a 
#' threshold of change is reached: 
#' \itemize{
#' \item{
#' `variation_drop` is based on the amount of variation from the 
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
#' `sites_drop`is based on how many sites the group has lost. 
#' It triggers new randomizations only when a certain number of sites have been
#' excluded since the last randomization.}
#' }
#' 
#' These two parameters with their defaults (`variation_drop = 0.2`, 
#' `sites_drop = 10`) result in marginal
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
#' \item{`sites_drop = 1` (or `0`) means new randomizations every time 
#' a site is treated (so randomizations at every division, and 
#' `variation_drop` is never used).}
#' \item{`variation_drop = 0` likewise means new randomizations at 
#' every division, and
#' `sites_drop` is then never used.}
#' \item{`variation_drop = 1` never triggers new randomizations, leaving
#' `sites_drop` to decide on its own, and `sites_drop = Inf` never
#' triggers randomizations, leaving `variation_drop` to decide on its own.}
#' \item{Both switched off (`variation_drop = 1` and
#' `sites_drop = Inf`) randomizes once, at the first division, and reuses
#' those trees for the whole tree. This is the fastest setting and the one that
#' fits the data least well.}}
#' 
#' ## --- 5. Groups where all distances are equal ---
#'
#' Some groups have all their dissimilarities equal to one single value: every
#' pair of sites inside the group is exactly as different as every other pair.
#' These are tied blocks, and they are common in presence-absence data, where
#' indices such as Simpson saturate at 1 as soon as two sites share no species.
#' A tied block means randomization brings nothing useful
#' and the `n_runs` trees are wasted computation time.
#'
#' To avoid this, the argument `tie_block_resolution = TRUE` (the default)
#' recognises such a group and
#' resolves it directly, peeling its sites off one at a time with every
#' division sitting at the common value. Setting it to `FALSE` divides tied 
#' blocks from randomized trees
#' like any other group (pre-1.4.0 behaviour).
#'
#'
#' ## --- 6. Node heights ---
#'
#' The height of nodes in a tree must be monotonous, i.e. a child node cannot
#'  be have a higher height than its parents. However, this situation can 
#' happen when building the tree, which is why all tree construction algorithms
#'  have a monotonicity section where node height is recalculated.
#' 
#' IHCT corrects the
#' heights once the whole tree is built, with a method that depends on 
#' `height_rule`:
#' 
#' `"max_child"` raises every division to the highest division inside it. It is
#' simple, and it is what bioregion 1.4.0 and earlier did, but a single high
#' division buried deep inside a group drags all of its parents up with it,
#' well above the dissimilarities those divisions actually summarize. We
#' found out it provides lower quality (lower cophenetic correlation 
#' coefficient) than `"least_squares"`, so we changed the default after 1.4.0.
#'
#' `"least_squares"` (the default) instead moves the heights as little as it
#' can: among all sets of heights with no branch doubling back, it takes the
#' one that stays closest to the heights the divisions were given, weighting
#' each division by the number of site pairs it stands for. With
#' `method = "average"` these are provably the heights that fit the
#' dissimilarities best on the tree shape at hand, so the cophenetic
#' correlation is never below what `"max_child"` gives and is usually a little
#' above it.
#'
#' ## --- 7. Computation time and parallelization ---
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
#' ## --- 8. Reproducing the trees of bioregion 1.4.0 ---
#'
#' For the same seed, the following settings give the tree that bioregion 1.4.0
#' and earlier produced:
#'
#' \preformatted{
#' ihct(dist_mat,
#'      method = "average",
#'      n_runs = 100,
#'      top_n_trees = 2,
#'      variation_drop = 0,
#'      height_rule = "max_child",
#'      tie_block_resolution = FALSE)
#' }
#'
#' `sites_drop` may be left at any value here, since `variation_drop = 0`
#' already rebuilds the trees at every division.
#'
#' @references
#' Dapporto L, Ramazzotti M, Fattorini S, Talavera G, Vila R & Dennis RLH
#' (2013) Recluster: an unbiased clustering procedure for beta-diversity
#' turnover. \emph{Ecography} 36, 1070--1075.
#'
#' Dapporto L, Ciolli G, Dennis RLH, Fox R & Shreeve TG (2015) A new procedure
#' for extrapolating turnover regionalization at mid-small spatial scales,
#' tested on British butterflies. \emph{Methods in Ecology and Evolution} 6,
#' 1287--1297.
#'
#' Kreft H & Jetz W (2010) A framework for delineating biogeographical regions
#' based on species distributions. \emph{Journal of Biogeography} 37,
#' 2029--2053.
#'
#' @seealso
#' For more details illustrated with a practical example,
#' see the vignette:
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#'
#' Associated functions:
#' [hclu_hierarclust] [cut_tree]
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
#' dist_mat <- net_to_mat(dissim[, 1:3],
#'                        weight = TRUE,
#'                        squared = TRUE,
#'                        symmetrical = TRUE)
#'
#' set.seed(1)
#' tree <- ihct(dist_mat,
#'              n_runs = 20,
#'              verbose = FALSE)
#' plot(tree)
#' cut_tree(tree, n_clust = 3)
#'
#' # Fastest setting: build the randomized trees once, at the first division,
#' # and reuse them all the way down
#' set.seed(1)
#' fast_tree <- ihct(dist_mat,
#'                   n_runs = 20,
#'                   variation_drop = 1,
#'                   sites_drop = Inf,
#'                   verbose = FALSE)
#'
#' @export
ihct <- function(dist_mat,
                 method = "average",
                 n_runs = 100,
                 top_n_trees = 2,
                 variation_drop = 0.2,
                 sites_drop = 10,
                 height_rule = c("least_squares", "max_child"),
                 tie_block_resolution = TRUE,
                 n_workers = 1,
                 size_parallel = 200,
                 verbose = TRUE) {

  # checked here rather than where it is used, at the end of the run
  height_rule <- match.arg(height_rule)
  if (inherits(dist_mat, "dist")) dist_mat <- as.matrix(dist_mat)
  if (!is.matrix(dist_mat) || !is.numeric(dist_mat) ||
      nrow(dist_mat) != ncol(dist_mat) || nrow(dist_mat) == 0) {
    stop("dist_mat must be a square numeric matrix of dissimilarities between ",
         "sites, or a dist object.", call. = FALSE)
  }
  n <- nrow(dist_mat)
  site_names <- rownames(dist_mat)
  if (is.null(site_names)) {
    site_names <- as.character(seq_len(n))
    rownames(dist_mat) <- site_names
    colnames(dist_mat) <- site_names
    if (verbose) {
      message("No labels detected, they have been assigned automatically.")
    }
  }
  if (anyDuplicated(site_names) > 0) {
    stop("The row names of dist_mat must be unique site names.", call. = FALSE)
  }
  if (top_n_trees > n_runs) top_n_trees <- n_runs
  if (n == 1) stop("At least two sites are needed to build a tree.")
  if (length(variation_drop) != 1 || is.na(variation_drop) ||
      variation_drop < 0 || variation_drop > 1) {
    stop("variation_drop must be a single number between 0 and 1.")
  }
  if (length(sites_drop) != 1 || is.na(sites_drop) || sites_drop < 0) {
    stop("sites_drop must be a single number of sites, 0 or more (0 or 1 to ",
         "rebuild the trees at every division, Inf to leave the rebuilds to ",
         "variation_drop).")
  }
  if (length(tie_block_resolution) != 1 || !is.logical(tie_block_resolution) ||
      is.na(tie_block_resolution)) {
    stop("tie_block_resolution must be TRUE or FALSE.")
  }
  if (length(n_workers) != 1 || is.na(n_workers) || !is.numeric(n_workers) ||
      n_workers < 1) {
    stop("n_workers must be a single number of processes, 1 or more.")
  }
  if (length(size_parallel) != 1 || is.na(size_parallel) ||
      !is.numeric(size_parallel) || size_parallel < 2) {
    stop("size_parallel must be a single number of sites, 2 or more.")
  }
  n_workers <- as.integer(n_workers)
  # Pruning a tree gives every node that lost sites the mean dissimilarity
  # between the two groups it still joins, which is the height UPGMA gives it;
  # the other linkage methods define their heights otherwise, so with them the
  # trees are built afresh at every group.
  # Either criterion at its "rebuild always" value makes the other idle and
  # leaves nothing to inherit, so the bookkeeping the reuse needs is skipped
  # altogether. Every division loses at least one site, so sites_drop of 1 or
  # less fires everywhere just as variation_drop = 0 does.
  pass_trees_down <- variation_drop > 0 && sites_drop > 1 && method == "average"
  if (variation_drop > 0 && sites_drop > 1 && !pass_trees_down && verbose) {
    message("Randomized trees are rebuilt at every division: reusing them ",
            "is only available with method = \"average\".")
  }
  # A tied block is resolved on its common dissimilarity, which is the height
  # min, max and mean all give it; the centroid-based linkages define their
  # heights otherwise, so they build trees on tied blocks like any other group.
  resolve_ties <- tie_block_resolution &&
    method %in% c("single", "complete", "average", "mcquitty")
  if (tie_block_resolution && !resolve_ties && verbose) {
    message("Groups whose dissimilarities are all equal are divided from ",
            "randomized trees like any other: resolving them directly is ",
            "only available with method = \"single\", \"complete\", ",
            "\"average\" or \"mcquitty\".")
  }
  keep_share <- 1 - variation_drop
  # A group gets fresh trees when either rule says so: when it has lost
  # `variation_drop` of the share of the variation its trees were built for,
  # or when it has lost `sites_drop` sites since then (a count of sites; Inf
  # switches the second rule off). The two catch different shapes of tree: the
  # share falls
  # slowly at the top of a tree that peels one site off at a time, where the
  # count then forces the rebuild, and quickly on a balanced tree, where the
  # count rarely fires before the share does.
  needs_fresh <- function(m, share, mark_share, mark_sites) {
    share <= keep_share * mark_share || (mark_sites - m) >= sites_drop
  }

  # Worker processes for the randomized runs of the large groups, if asked for
  # and if they can be had; NULL means every run is made here, one at a time.
  workers <- start_workers(n_workers, dist_mat, size_parallel, verbose)
  on.exit(stop_workers(workers), add = TRUE)

  # --- The tree under construction -------------------------------------------
  # One row per node (groups of sites and single sites alike), at most 2n - 1.
  n_max <- 2 * n - 1
  node_height <- numeric(n_max)        # height of the division (0 for sites)
  node_left <- integer(n_max)          # child nodes (0 for sites)
  node_right <- integer(n_max)
  node_site <- integer(n_max)          # site number for single-site nodes, 0 otherwise
  node_pairs <- numeric(n_max)         # site pairs joined at the division (0 for sites)
  node_sites <- vector("list", n_max)  # sites of groups still to be divided
  node_trees <- vector("list", n_max)  # candidate trees a group takes from its parent
  node_variation <- numeric(n_max)     # share of total variation a group carries
  node_mark <- numeric(n_max)          # the share when its trees were last built
  node_mark_n <- numeric(n_max)        # the sites of the group they were built for
  n_nodes <- 1L
  node_sites[[1]] <- seq_len(n)        # the root holds every site, in matrix order
  node_variation[1] <- 1               # every pair is still inside it
  node_mark[1] <- 1
  node_mark_n[1] <- n
  work <- 1L                           # groups waiting to be divided (last in, first out)
  n_placed <- 0L                       # progress: sites turned into leaves

  # the mean and the total sum of squares of the dissimilarities, read one row
  # at a time so that no copy of the matrix is made
  grand_mean <- 0; total_variation <- 0
  if (pass_trees_down) {
    sum_d <- 0; sum_d2 <- 0
    for (i in seq_len(n)) {
      row <- dist_mat[i, ]
      sum_d <- sum_d + sum(row); sum_d2 <- sum_d2 + sum(row * row)
    }
    n_pairs <- n * (n - 1) / 2
    grand_mean <- sum_d / 2 / n_pairs
    total_variation <- sum_d2 / 2 - n_pairs * grand_mean^2
    if (total_variation <= 0) pass_trees_down <- FALSE   # every pair is equal
  }

  add_node <- function() {
    n_nodes <<- n_nodes + 1L
    n_nodes
  }
  new_node <- function(sites) {
    id <- add_node()
    node_sites[[id]] <<- sites
    id
  }
  # A tied block: all the dissimilarities inside the group are the same value,
  # so every tree on it is equally good and there is nothing to decide. The
  # sites are peeled off one by one into a rake, every division at that
  # value, which reproduces the block's dissimilarities exactly.
  resolve_tied_block <- function(id, sites, height) {
    current <- id
    repeat {
      m <- length(sites)
      node_height[current] <<- height
      node_pairs[current] <<- m - 1
      first <- add_node(); node_site[first] <<- sites[1]
      rest <- add_node()
      node_left[current] <<- first
      node_right[current] <<- rest
      if (m == 2) {
        node_site[rest] <<- sites[2]
        break
      }
      sites <- sites[-1]
      current <- rest
    }
  }

  while (length(work) > 0) {
    id <- work[length(work)]
    work <- work[-length(work)]
    sites <- node_sites[[id]]
    node_sites[id] <- list(NULL)
    trees <- node_trees[[id]]          # taken from the parent, NULL at the root
    node_trees[id] <- list(NULL)
    m <- length(sites)

    if (m == 1) {
      node_site[id] <- sites
      n_placed <- n_placed + 1L
      next
    }

    if (m == 2) {
      groups <- list(sites[1], sites[2])
    } else {
      # fresh trees once the group has lost `variation_drop` of the variation
      # the last group they were built for carried
      build_fresh <- is.null(trees) ||
        needs_fresh(m, node_variation[id], node_mark[id], node_mark_n[id])
      # A tied block has nothing to decide. Looking costs a pass over its
      # sub-matrix, which is nothing next to building trees but not next to
      # reusing them; when trees are reused the check is therefore only made
      # if they are flat, which they have to be on a tied block (their heights
      # are means of its dissimilarities).
      if (resolve_ties &&
          (build_fresh || diff(range(trees[[1]]$height)) == 0) &&
          is_tied_block(dist_mat, sites)) {
        resolve_tied_block(id, sites, dist_mat[sites[2], sites[1]])
        n_placed <- n_placed + m
        next
      }
      division <- divide_sites(dist_mat, sites, site_names, method, n_runs,
                               top_n_trees, if (build_fresh) NULL else trees,
                               workers)
      groups <- division$groups
      trees <- division$trees
    }
    # smaller group first (same convention as hclust)
    if (length(groups[[1]]) > length(groups[[2]])) groups <- rev(groups)

    node_height[id] <- division_height(dist_mat, groups[[1]], groups[[2]], sites, method)
    node_pairs[id] <- length(groups[[1]]) * length(groups[[2]])
    node_left[id] <- new_node(groups[[1]])
    node_right[id] <- new_node(groups[[2]])

    # A group takes its parent's trees unless it still carries enough of the
    # variation to be worth fresh ones. They are pruned here, before the
    # parent's are released, and the smaller group is divided first, so that
    # only one large set of trees is ever waiting in the work list.
    if (pass_trees_down && max(lengths(groups)) > 2) {
      node_variation[c(node_left[id], node_right[id])] <-
        variation_share(dist_mat, sites, groups[[1]], groups[[2]],
                        node_variation[id], grand_mean, total_variation)
      # the mark the children are measured against is where these trees were
      # built, which is here when they are fresh and further up when they are not
      mark <- if (build_fresh) node_variation[id] else node_mark[id]
      mark_n <- if (build_fresh) m else node_mark_n[id]
      node_mark[c(node_left[id], node_right[id])] <- mark
      node_mark_n[c(node_left[id], node_right[id])] <- mark_n
      for (child in c(node_left[id], node_right[id])) {
        child_sites <- node_sites[[child]]
        if (length(child_sites) > 2 && !is.null(trees) &&
            !needs_fresh(length(child_sites), node_variation[child], mark, mark_n)) {
          keep_sites <- logical(n)
          keep_sites[child_sites] <- TRUE
          node_trees[[child]] <- prune_trees(trees, keep_sites, dist_mat)
        }
      }
    }
    trees <- NULL
    # the smaller group is divided first: put it last in the work list
    work <- c(work, node_right[id], node_left[id])

    if (verbose && interactive()) {
      cat(sprintf("\rIHCT: %d of %d sites placed, current group: %d sites     ",
                  n_placed, n, m))
      utils::flush.console()
    }
  }
  if (verbose && interactive()) cat("\n")

  # --- Monotone heights ------------------------------------------------------
  keep <- seq_len(n_nodes)
  height <- make_heights_monotone(node_height[keep], node_left[keep],
                                  node_right[keep], node_pairs[keep],
                                  rule = height_rule)

  nodes_to_hclust(height, node_left[keep], node_right[keep], node_site[keep],
                  site_names)
}

# The share of the dataset's total variation that each of the two groups of a
# division still carries.
#
# Two sites are separated at exactly one division of the tree, so the pairs
# inside a group of sites are the business of that group's own sub-tree and of
# no other part of the tree. The share of the total sum of squares those pairs
# carry,
#     sum over pairs inside the group of (d - mean of all d)^2 / total,
# is therefore a ceiling on how much fit the divisions inside the group can
# still win or lose: it is 1 for the whole set of sites, and falls as the tree
# is built. It is what decides where fresh randomized trees are worth building.
#
# The three pieces of a division add up to the share of the group being
# divided: what is inside the first group, what is inside the second, and what
# lies between them. Only the smaller group and the pairs between the two are
# read, and the larger group keeps the rest, so a division costs the size of
# its smaller group times the size of the group -- the same order as the
# dissimilarities already read to give the division its height.
variation_share <- function(dist_mat, sites, small, large, share, mean_d, total) {
  squared_sum <- function(block) sum((block - mean_d)^2)
  # the diagonal of the block below holds the distance of a site to itself
  between <- squared_sum(dist_mat[small, large, drop = FALSE])
  rows <- squared_sum(dist_mat[small, sites, drop = FALSE]) -
    length(small) * mean_d^2
  inside_small <- (rows - between) / 2
  small_share <- max(inside_small / total, 0)
  c(small_share, max(share - small_share - between / total, 0))
}

# TRUE when a group of sites is a tied block, i.e. when all the dissimilarities
# inside it are the same value. The rows are read one at a time and the first
# one that disagrees ends the check, so nothing of the size of the sub-matrix
# is ever built.
is_tied_block <- function(dist_mat, sites) {
  value <- dist_mat[sites[2], sites[1]]
  for (i in seq_len(length(sites) - 1L)) {
    if (any(dist_mat[sites[i], sites[-seq_len(i)]] != value)) return(FALSE)
  }
  TRUE
}

# One randomized tree: the group's dissimilarities with its sites in the order
# `shuffled` gives, clustered, and scored by its fit to the dissimilarities.
#
# The shuffled dissimilarities are read straight out of `dist_mat` into the
# vector fastcluster expects (ihct_shuffled_dist), and the score reads out of
# `dist_mat` too, so the square sub-matrix of the group is never built.
one_tree <- function(shuffled, dist_mat, method) {
  tree <- fastcluster::hclust(ihct_shuffled_dist(dist_mat, shuffled),
                              method = method)
  list(merge = tree$merge, height = tree$height,
       pairs = ihct_node_sizes(tree$merge)$pairs,
       leaf_site = shuffled,
       score = tree_fit_score(tree, dist_mat, method, shuffled))
}

# n_runs trees for a group of sites (integer positions), each built on the
# group's sites shuffled at random.
#
# The shuffles are drawn here, one after another, exactly as a single process
# would draw them, and only the building of the trees is handed to the workers.
# That is what makes the result independent of `n_workers`: the same seed gives
# the same shuffles, hence the same trees, hence the same tree in the end.
fresh_trees <- function(dist_mat, sites, method, n_runs, workers = NULL) {
  m <- length(sites)
  shuffles <- lapply(seq_len(n_runs), function(run) sites[sample.int(m)])
  if (is.null(workers) || m < workers$size_parallel) {
    return(lapply(shuffles, one_tree, dist_mat = dist_mat, method = method))
  }
  # A worker must not disturb the random numbers of this process: the next
  # group has to be drawn where this one left off, whoever built the trees.
  seed_before <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  trees <- if (workers$kind == "fork") {
    # a forked worker already sees this process's memory, matrix included
    parallel::mclapply(shuffles, one_tree, dist_mat = dist_mat, method = method,
                       mc.cores = workers$n)
  } else {
    # a worker of the socket pool was given the matrix once, when the pool was
    # started; only the site numbers of each run travel now
    parallel::parLapply(workers$cl, shuffles, one_tree_on_worker, method = method)
  }
  if (!is.null(seed_before)) assign(".Random.seed", seed_before, globalenv())
  # a worker that failed returns the error instead of a tree
  failed <- vapply(trees, inherits, TRUE, what = "try-error")
  if (any(failed)) {
    stop("Building the randomized trees on ", workers$n, " processes failed: ",
         conditionMessage(attr(trees[[which(failed)[1]]], "condition")),
         "\nRun again with n_workers = 1.", call. = FALSE)
  }
  trees
}

# What a worker of the socket pool runs. It lives here, at the top level of the
# package, on purpose: a function defined inside fresh_trees() would carry that
# call's variables -- the whole dissimilarity matrix among them -- to the
# workers at every group, which is the cost the pool exists to avoid. Written
# this way, only the function's name travels, and the matrix is the copy the
# worker was given once (see start_workers).
one_tree_on_worker <- function(shuffled, method) {
  one_tree(shuffled, get(".ihct_dist_mat", envir = globalenv()), method)
}

# A pool of worker processes for the randomized runs, or NULL when they are all
# to be made on this process.
#
# The two families of operating systems need different pools. On Unix a worker
# is a fork of this process: it already sees the dissimilarity matrix and
# nothing has to be set up. On Windows there is no fork, so the workers are
# fresh R processes, and each is given the matrix once, here, rather than at
# every group -- which is why a large matrix costs as much memory again per
# worker.
#
# Anything that goes wrong (the parallel package missing, sockets refused,
# which happens on managed installations, the package not loadable in a worker)
# is reported once and the runs are simply made here instead. The trees do not
# depend on this, only the time they take.
start_workers <- function(n_workers, dist_mat, size_parallel, verbose) {
  if (n_workers <= 1) return(NULL)
  give_up <- function(why) {
    warning("The randomized trees are built one after another on this ",
            "process: ", why, " Set n_workers = 1 to silence this.",
            call. = FALSE)
    NULL
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    return(give_up("the parallel package is not installed."))
  }
  if (.Platform$OS.type != "windows") {
    return(list(kind = "fork", n = n_workers, cl = NULL,
                size_parallel = size_parallel))
  }
  cl <- tryCatch(parallel::makePSOCKcluster(n_workers), error = function(e) e)
  if (inherits(cl, "error")) {
    return(give_up(paste0("the worker processes could not be started (",
                          conditionMessage(cl), ").")))
  }
  ready <- tryCatch({
    parallel::clusterEvalQ(cl, library(bioregion))
    # the matrix goes to the workers once, under a name of our own
    holder <- new.env(parent = emptyenv())
    assign(".ihct_dist_mat", dist_mat, envir = holder)
    parallel::clusterExport(cl, ".ihct_dist_mat", envir = holder)
    TRUE
  }, error = function(e) e)
  if (inherits(ready, "error")) {
    parallel::stopCluster(cl)
    return(give_up(paste0("the workers could not be prepared (",
                          conditionMessage(ready), ").")))
  }
  if (verbose) {
    message("Randomized trees for groups of ", size_parallel, " sites or more ",
            "are built on ", n_workers, " processes.")
  }
  list(kind = "socket", n = n_workers, cl = cl, size_parallel = size_parallel)
}

stop_workers <- function(workers) {
  if (!is.null(workers) && workers$kind == "socket") {
    try(parallel::stopCluster(workers$cl), silent = TRUE)
  }
  invisible(NULL)
}

# Pass a group's trees down to one of its two halves: every tree is pruned to
# the sites that half keeps (`keep` is TRUE for them, over all the sites of
# dist_mat) and scored again. Pruning only moves the heights of the nodes that
# lost sites, and those heights stay mean dissimilarities between the two
# groups a node joins, so the score is the same sum over nodes as for a fresh
# UPGMA tree (see tree_fit_score).
prune_trees <- function(trees, keep, dist_mat) {
  lapply(trees, function(tree) {
    pruned <- ihct_prune_tree(tree$merge, tree$height, tree$pairs,
                              tree$leaf_site, keep, dist_mat)
    pruned$score <- sum(pruned$pairs * pruned$height^2)
    pruned
  })
}

# Divide a group of sites (integer positions) in two, from n_runs candidate
# trees: the ones given in `trees`, or fresh ones when it is NULL, built on the
# worker processes of `workers` if the group is large enough for that to pay.
# Returns the two groups, each sorted by site name, and the trees the division
# was read from, which the two groups may inherit.
#
# The random draw of a division happens here, inside this call, and not in the
# loop that decides when to rebuild: that is what lets a caller record the
# state of the random numbers before a division and replay it afterwards, which
# is how the tests recover which randomized trees a division was read from.
divide_sites <- function(dist_mat, sites, site_names, method, n_runs,
                         top_n_trees, trees = NULL, workers = NULL) {
  m <- length(sites)
  if (is.null(trees)) {
    trees <- fresh_trees(dist_mat, sites, method, n_runs, workers)
  }

  # 1. keep the best trees (in case of equal scores, the first runs)
  scores <- vapply(trees, function(tree) tree$score, 0)
  best <- rank_by_score(scores)[seq_len(top_n_trees)]

  # 2. two groups from the best trees, working on sites sorted by name. The
  #    division a tree stands for is the one at its root, the two groups it
  #    joins last.
  sorted <- sites[order(site_names[sites])]
  side <- vapply(trees[best], function(tree) {
    as.integer(sorted %in% ihct_top_division(tree$merge, tree$leaf_site)[[1]])
  }, integer(m))
  if (top_n_trees == 1) {
    membership <- side[, 1]
  } else {
    # co-assignment: how often two sites fall on the same side, turned into a
    # dissimilarity, clustered with the same linkage and cut in two. Two sites
    # agree in a tree when they are both in its first group or both in its
    # second, so counting the agreements is two products of the 0/1 matrix of
    # the sides with itself.
    together <- tcrossprod(side) + tcrossprod(1 - side)
    coassign <- 1 - together / n_runs
    dimnames(coassign) <- list(site_names[sorted], site_names[sorted])
    membership <- stats::cutree(fastcluster::hclust(stats::as.dist(coassign), method = method), k = 2)
  }
  # group 1 is the one containing the first site in name order (cutree convention)
  list(groups = list(sorted[membership == membership[1]],
                     sorted[membership != membership[1]]),
       trees = trees)
}

# Score of a tree = its cophenetic correlation with the dissimilarities, or a
# quantity that ranks trees in the same order.
#
# For UPGMA ("average"), the height of a node is the mean dissimilarity between
# the two groups it joins, so the cophenetic distances are group means and the
# usual analysis-of-variance identity gives
#     sum of squared errors = sum(d^2) - sum over nodes of (pairs_k * height_k^2)
# where pairs_k is the number of site pairs first joined at node k. Since
# sum(d^2) is the same for all trees built on the same sites, trees can be
# ranked by sum(pairs_k * height_k^2) alone, which only needs the merge sizes.
#
# For other linkages the cophenetic correlation is computed exactly, by
# listing all pairs of sites (in C++).
#
# `d` holds the dissimilarities in the order of the tree's own sites, unless
# `leaf_site` says which row of `d` each site of the tree stands for, in which
# case `d` can be the whole matrix.
tree_fit_score <- function(tree, d, method, leaf_site = NULL) {
  if (method == "average") {
    sizes <- ihct_node_sizes(tree$merge)
    sum(sizes$pairs * tree$height^2)
  } else {
    ihct_cophenetic_correlation(tree$merge, tree$height, d, leaf_site)
  }
}

# Order trees from best to worst score. Scores that differ by less than a
# small tolerance are treated as equal, in which case earlier runs come first
# (this keeps the choice independent of rounding in the last digits).
rank_by_score <- function(scores, tolerance = 1e-10) {
  o <- order(scores, decreasing = TRUE)
  s <- scores[o]
  same_as_previous <- c(FALSE, abs(diff(s)) <= tolerance * max(abs(s[1]), 1))
  group <- cumsum(!same_as_previous)
  o[order(group, o)]
}

# Height of the division between two groups, according to the linkage method.
# `sites` is the whole group being divided (used by the centroid-based methods).
division_height <- function(dist_mat, group1, group2, sites, method) {
  between <- dist_mat[group1, group2]
  n1 <- length(group1); n2 <- length(group2)
  if (method %in% c("ward.D", "ward.D2", "centroid", "median")) {
    # Note: this is the historical approximation of these linkages in IHCT:
    # the "centroid" of a group is the mean of its rows of dissimilarities to
    # all sites of the group being divided.
    centroid1 <- colMeans(dist_mat[group1, sites, drop = FALSE])
    centroid2 <- colMeans(dist_mat[group2, sites, drop = FALSE])
    centroid_distance <- sum((centroid1 - centroid2)^2)
  }
  h <- switch(method,
              "single" = min(between),
              "complete" = max(between),
              "average" = mean(between),
              "mcquitty" = mean(between),
              "ward.D" = (n1 * n2) / (n1 + n2) * sqrt(centroid_distance),
              "ward.D2" = (n1 * n2) / (n1 + n2) * centroid_distance,
              "centroid" = centroid_distance,
              "median" = 0.5 * centroid_distance,
              stop("method argument is not valid"))
  max(h, 0)
}

# Remove the inversions of a tree, i.e. the divisions that are lower than a
# division they contain, and return the corrected heights.
#
# The tree is given as a node table: for every node, `height` (0 for a single
# site), `left` and `right` (the child nodes, 0 for a single site) and `pairs`
# (the number of site pairs joined at the division, i.e. the number of sites on
# the left times the number of sites on the right). Children always have larger
# node numbers than their parent, so decreasing node numbers go from the bottom
# of the tree to the top.
#
# Two rules:
#
# "max_child" raises every division to the highest division it contains. Simple,
# but it moves a division far above the dissimilarities it summarizes whenever a
# single small group deep in the tree is high.
#
# "least_squares" instead moves the heights as little as possible: it returns the
# monotone heights closest to the divisions' own heights, each weighted by the
# number of site pairs it summarizes. For UPGMA ("average", and "mcquitty" which
# uses the same height here) a division's height is the mean dissimilarity
# between the two groups, cophenetic distances are then group means of the
# dissimilarities, and the sum of squared differences between dissimilarities and
# cophenetic distances is
#     constant + sum over divisions of pairs_k * (height_k - mean_k)^2,
# so these heights are the ones that fit the dissimilarities best on the tree at
# hand: the cophenetic correlation is never below what "max_child" gives, and
# usually above. For the other linkages the same pooling is applied to the
# heights the linkage defines, which keeps them as close to the linkage as
# monotonicity allows but carries no such guarantee.
#
# How it works (isotonic regression on a tree, Pardalos & Xue 1999): heights are
# handled in blocks of divisions that share a common height, starting with one
# block per division. Nodes are visited from the bottom of the tree to the top;
# at each node, as long as a block just below is higher than the node's own
# block, the two are merged and the merged block takes the weighted mean of their
# heights, which may bring further blocks below into conflict. A node is left
# once no block below it is higher, so the heights are monotone at the end.
make_heights_monotone <- function(height, left, right, pairs,
                                  rule = c("least_squares", "max_child")) {
  rule <- match.arg(rule)
  is_group <- left > 0
  bottom_up <- rev(which(is_group))

  if (rule == "max_child") {
    for (id in bottom_up) {
      height[id] <- max(height[id], height[left[id]], height[right[id]])
    }
    return(height)
  }

  # weight and total height of each block, blocks named after their highest node
  weight <- pairs
  total <- pairs * height
  below <- vector("list", length(height))   # blocks directly below each block
  absorbed_by <- integer(length(height))    # 0 while a block still exists

  for (id in bottom_up) {
    # single sites have no height and so never conflict with a division
    conflicting <- c(left[id], right[id])
    conflicting <- conflicting[is_group[conflicting]]
    while (length(conflicting) > 0) {
      heights_below <- total[conflicting] / weight[conflicting]
      highest <- which.max(heights_below)
      if (heights_below[highest] <= total[id] / weight[id]) break
      block <- conflicting[highest]
      weight[id] <- weight[id] + weight[block]
      total[id] <- total[id] + total[block]
      absorbed_by[block] <- id
      conflicting <- c(conflicting[-highest], below[[block]])
      below[block] <- list(NULL)
    }
    below[[id]] <- conflicting
  }

  # a block absorbed into another always has a larger node number than the
  # block that absorbed it, so going top-down resolves the chains in one pass
  final <- integer(length(height))
  for (id in which(is_group)) {
    final[id] <- if (absorbed_by[id] == 0L) id else final[absorbed_by[id]]
  }
  height[is_group] <- (total / weight)[final[is_group]]
  height
}

# Turn the node table into an hclust object. Internal nodes are numbered in
# order of increasing height (ties keep the order in which nodes are met when
# walking the tree, left branch first, children before parents); sites are
# numbered in alphabetical order of their names, as hclust does.
nodes_to_hclust <- function(height, left, right, site, site_names) {
  n_nodes <- length(height)
  is_group <- left > 0
  n_groups <- sum(is_group)

  # nodes met children-first, left branch first (iterative walk from the root)
  visit_order <- integer(n_groups); n_visited <- 0L
  stack_id <- 1L; stack_state <- 0L
  while (length(stack_id) > 0) {
    k <- length(stack_id); id <- stack_id[k]
    if (stack_state[k] == 0L) {
      stack_state[k] <- 1L
      child <- left[id]
      if (is_group[child]) { stack_id <- c(stack_id, child); stack_state <- c(stack_state, 0L) }
    } else if (stack_state[k] == 1L) {
      stack_state[k] <- 2L
      child <- right[id]
      if (is_group[child]) { stack_id <- c(stack_id, child); stack_state <- c(stack_state, 0L) }
    } else {
      n_visited <- n_visited + 1L; visit_order[n_visited] <- id
      stack_id <- stack_id[-k]; stack_state <- stack_state[-k]
    }
  }
  visit_order <- visit_order[order(height[visit_order])]   # stable: ties keep visit order
  new_id <- integer(n_nodes); new_id[visit_order] <- seq_len(n_groups)

  labels <- sort(site_names)
  site_rank <- match(site_names, labels)
  ref <- function(child) if (is_group[child]) new_id[child] else -site_rank[site[child]]
  merge <- cbind(vapply(visit_order, function(id) ref(left[id]), 1L),
                 vapply(visit_order, function(id) ref(right[id]), 1L))

  # left-to-right order of the sites (iterative walk, left branch first)
  leaf_order <- integer(n_nodes - n_groups); n_leaves <- 0L
  stack <- 1L
  while (length(stack) > 0) {
    id <- stack[length(stack)]; stack <- stack[-length(stack)]
    if (is_group[id]) {
      stack <- c(stack, right[id], left[id])
    } else {
      n_leaves <- n_leaves + 1L; leaf_order[n_leaves] <- site_rank[site[id]]
    }
  }

  structure(list(merge = merge, height = height[visit_order], order = leaf_order,
                 labels = labels, method = "Iterative Hierarchical Consensus Tree"),
            class = "hclust")
}
