# Iterative Hierarchical Consensus Tree (IHCT)
#
# Builds a hierarchical tree of sites from a dissimilarity matrix, from the top
# down: the whole set of sites is divided in two groups, then each group is
# divided in two, and so on until every group is a single site.
#
# Why not simply use hclust? Dissimilarity matrices of species composition
# contain many tied values, and hierarchical clustering breaks ties according
# to the order of sites in the matrix, so the tree depends on that arbitrary
# order (Dapporto et al. 2013). IHCT therefore decides every division from
# many trees built on randomized versions of the matrix, and re-randomizes
# inside each group so that the decision taken at one level does not inherit
# the tie-breaking of the level above.
#
# What happens at one group of sites:
#   1. n_runs trees are built with fastcluster::hclust, each on a copy of the
#      dissimilarity sub-matrix whose sites were shuffled at random.
#   2. Each tree is scored by how well its cophenetic distances reproduce the
#      dissimilarities (cophenetic correlation). The top_n_trees best trees
#      are kept.
#   3. The two groups are the two branches of the best tree (top_n_trees = 1),
#      or the two groups obtained by clustering the sites according to how
#      often they fall on the same side in the best trees (top_n_trees > 1,
#      "co-assignment").
#   4. The height of the division is computed from the dissimilarities between
#      the two groups according to the linkage method (mean for UPGMA).
# The two groups are then processed in turn, the smaller one first. When all
# groups are single sites, heights are made monotone (a group is never lower
# than the groups it contains) and the tree is returned as an hclust object.
#
# Divisions are decided independently at each group, from different randomized
# trees, so a division can come out lower than a division inside one of its
# groups. Two rules remove these inversions (see make_heights_monotone):
# "least_squares" (default) moves the offending heights as little as possible,
# "max_child" raises a division to the highest division it contains (the rule
# used up to bioregion 1.4.0).
#
# Implementation notes (2026-09):
# - Sites are handled as integer positions in dist_mat; names are only used to
#   sort sites (the co-assignment step works on sites sorted by name, as we
#   originally implementated) and to label the final tree.
# - The groups still to be divided are kept in a work list instead of using R
#   recursion, so that deep chain-like trees do not exhaust R's call stack.
# - For UPGMA, the cophenetic correlation of a tree is obtained from the merge
#   sizes and heights only (see tree_fit_score), which avoids computing the
#   full cophenetic matrix of every randomized tree.

IHCT <- function(dist_mat,
                 method = "average",
                 n_runs = 100,
                 top_n_trees = 2,
                 height_rule = c("least_squares", "max_child"),
                 verbose = TRUE) {

  # checked here rather than where it is used, at the end of the run
  height_rule <- match.arg(height_rule)
  n <- nrow(dist_mat)
  site_names <- rownames(dist_mat)
  if (top_n_trees > n_runs) top_n_trees <- n_runs
  if (n == 1) stop("At least two sites are needed to build a tree.")

  # --- The tree under construction -------------------------------------------
  # One row per node (groups of sites and single sites alike), at most 2n - 1.
  n_max <- 2 * n - 1
  node_height <- numeric(n_max)        # height of the division (0 for sites)
  node_left <- integer(n_max)          # child nodes (0 for sites)
  node_right <- integer(n_max)
  node_site <- integer(n_max)          # site number for single-site nodes, 0 otherwise
  node_pairs <- numeric(n_max)         # site pairs joined at the division (0 for sites)
  node_sites <- vector("list", n_max)  # sites of groups still to be divided
  n_nodes <- 1L
  node_sites[[1]] <- seq_len(n)        # the root holds every site, in matrix order
  work <- 1L                           # groups waiting to be divided (last in, first out)
  n_placed <- 0L                       # progress: sites turned into leaves

  new_node <- function(sites) {
    n_nodes <<- n_nodes + 1L
    node_sites[[n_nodes]] <<- sites
    n_nodes
  }

  while (length(work) > 0) {
    id <- work[length(work)]
    work <- work[-length(work)]
    sites <- node_sites[[id]]
    node_sites[id] <- list(NULL)
    m <- length(sites)

    if (m == 1) {
      node_site[id] <- sites
      n_placed <- n_placed + 1L
      next
    }

    if (m == 2) {
      groups <- list(sites[1], sites[2])
    } else {
      groups <- divide_sites(dist_mat, sites, site_names, method, n_runs, top_n_trees)
    }
    # smaller group first (same convention as hclust)
    if (length(groups[[1]]) > length(groups[[2]])) groups <- rev(groups)

    node_height[id] <- division_height(dist_mat, groups[[1]], groups[[2]], sites, method)
    node_pairs[id] <- length(groups[[1]]) * length(groups[[2]])
    node_left[id] <- new_node(groups[[1]])
    node_right[id] <- new_node(groups[[2]])
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

# Divide a group of sites (integer positions) in two, from n_runs randomized
# trees. Returns a list of two integer vectors, each sorted by site name.
divide_sites <- function(dist_mat, sites, site_names, method, n_runs, top_n_trees) {
  m <- length(sites)

  # 1. randomized trees, each scored by its fit to the dissimilarities
  trees <- vector("list", n_runs)
  scores <- numeric(n_runs)
  for (run in seq_len(n_runs)) {
    shuffled <- sites[sample.int(m)]
    d_run <- dist_mat[shuffled, shuffled]
    trees[[run]] <- fastcluster::hclust(stats::as.dist(d_run), method = method)
    scores[run] <- tree_fit_score(trees[[run]], d_run, method)
  }

  # 2. keep the best trees (in case of equal scores, the first runs)
  best <- rank_by_score(scores)[seq_len(top_n_trees)]

  # 3. two groups from the best trees, working on sites sorted by name
  sorted <- sites[order(site_names[sites])]
  side <- vapply(trees[best], function(tree) {
    stats::cutree(tree, k = 2)[site_names[sorted]]
  }, integer(m))
  if (top_n_trees == 1) {
    membership <- side[, 1]
  } else {
    # co-assignment: how often two sites fall on the same side, turned into a
    # dissimilarity, clustered with the same linkage and cut in two
    together <- matrix(0, m, m)
    for (r in seq_len(ncol(side))) together <- together + outer(side[, r], side[, r], "==")
    coassign <- 1 - together / n_runs
    dimnames(coassign) <- list(site_names[sorted], site_names[sorted])
    membership <- stats::cutree(fastcluster::hclust(stats::as.dist(coassign), method = method), k = 2)
  }
  # group 1 is the one containing the first site in name order (cutree convention)
  list(sorted[membership == membership[1]], sorted[membership != membership[1]])
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
tree_fit_score <- function(tree, d, method) {
  if (method == "average") {
    sizes <- ihct_node_sizes(tree$merge)
    sum(sizes$pairs * tree$height^2)
  } else {
    ihct_cophenetic_correlation(tree$merge, tree$height, d)
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
