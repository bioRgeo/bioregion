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
#   1. n_runs candidate trees are obtained (see below).
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
# Where the candidate trees come from. Building n_runs trees at every group is
# what makes IHCT slow: with species composition data most divisions peel one
# site off the rest, so almost every group is nearly as large as the whole
# dataset and the cost grows with the cube of the number of sites. A group can
# take its parent's trees instead, pruned to its own sites: the sites of the
# other half are dropped, the nodes left with a single child are contracted,
# and every node that lost sites gets the mean dissimilarity between the two
# groups it still joins (see ihct_prune_tree). Pruned trees keep n_runs
# different tie-breakings alive for a fraction of the cost, but they are
# reshaped a little more at every level, so they must not be used where the
# division still matters much.
#
# Where it still matters is measured by the share of total variation of a group
# (see variation_share): two sites are separated at exactly one division, so
# the pairs inside a group are the business of that group's own sub-tree and of
# no other part of the tree, and the share of the dataset's total sum of
# squares those pairs carry bounds what the divisions inside the group can
# still win or lose. It is 1 at the whole set of sites and falls as the tree is
# built.
#
# Trees are rebuilt when a group has lost `variation_drop` of the variation the
# last group they were built for still carried, and are inherited in between.
# The comparison is against the last rebuild and not against a fixed level, so
# that rebuilds keep happening all the way down: a fixed level would be crossed
# once on any path -- the share only ever falls -- and everything below it would
# inherit for the rest of its depth however often its trees were pruned. Each
# rebuild moves the mark down with it, so the rebuilds on a path are spaced
# geometrically and there are at most log(smallest share) / log(1 - drop) of
# them.
#
# The share alone is not enough on a tree that peels a single site off at a
# time, which is what species composition data usually gives: the group barely
# changes from one division to the next, so its share falls very slowly and the
# trees are inherited through dozens of divisions near the root, which are
# exactly the divisions that own most of the pairs. `sites_drop` rebuilds them
# on a second criterion -- the number of sites lost since the last rebuild --
# which on such a tree fires every `sites_drop` divisions wherever it is, and
# on a balanced tree rarely fires at all, since a division there drops many
# sites at once and the share has already fallen. The two together were
# measured to give up less cophenetic correlation than either alone at the same
# cost.
#
# Tied blocks. A group whose dissimilarities are all the same value is a tied
# block: every tree on it has the same cophenetic distances, so all n_runs
# candidates score alike and the randomization has nothing to choose between.
# `tie_block_resolution` (TRUE by default) resolves such a group directly, as a
# rake whose divisions all sit at that common value (see resolve_tied_block),
# which costs no cophenetic correlation and saves the trees. It is independent
# of the reuse machinery -- a tied block needs no trees whether or not trees
# are being passed down -- but it is not valid for every linkage: the height it
# writes is the block's common value, which is what min, max and mean all
# return, so it applies to "single", "complete", "average" and "mcquitty" and
# not to the centroid-based linkages, whose heights are not dissimilarities.
# Note that the rake is one arbitrary resolution among many equally good
# ones; the randomized trees of bioregion 1.4.0 picked another (a balanced
# split), so the two differ in topology inside a tied block while fitting the
# dissimilarities identically.
#
# Switching the criteria off. Both are dials that run the same way: the lower
# the value the more often trees are rebuilt, and the top of each range turns
# that criterion off. `variation_drop = 1` leaves the rebuilds to `sites_drop`
# alone, `sites_drop = Inf` leaves them to `variation_drop` alone, and both at
# once builds the trees once at the root and prunes them all the way down. At
# the other end, either criterion at its "rebuild always" value (`variation_drop
# = 0`, or `sites_drop` of 0 or 1, since every division loses at least one site)
# forces fresh trees at every group and so makes the other criterion idle; that
# is why `pass_trees_down` is false for either of them. With the tied blocks
# left alone as well, the algorithm is the one of bioregion 1.4.0 and earlier:
# `variation_drop = 0, top_n_trees = 2, height_rule = "max_child",
# tie_block_resolution = FALSE` reproduces the trees of those versions seed for
# seed.
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
# - A candidate tree is a list(merge, height, pairs, leaf_site, score): the
#   merge matrix and the heights of an hclust object, the number of site pairs
#   joined at each node, the site of dist_mat every leaf stands for, and the
#   score of step 2. Fresh and inherited trees have the same form, so the rest
#   of the code does not have to know where a tree comes from.
# - For UPGMA, the cophenetic correlation of a tree is obtained from the merge
#   sizes and heights only (see tree_fit_score), which avoids computing the
#   full cophenetic matrix of every randomized tree.

IHCT <- function(dist_mat,
                 method = "average",
                 n_runs = 100,
                 top_n_trees = 2,
                 variation_drop = 0.2,
                 sites_drop = 10,
                 height_rule = c("least_squares", "max_child"),
                 tie_block_resolution = TRUE,
                 verbose = TRUE) {

  # checked here rather than where it is used, at the end of the run
  height_rule <- match.arg(height_rule)
  n <- nrow(dist_mat)
  site_names <- rownames(dist_mat)
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
                               top_n_trees, if (build_fresh) NULL else trees)
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

# n_runs trees for a group of sites (integer positions), each built by
# fastcluster on the dissimilarities of the group with its sites shuffled at
# random, and scored by its fit to the dissimilarities.
#
# The shuffled dissimilarities are read straight out of `dist_mat` into the
# vector fastcluster expects (ihct_shuffled_dist), and the score reads out of
# `dist_mat` too, so the square sub-matrix of the group is never built. The
# random draw, and therefore the trees, are the same as when it was.
fresh_trees <- function(dist_mat, sites, method, n_runs) {
  m <- length(sites)
  trees <- vector("list", n_runs)
  for (run in seq_len(n_runs)) {
    shuffled <- sites[sample.int(m)]
    tree <- fastcluster::hclust(ihct_shuffled_dist(dist_mat, shuffled),
                                method = method)
    trees[[run]] <- list(merge = tree$merge, height = tree$height,
                         pairs = ihct_node_sizes(tree$merge)$pairs,
                         leaf_site = shuffled,
                         score = tree_fit_score(tree, dist_mat, method, shuffled))
  }
  trees
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
# trees: the ones given in `trees`, or fresh ones when it is NULL. Returns the
# two groups, each sorted by site name, and the trees the division was read
# from, which the two groups may inherit.
divide_sites <- function(dist_mat, sites, site_names, method, n_runs,
                         top_n_trees, trees = NULL) {
  m <- length(sites)
  if (is.null(trees)) trees <- fresh_trees(dist_mat, sites, method, n_runs)

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
