# Tests of the Iterative Hierarchical Consensus Tree (R/IHCT.R)

# small dissimilarity matrix with many ties, as in real presence/absence data
make_matrix <- function(n_sites = 30, n_species = 40, seed = 1) {
  set.seed(seed)
  comat <- matrix(rbinom(n_sites * n_species, 1, 0.15), n_sites, n_species)
  comat[rowSums(comat) == 0, 1] <- 1
  rownames(comat) <- paste0("site", seq_len(n_sites))
  colnames(comat) <- paste0("sp", seq_len(n_species))
  dis <- dissimilarity(comat, metric = "Simpson")
  as.matrix(stats::as.dist(net_to_mat(dis[, 1:3], weight = TRUE,
                                      squared = TRUE, symmetrical = TRUE)))
}

# structural checks an hclust object must satisfy
expect_valid_tree <- function(hc, d) {
  n <- nrow(d)
  expect_s3_class(hc, "hclust")
  expect_equal(dim(hc$merge), c(n - 1, 2))
  expect_length(hc$height, n - 1)
  expect_setequal(hc$order, seq_len(n))
  expect_setequal(hc$labels, rownames(d))
  expect_setequal(-hc$merge[hc$merge < 0], seq_len(n))          # every site once
  expect_setequal(hc$merge[hc$merge > 0], seq_len(n - 2))       # every node once but the root
  expect_true(all(vapply(seq_len(n - 1), function(k) all(which(hc$merge == k, arr.ind = TRUE)[, 1] > k), TRUE)))
  expect_true(all(is.finite(hc$height)) && all(hc$height >= 0))
  expect_true(all(diff(hc$height) >= -1e-12))                   # no overlapping branches
  # sites of every node form a contiguous block of the leaf order
  pos <- integer(n); pos[hc$order] <- seq_len(n)
  members <- vector("list", n - 1)
  for (k in seq_len(n - 1)) {
    a <- if (hc$merge[k, 1] < 0) -hc$merge[k, 1] else members[[hc$merge[k, 1]]]
    b <- if (hc$merge[k, 2] < 0) -hc$merge[k, 2] else members[[hc$merge[k, 2]]]
    members[[k]] <- c(a, b)
    expect_equal(diff(range(pos[members[[k]]])) + 1, length(members[[k]]))
  }
  expect_silent(stats::cutree(hc, k = 1:min(n, 10)))
  expect_silent(stats::cophenetic(hc))
  expect_silent(stats::as.dendrogram(hc))
}

test_that("closed-form score ranks UPGMA trees like the cophenetic correlation", {
  d <- make_matrix()
  set.seed(2)
  for (i in 1:20) {
    perm <- sample(nrow(d)); dp <- d[perm, perm]
    hc <- fastcluster::hclust(stats::as.dist(dp), "average")
    lt <- dp[lower.tri(dp)]
    sizes <- ihct_node_sizes(hc$merge)
    sse <- sum(lt^2) - sum(sizes$pairs * hc$height^2)
    sst <- sum(lt^2) - sum(lt)^2 / length(lt)
    expect_equal(sqrt(1 - sse / sst), tree_eval(hc, dp)$cophcor, tolerance = 1e-12)
    expect_equal(sizes$size[nrow(hc$merge)], nrow(d))
  }
})

test_that("pair-enumeration cophenetic correlation equals tree_eval for any linkage", {
  d <- make_matrix()
  for (method in c("complete", "single", "mcquitty", "ward.D2")) {
    hc <- fastcluster::hclust(stats::as.dist(d), method)
    expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, d),
                 tree_eval(hc, d)$cophcor, tolerance = 1e-12)
  }
})

test_that("ihct_shuffled_dist gives what as.dist of the shuffled sub-matrix gives", {
  d <- make_matrix(40)
  n <- nrow(d)
  set.seed(3)
  for (i in 1:20) {
    # the whole matrix shuffled, and a group of it shuffled
    sites <- if (i %% 2) sample.int(n) else sample.int(n, sample(2:n, 1))
    expected <- stats::as.dist(d[sites, sites])
    got <- ihct_shuffled_dist(d, sites)
    expect_s3_class(got, "dist")
    expect_equal(as.numeric(got), as.numeric(expected), tolerance = 0)
    expect_equal(attr(got, "Size"), length(sites))
    expect_false(attr(got, "Diag"))
    expect_false(attr(got, "Upper"))
    # what it is there for: the same tree, without the square copy
    expect_equal(fastcluster::hclust(got, "average")$merge,
                 fastcluster::hclust(expected, "average")$merge)
  }
  # the smallest group a tree is built on
  expect_equal(as.numeric(ihct_shuffled_dist(d, c(7L, 2L))), d[7, 2], tolerance = 0)
})

test_that("cophenetic correlation reads the whole matrix when given leaf_site", {
  d <- make_matrix()
  n <- nrow(d)
  set.seed(4)
  for (method in c("average", "complete", "single", "mcquitty")) {
    for (i in 1:5) {
      sites <- sample.int(n, sample(5:n, 1))
      sub <- d[sites, sites]
      hc <- fastcluster::hclust(stats::as.dist(sub), method)
      # site i of the tree is row i of `sub`, or row sites[i] of `d`
      expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, d, sites),
                   ihct_cophenetic_correlation(hc$merge, hc$height, sub),
                   tolerance = 1e-12)
      expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, d, sites),
                   tree_eval(hc, sub)$cophcor, tolerance = 1e-12)
    }
  }
  # leaf_site left out is the same as the identity
  hc <- fastcluster::hclust(stats::as.dist(d), "complete")
  expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, d, seq_len(n)),
               ihct_cophenetic_correlation(hc$merge, hc$height, d))
})

test_that("tree_fit_score scores a tree the same way from either matrix", {
  d <- make_matrix()
  n <- nrow(d)
  set.seed(5)
  for (method in c("average", "complete", "ward.D2")) {
    sites <- sample.int(n)
    hc <- fastcluster::hclust(ihct_shuffled_dist(d, sites), method)
    expect_equal(tree_fit_score(hc, d, method, sites),
                 tree_fit_score(hc, d[sites, sites], method), tolerance = 1e-12)
  }
})

test_that("rank_by_score treats nearly equal scores as ties, earlier runs first", {
  expect_equal(rank_by_score(c(1, 3, 3 + 1e-14, 2)), c(2, 3, 4, 1))
  expect_equal(rank_by_score(c(0, 0, 0)), 1:3)
  expect_equal(rank_by_score(c(5, 4, 6)), c(3, 1, 2))
})

# a random binary tree in the node table format make_heights_monotone() takes:
# heights drawn at random, so parents are routinely lower than their children
random_node_table <- function(n) {
  height <- 0; left <- 0L; right <- 0L; site <- 0L; pairs <- 0
  sites_of <- list(seq_len(n)); n_nodes <- 1L; work <- 1L
  while (length(work)) {
    id <- work[1]; work <- work[-1]
    s <- sites_of[[id]]
    if (length(s) == 1) { site[id] <- s; next }
    k <- sample(seq_len(length(s) - 1), 1)
    for (g in list(s[seq_len(k)], s[-seq_len(k)])) {
      n_nodes <- n_nodes + 1L
      sites_of[[n_nodes]] <- g
      height[n_nodes] <- 0; pairs[n_nodes] <- 0
      left[n_nodes] <- 0L; right[n_nodes] <- 0L; site[n_nodes] <- 0L
      if (left[id] == 0L) left[id] <- n_nodes else right[id] <- n_nodes
    }
    height[id] <- runif(1)
    pairs[id] <- k * (length(s) - k)
    work <- c(work, left[id], right[id])
  }
  list(height = height, left = as.integer(left), right = as.integer(right),
       site = as.integer(site), pairs = pairs)
}

# weighted distance to the divisions' own heights, the quantity least_squares
# minimizes (and, for UPGMA, the fit of the tree to the dissimilarities up to a
# constant)
height_sse <- function(h, tr) sum(tr$pairs * (h - tr$height)^2)

expect_monotone <- function(h, tr) {
  groups <- which(tr$left > 0)
  expect_true(all(h[groups] >= h[tr$left[groups]] - 1e-12))
  expect_true(all(h[groups] >= h[tr$right[groups]] - 1e-12))
  expect_true(all(h >= 0))
}

test_that("least-squares heights are the known optimum on hand-made trees", {
  # a 3-site tree: the root separates one site from two (2 pairs) at 0.2, the
  # division inside the pair (1 pair) is at 0.5. Minimizing
  # 2 (h_root - 0.2)^2 + (h_pair - 0.5)^2 under h_root >= h_pair pools both at
  # (2 x 0.2 + 0.5) / 3 = 0.3
  tr <- list(height = c(0.2, 0, 0.5, 0, 0), left = c(2L, 0L, 4L, 0L, 0L),
             right = c(3L, 0L, 5L, 0L, 0L), pairs = c(2, 0, 1, 0, 0))
  expect_equal(make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                                     "least_squares"),
               c(0.3, 0, 0.3, 0, 0))
  expect_equal(make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                                     "max_child"),
               c(0.5, 0, 0.5, 0, 0))

  # heights already monotone: nothing moves, whatever the rule
  tr$height <- c(0.6, 0, 0.5, 0, 0)
  for (rule in c("least_squares", "max_child")) {
    expect_equal(make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                                       rule), tr$height)
  }

  # a chain of three divisions, one pair each: correcting the middle one (0.1)
  # against its child (0.5) puts it at 0.3, which is now above its own parent
  # (0.2), so the correction has to cascade and the three pool together at the
  # mean of 0.2, 0.1 and 0.5
  tr <- list(height = c(0.2, 0, 0.1, 0, 0.5, 0, 0),
             left  = c(2L, 0L, 4L, 0L, 6L, 0L, 0L),
             right = c(3L, 0L, 5L, 0L, 7L, 0L, 0L),
             pairs = c(1, 0, 1, 0, 1, 0, 0))
  h <- make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                             "least_squares")
  expect_equal(h[c(1, 3, 5)], rep(0.8 / 3, 3))
  expect_monotone(h, tr)
  # one violation only: the child is pooled with its parent, the root stays put
  tr$height <- c(0.9, 0, 0.1, 0, 0.5, 0, 0)
  h <- make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                             "least_squares")
  expect_equal(h[c(1, 3, 5)], c(0.9, 0.3, 0.3))
  expect_monotone(h, tr)
})

test_that("least-squares heights are monotone and fit better than max_child", {
  set.seed(4)
  for (i in 1:20) {
    tr <- random_node_table(sample(2:40, 1))
    h_ls <- make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                                  "least_squares")
    h_mc <- make_heights_monotone(tr$height, tr$left, tr$right, tr$pairs,
                                  "max_child")
    expect_monotone(h_ls, tr)
    expect_monotone(h_mc, tr)
    expect_lte(height_sse(h_ls, tr), height_sse(h_mc, tr) + 1e-12)
    # nothing is moved above the highest or below the lowest division
    groups <- which(tr$left > 0)
    expect_true(all(h_ls[groups] <= max(tr$height[groups]) + 1e-12))
    expect_true(all(h_ls[groups] >= min(tr$height[groups]) - 1e-12))
  }
})

# the groups of sites of every node, as a set; two trees with the same shape
# have the same one, whatever the heights (which decide how nodes are numbered)
clades <- function(hc) {
  members <- vector("list", nrow(hc$merge))
  for (k in seq_len(nrow(hc$merge))) {
    side <- lapply(hc$merge[k, ], function(j) {
      if (j < 0) hc$labels[-j] else members[[j]]
    })
    members[[k]] <- sort(unlist(side))
  }
  sort(vapply(members, paste, "", collapse = ","))
}

test_that("the height rule changes the heights of an IHCT tree but not its shape", {
  d <- make_matrix()
  set.seed(7); ls <- IHCT(d, n_runs = 20, height_rule = "least_squares", verbose = FALSE)
  set.seed(7); mc <- IHCT(d, n_runs = 20, height_rule = "max_child", verbose = FALSE)
  expect_valid_tree(ls, d)
  expect_valid_tree(mc, d)
  # the rule is applied once the topology is decided, so both trees have the
  # same nodes; only the heights, and hence the order the nodes are numbered
  # in, can differ
  expect_identical(clades(ls), clades(mc))
  # a division is never moved above the highest division it contains, so the
  # least-squares heights are never above the max_child ones, and here some
  # division is strictly lower, i.e. inversions did occur on this matrix
  expect_true(all(ls$height <= mc$height + 1e-12))
  expect_true(any(ls$height < mc$height))
  expect_gte(tree_eval(ls, d)$cophcor, tree_eval(mc, d)$cophcor)
  expect_error(IHCT(d, n_runs = 5, height_rule = "highest", verbose = FALSE))
})

test_that("IHCT returns a valid hclust tree, reproducible with a seed", {
  d <- make_matrix()
  set.seed(10); hc1 <- IHCT(d, method = "average", n_runs = 20, top_n_trees = 2, verbose = FALSE)
  set.seed(10); hc2 <- IHCT(d, method = "average", n_runs = 20, top_n_trees = 2, verbose = FALSE)
  expect_valid_tree(hc1, d)
  expect_identical(hc1, hc2)
  # the seed matters. Checked with every source of randomness on: the trees
  # rebuilt at every division (variation_drop = 0) and the tied blocks divided
  # from them rather than resolved directly (tie_block_resolution = FALSE).
  # Both shortcuts draw far fewer random numbers, and on a matrix this small
  # two seeds can then end on the same tree.
  set.seed(10); seed10 <- IHCT(d, n_runs = 20, variation_drop = 0,
                               tie_block_resolution = FALSE, verbose = FALSE)
  set.seed(11); seed11 <- IHCT(d, n_runs = 20, variation_drop = 0,
                               tie_block_resolution = FALSE, verbose = FALSE)
  expect_false(identical(seed10$merge, seed11$merge) &&
                 identical(seed10$height, seed11$height))
  # the tree fits the dissimilarities at least as well as a single UPGMA tree
  single <- fastcluster::hclust(stats::as.dist(d), "average")
  expect_gte(tree_eval(hc1, d)$cophcor, tree_eval(single, d)$cophcor - 1e-8)
})

test_that("IHCT works with the division of the best tree and with other linkages", {
  d <- make_matrix(seed = 3)
  set.seed(1); hc <- IHCT(d, method = "average", n_runs = 10, top_n_trees = 1, verbose = FALSE)
  expect_valid_tree(hc, d)
  for (method in c("complete", "single", "ward.D2")) {
    set.seed(1); hc <- IHCT(d, method = method, n_runs = 5, top_n_trees = 2, verbose = FALSE)
    expect_valid_tree(hc, d)
  }
})

test_that("IHCT handles tiny and degenerate matrices", {
  d2 <- matrix(c(0, 0.3, 0.3, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  hc <- IHCT(d2, n_runs = 5, verbose = FALSE)
  expect_valid_tree(hc, d2); expect_equal(hc$height, 0.3)
  d3 <- matrix(c(0, 0.2, 0.9, 0.2, 0, 0.8, 0.9, 0.8, 0), 3, dimnames = list(letters[1:3], letters[1:3]))
  set.seed(1); hc <- IHCT(d3, n_runs = 5, verbose = FALSE)
  expect_valid_tree(hc, d3)
  expect_equal(stats::cutree(hc, 2)[c("a", "b")], c(a = 1, b = 1))
  # identical sites (all dissimilarities 0) and saturated sites (all 1)
  d0 <- matrix(0, 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  set.seed(1); expect_valid_tree(IHCT(d0, n_runs = 5, verbose = FALSE), d0)
  d1 <- matrix(1, 5, 5, dimnames = list(letters[1:5], letters[1:5])); diag(d1) <- 0
  set.seed(1); expect_valid_tree(IHCT(d1, n_runs = 5, verbose = FALSE), d1)
})

test_that("a deep chain-shaped tree does not exhaust the call stack", {
  n <- 600
  d <- outer(seq_len(n), seq_len(n), pmax) / n; diag(d) <- 0
  dimnames(d) <- list(paste0("s", seq_len(n)), paste0("s", seq_len(n)))
  set.seed(1)
  hc <- IHCT(d, n_runs = 2, verbose = FALSE)
  expect_valid_tree(hc, d)
  expect_equal(tree_eval(hc, d)$cophcor, 1, tolerance = 1e-8)
})

# --- inherited trees (pruning) and the variation cutoff --------------------------

test_that("pruning a tree gives the induced tree with mean-dissimilarity heights", {
  d <- make_matrix(n_sites = 40, seed = 2)
  set.seed(5)
  for (i in 1:30) {
    sites <- sample(nrow(d), sample(8:nrow(d), 1))
    tree <- fresh_trees(d, sites, "average", 1)[[1]]
    keep <- logical(nrow(d))
    keep[sample(sites, sample(3:length(sites), 1))] <- TRUE
    pruned <- ihct_prune_tree(tree$merge, tree$height, tree$pairs,
                              tree$leaf_site, keep, d)
    expect_setequal(pruned$leaf_site, sites[keep[sites]])
    # the sites, height and number of pairs of every node, against a pure-R
    # pruning of the same tree (helper-ihct.R)
    expect_identical(ref_tree_key(pruned$merge, pruned$height, pruned$pairs,
                                  pruned$leaf_site),
                     ref_prune_tree(tree$merge, tree$leaf_site, keep, d))
    # pruning again is pruning once to the smaller set
    keep2 <- keep
    keep2[sample(pruned$leaf_site, 2)] <- FALSE
    twice <- ihct_prune_tree(pruned$merge, pruned$height, pruned$pairs,
                             pruned$leaf_site, keep2, d)
    expect_identical(ref_tree_key(twice$merge, twice$height, twice$pairs,
                                  twice$leaf_site),
                     ref_prune_tree(tree$merge, tree$leaf_site, keep2, d))
  }
})

test_that("pruning to a clade gives the UPGMA tree of the sub-matrix", {
  # the sites of a clade of a UPGMA tree only ever join each other, so the
  # sub-tree on them is the UPGMA tree of their sub-matrix. Ties would let the
  # two runs break them differently, so the matrix is jittered to remove them.
  d <- make_matrix(n_sites = 50, seed = 4)
  set.seed(6)
  noise <- matrix(runif(nrow(d)^2, 0, 1e-6), nrow(d))
  noise[lower.tri(noise)] <- t(noise)[lower.tri(noise)]
  d <- d + noise
  diag(d) <- 0

  tree <- fresh_trees(d, seq_len(nrow(d)), "average", 1)[[1]]
  members <- ref_node_members(tree$merge, tree$leaf_site)
  clade_nodes <- which(lengths(members) >= 4 & lengths(members) <= nrow(d) - 2)
  expect_gt(length(clade_nodes), 5)
  for (k in clade_nodes) {
    keep <- logical(nrow(d))
    keep[members[[k]]] <- TRUE
    pruned <- ihct_prune_tree(tree$merge, tree$height, tree$pairs,
                              tree$leaf_site, keep, d)
    sites <- members[[k]]
    upgma <- fastcluster::hclust(stats::as.dist(d[sites, sites]), "average")
    expect_identical(ref_tree_key(pruned$merge, pruned$height, pruned$pairs,
                                  pruned$leaf_site),
                     ref_tree_key(upgma$merge, upgma$height,
                                  ihct_node_sizes(upgma$merge)$pairs, sites))
  }
})

test_that("ihct_top_division gives the two branches of the root, as cutree does", {
  d <- make_matrix(n_sites = 30, seed = 5)
  set.seed(8)
  for (i in 1:20) {
    sites <- sample(nrow(d), sample(4:nrow(d), 1))
    tree <- fresh_trees(d, sites, "average", 1)[[1]]
    division <- ihct_top_division(tree$merge, tree$leaf_site)
    expect_setequal(unlist(division), sites)
    expect_equal(lapply(division, sort),
                 lapply(ref_top_division(tree$merge, tree$leaf_site), sort))
    # cutree(k = 2) cuts a tree above its last merge, i.e. at its root
    hc <- structure(list(merge = tree$merge, height = sort(tree$height),
                         order = seq_along(tree$leaf_site),
                         labels = as.character(tree$leaf_site),
                         method = "average"), class = "hclust")
    groups <- stats::cutree(hc, k = 2)
    expect_setequal(as.integer(names(groups)[groups == groups[1]]),
                    if (as.integer(names(groups)[1]) %in% division[[1]])
                      division[[1]] else division[[2]])
  }
})

test_that("a group whose dissimilarities are all equal is recognized", {
  d <- make_matrix(n_sites = 20, seed = 6)
  flat <- matrix(0.4, 5, 5); diag(flat) <- 0
  expect_true(is_tied_block(flat, 1:5))
  expect_true(is_tied_block(flat, c(2, 4, 5)))
  flat[1, 3] <- flat[3, 1] <- 0.6
  expect_false(is_tied_block(flat, 1:5))
  expect_true(is_tied_block(flat, c(2, 4, 5)))
  expect_false(is_tied_block(d, seq_len(nrow(d))))
})

test_that("groups with all dissimilarities equal become a chain at that height", {
  # sites that share no species at all are all at dissimilarity 1 from each
  # other: every tree on them is as good as any other, and IHCT peels them off
  # one by one instead of randomizing
  d <- matrix(1, 8, 8, dimnames = list(letters[1:8], letters[1:8]))
  diag(d) <- 0
  set.seed(1)
  hc <- IHCT(d, n_runs = 20, variation_drop = 0.2, verbose = FALSE)
  expect_valid_tree(hc, d)
  expect_equal(hc$height, rep(1, 7))
  # a chain: every merge but the first joins a single site to what came before
  expect_equal(sum(hc$merge[, 1] < 0 & hc$merge[, 2] < 0), 1)
  # the same with identical sites, which are all at dissimilarity 0
  d0 <- matrix(0, 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  set.seed(1)
  hc <- IHCT(d0, n_runs = 20, variation_drop = 0.2, verbose = FALSE)
  expect_valid_tree(hc, d0)
  expect_equal(hc$height, rep(0, 5))
})

test_that("inherited trees give valid trees for every division rule", {
  d <- make_matrix(n_sites = 50, seed = 7)
  for (cutoff in c(0, 0.1, 0.5, 1)) {
    for (top_n_trees in c(1, 2, 3)) {
      set.seed(3)
      hc <- IHCT(d, n_runs = 20, top_n_trees = top_n_trees, variation_drop = cutoff,
                 verbose = FALSE)
      expect_valid_tree(hc, d)
      # inheriting trees must not undo what the randomization is for: the tree
      # still has to fit the dissimilarities at least as well as a single one
      expect_gte(tree_eval(hc, d)$cophcor,
                 tree_eval(fastcluster::hclust(stats::as.dist(d), "average"),
                           d)$cophcor - 1e-8)
    }
  }
})

test_that("the variation cutoff changes the trees, saves runs, and is reproducible", {
  d <- make_matrix(n_sites = 60, seed = 8)
  set.seed(3); a <- IHCT(d, n_runs = 20, variation_drop = 0.3, verbose = FALSE)
  set.seed(3); b <- IHCT(d, n_runs = 20, variation_drop = 0.3, verbose = FALSE)
  set.seed(3); none <- IHCT(d, n_runs = 20, variation_drop = 0, verbose = FALSE)
  expect_identical(a, b)
  expect_false(identical(clades(a), clades(none)))

  # count the divisions that built their own trees instead of taking their
  # parent's, by watching the only function that builds them
  count_fresh <- function(cutoff) {
    n <- 0L
    build <- fresh_trees
    local_mocked_bindings(
      fresh_trees = function(...) { n <<- n + 1L; build(...) },
      .package = "bioregion")
    set.seed(3)
    IHCT(d, n_runs = 20, variation_drop = cutoff, verbose = FALSE)
    n
  }
  expect_lt(count_fresh(0.3), count_fresh(0))
})

test_that("only average linkage inherits trees, and the cutoff is checked", {
  d <- make_matrix(n_sites = 30, seed = 9)
  # the heights a pruned tree gets are mean dissimilarities, which is what
  # UPGMA defines; the other methods keep building their trees afresh, so
  # the cutoff makes no difference to them
  for (method in c("complete", "ward.D2")) {
    set.seed(2)
    with_cutoff <- suppressMessages(
      IHCT(d, method = method, n_runs = 10, variation_drop = 1, verbose = FALSE))
    set.seed(2)
    without <- IHCT(d, method = method, n_runs = 10, variation_drop = 0, verbose = FALSE)
    expect_identical(with_cutoff, without)
  }
  expect_message(IHCT(d, method = "complete", n_runs = 5, variation_drop = 0.5,
                      verbose = TRUE),
                 "only available with method")
  expect_error(IHCT(d, n_runs = 5, variation_drop = 2, verbose = FALSE),
               "between 0 and 1")
  expect_error(IHCT(d, n_runs = 5, variation_drop = -0.1, verbose = FALSE),
               "between 0 and 1")
  expect_error(IHCT(d, n_runs = 5, variation_drop = NA, verbose = FALSE),
               "between 0 and 1")
})

test_that("sites_drop builds trees again after a number of sites is lost", {
  # A matrix whose tree peels a single site off at a time: site n is further
  # from everything than site n - 1, and so on, so every division separates the
  # highest site from the rest. This is the shape sites_drop is for -- the group
  # barely changes from one division to the next, so its share of the variation
  # falls too slowly for variation_drop to rebuild anything near the top.
  n <- 60
  d <- outer(seq_len(n), seq_len(n), function(i, j) pmax(i, j) / n)
  diag(d) <- 0
  dimnames(d) <- list(sprintf("s%02d", seq_len(n)), sprintf("s%02d", seq_len(n)))
  count_fresh <- function(...) {
    k <- 0L
    build <- fresh_trees
    local_mocked_bindings(
      fresh_trees = function(...) { k <<- k + 1L; build(...) },
      .package = "bioregion")
    set.seed(3)
    IHCT(d, n_runs = 20, verbose = FALSE, ...)
    k
  }
  # the second criterion can only add rebuilds, never remove any
  variation_only <- count_fresh(variation_drop = 0.2, sites_drop = Inf)
  expect_gt(count_fresh(variation_drop = 0.2, sites_drop = 5), variation_only)
  # and the smaller the count, the more often trees are built again
  expect_gt(count_fresh(variation_drop = 0.2, sites_drop = 2),
            count_fresh(variation_drop = 0.2, sites_drop = 10))
  # a plain count of sites: a division always removes at least one site from a
  # group, so 1 or less asks for fresh trees everywhere, exactly as
  # variation_drop = 0 does
  expect_equal(count_fresh(variation_drop = 0.2, sites_drop = 1),
               count_fresh(variation_drop = 0, sites_drop = 10))
  expect_equal(count_fresh(variation_drop = 0.2, sites_drop = 0),
               count_fresh(variation_drop = 0, sites_drop = 10))
  # never more than building them at every division
  expect_lte(count_fresh(variation_drop = 0.2, sites_drop = 2),
             count_fresh(variation_drop = 0))

  # variation_drop = 0 means "build them everywhere", whatever sites_drop says
  set.seed(3); with_count <- IHCT(d, n_runs = 20, variation_drop = 0,
                                  sites_drop = 5, verbose = FALSE)
  set.seed(3); without <- IHCT(d, n_runs = 20, variation_drop = 0,
                               sites_drop = Inf, verbose = FALSE)
  expect_identical(with_count, without)

  # reproducible, and still a valid tree
  set.seed(4); a <- IHCT(d, n_runs = 20, sites_drop = 5, verbose = FALSE)
  set.seed(4); b <- IHCT(d, n_runs = 20, sites_drop = 5, verbose = FALSE)
  expect_identical(a, b)
  expect_valid_tree(a, d)

  expect_error(IHCT(d, n_runs = 5, sites_drop = -1, verbose = FALSE),
               "sites_drop must be a single number")
  expect_error(IHCT(d, n_runs = 5, sites_drop = NA, verbose = FALSE),
               "sites_drop must be a single number")
  expect_error(IHCT(d, n_runs = 5, sites_drop = c(2, 3), verbose = FALSE),
               "sites_drop must be a single number")
})

test_that("sites_drop rarely fires on a tree that splits into balanced halves", {
  # a group there loses many sites at a time, so its share of the variation has
  # already fallen and variation_drop has rebuilt the trees anyway: the count
  # costs nothing on this shape of data, which is why it can be on by default
  d <- make_matrix(n_sites = 60, seed = 8)
  count_fresh <- function(...) {
    k <- 0L
    build <- fresh_trees
    local_mocked_bindings(
      fresh_trees = function(...) { k <<- k + 1L; build(...) },
      .package = "bioregion")
    set.seed(3)
    IHCT(d, n_runs = 20, verbose = FALSE, ...)
    k
  }
  expect_equal(count_fresh(variation_drop = 0.2, sites_drop = 5),
               count_fresh(variation_drop = 0.2, sites_drop = Inf))
})

test_that("the defaults reuse trees and fit the dissimilarities about as well", {
  d <- make_matrix(n_sites = 60, seed = 8)
  expect_equal(formals(IHCT)$variation_drop, 0.2)
  expect_equal(formals(IHCT)$sites_drop, 10)
  set.seed(5); default <- IHCT(d, n_runs = 20, verbose = FALSE)
  set.seed(5); everywhere <- IHCT(d, n_runs = 20, variation_drop = 0, verbose = FALSE)
  expect_valid_tree(default, d)
  # the whole point of the defaults: nearly the same fit for less work
  expect_gt(tree_eval(default, d)$cophcor,
            tree_eval(everywhere, d)$cophcor - 0.01)
})

# Tests for tie_block_resolution ----------------------------------------------
test_that("a tied block is resolved as a rake at its own dissimilarity", {
  # three sites all at 0.4 from each other, joined to a fourth further away
  d <- matrix(0.9, 4, 4); diag(d) <- 0
  d[1:3, 1:3] <- 0.4; diag(d) <- 0
  rownames(d) <- colnames(d) <- paste0("s", 1:4)

  expect_true(is_tied_block(d, 1:3))
  expect_false(is_tied_block(d, 1:4))

  set.seed(1); hc <- IHCT(d, n_runs = 10, verbose = FALSE)
  expect_valid_tree(hc, d)
  # the block's two divisions both sit at its common value
  expect_equal(sort(hc$height), c(0.4, 0.4, 0.9))
})

test_that("tie_block_resolution = FALSE divides tied blocks from trees instead", {
  d <- make_matrix(n_sites = 40, seed = 4)
  set.seed(6); on  <- IHCT(d, n_runs = 20, verbose = FALSE)
  set.seed(6); off <- IHCT(d, n_runs = 20, tie_block_resolution = FALSE,
                           verbose = FALSE)
  expect_valid_tree(on, d)
  expect_valid_tree(off, d)
  # it changes the shape of the tree inside the tied blocks ...
  expect_false(identical(on$merge, off$merge))
  # ... but every tree on a tied block reproduces its dissimilarities exactly,
  # so it is not a way of fitting the data better or worse
  expect_equal(tree_eval(on, d)$cophcor, tree_eval(off, d)$cophcor,
               tolerance = 0.01)

  # it is the switch that has to be off to get the bioregion 1.4.0 algorithm
  set.seed(6); legacy_settings <- IHCT(d, n_runs = 20, variation_drop = 0,
                                       height_rule = "max_child",
                                       tie_block_resolution = FALSE,
                                       verbose = FALSE)
  expect_valid_tree(legacy_settings, d)

  expect_error(IHCT(d, n_runs = 5, tie_block_resolution = NA, verbose = FALSE),
               "tie_block_resolution must be TRUE or FALSE")
  expect_error(IHCT(d, n_runs = 5, tie_block_resolution = "yes", verbose = FALSE),
               "tie_block_resolution must be TRUE or FALSE")
})

test_that("tied blocks are only resolved where the linkage gives them their own height", {
  d <- matrix(0.9, 6, 6); d[1:4, 1:4] <- 0.4; diag(d) <- 0
  rownames(d) <- colnames(d) <- paste0("s", 1:6)
  # min, max and mean of a tied block all return its common value
  for (meth in c("single", "complete", "average", "mcquitty")) {
    set.seed(2); hc <- IHCT(d, method = meth, n_runs = 10, verbose = FALSE)
    expect_true(sum(abs(hc$height - 0.4) < 1e-12) >= 3,
                label = paste0(meth, ": divisions of the tied block at 0.4"))
  }
  # the centroid-based linkages define their heights otherwise, so the block is
  # divided from randomised trees like any other group and the user is told
  expect_message(IHCT(d, method = "ward.D2", n_runs = 5, verbose = TRUE),
                 "only available with method")
  set.seed(2)
  expect_silent(IHCT(d, method = "ward.D2", n_runs = 5,
                     tie_block_resolution = FALSE, verbose = FALSE))
})

test_that("either criterion at its 'randomize always' value switches reuse off", {
  d <- make_matrix(n_sites = 40, seed = 9)
  # no trees are ever passed down, so the two shortcuts cannot differ
  set.seed(7); by_variation <- IHCT(d, n_runs = 20, variation_drop = 0,
                                    sites_drop = 10, verbose = FALSE)
  set.seed(7); by_sites <- IHCT(d, n_runs = 20, variation_drop = 0.2,
                                sites_drop = 1, verbose = FALSE)
  set.seed(7); by_sites0 <- IHCT(d, n_runs = 20, variation_drop = 0.2,
                                 sites_drop = 0, verbose = FALSE)
  expect_identical(by_variation, by_sites)
  expect_identical(by_variation, by_sites0)

  # and at the other end, both switched off builds the trees once, at the root
  k <- 0L
  build <- fresh_trees
  local_mocked_bindings(
    fresh_trees = function(...) { k <<- k + 1L; build(...) },
    .package = "bioregion")
  set.seed(7); once <- IHCT(d, n_runs = 20, variation_drop = 1,
                            sites_drop = Inf, verbose = FALSE)
  expect_equal(k, 1L)
  expect_valid_tree(once, d)
})
