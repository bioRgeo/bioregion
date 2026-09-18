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
  set.seed(11); hc3 <- IHCT(d, method = "average", n_runs = 20, top_n_trees = 2, verbose = FALSE)
  expect_valid_tree(hc1, d)
  expect_identical(hc1, hc2)
  expect_false(identical(hc1$merge, hc3$merge) && identical(hc1$height, hc3$height))
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
