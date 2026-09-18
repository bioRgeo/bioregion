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
