# Tests of tree_eval() (R/utils.R) and tree_eval_cpp() (src/ihct_trees.cpp)
#
# tree_eval() used to build the full cophenetic matrix and call stats::cor();
# it now reads the fit of the tree straight from its merges in C++. These
# tests check that the two give the same cophenetic correlation and msd, to
# within 1e-15 (a few units in the last place), on every kind of tree the
# package evaluates.

# the previous tree_eval(), verbatim, as the reference
tree_eval_reference <- function(tree, dist_mat, method = "pearson") {
  if(inherits(dist_mat, "dist")){
    dist_mat <- as.matrix(dist_mat)
  }
  coph <- as.matrix(stats::cophenetic(tree))
  coph <- coph[match(rownames(dist_mat), rownames(coph)),
               match(rownames(dist_mat), colnames(coph))]
  lower_tri_idx <- lower.tri(dist_mat)
  cophcor <- stats::cor(dist_mat[lower_tri_idx],
                        coph[lower_tri_idx],
                        method = method)
  diff_matrix <- dist_mat - coph
  msd <- mean(diff_matrix[lower_tri_idx]^2)
  return(list(cophcor = cophcor,
              msd = msd))
}

# cophcor within 1e-15, msd within 1e-15 of its own size
expect_same_fit <- function(new, ref, label = "") {
  expect_lte(abs(new$cophcor - ref$cophcor), 1e-15,
             label = paste(label, "cophcor difference"))
  expect_lte(abs(new$msd - ref$msd), 1e-15 * max(ref$msd, 1e-300),
             label = paste(label, "msd difference"))
}

# Simpson dissimilarities of a sample of sites, as a named square matrix
sample_dissim <- function(comat, n_sites, metric = "Simpson", seed = 1) {
  set.seed(seed)
  comat <- comat[sample(nrow(comat), min(n_sites, nrow(comat))), ]
  dis <- dissimilarity(comat, metric = metric)
  as.matrix(stats::as.dist(net_to_mat(dis[, 1:3], weight = TRUE,
                                      squared = TRUE, symmetrical = TRUE)))
}

LINKAGES <- c("average", "complete", "single", "mcquitty",
              "ward.D", "ward.D2", "centroid", "median")

test_that("tree_eval matches the previous version on hclust trees", {
  mats <- list(fish = sample_dissim(fishmat, 300),
               vege = sample_dissim(vegemat, 300),
               fish_sorensen = sample_dissim(fishmat, 300, "Sorensen"))
  set.seed(2)
  x <- matrix(stats::runif(300 * 10), 300)
  rownames(x) <- paste0("site", seq_len(300))
  mats$euclidean <- as.matrix(stats::dist(x))

  for (nm in names(mats)) {
    D <- mats[[nm]]
    for (method in LINKAGES) {
      # a tree built on a shuffled copy, so that its sites are not in the
      # order of the matrix it is evaluated against
      p <- sample(nrow(D))
      hc <- fastcluster::hclust(stats::as.dist(D[p, p]), method)
      expect_same_fit(tree_eval(hc, D), tree_eval_reference(hc, D),
                      paste(nm, method))
      # a dist object or the shuffled matrix give the same fit
      expect_same_fit(tree_eval(hc, stats::as.dist(D)),
                      tree_eval_reference(hc, D), paste(nm, method, "dist"))
      expect_same_fit(tree_eval(hc, D[p, p]), tree_eval_reference(hc, D),
                      paste(nm, method, "shuffled"))
    }
  }
})

test_that("tree_eval stays exact on dissimilarities packed near 1", {
  # The case where the one-pass formula loses its digits: values whose spread
  # is tiny compared to their mean (here 0.999 to 1, as with high turnover).
  set.seed(3)
  x <- matrix(stats::runif(300 * 10), 300)
  rownames(x) <- paste0("site", seq_len(300))
  E <- as.matrix(stats::dist(x))
  D <- 0.999 + 0.001 * E / max(E)
  diag(D) <- 0
  for (method in c("average", "complete", "ward.D2")) {
    hc <- fastcluster::hclust(stats::as.dist(D), method)
    expect_same_fit(tree_eval(hc, D), tree_eval_reference(hc, D),
                    paste("near 1", method))
  }
})

test_that("tree_eval matches the previous version on diana and phylo trees", {
  D <- sample_dissim(fishmat, 200)

  # diana (hclu_diana)
  dv <- cluster::diana(stats::as.dist(D), diss = TRUE)
  expect_same_fit(tree_eval(dv, stats::as.dist(D)),
                  tree_eval_reference(dv, stats::as.dist(D)), "diana")

  # phylo: the previous version read the patristic distances of the tree
  # (ape's cophenetic), which as.hclust() turns into heights
  hc <- fastcluster::hclust(stats::as.dist(D), "average")
  ph <- ape::as.phylo(hc)
  expect_same_fit(tree_eval(ph, D), tree_eval_reference(ph, D), "phylo")

  # the consensus tree of hclu_hierarclust: least-squares branch lengths, then
  # polytomies resolved with zero-length branches
  set.seed(4)
  trees <- lapply(1:10, function(i) {
    p <- sample(nrow(D))
    ape::as.phylo(fastcluster::hclust(stats::as.dist(D[p, p]), "average"))
  })
  cons <- ape::consensus(trees, p = 0.5)
  cons <- phangorn::nnls.tree(stats::as.dist(D), cons,
                              method = "ultrametric", trace = 0)
  cons <- ape::multi2di(cons)
  expect_same_fit(tree_eval(cons, D), tree_eval_reference(cons, D),
                  "consensus phylo")
})

test_that("tree_eval handles the edge cases as the previous version did", {
  # all dissimilarities equal: no correlation, with stats::cor()'s warning
  D <- matrix(0.5, 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  diag(D) <- 0
  hc <- fastcluster::hclust(stats::as.dist(D), "average")
  expect_warning(new <- tree_eval(hc, D), "standard deviation is zero")
  ref <- suppressWarnings(tree_eval_reference(hc, D))
  expect_true(is.na(new$cophcor) && is.na(ref$cophcor))
  expect_equal(new$msd, ref$msd)

  # the same with enough sites and a value such as 0.3, for which the sum over
  # all pairs divided by their number is not exactly 0.3
  for (n in c(100, 300)) {
    for (v in c(0.3, 0.1, 1/3)) {
      D <- matrix(v, n, n, dimnames = list(paste0("s", 1:n), paste0("s", 1:n)))
      diag(D) <- 0
      hc <- fastcluster::hclust(stats::as.dist(D), "average")
      expect_warning(new <- tree_eval(hc, D), "standard deviation is zero")
      expect_true(is.na(new$cophcor), label = paste("constant", n, v))
      expect_equal(new$msd, suppressWarnings(tree_eval_reference(hc, D))$msd)
    }
  }

  # all heights equal (a rake), the dissimilarities not
  D <- sample_dissim(fishmat, 300)
  hc <- fastcluster::hclust(stats::as.dist(D), "average")
  hc$height[] <- 0.3
  expect_warning(new <- tree_eval(hc, D), "standard deviation is zero")
  ref <- suppressWarnings(tree_eval_reference(hc, D))
  expect_true(is.na(new$cophcor) && is.na(ref$cophcor))
  expect_equal(new$msd, ref$msd, tolerance = 1e-15)

  # a perfect fit gives exactly 1, never above
  hc <- fastcluster::hclust(stats::as.dist(stats::cophenetic(
    fastcluster::hclust(stats::dist(1:20), "average"))), "average")
  D <- as.matrix(stats::cophenetic(hc))
  fit <- tree_eval(hc, D)
  expect_lte(fit$cophcor, 1)
  expect_equal(fit$cophcor, 1, tolerance = 1e-15)
  expect_equal(fit$msd, 0, tolerance = 1e-15)

  # two sites: a single pair, so no correlation
  D <- matrix(c(0, 0.3, 0.3, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  hc <- fastcluster::hclust(stats::as.dist(D), "average")
  expect_warning(fit <- tree_eval(hc, D), "standard deviation is zero")
  expect_true(is.na(fit$cophcor))
  expect_equal(fit$msd, 0)

  # without names, sites are matched by position
  D <- sample_dissim(fishmat, 50)
  hc <- fastcluster::hclust(stats::as.dist(D), "complete")
  D_unnamed <- unname(D)
  hc_unnamed <- hc
  hc_unnamed$labels <- NULL
  expect_equal(tree_eval(hc_unnamed, D_unnamed), tree_eval(hc, D))

  # a tree that does not match the matrix is an error, not a wrong number
  expect_error(tree_eval(hc, D[-1, -1]), "do not match")
  D_renamed <- D
  rownames(D_renamed)[1] <- "elsewhere"
  expect_error(tree_eval(hc, D_renamed), "do not match")

  # merges out of order (a node used before it is created, as the consensus
  # tree of hclu_hierarclust could give) are read correctly as long as the
  # last merge is the root...
  D <- matrix(0, 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  set.seed(6)
  D[lower.tri(D)] <- stats::runif(10)
  D <- D + t(D)
  in_order <- structure(list(merge = rbind(c(-1, -2), c(-4, -5), c(1, -3),
                                           c(2, 3)),
                             height = c(0.1, 0.2, 0.5, 0.9),
                             order = c(1, 2, 3, 4, 5), labels = letters[1:5]),
                        class = "hclust")
  out_of_order <- in_order
  out_of_order$merge <- rbind(c(-4, -5), c(3, -3), c(-1, -2), c(1, 2))
  out_of_order$height <- c(0.2, 0.5, 0.1, 0.9)
  expect_equal(tree_eval(out_of_order, D), tree_eval(in_order, D))
  expect_same_fit(tree_eval(out_of_order, D),
                  tree_eval_reference(in_order, D), "out of order")
  # ... and an error otherwise, rather than a wrong fit
  root_not_last <- in_order
  root_not_last$merge <- rbind(c(-1, -2), c(-4, -5), c(2, 4), c(1, -3))
  root_not_last$height <- c(0.1, 0.2, 0.9, 0.5)
  expect_error(tree_eval(root_not_last, D), "not its root")
})

test_that("tree_eval_cpp reads the same fit from the whole matrix", {
  D <- sample_dissim(vegemat, 200)
  set.seed(5)
  for (method in c("average", "complete", "ward.D2")) {
    sites <- sample.int(nrow(D), 80)
    sub <- D[sites, sites]
    hc <- fastcluster::hclust(stats::as.dist(sub), method)
    expect_equal(tree_eval_cpp(hc$merge, hc$height, D, sites),
                 tree_eval_cpp(hc$merge, hc$height, sub), tolerance = 1e-15)
    expect_same_fit(as.list(tree_eval_cpp(hc$merge, hc$height, D, sites)),
                    tree_eval_reference(hc, sub), paste("leaf_site", method))
  }
})

test_that("hclu_hierarclust and hclu_diana report the same fit as before", {
  # Every call to tree_eval() made by the clustering functions is checked
  # against the previous version on the very tree and matrix it was given.
  tree_eval_new <- tree_eval
  n_calls <- 0
  local_mocked_bindings(tree_eval = function(tree, dist_mat) {
    n_calls <<- n_calls + 1
    new <- tree_eval_new(tree, dist_mat)
    expect_same_fit(new, tree_eval_reference(tree, dist_mat),
                    class(tree)[1])
    new
  })
  dissim <- dissimilarity(fishmat[1:150, ], metric = "Simpson")
  reported_fit <- function(clust) {
    list(cophcor = clust$algorithm$final.tree.coph.cor,
         msd = clust$algorithm$final.tree.msd)
  }

  # a single tree
  clust <- hclu_hierarclust(dissim, randomize = FALSE, verbose = FALSE)
  expect_equal(n_calls, 1)

  # best of the randomized trees: one call per trial, and the best one is
  # reported
  clust <- hclu_hierarclust(dissim, n_runs = 10, seed = 1,
                            optimal_tree_method = "best",
                            keep_trials = "all", verbose = FALSE)
  expect_equal(n_calls, 1 + 10)
  best <- which.max(vapply(clust$algorithm$trials,
                           function(trial) trial$cophcor, numeric(1)))
  expect_equal(reported_fit(clust),
               clust$algorithm$trials[[best]][c("cophcor", "msd")])

  # consensus tree
  clust <- hclu_hierarclust(dissim, n_runs = 10, seed = 1,
                            optimal_tree_method = "consensus",
                            verbose = FALSE)
  expect_equal(n_calls, 1 + 10 + 11)

  # iterative hierarchical consensus tree, for a linkage scored by the
  # closed form and one scored by tree_eval_cpp
  for (method in c("average", "complete")) {
    clust <- hclu_hierarclust(dissim, method = method, n_runs = 10, seed = 1,
                              optimal_tree_method = "ihct", verbose = FALSE)
  }
  expect_equal(n_calls, 1 + 10 + 11 + 2)

  # diana
  clust <- hclu_diana(dissim, verbose = FALSE)
  expect_equal(n_calls, 1 + 10 + 11 + 2 + 1)
})
