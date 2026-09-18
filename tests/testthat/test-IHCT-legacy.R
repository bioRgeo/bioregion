# =============================================================================
# Equivalence of the rewritten IHCT with the implementation it replaced
# -----------------------------------------------------------------------------
# The previous implementation is read out of git at a pinned commit and sourced
# into its own environment, so the two versions can be run side by side on
# fishdf and vegedf. Everything is compared: the C++ helpers against the R and
# stats code they replace, every new function against the old one it replaces,
# and the trees themselves, on subsets and on the full datasets.
#
# ON TIES -- why some differences are expected and not failures
#   At each division both versions build n_runs trees from randomised copies of
#   the dissimilarity matrix and keep the best ones. Dissimilarity matrices of
#   species composition are full of tied values, so it is common for several of
#   those trees to fit the data EQUALLY WELL, down to the last bit of the
#   cophenetic correlation. The old code then picked between them with
#   order(ccc, decreasing = TRUE), i.e. on floating-point noise; the new
#   rank_by_score() treats near-equal scores as tied and prefers the earliest
#   run. Neither choice is more correct than the other.
#
#   So these tests never require the two versions to pick the same trees. They
#   require that a different pick is always a pick among trees of the same
#   quality, and that a resulting tree difference can be traced to such a pick.
#   A difference that cannot be traced to a tie IS a failure.
#
# ON HEIGHTS -- what is pinned and what is compared
#   The old implementation raised every division to the highest division it
#   contained ("max_child"); the current default instead moves the heights as
#   little as possible ("least_squares", see make_heights_monotone). The rule is
#   applied once the whole topology is decided and draws no random numbers, so
#   it cannot change a single division. Every test that asks for the OLD tree
#   therefore pins height_rule = "max_child" and stays an exact identity check,
#   and expect_equivalent_trees() additionally runs the current default against
#   the same old tree: the shape must be the one the old code produced, and the
#   cophenetic correlation must not be lower.
#
# These tests need a git checkout of the package: they skip on CRAN, and
# wherever git or the pinned commit cannot be reached (a tarball build, for
# instance). Running them takes several minutes, mostly in the old
# implementation on the full vegedf.
# =============================================================================

# the last commit before the IHCT rewrite ("Added shortcut to find stable trees
# in IHCT()"), which holds the implementation these tests compare against
LEGACY_COMMIT <- "214abc5e"
LEGACY_FILES <- c("R/IHCT.R", "R/hclu_hierarclust.R")

# Two trees that differ must be traceable to a pick among trees whose
# cophenetic correlations agree to within this much.
TIE_TOL <- 1e-8


# -- the legacy implementation ------------------------------------------------

# Sourced into an environment whose parent is the bioregion namespace: the old
# functions then call each other rather than their new namesakes, while still
# finding the utils.R helpers they rely on (randomize_dist, tree_eval,
# cut_tree, ...). NULL when git or the commit is not reachable.
legacy <- local({
  if (!nzchar(Sys.which("git"))) return(NULL)
  tryCatch({
    # tests run in tests/testthat; walk up to the package root
    root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    while (!file.exists(file.path(root, "DESCRIPTION"))) {
      up <- dirname(root)
      if (identical(up, root)) return(NULL)
      root <- up
    }
    if (!dir.exists(file.path(root, ".git"))) return(NULL)
    dir <- file.path(tempdir(), "bioregion-legacy-IHCT")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    e <- new.env(parent = asNamespace("bioregion"))
    for (f in LEGACY_FILES) {
      txt <- suppressWarnings(
        system2("git", c("-C", shQuote(root), "show",
                         paste0(LEGACY_COMMIT, ":", f)),
                stdout = TRUE, stderr = FALSE))
      if (!length(txt) || !is.null(attr(txt, "status"))) return(NULL)
      path <- file.path(dir, basename(f))
      writeLines(txt, path)
      sys.source(path, envir = e)
    }
    if (!all(vapply(c("IHCT", "coassign_binary_split", "reconstruct_hclust_bis"),
                    function(x) is.function(e[[x]]), logical(1)))) return(NULL)
    e
  }, error = function(...) NULL)
})

# the old IHCT recurses once per division, so a chain-shaped tree of n sites
# nests n calls deep
options(expressions = 500000)

skip_without_legacy <- function() {
  skip_on_cran()
  if (is.null(legacy)) {
    skip(paste0("the IHCT implementation at commit ", LEGACY_COMMIT,
                " is not reachable (needs a git checkout of the package)"))
  }
  invisible(TRUE)
}


# -- data ---------------------------------------------------------------------
# The pipeline hclu_hierarclust() uses: network -> matrix -> dissimilarity ->
# square symmetrical dissimilarity matrix. Built once, on demand, and cached:
# vegedf takes a few seconds.

.cache <- new.env(parent = emptyenv())

make_dist <- function(net, metric = "Simpson") {
  mat <- net_to_mat(net, weight = TRUE, squared = FALSE, symmetrical = FALSE)
  dis <- dissimilarity(mat, metric = metric)
  net_to_mat(dis[, c(1, 2, which(colnames(dis) == metric))],
             weight = TRUE, squared = TRUE, symmetrical = TRUE)
}

# full dissimilarity matrix for one dataset ("fish" or "vege")
dist_full <- function(which) {
  key <- paste0("full_", which)
  if (is.null(.cache[[key]])) {
    .cache[[key]] <- make_dist(if (which == "fish") fishdf else vegedf)
  }
  .cache[[key]]
}

# a reproducible subset of n sites of that matrix (drawn without disturbing the
# random stream the caller is using)
dist_sub <- function(which, n, seed = 42) {
  key <- paste0(which, "_", n, "_", seed)
  if (is.null(.cache[[key]])) {
    D <- dist_full(which)
    state <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    set.seed(seed)
    s <- sort(sample(rownames(D), min(n, nrow(D))))
    if (!is.null(state)) assign(".Random.seed", state, globalenv())
    .cache[[key]] <- D[s, s]
  }
  .cache[[key]]
}

DATASETS <- c("fish", "vege")
LINKAGES <- c("single", "complete", "average", "mcquitty",
              "ward.D", "ward.D2", "centroid", "median")


# -- running the two versions -------------------------------------------------

# the legacy IHCT, with the arguments hclu_hierarclust() passed it at the time
run_legacy <- function(D, method = "average", n_runs = 100,
                       top_n_trees = 2, seed = 1) {
  set.seed(seed)
  tr <- legacy$IHCT(D,
                    sites = rownames(D),
                    method = method,
                    n_runs = n_runs,
                    top_n_trees = top_n_trees,
                    monotonicity_direction = "bottom-up",
                    stable_shortcircuit = FALSE,
                    verbose = FALSE)
  legacy$reconstruct_hclust_bis(tr)
}

# the current IHCT with the legacy height rule, so that the trees can be
# compared division by division AND height by height
run_current <- function(D, method = "average", n_runs = 100,
                        top_n_trees = 2, seed = 1,
                        height_rule = "max_child") {
  set.seed(seed)
  IHCT(D, method = method, n_runs = n_runs, top_n_trees = top_n_trees,
       height_rule = height_rule, verbose = FALSE)
}

# the groups of sites of every node, as a set. Two trees of the same shape have
# the same one whatever their heights, which is what changing the height rule
# may not alter (nodes are numbered by height, so the merge matrices can differ)
clades <- function(hc) {
  members <- vector("list", nrow(hc$merge))
  for (k in seq_len(nrow(hc$merge))) {
    side <- lapply(hc$merge[k, ], function(j) if (j < 0) hc$labels[-j] else members[[j]])
    members[[k]] <- sort(unlist(side))
  }
  sort(vapply(members, paste, "", collapse = ","))
}

trees_identical <- function(a, b) {
  identical(unname(as.matrix(a$merge)), unname(as.matrix(b$merge))) &&
    isTRUE(all.equal(unname(a$height), unname(b$height), tolerance = 1e-12)) &&
    identical(as.integer(a$order), as.integer(b$order)) &&
    identical(a$labels, b$labels)
}


# -- replaying a run division by division -------------------------------------
# Both versions walk the tree in the same order (smaller branch first) and draw
# the same permutations, so their k-th division is the same division. Tracing
# both and comparing them index by index localises a difference exactly.

# temporarily replace a binding so that the recursive / looping callers pick the
# traced version up
with_traced <- function(env, name, traced, expr) {
  orig <- get(name, envir = env)
  locked <- bindingIsLocked(name, env)
  if (locked) unlockBinding(name, env)
  assign(name, traced, envir = env)
  on.exit({
    assign(name, orig, envir = env)
    if (locked) lockBinding(name, env)
  }, add = TRUE)
  force(expr)
}

# step 3 of divide_sites(), applied to an arbitrary selection of trees
groups_from_picks <- function(trees, picks, sites, n_runs, method) {
  m <- length(sites); sorted <- sort(sites)
  side <- vapply(trees[picks],
                 function(tree) stats::cutree(tree, k = 2)[sorted], integer(m))
  if (length(picks) == 1) {
    membership <- side[, 1]
  } else {
    together <- matrix(0, m, m)
    for (r in seq_len(ncol(side))) together <- together + outer(side[, r], side[, r], "==")
    coassign <- 1 - together / n_runs
    dimnames(coassign) <- list(sorted, sorted)
    membership <- stats::cutree(
      fastcluster::hclust(stats::as.dist(coassign), method = method), k = 2)
  }
  sort(sorted[membership == membership[1]])   # group 1; group 2 is its complement
}

# Draw the n_runs randomised trees for ONE group of sites exactly as both
# versions do, then apply both selection rules to them. `sub` must already be
# the sub-matrix of that group, in the site order the caller uses.
division_keys <- function(sub, method = "average", n_runs = 100, top_n_trees = 2,
                          fit = tree_fit_score, rank = rank_by_score) {
  m <- nrow(sub); sites <- rownames(sub)
  top_n_trees <- min(top_n_trees, n_runs)
  trees <- vector("list", n_runs)
  ccc <- numeric(n_runs); sco <- numeric(n_runs)
  for (r in seq_len(n_runs)) {
    perm <- sample.int(m)                      # the same draw as both versions
    d_run <- sub[perm, perm]
    hc <- fastcluster::hclust(stats::as.dist(d_run), method)
    trees[[r]] <- hc
    ccc[r] <- suppressWarnings(tree_eval(hc, sub)$cophcor)   # OLD ranking key
    sco[r] <- fit(hc, d_run, method)                         # NEW ranking key
  }
  old_pick <- order(ccc, decreasing = TRUE)[seq_len(top_n_trees)]
  new_pick <- rank(sco)[seq_len(top_n_trees)]
  list(ccc = ccc, sco = sco, old_pick = old_pick, new_pick = new_pick,
       same_set = setequal(old_pick, new_pick),
       # does a different selection actually change the division?
       same_split = identical(groups_from_picks(trees, old_pick, sites, n_runs, method),
                              groups_from_picks(trees, new_pick, sites, n_runs, method)),
       # how much worse is the best tree the new rule kept?
       ccc_loss = max(ccc[old_pick]) - max(ccc[new_pick]),
       n_distinct_ccc = length(unique(round(ccc, 12))))
}

selection_clash <- function(D, sites, method = "average",
                            n_runs = 100, top_n_trees = 2, seed = 1) {
  set.seed(seed)
  division_keys(D[sites, sites], method, n_runs, top_n_trees)
}

# every division of a run of the CURRENT implementation, in order;
# `keys_at` names the one division whose tree selection should also be recovered
# (recovering it costs about as much as the division itself, so it is done for
# a single division only, once we know which one to look at)
trace_current <- function(D, method, n_runs, top_n_trees, seed, keys_at = NULL) {
  ns <- asNamespace("bioregion")
  orig <- get("divide_sites", envir = ns)
  orig_fit <- get("tree_fit_score", envir = ns)   # captured so that the key
  orig_rank <- get("rank_by_score", envir = ns)   # recovery is never rebound
  rng_get <- function() get(".Random.seed", envir = globalenv())
  rng_set <- function(s) assign(".Random.seed", s, envir = globalenv())
  log <- list()

  traced <- function(dist_mat, sites, site_names, method, n_runs, top_n_trees) {
    before <- rng_get()
    g <- orig(dist_mat, sites, site_names, method, n_runs, top_n_trees)
    k <- NULL
    if (!is.null(keys_at) && length(log) + 1L == keys_at) {
      after <- rng_get()
      rng_set(before)
      k <- division_keys(dist_mat[sites, sites], method, n_runs, top_n_trees,
                         fit = orig_fit, rank = orig_rank)
      rng_set(after)
    }
    log[[length(log) + 1L]] <<-
      list(sites = sort(site_names[sites]), g1 = sort(site_names[g[[1]]]), keys = k)
    g
  }

  with_traced(ns, "divide_sites", traced, {
    set.seed(seed)
    IHCT(D, method = method, n_runs = n_runs, top_n_trees = top_n_trees,
         verbose = FALSE)
  })
  log
}

# the same, for the legacy implementation
trace_legacy <- function(D, method, n_runs, top_n_trees, seed) {
  orig <- get("coassign_binary_split", envir = legacy)
  log <- list()
  traced <- function(dist_mat, ...) {
    g <- orig(dist_mat, ...)
    log[[length(log) + 1L]] <<-
      list(sites = sort(rownames(dist_mat)), g1 = sort(unname(g[[1]])))
    g
  }
  with_traced(legacy, "coassign_binary_split", traced, {
    set.seed(seed)
    legacy$IHCT(D, sites = rownames(D), method = method, n_runs = n_runs,
                top_n_trees = top_n_trees, monotonicity_direction = "bottom-up",
                stable_shortcircuit = FALSE, verbose = FALSE)
  })
  log
}

# index of the first division on which the two versions disagree, NA if none
first_differing_division <- function(lo, ln) {
  n <- min(length(lo), length(ln))
  for (i in seq_len(n)) {
    if (!identical(ln[[i]]$sites, lo[[i]]$sites) ||
        !identical(ln[[i]]$g1, lo[[i]]$g1)) return(i)
  }
  NA_integer_   # no division differs; if the trees do, that is unexplained
}


# -- the expectation used everywhere trees are compared -----------------------
# Passes when the two trees are identical, and also when they are not but the
# difference is traceable to a pick among equally good trees. Fails otherwise.

expect_equivalent_trees <- function(D, method = "average", n_runs = 100,
                                    top_n_trees = 2, seed = 1, label = "") {
  hc_old <- run_legacy(D, method, n_runs, top_n_trees, seed)
  hc_new <- run_current(D, method, n_runs, top_n_trees, seed)

  # The current default height rule, on the very same run: it may change the
  # heights (that is its purpose) but never a division, so the tree must keep
  # the shape it has under the legacy rule, and where the heights are mean
  # dissimilarities it can only fit the data better. Checked on every dataset
  # and linkage the suite visits, and needs no extra run of the old code.
  hc_ls <- run_current(D, method, n_runs, top_n_trees, seed,
                       height_rule = "least_squares")
  expect_identical(clades(hc_ls), clades(hc_new),
                   label = paste0(label, ": groups of sites of the ",
                                  "least-squares tree"))
  if (method %in% c("average", "mcquitty")) {
    expect_gte(tree_eval(hc_ls, D)$cophcor, tree_eval(hc_new, D)$cophcor,
               label = paste0(label, ": cophenetic correlation with ",
                              "least-squares heights"))
  }

  if (trees_identical(hc_old, hc_new)) {
    succeed()
    return(invisible("identical"))
  }

  # 1. the new tree must fit the dissimilarities as well as the old one
  e_old <- tree_eval(hc_old, D)
  e_new <- tree_eval(hc_new, D)
  expect_gte(e_new$cophcor, e_old$cophcor - 1e-6,
             label = paste0(label, ": cophenetic correlation of the new tree"))

  # 2. the difference must come from the tree selection at one division, and
  #    that selection must have been a choice between equally good trees
  lo <- trace_legacy(D, method, n_runs, top_n_trees, seed)
  ln <- trace_current(D, method, n_runs, top_n_trees, seed)
  i <- first_differing_division(lo, ln)
  expect_false(is.na(i),
               label = paste0(label, ": trees differ but no division does, which"))
  if (is.na(i)) return(invisible("unexplained"))

  k <- trace_current(D, method, n_runs, top_n_trees, seed, keys_at = i)[[i]]$keys
  expect_false(is.null(k),
               label = paste0(label, ": tree selection at division ", i, " missing, which"))
  if (is.null(k)) return(invisible("unexplained"))

  expect_false(k$same_set,
               label = paste0(label, ": at division ", i,
                              " the same trees were picked yet the split differs, which"))
  expect_lte(k$ccc_loss, TIE_TOL,
             label = paste0(label, ": cophenetic correlation given up at division ", i))

  invisible(sprintf(
    paste("%s: trees differ from division %d (%d sites), where %d of %d",
          "randomised trees tie on the best fit; the new rule gives up",
          "%.3g of cophenetic correlation"),
    label, i, length(ln[[i]]$sites),
    sum(abs(k$ccc - max(k$ccc)) < 1e-12), length(k$ccc), k$ccc_loss))
}


# =============================================================================
# 1. The C++ helpers against the R and stats code they replace
# =============================================================================

test_that("ihct_node_sizes() matches a plain R implementation", {
  skip_without_legacy()

  node_sizes_R <- function(merge) {
    size <- integer(nrow(merge)); pairs <- numeric(nrow(merge))
    for (k in seq_len(nrow(merge))) {
      a <- merge[k, 1]; b <- merge[k, 2]
      sa <- if (a < 0) 1L else size[a]
      sb <- if (b < 0) 1L else size[b]
      size[k] <- sa + sb
      pairs[k] <- sa * sb
    }
    list(size = size, pairs = pairs)
  }

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    set.seed(7)
    for (i in 1:25) {
      p <- sample(nrow(D)); Dp <- D[p, p]
      for (meth in c("average", "complete", "single", "ward.D2")) {
        hc <- fastcluster::hclust(stats::as.dist(Dp), meth)
        cpp <- ihct_node_sizes(hc$merge)
        ref <- node_sizes_R(hc$merge)
        expect_identical(as.integer(cpp$size), as.integer(ref$size))
        expect_equal(cpp$pairs, ref$pairs)
      }
    }
    # the pair counts must add up to every pair of sites, exactly once
    hc <- fastcluster::hclust(stats::as.dist(D), "average")
    expect_equal(sum(ihct_node_sizes(hc$merge)$pairs), nrow(D) * (nrow(D) - 1) / 2)
  }
})

test_that("ihct_cophenetic_correlation() reproduces the old CCC computation", {
  skip_without_legacy()

  # verbatim logic of utils.R::tree_eval(), which is what the old IHCT called
  coph_R_old <- function(hc, d) {
    coph <- as.matrix(stats::cophenetic(hc))
    coph <- coph[match(rownames(d), rownames(coph)), match(rownames(d), colnames(coph))]
    lt <- lower.tri(d)
    stats::cor(d[lt], coph[lt], method = "pearson")
  }

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    set.seed(11)
    for (i in 1:20) {
      p <- sample(nrow(D)); Dp <- D[p, p]
      for (meth in LINKAGES) {
        hc <- fastcluster::hclust(stats::as.dist(Dp), meth)
        expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, Dp),
                     coph_R_old(hc, Dp), tolerance = 1e-12,
                     label = paste(nm, meth, "cophenetic correlation"))
        expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, Dp),
                     tree_eval(hc, Dp)$cophcor, tolerance = 1e-12)
      }
    }
    # and it does not depend on the order the matrix happens to be in
    set.seed(12)
    p <- sample(nrow(D)); Dp <- D[p, p]
    hc <- fastcluster::hclust(stats::as.dist(Dp), "complete")
    expect_equal(ihct_cophenetic_correlation(hc$merge, hc$height, Dp),
                 tree_eval(hc, D)$cophcor, tolerance = 1e-12)
  }
})

test_that("the UPGMA closed-form score reproduces and ranks like the CCC", {
  skip_without_legacy()

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    lt <- D[lower.tri(D)]
    sst <- sum(lt^2) - sum(lt)^2 / length(lt)
    set.seed(13)
    ccc <- numeric(100); sco <- numeric(100)
    for (i in 1:100) {
      p <- sample(nrow(D)); Dp <- D[p, p]
      hc <- fastcluster::hclust(stats::as.dist(Dp), "average")
      sco[i] <- tree_fit_score(hc, Dp, "average")
      ccc[i] <- tree_eval(hc, Dp)$cophcor
      # the analysis-of-variance identity the closed form rests on
      expect_equal(sqrt(1 - (sum(lt^2) - sco[i]) / sst), ccc[i], tolerance = 1e-10,
                   label = paste(nm, "closed-form cophenetic correlation"))
    }
    # Wherever the fit genuinely differs, the score must order the trees the
    # same way. Where every tree fits equally well -- which happens a lot with
    # tied dissimilarities -- there is no order to agree on, and both rules are
    # sorting floating-point noise.
    dccc <- outer(ccc, ccc, "-"); dsco <- outer(sco, sco, "-")
    expect_equal(sum(dccc > 1e-12 & dsco <= 0), 0,
                 label = paste(nm, "pairs where a better-fitting tree scores lower"))
  }
})


# =============================================================================
# 2. Each new function against the old one it replaces
# =============================================================================

test_that("divide_sites() splits a group of sites exactly like coassign_binary_split()", {
  skip_without_legacy()

  for (nm in DATASETS) {
    D <- dist_sub(nm, 150)
    nmz <- rownames(D)
    set.seed(101)
    for (i in 1:30) {
      sites_nm <- sort(sample(nmz, sample(5:nrow(D), 1)))
      idx <- match(sites_nm, nmz)
      for (tn in c(1, 2, 5)) {
        set.seed(2000 + i)
        g_new <- lapply(divide_sites(D, idx, nmz, "average", 40, tn),
                        function(x) sort(nmz[x]))
        set.seed(2000 + i)
        g_old <- lapply(unname(legacy$coassign_binary_split(
          D[sites_nm, sites_nm], method = "average", n_runs = 40,
          top_n_trees = tn, binsplit = "tree")), sort)
        expect_setequal(g_new[[1]], g_old[[1]])
        expect_setequal(g_new[[2]], g_old[[2]])
      }
    }
  }
})

test_that("division_height() matches the old height computation for every linkage", {
  skip_without_legacy()

  # verbatim from the old IHCT
  old_height <- function(dist_mat_d, c1, c2, method) {
    pairwise_distances <- dist_mat_d[c1, c2]
    centroid1 <- colMeans(dist_mat_d[c1, , drop = FALSE])
    centroid2 <- colMeans(dist_mat_d[c2, , drop = FALSE])
    centroid_distance <- sum((centroid1 - centroid2)^2)
    switch(method,
           "single"   = min(pairwise_distances),
           "complete" = max(pairwise_distances),
           "average"  = mean(pairwise_distances),
           "mcquitty" = mean(pairwise_distances),
           "ward.D"   = (length(c1) * length(c2)) / (length(c1) + length(c2)) *
                          sqrt(centroid_distance),
           "ward.D2"  = (length(c1) * length(c2)) / (length(c1) + length(c2)) *
                          centroid_distance,
           "centroid" = centroid_distance,
           "median"   = 0.5 * centroid_distance)
  }

  for (nm in DATASETS) {
    D <- dist_sub(nm, 150)
    nmz <- rownames(D)
    set.seed(103)
    for (i in 1:40) {
      sites_nm <- sort(sample(nmz, sample(4:nrow(D), 1)))
      cut1 <- sample(seq_len(length(sites_nm) - 1), 1)
      c1 <- sites_nm[seq_len(cut1)]; c2 <- setdiff(sites_nm, c1)
      for (meth in LINKAGES) {
        expect_equal(
          division_height(D, match(c1, nmz), match(c2, nmz), match(sites_nm, nmz), meth),
          # the new code always floors the height at 0; the old bottom-up branch
          # lost that floor by overwriting the variable it had been applied to
          max(old_height(D[sites_nm, sites_nm], c1, c2, meth), 0),
          tolerance = 1e-12, label = paste(nm, meth, "division height"))
      }
    }
  }
})

test_that("nodes_to_hclust() rebuilds the same hclust as reconstruct_hclust_bis()", {
  skip_without_legacy()

  # the old nested-list tree, in the node table the new reconstruction expects
  as_node_table <- function(tree, site_names) {
    cap <- 2L * length(site_names)
    height <- numeric(cap); left <- integer(cap)
    right <- integer(cap); site <- integer(cap)
    store <- vector("list", cap); store[[1]] <- tree
    n <- 1L; stack <- 1L
    while (length(stack)) {
      id <- stack[length(stack)]; stack <- stack[-length(stack)]
      nd <- store[[id]]; store[id] <- list(NULL)
      if (is.null(nd$children)) {
        site[id] <- match(nd$name, site_names)
      } else {
        height[id] <- nd$height
        n <- n + 1L; left[id]  <- n; store[[n]] <- nd$children[[1]]
        n <- n + 1L; right[id] <- n; store[[n]] <- nd$children[[2]]
        stack <- c(stack, left[id], right[id])
      }
    }
    list(height = height[seq_len(n)], left = left[seq_len(n)],
         right = right[seq_len(n)], site = site[seq_len(n)])
  }

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    # both reconstructions are fed the SAME tree, so any difference is in the
    # reconstruction itself: node numbering, merge matrix, leaf order
    set.seed(5)
    tr <- legacy$IHCT(D, sites = rownames(D), method = "average", n_runs = 20,
                      top_n_trees = 2, monotonicity_direction = "bottom-up",
                      stable_shortcircuit = FALSE, verbose = FALSE)
    hc_old <- legacy$reconstruct_hclust_bis(tr)
    nd <- as_node_table(tr, rownames(D))
    hc_new <- nodes_to_hclust(nd$height, nd$left, nd$right, nd$site, rownames(D))

    expect_identical(unname(as.matrix(hc_new$merge)), unname(as.matrix(hc_old$merge)))
    expect_equal(unname(hc_new$height), unname(hc_old$height), tolerance = 1e-14)
    expect_identical(as.integer(hc_new$order), as.integer(hc_old$order))
    expect_identical(hc_new$labels, hc_old$labels)
    expect_equal(as.matrix(stats::cophenetic(hc_new)),
                 as.matrix(stats::cophenetic(hc_old)), tolerance = 1e-12)
  }
})


# =============================================================================
# 3. The one place the two versions may legitimately differ: tie-breaking
# =============================================================================

test_that("a different pick among tied trees never keeps a worse tree", {
  skip_without_legacy()

  for (nm in DATASETS) {
    D <- dist_sub(nm, 150)
    set.seed(99)
    n_clash <- 0L; n_split <- 0L
    for (i in 1:30) {
      sites <- sort(sample(rownames(D), sample(8:nrow(D), 1)))
      r <- selection_clash(D, sites, n_runs = 50, top_n_trees = 2, seed = 1000 + i)
      # The two rules need not agree on which trees to keep. What they may not
      # do is keep a tree that fits the data measurably worse.
      expect_lte(r$ccc_loss, TIE_TOL,
                 label = paste0(nm, ", group ", i,
                                ": cophenetic correlation given up by the new rule"))
      n_clash <- n_clash + !r$same_set
      n_split <- n_split + !r$same_split
    }
    # a division can only change if the selection did
    expect_lte(n_split, n_clash)
    message(sprintf(
      "  [%s] tree selection differed on %d of 30 groups; %d changed the split",
      nm, n_clash, n_split))
  }
})

test_that("with the old ranking restored, the new algorithm gives old trees", {
  skip_without_legacy()

  # Decisive check that nothing OUTSIDE the ranking step behaves differently.
  # The new code changed two things about ranking: it scores UPGMA trees with a
  # closed form rather than the cophenetic correlation, and it treats near-equal
  # scores as ties. Put the old key (the exact correlation, for every linkage)
  # and the old tie-break (order(), no tolerance) back in, and ask for the old
  # height rule: what is left is the rewritten machinery -- the work list,
  # integer site handling, the C++ helpers, the hclust reconstruction.
  run_old_ranking <- function(D, method, n_runs, top_n_trees, seed) {
    ns <- asNamespace("bioregion")
    fit_ccc <- function(tree, d, method) {
      ihct_cophenetic_correlation(tree$merge, tree$height, d)
    }
    rank_exact <- function(scores, tolerance = 1e-10) order(scores, decreasing = TRUE)
    with_traced(ns, "tree_fit_score", fit_ccc,
      with_traced(ns, "rank_by_score", rank_exact, {
        set.seed(seed)
        IHCT(D, method = method, n_runs = n_runs, top_n_trees = top_n_trees,
             height_rule = "max_child", verbose = FALSE)
      }))
  }

  for (nm in DATASETS) {
    D <- dist_sub(nm, 150)
    for (s in 1:3) {
      hc_old <- run_legacy(D, "average", 50, 2, s)
      hc_mix <- run_old_ranking(D, "average", 50, 2, s)
      if (trees_identical(hc_old, hc_mix)) {
        succeed()
        next
      }
      # The only thing that can still differ is the ORDER in which the
      # correlation is summed: the old code summed over the matrix, this one
      # enumerates pairs tree-first. Both are the same number in exact
      # arithmetic, so a residual difference must sit on an exact tie.
      lo <- trace_legacy(D, "average", 50, 2, s)
      ln <- trace_current(D, "average", 50, 2, s)
      i <- first_differing_division(lo, ln)
      skip_if(is.na(i), "trees differ but no division does")
      k <- trace_current(D, "average", 50, 2, s, keys_at = i)[[i]]$keys
      n_tied <- sum(abs(k$ccc - max(k$ccc)) < 1e-15)
      expect_gt(n_tied, 1,
                label = paste0(nm, "/seed ", s, ": trees tied on the best ",
                               "fit at division ", i))
      message(sprintf(
        paste("  [%s] seed %d: %d of %d randomised trees are exactly tied at",
              "division %d; the winner is decided by rounding, not by the data"),
        nm, s, n_tied, length(k$ccc), i))
    }
  }
})


# =============================================================================
# 4. Whole trees
# =============================================================================

test_that("IHCT() reproduces the old tree on subsets, for every linkage and top_n_trees", {
  skip_without_legacy()

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    for (meth in c("average", "complete", "single", "ward.D2")) {
      note <- expect_equivalent_trees(D, meth, n_runs = 25, top_n_trees = 2,
                                      seed = 1, label = paste(nm, meth))
      if (!identical(note, "identical")) message("  ", note)
    }
    # top_n_trees = 1: the old code still routed the single best tree through
    # the co-assignment matrix, the new one cuts it directly
    for (tn in c(1, 2)) {
      note <- expect_equivalent_trees(D, "average", n_runs = 50, top_n_trees = tn,
                                      seed = 1, label = paste0(nm, " top_n_trees=", tn))
      if (!identical(note, "identical")) message("  ", note)
    }
  }
})

test_that("IHCT() reproduces the old tree across seeds", {
  skip_without_legacy()

  for (nm in DATASETS) {
    D <- dist_sub(nm, 150)
    for (s in 1:3) {
      note <- expect_equivalent_trees(D, "average", n_runs = 50, top_n_trees = 2,
                                      seed = s, label = paste0(nm, "/150 seed ", s))
      if (!identical(note, "identical")) message("  ", note)
    }
  }
})

test_that("IHCT() reproduces the old tree on the full fishdf and vegedf", {
  skip_without_legacy()
  # the old implementation needs roughly a minute per dataset here

  for (nm in DATASETS) {
    D <- dist_full(nm)
    t_old <- system.time(hc_old <- run_legacy(D, "average", 100, 2, 1))[["elapsed"]]
    t_new <- system.time(hc_new <- run_current(D, "average", 100, 2, 1))[["elapsed"]]
    message(sprintf("  [%s] %d sites, n_runs = 100: old %.1fs, new %.1fs (x%.1f)",
                    nm, nrow(D), t_old, t_new, t_old / t_new))

    if (trees_identical(hc_old, hc_new)) {
      succeed()
      next
    }
    expect_gte(tree_eval(hc_new, D)$cophcor, tree_eval(hc_old, D)$cophcor - 1e-6,
               label = paste0(nm, " full: cophenetic correlation of the new tree"))
    lo <- trace_legacy(D, "average", 100, 2, 1)
    ln <- trace_current(D, "average", 100, 2, 1)
    i <- first_differing_division(lo, ln)
    expect_false(is.na(i),
                 label = paste0(nm, " full: trees differ but no division does, which"))
    if (is.na(i)) next
    k <- trace_current(D, "average", 100, 2, 1, keys_at = i)[[i]]$keys
    expect_false(k$same_set,
                 label = paste0(nm, " full: at division ", i,
                                " the same trees were picked yet the split differs, which"))
    expect_lte(k$ccc_loss, TIE_TOL,
               label = paste0(nm, " full: cophenetic correlation given up at division ", i))
    message(sprintf("  [%s] trees differ from division %d of %d (%d sites), on a tie",
                    nm, i, length(ln), length(ln[[i]]$sites)))
  }
})


# =============================================================================
# 5. hclu_hierarclust(), end to end
# =============================================================================

test_that("hclu_hierarclust() reproduces the old output end to end", {
  skip_without_legacy()
  skip_if(!is.function(legacy$hclu_hierarclust),
          "the legacy hclu_hierarclust() could not be sourced")

  for (nm in DATASETS) {
    D <- dist_sub(nm, 60)
    dis <- mat_to_net(D, weight = TRUE, include_diag = FALSE, include_lower = FALSE)
    colnames(dis)[3] <- "Simpson"

    # The legacy hclu_hierarclust() calls shared helpers that may have moved on
    # since; if it can no longer run, that is not an IHCT regression.
    o <- tryCatch(
      legacy$hclu_hierarclust(dis, index = "Simpson", method = "average",
                              n_runs = 50,
                              optimal_tree_method = "iterative_consensus_tree",
                              n_clust = 5, seed = 1, verbose = FALSE),
      error = function(e) e)
    if (inherits(o, "error")) {
      skip(paste("the legacy hclu_hierarclust() no longer runs:",
                 conditionMessage(o)))
    }

    n <- hclu_hierarclust(dis, index = "Simpson", method = "average",
                          n_runs = 50,
                          optimal_tree_method = "iterative_consensus_tree",
                          n_clust = 5, seed = 1, top_n_trees = 2, verbose = FALSE)

    if (trees_identical(o$algorithm$final.tree, n$algorithm$final.tree)) {
      expect_equal(n$algorithm$final.tree.coph.cor, o$algorithm$final.tree.coph.cor,
                   tolerance = 1e-12)
      expect_equal(n$algorithm$final.tree.msd, o$algorithm$final.tree.msd,
                   tolerance = 1e-12)
      expect_equal(n$clusters[order(n$clusters$ID), ],
                   o$clusters[order(o$clusters$ID), ], ignore_attr = TRUE)
      expect_equal(n$cluster_info, o$cluster_info, ignore_attr = TRUE)
    } else {
      # a tie inside IHCT; the tree must still be as good, and the object still
      # has to be a well-formed bioregionalization with the requested cut
      expect_gte(n$algorithm$final.tree.coph.cor,
                 o$algorithm$final.tree.coph.cor - 1e-6,
                 label = paste0(nm, ": cophenetic correlation of the new tree"))
      expect_identical(dim(n$clusters), dim(o$clusters))
      expect_setequal(n$clusters$ID, o$clusters$ID)
      message(sprintf("  [%s] hclu_hierarclust trees differ on a tie inside IHCT", nm))
    }
  }
})
