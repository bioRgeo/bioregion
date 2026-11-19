# Inputs -----------------------------------------------------------------------
tree <- hclu_hierarclust(wf_dissim, optimal_tree_method = "best",
                         verbose = FALSE)
tree_1 <- cut_tree(tree, n_clust = 5, verbose = FALSE)
tree_3 <- cut_tree(tree, n_clust = c(3, 5, 10), verbose = FALSE)
dyn <- cut_tree(tree, dynamic_tree_cut = TRUE, verbose = FALSE)
diana <- suppressWarnings(hclu_diana(wf_dissim, n_clust = c(3, 5),
                                     verbose = FALSE))
pam <- nhclu_pam(wf_dissim, n_clust = c(3, 6), seed = 1)
greedy <- netclu_greedy(wf_sim)
greedy_bip <- netclu_greedy(mat_to_net(wf_comat), weight = FALSE,
                            bipartite = TRUE)

# Nested partitions without a tree, as the levels of netclu_infomap()
nested <- tree_3
nested$name <- "netclu_infomap"
nested$algorithm <- list()
nested$inputs$hierarchical <- TRUE

# Nested partitions of sites and species
nested_bip <- greedy_bip
nested_bip$clusters$K_coarse <- ifelse(
  greedy_bip$clusters[, 2] %in% unique(greedy_bip$clusters[, 2])[1:2],
  "1", "2")
nested_bip$clusters <- nested_bip$clusters[, c(1, 3, 2)]
attr(nested_bip$clusters, "node_type") <- attr(greedy_bip$clusters,
                                               "node_type")
nested_bip$inputs$hierarchical <- TRUE

# Tests for valid outputs ------------------------------------------------------
test_that("bars", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_null(plot(pam))
  expect_no_error(plot(pam, partitions = 2))
  expect_no_error(plot(pam, partitions = "K_6", horizontal = TRUE))
  expect_no_error(plot(pam, flip = TRUE, cex = 0.5))
  expect_no_error(plot(pam, items_per_unit = 3))
  expect_no_error(plot(pam, max_clusters = 3))

  # Axes, titles and scales of the panels
  expect_no_error(plot(pam, scales = "free"))
  expect_no_error(plot(pam, scales = "free", horizontal = TRUE))
  expect_no_error(plot(pam, axes = FALSE, labels = FALSE))
  expect_no_error(plot(pam, main = c("3 clusters", "6 clusters"),
                       xlab = "Clusters", ylab = "Sites"))
  expect_no_error(plot(pam, partitions = 1, main = "PAM", ylab = ""))

  # Sites without a cluster are left out
  expect_true(anyNA(greedy$clusters[, 2]))
  expect_no_error(plot(greedy))

  # Colors of a partition with more or with fewer clusters
  expect_no_error(plot(tree_3, plot_type = "bars", partitions = 1,
                       color_partition = 3))
  expect_no_error(plot(tree_3, plot_type = "bars", partitions = 3,
                       color_partition = "K_3"))

  # Sites and species
  expect_no_error(plot(greedy_bip))
  expect_no_error(plot(greedy_bip, bipartite_layout = "stack",
                       horizontal = TRUE))
  expect_no_error(plot(greedy_bip, bipartite_layout = "side"))
  expect_no_error(plot(greedy_bip, items_per_unit = c(2, 5),
                       max_clusters = 2))
  expect_no_error(plot(greedy_bip, items_per_unit = c(1, 3),
                       horizontal = TRUE))
})

test_that("a cluster named 'other' is kept apart from the 'other' bar", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  # The largest cluster is named "other", and the smallest clusters are
  # grouped in the last bar
  pam_other <- nhclu_pam(wf_dissim, n_clust = 10, seed = 1)
  cl <- pam_other$clusters[, 2]
  largest <- names(which.max(table(cl)))
  pam_other$clusters[cl == largest, 2] <- "other"
  expect_no_error(plot(pam_other, max_clusters = 4))

  cl <- pam_other$clusters[, 2]
  cnt <- bars_content(cl, cl, attr(pam_other$clusters, "node_type"),
                      cluster_colors(pam_other)[[1]],
                      k_unit = c(site = 1, species = 1), max_clusters = 4,
                      bipartite_layout = "stack", flip = FALSE)
  # 3 clusters, including the cluster "other", then the bar of the 7
  # smallest clusters
  expect_length(cnt$clusters, 4)
  expect_equal(cnt$is_other, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(cnt$n_other, 7)
  expect_true("other" %in% cnt$clusters[!cnt$is_other])
  expect_equal(sum(cnt$n_items[[4]]), sum(!(cl %in% cnt$clusters[1:3])))
})

test_that("trees", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  # Trees without clusters, with one or several partitions
  expect_no_error(plot(tree))
  expect_no_error(plot(tree_1))
  expect_no_error(plot(tree_3, labels = FALSE))
  expect_no_error(plot(tree_3, cex = 0.5))
  expect_null(plot(tree_1))

  # Dynamic tree cut, including clusters split across branches ("hybrid")
  expect_no_error(plot(dyn))
  dyn_hybrid <- cut_tree(tree, dynamic_tree_cut = TRUE,
                         dynamic_method = "hybrid",
                         dissimilarity = wf_dissim, verbose = FALSE)
  expect_no_error(plot(dyn_hybrid))
  expect_no_error(plot(dyn_hybrid, horizontal = TRUE))

  expect_no_error(plot(diana))

  # Branches colored by a partition chosen by its number or its name
  expect_no_error(plot(tree_3, branch_partition = 1))
  expect_no_error(plot(tree_3, branch_partition = "K_5"))

  # Selected partitions, mirrored and horizontal trees
  expect_no_error(plot(tree_3, partitions = 2:3))
  expect_no_error(plot(tree_3, flip = TRUE))
  expect_no_error(plot(tree_3, horizontal = TRUE))
  expect_no_error(plot(tree_3, horizontal = TRUE, flip = TRUE,
                       labels = FALSE))
  expect_no_error(plot(tree, horizontal = TRUE))

  # Axes and titles
  expect_no_error(plot(tree_3, main = "Tree", ylab = "Height"))
  expect_no_error(plot(tree_3, horizontal = TRUE, xlab = "Height",
                       axes = FALSE))

  # Collapsed clusters
  expect_no_error(plot(tree_3, collapse = TRUE))
  expect_no_error(plot(tree_3, partitions = 1:2, collapse = TRUE,
                       horizontal = TRUE))
  expect_no_error(plot(diana, collapse = TRUE))

  # One step per partition, with or without a tree
  expect_no_error(plot(tree_3, heights = "levels"))
  expect_no_error(plot(nested, plot_type = "tree"))
  expect_no_error(plot(nested, plot_type = "tree", collapse = TRUE))
  expect_no_error(plot(nested_bip, plot_type = "tree", labels = FALSE))
  expect_no_error(plot(greedy, plot_type = "tree", labels = FALSE))
})

test_that("trees of clusters", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_no_error(plot(tree_3, plot_type = "cluster_tree"))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree",
                       horizontal = TRUE, flip = TRUE))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree", heights = "rank"))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree",
                       heights = "levels"))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree", partitions = 1))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree",
                       items_per_unit = 2))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree",
                       color_partition = 1))
  expect_no_error(plot(tree_3, plot_type = "cluster_tree", labels = FALSE,
                       axes = FALSE, main = "Clusters"))

  # hclu_diana() has no cut heights
  expect_no_error(plot(diana, plot_type = "cluster_tree"))

  # Default plot for nested partitions without a tree
  expect_no_error(plot(nested))
  expect_no_error(plot(nested_bip, horizontal = TRUE))
})

test_that("nested partitions are rebuilt as a tree", {
  out <- nested_tree(tree_3$clusters, c("K_3", "K_5", "K_10"))

  expect_s3_class(out, "hclust")
  expect_setequal(out$labels, tree_3$clusters$ID)
  expect_equal(nrow(out$merge), nrow(tree_3$clusters) - 1)
  expect_setequal(out$order, seq_along(out$labels))
  # One step of height per partition, root at the top
  expect_setequal(unique(out$height), 1:4)

  # Each cluster is a single branch of the rebuilt tree
  cl <- as.character(tree_3$clusters$K_10[match(out$labels,
                                                tree_3$clusters$ID)])
  expect_true(clusters_are_branches(out, cl, merge_clusters(out, cl)))

  expect_true(partitions_nested(tree_3$clusters, c("K_3", "K_5", "K_10")))
  expect_false(partitions_nested(pam$clusters, c("K_3", "K_6")))
})

test_that("plot uses the colors stored by bioregion_colors()", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  colored <- bioregion_colors(tree_3)
  recut <- cut_tree(colored, n_clust = 8, verbose = FALSE)

  # bioregion_colors() is only called when the colors are missing or do not
  # match the partitions
  local_mocked_bindings(bioregion_colors = function(...) {
    stop("colors recomputed")
  })
  expect_no_error(plot(colored))
  expect_no_error(plot(colored, plot_type = "bars"))
  expect_error(plot(tree_3), "colors recomputed")
  expect_error(plot(recut), "colors recomputed")
})

test_that("arguments that do not apply to the plot are ignored with a warning", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  # Arguments of trees for bars, and of bars for trees
  expect_warning(plot(pam, collapse = TRUE),
                 "collapse does not apply to plot_type = 'bars' and is ignored.",
                 fixed = TRUE)
  expect_warning(plot(pam, heights = "levels", branch_partition = 1),
                 paste0("heights, branch_partition do not apply to ",
                        "plot_type = 'bars' and are ignored."),
                 fixed = TRUE)
  expect_warning(plot(tree_3, max_clusters = 5, scales = "free"),
                 paste0("max_clusters, scales do not apply to ",
                        "plot_type = 'tree' and are ignored."),
                 fixed = TRUE)
  expect_warning(plot(tree_3, color_partition = 1),
                 "color_partition does not apply to plot_type = 'tree'",
                 fixed = TRUE)
  expect_warning(plot(tree_3, items_per_unit = 2),
                 "items_per_unit does not apply to plot_type = 'tree'",
                 fixed = TRUE)
  expect_warning(plot(tree_3, plot_type = "cluster_tree", collapse = TRUE),
                 "collapse does not apply to plot_type = 'cluster_tree'",
                 fixed = TRUE)
  expect_warning(plot(nested, bipartite_layout = "side"),
                 "bipartite_layout does not apply to plot_type = 'cluster_tree'",
                 fixed = TRUE)

  # Bars: layout of sites and species without species, scales of one panel
  expect_warning(plot(pam, bipartite_layout = "side"),
                 "bipartite_layout only applies to bipartite clusterings",
                 fixed = TRUE)
  expect_warning(plot(pam, partitions = 1, scales = "free"),
                 "scales only applies to bars of several partitions",
                 fixed = TRUE)

  # No warning for the arguments that apply, or left to their default value
  expect_no_warning(plot(pam, scales = "free", max_clusters = 3,
                         items_per_unit = 2, color_partition = 1))
  expect_no_warning(plot(greedy_bip, bipartite_layout = "side"))
  expect_no_warning(plot(tree_3, heights = "levels", collapse = TRUE,
                         branch_partition = 1))
  expect_no_warning(plot(tree_3, plot_type = "cluster_tree",
                         color_partition = 1, items_per_unit = 2,
                         branch_partition = 2, heights = "rank"))
  expect_no_warning(plot(pam, collapse = FALSE, heights = "auto",
                         bipartite_layout = "mirror", scales = "fixed",
                         max_clusters = 25))
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_error(plot(tree_3, plot_type = "map"),
               "^plot_type must be one of")
  expect_error(plot(tree_3, heights = "log"),
               "^heights must be one of")
  expect_error(plot(tree_3, collapse = 1),
               "collapse must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, labels = 1),
               "labels must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, flip = "yes"),
               "flip must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, horizontal = NA),
               "horizontal must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, axes = "no"),
               "axes must be a boolean.", fixed = TRUE)
  expect_error(plot(pam, scales = "same"),
               "scales must be one of 'fixed' or 'free'.", fixed = TRUE)
  expect_error(plot(pam, main = 1),
               "main must be a character.", fixed = TRUE)
  expect_error(plot(pam, xlab = c("a", "b")),
               "xlab must be of length 1.", fixed = TRUE)
  expect_error(plot(pam, ylab = 1),
               "ylab must be a character.", fixed = TRUE)
  expect_error(plot(greedy_bip, bipartite_layout = "top"),
               "^bipartite_layout must be one of")
  expect_error(plot(pam, max_clusters = 1),
               "max_clusters must be higher than 1.", fixed = TRUE)
  expect_error(plot(pam, items_per_unit = "all"),
               "^items_per_unit must be 'auto' or")
  expect_error(plot(pam, items_per_unit = c(1, 2, 3)),
               "items_per_unit must be of length 1 or 2.", fixed = TRUE)
  expect_error(plot(pam, items_per_unit = 0),
               "^items_per_unit must be composed of values strictly higher")

  expect_error(plot(tree_3, partitions = 0),
               "partitions should be strictly higher than 0.", fixed = TRUE)
  expect_error(plot(tree_3, partitions = 4),
               "^partitions should be lower or equal to 3")
  expect_error(plot(tree_3, partitions = "K_4"),
               "^If partitions is a character, it should be a")
  expect_error(plot(tree_3, branch_partition = c(2, 3)),
               "branch_partition must be of length 1.", fixed = TRUE)
  expect_error(plot(tree_3, branch_partition = "zz"),
               "^If branch_partition is a character, it should be a")
  expect_error(plot(tree_3, branch_partition = 10),
               "^branch_partition should be lower or equal to")
  expect_error(plot(pam, color_partition = 1:2),
               "color_partition must be of length 1.", fixed = TRUE)

  # Uncut trees: only the tree itself can be plotted
  expect_error(plot(tree, plot_type = "bars"),
               "^No clusters have been generated for your hierarchical tree")
  expect_error(plot(tree, plot_type = "cluster_tree"),
               "^No clusters have been generated for your hierarchical tree")
  expect_error(plot(tree, partitions = 1),
               "partitions cannot be used before the tree is cut.",
               fixed = TRUE)
  expect_error(plot(tree, branch_partition = 1),
               "branch_partition cannot be used before the tree is cut.",
               fixed = TRUE)
  expect_error(plot(tree, color_partition = "K_3"),
               "color_partition cannot be used before the tree is cut.",
               fixed = TRUE)
  expect_error(plot(tree, heights = "levels"),
               "heights = 'levels' cannot be used before the tree is cut.",
               fixed = TRUE)
  expect_error(plot(tree, partitions = 1, heights = "levels"),
               "partitions, heights = 'levels' cannot be used before",
               fixed = TRUE)

  # x must be a valid bioregion.clusters object, as in map_bioregions()
  expect_error(plot.bioregion.clusters(pam$clusters),
               "x must be a bioregion.clusters object.", fixed = TRUE)
  altered <- pam
  altered$name <- NULL
  expect_error(plot(altered), "x is a bioregion.cluster object but it has ",
               fixed = TRUE)
  altered <- pam
  attr(altered$clusters, "node_type") <- NULL
  expect_error(plot(altered), "x is a bioregion.cluster object but it has ",
               fixed = TRUE)
  altered <- pam
  altered$clusters <- NA
  expect_error(plot(altered),
               "x does not have the expected type of 'clusters' slot.",
               fixed = TRUE)

  # Heights
  expect_error(plot(tree_3, heights = "rank"),
               "heights = 'rank' is only available for plot_type = ",
               fixed = TRUE)
  expect_error(plot(nested, heights = "tree"),
               "heights = 'tree' needs the tree of hclu_hierarclust()",
               fixed = TRUE)
  expect_error(plot(pam, plot_type = "tree"),
               "^The partitions are not nested")

  # Clusters split across the tree
  expect_error(plot(dyn, collapse = TRUE),
               "span distinct branches across the tree", fixed = TRUE)
  expect_error(plot(dyn, plot_type = "cluster_tree"),
               "span distinct branches across the tree", fixed = TRUE)
})
