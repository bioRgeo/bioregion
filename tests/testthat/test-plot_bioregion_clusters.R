# Inputs -----------------------------------------------------------------------
tree <- hclu_hierarclust(wf_dissim, optimal_tree_method = "best",
                         verbose = FALSE)
tree_1 <- cut_tree(tree, n_clust = 5, verbose = FALSE)
tree_3 <- cut_tree(tree, n_clust = c(3, 5, 10), verbose = FALSE)

# Tests for valid outputs ------------------------------------------------------
test_that("plot of hierarchical trees", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  # Trees without clusters, with one or several partitions
  expect_no_error(plot(tree))
  expect_no_error(plot(tree_1))
  expect_no_error(plot(tree_3, tree_labels = FALSE))
  expect_no_error(plot(tree_3, cex = 0.5))
  expect_null(plot(tree_1))

  # Dynamic tree cut, including clusters split across branches ("hybrid")
  dyn <- cut_tree(tree, dynamic_tree_cut = TRUE, verbose = FALSE)
  expect_no_error(plot(dyn))
  dyn_hybrid <- cut_tree(tree, dynamic_tree_cut = TRUE,
                         dynamic_method = "hybrid",
                         dissimilarity = wf_dissim, verbose = FALSE)
  expect_no_error(plot(dyn_hybrid))

  # DIANA
  diana <- hclu_diana(wf_dissim, n_clust = c(3, 5), verbose = FALSE)
  expect_no_error(suppressWarnings(plot(diana)))

  # Branches coloured by a partition chosen by its number or its name
  expect_no_error(plot(tree_3, tree_branch_partition = 2))
  expect_no_error(plot(tree_3,
                       tree_branch_partition = colnames(tree_3$clusters)[3]))

  # Mirrored and horizontal trees, with and without clusters
  expect_no_error(plot(tree_3, tree_flip = TRUE))
  expect_no_error(plot(tree_3, tree_horizontal = TRUE))
  expect_no_error(plot(tree_3, tree_horizontal = TRUE, tree_flip = TRUE,
                       tree_labels = FALSE))
  expect_no_error(plot(tree, tree_horizontal = TRUE))
  expect_no_error(plot(dyn_hybrid, tree_horizontal = TRUE))
})

test_that("plot uses the colours stored by bioregion_colors()", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  colored <- bioregion_colors(tree_3)
  recut <- cut_tree(colored, n_clust = 8, verbose = FALSE)

  # bioregion_colors() is only called when the colours are missing or do not
  # match the partitions
  local_mocked_bindings(bioregion_colors = function(...) {
    stop("colours recomputed")
  })
  expect_no_error(plot(colored))
  expect_error(plot(tree_3), "colours recomputed")
  expect_error(plot(recut), "colours recomputed")
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_error(plot(tree_3, tree_branch_partition = c(2, 3)),
               "tree_branch_partition must be of length 1.", fixed = TRUE)
  expect_error(plot(tree_3, tree_branch_partition = "zz"),
               "^If tree_branch_partition is a character, it should be a")
  expect_error(plot(tree_3, tree_branch_partition = 0),
               "tree_branch_partition should be strictly higher than 0.",
               fixed = TRUE)
  expect_error(plot(tree_3, tree_branch_partition = 10),
               "^tree_branch_partition should be lower or equal to")
  expect_error(plot(tree_3, tree_labels = 1),
               "tree_labels must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, tree_flip = "yes"),
               "tree_flip must be a boolean.", fixed = TRUE)
  expect_error(plot(tree_3, tree_horizontal = NA),
               "tree_horizontal must be a boolean.", fixed = TRUE)

  other <- tree_1
  other$name <- "nhclu_pam"
  expect_error(plot(other),
               "plot() is not available yet for the outputs of nhclu_pam().",
               fixed = TRUE)
})
