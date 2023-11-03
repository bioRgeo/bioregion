# Preamble code ----------------------------------------------------------------
dissim <- dissimilarity(fishmat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  clust1 <- hclu_diana(dissim, n_clust = 1, index = "Simpson")
  
  expect_equal(clust1$name, "divisive_hierarchical_clustering")
  expect_equal(ncol(clust1$clusters), 2L)
  # expect_equal(length(unique(clust1$clusters$K_6)), 6L)
})

# Tests for invalid inputs ----------------------------------------------------
# test_that("error messages with wrong inputs", {
#   expect_error(
#     nhclu_clara(dissimilarity = "zz"),
#     "dissimilarity is not a bioregion.pairwise.metric object, a
#            dissimilarity matrix (class dist) or a data.frame with at least 3
#            columns (site1, site2, and your dissimilarity index)",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, index = "zz"),
#     "Argument index should be one of the column names of dissimilarity",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = "zz"),
#     "n_clust must an integer or a vector of integers determining the
#            number of clusters.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, maxiter = "zz"),
#     "maxiter must be numeric.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, initializer = "zz"),
#     "initializer must be either 'BUILD' (used in classic PAM) or 'LAB'
#          (linear approximative BUILD).",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, fasttol = "zz"),
#     "fasttol must be numeric.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, numsamples = "zz"),
#     "numsamples must be numeric.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, sampling = "zz"),
#     "sampling must be numeric.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, independent = "zz"),
#     "independent must be a boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     nhclu_clara(dissim, n_clust = 2, seed = "zz"),
#     "seed must be numeric.",
#     fixed = TRUE)
#   
# })
