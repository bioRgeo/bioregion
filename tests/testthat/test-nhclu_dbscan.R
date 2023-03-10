# Preamble code ----------------------------------------------------------------
dissim <- dissimilarity(vegemat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  clust1 <- nhclu_dbscan(dissim, index = "Simpson", plot = FALSE)
  
  expect_equal(clust1$name, "dbscan")
  expect_equal(ncol(clust1$clusters), 2L)
  expect_equal(length(unique(clust1$clusters$K_6)), 6L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    nhclu_dbscan(dissimilarity = "zz"),
    "dissimilarity is not a bioregion.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index)",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, index = "zz", plot = FALSE),
    "Argument index should be one of the column names of dissimilarity",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, minPts = "zz", plot = FALSE),
    "minPts must be a numeric.",
    fixed = TRUE)

  expect_error(
    nhclu_dbscan(dissim, eps = "zz", plot = FALSE),
    "eps must be a numeric.",
    fixed = TRUE)

  expect_error(
    nhclu_dbscan(dissim, plot = "zz"),
    "plot must be a Boolean.",
    fixed = TRUE)
})
