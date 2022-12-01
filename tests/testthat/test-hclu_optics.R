# Preamble code ----------------------------------------------------------------
dissim <- dissimilarity(vegemat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  tree1 <- hclu_optics(dissim, index = "Simpson")
  
  expect_equal(nrow(tree1$clusters), 715L)
  expect_equal(ncol(tree1$clusters), 2L)
  expect_equal(length(unique(tree1$clusters)), 2L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    hclu_optics(dissim, index = "zz"),
    "Argument index should be one of the column names of dissimilarity",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, minPts = "zz"),
    "minPts must be a positive integer, indicating the number of points
      to form a dense region (see dbscan::dbscan()).",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, eps = "zz"),
    "eps must be a positive integer, indicating the upper limit of the
         size of the epsilon neighborhood (see dbscan::optics()).",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, xi = "zz"),
    "xi must be a numeric in the ]0, 1[ interval
           (see dbscan::optics()).",
    fixed = TRUE)
})
