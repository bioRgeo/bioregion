# Preamble code ----------------------------------------------------------------
dissim <- dissimilarity(fishmat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  clust1 <- nhclu_pam(dissim, n_clust = 2:10, index = "Simpson")
  
  expect_equal(clust1$name, "pam")
  expect_equal(ncol(clust1$clusters), 10L)
  expect_equal(length(unique(clust1$clusters$K_6)), 6L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    nhclu_pam(dissimilarity = "zz"),
    "dissimilarity is not a bioRgeo.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index)",
    fixed = TRUE)
  
  expect_error(
    nhclu_pam(dissim, index = "zz"),
    "Argument index should be one of the column names of dissimilarity",
    fixed = TRUE)
  
  expect_error(
    nhclu_pam(dissim, n_clust = "zz"),
    "n_clust must an integer or a vector of integers determining the
           number of clusters.",
    fixed = TRUE)
  
  expect_error(
    nhclu_pam(dissim, variant = "zz"),
    "variant is a character string indicating. Available options are
         original, o_1, o_2, f_3, f_4, f_5 or faster.",
    fixed = TRUE)
  
  expect_error(
    nhclu_pam(dissim, nstart = "zz"),
    "nstart must be numeric.",
    fixed = TRUE)
  
  expect_error(
    nhclu_pam(dissim, cluster_only = "zz"),
    "cluster_only must be a boolean.",
    fixed = TRUE)
  
})
