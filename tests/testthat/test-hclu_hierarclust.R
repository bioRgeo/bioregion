# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  tree1 <- hclu_hierarclust(dissim, n_clust = 5)
  
  expect_equal(nrow(tree1$clusters), 20L)
  expect_equal(ncol(tree1$clusters), 2L)
  expect_equal(length(unique(tree1$clusters$K_5)), 5L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    hclu_hierarclust(dissim, n_clust = "5"),
    "n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from partition_metrics()",
    fixed = TRUE)
})
