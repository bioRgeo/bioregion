# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

simil <- similarity(comat, metric = "all")
dissimilarity <- similarity_to_dissimilarity(simil)

tree1 <- hclu_hierarclust(dissimilarity, n_clust = 5)

# Tests for valid outputs -----------------------------------------------------
test_that("class list and bioregion.clusters", {
  tree2 <- cut_tree(tree1, cut_height = 0.05)
  
  expect_identical(class(tree2)[1], "bioregion.clusters")
  expect_identical(class(tree2)[2], "list")
  
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    cut_tree(tree1, n_clust = 1.5),
    "n_clust must an integer or a vector of integers determining the
             number of clusters.", fixed = TRUE)
})

test_that("error messages with wrong inputs 2", {
  expect_error(
    cut_tree(tree1, cut_height = 0.05, n_clust = 3),
    "Please provide either n_clust or cut_height, but not both at the
           same time.", fixed = TRUE)
})
