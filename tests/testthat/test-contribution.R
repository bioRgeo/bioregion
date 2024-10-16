# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

dissim <- dissimilarity(comat, metric = "Simpson")
clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

net <- similarity(comat, metric = "Simpson")
com <- netclu_greedy(net)

clust2 <- nhclu_kmeans(dissim, n_clust = 3:4, index = "Simpson")

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  contrib <- contribution(cluster_object = clust1, comat = comat,
                          indices = "contribution")
  
  expect_equal(inherits(contrib, "data.frame"), TRUE)
  expect_equal(dim(contrib)[1], 75)
  expect_equal(dim(contrib)[2], 3)
  
  contrib2 <- contribution(cluster_object = com, comat = comat,
                           indices = "contribution")
  
  expect_equal(inherits(contrib2, "data.frame"), TRUE)
  expect_equal(dim(contrib2)[1], 50)
  expect_equal(dim(contrib2)[2], 3)
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    contribution("zz"),
    "This function is designed to work on bioregion.clusters objects and
         on a site x species matrix.",
    fixed = TRUE)
  
  expect_error(
    contribution(clust2, comat = comat, indices = "contribution"),
    "This function is designed to be applied on a single partition.Your cluster_object has multiple partitions (select only one).",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = comat, indices = "zz"),
    "Please choose algorithm among the followings values:
    contribution or Cz.",
    fixed = TRUE)
  
})
