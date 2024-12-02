# Inputs -----------------------------------------------------------------------
comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
                         prob = 1/1:1001), 10, 12)
rownames(comat_1) <- paste0("Site", 1:10)
colnames(comat_1) <- paste0("Species", 1:12)
comat_1 <- cbind(comat_1,
                 matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 1:10),
                                        paste0("Species", 13:20))))

comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
                         prob = 1/1:1001), 10, 12)
rownames(comat_2) <- paste0("Site", 11:20)
colnames(comat_2) <- paste0("Species", 9:20)
comat_2 <- cbind(matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 11:20),
                                        paste0("Species", 1:8))),
                 comat_2)

comat <- rbind(comat_1, comat_2)

dissim <- dissimilarity(comat, metric = "Simpson")
clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

net <- similarity(comat, metric = "Simpson")
com <- netclu_greedy(net)

multi_clust <- nhclu_kmeans(dissim, n_clust = 3:4, index = "Simpson")

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  test_output <- bioregion_metrics(cluster_object = clust1, comat = comat)
  
  expect_equal(inherits(test_output, "data.frame"), TRUE)
  expect_equal(dim(test_output)[1], 3)
  expect_equal(dim(test_output)[2], 5)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    site_species_metrics("zz"),
    "This function is designed to work on bioregion.clusters objects and
         on a site x species matrix.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(multi_clust, comat = comat),
    "This function is designed to be applied on a single partition.Your cluster_object has multiple partitions (select only one).",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(com, comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)
  
})
