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

multi_clust <- nhclu_kmeans(dissim, n_clust = 3:4, index = "Simpson")

net_bip <- mat_to_net(comat, weight = TRUE)
clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  rho <- contribution(cluster_object = clust1, comat = comat,
                      indices = "rho")
  
  expect_equal(inherits(rho, "data.frame"), TRUE)
  expect_equal(dim(rho)[1], 75)
  expect_equal(dim(rho)[2], 3)
  
  rho2 <- contribution(cluster_object = com, comat = comat,
                       indices = "rho")
  
  expect_equal(inherits(rho2, "data.frame"), TRUE)
  expect_equal(dim(rho2)[1], 50)
  expect_equal(dim(rho2)[2], 3)
  
  suppressWarnings({
    rho3 <- contribution(cluster_object = clust_bip, comat = comat,
                         bipartite_link = net_bip,
                         indices = c("rho", "Cz"))
  })
  expect_equal(inherits(rho3, "list"), TRUE)
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    contribution("zz"),
    "This function is designed to work on bioregion.clusters objects and
         on a site x species matrix.",
    fixed = TRUE)
  
  expect_error(
    contribution(multi_clust, comat = comat, indices = "rho"),
    "This function is designed to be applied on a single partition.Your cluster_object has multiple partitions (select only one).",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = comat, indices = "zz"),
    "Please choose algorithm among the followings values:
    rho or Cz.",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = comat, bipartite_link = "zz", indices = "Cz"),
    "Cz metrics can only be computed for a bipartite partition (where
         both sites and species are assigned to a bioregion.",
    fixed = TRUE)
  
  expect_error(
    contribution(com, comat = comat, bipartite_link = NULL, indices = "Cz"),
    "bipartite_link is needed to compute Cz indices.",
    fixed = TRUE)
  
})
