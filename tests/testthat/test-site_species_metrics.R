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

clust_h <- hclu_hierarclust(dissim,
                            optimal_tree_method = "best",
                            n_clust = NULL,
                            cut_height = NULL)

simil <- dissimilarity_to_similarity(dissim)
clust_louv <- netclu_louvain(simil)
clust_louv$clusters <- NULL

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  rho <- site_species_metrics(cluster_object = clust1, 
                              comat = comat,
                              indices = "rho")
  expect_equal(inherits(rho, "data.frame"), TRUE)
  expect_equal(dim(rho)[1], 75)
  expect_equal(dim(rho)[2], 3)
  
  rho2 <- site_species_metrics(cluster_object = com, 
                               comat = comat,
                               indices = "rho")
  expect_equal(inherits(rho2, "data.frame"), TRUE)
  expect_equal(dim(rho2)[1], 50)
  expect_equal(dim(rho2)[2], 3)
  
  suppressWarnings({
    rho3 <- site_species_metrics(cluster_object = clust_bip, 
                                 comat = comat,
                                 net = net_bip,
                                 indices = c("rho", "Cz"))
  })
  expect_equal(inherits(rho3, "list"), TRUE)
  
  #aff <- site_species_metrics(cluster_object = clust_bip, 
  #                               comat = comat,
  #                               indices = "affinity")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    site_species_metrics("zz"),
    "^This function is designed to work on bioregion.clusters objects ")
  
  expect_error(
    site_species_metrics(multi_clust, 
                         comat = comat, 
                         indices = "rho"),
    "^This function is designed to be applied on a single partition.")
  
  expect_error(
    site_species_metrics(clust_h),
    "^No clusters have been generated for your hierarchical")
  
  expect_error(
    site_species_metrics(clust_louv),
    "^cluster_object does not have the expected type of")
  
  expect_error(
    site_species_metrics(com, 
                         comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat,
                         indices = c(1,2)),
    "indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         indices = "zz"),
    "^Please choose indices from the following")
  
  expect_error(
    site_species_metrics(com, comat = comat, 
                         net = "zz",
                         indices = "Cz"),
    "^Cz metrics can only be computed for a bipartite partition")
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = NULL,
                         indices = "Cz"),
    "net is needed to compute Cz indices.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = "zz"),
    "net should be a data.frame with at least two columns,
           corresponding to the sites and species. By default, sites are
           considered to be in the first column, and species in the second.
           This can be changed with the arguments 'site_col' and
           'species_col'.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = net_bip,
                         site_col = "zz"),
    "site_col must be numeric.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = net_bip,
                         species_col = "zz"),
    "species_col must be numeric.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(clust_bip, 
                         comat = comat, 
                         net = net_bip,
                         indices = "Cz",
                         site_col = 30),
    "The site column ('site_col') is incorrect.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(clust_bip, 
                         comat = comat, 
                         net = net_bip,
                         indices = "Cz",
                         species_col = 30),
    "The species column ('species_col') is incorrect.",
    fixed = TRUE)
  
})
