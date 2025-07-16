# Inputs -----------------------------------------------------------------------
install_binaries()

net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

clu1 <- netclu_louvain(net, 
                       lang = "igraph", 
                       bipartite = TRUE)
clu2 <- hclu_hierarclust(net,
                         optimal_tree_method = "best",
                         verbose = FALSE)
clu3 <- netclu_louvain(net, 
                       lang = "igraph")
clu4 <- netclu_louvain(net,
                       lang = "igraph",
                       bipartite = TRUE,
                       return_node_type = "site")
clu5 <- netclu_beckett(net)
clu6 <- netclu_infomap(net,
                       bipartite = TRUE)

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {

  sub <- site_species_subset(clu1, 
                             node_type = "site")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "site")

  sub <- site_species_subset(clu1, 
                             node_type = "species")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "species")
  
  sub <- site_species_subset(clu5, 
                             node_type = "species")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "species")
  
  sub <- site_species_subset(clu6, 
                             node_type = "site")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "site")

})

# Tests for invalid inputs -----------------------------------------------------
test_that("indalid inputs", {

  expect_error(
    site_species_subset(clu1, 
                        node_type = 1),
    "node_type must be a character.",
    fixed = TRUE)

  expect_error(
    site_species_subset(clu1, 
                        node_type = c(1,1)),
    "node_type must be of length 1.",
    fixed = TRUE)

  expect_error(
    site_species_subset(clu1, 
                        node_type = "1"),
    "^Please choose node_type from the following")

  expect_error(
    site_species_subset("1", 
                        node_type = "site"),
    "clusters must be a bioregion.clusters object.",
    fixed = TRUE)

  expect_error(
    site_species_subset(clu2, 
                        node_type = "site"),
    "clusters must be an output of a 'netclu_' function.",
    fixed = TRUE)

  expect_error(
    site_species_subset(clu3, 
                        node_type = "site"),
    "clusters must be based on a bipartite network.",
    fixed = TRUE)

  expect_error(
    site_species_subset(clu4, 
                        node_type = "site"),
    "clusters must contain both types of node.",
    fixed = TRUE)

})