# Inputs -----------------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

clu1 <- netclu_louvain(net, lang = "igraph", bipartite = TRUE)
clu2 <- hclu_hierarclust(net)
clu3 <- netclu_louvain(net, lang = "igraph")
clu4 <- netclu_louvain(net, 
                       lang = "igraph", 
                       bipartite = TRUE, 
                       return_node_type = "sites")

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {
  
  sub <- subset_node(clu1, node_type = "sites")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "sites")
  
  sub <- subset_node(clu1, node_type = "species")
  expect_equal(inherits(sub, "bioregion.clusters"), TRUE)
  expect_equal(sub$args$return_node_type, "species")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("indalid inputs", {
  
  expect_error(
    subset_node(clu1, node_type = 1),
    "node_type must be a character.", 
    fixed = TRUE)
  
  expect_error(
    subset_node(clu1, node_type = "1"),
    "Please choose node_type among the followings values:
sites and species", 
    fixed = TRUE)
  
  expect_error(
    subset_node("1", node_type = "sites"),
    "clusters must be a bioregion.clusters object.", 
    fixed = TRUE)
  
  expect_error(
    subset_node(clu2, node_type = "sites"),
    "clusters must be an output of a 'netclu_' function.", 
    fixed = TRUE)
  
  expect_error(
    subset_node(clu3, node_type = "sites"),
    "clusters must be based on a bipartite network.", 
    fixed = TRUE)
  
  expect_error(
    subset_node(clu4, node_type = "sites"),
    "clusters must contain both types of node.", 
    fixed = TRUE)
  
})