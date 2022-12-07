# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

net <- similarity(comat, metric = "Simpson")

net_bip <- mat_to_net(comat, weight = TRUE)

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  clust1 <- netclu_labelprop(net)
  
  expect_equal(clust1$name, "netclu_labelprop")
  expect_equal(ncol(clust1$clusters), 2L)
  expect_equal(nrow(clust1$clusters), 5L)
  expect_equal(length(unique(clust1$clusters)), 2L)
  
  clust2 <- netclu_labelprop(net_bip, bipartite = TRUE)
  
  expect_equal(clust2$name, "netclu_labelprop")
  expect_equal(ncol(clust2$clusters), 2L)
  expect_equal(nrow(clust2$clusters), 15L)
  expect_equal(length(unique(clust2$clusters)), 2L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    netclu_labelprop(net = "zz"),
    "net must be a data.frame.",
    fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net, weight = "zz"),
    "weight must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net, index = "zz"),
    "index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net, bipartite = "zz"),
    "bipartite must be a boolean",
    fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net_bip, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net_bip, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net_bip, bipartite = TRUE, return_node_type = "zz"),
    "Please choose return_node_type among the followings values:
both, sites and species",
    fixed = TRUE)
  
  expect_error(
    netclu_labelprop(net = net, algorithm_in_output = "zz"),
    "algorithm_in_output must be a boolean",
    fixed = TRUE)
})
