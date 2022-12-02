# Preamble code ----------------------------------------------------------------
comat <- matrix(c(1, 1, 1, 0, 0,
                  1, 1, 1, 0, 0,
                  1, 0, 0, 0, 0,
                  0, 0, 0, 1, 1,
                  0, 1, 1, 1, 0),
                byrow = TRUE, nrow = 5, ncol = 5,
                dimnames = list(paste0("Site", 1:5),
                                paste0("Species", 1:5)))

net <- similarity(comat, metric = "Simpson")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  com <- netclu_infomap(net)
  
  expect_equal(nrow(com$clusters), 5L)
  expect_equal(ncol(com$clusters), 2L)
  expect_equal(length(unique(com$clusters)), 2L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(netclu_infomap("zz"), "net must be a data.frame.",
               fixed = TRUE)
  
  expect_error(netclu_infomap(net, weight = "zz"),
               "weight must be a boolean.", fixed = TRUE)
  
  expect_error(netclu_infomap(net, index = "zz"),
               "index is a character, it should be a column name (and not the
                    first or second column).", fixed = TRUE)
  
  expect_error(netclu_infomap(net, nbmod = "zz"),
               "nbmod must be numeric.", fixed = TRUE)
  
  expect_error(netclu_infomap(net, markovtime = "zz"),
               "markovtime must be numeric.", fixed = TRUE)
  
  expect_error(netclu_infomap(net, seed = "zz"),
               "seed must be numeric.", fixed = TRUE)
  
  expect_error(netclu_infomap(net, numtrials = "zz"),
               "numtrials must be numeric.", fixed = TRUE)
  
  expect_error(netclu_infomap(net, twolevel = "zz"),
               "twolevel must be a boolean", fixed = TRUE)
  
  expect_error(netclu_infomap(net, show_hierarchy = "zz"),
               "show_hierarchy must be a boolean", fixed = TRUE)
  
  expect_error(netclu_infomap(net, directed = "zz"),
               "directed must be a boolean", fixed = TRUE)
  
  expect_error(netclu_infomap(net, bipartite_version = "zz"),
               "bipartite_version must be a boolean", fixed = TRUE)
  
  expect_error(netclu_infomap(net, bipartite = "zz"),
               "bipartite must be a boolean", fixed = TRUE)
  
  # expect_error(netclu_infomap(net, site_col = "zz"),
  #              "site_col must be a boolean", fixed = TRUE)
  
  # expect_error(netclu_infomap(net, species_col = "zz"),
  #              "site_col must be a boolean", fixed = TRUE)
  
  # expect_error(netclu_infomap(net, return_node_type = "zz"),
  #              "return_node_type must be a boolean", fixed = TRUE)
  
  expect_error(netclu_infomap(net, delete_temp = "zz"),
               "delete_temp must be a boolean", fixed = TRUE)
  
  # expect_error(netclu_infomap(net, path_temp = "zz"),
  #              "delete_temp must be a boolean", fixed = TRUE)
  
  # expect_error(netclu_infomap(net, binpath = "zz"),
  #              "binpath must be a boolean", fixed = TRUE)
})
