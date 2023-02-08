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
  com <- netclu_leiden(net)
  
  expect_equal(nrow(com$clusters), 5L)
  expect_equal(ncol(com$clusters), 2L)
  expect_equal(length(unique(com$clusters)), 2L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(netclu_leiden("zz"), "net must be a data.frame.",
               fixed = TRUE)
  
  expect_error(netclu_leiden(net, weight = "zz"),
               "weight must be a boolean.", fixed = TRUE)
  
  expect_error(netclu_leiden(net, index = "zz"),
               "index is a character, it should be a column name (and not the
                    first or second column).", fixed = TRUE)

  expect_error(netclu_leiden(net, bipartite = "zz"),
               "bipartite must be a boolean", fixed = TRUE)
  
  expect_error(netclu_leiden(net, algorithm_in_output = "zz"),
               "algorithm_in_output must be a boolean", fixed = TRUE)
  
})
