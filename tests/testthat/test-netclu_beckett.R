# Preamble code ----------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  clust1 <- netclu_beckett(net)
  
  expect_equal(clust1$name, "netclu_beckett")
  expect_equal(ncol(clust1$clusters), 2L)
  expect_equal(length(unique(clust1$clusters)), 2L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    netclu_beckett(net = "zz"),
    "net must be a data.frame.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, weight = "zz"),
    "weight must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, index = "zz"),
    "index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, return_node_type = "zz"),
    "Please choose return_node_type among the followings values:
both, sites and species",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, forceLPA = "zz"),
    "forceLPA must be a boolean",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net = net, algorithm_in_output = "zz"),
    "algorithm_in_output must be a boolean",
    fixed = TRUE)
})
