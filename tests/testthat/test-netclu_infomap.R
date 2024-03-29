# Inputs -----------------------------------------------------------------------
install_binaries()

net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))


# Tests for valid outputs ------------------------------------------------------
#test_that("valid output", {
#  
#})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {

  expect_error(
    netclu_infomap(net, binpath = 1),
    "binpath must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, binpath = c("zz","zz")),
    "binpath must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, version = 1),
    "version must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, version = c("zz","zz")),
    "version must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite_version = 1),
    "bipartite_version must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite_version = c("zz","zz")),
    "bipartite_version must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
})
