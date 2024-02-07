# Preamble code ----------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

net2 <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c("a", "b", "a", "c", "d", "b", "d")
)

# Tests for valid outputs -----------------------------------------------------
test_that("dimension of output", {
  mat <- net_to_mat(net, weight = TRUE)
  
  expect_equal(nrow(mat), 3L)
  expect_equal(ncol(mat), 4L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    net_to_mat(net, weight = "zz"),
    "weight must be a boolean.", fixed = TRUE)
  
  expect_error(
    net_to_mat(net, squared = "zz"),
    "squared must be a boolean", fixed = TRUE)
  
  expect_error(
    net_to_mat(net, symmetrical = "zz"),
    "symmetrical must be a boolean", fixed = TRUE)
  
  expect_error(
    net_to_mat(net, missing_value = "zz"),
    "missing_value must be numeric.", fixed = TRUE)
  
  expect_error(
    net_to_mat(net2, weight = TRUE),
    "The third column of net must be numeric", fixed = TRUE)
  
  
})
