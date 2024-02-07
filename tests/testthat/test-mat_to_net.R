# Preamble code ----------------------------------------------------------------
mat <- matrix(sample(1000, 50), 5, 10)
rownames(mat) <- paste0("Site", 1:5)
colnames(mat) <- paste0("Species", 1:10)

# Tests for valid outputs -----------------------------------------------------
test_that("dimension of output", {
  net <- mat_to_net(mat, weight = TRUE)
  
  expect_equal(nrow(net), 50L)
  expect_equal(ncol(net), 3L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    mat_to_net(mat, weight = "zz"),
    "weight must be a boolean", fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, remove_zeroes = "zz"),
    "remove_zeroes must be a boolean", fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, include_diag = "zz"),
    "include_diag must be a boolean", fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, include_lower = "zz"),
    "include_lower must be a boolean", fixed = TRUE)
  
  expect_message(
    mat_to_net(mat, include_diag = FALSE),
    "include_diag is only used with squared matrix.", fixed = TRUE)
  
  expect_message(
    mat_to_net(mat, include_lower = FALSE),
    "include_lower is only used with squared matrix.", fixed = TRUE)

})
