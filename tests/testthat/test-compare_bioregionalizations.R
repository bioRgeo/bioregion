# Inputs -----------------------------------------------------------------------
bioregionalizations <- data.frame(matrix(nr = 4, nc = 4, 
                                         c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
                                         byrow = TRUE))

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  test_output <- compare_bioregionalizations(bioregionalizations)
  
  expect_equal(inherits(test_output, "list"), TRUE)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    compare_bioregionalizations(bioregionalizations, indices = 1),
    "indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    compare_bioregionalizations(bioregionalizations, indices = "zz"),
    "Please choose algorithm among the followings values:
    rand or jaccard",
    fixed = TRUE)
  
})
