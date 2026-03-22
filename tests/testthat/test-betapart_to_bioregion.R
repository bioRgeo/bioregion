# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    betapart_to_bioregion(1),
    "^This function is deprecated, please use as_bioregion_pairwise instead.")

})