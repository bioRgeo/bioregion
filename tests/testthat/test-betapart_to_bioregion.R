# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_warning(
    expect_error(
      betapart_to_bioregion(1),
      "betapart_result must be a valid object from the betapart package.",
      fixed = TRUE),
    "^'betapart_to_bioregion' is deprecated."
  )

  
})