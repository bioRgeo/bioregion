# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

beta_div <- betapart::beta.pair.abund(comat)

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  dissimil <- betapart_to_bioregion(beta_div)
  expect_equal(inherits(dissimil, "bioregion.pairwise.metric"), TRUE)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  test <- 1
  expect_error(
    betapart_to_bioregion(test),
    "betapart_result must be a valid object from the betapart package.",
    fixed = TRUE)
  
})