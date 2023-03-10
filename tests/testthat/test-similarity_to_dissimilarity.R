# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001),
                5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

simil <- similarity(comat, metric = "all")

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  dissimilarity <- similarity_to_dissimilarity(simil)
  
  expect_equal(nrow(simil), nrow(dissimilarity))
  expect_equal(ncol(simil), ncol(dissimilarity))
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    similarity_to_dissimilarity(NULL),
    "similarity should be a bioregion object created by similarity() or
         dissimilarity_to_similarity()",
    fixed = TRUE)
})
