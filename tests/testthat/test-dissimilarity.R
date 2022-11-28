# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  dissim <- dissimilarity(comat, metric = "all")
  
  expect_equal(ncol(dissim), 15L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    dissim <- dissimilarity(comat, metric = "zzz"),
    "One or several similarity metric(s) chosen is not available.
     Please chose among the followings:
         abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, Bray, Brayturn or
         Euclidean.",
    fixed = TRUE)
})

test_that("error messages with wrong inputs 2", {
  expect_error(
    dissim <- dissimilarity(comat, method = "zzz"),
    "The method is not available.
     Please chose among the followings:
         prodmat, loops",
    fixed = TRUE)
})