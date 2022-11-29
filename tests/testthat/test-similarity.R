# Preamble code ----------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

# Tests for valid outputs -----------------------------------------------------
test_that("number of columns in output", {
  simil <- similarity(comat, metric = c("abc", "ABC", "Simpson", "Brayturn"))
  
  expect_equal(nrow(simil), 10L)
  expect_equal(ncol(simil), 10L)
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    similarity(comat, metric = "zzz"),
    "One or several similarity metric(s) chosen is not available.
     Please chose among the followings:
         abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, Bray, Brayturn or
         Euclidean.",
    fixed = TRUE)
})

test_that("error messages with wrong inputs 2", {
  expect_error(
    similarity(comat, method = "zzz"),
    "The method is not available.
     Please chose among the followings: prodmat, loops",
    fixed = TRUE)
})