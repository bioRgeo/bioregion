# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001),
                5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

simil <- similarity(comat, metric = "all", formula = "a + b")
dissimil <- similarity_to_dissimilarity(simil, include_formula = TRUE)
dissimil2 <- similarity_to_dissimilarity(simil, include_formula = FALSE)

simil2 <- simil
attr(simil2, "type") <- NULL

simil3 <- dissimilarity(comat, metric = "all")

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  expect_equal(dissimil$Jaccard, 1-simil$Jaccard)
  expect_equal(dissimil2$Jaccard, 1-simil$Jaccard)
  expect_equal(dissimil$Jaccardturn, 1-simil$Jaccardturn)
  expect_equal(dissimil$Sorensen, 1-simil$Sorensen)
  expect_equal(dissimil$Simpson, 1-simil$Simpson)
  expect_equal(dissimil$Bray, 1-simil$Bray)
  expect_equal(dissimil$Brayturn, 1-simil$Brayturn)
  expect_equal(dissimil$Euclidean, (1/simil$Euclidean)-1)
  expect_equal(dissimil[,dim(simil)[2]], 1-simil[,dim(simil)[2]])
  expect_equal(dissimil2[,dim(simil)[2]], simil[,dim(simil)[2]])
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    similarity_to_dissimilarity("2", 
                                include_formula = TRUE),
    "^similarity should be a bioregion.pairwise object created by")
  
  expect_error(
    similarity_to_dissimilarity(simil2, 
                                include_formula = TRUE),
    "^similarity is a bioregion.pairwise object but it has not")
  
  expect_error(
    similarity_to_dissimilarity(simil3, 
                                include_formula = TRUE),
    "^similarity is already composed of")
  
  expect_error(
    similarity_to_dissimilarity(simil, include_formula = "zz"),
    "include_formula must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    similarity_to_dissimilarity(simil, include_formula = c(TRUE,TRUE)),
    "include_formula must be of length 1.",
    fixed = TRUE)
  
})

