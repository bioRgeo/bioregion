# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001),
                5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

dissimil <- dissimilarity(comat, metric = "all", formula = "a + b")
simil <- dissimilarity_to_similarity(dissimil, include_formula = TRUE)
simil2 <- dissimilarity_to_similarity(dissimil, include_formula = FALSE)

dissimil2 <- dissimil
attr(dissimil2, "type") <- NULL

dissimil3 <- similarity(comat, metric = "all")

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  expect_equal(simil$Jaccard, 1-dissimil$Jaccard)
  expect_equal(simil2$Jaccard, 1-dissimil$Jaccard)
  expect_equal(simil$Jaccardturn, 1-dissimil$Jaccardturn)
  expect_equal(simil$Sorensen, 1-dissimil$Sorensen)
  expect_equal(simil$Simpson, 1-dissimil$Simpson)
  expect_equal(simil$Bray, 1-dissimil$Bray)
  expect_equal(simil$Brayturn, 1-dissimil$Brayturn)
  expect_equal(simil$Euclidean, 1/(1+dissimil$Euclidean))
  expect_equal(simil[,dim(simil)[2]], 1-dissimil[,dim(simil)[2]])
  expect_equal(simil2[,dim(simil)[2]], dissimil[,dim(simil)[2]])
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    dissimilarity_to_similarity("2", include_formula = TRUE),
    "dissimilarity should be a bioregion.pairwise.metric object created by 
dissimilarity() or similarity_to_dissimilarity().",
    fixed = TRUE)
  
  expect_error(
    dissimilarity_to_similarity(dissimil2, include_formula = TRUE),
    "dissimilarity is a bioregion.pairwise.metric object but it has not 
been possible to identify the object's type (similarity or dissimilarity) 
probably because the bioregion.pairwise.metric object has been altered.",
    fixed = TRUE)
  
  expect_error(
    dissimilarity_to_similarity(dissimil3, include_formula = TRUE),
    "dissimilarity is already composed of 
similarity metrics. If you want to convert it to dissimilarity, use 
similarity_to_dissimilarity().",
    fixed = TRUE)
  
  expect_error(
    dissimilarity_to_similarity(dissimil, include_formula = "zz"),
    "include_formula must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    dissimilarity_to_similarity(dissimil, include_formula = c(TRUE,TRUE)),
    "include_formula must be of length 1.",
    fixed = TRUE)
  
})