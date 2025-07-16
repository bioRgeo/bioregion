# Inputs -----------------------------------------------------------------------
clu1 <- data.frame(matrix(nr = 4, nc = 4, 
                          c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
                          byrow = TRUE))

dissim <- dissimilarity(fishmat, metric = "all")
clu2 <- hclu_hierarclust(dissim,
                         optimal_tree_method = "best",
                         n_clust = NULL,
                         cut_height = NULL,
                         verbose = FALSE)

sim <- similarity(fishmat, metric = "all")
clu3 <- netclu_louvain(sim)
clu3$clusters <- NULL

df <- clu1
df[1,1] <- NA

ldf <- list(df, 1)

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  comp <- compare_bioregionalizations(clu1)
  expect_equal(inherits(comp, "list"), TRUE)
  
  comp <- compare_bioregionalizations(clu1,
                                      cor_frequency = TRUE)
  expect_equal(as.numeric(comp$bioregionalization_freq_cor[1]) > 0, TRUE)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    compare_bioregionalizations(clu2),
    "^No clusters have been generated")
  
  expect_error(
    compare_bioregionalizations(clu3),
    "^bioregionalizations does not have the expected type of ")
  
  expect_error(
    compare_bioregionalizations(df),
    "^NA")
  
  expect_error(
    compare_bioregionalizations(ldf),
    "^All elements in bioregionalizations should be")
  
  expect_error(
    compare_bioregionalizations(1),
    "^This function is designed to work either on")
  
  expect_error(
    compare_bioregionalizations(data.frame(ID = 1:10)),
    "^This function is designed to be applied on multiple ")
  
  expect_error(
    compare_bioregionalizations(clu1, 
                                indices = 1),
    "indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    compare_bioregionalizations(clu1, 
                                indices = "zz"),
    "^Please choose indices from the following")
  
  expect_error(
    compare_bioregionalizations(clu1, 
                                indices = c("rand","zz")),
    "^Please choose indices from the following")
  
})
