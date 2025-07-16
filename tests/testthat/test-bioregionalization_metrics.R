# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

d <- dist(comat)
simil <- similarity(comat, metric = "all")
dissim <- similarity_to_dissimilarity(simil)
comat_df <- mat_to_net(comat, weight = TRUE, remove_zeroes = TRUE)

comat2 <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                10, 50)
rownames(comat2) <- paste0("Site", 31:40)
colnames(comat2) <- paste0("Species", 51:100)

d2 <- dist(comat2)
simil2 <- similarity(comat2, metric = "all")
dissim2 <- similarity_to_dissimilarity(simil2)
comat_df2 <- mat_to_net(comat2, weight = TRUE, remove_zeroes = TRUE)

clu1 <- hclu_hierarclust(dissim, 
                         n_clust = 5,
                         index = "Simpson",
                         optimal_tree_method = "best",
                         verbose = FALSE)

clu2 <- hclu_hierarclust(dissim,
                         optimal_tree_method = "best",
                         n_clust = NULL,
                         cut_height = NULL,
                         verbose = FALSE)

clu3 <- netclu_louvain(simil)
clu3$clusters <- NULL

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  a <- bioregionalization_metrics(clu1,
                         dissimilarity = dissim,
                         net = comat_df,
                         site_col = "Node1",
                         species_col = "Node2",
                         eval_metric = c("tot_endemism",
                                         "avg_endemism",
                                         "pc_distance",
                                         "anosim"))
  expect_identical(class(a)[1], "bioregion.bioregionalization.metrics")
  expect_identical(class(a)[2], "list")

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    bioregionalization_metrics(clu2),
    "^No clusters have been generated")
  
  expect_error(
    bioregionalization_metrics(clu3),
    "^bioregionalization does not have the expected type of ")
  
  expect_error(
    bioregionalization_metrics(NULL),
    "^This function is designed to work on bioregion.clusters objects")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               eval_metric = c(1,2)),
    "eval_metric must be a character.")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               eval_metric = c("yy","zz")),
    "^One or several")
  
  expect_error(
    expect_warning(
      bioregionalization_metrics(clu1,
                                 eval_metric = "pc_distance"),
      "^No dissimilarity oject provided"),
    "^No evaluation metric")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = simil,
                               eval_metric = "all"),
    "^dissimilarity must be an object containing dissimilarity")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = c("z","z"),
                               eval_metric = "all"),
    "dissimilarity_index must be of length 1.")
  
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = TRUE,
                               eval_metric = "all"),
    "dissimilarity_index must be a character.")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "shalala",
                               eval_metric = "all"),
    "^dissimilarity_index does not exist")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = 1,
                               dissimilarity_index = NULL,
                               eval_metric = "all"),
    "^dissimilarity must be a bioregion.pairwise object or")

  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim2,
                               dissimilarity_index = "Sorensen",
                               eval_metric = "all"),
    "bioregionalization and dissimilarity have different number of sites.")
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                               site_col = "Weight", 
                               species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                               site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                               species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                              site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df, 
                               site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                               site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df,
                               species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(clu1,
                               dissimilarity = dissim,
                               dissimilarity_index = "Sorensen",
                               net = comat_df2,
                               eval_metric = "all"),
    "^Some elements of the cluster table")
  
})
