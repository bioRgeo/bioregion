# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

simil <- similarity(comat, metric = "all")
dissim<- similarity_to_dissimilarity(simil)

tree <- hclu_hierarclust(dissim,
                         optimal_tree_method = "best",
                         n_clust = 5)

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  expect_identical(class(tree)[1], "bioregion.clusters")
  expect_identical(class(tree)[2], "list")

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    cut_tree(tree, 
             n_clust = 1.5),
    "^n_clust must an integer")
  
  expect_error(
    cut_tree(tree, 
             n_clust = "a"),
    "^n_clust must be one of those")
  
  expect_error(
    cut_tree(tree, 
             n_clust = 1,
             cut_height = 1),
    "^Please provide either n_clust or cut_height,")
  
  expect_error(
    cut_tree(tree, 
             n_clust = NULL, 
             cut_height = "zz"),
    "cut_height must be numeric.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             n_clust = NULL, 
             cut_height = -1),
    "cut_height must be composed of values higher than 0.",
    fixed = TRUE)
  
  
  expect_error(
    cut_tree(tree, 
             find_h = 1),
    "find_h must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             find_h = c("zz","zz")),
    "find_h must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,  
             h_min =  c("zz","zz")),
    "h_min must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             h_min = "zz"),
    "h_min must be numeric.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,  
             h_min = -1),
    "h_min must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             h_max =  c("zz","zz")),
    "h_max must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,  
             h_max = "zz"),
    "h_max must be numeric.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             h_max = -1),
    "h_max must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             h_min = 1,
             h_max = 0),
    "h_min must be inferior to h_max.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = 1),
    "dynamic_tree_cut must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = c("zz","zz")),
    "dynamic_tree_cut must be of length 1.",
    fixed = TRUE)
  
  expect_message(
    cut_tree(tree, 
             n_clust = 1,
             dynamic_tree_cut = TRUE),
    "^The dynamic tree cut method was requested,")
  
  expect_message(
    cut_tree(tree, 
             cut_height = 1,
             dynamic_tree_cut = TRUE),
    "^The dynamic tree cut method was requested,")
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = TRUE,
             dynamic_method = TRUE),
    "dynamic_method must be a character.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = TRUE,
             dynamic_method = c(1, "1")),
    "dynamic_method must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = TRUE,
             dynamic_method = "a"),
    "^Please choose dynamic_method from the following")
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = TRUE,
             dynamic_method = "hybrid",
             dissimilarity = 1),
    "^dissimilarity is not a")
  
  expect_error(
    cut_tree(tree, 
             dynamic_tree_cut = TRUE,
             dynamic_method = "hybrid",
             dissimilarity = data.frame(ID = 1)),
    "^dissimilarity is not a")
  
  
  expect_error(
    cut_tree(tree = 1),
    "^This function is designed to work either on outputs")
  
  
})


