# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

simil <- similarity(comat, metric = "all")
dissim<- similarity_to_dissimilarity(simil)

tree <- hclu_hierarclust(dissim,
                         optimal_tree_method = "best",
                         n_clust = 5,
                         verbose = FALSE)

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
  
  expect_error(
    cut_tree(tree,
             show_hierarchy = 1),
    "show_hierarchy must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,
             show_hierarchy = c(TRUE, FALSE)),
    "show_hierarchy must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,
             verbose = 1),
    "verbose must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    cut_tree(tree,
             verbose = c(TRUE, FALSE)),
    "verbose must be of length 1.",
    fixed = TRUE)
  
})

# Tests for show_hierarchy argument --------------------------------------------
test_that("show_hierarchy argument works in cut_tree", {
  
  tree_no_cut <- hclu_hierarclust(dissim,
                                  optimal_tree_method = "best",
                                  verbose = FALSE)
  
  cut1 <- cut_tree(tree_no_cut,
                   n_clust = c(3, 5),
                   show_hierarchy = FALSE,
                   verbose = FALSE)
  expect_equal(cut1$args$show_hierarchy, FALSE)
  expect_equal(inherits(cut1, "bioregion.clusters"), TRUE)
  
  cut2 <- cut_tree(tree_no_cut,
                   n_clust = c(3, 5),
                   show_hierarchy = TRUE,
                   verbose = FALSE)
  expect_equal(cut2$args$show_hierarchy, TRUE)
  expect_equal(inherits(cut2, "bioregion.clusters"), TRUE)
  
  # Both should have same structure
  expect_equal(dim(cut1$clusters), dim(cut2$clusters))
  
  # Test that summary works with both show_hierarchy settings
  quietly(expect_no_error(summary(cut1)))
  quietly(expect_no_error(summary(cut2)))
  
  # Verify hierarchical status is properly set when multiple cuts
  expect_equal(cut1$inputs$hierarchical, TRUE)
  expect_equal(cut2$inputs$hierarchical, TRUE)
  
  # Verify cluster_info is present
  expect_true(!is.null(cut1$cluster_info))
  expect_true(!is.null(cut2$cluster_info))
  expect_equal(nrow(cut1$cluster_info), 2)
  expect_equal(nrow(cut2$cluster_info), 2)
  
  # Test with single cut
  cut3 <- cut_tree(tree_no_cut,
                   n_clust = 5,
                   show_hierarchy = FALSE,
                   verbose = FALSE)
  expect_equal(cut3$args$show_hierarchy, FALSE)
  quietly(expect_no_error(summary(cut3)))
  
  cut4 <- cut_tree(tree_no_cut,
                   n_clust = 5,
                   show_hierarchy = TRUE,
                   verbose = FALSE)
  expect_equal(cut4$args$show_hierarchy, TRUE)
  quietly(expect_no_error(summary(cut4)))
  
  # With cut_height instead of n_clust
  cut5 <- cut_tree(tree_no_cut,
                   cut_height = c(0.3, 0.5),
                   show_hierarchy = FALSE,
                   verbose = FALSE)
  expect_equal(cut5$args$show_hierarchy, FALSE)
  quietly(expect_no_error(summary(cut5)))
  
  cut6 <- cut_tree(tree_no_cut,
                   cut_height = c(0.3, 0.5),
                   show_hierarchy = TRUE,
                   verbose = FALSE)
  expect_equal(cut6$args$show_hierarchy, TRUE)
  quietly(expect_no_error(summary(cut6)))
  
})


