# # Preamble code ----------------------------------------------------------------
# comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#                 20, 25)
# rownames(comat) <- paste0("Site",1:20)
# colnames(comat) <- paste0("Species",1:25)
# 
# dissim <- dissimilarity(comat, metric = "all")
# 
# # Tests for valid outputs -----------------------------------------------------
# test_that("number of columns in output", {
#   tree1 <- hclu_hierarclust(dissim, n_clust = 5)
#   
#   expect_equal(nrow(tree1$clusters), 20L)
#   expect_equal(ncol(tree1$clusters), 2L)
#   expect_equal(length(unique(tree1$clusters$K_5)), 5L)
# })
# 
# # Tests for invalid inputs ----------------------------------------------------
# test_that("error messages with wrong inputs", {
#   expect_error(
#     hclu_hierarclust(dissimilarity = "zz"),
#     "dissimilarity is not a bioregion.pairwise.metric object, a
#            dissimilarity matrix (class dist) or a data.frame with at least 3
#            columns (site1, site2, and your dissimilarity index)",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = "5"),
#     "n_clust must be one of those:
#         * an integer determining the number of clusters
#         * a vector of integers determining the numbers of clusters for each cut
#         * the output from partition_metrics()",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, method = "zz"),
#     "method is a character string indicating what hierarchical
#          classification method to use. See help for available options",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, randomize = "zz"),
#     "randomize must be a Boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs = "zz"),
#     "n_runs must be a positive integer.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, keep_trials = "zz"),
#     "keep_trials must be a Boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, optimal_tree_method = "zz"),
#     "optimal_tree_method must be a character string. Only available
#          option at the moment is best.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = "zz"),
#     "n_clust must be one of those:
#         * an integer determining the number of clusters
#         * a vector of integers determining the numbers of clusters for each cut
#         * the output from partition_metrics()",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, cut_height = "zz"),
#     "cut_height must be a positive integer.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, find_h = "zz"),
#     "find_h must be a Boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_max = "zz"),
#     "h_max must be a positive integer.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_min = "zz"),
#     "h_min must be a positive integer.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_min = 1, h_max = 0),
#     "h_min must be inferior to h_max.",
#     fixed = TRUE)
# })
