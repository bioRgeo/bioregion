# # Preamble code ----------------------------------------------------------------
# comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#                 20, 25)
# rownames(comat) <- paste0("Site", 1:20)
# colnames(comat) <- paste0("Species", 1:25)
# 
# simil <- similarity(comat, metric = "all")
# dissimilarity <- similarity_to_dissimilarity(simil)
# 
# comat_df <- mat_to_net(comat, weight = TRUE, remove_zeroes = TRUE)
# 
# tree1 <- hclu_hierarclust(dissimilarity, n_clust = 2:20,
#                           index = "Simpson")
# 
# # Tests for valid outputs -----------------------------------------------------
# test_that("class list and bioregion.clusters", {
#   a <- partition_metrics(tree1,
#                          dissimilarity = dissimilarity,
#                          net = comat_df,
#                          site_col = "Node1",
#                          species_col = "Node2",
#                          eval_metric = c("tot_endemism",
#                                          "avg_endemism",
#                                          "pc_distance",
#                                          "anosim"))
#   
#   
#   expect_identical(class(a)[1], "bioregion.partition.metrics")
#   expect_identical(class(a)[2], "list")
#   
# })
# 
# # Tests for invalid inputs ----------------------------------------------------
# test_that("error messages with wrong inputs", {
#   expect_error(
#     partition_metrics(NULL,
#                       dissimilarity = dissimilarity,
#                       net = comat_df,
#                       site_col = "Node1",
#                       species_col = "Node2",
#                       eval_metric = c("tot_endemism",
#                                       "avg_endemism",
#                                       "pc_distance",
#                                       "anosim")),
#     "This function is designed to work either on bioregion.clusters objects
#          (outputs from clustering functions)", fixed = TRUE)
# })
