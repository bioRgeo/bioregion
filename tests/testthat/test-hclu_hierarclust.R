# # Inputs -----------------------------------------------------------------------
# dissim <- dissimilarity(fishmat, metric = "all")
# dissim2 <- dissim[,1:2]
# 
# df <- data.frame(ID1=dissim$Site1, ID2=dissim$Site2, W=dissim$Eucl)
# 
# d <- dist(fishmat)
# mat <- fishmat
# rownames(mat) <- NULL
# colnames(mat) <- NULL
# d2 <- dist(mat)
# d3 <- d
# d3[1] <- "1"
# d4 <- d
# d4[1] <- NA
# 
# 
# uni <- data.frame(
#   Site1 = c("c", "b", "a"),
#   Site2 = c("a", "c", "b"),
#   Weight = c(10, 100, 1))
# 
# unina1 <- uni
# unina1[1,1] <- NA
# 
# unina2 <- uni
# unina2[1,3] <- NA
# 
# unichar <- uni
# unichar$Weight <- unichar$Site1
# 
# uni2 <- data.frame(
#   Site1 = c("a", "c", "a"),
#   Site2 = c("a", "a", "b"),
#   Weight = c(10, 100, 1))
# 
# uni3 <- data.frame(
#   Site1 = c("c", "b", "a"),
#   Site2 = c("a", "a", "b"),
#   Weight = c(10, 100, 1))
# 
# uni4 <- data.frame(
#   Site1 = c("c", "a", "a"),
#   Site2 = c("a", "b", "b"),
#   Weight = c(10, 100, 1))
# 
# # Tests for valid outputs ------------------------------------------------------
# test_that("valid output", {
#   
#   clust <- hclu_hierarclust(dissim,
#                        index = "Simpson",
#                        method = "average",
#                        randomize = TRUE,
#                        n_runs = 30,
#                        keep_trials = FALSE,
#                        optimal_tree_method = "best",
#                        n_clust = c(1,2,3),
#                        cut_height = NULL,
#                        find_h = TRUE,
#                        h_max = 1,
#                        h_min = 0)
#   expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
#   expect_equal(clust$name, "hclu_hierarclust")
#   expect_equal(clust$args$index, "Simpson")
#   expect_equal(clust$args$method, "average")
#   expect_equal(clust$args$randomize, TRUE)
#   expect_equal(clust$args$n_runs, 30)
#   expect_equal(clust$args$keep_trials, FALSE)
#   expect_equal(clust$args$optimal_tree_method, "best")
#   expect_equal(clust$args$n_clust, c(1,2,3))
#   expect_equal(clust$args$cut_height, NULL)
#   expect_equal(clust$args$find_h, TRUE)
#   expect_equal(clust$args$h_max, 1)
#   expect_equal(clust$args$h_min, 0)
#   expect_equal(clust$inputs$bipartite, FALSE)
#   expect_equal(clust$inputs$weight, TRUE)
#   expect_equal(clust$inputs$pairwise, TRUE)
#   expect_equal(clust$inputs$pairwise_metric, "Simpson")
#   expect_equal(clust$inputs$dissimilarity, TRUE)
#   expect_equal(clust$inputs$nb_sites, 338)
#   expect_equal(clust$inputs$hierarchical, TRUE)
#   expect_equal(dim(clust$clusters)[2], 4)
#   
#   clust <- hclu_hierarclust(dissim,
#                        index = "Euclidean",
#                        n_clust = c(5,10),
#                        find_h = FALSE)
#   expect_equal(colnames(clust$clusters)[2], "K_5")
#   expect_equal(colnames(clust$clusters)[3], "K_10")
#   
#   clust <- hclu_hierarclust(dissim,
#                        index = "Euclidean",
#                        n_clust = c(10,5),
#                        find_h = FALSE)
#   expect_equal(colnames(clust$clusters)[2], "K_5")
#   expect_equal(colnames(clust$clusters)[3], "K_10")
#   
# })
# 
# # Tests for invalid inputs -----------------------------------------------------
# test_that("invalid inputs", {
#   
#   expect_error(
#     hclu_hierarclust(dissimilarity = "zz"),
#     "dissimilarity is not a bioregion.pairwise.metric object, 
# a dissimilarity matrix (class dist) or 
# a data.frame with at least 3 columns (site1, site2 and your dissimilarity index).",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim2),
#     "dissimilarity must be a data.frame with at least three columns.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(uni[,-3]),
#     "dissimilarity must be a data.frame with at least three columns.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(unina1),
#     "NA(s) detected in dissimilarity.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(unina2),
#     "NA(s) detected in dissimilarity.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(uni4),
#     "The first two columns of dissimilarity contain duplicated pairs of sites!",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(uni3),
#     "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(uni2),
#     "dissimilarity contains rows with the same site on both columns!",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(unichar),
#     "The weight column must be numeric.",
#     fixed = TRUE)  
#   
#   expect_message(
#     hclu_hierarclust(d2),
#     "No labels detected, they have been assigned automatically.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(d3),
#     "dissimilarity must be numeric.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(d4),
#     "NA(s) detected in dissimilarity.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(dissim, index = c("zz",1)),
#     "index must be of length 1.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, index = "zz"),
#     "If index is a character, it should be a column name (and not the
#                     first or second column).",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, index = "Site1"),
#     "If index is a character, it should be a column name (and not the
#                     first or second column).",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, index = 0.1),
#     "If index is numeric, it should be an integer.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, index = 2),
#     "index should be stricltly higher than 2.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(uni, index = 4),
#     "index should be lower or equal to 3.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, method = c(1,"z")),
#     "method must be of length 1.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, method = 1),
#     "method must be a character.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, method = "zz"),
#     "Please choose method among the followings values:
# ward.D, ward.D2, single, complete, average, mcquitty, median or centroid",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, randomize = 1),
#     "randomize must be a boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, randomize = c("zz","zz")),
#     "randomize must be of length 1.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs =  c("zz","zz")),
#     "n_runs must be of length 1.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs = "zz"),
#     "n_runs must be numeric.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs = 1.1),
#     "n_runs must be an integer.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs = -1),
#     "n_runs must be strictly higher than 0.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_runs = 0),
#     "n_runs must be strictly higher than 0.",
#     fixed = TRUE)   
#   
#   expect_error(
#     hclu_hierarclust(dissim, keep_trials = 1),
#     "keep_trials must be a boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, keep_trials = c("zz","zz")),
#     "keep_trials must be of length 1.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(dissim, optimal_tree_method = 1),
#     "optimal_tree_method must be a character.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, optimal_tree_method = c("1","1")),
#     "optimal_tree_method must be of length 1.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, optimal_tree_method = "zz"),
#     "Please choose optimal_tree_method among the followings values:
# best",
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
#     hclu_hierarclust(dissim, n_clust = c(1.1,2)),
#     "n_clust must be composed of integer(s).",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = -1),
#     "n_clust must be composed of value(s) strictly higher than 0.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = c(1,-1)),
#     "n_clust must be composed of value(s) strictly higher than 0.",
#     fixed = TRUE) 
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = 0),
#     "n_clust must be composed of value(s) strictly higher than 0.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = 1, cut_height = 1),
#     "Please provide either n_clust or cut_height, but not both at the
#            same time.",
#     fixed = TRUE)   
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = NULL, cut_height = "zz"),
#     "cut_height must be numeric.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, n_clust = NULL, cut_height = -1),
#     "cut_height must be composed of value(s) higher than 0.",
#     fixed = TRUE)   
#   
#   expect_error(
#     hclu_hierarclust(dissim, find_h = 1),
#     "find_h must be a boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     hclu_hierarclust(dissim, find_h = c("zz","zz")),
#     "find_h must be of length 1.",
#     fixed = TRUE)   
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_min =  c("zz","zz")),
#     "h_min must be of length 1.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_min = "zz"),
#     "h_min must be numeric.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_min = -1),
#     "h_min must be higher than 0.",
#     fixed = TRUE)   
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_max =  c("zz","zz")),
#     "h_max must be of length 1.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_max = "zz"),
#     "h_max must be numeric.",
#     fixed = TRUE)  
#   
#   expect_error(
#     hclu_hierarclust(dissim, h_max = -1),
#     "h_max must be higher than 0.",
#     fixed = TRUE)    
#  
# })
