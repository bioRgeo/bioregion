# Inputs -----------------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

simil <- similarity(comat, metric = "all")
dissimil <- dissimilarity(comat, metric = "all")

uni <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "c", "b"),
  Weight = c(10, 100, 1))

uni2 <- data.frame(
  Site1 = c("a", "b", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

uni3 <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

net2 <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c("a", "b", "a", "c", "d", "b", "d")
)

net3 <- rbind(net, net)

net4 <- net
net4[1,1] <- NA

net5 <- net
net5[1,3] <- NA

net6 <- data.frame(
  Obj1 = c("a","b","c","a"),
  Obj2 = c("A","B","C","V"),
  Weight = c(1,-1,0,2),
  Weight2 = c(1,1,0,2),
  Weight3 = c(1,-1,0,2)
)

fdf <- fishdf[1:1000,]
vdf <- vegedf[1:1000,]
simf <- similarity(fishmat, metric = "all")


# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  clust <- netclu_leiden(simil,
                         weight = TRUE,
                         cut_weight = 0,
                         index = 3,      
                         seed = NULL,
                         objective_function = "CPM",
                         resolution_parameter = 1,
                         beta = 0.01,
                         n_iterations = 2,
                         vertex_weights = NULL,
                         bipartite = FALSE,
                         site_col = 1,
                         species_col = 2,
                         return_node_type = "both",
                         algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "netclu_leiden")
  expect_equal(clust$args$weight, TRUE)
  expect_equal(clust$args$cut_weight, 0)
  expect_equal(clust$args$index, 3)
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$objective_function, "CPM")
  expect_equal(clust$args$resolution_parameter, 1)
  expect_equal(clust$args$beta, 0.01)
  expect_equal(clust$args$n_iterations, 2)
  expect_equal(clust$args$vertex_weights, NULL)
  expect_equal(clust$args$bipartite, FALSE)
  expect_equal(clust$args$site_col, 1)
  expect_equal(clust$args$species_col, 2)
  expect_equal(clust$args$return_node_type, "both")
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Jaccard")
  expect_equal(clust$inputs$dissimilarity, FALSE)
  expect_equal(clust$inputs$nb_sites, 5)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(dim(clust$clusters)[1], 5)
  
  clust <- netclu_leiden(simil,
                         weight = FALSE,
                         index = 3,
                         bipartite = FALSE,
                         site_col = 1,
                         species_col = 2,
                         return_node_type = "both",
                         algorithm_in_output = TRUE)
  expect_equal(clust$args$weight, FALSE)
  expect_equal(clust$inputs$weight, FALSE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, NA)
  expect_equal(clust$inputs$dissimilarity, FALSE)
  
  clust <- netclu_leiden(net, bipartite = TRUE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(clust$args$return_node_type, "both")
  
  clust <- netclu_leiden(net, 
                         bipartite = TRUE, 
                         return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(clust$args$return_node_type, "species")
  
  clust <- netclu_leiden(net, 
                         bipartite = TRUE, 
                         return_node_type = "sites")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(clust$args$return_node_type, "sites")
  
  clust <- netclu_leiden(net, cut_weight = 100)
  expect_equal(colnames(clust$clusters), c("ID","K_0"))
  expect_equal(length(table(clust$clusters$K_0)), 0)
  expect_equal(clust$cluster_info[1,1], "K_0")
  expect_equal(clust$cluster_info[1,2], 0)
  
  clust1 <- netclu_leiden(fdf, seed = 1)
  clust2 <- netclu_leiden(fdf, seed = 1)
  expect_equal(sum(clust1$clusters$K_266==clust2$clusters$K_266), 266)
  
  clust1 <- netclu_leiden(vdf, seed = 2)
  clust2 <- netclu_leiden(vdf, seed = 2)
  expect_equal(sum(clust1$clusters$K_763==clust2$clusters$K_763), 873)
  
  clust1 <- netclu_leiden(simf, seed = 3)
  clust2 <- netclu_leiden(simf, seed = 3)
  expect_equal(sum(clust1$clusters$K_338==clust2$clusters$K_338), 338)
  
  r1 <- runif(1)
  clust1 <- netclu_leiden(vdf)
  r2 <- runif(1)
  clust2 <- netclu_leiden(vdf)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_leiden(vdf, seed = 1)
  r2 <- runif(1)
  clust2 <- netclu_leiden(vdf, seed = 1)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_leiden(vdf, seed = 1000)
  r2 <- runif(1)
  clust2 <- netclu_leiden(vdf, seed = 1000)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    netclu_leiden(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
#   expect_message(
#     netclu_leiden(net, bipartite = FALSE),
#     "net is not a bioregion.pairwise.metric object. 
# Note that some functions required dissimilarity metrics (hclu_ & nhclu_) and
# others similarity metrics (netclu_). 
# Please carefully check your data before using the clustering functions.",
#     fixed = TRUE)
  
  expect_error(
    netclu_leiden(dissimil),
    "net seems to be a dissimilarity object. 
This function should be applied on similarities, not dissimilarities. 
Use dissimilarity_to_similarity() before using this function.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, cut_weight =  c("zz","zz")),
    "cut_weight must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, cut_weight = "zz"),
    "cut_weight must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, cut_weight = -1),
    "cut_weight must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net[,-3], weight = TRUE),
    "net must be a data.frame with at least three columns if weight equal 
        TRUE.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = "zz"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = "Site1"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = 2),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(simil, index = 1),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(simil, index = 16),
    "index should be lower or equal to 15.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, seed = -1),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_leiden(net, seed = 0),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_leiden(uni, bipartite = TRUE),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, return_node_type = "zz"),
    "Please choose return_node_type among the followings values:
both, sites or species",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, site_col = "Weight", 
                        species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(uni2),
    "The network contains self-loop(s)!"
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(uni3),
    "The network is directed, this function is designed for undirected networks!"
    , fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, objective_function = 1),
    "objective_function must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, objective_function = c("zz","zz")),
    "objective_function must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, objective_function = "zz"),
    "Please choose objective_function among the following values: 
CPM or modularity",
    fixed = TRUE)
  
  expect_error(
    netclu_leiden(net, resolution_parameter =  c("zz","zz")),
    "resolution_parameter must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, resolution_parameter = "zz"),
    "resolution_parameter must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, resolution_parameter = 1.1),
    "resolution_parameter must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, resolution_parameter = -1),
    "resolution_parameter must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_leiden(net, beta =  c("zz","zz")),
    "beta must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, beta = "zz"),
    "beta must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, n_iterations =  c("zz","zz")),
    "n_iterations must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, n_iterations = "zz"),
    "n_iterations must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, n_iterations = 1.1),
    "n_iterations must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_leiden(net, n_iterations = -1),
    "n_iterations must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_leiden(net, n_iterations = 0),
    "n_iterations must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_leiden(net, vertex_weights = "zz"),
    "vertex_weights must be numeric.",
    fixed = TRUE) 
  
})
