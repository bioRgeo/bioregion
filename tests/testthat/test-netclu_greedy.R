# Inputs -----------------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

#comat[1,]=0

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
  
  clust <- netclu_greedy(simil,
                         weight = TRUE,
                         cut_weight = 0,
                         index = 3,
                         bipartite = FALSE,
                         site_col = 1,
                         species_col = 2,
                         return_node_type = "both",
                         algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "netclu_greedy")
  expect_equal(clust$args$weight, TRUE)
  expect_equal(clust$args$cut_weight, 0)
  expect_equal(clust$args$index, 3)
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
  
  clust <- netclu_greedy(simil,
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
  
  clust2 <- netclu_greedy(simil,
                         weight = FALSE,
                         index = 3,
                         bipartite = FALSE,
                         site_col = 1,
                         species_col = 2,
                         return_node_type = "both",
                         algorithm_in_output = TRUE)
  expect_equal(sum(clust$clusters$K_2==clust2$clusters$K_2), 5)
  
  clust <- netclu_greedy(net, bipartite = TRUE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(clust$args$return_node_type, "both")
  
  clust <- netclu_greedy(net, 
                          bipartite = TRUE, 
                          return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(clust$args$return_node_type, "species")
  
  clust <- netclu_greedy(net, 
                          bipartite = TRUE, 
                          return_node_type = "site")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(clust$args$return_node_type, "site")
  
  clust <- netclu_greedy(net, cut_weight = 100)
  expect_equal(colnames(clust$clusters), c("ID","K_0"))
  expect_equal(length(table(clust$clusters$K_0)), 0)
  expect_equal(clust$cluster_info[1,1], "K_0")
  expect_equal(clust$cluster_info[1,2], 0)
  
  clust1 <- netclu_greedy(fdf)
  clust2 <- netclu_greedy(fdf)
  expect_equal(sum(clust1$clusters$K_6==clust2$clusters$K_6), 266)
  
  clust1 <- netclu_greedy(vdf)
  clust2 <- netclu_greedy(vdf)
  expect_equal(sum(clust1$clusters$K_3==clust2$clusters$K_3), 873)
  
  clust1 <- netclu_greedy(simf)
  clust2 <- netclu_greedy(simf)
  expect_equal(sum(clust1$clusters$K_3==clust2$clusters$K_3), 338)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    netclu_greedy(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
#   expect_message(
#     netclu_greedy(net, bipartite = FALSE),
#     "net is not a bioregion.pairwise.metric object. 
# Note that some functions required dissimilarity metrics (hclu_ & nhclu_) and
# others similarity metrics (netclu_). 
# Please carefully check your data before using the clustering functions.",
#     fixed = TRUE)
  
  expect_error(
    netclu_greedy(dissimil),
    "^net seems to be a dissimilarity object")
  
  expect_error(
    netclu_greedy("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, cut_weight =  c("zz","zz")),
    "cut_weight must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_greedy(net, cut_weight = "zz"),
    "cut_weight must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_greedy(net, cut_weight = -1),
    "cut_weight must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net[,-3], weight = TRUE),
    "^net must be a data.frame with at least three columns ")
  
  expect_error(
    netclu_greedy(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, index = "zz"),
    "^If index is a character, it should be a column ")
  
  expect_error(
    netclu_greedy(net, index = "Site1"),
    "^If index is a character, it should be a column ")
  
  expect_error(
    netclu_greedy(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(simil, index = 1),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(simil, index = 16),
    "index should be lower or equal to 15.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net6, weight = TRUE),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(uni, bipartite = TRUE),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, return_node_type = "zz"),
    "^Please choose return_node_type from the following:")
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, site_col = "Weight", 
                   species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(uni2),
    "The network contains self-loop(s)!"
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(uni3),
    "The network is directed, this function is designed for undirected networks!"
    , fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_greedy(net, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
})
