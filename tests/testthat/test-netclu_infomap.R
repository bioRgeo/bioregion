# Inputs -----------------------------------------------------------------------
quietly(install_binaries(verbose = FALSE))

net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

simil <- similarity(comat, metric = "all")
similnull <- simil
colnames(similnull) <- NULL
similna <- simil
colnames(similna) <- NA

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
  
  clust <- netclu_infomap(simil,
                          weight = TRUE,
                          cut_weight = 0,
                          index = 3,
                          seed = NULL,
                          nbmod = 0,
                          markovtime = 1,
                          numtrials = 1,
                          twolevel = FALSE,
                          show_hierarchy = FALSE,
                          directed = FALSE,
                          bipartite_version = FALSE,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          version = "2.8.0",
                          binpath = "tempdir",
                          check_install = TRUE,
                          path_temp = "infomap_temp",
                          delete_temp = TRUE)
  expect_equal(clust$args$index, 3)
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$nbmod, 0)
  expect_equal(clust$args$markovtime, 1)
  expect_equal(clust$args$numtrials, 1)
  expect_equal(clust$args$twolevel, FALSE)
  expect_equal(clust$args$show_hierarchy, FALSE)
  expect_equal(clust$args$directed, FALSE)
  expect_equal(clust$args$bipartite_version, FALSE)
  expect_equal(clust$args$bipartite, FALSE)
  expect_equal(clust$args$site_col, 1)
  expect_equal(clust$args$species_col, 2)
  expect_equal(clust$args$return_node_type, "both")
  expect_equal(clust$args$version, "2.8.0")
  #expect_equal(clust$args$binpath, "tempdir")
  expect_equal(clust$args$check_install, TRUE)
  #expect_equal(clust$args$path_temp, "infomap_temp")
  expect_equal(clust$args$delete_temp, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Jaccard")
  expect_equal(clust$inputs$dissimilarity, FALSE)
  expect_equal(clust$inputs$nb_sites, 5)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  expect_equal(clust$inputs$node_type, "site")
  expect_equal(sum(attr(clust$clusters, "node_type")=="site"), 
               dim(clust$clusters)[1])
  expect_equal(clust$algorithm$version, "2.8.0")
  expect_equal(dim(clust$clusters)[1], 5)
  
  clust <- netclu_infomap(simil,
                          weight = TRUE,
                          index = 5)
  expect_equal(clust$args$index, 5)
  expect_equal(clust$inputs$pairwise_metric, "Sorensen")
  
  clust <- netclu_infomap(simil,
                          weight = TRUE,
                          index = "Simpson")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  
  clust <- netclu_infomap(similna,
                          weight = TRUE,
                          index = 3)
  expect_equal(clust$args$index, 3)
  expect_equal(is.na(clust$inputs$pairwise_metric), TRUE)
  
  clust <- netclu_infomap(similna,
                          weight = TRUE,
                          index = 7)
  expect_equal(clust$args$index, 7)
  expect_equal(is.na(clust$inputs$pairwise_metric), TRUE)
  
  clust <- netclu_infomap(simil,
                          weight = TRUE,
                          index = 3,
                          seed = 1,
                          nbmod = 0,
                          markovtime = 1,
                          numtrials = 1,
                          twolevel = FALSE,
                          show_hierarchy = FALSE,
                          directed = FALSE,
                          bipartite_version = FALSE,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          version = "2.8.0",
                          binpath = "tempdir",
                          check_install = TRUE,
                          path_temp = "infomap_temp",
                          delete_temp = TRUE)
  expect_equal(clust$args$seed, 1)
  
  clust2 <- netclu_infomap(simil,
                          weight = TRUE,
                          index = 3,
                          nbmod = 0,
                          markovtime = 1,
                          seed = 1,
                          numtrials = 1,
                          twolevel = FALSE,
                          show_hierarchy = FALSE,
                          directed = FALSE,
                          bipartite_version = FALSE,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          version = "2.8.0",
                          binpath = "tempdir",
                          check_install = TRUE,
                          path_temp = "infomap_temp",
                          delete_temp = TRUE)
  expect_equal(clust2$args$seed, 1)
  expect_equal(sum(clust$clusters$K_1==clust2$clusters$K_1), 5)
  
  clust <- netclu_infomap(net)
  expect_equal(dim(clust$clusters)[1], 7)
  
  clust <- netclu_infomap(net, 
                          bipartite = TRUE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(clust$args$return_node_type, "both")
  
  clust <- netclu_infomap(net, 
                          bipartite = TRUE, 
                          return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(clust$args$return_node_type, "species")
  expect_equal(clust$inputs$node_type, "species")
  
  clust <- netclu_infomap(net, 
                          bipartite = TRUE, 
                          return_node_type = "site")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(clust$args$return_node_type, "site")
  expect_equal(clust$inputs$node_type, "site")
  
  clust <- netclu_infomap(net, cut_weight = 40, seed = 1)
  expect_equal(colnames(clust$clusters), c("ID","K_2"))
  expect_equal(length(table(clust$clusters$K_2)), 2)
  expect_equal(clust$cluster_info[1,1], "K_2")
  expect_equal(clust$cluster_info[1,2], 2)
  
  clust <- netclu_infomap(net, cut_weight = 60, seed = 1)
  expect_equal(colnames(clust$clusters), c("ID","K_1"))
  expect_equal(length(table(clust$clusters$K_1)), 1)
  expect_equal(clust$cluster_info[1,1], "K_1")
  expect_equal(clust$cluster_info[1,2], 1)
  
  clust1 <- netclu_infomap(fdf, seed = 1)
  clust2 <- netclu_infomap(fdf, seed = 1)
  expect_equal(clust1$args$seed==clust2$args$seed, TRUE)
  
  clust1 <- netclu_infomap(fdf, seed = 1)
  clust2 <- netclu_infomap(fdf, seed = 1)
  expect_equal(sum(clust1$clusters$K_8==clust2$clusters$K_8), 266)
  expect_equal(sum(clust1$clusters$K_10==clust2$clusters$K_10), 266)
  
  clust1 <- netclu_infomap(vdf, seed = 1)
  clust2 <- netclu_infomap(vdf, seed = 1)
  expect_equal(sum(clust1$clusters$K_2==clust2$clusters$K_2), 873)
  
  r1 <- runif(1)
  clust1 <- netclu_infomap(vdf, seed = NULL)
  r2 <- runif(1)
  clust2 <- netclu_infomap(vdf, seed = NULL)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_infomap(vdf, seed = 1)
  r2 <- runif(1)
  clust2 <- netclu_infomap(vdf, seed = 1)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_infomap(vdf, seed = 1000)
  r2 <- runif(1)
  clust2 <- netclu_infomap(vdf, seed = 1000)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  clust <- netclu_infomap(simf, seed = 1, show_hierarchy = FALSE)
  expect_equal(colnames(clust$clusters), c("ID","K_6","K_7"))
  expect_equal(length(table(clust$clusters$K_6)), 6)
  expect_equal(length(table(clust$clusters$K_7)), 7)
  expect_equal(clust$cluster_info[1,1], "K_6")
  expect_equal(clust$cluster_info[1,2], 6)
  expect_equal(clust$cluster_info[1,3], 1)
  expect_equal(clust$cluster_info[2,1], "K_7")
  expect_equal(clust$cluster_info[2,2], 7)
  expect_equal(clust$cluster_info[2,3], 2)
  expect_equal(clust$inputs$hierarchical, TRUE)
  
  # Test data_type with bipartite network (weighted)
  clust <- netclu_infomap(net, bipartite = TRUE, weight = TRUE)
  expect_equal(clust$inputs$data_type, "abundance")
  
  # Test data_type with bipartite network (unweighted)
  clust <- netclu_infomap(net, bipartite = TRUE, weight = FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  # Test data_type with similarity metrics (occurrence-based)
  clust <- netclu_infomap(simil, index = "Jaccard")
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- netclu_infomap(simil, index = "Simpson")
  expect_equal(clust$inputs$data_type, "occurrence")
  
  # Test data_type with similarity metrics (abundance-based)
  clust <- netclu_infomap(simil, index = "Bray")
  expect_equal(clust$inputs$data_type, "abundance")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {

  expect_error(
    netclu_infomap(net, binpath = 1),
    "binpath must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, binpath = c("zz","zz")),
    "binpath must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, check_install = 1),
    "check_install must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, check_install = c(TRUE,FALSE)),
    "check_install must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, path_temp = 1),
    "path_temp must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, path_temp = c("zz","zz")),
    "path_temp must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, delete_temp = c("zz","zz")),
    "delete_temp must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, delete_temp =  1),
    "delete_temp must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, version = 1),
    "version must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, version = c("zz","zz")),
    "version must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, nbmod = c("zz","zz")),
    "nbmod must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, nbmod = "zz"),
    "nbmod must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, nbmod = 1.1),
    "nbmod must be an integer.",
    fixed = TRUE)  

  expect_error(
    netclu_infomap(net, nbmod = -1),
    "nbmod must be higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, markovtime =  c("zz","zz")),
    "markovtime must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, markovtime = "zz"),
    "markovtime must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, markovtime = -1.1),
    "markovtime must be strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, markovtime = 0),
    "markovtime must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_infomap(net, seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, seed = -1),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_infomap(net, seed = 0),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_infomap(net, numtrials =  c("zz","zz")),
    "numtrials must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, numtrials = "zz"),
    "numtrials must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, numtrials = 1.1),
    "numtrials must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, numtrials = -1),
    "numtrials must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_infomap(net, numtrials = 0),
    "numtrials must be strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, twolevel = 1),
    "twolevel must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, twolevel = c("zz","zz")),
    "twolevel must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, show_hierarchy = 1),
    "show_hierarchy must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, show_hierarchy = c("zz","zz")),
    "show_hierarchy must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite_version = 1),
    "bipartite_version must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite_version = c("zz","zz")),
    "bipartite_version must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
#   expect_message(
#     netclu_infomap(net, bipartite = FALSE),
#     "net is not a bioregion.pairwise object. 
# Note that some functions required dissimilarity metrics (hclu_ & nhclu_) and
# others similarity metrics (netclu_). 
# Please carefully check your data before using the clustering functions.",
#     fixed = TRUE)
  
  expect_error(
    netclu_infomap(dissimil),
    "^net seems to be a dissimilarity object")
  
  expect_error(
    netclu_infomap("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, cut_weight =  c("zz","zz")),
    "cut_weight must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, cut_weight = "zz"),
    "cut_weight must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_infomap(net, cut_weight = -1),
    "cut_weight must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net[,-3], weight = TRUE),
    "^net must be a data.frame with at least three columns if ")
  
  expect_error(
    netclu_infomap(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, index = "zz"),
    "^If index is a character, it should be a column ")
  
  expect_error(
    netclu_infomap(net, index = "Site1"),
    "^If index is a character, it should be a column ")
  
  expect_error(
    netclu_infomap(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(simil, index = 1),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(simil, index = 16),
    "index should be lower or equal to 15.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net6, weight = TRUE),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(uni, bipartite = TRUE),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, return_node_type = "zz"),
    "^Please choose return_node_type from the following:")
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, site_col = "Weight", 
                   species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, directed = "zz"),
    "directed must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, directed = c("zz",1)),
    "directed must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(uni2),
    "The network contains self-loop(s)!"
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(uni3, directed = FALSE),
    "net should not be directed if directed = FALSE."
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, bipartite = TRUE, directed = TRUE),
    "directed cannot be set to TRUE if the network is bipartite!"
    , fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, cut_weight = 100),
    "^The network is empty. ")
  
})
