# Inputs -----------------------------------------------------------------------
install_binaries()

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
  
  clust <- netclu_louvain(simil,
                          weight = TRUE,
                          cut_weight = 0,
                          index = 3,
                          lang = "igraph",
                          resolution = 1,
                          seed = NULL,
                          q = 0,
                          c = 0.5,
                          k = 1,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          binpath = "tempdir",
                          check_install = TRUE,
                          path_temp = "louvain_temp",
                          delete_temp = TRUE,
                          algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "netclu_louvain")
  expect_equal(clust$args$weight, TRUE)
  expect_equal(clust$args$cut_weight, 0)
  expect_equal(clust$args$index, 3)
  expect_equal(clust$args$lang, "igraph")
  expect_equal(clust$args$resolution, 1)
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$q, 0)
  expect_equal(clust$args$k, 1)
  expect_equal(clust$args$bipartite, FALSE)
  expect_equal(clust$args$site_col, 1)
  expect_equal(clust$args$species_col, 2)
  expect_equal(clust$args$return_node_type, "both")
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
  expect_equal(dim(clust$clusters)[1], 5)
  
  clust <- netclu_louvain(net)
  expect_equal(dim(clust$clusters)[1], 7)
  
  clust <- netclu_louvain(net, bipartite = TRUE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(clust$args$return_node_type, "both")
  
  clust <- netclu_louvain(net, 
                          bipartite = TRUE, 
                          return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(clust$args$return_node_type, "species")
  
  clust <- netclu_louvain(net, 
                          bipartite = TRUE, 
                          return_node_type = "site")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(clust$args$return_node_type, "site")
  
  clust <- netclu_louvain(net, cut_weight = 100, lang = "igraph")
  expect_equal(colnames(clust$clusters), c("ID","K_0"))
  expect_equal(length(table(clust$clusters$K_0)), 0)
  expect_equal(clust$cluster_info[1,1], "K_0")
  expect_equal(clust$cluster_info[1,2], 0)
  
  clust <- netclu_louvain(net, cut_weight = 60, lang = "igraph")
  expect_equal(colnames(clust$clusters), c("ID","K_1"))
  expect_equal(length(table(clust$clusters$K_1)), 1)
  expect_equal(clust$cluster_info[1,1], "K_1")
  expect_equal(clust$cluster_info[1,2], 1)
  
  clust1 <- netclu_louvain(fdf, lang = "igraph", seed=1)
  clust2 <- netclu_louvain(fdf, lang = "igraph", seed=1)
  expect_equal(sum(clust1$clusters$K_6==clust2$clusters$K_6), 266)
  
  clust1 <- netclu_louvain(vdf, lang = "igraph", seed=2)
  clust2 <- netclu_louvain(vdf, lang = "igraph", seed=2)
  expect_equal(sum(clust1$clusters$K_3==clust2$clusters$K_3), 873)
  
  clust1 <- netclu_louvain(simf, lang = "igraph", seed=1)
  clust2 <- netclu_louvain(simf, lang = "igraph", seed=1)
  expect_equal(sum(clust1$clusters$K_4==clust2$clusters$K_4), 338)
  
  r1 <- runif(1)
  clust1 <- netclu_louvain(vdf, lang = "igraph")
  r2 <- runif(1)
  clust2 <- netclu_louvain(vdf, lang = "igraph")
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_louvain(vdf, lang = "igraph", seed = 1)
  r2 <- runif(1)
  clust2 <- netclu_louvain(vdf, lang = "igraph", seed = 1)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  r1 <- runif(1)
  clust1 <- netclu_louvain(vdf, lang = "igraph", seed = 1000)
  r2 <- runif(1)
  clust2 <- netclu_louvain(vdf, lang = "igraph", seed = 1000)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  clust <- netclu_louvain(net, cut_weight = 60, lang = "cpp")
  expect_equal(colnames(clust$clusters), c("ID","K_1"))
  expect_equal(length(table(clust$clusters$K_1)), 1)
  expect_equal(clust$cluster_info[1,1], "K_1")
  expect_equal(clust$cluster_info[1,2], 1)
  
  clust <- netclu_louvain(fdf, lang = "cpp")
  expect_equal(clust$inputs$hierarchical, TRUE)

  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    netclu_louvain(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
#   expect_message(
#     netclu_louvain(net, bipartite = FALSE),
#     "net is not a bioregion.pairwise object. 
# Note that some functions required dissimilarity metrics (hclu_ & nhclu_) and
# others similarity metrics (netclu_). 
# Please carefully check your data before using the clustering functions.",
#     fixed = TRUE)
  
  expect_error(
    netclu_louvain(dissimil),
    "^net seems to be a dissimilarity object")
  
  expect_error(
    netclu_louvain("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, cut_weight =  c("zz","zz")),
    "cut_weight must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, cut_weight = "zz"),
    "cut_weight must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, cut_weight = -1),
    "cut_weight must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net[,-3], weight = TRUE),
    "^net must be a data.frame with at least three columns if")
  
  expect_error(
    netclu_louvain(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, index = "zz"),
    "^If index is a character, it should be a ")
  
  expect_error(
    netclu_louvain(net, index = "Site1"),
    "^If index is a character, it should ")
  
  expect_error(
    netclu_louvain(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(simil, index = 1),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(simil, index = 16),
    "index should be lower or equal to 15.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net6, weight = TRUE),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive values.", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(uni, bipartite = TRUE),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, return_node_type = "zz"),
    "^Please choose return_node_type from the following:")
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, site_col = "Weight", 
                 species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(uni2),
    "The network contains self-loop(s)!"
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(uni3),
    "The network is directed, this function is designed for undirected networks!"
    , fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang = 1),
    "lang must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang = c("zz","zz")),
    "lang must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang = "zz"),
    "^Please choose lang from the following:")
  
  expect_error(
    netclu_louvain(net, resolution =  c("zz","zz")),
    "resolution must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, resolution = "zz"),
    "resolution must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, resolution = -1.1),
    "resolution must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, resolution = 0),
    "resolution must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, seed = -1),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, seed = 0),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, q =  c("zz","zz")),
    "q must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, q = "zz"),
    "q must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, q = 1.1),
    "q must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, q = -1),
    "q must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, c =  c("zz","zz")),
    "c must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, c = "zz"),
    "c must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, c = -1),
    "c must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, c = 0),
    "c must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, c = 1),
    "c must be in the interval (0,1)!",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, k =  c("zz","zz")),
    "k must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, k = "zz"),
    "k must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_louvain(net, k = -1),
    "k must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, k = 0),
    "k must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_louvain(net, algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, cut_weight = 100, lang = "cpp"),
    "^The network is empty.")
  
  expect_error(
    netclu_louvain(net, lang="cpp", binpath = 1),
    "binpath must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", binpath = c("zz","zz")),
    "binpath must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", check_install = 1),
    "check_install must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", check_install = c(TRUE,FALSE)),
    "check_install must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", path_temp = 1),
    "path_temp must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", path_temp = c("zz","zz")),
    "path_temp must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", delete_temp = c("zz","zz")),
    "delete_temp must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_louvain(net, lang="cpp", delete_temp =  1),
    "delete_temp must be a boolean.",
    fixed = TRUE)
  
})
