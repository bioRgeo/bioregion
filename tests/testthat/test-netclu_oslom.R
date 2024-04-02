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


# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  clust <- netclu_oslom(simil,
                        weight = TRUE,
                        index = 3,
                        reassign = "no",
                        r = 10,
                        hr = 50,
                        seed = 0,
                        t = 0.1,
                        cp = 0.5,
                        directed = FALSE,
                        bipartite = FALSE,
                        site_col = 1,
                        species_col = 2,
                        return_node_type = "both",
                        binpath = "tempdir",
                        path_temp = "oslom_temp",
                        delete_temp = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "netclu_oslom")
  expect_equal(clust$args$weight, TRUE)
  expect_equal(clust$args$index, 3)
  expect_equal(clust$args$reassign, "no")
  expect_equal(clust$args$r, 10)
  expect_equal(clust$args$hr, 50)
  #expect_equal(clust$args$seed, 0)
  expect_equal(clust$args$t, 0.1)
  expect_equal(clust$args$cp, 0.5)
  expect_equal(clust$args$directed, FALSE)
  expect_equal(clust$args$bipartite, FALSE)
  expect_equal(clust$args$site_col, 1)
  expect_equal(clust$args$species_col, 2)
  expect_equal(clust$args$return_node_type, "both")
  #expect_equal(clust$args$binpath, "tempdir")
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
  
  clust <- netclu_oslom(simil,
                        weight = TRUE,
                        index = 3,
                        reassign = "no",
                        r = 10,
                        hr = 50,
                        seed = 1,
                        t = 0.1,
                        cp = 0.5,
                        directed = FALSE,
                        bipartite = FALSE,
                        site_col = 1,
                        species_col = 2,
                        return_node_type = "both",
                        binpath = "tempdir",
                        path_temp = "oslom_temp",
                        delete_temp = TRUE)
  expect_equal(clust$args$seed, 1)
  
  clust2 <- netclu_oslom(simil,
                        weight = TRUE,
                        index = 3,
                        reassign = "no",
                        r = 10,
                        hr = 50,
                        seed = 1,
                        t = 0.1,
                        cp = 0.5,
                        directed = FALSE,
                        bipartite = FALSE,
                        site_col = 1,
                        species_col = 2,
                        return_node_type = "both",
                        binpath = "tempdir",
                        path_temp = "oslom_temp",
                        delete_temp = TRUE)
  expect_equal(clust2$args$seed, 1)
  expect_equal(sum(clust$clusters$K_5==clust2$clusters$K_5), 5)
  
  clust <- netclu_oslom(net)
  expect_equal(dim(clust$clusters)[1], 7)
  
  clust <- netclu_oslom(net, bipartite = TRUE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(clust$args$return_node_type, "both")
  
  clust <- netclu_oslom(net, 
                        bipartite = TRUE, 
                        return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(clust$args$return_node_type, "species")
  
  clust <- netclu_oslom(net, 
                        bipartite = TRUE, 
                        return_node_type = "sites")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(clust$args$return_node_type, "sites")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    netclu_oslom(net, binpath = 1),
    "binpath must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, binpath = c("zz","zz")),
    "binpath must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, path_temp = 1),
    "path_temp must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, delete_temp = c("zz","zz")),
    "delete_temp must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, delete_temp =  1),
    "delete_temp must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, binpath = c("zz","zz")),
    "binpath must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, directed = "zz"),
    "directed must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_infomap(net, directed = c("zz",1)),
    "directed must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, reassign = 1),
    "reassign must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, reassign = c("zz","zz")),
    "reassign must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, reassign = "zz"),
    "Please choose reassign among the following values: 
no, random or simil",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, r =  c("zz","zz")),
    "r must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, r = "zz"),
    "r must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, r = 1.1),
    "r must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, r = -1),
    "r must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, r = 0),
    "r must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, hr =  c("zz","zz")),
    "hr must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, hr = "zz"),
    "hr must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, hr = 1.1),
    "hr must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, hr = -1),
    "hr must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, seed = -1),
    "seed must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, t =  c("zz","zz")),
    "t must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, t = "zz"),
    "t must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, t = -1),
    "t must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, t = 0),
    "t must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, t = 1),
    "t must be in the interval (0,1)!",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, cp =  c("zz","zz")),
    "cp must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    netclu_oslom(net, cp = "zz"),
    "cp must be numeric.",
    fixed = TRUE)  

  expect_error(
    netclu_oslom(net, cp = -1),
    "cp must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, cp = 0),
    "cp must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, cp = 1),
    "cp must be in the interval (0,1)!",
    fixed = TRUE) 
  
  expect_error(
    netclu_oslom(net, bipartite = 1),
    "bipartite must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = c("zz","zz")),
    "bipartite must be of length 1.",
    fixed = TRUE)
  
#   expect_message(
#     netclu_oslom(net, bipartite = FALSE),
#     "net is not a bioregion.pairwise.metric object. 
# Note that some functions required dissimilarity metrics (hclu_ & nhclu_) and
# others similarity metrics (netclu_). 
# Please carefully check your data before using the clustering functions.",
#     fixed = TRUE)
  
  expect_error(
    netclu_oslom(dissimil),
    "net seems to be a dissimilarity object. 
This function should be applied on similarities, not dissimilarities. 
Use dissimilarity_to_similarity() before using this function.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net[,-3], weight = TRUE),
    "net must be a data.frame with at least three columns if weight equal 
        TRUE.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = "zz"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = "Site1"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = 2),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(simil, index = 1),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(simil, index = 16),
    "index should be lower or equal to 15.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, weight = FALSE, reassign = "simil"),
    "A reassignement based on similarity should not be use when weight
           equal FALSE", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(uni, bipartite = TRUE),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, return_node_type = "zz"),
    "Please choose return_node_type among the followings values:
both, sites or species",
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, site_col = "Weight", 
                   species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, directed = "zz"),
    "directed must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, directed = c("zz",1)),
    "directed must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_oslom(uni2),
    "The network contains self-loop(s)!"
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(uni3, directed = FALSE),
    "net should not be directed if directed = FALSE."
    , fixed = TRUE)
  
  expect_error(
    netclu_oslom(net, bipartite = TRUE, directed = TRUE),
    "directed cannot be set to TRUE if the network is bipartite!"
    , fixed = TRUE)
  
})
