# Inputs -----------------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

uni <- data.frame(
  Site1 = c("b", "a"),
  Site2 = c("a", "b"),
  Weight = c(10, 100))

simil <- similarity(net_to_mat(net), metric = c("abc"))

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
  
  clust <- netclu_beckett(net,
                          weight = TRUE,
                          index = names(net)[3],
                          forceLPA = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "netclu_beckett")
  expect_equal(clust$args$weight, TRUE)
  expect_equal(clust$args$index, "Weight")
  expect_equal(clust$args$forceLPA, FALSE)
  expect_equal(clust$args$site_col, 1)
  expect_equal(clust$args$species_col, 2)
  expect_equal(clust$args$return_node_type, "both")
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, TRUE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, FALSE)
  expect_equal(clust$inputs$pairwise_metric, NA)
  expect_equal(clust$inputs$dissimilarity, FALSE)
  expect_equal(clust$inputs$nb_sites, 3)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(dim(clust$clusters)[1], 7)
  expect_equal(dim(clust$clusters)[2], 2)
  
  clust2 <- netclu_beckett(net,
                           weight = TRUE,
                           index = names(net)[3],
                           forceLPA = FALSE,
                           site_col = 1,
                           species_col = 2,
                           return_node_type = "both",
                           algorithm_in_output = TRUE)
  expect_equal(sum(clust$clusters$K_2==clust2$clusters$K_2), 7)

  clust <- netclu_beckett(net, return_node_type = "species")
  expect_equal(dim(clust$clusters)[1], 4)
  expect_equal(dim(clust$clusters)[2], 2)
  expect_equal(clust$args$return_node_type, "species")
  
  clust <- netclu_beckett(net, return_node_type = "sites")
  expect_equal(dim(clust$clusters)[1], 3)
  expect_equal(dim(clust$clusters)[2], 2)
  expect_equal(clust$args$return_node_type, "sites")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    netclu_beckett(net, forceLPA = "zz"),
    "forceLPA must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, forceLPA = c("zz",1)),
    "forceLPA must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, algorithm_in_output = "zz"),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, algorithm_in_output = c("zz",1)),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(uni),
    "The network is not bipartite!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, return_node_type = 1),
    "return_node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, return_node_type = "zz"),
    "Please choose return_node_type among the followings values:
both, sites or species",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, site_col = "Weight", species_col = "Weight"),
    "site_col and species_col should not be the same."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, site_col = "zz"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, species_col = "zz"),
    "If species_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, site_col = "Weight"),
    "If site_col is a character, it should be the first or second column name."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, site_col = 3),
    "If site_col is numeric, it should be equal to 1 or 2."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, site_col = FALSE),
    "site_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, species_col = FALSE),
    "species_col should be numeric or character."
    , fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net[,-3], weight = TRUE),
    "net must be a data.frame with at least three columns if weight equal 
        TRUE.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)

  expect_error(
    netclu_beckett(net, index = "zz"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, index = "Site1"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, index = 2),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE, index = 3),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE, index = "Weight"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE, index = 5),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
  expect_error(
    netclu_beckett(net6, weight = TRUE, index = "Weight3"),
    "The weight column should contain only positive reals:
          negative value(s) detected!", 
    fixed = TRUE)
  
})
