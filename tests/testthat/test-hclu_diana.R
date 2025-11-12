# Inputs -----------------------------------------------------------------------
dissim <- dissimilarity(fishmat, metric = "all")
dissim2 <- dissim[,1:2]

df <- data.frame(ID1=dissim$Site1, ID2=dissim$Site2, W=dissim$Eucl)

d <- dist(fishmat)
mat <- fishmat
rownames(mat) <- NULL
colnames(mat) <- NULL
d2 <- dist(mat)
d3 <- d
d3[1] <- "1"
d4 <- d
d4[1] <- NA


uni <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "c", "b"),
  Weight = c(10, 100, 1))

unina1 <- uni
unina1[1,1] <- NA

unina2 <- uni
unina2[1,3] <- NA

unichar <- uni
unichar$Weight <- unichar$Site1

uni2 <- data.frame(
  Site1 = c("a", "c", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

uni3 <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

uni4 <- data.frame(
  Site1 = c("c", "a", "a"),
  Site2 = c("a", "b", "b"),
  Weight = c(10, 100, 1))

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  clust <- hclu_diana(dissim,
                      index = "Simpson",
                      n_clust = c(1,2,3),
                      cut_height = NULL,
                      find_h = FALSE,
                      h_max = 1,
                      h_min = 0,
                      verbose = FALSE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "hclu_diana")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$n_clust, c(1,2,3))
  expect_equal(clust$args$cut_height, NULL)
  expect_equal(clust$args$find_h, FALSE)
  expect_equal(clust$args$h_max, 1)
  expect_equal(clust$args$h_min, 0)
  expect_equal(clust$args$verbose, FALSE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, TRUE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, TRUE)
  expect_equal(clust$inputs$data_type, "occurrence")
  expect_equal(clust$inputs$node_type, "site")
  expect_equal(sum(attr(clust$clusters, "node_type")=="site"), 
               dim(clust$clusters)[1])
  expect_equal(dim(clust$clusters)[2], 4)

  clust <- hclu_diana(dissim,
                      index = 7,
                      n_clust = c(1,2,3),
                      cut_height = NULL,
                      find_h = FALSE,
                      h_max = 1,
                      h_min = 0,
                      verbose = FALSE)
  expect_equal(clust$args$index, 7)
  expect_equal(clust$inputs$pairwise_metric, "Bray")
  expect_equal(clust$inputs$node_type, "site")
  expect_equal(sum(attr(clust$clusters, "node_type")=="site"), 
               dim(clust$clusters)[1])
  expect_equal(dim(clust$clusters)[2], 4)  
  
  clust <- hclu_diana(dissim,
                      index = "Euclidean",
                      n_clust = c(5,10),
                      find_h = FALSE,
                      verbose = FALSE)
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")
  
  clust <- hclu_diana(dissim,
                      index = "Euclidean",
                      n_clust = c(10,5),
                      find_h = FALSE,
                      verbose = FALSE)
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")
  
  # Test data_type with different dissimilarity metrics
  # We expect warnings here as the algo cannot find 3 clusters with those metrics
  # on this dataset
  expect_warning(clust <- hclu_diana(dissim, 
                                     index = "Simpson", 
                                     n_clust = 3,
                                     verbose = FALSE),)
  expect_equal(clust$inputs$data_type, "occurrence")

  expect_warning(clust <- hclu_diana(dissim, 
                                     index = "Jaccard", 
                                     n_clust = 3,
                                     verbose = FALSE))
  expect_equal(clust$inputs$data_type, "occurrence")

  expect_warning(clust <- hclu_diana(dissim, 
                                     index = "Bray", 
                                     n_clust = 3,
                                     verbose = FALSE))
  expect_equal(clust$inputs$data_type, "abundance")

  expect_warning(clust <- hclu_diana(dissim, 
                                     index = "Euclidean", 
                                     n_clust = 3,
                                     verbose = FALSE))
  expect_equal(clust$inputs$data_type, NA)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    hclu_diana(dissimilarity = "zz"),
    "^dissimilarity is not a bioregion.pairwise object")
  
  expect_error(
    hclu_diana(dissim2),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(uni[,-3]),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(unina1),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(unina2),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(uni4),
    "The first two columns of dissimilarity contain duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(uni3),
    "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(uni2),
    "dissimilarity contains rows with the same site on both columns!",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(unichar),
    "The weight column must be numeric.",
    fixed = TRUE)  
  
  expect_message(
    hclu_diana(d2, verbose = FALSE),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(d3),
    "dissimilarity must be numeric.",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(d4),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(dissim, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(dissim, index = "zz"),
    "^If index is a character, it should be ")
  
  expect_error(
    hclu_diana(dissim, index = "Site1"),
    "^If index is a character, it should be ")
  
  expect_error(
    hclu_diana(dissim, index = 0.1),
    "^If index is numeric, it should be an integer.")
  
  expect_error(
    hclu_diana(dissim, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(uni, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(dissim, n_clust = "zz"),
    "n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from bioregionalization_metrics()",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, n_clust = c(1.1,2)),
    "n_clust must be composed of integers.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, n_clust = -1),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(dissim, n_clust = c(1,-1)),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_diana(dissim, n_clust = 0),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, 
               n_clust = 1, 
               cut_height = 1),
    "^Please provide either n_clust or cut_height,")   
  
  expect_error(
    hclu_diana(dissim, n_clust = NULL, cut_height = "zz"),
    "cut_height must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, n_clust = NULL, cut_height = -1),
    "cut_height must be composed of values higher than 0.",
    fixed = TRUE)   
  
  expect_error(
    hclu_diana(dissim, find_h = 1),
    "find_h must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    hclu_diana(dissim, find_h = c("zz","zz")),
    "find_h must be of length 1.",
    fixed = TRUE)   
  
  expect_error(
    hclu_diana(dissim, h_min =  c("zz","zz")),
    "h_min must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, h_min = "zz"),
    "h_min must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, h_min = -1),
    "h_min must be higher than 0.",
    fixed = TRUE)   
  
  expect_error(
    hclu_diana(dissim, 
               h_max =  c("zz","zz")),
    "h_max must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, 
               h_max = "zz"),
    "h_max must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_diana(dissim, 
               h_max = -1),
    "h_max must be higher than 0.",
    fixed = TRUE)    
  
  expect_error(
    hclu_diana(dissim, 
               h_min = 1,
               h_max = 0),
    "h_min must be inferior to h_max.",
    fixed = TRUE)    
  
  expect_error(
    hclu_diana(dissim, 
               verbose = 1),
    "verbose must be a boolean.",
    fixed = TRUE)   
  
  expect_error(
    hclu_diana(dissim, 
               verbose = c(TRUE, FALSE)),
    "verbose must be of length 1.",
    fixed = TRUE)   
  
})
