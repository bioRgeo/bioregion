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
  
  clust <- nhclu_kmeans(dissim,
                       index = "Simpson",
                       seed = NULL,
                       n_clust = c(1,2,3),
                       iter_max = 10,
                       nstart = 10,
                       algorithm = "Hartigan-Wong",
                       algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "nhclu_kmeans")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$n_clust, c(1,2,3))
  expect_equal(clust$args$iter_max, 10)
  expect_equal(clust$args$nstart, 10)
  expect_equal(clust$args$algorithm, "Hartigan-Wong")
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, TRUE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  expect_equal(dim(clust$clusters)[2], 4)
  
  clust <- nhclu_kmeans(dissim,
                       index = "Euclidean",
                       n_clust = c(5,10))
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")
  
  clust <- nhclu_kmeans(dissim,
                       index = "Euclidean",
                       n_clust = c(10,5))
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")
  
  clust1 <- nhclu_kmeans(dissim,
                      index = "Euclidean",
                      n_clust = 5,
                      seed = 1)
  clust2 <- nhclu_kmeans(dissim,
                      index = "Euclidean",
                      n_clust = 5,
                      seed =1)
  expect_equal(sum(clust1$clusters$K_5==clust2$clusters$K_5), 338)
  
  r1 <- runif(1)
  clust1 <- nhclu_kmeans(dissim,
                      index = "Euclidean",
                      n_clust = 5,
                      seed = 1)
  r2 <- runif(1)
  clust2 <- nhclu_kmeans(dissim,
                      index = "Euclidean",
                      n_clust = 5,
                      seed = 1)
  r3 <- runif(1)
  expect_equal(r1!=r2, TRUE)
  expect_equal(r2!=r3, TRUE)
  expect_equal(r1!=r3, TRUE)
  
  # Test data_type with different dissimilarity metrics
  clust <- nhclu_kmeans(dissim, index = "Simpson", n_clust = 3)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_kmeans(dissim, index = "Jaccard", n_clust = 3)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_kmeans(dissim, index = "Bray", n_clust = 3)
  expect_equal(clust$inputs$data_type, "abundance")
  
  clust <- nhclu_kmeans(dissim, index = "Euclidean", n_clust = 3)
  expect_equal(clust$inputs$data_type, "unknown")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    nhclu_kmeans(dissimilarity = "zz"),
    "^dissimilarity is not a bioregion.pairwise object")
  
  expect_error(
    nhclu_kmeans(dissim2),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(uni[,-3]),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(unina1),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(unina2),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(uni4),
    "The first two columns of dissimilarity contain duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(uni3),
    "^The first two columns of dissimilarity contain") 
  
  expect_error(
    nhclu_kmeans(uni2),
    "dissimilarity contains rows with the same site on both columns!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(unichar),
    "The weight column must be numeric.",
    fixed = TRUE)  
  
  expect_message(
    nhclu_kmeans(d2),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(d3),
    "dissimilarity must be numeric.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(d4),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, index = "zz"),
    "^If index is a character, it should be a column name")
  
  expect_error(
    nhclu_kmeans(dissim, index = "Site1"),
    "^If index is a character, it should be a column name")
  
  expect_error(
    nhclu_kmeans(dissim, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(uni, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, seed = -1),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, seed = 0),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = "zz"),
    "n_clust must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = c(1.1,2)),
    "n_clust must be composed of integers.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = -1),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = c(1,-1)),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = 0),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, iter_max =  c("zz","zz")),
    "iter_max must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, iter_max = "zz"),
    "iter_max must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, iter_max = 1.1),
    "iter_max must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, iter_max = -1),
    "iter_max must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, algorithm = 1),
    "algorithm must be a character.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, algorithm = "zz"),
    "^Please choose algorithm from the following:")
  
  expect_error(
    nhclu_kmeans(dissim, nstart =  c("zz","zz")),
    "nstart must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, nstart = "zz"),
    "nstart must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, nstart = 1.1),
    "nstart must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, nstart = -1),
    "nstart must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
})
