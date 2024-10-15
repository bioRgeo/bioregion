# Inputs -----------------------------------------------------------------------
dissim <- dissimilarity(comat, metric = "Simpson")
sim <- dissimilarity_to_similarity(dissim)

df <- data.frame(ID1 = sim$Site1, ID2 = sim$Site2, W = sim$Simpson)

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
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         seed = NULL,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "nhclu_affprop")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$p, NA)
  expect_equal(clust$args$q, NA)
  expect_equal(clust$args$maxits, 1000)
  expect_equal(clust$args$convits, 100)
  expect_equal(clust$args$lam, 0.9)
  expect_equal(clust$args$details, FALSE)
  expect_equal(clust$args$nonoise, FALSE)
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, FALSE)
  expect_equal(clust$inputs$nb_sites, 20)
  expect_equal(clust$inputs$hierarchical, FALSE)
  # expect_equal(dim(clust$clusters)[2], 4)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    nhclu_affprop("zz"),
    "similarity is not a bioregion.pairwise.metric object, 
a dissimilarity matrix (class dist) or 
a data.frame with at least 3 columns (site1, site2 and your dissimilarity index).",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(dissim2),
    "similarity seems to be a dissimilarity object. 
This function should be applied on similarities, not dissimilarities. 
Use dissimilarity_to_similarity() before using this function.",
    fixed = TRUE)
  
})
