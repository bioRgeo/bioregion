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
  
  clust <- nhclu_dbscan(dissim,
                       index = "Simpson",
                       minPts = NULL,
                       eps = NULL,
                       plot = FALSE,
                       algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "nhclu_dbscan")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$minPts, NULL)
  expect_equal(clust$args$eps, NULL)
  expect_equal(clust$args$plot, FALSE)
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, TRUE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  expect_equal(dim(clust$clusters)[2], 2)
  
  clust <- nhclu_dbscan(dissim,
                        index = "Euclidean",
                        minPts = c(3,4),
                        eps = c(0,0.1,0.2),
                        plot = FALSE)
  expect_equal(colnames(clust$clusters)[2], "K_12_1")
  expect_equal(colnames(clust$clusters)[3], "K_12_2")
  expect_equal(colnames(clust$clusters)[4], "K_12_3")
  expect_equal(colnames(clust$clusters)[5], "K_19_1")
  expect_equal(colnames(clust$clusters)[6], "K_19_2")
  expect_equal(colnames(clust$clusters)[7], "K_19_3")
  
  clust1 <- nhclu_dbscan(dissim,
                        index = "Euclidean",
                        minPts = c(3,4),
                        eps = c(0,0.1,0.2),
                        plot = FALSE)
  clust2 <- nhclu_dbscan(dissim,
                         index = "Euclidean",
                         minPts = c(3,4),
                         eps = c(0,0.1,0.2),
                         plot = FALSE)
  expect_equal(sum(clust1$clusters$K_19_3==clust2$clusters$K_19_3, 
                   na.rm = TRUE), 136)
  
  # Test data_type with different dissimilarity metrics
  clust <- nhclu_dbscan(dissim, index = "Simpson", minPts = 3, eps = 0.1)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_dbscan(dissim, index = "Jaccard", minPts = 3, eps = 0.1)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_dbscan(dissim, index = "Bray", minPts = 3, eps = 0.1)
  expect_equal(clust$inputs$data_type, "abundance")
  
  clust <- nhclu_dbscan(dissim, index = "Euclidean", minPts = 3, eps = 0.1)
  expect_equal(clust$inputs$data_type, "unknown")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    nhclu_dbscan(dissimilarity = "zz"),
    "^dissimilarity is not a bioregion.pairwise object")
  
  expect_error(
    nhclu_dbscan(dissim2),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(uni[,-3]),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(unina1),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_dbscan(unina2),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_dbscan(uni4),
    "The first two columns of dissimilarity contain duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(uni3),
    "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(uni2),
    "dissimilarity contains rows with the same site on both columns!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(unichar),
    "The weight column must be numeric.",
    fixed = TRUE)  
  
  expect_message(
    nhclu_dbscan(d2),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(d3),
    "dissimilarity must be numeric.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(d4),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(dissim, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, index = "zz"),
    "^If index is a character, it should be a column name")
  
  expect_error(
    nhclu_dbscan(dissim, index = "Site1"),
    "^If index is a character, it should be a column name")
  
  expect_error(
    nhclu_dbscan(dissim, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(uni, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, minPts = "zz"),
    "minPts must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_dbscan(dissim, minPts = c(1.1,2)),
    "minPts must be composed of integers.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_dbscan(dissim, minPts = -1),
    "minPts must be strictly higher than 2.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(dissim, minPts = 2),
    "minPts must be strictly higher than 2.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(dissim, eps = "zz"),
    "eps must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_dbscan(dissim, eps = -1),
    "eps must be composed of values higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_dbscan(dissim, plot = 1),
    "plot must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, plot = c(TRUE,FALSE)),
    "plot must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_dbscan(dissim, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
})
