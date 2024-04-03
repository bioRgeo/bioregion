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
                       n_clust = c(1,2,3),
                       iter_max = 10,
                       nstart = 10,
                       algorithm = "Hartigan-Wong",
                       algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "nhclu_kmeans")
  expect_equal(clust$args$index, "Simpson")
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
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    nhclu_kmeans(dissimilarity = "zz"),
    "dissimilarity is not a bioregion.pairwise.metric object, 
a dissimilarity matrix (class dist) or 
a data.frame with at least 3 columns (site1, site2 and your dissimilarity index).",
    fixed = TRUE)
  
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
    "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
    fixed = TRUE) 
  
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
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, index = "Site1"),
    "If index is a character, it should be a column name (and not the
                    first or second column).",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, index = 2),
    "index should be stricltly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(uni, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = "zz"),
    "n_clust must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = c(1.1,2)),
    "n_clust must be composed of integer(s).",
    fixed = TRUE)  
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = -1),
    "n_clust must be composed of value(s) strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = c(1,-1)),
    "n_clust must be composed of value(s) strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_kmeans(dissim, n_clust = 0),
    "n_clust must be composed of value(s) strictly higher than 0.",
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
    "Please choose algorithm among the followings values:
Hartigan-Wong, Lloyd, Forgy or MacQueen",
    fixed = TRUE)
  
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
