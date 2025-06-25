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

  clust <- hclu_hierarclust(dissim,
                       index = "Simpson",
                       method = "average",
                       randomize = TRUE,
                       n_runs = 30,
                       keep_trials = FALSE,
                       optimal_tree_method = "best",
                       n_clust = c(1,2,3),
                       cut_height = NULL,
                       find_h = TRUE,
                       h_max = 1,
                       h_min = 0)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "hclu_hierarclust")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$method, "average")
  expect_equal(clust$args$randomize, TRUE)
  expect_equal(clust$args$n_runs, 30)
  expect_equal(clust$args$keep_trials, FALSE)
  expect_equal(clust$args$optimal_tree_method, "best")
  expect_equal(clust$args$n_clust, c(1,2,3))
  expect_equal(clust$args$cut_height, NULL)
  expect_equal(clust$args$find_h, TRUE)
  expect_equal(clust$args$h_max, 1)
  expect_equal(clust$args$h_min, 0)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, TRUE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, TRUE)
  expect_equal(dim(clust$clusters)[2], 4)
  
  clust <- hclu_hierarclust(d,
                            optimal_tree_method = "best",)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "hclu_hierarclust")

  clust <- hclu_hierarclust(dissim,
                            optimal_tree_method = "best",
                            index = "Euclidean",
                            n_clust = c(5,10),
                            find_h = FALSE)
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")

  clust <- hclu_hierarclust(dissim,
                            optimal_tree_method = "best",
                            index = "Euclidean",
                            n_clust = c(10,5),
                            find_h = FALSE)
  expect_equal(colnames(clust$clusters)[2], "K_5")
  expect_equal(colnames(clust$clusters)[3], "K_10")

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    hclu_hierarclust(dissimilarity = "zz",
                     optimal_tree_method = "best"),
    "^dissimilarity is not a bioregion.pairwise object")
  
  expect_error(
    hclu_hierarclust(dissim2,
                     optimal_tree_method = "best"),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(uni[,-3],
                     optimal_tree_method = "best"),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(unina1,
                     optimal_tree_method = "best"),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(unina2,
                     optimal_tree_method = "best"),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(uni4,
                     optimal_tree_method = "best"),
    "The first two columns of dissimilarity contain duplicated pairs of sites!",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(uni3,
                     optimal_tree_method = "best"),
    "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(uni2,
                     optimal_tree_method = "best"),
    "dissimilarity contains rows with the same site on both columns!",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(unichar,
                     optimal_tree_method = "best"),
    "The weight column must be numeric.",
    fixed = TRUE)

  expect_message(
    hclu_hierarclust(d2,
                     optimal_tree_method = "best"),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(d3,
                     optimal_tree_method = "best"),
    "dissimilarity must be numeric.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(d4,
                     optimal_tree_method = "best"),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     index = c("zz",1),
                     optimal_tree_method = "best"),
    "index must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     index = "zz",
                     optimal_tree_method = "best"),
    "^If index is a character, it should be ")

  expect_error(
    hclu_hierarclust(dissim, 
                     index = "Site1",
                     optimal_tree_method = "best"),
    "^If index is a character, it should be ")

  expect_error(
    hclu_hierarclust(dissim, 
                     index = 0.1,
                     optimal_tree_method = "best"),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     index = 2,
                     optimal_tree_method = "best"),
    "index should be strictly higher than 2.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(uni, 
                     index = 4,
                     optimal_tree_method = "best"),
    "index should be lower or equal to 3.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     method = c(1,"z"),
                     optimal_tree_method = "best"),
    "method must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     method = 1,
                     optimal_tree_method = "best"),
    "method must be a character.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     method = "zz",
                     optimal_tree_method = "best"),
    "^Please choose method from")

  expect_error(
    hclu_hierarclust(dissim, 
                     randomize = 1,
                     optimal_tree_method = "best"),
    "randomize must be a boolean.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     randomize = c("zz","zz"),
                     optimal_tree_method = "best"),
    "randomize must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_runs =  c("zz","zz"),
                     optimal_tree_method = "best"),
    "n_runs must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, ,
                     optimal_tree_method = "best",
                     n_runs = "zz"),
    "n_runs must be numeric.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_runs = 1.1,
                     optimal_tree_method = "best"),
    "n_runs must be an integer.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_runs = -1,
                     optimal_tree_method = "best"),
    "n_runs must be strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_runs = 0,
                     optimal_tree_method = "best"),
    "n_runs must be strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     keep_trials = 1,
                     optimal_tree_method = "best"),
    "keep_trials must be a boolean.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     keep_trials = c("zz","zz"),
                     optimal_tree_method = "best"),
    "keep_trials must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, optimal_tree_method = 1),
    "optimal_tree_method must be a character.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, optimal_tree_method = c("1","1")),
    "optimal_tree_method must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, optimal_tree_method = "zz"),
    "^Please choose optimal_tree_method from the following")

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = "zz",
                     optimal_tree_method = "best"),
    "^n_clust must be one of those:")

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = c(1.1,2),
                     optimal_tree_method = "best"),
    "n_clust must be composed of integers.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = -1,
                     optimal_tree_method = "best"),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = c(1,-1),
                     optimal_tree_method = "best"),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = 0,
                     optimal_tree_method = "best"),
    "n_clust must be composed of values strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = 1, 
                     cut_height = 1,
                     optimal_tree_method = "best"),
    "^Please provide either n_clust or cut_height")

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = NULL, 
                     cut_height = "zz",
                     optimal_tree_method = "best"),
    "cut_height must be numeric.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     n_clust = NULL, 
                     cut_height = -1,
                     optimal_tree_method = "best"),
    "cut_height must be composed of values higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     find_h = 1,
                     optimal_tree_method = "best"),
    "find_h must be a boolean.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     find_h = c("zz","zz"),
                     optimal_tree_method = "best"),
    "find_h must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_min =  c("zz","zz"),
                     optimal_tree_method = "best"),
    "h_min must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_min = "zz",
                     optimal_tree_method = "best"),
    "h_min must be numeric.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_min = -1,
                     optimal_tree_method = "best"),
    "h_min must be higher than 0.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_max =  c("zz","zz"),
                     optimal_tree_method = "best"),
    "h_max must be of length 1.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_max = "zz",
                     optimal_tree_method = "best"),
    "h_max must be numeric.",
    fixed = TRUE)

  expect_error(
    hclu_hierarclust(dissim, 
                     h_max = -1,
                     optimal_tree_method = "best"),
    "h_max must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     h_min = 1,
                     h_max = 0,
                     optimal_tree_method = "best"),
    "h_min must be inferior to h_max.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     consensus_p =  c("zz","zz"),
                     optimal_tree_method = "best"),
    "consensus_p must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     consensus_p = "zz",
                     optimal_tree_method = "best"),
    "consensus_p must be numeric.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     consensus_p = -1,
                     optimal_tree_method = "best"),
    "consensus_p must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     consensus_p = 0.4,
                     optimal_tree_method = "best"),
    "consensus_p must be between 0.5 and 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_hierarclust(dissim, 
                     consensus_p = 1.1,
                     optimal_tree_method = "best"),
    "consensus_p must be between 0.5 and 1.",
    fixed = TRUE)
  
  expect_message(
    hclu_hierarclust(dissim, 
                     optimal_tree_method = "iterative_consensus_tree"),
    "^Building the iterative")
  
  expect_warning(
    hclu_hierarclust(dissim, 
                     method = "mcquitty",
                     optimal_tree_method = "iterative_consensus_tree"),
    "^mcquitty")
  
  expect_error(
    hclu_hierarclust(dissim, 
                     optimal_tree_method = "consensus",
                     n_runs = 1),
    "At least two trees are required to calculate a consensus.",
    fixed = TRUE)

})
