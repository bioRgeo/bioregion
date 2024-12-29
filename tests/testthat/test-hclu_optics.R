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
  
  clust <- hclu_optics(dissim,
                       index = "Simpson",
                       minPts = NULL,
                       eps = NULL,
                       xi = 0.05,
                       minimum = FALSE,
                       show_hierarchy = FALSE,
                       algorithm_in_output = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "hclu_optics")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$minPts, NULL)
  expect_equal(clust$args$eps, NULL)
  expect_equal(clust$args$xi, 0.05)
  expect_equal(clust$args$minimum, FALSE)
  expect_equal(clust$args$show_hierarchy, FALSE)
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, TRUE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(dim(clust$clusters)[2], 2)
  
  clust1 <- hclu_optics(dissim,
                       index = "Simpson",
                       minPts = NULL,
                       eps = NULL,
                       xi = 0.05,
                       minimum = FALSE,
                       show_hierarchy = FALSE,
                       algorithm_in_output = TRUE)
  clust2 <- hclu_optics(dissim,
                       index = "Simpson",
                       minPts = NULL,
                       eps = NULL,
                       xi = 0.05,
                       minimum = FALSE,
                       show_hierarchy = FALSE,
                       algorithm_in_output = TRUE)
  expect_equal(sum(clust1$clusters$K_3==clust2$clusters$K_3, na.rm = TRUE), 337)
  
  clust1 <- hclu_optics(dissim,
                        index = 5,
                        show_hierarchy = FALSE)
  expect_equal(clust$inputs$hierarchical, FALSE)
  
  clust2 <- hclu_optics(dissim,
                        index = 5,
                        show_hierarchy = TRUE)
  expect_equal(clust2$inputs$hierarchical, TRUE)
  tab12 <- table(clust1$clusters$K_18,clust2$clusters$K_18)
  expect_equal(sum(apply(tab12==0,1,sum)==17),18)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    hclu_optics(dissimilarity = "zz"),
    "^dissimilarity is not a bioregion.pairwise.metric object")
  
  expect_error(
    hclu_optics(dissim2),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(uni[,-3]),
    "dissimilarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(unina1),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(unina2),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(uni4),
    "The first two columns of dissimilarity contain duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(uni3),
    "The first two columns of dissimilarity contain (unordered) duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(uni2),
    "dissimilarity contains rows with the same site on both columns!",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(unichar),
    "The weight column must be numeric.",
    fixed = TRUE)  
  
  expect_message(
    hclu_optics(d2),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(d3),
    "dissimilarity must be numeric.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(d4),
    "NA(s) detected in dissimilarity.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, index = c("zz",1)),
    "index must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, index = "zz"),
    "^If index is a character, it should be ")
  
  expect_error(
    hclu_optics(dissim, index = "Site1"),
    "^If index is a character, it should be ")
  
  expect_error(
    hclu_optics(dissim, index = 0.1),
    "If index is numeric, it should be an integer.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, index = 2),
    "index should be strictly higher than 2.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(uni, index = 4),
    "index should be lower or equal to 3.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, minPts =  c("zz","zz")),
    "minPts must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, minPts = "zz"),
    "minPts must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, minPts = 1.1),
    "minPts must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, minPts = -1),
    "minPts must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, eps =  c("zz","zz")),
    "eps must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, eps = "zz"),
    "eps must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, eps = 1.1),
    "eps must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, eps = -1),
    "eps must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, xi =  c("zz","zz")),
    "xi must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, xi = "zz"),
    "xi must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    hclu_optics(dissim, xi = -1),
    "xi must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, xi = 0),
    "xi must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, xi = 1),
    "xi must be in the interval (0,1), (see dbscan::optics())",
    fixed = TRUE) 
  
  expect_error(
    hclu_optics(dissim, minimum = 1),
    "minimum must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, minimum = c(TRUE,FALSE)),
    "minimum must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, minimum = 1),
    "minimum must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, show_hierarchy = c(TRUE,FALSE)),
    "show_hierarchy must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, show_hierarchy = 1),
    "show_hierarchy must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    hclu_optics(dissim, algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
})
