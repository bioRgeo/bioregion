# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001),
                5, 10)
rownames(comat) <- paste0("s", 1:5)
colnames(comat) <- paste0("sp", 1:10)
comat1 <- matrix(sample(0:1000, size = 100, replace = TRUE, prob = 1 / 1:1001),
                10, 10)
rownames(comat1) <- paste0("s", 1:10)
colnames(comat1) <- paste0("sp", 1:10)
comat2 <- matrix(sample(0:1000, size = 100, replace = TRUE, prob = 1 / 1:1001),
                5, 100)
rownames(comat2) <- paste0("s", 1:5)
colnames(comat2) <- paste0("sp", 1:100)

comat3 <- matrix(sample(0:1000, size = 100, replace = TRUE, prob = 1 / 1:1001),
                 5, 100)
rownames(comat3) <- paste0("x", 1:5)
colnames(comat3) <- paste0("sp", 1:100)

simil0 <- similarity(comat, metric = "all", formula = "a + b")
attr(simil0, "type") <- NULL

simil00 <- similarity(comat, metric = "all", formula = "a + b")
attr(simil00, "type") <- "umstinetnu"

simil1 <- similarity(comat, metric = "all", formula = "a + b")

dissimil1 <- dissimilarity(comat, metric = "all", formula = "a + b")

simil2 <- similarity(comat1, metric = "all", formula = "a + b")

simil3 <- similarity(comat2, metric = "all", formula = "a + b")

simil4 <- similarity(comat3, metric = "all", formula = "a + b")


# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  pair <- combine_bioregion_pairwise(primary_metrics = simil1, 
                                     secondary_metrics = simil3,
                                     new_metrics = NULL)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 5)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[1], 10)
  expect_equal(dim(pair)[2], 30)
  expect_equal(pair[1,1], "s1")
  
  pair <- combine_bioregion_pairwise(primary_metrics = similarity(comat, metric = "abc"), 
                                     secondary_metrics = NULL,
                                     new_metrics = "a + b")
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 5)
  expect_equal(attr(pair, "nb_species"), 10)
  expect_equal(dim(pair)[1], 10)
  expect_equal(dim(pair)[2], 6)
  expect_equal(pair[1,1], "s1")
  
  pair <- combine_bioregion_pairwise(similarity(comat, metric = "abc"), 
                                     similarity(comat, metric = "Simpson"),
                                     new_metrics = c("a + b", "Simpson*Simpson"))
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 5)
  expect_equal(attr(pair, "nb_species"), 10)
  expect_equal(dim(pair)[1], 10)
  expect_equal(dim(pair)[2], 8)
  expect_equal(pair[1,1], "s1")

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = 1,
                               secondary_metrics = 1,
                               new_metrics = NULL),
    "^primary_metrics should be a bioregion.pairwise.metric object")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil0,
                               secondary_metrics = 1,
                               new_metrics = NULL),
    "^primary_metrics is a bioregion.pairwise.metric object but it has not")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil00,
                               secondary_metrics = 1,
                               new_metrics = NULL),
    "^primary_metrics is a bioregion.pairwise.metric object but it has not")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = 1,
                               new_metrics = NULL),
    "^secondary_metrics should be a bioregion.pairwise.metric object")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil0,
                               new_metrics = NULL),
    "^secondary_metrics is a bioregion.pairwise.metric object but it has not")
  
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil00,
                               new_metrics = NULL),
    "^secondary_metrics is a bioregion.pairwise.metric object but it has not")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = dissimil1,
                               new_metrics = NULL),
    "^primary_metrics and secondary_metrics should have the same type")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil2,
                               new_metrics = NULL),
    "^primary_metrics and secondary_metrics should have the same number")
  
  expect_message(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil3,
                               new_metrics = NULL),
    "^primary_metrics and secondary_metrics are based")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil4,
                               new_metrics = NULL),
    "^primary_metrics and secondary_metrics should have the same sites")
  
  expect_error(
    combine_bioregion_pairwise(primary_metrics = simil1,
                               secondary_metrics = simil3,
                               new_metrics = 1),
    "new_metrics must be a character.",
    fixed = TRUE)
  

})