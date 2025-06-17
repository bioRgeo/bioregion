# Inputs -----------------------------------------------------------------------
mat0 <- matrix(runif(100), 10, 10)
mat1 <- mat0
mat1[1,1] <- NA
mat2 <- mat0
rownames(mat2) <- paste0("s",1:10)

nbsite <- 10
nbsp <- 10
set.seed(1)
comat <- matrix(runif(nbsite*nbsp), nbsite, nbsp)
rownames(comat) <- paste0("s", 1:nbsite)
colnames(comat) <- paste0("sp", 1:nbsp)
comatbin <- comat
comatbin[comat > 0.7] <- 1
comatbin[comat <= 0.7] <- 0

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  pair <- as_bioregion_pairwise(list(mat0,mat0,mat0), 
                                metric_name = c("abc", "ABC", "Euclidean"),
                                pkg = NULL,
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  
  pair <- as_bioregion_pairwise(list(mat0,mat0,mat0), 
                                metric_name = c("abc", "ABC", "Euclidean"),
                                pkg = NULL,
                                is_similarity = TRUE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  
  pair <- as_bioregion_pairwise(list(mat0=mat0,mat1=mat0,mat2=mat0), 
                                metric_name = NULL,
                                pkg = NULL,
                                is_similarity = TRUE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(colnames(pair)[3], "mat0")
  
  pair <- as_bioregion_pairwise(list(mat2,mat0,mat0), 
                                metric_name = NULL,
                                pkg = NULL,
                                is_similarity = TRUE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(mat0, 
                                metric_name = NULL,
                                pkg = NULL,
                                is_similarity = TRUE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(colnames(pair)[3], "Metric")
  
  pair <- as_bioregion_pairwise(mat0, 
                                metric_name = "Mettric",
                                pkg = NULL,
                                is_similarity = TRUE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "similarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(colnames(pair)[3], "Mettric")
  

})

# Tests for external packages --------------------------------------------------
test_that("adespatial", {
  skip_if_not_installed("adespatial")
  library(adespatial)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = "adespatial",
                          is_similarity = FALSE),
    "mat does not seem to be an output from the adespatial package.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(beta.div(comatbin),
                          metric_name = NULL,
                          pkg = "adespatial",
                          is_similarity = FALSE),
    "D is NULL. Check that save.D=TRUE.",
    fixed = TRUE)
  
  pair <- as_bioregion_pairwise(beta.div(comatbin,
                                         save.D = TRUE), 
                                metric_name = NULL,
                                pkg = "adespatial",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(colnames(pair)[3], "hellinger")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(beta.div.comp(comatbin), 
                                metric_name = NULL,
                                pkg = "adespatial",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(colnames(pair)[3], "Podani family, Jaccard")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(beta.div.comp(comatbin,
                                              save.abc = TRUE), 
                                metric_name = NULL,
                                pkg = "adespatial",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 6)
  expect_equal(colnames(pair)[3], "Podani family, Jaccard")
  expect_equal(colnames(pair)[4], "a")
  expect_equal(pair[1,1], "s1")

})

test_that("betapart", {
  skip_if_not_installed("betapart")
  library(betapart)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = "betapart",
                          is_similarity = FALSE),
    "mat does not seem to be an output from the betapart package.",
    fixed = TRUE)
  
  pair <- as_bioregion_pairwise(beta.pair(comatbin), 
                                metric_name = NULL,
                                pkg = "betapart",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 5)
  expect_equal(colnames(pair)[3], "beta.sim")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(beta.pair(comatbin,
                                          index.family = "jaccard"), 
                                metric_name = NULL,
                                pkg = "betapart",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 5)
  expect_equal(colnames(pair)[3], "beta.jtu")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(beta.pair.abund(comat), 
                                metric_name = NULL,
                                pkg = "betapart",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 5)
  expect_equal(colnames(pair)[3], "beta.bray.bal")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(betapart.core(comatbin), 
                                metric_name = NULL,
                                pkg = "betapart",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 5)
  expect_equal(colnames(pair)[3], "a")
  expect_equal(pair[1,1], "s1")
  
})

test_that("ecodist", {
  skip_if_not_installed("ecodist")
  library(ecodist)
  registerS3method("dim", "dist", get("dim.dist", envir = asNamespace("proxy")))
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = "ecodist",
                          is_similarity = FALSE),
    "mat does not seem to be an output from the ecodist package.",
    fixed = TRUE)
  
  pair <- as_bioregion_pairwise(distance(comatbin,
                                         method = "jaccard"), 
                                metric_name = NULL,
                                pkg = "ecodist",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 3)
  expect_equal(colnames(pair)[3], "jaccard")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(bcdist(comat), 
                                metric_name = NULL,
                                pkg = "ecodist",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 3)
  expect_equal(colnames(pair)[3], "bray-curtis")
  expect_equal(pair[1,1], "s1")

})

test_that("vegan", {
  skip_if_not_installed("vegan")
  library(vegan)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = "vegan",
                          is_similarity = FALSE),
    "mat does not seem to be an output from the vegan package.",
    fixed = TRUE)
  
  pair <- as_bioregion_pairwise(vegdist(comatbin,
                                        method = "jaccard"), 
                                metric_name = NULL,
                                pkg = "vegan",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 3)
  expect_equal(colnames(pair)[3], "jaccard")
  expect_equal(pair[1,1], "s1")
  
  pair <- as_bioregion_pairwise(designdist(comatbin,
                                           method = "(A+B-2*J)/(A+B-J)",
                                           terms = "binary"), 
                                metric_name = NULL,
                                pkg = "vegan",
                                is_similarity = FALSE)
  expect_equal(inherits(pair, "bioregion.pairwise.metric"), TRUE)
  expect_equal(attr(pair, "type"), "dissimilarity")
  expect_equal(attr(pair, "nb_sites"), 10)
  expect_equal(attr(pair, "nb_species"), NA)
  expect_equal(dim(pair)[2], 3)
  expect_equal(colnames(pair)[3], "binary (A+B-2*J)/(A+B-J)")
  expect_equal(pair[1,1], "s1")
  
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = 1,
                          pkg = NULL,
                          is_similarity = FALSE),
    "metric_name must be a character.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = c("i","i"),
                          pkg = NULL,
                          is_similarity = FALSE),
    "metric_name should have the same length as mat.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,mat0),
                          metric_name = c("i","i"),
                          pkg = NULL,
                          is_similarity = FALSE),
    "metric_name should have the same length as mat.",
    fixed = TRUE)
  
  expect_message(
    expect_error(    
      as_bioregion_pairwise(mat0,
                            metric_name = "1",
                            pkg = "nistn",
                            is_similarity = FALSE),
      "^Please choose pkg from the following:"),
    "metric_name will be ignored when pkg is not NULL.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = 1),
    "is_similarity must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(mat0,
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = c(FALSE,FALSE)),
    "is_similarity must be of length 1.",
    fixed = TRUE)
  
  expect_message(
    expect_error(    
      as_bioregion_pairwise(mat0,
                            metric_name = NULL,
                            pkg = "nistn",
                            is_similarity = TRUE),
      "^Please choose pkg from the following:"),
    "is_similarity will be ignored when pkg is not NULL.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(1,
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "mat must be a matrix, a dist object, or a list of these.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(list(mat0,1,mat0),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "mat must be a matrix, a dist object, or a list of these.",
    fixed = TRUE)
  
  expect_error(
    as_bioregion_pairwise(mat1,
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat1,mat1),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,matrix(1,1,1)),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,matrix("1",10,10)),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,matrix(1,1,1)),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,matrix(1,10,100)),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "^mat must be or contain only numeric")
  
  expect_error(
    as_bioregion_pairwise(list(mat0,mat0,matrix(1,3,3)),
                          metric_name = NULL,
                          pkg = NULL,
                          is_similarity = FALSE),
    "mat must contain only square matrices with the same number sites.",
    fixed = TRUE)
  
})