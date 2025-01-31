# Inputs -----------------------------------------------------------------------
comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
                         prob = 1/1:1001), 10, 12)
rownames(comat_1) <- paste0("Site", 1:10)
colnames(comat_1) <- paste0("Species", 1:12)
comat_1 <- cbind(comat_1,
                 matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 1:10),
                                        paste0("Species", 13:20))))

comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
                         prob = 1/1:1001), 10, 12)
rownames(comat_2) <- paste0("Site", 11:20)
colnames(comat_2) <- paste0("Species", 9:20)
comat_2 <- cbind(matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 11:20),
                                        paste0("Species", 1:8))),
                 comat_2)

comat <- rbind(comat_1, comat_2)

dissim <- dissimilarity(comat, metric = "Simpson")
clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

net <- similarity(comat, metric = "Simpson")
com <- netclu_greedy(net)

multi_clust <- nhclu_kmeans(dissim, n_clust = 3:4, index = "Simpson")

clust_h <- hclu_hierarclust(dissim,
                            optimal_tree_method = "best",
                            n_clust = NULL,
                            cut_height = NULL)

simil <- dissimilarity_to_similarity(dissim)
clust_louv <- netclu_louvain(simil)
clust_louv$clusters <- NULL


# Spatial coherence
vegedissim <- dissimilarity(vegemat)
hclu <- nhclu_kmeans(dissimilarity = vegedissim, n_clust = 4)
vegemap <- map_bioregions(hclu, 
                          vegesf, 
                          write_clusters = TRUE, 
                          plot = FALSE)

vegemap_error <- vegemap
vegemap_error$K_4 <- NA

vegenodf <- vegemap
class(vegenodf) <- "sf"

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {

  #test_output <- bioregion_metrics(bioregionalization = clust1, 
  #                                 comat = comat)
  #expect_equal(inherits(test_output, "data.frame"), TRUE)
  #expect_equal(dim(test_output)[1], 3)
  #expect_equal(dim(test_output)[2], 5)

  test_output2 <- bioregion_metrics(bioregionalization = hclu, 
                                    comat = vegemat,
                                    map = vegemap, 
                                    col_bioregion = 2)
  expect_equal(inherits(test_output2, "data.frame"), TRUE)
  expect_equal(dim(test_output2)[1], 4)
  expect_equal(dim(test_output2)[2], 6)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {

  expect_error(
    bioregion_metrics("zz"),
    "^This function is designed to work on bioregion.clusters objects ")

  expect_error(
    bioregion_metrics(multi_clust, comat = comat),
    "^This function is designed to be applied on a single")
  
  expect_error(
    bioregion_metrics(clust_h),
    "^No clusters have been generated for your hierarchical")
  
  expect_error(
    bioregion_metrics(clust_louv),
    "^bioregionalization does not have the expected type of")

  expect_error(
    bioregion_metrics(com, 
                      comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = "zz"),
    "map must be a 'sf' spatial data.frame with bioregions and sites.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegenodf),
    "map must be a 'sf' spatial data.frame with bioregions and sites.",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap[, 2, drop = FALSE],
                      col_bioregion = 2),
    "map must have at least 3 columns: sites, bioregions and geometry.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap,
                      col_bioregion = NULL),
    "^col_bioregion must be defined")
  
  expect_warning(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = NULL,
                      col_bioregion = 1),
    "^col_bioregion is defined but is not considered since")
  
  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap,
                      col_bioregion = c("zz",1)),
    "col_bioregion must be of length 1.",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap,
                      col_bioregion = "zz"),
    "col_bioregion must be numeric.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap,
                      col_bioregion = 0.1),
    "col_bioregion must be an integer.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(hclu, 
                      comat = vegemat, 
                      map = vegemap,
                      col_bioregion = 0),
    "col_bioregion must be strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(hclu, comat = vegemat, map = vegemap_error,
                      col_bioregion = 2),
    "There is no bioregion in the Bioregion column.",
    fixed = TRUE)

})
