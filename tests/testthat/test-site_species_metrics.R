# Tests for site_species_metrics()
# 
# This file contains unit tests for input validation, output structure, 
# and error handling.
#
# For NUMERICAL VALIDATION of calculations (comparing against manual 
# calculations), see: tests/validation/validate_site_species_metrics.R
#
# The validation tests verify that formulas are correctly implemented by
# comparing function outputs with hand-calculated expected values for
# simple test cases.
#


# Inputs -----------------------------------------------------------------------
set.seed(123)  # For reproducibility
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
                20, 25)
rownames(comat) <- paste0("Site",1:20)
colnames(comat) <- paste0("Species",1:25)

# Simple test matrix with known values for validation
comat_simple <- matrix(c(
  1, 1, 0, 0,
  1, 1, 0, 0,
  0, 0, 1, 1,
  0, 0, 1, 1
), nrow = 4, byrow = TRUE)
rownames(comat_simple) <- paste0("Site", 1:4)
colnames(comat_simple) <- paste0("Sp", 1:4)

dissim <- dissimilarity(comat, metric = "Simpson")
clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

# Create multiple bioregionalizations for testing
clust_k3 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
clust_k4 <- nhclu_kmeans(dissim, n_clust = 4, index = "Simpson")
multi_clust_manual <- clust_k3
multi_clust_manual$clusters <- cbind(clust_k3$clusters, 
                                     K_4 = clust_k4$clusters[, 2])

net <- similarity(comat, metric = "Simpson")
com <- netclu_greedy(net)

multi_clust <- nhclu_kmeans(dissim, n_clust = 3:4, index = "Simpson")

net_bip <- mat_to_net(comat, weight = TRUE)
clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)

clust_h <- hclu_hierarclust(dissim,
                            optimal_tree_method = "best",
                            n_clust = NULL,
                            cut_height = NULL,
                            verbose = FALSE)

simil <- dissimilarity_to_similarity(dissim)
clust_louv <- netclu_louvain(simil)
clust_louv$clusters <- NULL

# Tests for valid outputs ------------------------------------------------------
test_that("output structure for single bioregionalization", {
  
  rho <- site_species_metrics(bioregionalization = clust1, 
                              comat = comat,
                              indices = "rho")
  
  # Test new list structure
  expect_true(inherits(rho, "bioregion.site.species.metrics"))
  expect_true(inherits(rho, "list"))
  expect_equal(length(rho), 1)
  expect_equal(attr(rho, "n_bioregionalizations"), 1)
  expect_equal(attr(rho, "indices"), "rho")
  
  # Test structure of first element
  expect_true(is.list(rho[[1]]))
  expect_true("name" %in% names(rho[[1]]))
  expect_true("metrics" %in% names(rho[[1]]))
  expect_true("indices" %in% names(rho[[1]]))
  expect_true("args" %in% names(rho[[1]]))
  
  # Test metrics data.frame
  expect_true(inherits(rho[[1]]$metrics, "data.frame"))
  expect_equal(nrow(rho[[1]]$metrics), 75)  # 3 bioregions * 25 species
  expect_equal(ncol(rho[[1]]$metrics), 3)
  expect_true(all(c("Bioregion", "Species", "rho") %in% colnames(rho[[1]]$metrics)))
  
  # Test access by name
  expect_equal(names(rho), "K_3")
  expect_identical(rho$K_3$metrics, rho[[1]]$metrics)
})

test_that("output structure for multiple bioregionalizations", {
  
  result <- site_species_metrics(bioregionalization = multi_clust_manual, 
                                 comat = comat,
                                 indices = "rho")
  
  # Test list structure
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_true(inherits(result, "list"))
  expect_equal(length(result), 2)
  expect_equal(attr(result, "n_bioregionalizations"), 2)
  expect_equal(names(result), c("K_3", "K_4"))
  
  # Test first bioregionalization
  expect_equal(result[[1]]$name, "K_3")
  expect_equal(nrow(result[[1]]$metrics), 75)  # 3 * 25
  expect_equal(length(unique(result[[1]]$metrics$Bioregion)), 3)
  
  # Test second bioregionalization
  expect_equal(result[[2]]$name, "K_4")
  expect_equal(nrow(result[[2]]$metrics), 100)  # 4 * 25
  expect_equal(length(unique(result[[2]]$metrics$Bioregion)), 4)
  
  # Test access by name
  expect_identical(result$K_3$metrics, result[[1]]$metrics)
  expect_identical(result$K_4$metrics, result[[2]]$metrics)
})

test_that("all indices work correctly", {
  
  # Single index
  rho <- site_species_metrics(clust1, comat, indices = "rho")
  expect_true("rho" %in% colnames(rho[[1]]$metrics))
  expect_false("affinity" %in% colnames(rho[[1]]$metrics))
  
  # Multiple indices
  all_metrics <- site_species_metrics(
    clust1, comat, 
    indices = c("rho", "affinity", "fidelity", "indicator_value")
  )
  expect_true(all(c("rho", "affinity", "fidelity", "indval") %in% 
                    colnames(all_metrics[[1]]$metrics)))
  expect_equal(ncol(all_metrics[[1]]$metrics), 6)  # Bioregion, Species, + 4 metrics
})

test_that("Cz indices work for bipartite", {
  
  suppressWarnings({
    cz_only <- site_species_metrics(
      bioregionalization = clust_bip, 
      comat = comat,
      net = net_bip,
      indices = "Cz"
    )
  })
  
  expect_true(inherits(cz_only, "bioregion.site.species.metrics"))
  expect_true("cz_metrics" %in% names(cz_only[[1]]))
  expect_true(inherits(cz_only[[1]]$cz_metrics, "data.frame"))
  expect_true(all(c("Node", "Bioregion", "Category", "C", "z") %in% 
                    colnames(cz_only[[1]]$cz_metrics)))
  
  # Cz with other indices
  suppressWarnings({
    cz_mixed <- site_species_metrics(
      bioregionalization = clust_bip, 
      comat = comat,
      net = net_bip,
      indices = c("Cz", "rho")
    )
  })
  
  expect_true("metrics" %in% names(cz_mixed[[1]]))
  expect_true("cz_metrics" %in% names(cz_mixed[[1]]))
  expect_true("rho" %in% colnames(cz_mixed[[1]]$metrics))
})

test_that("print method works", {
  
  result <- site_species_metrics(clust1, comat, indices = "rho")
  
  # Should not error
  expect_output(print(result), "Site and species contribution metrics")
  expect_output(print(result), "Number of bioregionalizations")
  expect_output(print(result), "Computed indices")
  
  # Multiple bioregionalizations
  result_multi <- site_species_metrics(multi_clust_manual, comat, indices = "rho")
  expect_output(print(result_multi), "K_3, K_4")
})

test_that("str method works", {
  
  result <- site_species_metrics(clust1, comat, indices = "rho")
  
  # Should not error and show limited depth
  expect_output(str(result), "bioregion.site.species.metrics")
  expect_output(str(result), "\\$ K_3")
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    site_species_metrics("zz"),
    "^This function is designed to work on bioregion.clusters objects ")
  
  # REMOVED: No longer expect error for multiple bioregionalizations
  # This is now a valid use case!
  
  expect_error(
    site_species_metrics(clust_h),
    "^No clusters have been generated for your hierarchical")
  
  expect_error(
    site_species_metrics(clust_louv),
    "^bioregionalization does not have the expected type of")
  
  expect_error(
    site_species_metrics(com, 
                         comat = "zz"),
    "comat must be a matrix.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat,
                         indices = c(1,2)),
    "indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         indices = "zz"),
    "^Please choose indices from the following")
  
  expect_error(
    site_species_metrics(com, comat = comat, 
                         net = "zz",
                         indices = "Cz"),
    "^Cz metrics can only be computed for a bipartite")
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = NULL,
                         indices = "Cz"),
    "net is needed to compute Cz indices.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = "zz"),
    "net should be a data.frame with at least two columns,
           corresponding to the sites and species. By default, sites are
           considered to be in the first column, and species in the second.
           This can be changed with the arguments 'site_col' and
           'species_col'.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = net_bip,
                         site_col = "zz"),
    "site_col must be numeric.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(com, 
                         comat = comat, 
                         net = net_bip,
                         species_col = "zz"),
    "species_col must be numeric.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(clust_bip, 
                         comat = comat, 
                         net = net_bip,
                         indices = "Cz",
                         site_col = 30),
    "The site column ('site_col') is incorrect.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(clust_bip, 
                         comat = comat, 
                         net = net_bip,
                         indices = "Cz",
                         species_col = 30),
    "The species column ('species_col') is incorrect.",
    fixed = TRUE)
  
})

# Tests with known expected values ---------------------------------------------
test_that("rho calculation is correct with known values", {
  
  # Simple case: 4 sites, 4 species, 2 bioregions
  # Sites 1-2 in bioregion 1, Sites 3-4 in bioregion 2
  # Species 1-2 only in bioregion 1, Species 3-4 only in bioregion 2
  
  # Create simple clustering
  simple_clusters <- data.frame(
    ID = paste0("Site", 1:4),
    K_2 = c(1, 1, 2, 2)
  )
  class(simple_clusters) <- c("bioregion.clusters", "data.frame")
  
  simple_bioregion <- list(
    name = "test",
    args = list(),
    inputs = list(bipartite = FALSE),
    algorithm = list(),
    clusters = simple_clusters,
    cluster_info = list(n_clust = 2)
  )
  class(simple_bioregion) <- "bioregion.clusters"
  
  result <- site_species_metrics(
    bioregionalization = simple_bioregion,
    comat = comat_simple,
    indices = "rho"
  )
  
  # Species 1 and 2 should have perfect association with bioregion 1
  # Species 3 and 4 should have perfect association with bioregion 2
  sp1_br1 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp1" & 
                                   result[[1]]$metrics$Bioregion == 1, "rho"]
  sp1_br2 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp1" & 
                                   result[[1]]$metrics$Bioregion == 2, "rho"]
  
  # Sp1 should be strongly associated with bioregion 1
  expect_true(sp1_br1 > 0)
  # Sp1 should be negatively associated with bioregion 2
  expect_true(sp1_br2 < 0)
  
  # Check symmetry: Sp3 should behave opposite to Sp1
  sp3_br1 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp3" & 
                                   result[[1]]$metrics$Bioregion == 1, "rho"]
  sp3_br2 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp3" & 
                                   result[[1]]$metrics$Bioregion == 2, "rho"]
  
  expect_true(sp3_br1 < 0)
  expect_true(sp3_br2 > 0)
})

test_that("affinity and fidelity calculations are correct", {
  
  # Simple case with known values
  simple_clusters <- data.frame(
    ID = paste0("Site", 1:4),
    K_2 = c(1, 1, 2, 2)
  )
  class(simple_clusters) <- c("bioregion.clusters", "data.frame")
  
  simple_bioregion <- list(
    name = "test",
    args = list(),
    inputs = list(bipartite = FALSE),
    algorithm = list(),
    clusters = simple_clusters,
    cluster_info = list(n_clust = 2)
  )
  class(simple_bioregion) <- "bioregion.clusters"
  
  result <- site_species_metrics(
    bioregionalization = simple_bioregion,
    comat = comat_simple,
    indices = c("affinity", "fidelity", "indicator_value")
  )
  
  # Species 1 appears in sites 1-2 (both in bioregion 1)
  # Affinity in bioregion 1 = 2 occurrences / 2 sites in bioregion = 1.0
  sp1_br1 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp1" & 
                                   result[[1]]$metrics$Bioregion == 1, ]
  
  expect_equal(sp1_br1$affinity, 1.0)
  expect_equal(sp1_br1$fidelity, 1.0)  # All occurrences in this bioregion
  expect_equal(sp1_br1$indval, 1.0)  # Perfect indicator
  
  # Species 1 in bioregion 2 = 0 occurrences / 2 sites = 0
  sp1_br2 <- result[[1]]$metrics[result[[1]]$metrics$Species == "Sp1" & 
                                   result[[1]]$metrics$Bioregion == 2, ]
  
  expect_equal(sp1_br2$affinity, 0.0)
  expect_equal(sp1_br2$fidelity, 0.0)
  expect_equal(sp1_br2$indval, 0.0)
})

test_that("multiple bioregionalizations produce independent results", {
  
  # Create two different clustering solutions
  result <- site_species_metrics(
    bioregionalization = multi_clust_manual,
    comat = comat,
    indices = "rho"
  )
  
  # Each bioregionalization should have different number of rows
  expect_equal(nrow(result$K_3$metrics), 75)  # 3 bioregions * 25 species
  expect_equal(nrow(result$K_4$metrics), 100)  # 4 bioregions * 25 species
  
  # Check that each has correct number of unique bioregions
  expect_equal(length(unique(result$K_3$metrics$Bioregion)), 3)
  expect_equal(length(unique(result$K_4$metrics$Bioregion)), 4)
  
  # Both should have same species
  expect_equal(
    sort(unique(result$K_3$metrics$Species)),
    sort(unique(result$K_4$metrics$Species))
  )
})

test_that("all species and bioregions are present in output", {
  
  result <- site_species_metrics(clust1, comat, indices = "rho")
  
  # All species should be present
  expect_equal(length(unique(result[[1]]$metrics$Species)), ncol(comat))
  expect_setequal(unique(result[[1]]$metrics$Species), colnames(comat))
  
  # All bioregions should be present
  expect_equal(length(unique(result[[1]]$metrics$Bioregion)), 3)
  
  # Each species should appear once per bioregion
  species_counts <- table(result[[1]]$metrics$Species)
  expect_true(all(species_counts == 3))
})

test_that("metrics have valid ranges", {
  
  result <- site_species_metrics(
    clust1, comat, 
    indices = c("rho", "affinity", "fidelity", "indicator_value")
  )
  
  metrics <- result[[1]]$metrics
  
  # Affinity: proportion, should be [0, 1]
  expect_true(all(metrics$affinity >= 0 & metrics$affinity <= 1))
  
  # Fidelity: proportion, should be [0, 1]
  expect_true(all(metrics$fidelity >= 0 & metrics$fidelity <= 1))
  
  # Indicator value: product of proportions, should be [0, 1]
  expect_true(all(metrics$indval >= 0 & metrics$indval <= 1))
  
  # Rho: standardized, can be negative or positive
  # Check for finite values (no NaN, Inf)
  expect_true(all(is.finite(metrics$rho)))
})

test_that("truncated print for many bioregionalizations", {
  
  # Create many bioregionalizations
  set.seed(456)
  many_clust <- clust1
  for (k in 4:10) {
    temp <- nhclu_kmeans(dissim, n_clust = k, index = "Simpson")
    many_clust$clusters <- cbind(many_clust$clusters, 
                                 temp$clusters[, 2])
  }
  colnames(many_clust$clusters) <- c("ID", paste0("K_", 3:10))
  
  result <- site_species_metrics(many_clust, comat, indices = "rho")
  
  # Should have 8 bioregionalizations
  expect_equal(length(result), 8)
  expect_equal(attr(result, "n_bioregionalizations"), 8)
  
  # Print should mention truncation
  output <- capture.output(print(result))
  expect_true(any(grepl("Showing first 3 out of", output)))
  expect_true(any(grepl("and 5 more bioregionalization", output)))
})

test_that("consistency between single and multiple bioregionalization results", {
  
  # Compute separately
  result_k3_only <- site_species_metrics(clust_k3, comat, indices = "rho")
  result_k4_only <- site_species_metrics(clust_k4, comat, indices = "rho")
  
  # Compute together
  result_both <- site_species_metrics(multi_clust_manual, comat, indices = "rho")
  
  # Results should be identical
  expect_equal(
    result_k3_only[[1]]$metrics[order(result_k3_only[[1]]$metrics$Species,
                                       result_k3_only[[1]]$metrics$Bioregion), ],
    result_both$K_3$metrics[order(result_both$K_3$metrics$Species,
                                   result_both$K_3$metrics$Bioregion), ],
    tolerance = 1e-10
  )
  
  expect_equal(
    result_k4_only[[1]]$metrics[order(result_k4_only[[1]]$metrics$Species,
                                       result_k4_only[[1]]$metrics$Bioregion), ],
    result_both$K_4$metrics[order(result_both$K_4$metrics$Species,
                                   result_both$K_4$metrics$Bioregion), ],
    tolerance = 1e-10
  )
})

