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
  

  
  # Test that net is required for Cz when no comat provided
  expect_error(
    site_species_metrics(com, 
                         net = NULL,
                         indices = "Cz"),
    "Either 'comat' or 'net' must be provided.",
    fixed = TRUE)
  
  # Test invalid net format (when comat not provided)
  expect_error(
    site_species_metrics(com, 
                         net = "zz"),
    "net should be a data.frame with at least two columns",
    fixed = TRUE)
  
  # Test site_col validation when using net (no comat)
  # NOTE: site_col can now be character or numeric, so test invalid character name
  expect_error(
    site_species_metrics(com, 
                         net = net_bip,
                         site_col = "NonExistentColumn"),
    "site_col 'NonExistentColumn' not found in net column names.",
    fixed = TRUE)
  
  # Test species_col validation when using net (no comat)
  expect_error(
    site_species_metrics(com, 
                         net = net_bip,
                         species_col = "NonExistentColumn"),
    "species_col 'NonExistentColumn' not found in net column names.",
    fixed = TRUE)
  
  # Test site_col index out of bounds
  expect_error(
    site_species_metrics(clust_bip, 
                         net = net_bip,
                         indices = "Cz",
                         site_col = 30),
    "The site column ('site_col') is incorrect.",
    fixed = TRUE)
  
  # Test species_col index out of bounds
  expect_error(
    site_species_metrics(clust_bip, 
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

# Tests for Feature 1: Interchangeable comat/net -------------------------------
test_that("only comat provided works (current behavior)", {
  
  result <- site_species_metrics(
    bioregionalization = clust1,
    comat = comat,
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)  # 3 bioregions * 25 species
})

test_that("only net provided works (new behavior)", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  
  result <- site_species_metrics(
    bioregionalization = clust1,
    net = net_test,
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)
})

test_that("both comat and net provided uses comat with warning", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  
  expect_warning(
    result <- site_species_metrics(
      bioregionalization = clust1,
      comat = comat,
      net = net_test,
      indices = "rho"
    ),
    "Both 'comat' and 'net' provided. Using 'comat' and ignoring 'net'."
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
})

test_that("neither comat nor net provided errors", {
  
  expect_error(
    site_species_metrics(
      bioregionalization = clust1,
      indices = "rho"
    ),
    "Either 'comat' or 'net' must be provided.",
    fixed = TRUE
  )
})

test_that("weight parameter affects Cz indices but not other metrics", {
  
  # Create test data for unipartite network
  comat_test <- matrix(sample(0:10, 100, replace = TRUE), 10, 10)
  rownames(comat_test) <- paste0("Site", 1:10)
  colnames(comat_test) <- paste0("Species", 1:10)
  
  dissim_test <- dissimilarity(comat_test, metric = "Simpson")
  clust_test <- nhclu_kmeans(dissim_test, n_clust = 3)
  
  net_weighted <- mat_to_net(comat_test, weight = TRUE)
  net_unweighted <- mat_to_net(comat_test, weight = FALSE)
  
  # Test with non-Cz indices - weight should not matter
  result_rho_weighted <- site_species_metrics(
    bioregionalization = clust_test,
    net = net_weighted,
    weight = TRUE,
    indices = "rho"
  )
  
  result_rho_unweighted <- site_species_metrics(
    bioregionalization = clust_test,
    net = net_unweighted,
    weight = FALSE,
    indices = "rho"
  )
  
  # Rho should be identical regardless of weight
  expect_equal(result_rho_weighted[[1]]$metrics, 
               result_rho_unweighted[[1]]$metrics)
  
  # Test with Cz indices - weight SHOULD matter
  # NOTE: Must provide ONLY net (not comat), otherwise comat takes precedence
  # When net is provided, the function auto-detects if it's weighted (3+ columns)
  # and uses the weights appropriately in Cz calculations
  suppressWarnings({
    result_cz_weighted <- site_species_metrics(
      bioregionalization = clust_test,
      net = net_weighted,
      indices = "Cz"
    )
    
    result_cz_unweighted <- site_species_metrics(
      bioregionalization = clust_test,
      net = net_unweighted,
      indices = "Cz"
    )
  })
  
  # Cz metrics should differ when weights differ
  expect_false(identical(result_cz_weighted[[1]]$cz_metrics$C, 
                        result_cz_unweighted[[1]]$cz_metrics$C))
})

test_that("reversed site_col/species_col order", {
  
  # Create net with reversed columns (species, site, weight)
  net_reversed <- mat_to_net(comat, weight = TRUE)
  net_reversed <- net_reversed[, c(2, 1, 3)]
  colnames(net_reversed) <- c("Species", "Site", "Weight")
  
  result <- site_species_metrics(
    bioregionalization = clust1,
    net = net_reversed,
    site_col = 2,
    species_col = 1,
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)
})

test_that("results equivalent when using comat vs equivalent net", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  
  result_comat <- site_species_metrics(
    bioregionalization = clust1,
    comat = comat,
    indices = "rho"
  )
  
  result_net <- site_species_metrics(
    bioregionalization = clust1,
    net = net_test,
    indices = "rho"
  )
  
  # Sort both for comparison
  metrics_comat <- result_comat[[1]]$metrics
  metrics_comat <- metrics_comat[order(metrics_comat$Species, 
                                       metrics_comat$Bioregion), ]
  
  metrics_net <- result_net[[1]]$metrics
  metrics_net <- metrics_net[order(metrics_net$Species, 
                                   metrics_net$Bioregion), ]
  
  # Remove rownames to avoid error in expect_equal
  rownames(metrics_comat) <- NULL
  rownames(metrics_net) <- NULL

  expect_equal(metrics_comat, metrics_net, tolerance = 1e-10)
})

test_that("site_col and species_col as integers (column indices)", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  
  result <- site_species_metrics(
    bioregionalization = clust1,
    net = net_test,
    site_col = 1,
    species_col = 2,
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)
})

test_that("site_col and species_col as characters (column names)", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  colnames(net_test) <- c("SiteCol", "SpeciesCol", "WeightCol")
  
  result <- site_species_metrics(
    bioregionalization = clust1,
    net = net_test,
    site_col = "SiteCol",
    species_col = "SpeciesCol",
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)
})

test_that("parameter inheritance from bioregionalization$args", {
  
  # Create bioregionalization with stored parameters
  net_test <- mat_to_net(comat, weight = TRUE)
  colnames(net_test) <- c("Site", "Species", "Weight")
  
  # Create a clustering that stores args
  clust_with_args <- clust1
  clust_with_args$args$weight <- TRUE
  clust_with_args$args$index <- "Weight"
  clust_with_args$args$site_col <- "Site"
  clust_with_args$args$species_col <- "Species"
  
  # Call without explicitly providing parameters
  result <- site_species_metrics(
    bioregionalization = clust_with_args,
    net = net_test,
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(nrow(result[[1]]$metrics), 75)
})

test_that("explicitly provided parameters override bioregionalization", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  
  # Create clustering with stored args
  clust_with_args <- clust1
  clust_with_args$args$site_col <- 1
  clust_with_args$args$species_col <- 2
  
  # Override with explicit parameters
  net_reversed <- net_test[, c(2, 1, 3)]
  
  result <- site_species_metrics(
    bioregionalization = clust_with_args,
    net = net_reversed,
    site_col = 2,  # Explicitly override
    species_col = 1,  # Explicitly override
    indices = "rho"
  )
  
  expect_true(inherits(result, "bioregion.site.species.metrics"))
})

test_that("works with different clustering functions (netclu_*, hclu_*, nhclu_*)", {
  
  # netclu_*
  result_netclu <- site_species_metrics(
    bioregionalization = com,
    comat = comat,
    indices = "rho"
  )
  expect_true(inherits(result_netclu, "bioregion.site.species.metrics"))
  
  # nhclu_*
  result_nhclu <- site_species_metrics(
    bioregionalization = clust1,
    comat = comat,
    indices = "rho"
  )
  expect_true(inherits(result_nhclu, "bioregion.site.species.metrics"))
  
  # hclu_* (when cut)
  dissim_test <- dissimilarity(comat, metric = "Simpson")
  clust_hclu <- hclu_hierarclust(dissim_test, n_clust = 3)
  result_hclu <- site_species_metrics(
    bioregionalization = clust_hclu,
    comat = comat,
    indices = "rho"
  )
  expect_true(inherits(result_hclu, "bioregion.site.species.metrics"))
})

test_that("invalid site_col name errors appropriately", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  colnames(net_test) <- c("Site", "Species", "Weight")
  
  expect_error(
    site_species_metrics(
      bioregionalization = clust1,
      net = net_test,
      site_col = "NonExistentColumn",
      species_col = "Species",
      indices = "rho"
    ),
    "site_col 'NonExistentColumn' not found in net column names.",
    fixed = TRUE
  )
})

test_that("invalid species_col name errors appropriately", {
  
  net_test <- mat_to_net(comat, weight = TRUE)
  colnames(net_test) <- c("Site", "Species", "Weight")
  
  expect_error(
    site_species_metrics(
      bioregionalization = clust1,
      net = net_test,
      site_col = "Site",
      species_col = "NonExistentColumn",
      indices = "rho"
    ),
    "species_col 'NonExistentColumn' not found in net column names.",
    fixed = TRUE
  )
})

# Tests for Feature 2: Cz for Unipartite Networks ------------------------------
test_that("Cz calculations work for unipartite networks", {
  
  # Create test data
  comat_uni <- matrix(sample(0:10, 100, replace = TRUE), 10, 10)
  rownames(comat_uni) <- paste0("Site", 1:10)
  colnames(comat_uni) <- paste0("Species", 1:10)
  
  # Create clustering
  dissim_uni <- dissimilarity(comat_uni, metric = "Simpson")
  clust_uni <- nhclu_kmeans(dissim_uni, n_clust = 3)
  
  # Convert to network format
  net_uni <- mat_to_net(comat_uni, weight = TRUE)
  
  # Compute metrics including Cz
  suppressWarnings({
    metrics <- site_species_metrics(
      bioregionalization = clust_uni,
      comat = comat_uni,
      net = net_uni,
      indices = c("Cz", "indicator_value")
    )
  })
  
  # Check output structure
  expect_true("cz_metrics" %in% names(metrics[[1]]))
  expect_true("species_clusters" %in% names(metrics[[1]]))
  expect_true("Category" %in% colnames(metrics[[1]]$cz_metrics))
  
  # Check that species are assigned
  cz_data <- metrics[[1]]$cz_metrics
  expect_true("species" %in% unique(cz_data$Category))
  expect_true("site" %in% unique(cz_data$Category))
  
  # Check species_clusters structure
  sp_clust <- metrics[[1]]$species_clusters
  expect_equal(nrow(sp_clust), ncol(comat_uni))
  expect_true(all(sp_clust$Species %in% colnames(comat_uni)))
})

test_that("species assignment handles ties correctly", {
  
  # Create matrix where some species have identical patterns
  comat_tied <- matrix(c(
    1, 1, 1, 1,
    1, 1, 1, 1,
    0, 0, 0, 0,
    0, 0, 0, 0,
    2, 2, 0, 0,
    2, 2, 0, 0
  ), nrow = 6, byrow = TRUE)
  rownames(comat_tied) <- paste0("Site", 1:6)
  colnames(comat_tied) <- paste0("Species", 1:4)
  
  # Create clustering - Sites 1-2 vs 3-4 vs 5-6
  clust_tied <- data.frame(
    ID = paste0("Site", 1:6),
    K_3 = c(1, 1, 2, 2, 3, 3)
  )
  class(clust_tied) <- c("bioregion.clusters", "data.frame")
  
  bioregion_tied <- list(
    name = "test",
    args = list(),
    inputs = list(bipartite = FALSE),
    algorithm = list(),
    clusters = clust_tied,
    cluster_info = list(n_clust = 3)
  )
  class(bioregion_tied) <- "bioregion.clusters"
  
  net_tied <- mat_to_net(comat_tied, weight = TRUE)
  
  # Compute Cz metrics
  suppressWarnings({
    result <- site_species_metrics(
      bioregionalization = bioregion_tied,
      comat = comat_tied,
      net = net_tied,
      indices = "Cz"
    )
  })
  
  # Check that tied_species is present if ties occurred
  # Species 1 and 2 have identical patterns, so may have ties
  expect_true("species_clusters" %in% names(result[[1]]))
  
  # All species should be assigned to some bioregion
  sp_clust <- result[[1]]$species_clusters
  expect_equal(nrow(sp_clust), 4)
  expect_true(all(sp_clust$Bioregion %in% c(1, 2, 3)))
})

test_that("Cz output structure for unipartite includes expected elements", {
  
  comat_uni <- matrix(sample(0:10, 100, replace = TRUE), 10, 10)
  rownames(comat_uni) <- paste0("Site", 1:10)
  colnames(comat_uni) <- paste0("Species", 1:10)
  
  dissim_uni <- dissimilarity(comat_uni, metric = "Simpson")
  clust_uni <- nhclu_kmeans(dissim_uni, n_clust = 3)
  net_uni <- mat_to_net(comat_uni, weight = TRUE)
  
  suppressWarnings({
    result <- site_species_metrics(
      bioregionalization = clust_uni,
      comat = comat_uni,
      net = net_uni,
      indices = "Cz"
    )
  })
  
  # Check all expected elements are present
  expect_true("name" %in% names(result[[1]]))
  expect_true("indices" %in% names(result[[1]]))
  expect_true("cz_metrics" %in% names(result[[1]]))
  expect_true("species_clusters" %in% names(result[[1]]))
  expect_true("args" %in% names(result[[1]]))
  
  # Check cz_metrics has required columns
  cz_data <- result[[1]]$cz_metrics
  expect_true(all(c("Node", "Bioregion", "Category", "C", "z") %in% 
                    colnames(cz_data)))
})

test_that("bipartite vs unipartite Cz metrics have similar structure", {
  
  # Unipartite
  comat_uni <- matrix(sample(0:10, 100, replace = TRUE), 10, 10)
  rownames(comat_uni) <- paste0("Site", 1:10)
  colnames(comat_uni) <- paste0("Species", 1:10)
  
  dissim_uni <- dissimilarity(comat_uni, metric = "Simpson")
  clust_uni <- nhclu_kmeans(dissim_uni, n_clust = 3)
  net_uni <- mat_to_net(comat_uni, weight = TRUE)
  
  suppressWarnings({
    result_uni <- site_species_metrics(
      bioregionalization = clust_uni,
      comat = comat_uni,
      net = net_uni,
      indices = "Cz"
    )
  })
  
  # Bipartite
  net_bip <- mat_to_net(comat, weight = TRUE)
  clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)
  
  suppressWarnings({
    result_bip <- site_species_metrics(
      bioregionalization = clust_bip,
      comat = comat,
      net = net_bip,
      indices = "Cz"
    )
  })
  
  # Both should have cz_metrics
  expect_true("cz_metrics" %in% names(result_uni[[1]]))
  expect_true("cz_metrics" %in% names(result_bip[[1]]))
  
  # Both should have same column structure
  expect_equal(
    names(result_uni[[1]]$cz_metrics),
    names(result_bip[[1]]$cz_metrics)
  )
  
  # Unipartite should have species_clusters, bipartite should not
  expect_true("species_clusters" %in% names(result_uni[[1]]))
  expect_false("species_clusters" %in% names(result_bip[[1]]))
})

test_that("Cz edge cases: single bioregion", {
  
  comat_single <- matrix(sample(0:10, 50, replace = TRUE), 5, 10)
  rownames(comat_single) <- paste0("Site", 1:5)
  colnames(comat_single) <- paste0("Species", 1:10)
  
  # Create clustering with single bioregion
  clust_single <- data.frame(
    ID = paste0("Site", 1:5),
    K_1 = rep(1, 5)
  )
  class(clust_single) <- c("bioregion.clusters", "data.frame")
  
  bioregion_single <- list(
    name = "test",
    args = list(),
    inputs = list(bipartite = FALSE),
    algorithm = list(),
    clusters = clust_single,
    cluster_info = list(n_clust = 1)
  )
  class(bioregion_single) <- "bioregion.clusters"
  
  net_single <- mat_to_net(comat_single, weight = TRUE)
  
  # Should work without errors
  suppressWarnings({
    result <- site_species_metrics(
      bioregionalization = bioregion_single,
      comat = comat_single,
      net = net_single,
      indices = "Cz"
    )
  })
  
  expect_true("cz_metrics" %in% names(result[[1]]))
  expect_true("species_clusters" %in% names(result[[1]]))
})

test_that("Cz edge cases: species present in only one bioregion", {
  
  # Create matrix where species are exclusive to bioregions
  comat_exclusive <- matrix(c(
    1, 1, 0, 0,
    1, 1, 0, 0,
    0, 0, 1, 1,
    0, 0, 1, 1
  ), nrow = 4, byrow = TRUE)
  rownames(comat_exclusive) <- paste0("Site", 1:4)
  colnames(comat_exclusive) <- paste0("Species", 1:4)
  
  clust_exclusive <- data.frame(
    ID = paste0("Site", 1:4),
    K_2 = c(1, 1, 2, 2)
  )
  class(clust_exclusive) <- c("bioregion.clusters", "data.frame")
  
  bioregion_exclusive <- list(
    name = "test",
    args = list(),
    inputs = list(bipartite = FALSE),
    algorithm = list(),
    clusters = clust_exclusive,
    cluster_info = list(n_clust = 2)
  )
  class(bioregion_exclusive) <- "bioregion.clusters"
  
  net_exclusive <- mat_to_net(comat_exclusive, weight = TRUE)
  
  suppressWarnings({
    result <- site_species_metrics(
      bioregionalization = bioregion_exclusive,
      comat = comat_exclusive,
      net = net_exclusive,
      indices = c("Cz", "indicator_value")
    )
  })
  
  # Species should be assigned to their exclusive bioregions
  sp_clust <- result[[1]]$species_clusters
  expect_equal(nrow(sp_clust), 4)
  
  # Species 1-2 should be in bioregion 1, Species 3-4 in bioregion 2
  expect_true(all(sp_clust$Species[sp_clust$Bioregion == 1] %in% c("Species1", "Species2")))
  expect_true(all(sp_clust$Species[sp_clust$Bioregion == 2] %in% c("Species3", "Species4")))
  
  # Check that indicator values are perfect (1.0) for exclusive species
  metrics <- result[[1]]$metrics
  sp1_br1 <- metrics[metrics$Species == "Species1" & metrics$Bioregion == 1, "indval"]
  expect_equal(sp1_br1, 1.0)
})

