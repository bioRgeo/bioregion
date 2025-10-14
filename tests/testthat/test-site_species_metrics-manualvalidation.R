# Validation Tests for site_species_metrics()
# 
# This script validates the calculations in site_species_metrics() by comparing 
# function outputs against manually calculated expected values for simple, 
# test cases. 

# Source helper functions and manual calculations
source(test_path("../helpers/helper_functions_validation.R"))
source(test_path("../helpers/manual_calculations_site_species_metrics.R"))

# Test 1: Perfect Partitioning - Affinity, Fidelity, IndVal ---------------
test_that("Affinity, Fidelity, IndVal - Perfect Partitioning", {
  
  # Get test case data with manual calculations
  tc1 <- calculate_test_case_1()
  
  # Create proper bioregion.clusters object
  mock_clust <- create_mock_clusters(tc1$comat, tc1$clusters)
  
  # Run site_species_metrics function
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc1$comat,
    indices = c("affinity", "fidelity", "indicator_value"),
    verbose = FALSE
  )
  
  # Extract metrics from result
  # Result is a list with one element (single bioregionalization)
  metrics_df <- result[[1]]$metrics
  
  # Test structure
  expect_true(inherits(result, "bioregion.site.species.metrics"))
  expect_equal(length(result), 1)
  expect_equal(attr(result, "n_bioregionalizations"), 1)
  
  # Compare with expected values
  # For each species in the expected data, find corresponding row in results
  for (i in 1:nrow(tc1$expected)) {
    species <- tc1$expected$Species[i]
    bioregion <- tc1$expected$Bioregion[i]
    
    # Find matching row in actual results
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    # Check that we found exactly one matching row
    expect_equal(nrow(result_row), 1,
                info = sprintf("Expected 1 row for %s in Bioregion %s", 
                              species, bioregion))
    
    # Compare affinity (within tolerance for floating-point arithmetic)
    expect_equal(result_row$affinity, 
                tc1$expected$affinity[i], 
                tolerance = 1e-10,
                label = sprintf("Affinity for %s in Bioregion %s", 
                               species, bioregion))
    
    # Compare fidelity
    expect_equal(result_row$fidelity, 
                tc1$expected$fidelity[i], 
                tolerance = 1e-10,
                label = sprintf("Fidelity for %s in Bioregion %s", 
                               species, bioregion))
    
    # Compare indicator value (IndVal)
    expect_equal(result_row$indval, 
                tc1$expected$indval[i], 
                tolerance = 1e-10,
                label = sprintf("IndVal for %s in Bioregion %s", 
                               species, bioregion))
  }
  
  # Additional check: verify that perfect partitioning gives perfect metrics
  # All species should have affinity = 1.0, fidelity = 1.0, indval = 1.0
  expect_true(all(tc1$expected$affinity == 1.0))
  expect_true(all(tc1$expected$fidelity == 1.0))
  expect_true(all(tc1$expected$indval == 1.0))
})


# Test 2: Partial Overlap - Affinity, Fidelity, IndVal --------------------
test_that("Affinity, Fidelity, IndVal - Partial Overlap", {
  
  tc2 <- calculate_test_case_2()
  
  mock_clust <- create_mock_clusters(tc2$comat, tc2$clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc2$comat,
    indices = c("affinity", "fidelity", "indicator_value"),
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Compare with expected values
  for (i in 1:nrow(tc2$expected)) {
    species <- tc2$expected$Species[i]
    bioregion <- tc2$expected$Bioregion[i]
    
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    expect_equal(nrow(result_row), 1)
    
    expect_equal(result_row$affinity, 
                tc2$expected$affinity[i], 
                tolerance = 1e-10,
                label = sprintf("Affinity for %s in Bioregion %s", 
                               species, bioregion))
    
    expect_equal(result_row$fidelity, 
                tc2$expected$fidelity[i], 
                tolerance = 1e-10,
                label = sprintf("Fidelity for %s in Bioregion %s", 
                               species, bioregion))
    
    expect_equal(result_row$indval, 
                tc2$expected$indval[i], 
                tolerance = 1e-10,
                label = sprintf("IndVal for %s in Bioregion %s", 
                               species, bioregion))
  }
})


# Test 3: Three Bioregions with Generalists --------------------------------
test_that("Affinity, Fidelity, IndVal - Three Bioregions with Generalists", {
  
  tc3 <- calculate_test_case_3()
  
  mock_clust <- create_mock_clusters(tc3$comat, tc3$clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc3$comat,
    indices = c("affinity", "fidelity", "indicator_value"),
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Compare with expected values
  for (i in 1:nrow(tc3$expected)) {
    species <- tc3$expected$Species[i]
    bioregion <- tc3$expected$Bioregion[i]
    
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    expect_equal(nrow(result_row), 1)
    
    expect_equal(result_row$affinity, 
                tc3$expected$affinity[i], 
                tolerance = 1e-10,
                label = sprintf("Affinity for %s in Bioregion %s", 
                               species, bioregion))
    
    expect_equal(result_row$fidelity, 
                tc3$expected$fidelity[i], 
                tolerance = 1e-10,
                label = sprintf("Fidelity for %s in Bioregion %s", 
                               species, bioregion))
    
    expect_equal(result_row$indval, 
                tc3$expected$indval[i], 
                tolerance = 1e-10,
                label = sprintf("IndVal for %s in Bioregion %s", 
                               species, bioregion))
  }
  
})


# Test 4: Rho Calculation - Perfect Partitioning ---------------------------
test_that("Rho calculation - Perfect Partitioning", {
  
  tc1 <- calculate_test_case_1()
  
  mock_clust <- create_mock_clusters(tc1$comat, tc1$clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc1$comat,
    indices = "rho",
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Rho should be calculated for ALL species-bioregion combinations
  # (not just the assigned bioregion)
  expect_equal(nrow(metrics_df), 8)  # 4 species * 2 bioregions
  
  # Compare with manually calculated rho values
  for (i in 1:nrow(tc1$expected_rho)) {
    species <- tc1$expected_rho$Species[i]
    bioregion <- tc1$expected_rho$Bioregion[i]
    
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    expect_equal(nrow(result_row), 1)
    
    # Rho values should match manual calculations
    # Using slightly larger tolerance due to sqrt operations
    expect_equal(result_row$rho, 
                tc1$expected_rho$rho[i], 
                tolerance = 1e-6,
                label = sprintf("Rho for %s in Bioregion %s", 
                               species, bioregion))
  }
  
  # Additional checks:
  # - Species in their own bioregion should have positive rho
  # - Species absent from a bioregion should have negative rho
  sp_a_br1 <- metrics_df[metrics_df$Species == "Sp_A" & 
                          metrics_df$Bioregion == 1, "rho"]
  sp_a_br2 <- metrics_df[metrics_df$Species == "Sp_A" & 
                          metrics_df$Bioregion == 2, "rho"]
  
  expect_true(sp_a_br1 > 0, label = "Sp_A should have positive rho in Bioregion 1")
  expect_true(sp_a_br2 < 0, label = "Sp_A should have negative rho in Bioregion 2")
})


# Test 5: Rho Calculation - Partial Overlap --------------------------------
test_that("Rho calculation - Partial Overlap", {
  
  tc2 <- calculate_test_case_2()
  
  mock_clust <- create_mock_clusters(tc2$comat, tc2$clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc2$comat,
    indices = "rho",
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Compare with manually calculated rho values
  for (i in 1:nrow(tc2$expected_rho)) {
    species <- tc2$expected_rho$Species[i]
    bioregion <- tc2$expected_rho$Bioregion[i]
    
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    expect_equal(result_row$rho, 
                tc2$expected_rho$rho[i], 
                tolerance = 1e-6,
                label = sprintf("Rho for %s in Bioregion %s", 
                               species, bioregion))
  }
  
  # Sp_B is present in both bioregions equally, so rho should be 0
  sp_b_br1 <- metrics_df[metrics_df$Species == "Sp_B" & 
                          metrics_df$Bioregion == 1, "rho"]
  sp_b_br2 <- metrics_df[metrics_df$Species == "Sp_B" & 
                          metrics_df$Bioregion == 2, "rho"]
  
  expect_equal(sp_b_br1, 0, tolerance = 1e-10,
              label = "Sp_B should have rho = 0 in Bioregion 1 (equal distribution)")
  expect_equal(sp_b_br2, 0, tolerance = 1e-10,
              label = "Sp_B should have rho = 0 in Bioregion 2 (equal distribution)")
})


# Test 6: Rho Calculation - Minimal Example --------------------------------
test_that("Rho calculation - Minimal Example (3 sites, 2 species)", {
  
  tc5 <- calculate_test_case_5()
  
  mock_clust <- create_mock_clusters(tc5$comat, tc5$clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc5$comat,
    indices = "rho",
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Compare with manually calculated rho values
  for (i in 1:nrow(tc5$expected_rho)) {
    species <- tc5$expected_rho$Species[i]
    bioregion <- tc5$expected_rho$Bioregion[i]
    
    result_row <- metrics_df[
      metrics_df$Species == species & 
      metrics_df$Bioregion == bioregion, 
    ]
    
    expect_equal(result_row$rho, 
                tc5$expected_rho$rho[i], 
                tolerance = 1e-6,
                label = sprintf("Rho for %s in Bioregion %s (minimal example)", 
                               species, bioregion))
  }
})


# Test 7: Cz Indices - Bipartite Network -----------------------------------
test_that("Cz indices - Bipartite network", {
  
  tc4 <- calculate_test_case_4()
  
  # Create bipartite cluster object with both site and species clusters
  mock_clust_bip <- create_mock_clusters(
    tc4$comat, 
    tc4$site_clusters,
    bipartite = TRUE,
    species_clusters = tc4$species_clusters,
    weight = TRUE,
    weight_index = 3  # Third column in net dataframe is the weight
  )
  
  # Create network data.frame for Cz calculations
  net <- mat_to_net(tc4$comat, weight = TRUE)
  
  # No warning expected anymore - function uses both comat and net
  result <- site_species_metrics(
    bioregionalization = mock_clust_bip,
    comat = tc4$comat,
    net = net,
    indices = "Cz",
    verbose = FALSE
  )
  
  # Cz metrics are stored separately from other metrics
  cz_df <- result[[1]]$cz_metrics
  
  # Check structure
  expect_true("Node" %in% colnames(cz_df))
  expect_true("Bioregion" %in% colnames(cz_df))
  expect_true("C" %in% colnames(cz_df))
  expect_true("z" %in% colnames(cz_df))
  
  # Compare with expected Cz values
  for (i in 1:nrow(tc4$expected_cz)) {
    node <- tc4$expected_cz$Node[i]
    
    result_row <- cz_df[cz_df$Node == node, ]
    
    expect_equal(nrow(result_row), 1,
                info = sprintf("Expected 1 row for node %s", node))
    
    # Compare participation coefficient (C)
    # C should be 0 for perfect within-bioregion connectivity
    expect_equal(result_row$C, 
                tc4$expected_cz$C[i], 
                tolerance = 1e-6,
                label = sprintf("Participation coefficient for %s", node))
    
    # Compare z-score
    # In this test case, all nodes have the same degree as their bioregion mean
    # so z-score should be 0 (or NA if sd = 0)
    if (!is.na(tc4$expected_cz$z[i])) {
      expect_equal(result_row$z, 
                  tc4$expected_cz$z[i], 
                  tolerance = 1e-6,
                  label = sprintf("Z-score for %s", node))
    }
  }

})


# Test 8: All Indices Combined ---------------------------------------------
test_that("All indices computed together", {
  
  tc2 <- calculate_test_case_2()
  
  mock_clust <- create_mock_clusters(tc2$comat, tc2$clusters)
  
  # Request all non-Cz indices at once
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = tc2$comat,
    indices = c("rho", "affinity", "fidelity", "indicator_value"),
    verbose = FALSE
  )
  
  metrics_df <- result[[1]]$metrics
  
  # Check that all columns are present
  expect_true(all(c("Bioregion", "Species", "rho", "affinity", 
                   "fidelity", "indval") %in% colnames(metrics_df)))
  
  # Verify that metrics match expected values for at least one species
  sp_a_br1 <- metrics_df[metrics_df$Species == "Sp_A" & 
                          metrics_df$Bioregion == 1, ]
  
  expect_equal(sp_a_br1$affinity, 1.0, tolerance = 1e-10)
  expect_equal(sp_a_br1$fidelity, 1.0, tolerance = 1e-10)
  expect_equal(sp_a_br1$indval, 1.0, tolerance = 1e-10)
  expect_equal(sp_a_br1$rho, 1.732051, tolerance = 1e-6)
})


# Test 9: Edge Case - Single Site in a Bioregion --------------------------
test_that("Edge case - Single site in a bioregion", {
  
  # Create a matrix with an uneven bioregionalization
  # Bioregion 1: 1 site, Bioregion 2: 3 sites
  comat <- matrix(c(
    1, 1, 0,  # Site1 (Br1)
    0, 1, 1,  # Site2 (Br2)
    0, 1, 1,  # Site3 (Br2)
    0, 1, 1   # Site4 (Br2)
  ), nrow = 4, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:4)
  colnames(comat) <- paste0("Sp_", LETTERS[1:3])
  
  clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c(1, 2, 2, 2),
    stringsAsFactors = FALSE
  )
  
  mock_clust <- create_mock_clusters(comat, clusters)
  
  # Should not error even with single-site bioregion
  expect_no_error({
    result <- site_species_metrics(
      bioregionalization = mock_clust,
      comat = comat,
      indices = c("rho", "affinity", "fidelity", "indicator_value"),
      verbose = FALSE
    )
  })
  
  metrics_df <- result[[1]]$metrics
  
  # Sp_A is only in the single-site bioregion
  sp_a_br1 <- metrics_df[metrics_df$Species == "Sp_A" & 
                          metrics_df$Bioregion == 1, ]
  
  # Affinity should be 1.0 (present in the only site)
  expect_equal(sp_a_br1$affinity, 1.0, tolerance = 1e-10)
  
  # Fidelity should be 1.0 (only in this bioregion)
  expect_equal(sp_a_br1$fidelity, 1.0, tolerance = 1e-10)
})


# Test 10: Multiple Bioregionalizations -----------------------------------
test_that("Multiple bioregionalizations (K_2 and K_3)", {
  
  # Create a matrix with multiple bioregionalizations
  comat <- matrix(c(
    1, 1, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0,
    0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 1, 1
  ), nrow = 6, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:6)
  colnames(comat) <- paste0("Sp_", LETTERS[1:6])
  
  # Two bioregionalizations: K_2 and K_3
  clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c(1, 1, 2, 2, 2, 2),  # 2 bioregions
    K_3 = c(1, 1, 2, 2, 3, 3),  # 3 bioregions
    stringsAsFactors = FALSE
  )
  
  mock_clust <- create_mock_clusters(comat, clusters)
  
  result <- site_species_metrics(
    bioregionalization = mock_clust,
    comat = comat,
    indices = c("affinity", "fidelity"),
    verbose = FALSE
  )
  
  # Should return a list with 2 elements (one per bioregionalization)
  expect_equal(length(result), 2)
  expect_equal(names(result), c("K_2", "K_3"))
  
  # Check K_2 results
  k2_metrics <- result$K_2$metrics
  expect_true(all(k2_metrics$Bioregion %in% c(1, 2)))
  
  # Check K_3 results
  k3_metrics <- result$K_3$metrics
  expect_true(all(k3_metrics$Bioregion %in% c(1, 2, 3)))
  
  # Verify that Sp_A has perfect metrics in both bioregionalizations
  sp_a_k2 <- k2_metrics[k2_metrics$Species == "Sp_A" & k2_metrics$Bioregion == 1, ]
  expect_equal(sp_a_k2$affinity, 1.0, tolerance = 1e-10)
  expect_equal(sp_a_k2$fidelity, 1.0, tolerance = 1e-10)
  

  sp_a_k3 <- k3_metrics[k3_metrics$Species == "Sp_A" & k3_metrics$Bioregion == 1, ]
  expect_equal(sp_a_k3$affinity, 1.0, tolerance = 1e-10)
  expect_equal(sp_a_k3$fidelity, 1.0, tolerance = 1e-10)
})
