# Test summary.bioregion.clusters function
# Testing with different clustering algorithms and configurations

# Install binaries for infomap tests
quietly(install_binaries())

# Test data setup --------------------------------------------------------------
dissim_fish <- dissimilarity(fishmat, metric = "Simpson")
sim_fish <- similarity(fishmat, metric = "Simpson")

# Test 1: Non-hierarchical clustering ------------------------------------------
test_that("summary works for non-hierarchical clustering (nhclu_pam)", {
  skip_on_cran()
  
  clust_pam <- nhclu_pam(dissim_fish, n_clust = 5, seed = 123)
  
  # Test that summary runs without error
  expect_no_error(quietly(summary(clust_pam)))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_pam)
  expect_identical(result, clust_pam)
  
  # Test with different parameters
  expect_no_error(quietly(summary(clust_pam, n_bioregionalizations = 1, n_top_clusters = 3)))
  expect_no_error(quietly(summary(clust_pam, n_bioregionalizations = 10, n_top_clusters = 20)))
})

# Test 2: Hierarchical clustering (hclu_hierarclust) --------------------------
test_that("summary works for hierarchical clustering (hclu_hierarclust)", {
  skip_on_cran()
  
  clust_hclust <- hclu_hierarclust(
    dissim_fish,
    n_clust = c(3, 10, 50),
    optimal_tree_method = "best",
    randomize = FALSE
  )
  
  # Test that summary runs without error
  expect_no_error(quietly(summary(clust_hclust)))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_hclust, n_bioregionalizations = 3, n_top_clusters = 10)
  expect_identical(result, clust_hclust)
  
  # Test with different parameters
  expect_no_error(quietly(summary(clust_hclust, n_bioregionalizations = 1, n_top_clusters = 5)))
  expect_no_error(quietly(summary(clust_hclust, n_bioregionalizations = 2, n_top_clusters = 15)))
})

# Test 3: Network clustering - unipartite, no hierarchy ------------------------
test_that("summary works for unipartite netclu_infomap without hierarchy", {
  skip_on_cran()
  skip_if_not(isTRUE(file.exists(
    system.file("bin", "infomap", package = "bioregion")
  )))
  
  clust_infomap_uni_flat <- netclu_infomap(
    sim_fish,
    bipartite = FALSE,
    bipartite_version = FALSE,
    show_hierarchy = FALSE,
    nbmod = 5,
    seed = 123
  )
  
  # Test that summary runs without error
  expect_no_error(quietly(summary(clust_infomap_uni_flat)))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_infomap_uni_flat, n_bioregionalizations = 3, n_top_clusters = 10)
  expect_identical(result, clust_infomap_uni_flat)
  
  # Check that hierarchical structure is not displayed
  output <- capture.output(quietly(summary(clust_infomap_uni_flat)))
  expect_false(any(grepl("Hierarchical structure", output)))
})

# Test 4: Network clustering - unipartite, with hierarchy ---------------------
test_that("summary works for unipartite netclu_infomap with hierarchy", {
  skip_on_cran()
  skip_if_not(isTRUE(file.exists(
    system.file("bin", "infomap", package = "bioregion")
  )))
  
  clust_infomap_uni_hier <- netclu_infomap(
    sim_fish,
    bipartite = FALSE,
    bipartite_version = FALSE,
    show_hierarchy = TRUE,
    seed = 123
  )
  
  # Test that summary runs without error
  expect_no_error(summary(clust_infomap_uni_hier))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_infomap_uni_hier, n_bioregionalizations = 3, n_top_clusters = 10)
  expect_identical(result, clust_infomap_uni_hier)
  
  # Check that hierarchical structure IS displayed
  output <- capture.output(summary(clust_infomap_uni_hier))
  expect_true(any(grepl("Hierarchical structure", output)))
  
  # Check for tree characters in hierarchy display
  expect_true(any(grepl("├─|└─", output)))
})

# Test 5: Network clustering - bipartite, no hierarchy ------------------------
test_that("summary works for bipartite netclu_infomap without hierarchy", {
  skip_on_cran()
  skip_if_not(isTRUE(file.exists(
    system.file("bin", "infomap", package = "bioregion")
  )))
  
  clust_infomap_bip_flat <- netclu_infomap(
    sim_fish,
    bipartite = TRUE,
    bipartite_version = FALSE,
    show_hierarchy = FALSE,
    nbmod = 5,
    seed = 123
  )
  
  # Test that summary runs without error
  expect_no_error(summary(clust_infomap_bip_flat))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_infomap_bip_flat, n_bioregionalizations = 3, n_top_clusters = 10)
  expect_identical(result, clust_infomap_bip_flat)
  
  # Check that site/species breakdown is shown
  output <- capture.output(summary(clust_infomap_bip_flat))
  expect_true(any(grepl("Site nodes:", output)))
  expect_true(any(grepl("Species nodes:", output)))
  expect_true(any(grepl("sites.*species", output)))
})

# Test 6: Network clustering - bipartite, with hierarchy ----------------------
test_that("summary works for bipartite netclu_infomap with hierarchy", {
  skip_on_cran()
  skip_if_not(isTRUE(file.exists(
    system.file("bin", "infomap", package = "bioregion")
  )))
  
  clust_infomap_bip_hier <- netclu_infomap(
    sim_fish,
    bipartite = TRUE,
    bipartite_version = FALSE,
    show_hierarchy = TRUE,
    seed = 123
  )
  
  # Test that summary runs without error
  expect_no_error(summary(clust_infomap_bip_hier))
  
  # Test that summary returns the object invisibly
  result <- summary(clust_infomap_bip_hier, n_bioregionalizations = 3, n_top_clusters = 10)
  expect_identical(result, clust_infomap_bip_hier)
  
  # Check that hierarchical structure IS displayed
  output <- capture.output(summary(clust_infomap_bip_hier))
  expect_true(any(grepl("Hierarchical structure", output)))
  
  # Check that site/species breakdown is shown
  expect_true(any(grepl("Site nodes:", output)))
  expect_true(any(grepl("Species nodes:", output)))
  
  # Check for tree characters in hierarchy display
  expect_true(any(grepl("├─|└─", output)))
})

# Test 7: Parameter validation -------------------------------------------------
test_that("summary respects n_bioregionalizations and n_top_clusters parameters", {
  skip_on_cran()
  
  clust_pam <- nhclu_pam(dissim_fish, n_clust = 5, seed = 123)
  
  # Test n_bioregionalizations = 1 shows only 1 bioregionalization
  output1 <- capture.output(summary(clust_pam, n_bioregionalizations = 1, n_top_clusters = 5))
  bioregionalization_lines1 <- grep("^Bioregionalization [0-9]+:", output1)
  expect_equal(length(bioregionalization_lines1), 1)
  
  # Test n_top_clusters limits the number of clusters shown
  output2 <- capture.output(summary(clust_pam, n_bioregionalizations = 1, n_top_clusters = 3))
  cluster_lines2 <- grep("^  Cluster", output2)
  expect_true(length(cluster_lines2) <= 3)
})

# Test 8: Edge cases -----------------------------------------------------------
test_that("summary handles edge cases correctly", {
  skip_on_cran()
  
  clust_pam <- nhclu_pam(dissim_fish, n_clust = 5, seed = 123)
  
  # Test with n_bioregionalizations larger than available
  expect_no_error(summary(clust_pam, n_bioregionalizations = 100))
  
  # Test with n_top_clusters larger than available
  expect_no_error(summary(clust_pam, n_top_clusters = 1000))
  
  # Test with n_bioregionalizations = 0 or negative (should show at least available)
  expect_no_error(summary(clust_pam, n_bioregionalizations = 0))
})

# Test 9: Display hierarchy function -------------------------------------------
test_that("display_hierarchy function works correctly", {
  skip_on_cran()
  skip_if_not(isTRUE(file.exists(
    system.file("bin", "infomap", package = "bioregion")
  )))
  
  clust_infomap_hier <- netclu_infomap(
    sim_fish,
    bipartite = FALSE,
    bipartite_version = FALSE,
    show_hierarchy = TRUE,
    seed = 123
  )
  
  # Test that display_hierarchy function is called properly
  output <- capture.output(summary(clust_infomap_hier, n_bioregionalizations = 3))
  
  # Check for proper tree structure
  has_tree_chars <- any(grepl("├─", output)) || any(grepl("└─", output))
  has_vertical_line <- any(grepl("│", output))
  
  expect_true(has_tree_chars)
  
  # Check for cluster sizes in hierarchy
  expect_true(any(grepl("\\(n=[0-9]+\\)", output)))
})

# Test 10: Output consistency --------------------------------------------------
test_that("summary output is consistent across runs", {
  skip_on_cran()
  
  clust_pam <- nhclu_pam(dissim_fish, n_clust = 5, seed = 123)
  
  # Capture output twice
  output1 <- capture.output(summary(clust_pam))
  output2 <- capture.output(summary(clust_pam))
  
  # Output should be identical
  expect_identical(output1, output2)
})
