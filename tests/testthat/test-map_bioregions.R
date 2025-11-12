# Inputs -----------------------------------------------------------------------
data("fishmat")
data("fishsf")
net <- similarity(fishmat, metric = "Simpson")
clu <- netclu_greedy(net)

cludf <- clu$clusters
cludfNA <- cludf
cludfNA[1,1] <- NA

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {

  map <- map_bioregions(clusters = clu, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  expect_equal(class(map)[2], "data.frame")
  
  map <- map_bioregions(clusters = clu, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = TRUE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  expect_equal(class(map)[2], "data.frame")
   
})
 
# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
   
  expect_error(
    map_bioregions(clusters = "zz"),
   "If not a bioregion.clusters's object, clusters must be a data.frame.",
   fixed = TRUE)
  
  expect_error(
    map_bioregions(clusters = data.frame(ID = rep(1,5))),
    "clusters must be a data.frame with at least two columns.",
    fixed = TRUE)
  
  expect_error(
    expect_message(
      map_bioregions(clusters = data.frame(ID = rep(1,5), 
                                           ID2 = rep(1,5)),
                     geometry = fishsf),
      "Duplicated site ID detected!",
      fixed = TRUE),
    "Some cluster sites are not found in the geometry")
  
  expect_error(
    map_bioregions(clusters = cludfNA, 
                   geometry = fishsf, 
                   write_clusters = TRUE, 
                   plot = FALSE),
    "^NA")
  
  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = "zz"),
    "It seems that the geometry used is not an sf object.",
    fixed = TRUE)

  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = fishsf, 
                   write_clusters = NA, 
                   plot = FALSE),
    "write_clusters must be a boolean.",
    fixed = TRUE)

  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = fishsf, 
                   write_clusters = TRUE, 
                   plot = NA),
    "plot must be a boolean.",
    fixed = TRUE)
  
})

# Tests for color integration --------------------------------------------------
test_that("color integration works", {
  
  # Add colors to clusters
  clu_colored <- bioregion_colors(clu)
  
  # Test that function runs without error with colors
  expect_no_error(
    map <- map_bioregions(clusters = clu_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  )
  
  # Test that output is still valid
  map <- map_bioregions(clusters = clu_colored, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  
  # Test with multiple partitions
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim,  
                               optimal_tree_method = "best", 
                               n_clust = c(3, 4, 8),
                               verbose = FALSE)
  clu_hier_colored <- bioregion_colors(clu_hier)
  
  expect_no_error(
    map <- map_bioregions(clusters = clu_hier_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  )
  
  map <- map_bioregions(clusters = clu_hier_colored, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 5L)  # ID + 3 partitions + geometry
  expect_equal(class(map)[1], "sf")
})

test_that("backwards compatibility without colors", {
  
  # Test that clusters without colors still work as before
  map <- map_bioregions(clusters = clu, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  
  # Test with data.frame input (no colors possible)
  expect_no_error(
    map <- map_bioregions(clusters = cludf, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  )
})

test_that("bipartite networks", {
  
  # Create bipartite network clusters
  net_bip <- mat_to_net(fishmat, weight = TRUE)
  clu_bip <- netclu_greedy(net_bip, bipartite = TRUE)
  
  # Check that node_type attribute exists
  expect_equal(length(attr(clu_bip$clusters, "node_type")), 533L)
  expect_equal(sum(attr(clu_bip$clusters, "node_type") == "site"), 338L)
  expect_equal(sum(attr(clu_bip$clusters, "node_type") == "species"), 195L)
  
  # Test mapping works (should filter to sites only)
  map <- map_bioregions(clusters = clu_bip, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)  # Only sites, not species
  expect_equal(class(map)[1], "sf")
  
  # Test with colors
  clu_bip_colored <- bioregion_colors(clu_bip)
  expect_no_error(
    map <- map_bioregions(clusters = clu_bip_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  )
  expect_equal(nrow(map), 338L)
})

# Tests for multi-panel data structure (no actual plotting to avoid issues)
test_that("multi-panel data structure", {
  
  # Create hierarchical clustering with multiple partitions
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = c(3, 4, 8),
                               verbose = FALSE)
  clu_hier_colored <- bioregion_colors(clu_hier)
  
  # Test that write_clusters returns correct structure with multiple partitions
  map <- map_bioregions(clusters = clu_hier_colored, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 5L)  # ID + 3 partitions + geometry
  expect_equal(class(map)[1], "sf")
  
  # Verify all partition columns exist
  expect_true("K_3" %in% colnames(map))
  expect_true("K_4" %in% colnames(map))
  expect_true("K_8" %in% colnames(map))
  
  # Test with multiple partitions (use values that work reliably)
  clu_many <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = c(2, 3, 4),
                               verbose = FALSE)
  clu_many_colored <- bioregion_colors(clu_many)
  
  map_many <- map_bioregions(clusters = clu_many_colored, 
                             geometry = fishsf, 
                             write_clusters = TRUE, 
                             plot = FALSE)
  expect_equal(nrow(map_many), 338L)
  expect_equal(ncol(map_many), 5L)  # ID + 3 partitions + geometry
  expect_true(all(c("K_2", "K_3", "K_4") %in% colnames(map_many)))
  
  # Verify colors are present
  expect_true(!is.null(clu_many_colored$colors))
  expect_equal(length(clu_many_colored$colors), 3L)
})

test_that("multi-panel unipartite with colors and plotting", {
  data(fishmat)
  data(fishsf)
  
  dissim <- similarity_to_dissimilarity(
    similarity(fishmat, metric = "Simpson")
  )
  
  # Create clusters with multiple partitions and colors
  clu <- hclu_hierarclust(dissim, 
                          optimal_tree_method = "best", 
                          n_clust = c(2, 3, 4),
                          verbose = FALSE)
  clu_colored <- bioregion_colors(clu)
  
  # Test that map_bioregions works with colors (structure check)
  map_result <- map_bioregions(clusters = clu_colored, 
                                geometry = fishsf, 
                                write_clusters = TRUE, 
                                plot = FALSE)
  
  # Verify structure
  expect_equal(nrow(map_result), 338L)
  expect_equal(ncol(map_result), 5L)  # ID + 3 partitions + geometry
  expect_equal(class(map_result)[1], "sf")
  expect_true(all(c("K_2", "K_3", "K_4") %in% colnames(map_result)))
  
  # Verify colors exist in cluster object
  expect_true(!is.null(clu_colored$colors))
  expect_equal(length(clu_colored$colors), 3L)
  expect_true(all(c("K_2", "K_3", "K_4") %in% names(clu_colored$colors)))
  
  # Verify each color data frame has correct structure
  for(partition_name in c("K_2", "K_3", "K_4")) {
    color_df <- clu_colored$colors[[partition_name]]
    expect_true(is.data.frame(color_df))
    expect_equal(ncol(color_df), 2L)
    expect_true(all(c("cluster", "color") %in% colnames(color_df)))
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", color_df$color)))  # Valid hex colors
  }
  
  # Test actual plotting with colors (multi-panel: 3 partitions in 1x3 layout)
  expect_no_error({
    png(tempfile(fileext = ".png"), width = 1800, height = 600)
    map_bioregions(clusters = clu_colored, 
                   geometry = fishsf, 
                   plot = TRUE)
    dev.off()
  })
})

test_that("multi-panel bipartite with colors and plotting", {
  data(fishdf)
  data(fishsf)
  
  # Install binaries for infomap
  quietly(install_binaries(verbose = FALSE))
  
  # Create bipartite clusters using netclu_infomap (creates hierarchical partitions)
  clu_bip <- netclu_infomap(fishdf, bipartite = TRUE, seed = 123)
  clu_bip_colored <- bioregion_colors(clu_bip)
  
  # Verify infomap created multiple partitions (hierarchical)
  expect_equal(nrow(clu_bip$cluster_info), 2L)
  
  # Test with write_clusters to verify structure
  map_result <- map_bioregions(clusters = clu_bip_colored, 
                                geometry = fishsf, 
                                write_clusters = TRUE, 
                                plot = FALSE)
  
  # Verify structure - should only have site rows (not species)
  expect_equal(nrow(map_result), 338L)
  expect_equal(class(map_result)[1], "sf")
  
  # Should have ID + 2 partitions + geometry (infomap creates 2 hierarchical levels)
  expect_equal(ncol(map_result), 4L)
  
  # Verify partition columns exist (exclude ID and geometry columns)
  all_cols <- colnames(map_result)
  partition_cols <- setdiff(all_cols, c(colnames(fishsf)[1], attr(map_result, "sf_column")))
  expect_equal(length(partition_cols), 2L)
  
  # Verify colors exist in cluster object
  expect_true(!is.null(clu_bip_colored$colors))
  expect_equal(length(clu_bip_colored$colors), 2L)
  
  # Verify each color data frame has correct structure
  for(partition_name in names(clu_bip_colored$colors)) {
    color_df <- clu_bip_colored$colors[[partition_name]]
    expect_true(is.data.frame(color_df))
    expect_equal(ncol(color_df), 2L)
    expect_true(all(c("cluster", "color") %in% colnames(color_df)))
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", color_df$color)))  # Valid hex colors
  }
  
  # Verify that node_type attribute is preserved
  expect_true(!is.null(attr(clu_bip$clusters, "node_type")))
  
  # Verify that only site nodes are in the map result
  # (the original bipartite has both sites and species, but map should only show sites)
  expect_true(nrow(clu_bip$clusters) > nrow(map_result))  # Original has sites + species
  
  # Test actual plotting with colors (multi-panel)
  expect_no_error({
    png(tempfile(fileext = ".png"), width = 1200, height = 600)
    map_bioregions(clusters = clu_bip_colored, 
                   geometry = fishsf, 
                   plot = TRUE)
    dev.off()
  })
})

test_that("multi-panel unipartite layout selection", {
  data(fishmat)
  data(fishsf)
  
  dissim <- similarity_to_dissimilarity(
    similarity(fishmat, metric = "Simpson")
  )
  
  # Test with 2 partitions (should use 1x2 layout)
  clu_2 <- hclu_hierarclust(dissim, 
                            optimal_tree_method = "best", 
                            n_clust = c(2, 3),
                            verbose = FALSE)
  clu_2_colored <- bioregion_colors(clu_2)
  
  map_2 <- map_bioregions(clusters = clu_2_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_2), 4L)  # ID + 2 partitions + geometry
  expect_true(all(c("K_2", "K_3") %in% colnames(map_2)))
  
  # Test with 4 partitions (should use 1x4 layout)
  clu_4 <- hclu_hierarclust(dissim, 
                            optimal_tree_method = "best", 
                            n_clust = c(2, 3, 4, 5),
                            verbose = FALSE)
  clu_4_colored <- bioregion_colors(clu_4)
  
  map_4 <- map_bioregions(clusters = clu_4_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_4), 6L)  # ID + 4 partitions + geometry
  expect_true(all(c("K_2", "K_3", "K_4", "K_5") %in% colnames(map_4)))
  
  # Test with 6 partitions (should use 2-column layout: 3x2)
  clu_6 <- hclu_hierarclust(dissim, 
                            optimal_tree_method = "best", 
                            n_clust = 2:7,
                            verbose = FALSE)
  clu_6_colored <- bioregion_colors(clu_6)
  
  map_6 <- map_bioregions(clusters = clu_6_colored, 
                          geometry = fishsf, 
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_6), 8L)  # ID + 6 partitions + geometry
  expect_true(all(c("K_2", "K_3", "K_4", "K_5", "K_6", "K_7") %in% colnames(map_6)))
})

# Tests for bioregionalization parameter ---------------------------------------
test_that("bioregionalization parameter with integer indices", {
  
  # Create multi-partition clusters
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:5,
                               verbose = FALSE)
  
  # Test single partition by index
  map_1 <- map_bioregions(clusters = clu_hier, 
                          geometry = fishsf, 
                          bioregionalization = 1,
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_1), 3L)  # ID + 1 partition + geometry
  expect_true("K_2" %in% colnames(map_1))
  expect_false("K_3" %in% colnames(map_1))
  
  # Test multiple partitions by index
  map_2 <- map_bioregions(clusters = clu_hier, 
                          geometry = fishsf, 
                          bioregionalization = c(1, 3),
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_2), 4L)  # ID + 2 partitions + geometry
  expect_true(all(c("K_2", "K_4") %in% colnames(map_2)))
  expect_false(any(c("K_3", "K_5") %in% colnames(map_2)))
  
  # Test all partitions explicitly
  map_all <- map_bioregions(clusters = clu_hier, 
                            geometry = fishsf, 
                            bioregionalization = 1:4,
                            write_clusters = TRUE, 
                            plot = FALSE)
  expect_equal(ncol(map_all), 6L)  # ID + 4 partitions + geometry
  expect_true(all(c("K_2", "K_3", "K_4", "K_5") %in% colnames(map_all)))
})

test_that("bioregionalization parameter with character names", {
  
  # Create multi-partition clusters
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:5,
                               verbose = FALSE)
  
  # Test single partition by name
  map_1 <- map_bioregions(clusters = clu_hier, 
                          geometry = fishsf, 
                          bioregionalization = "K_3",
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_1), 3L)  # ID + 1 partition + geometry
  expect_true("K_3" %in% colnames(map_1))
  expect_false("K_2" %in% colnames(map_1))
  
  # Test multiple partitions by name
  map_2 <- map_bioregions(clusters = clu_hier, 
                          geometry = fishsf, 
                          bioregionalization = c("K_2", "K_4"),
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_2), 4L)  # ID + 2 partitions + geometry
  expect_true(all(c("K_2", "K_4") %in% colnames(map_2)))
  expect_false(any(c("K_3", "K_5") %in% colnames(map_2)))
  
  # Test non-contiguous selection
  map_3 <- map_bioregions(clusters = clu_hier, 
                          geometry = fishsf, 
                          bioregionalization = c("K_5", "K_2"),
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_3), 4L)
  expect_true(all(c("K_2", "K_5") %in% colnames(map_3)))
})

test_that("bioregionalization parameter with colored clusters", {
  
  # Create multi-partition clusters with colors
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:5,
                               verbose = FALSE)
  clu_colored <- bioregion_colors(clu_hier)
  
  # Test that colors are preserved when selecting partitions by index
  map_1 <- map_bioregions(clusters = clu_colored, 
                          geometry = fishsf, 
                          bioregionalization = c(1, 3),
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_1), 4L)  # ID + 2 partitions + geometry
  expect_true(all(c("K_2", "K_4") %in% colnames(map_1)))
  
  # Test that colors are preserved when selecting partitions by name
  map_2 <- map_bioregions(clusters = clu_colored, 
                          geometry = fishsf, 
                          bioregionalization = c("K_3", "K_5"),
                          write_clusters = TRUE, 
                          plot = FALSE)
  expect_equal(ncol(map_2), 4L)
  expect_true(all(c("K_3", "K_5") %in% colnames(map_2)))
  
  # Test plotting with selected partitions (should use custom colors)
  expect_no_error({
    png(tempfile(fileext = ".png"), width = 1200, height = 600)
    map_bioregions(clusters = clu_colored, 
                   geometry = fishsf, 
                   bioregionalization = c(2, 4),
                   plot = TRUE)
    dev.off()
  })
})

test_that("bioregionalization parameter NULL default", {
  
  # Create multi-partition clusters
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:4,
                               verbose = FALSE)
  
  # Test that NULL (default) plots all partitions
  map_null <- map_bioregions(clusters = clu_hier, 
                             geometry = fishsf, 
                             bioregionalization = NULL,
                             write_clusters = TRUE, 
                             plot = FALSE)
  expect_equal(ncol(map_null), 5L)  # ID + 3 partitions + geometry
  expect_true(all(c("K_2", "K_3", "K_4") %in% colnames(map_null)))
  
  # Test that omitting parameter is same as NULL
  map_default <- map_bioregions(clusters = clu_hier, 
                                geometry = fishsf, 
                                write_clusters = TRUE, 
                                plot = FALSE)
  expect_equal(ncol(map_default), 5L)
  expect_true(all(c("K_2", "K_3", "K_4") %in% colnames(map_default)))
  
  # Verify both approaches give same result
  expect_equal(colnames(map_null), colnames(map_default))
})

test_that("bioregionalization parameter error handling", {
  
  # Create multi-partition clusters
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:5,
                               verbose = FALSE)
  
  # Test invalid index (out of range)
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = 5,
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization indices must be between 1 and 4"
  )
  
  # Test invalid index (zero)
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = 0,
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization indices must be between 1 and 4"
  )
  
  # Test invalid index (negative)
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = -1,
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization indices must be between 1 and 4"
  )
  
  # Test non-integer numeric
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = 1.5,
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization must be an integer or a vector of integers"
  )
  
  # Test invalid name
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = "K_10",
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization name\\(s\\) not found in clusters"
  )
  
  # Test partially invalid names
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = c("K_2", "K_99"),
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization name\\(s\\) not found in clusters"
  )
  
  # Test invalid type
  expect_error(
    map_bioregions(clusters = clu_hier, 
                   geometry = fishsf, 
                   bioregionalization = TRUE,
                   write_clusters = TRUE, 
                   plot = FALSE),
    "bioregionalization must be NULL, an integer"
  )
})

test_that("bioregionalization with bipartite networks", {
  
  # Create bipartite network clusters
  net_bip <- mat_to_net(fishmat, weight = TRUE)
  clu_bip <- netclu_greedy(net_bip, bipartite = TRUE)
  
  # Note: netclu_greedy with bipartite typically creates one partition
  # But we can still test that the parameter works
  
  # Test with single partition
  map_bip <- map_bioregions(clusters = clu_bip, 
                            geometry = fishsf, 
                            bioregionalization = 1,
                            write_clusters = TRUE, 
                            plot = FALSE)
  expect_equal(nrow(map_bip), 338L)  # Only sites
  expect_equal(class(map_bip)[1], "sf")
  
  # Test that filtering works correctly (sites only, not species)
  # Get partition name
  partition_name <- colnames(clu_bip$clusters)[-1][1]
  expect_true(partition_name %in% colnames(map_bip))
})

test_that("bioregionalization preserves row order", {
  
  # Create multi-partition clusters
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  clu_hier <- hclu_hierarclust(dissim, 
                               optimal_tree_method = "best", 
                               n_clust = 2:4,
                               verbose = FALSE)
  
  # Get full map
  map_full <- map_bioregions(clusters = clu_hier, 
                             geometry = fishsf, 
                             write_clusters = TRUE, 
                             plot = FALSE)
  
  # Get partial map
  map_partial <- map_bioregions(clusters = clu_hier, 
                                geometry = fishsf, 
                                bioregionalization = c(1, 3),
                                write_clusters = TRUE, 
                                plot = FALSE)
  
  # Check that IDs are in same order
  expect_equal(map_full[[1]], map_partial[[1]])
  
  # Check that cluster values match for selected partitions
  expect_equal(map_full$K_2, map_partial$K_2)
  expect_equal(map_full$K_4, map_partial$K_4)
})
