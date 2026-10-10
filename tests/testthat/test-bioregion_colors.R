# Preamble code ---------------------------------------------------------------
comat_small <- matrix(sample(0:100, 40, replace = TRUE), 5, 8)
rownames(comat_small) <- paste0("Site", 1:5)
colnames(comat_small) <- paste0("Species", 1:8)

net_small <- similarity(comat_small, metric = "Simpson")
clust_small <- netclu_greedy(net_small)

net_fish <- similarity(fishmat, metric = "Simpson")
clust_fish <- netclu_greedy(net_fish)

net_bip <- mat_to_net(vegemat, weight = TRUE)
clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)

# Tests for bioregion_colors ---------------------------------------------------

test_that("bioregion_colors works with few clusters", {
  clust_colored <- bioregion_colors(clust_small)
  
  # Check that colors element exists and is a list
  expect_true(!is.null(clust_colored$colors))
  expect_true(inherits(clust_colored$colors, "list"))
  
  # Check that there's one element per partition
  nb_partitions <- ncol(clust_small$clusters) - 1
  expect_equal(length(clust_colored$colors), nb_partitions)
  
  # Check first partition's color data frame structure
  first_partition_colors <- clust_colored$colors[[1]]
  expect_equal(ncol(first_partition_colors), 2)
  expect_equal(names(first_partition_colors), c("cluster", "color"))
  
  # Check that all colors are hex codes
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", first_partition_colors$color)))
  
  # Check that cluster column is character
  expect_true(is.character(first_partition_colors$cluster))
  
  # Check that clusters_colors element exists
  expect_true(!is.null(clust_colored$clusters_colors))
  expect_true(inherits(clust_colored$clusters_colors, "data.frame"))
  
  # Check clusters_colors has same structure as clusters
  expect_equal(ncol(clust_colored$clusters_colors), ncol(clust_colored$clusters))
  expect_equal(names(clust_colored$clusters_colors), names(clust_colored$clusters))
  
  # Check that ID column is preserved
  expect_equal(clust_colored$clusters_colors[[1]], clust_colored$clusters[[1]])
  
  # Check that cluster columns contain hex colors
  expect_true(all(grepl("^#", clust_colored$clusters_colors[[2]][!is.na(clust_colored$clusters_colors[[2]])])))
})

test_that("bioregion_colors assigns grey to clusters beyond top 12", {
  clust_colored <- bioregion_colors(clust_fish)
  
  # Get first partition colors
  first_partition_colors <- clust_colored$colors[[1]]
  
  # Count unique clusters
  all_clusters <- unique(as.character(clust_fish$clusters[[2]]))
  all_clusters <- all_clusters[!is.na(all_clusters)]
  nb_clusters <- length(all_clusters)
  
  if (nb_clusters > 12) {
    # Check that we have vivid colors and grey shades
    vivid_palette <- rcartocolor::carto_pal(n = 12, name = "Vivid")
    
    # Count how many colors match vivid palette
    nb_vivid <- sum(first_partition_colors$color %in% vivid_palette)
    
    # Count grey shades (colors with equal RGB values, excluding black)
    is_grey <- sapply(first_partition_colors$color, function(col) {
      if (col == "#000000") return(FALSE)
      r <- substr(col, 2, 3)
      g <- substr(col, 4, 5)
      b <- substr(col, 6, 7)
      return(r == g && g == b)
    })
    nb_grey <- sum(is_grey)
    
    # We should have vivid colors
    expect_true(nb_vivid > 0)
    
    # If we have more than 12 clusters, we should have grey shades
    if (nb_clusters > 12) {
      expect_true(nb_grey > 0)
      expect_equal(nb_vivid + nb_grey, nb_clusters)
    }
  }
})

test_that("bioregion_colors applies cutoff correctly", {
  # Calculate cluster sizes
  cluster_sizes <- table(as.character(clust_fish$clusters[[2]]))
  
  # Find a reasonable cutoff (e.g., median)
  cutoff_val <- 3
  
  clust_colored <- bioregion_colors(clust_fish, cutoff_insignificant = cutoff_val)
  
  # Get first partition colors
  first_partition_colors <- clust_colored$colors[[1]]
  
  # Check that small clusters are colored black
  small_clusters <- names(cluster_sizes)[cluster_sizes <= cutoff_val]
  
  if (length(small_clusters) > 0) {
    black_clusters <- first_partition_colors$cluster[first_partition_colors$color == "#000000"]
    
    # All small clusters should be black
    expect_true(all(small_clusters %in% black_clusters))
    
    # Count black clusters
    expect_equal(length(black_clusters), length(small_clusters))
  }
})

test_that("cluster_ordering parameter works correctly", {
  clust_colored <- bioregion_colors(clust_fish, cluster_ordering = "n_sites")
  
  # Check that it runs without error
  expect_true(!is.null(clust_colored$colors))
  
  # Get first partition colors
  first_partition_colors <- clust_colored$colors[[1]]
  
  # Get cluster sizes
  cluster_sizes <- table(as.character(clust_fish$clusters[[2]]))
  largest_cluster <- names(cluster_sizes)[which.max(cluster_sizes)]
  
  # Get vivid palette
  nb_clusters <- length(unique(as.character(clust_fish$clusters[[2]])))
  nb_colors <- min(12, nb_clusters)
  vivid_palette <- rcartocolor::carto_pal(n = nb_colors, name = "Vivid")
  
  # The largest cluster should get the first vivid color
  largest_color <- first_partition_colors$color[first_partition_colors$cluster == largest_cluster]
  
  # It should be a vivid color (in the palette)
  expect_true(largest_color %in% vivid_palette)
})

test_that("n_species and n_both ordering work for bipartite", {
  # Test n_species ordering
  expect_error({
    clust_species <- bioregion_colors(clust_bip, cluster_ordering = "n_species")
  }, NA)
  
  # Test n_both ordering
  expect_error({
    clust_both <- bioregion_colors(clust_bip, cluster_ordering = "n_both")
  }, NA)
  
  expect_true(!is.null(clust_both$colors))
})

test_that("error thrown for bipartite ordering on non-bipartite data", {
  expect_error(
    bioregion_colors(clust_fish, cluster_ordering = "n_species"),
    "can only be used with bipartite clustering"
  )
  
  expect_error(
    bioregion_colors(clust_fish, cluster_ordering = "n_both"),
    "can only be used with bipartite clustering"
  )
})

test_that("cluster IDs treated as character throughout", {
  clust_colored <- bioregion_colors(clust_fish)
  
  # Get first partition colors
  first_partition_colors <- clust_colored$colors[[1]]
  
  # Check that cluster column in colors is character
  expect_true(is.character(first_partition_colors$cluster))
  
  # Check that all cluster IDs are character (not factor or numeric)
  expect_false(is.factor(first_partition_colors$cluster))
  expect_false(is.numeric(first_partition_colors$cluster))
  
  # Check that clusters_colors contains hex codes (character)
  expect_true(is.character(clust_colored$clusters_colors[[2]]))
  
  # Check that hex codes are properly formatted
  non_na_colors <- clust_colored$clusters_colors[[2]][!is.na(clust_colored$clusters_colors[[2]])]
  expect_true(all(grepl("^#", non_na_colors)))
})

test_that("different rcartocolor palettes work", {
  palettes <- c("Bold", "Prism", "Safe", "Pastel")
  
  for (pal in palettes) {
    expect_error({
      clust_colored <- bioregion_colors(clust_small, palette = pal)
      expect_true(!is.null(clust_colored$colors))
    }, NA, info = paste("Failed with palette:", pal))
  }
})

test_that("grey shades have sufficient contrast", {
  # Create clustering with many clusters
  # If fishmat doesn't have enough clusters, this test may need adjustment
  clust_colored <- bioregion_colors(clust_fish)
  
  # Get first partition colors
  first_partition_colors <- clust_colored$colors[[1]]
  
  # Extract grey shades
  is_grey <- sapply(first_partition_colors$color, function(col) {
    if (col == "#000000") return(FALSE)
    r <- substr(col, 2, 3)
    g <- substr(col, 4, 5)
    b <- substr(col, 6, 7)
    return(r == g && g == b)
  })
  
  grey_colors <- first_partition_colors$color[is_grey]
  
  if (length(grey_colors) > 1) {
    # Convert to numeric values
    grey_values <- sapply(grey_colors, function(col) {
      strtoi(substr(col, 2, 3), base = 16)
    })
    
    # Check that there's variation in grey values
    expect_true(length(unique(grey_values)) > 1 || length(grey_values) == 1)
    
    # Check that greys are within expected range (64 to 204)
    expect_true(all(grey_values >= 64 & grey_values <= 204))
  }
})

test_that("print.bioregion.clusters displays color info", {
  clust_colored <- bioregion_colors(clust_fish)
  
  # Capture output
  output <- capture.output(print(clust_colored))
  output_text <- paste(output, collapse = "\n")
  
  # Check that color information is present
  expect_true(any(grepl("Color palette", output_text)))
})

test_that("input validation catches errors", {
  # Test with non-bioregion.clusters object
  expect_error(
    bioregion_colors(data.frame(x = 1:10)),
    "must be a bioregion.clusters object"
  )
  
  # Test with invalid palette
  expect_error(
    bioregion_colors(clust_fish, palette = "InvalidPalette"),
    "Invalid palette name"
  )
  
  # Test with invalid cluster_ordering
  expect_error(
    bioregion_colors(clust_fish, cluster_ordering = "invalid"),
    "must be one of"
  )
  
  # Test with invalid cutoff_insignificant
  expect_error(
    bioregion_colors(clust_fish, cutoff_insignificant = "text")
  )
})

test_that("handles multiple partitions correctly", {
  dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
  
  quietly(
    clust_hier <- hclu_hierarclust(dissim, n_clust = c(3, 4, 8),
                                   optimal_tree_method = "best")
  )
  
  clust_colored <- bioregion_colors(clust_hier)
  
  # Check that colors are assigned
  expect_true(!is.null(clust_colored$colors))
  
  # Check that colors is a list
  expect_true(inherits(clust_colored$colors, "list"))
  
  # Check that there's one color set per partition
  nb_partitions <- ncol(clust_hier$clusters) - 1
  expect_equal(length(clust_colored$colors), nb_partitions)
  
  # Check each partition has appropriate number of colors
  for (i in seq_len(nb_partitions)) {
    partition_name <- names(clust_hier$clusters)[i + 1]
    partition_clusters <- unique(as.character(clust_hier$clusters[[i + 1]]))
    partition_clusters <- partition_clusters[!is.na(partition_clusters)]
    
    expect_equal(nrow(clust_colored$colors[[partition_name]]), length(partition_clusters))
  }
  
  # Check that clusters_colors has same number of columns
  expect_equal(ncol(clust_colored$clusters_colors), ncol(clust_colored$clusters))
})

test_that("bioregion_colors preserves original cluster object structure", {
  clust_colored <- bioregion_colors(clust_fish)
  
  # Check that original elements are preserved
  expect_equal(clust_colored$name, clust_fish$name)
  expect_equal(clust_colored$args, clust_fish$args)
  expect_equal(clust_colored$inputs, clust_fish$inputs)
  expect_equal(clust_colored$clusters, clust_fish$clusters)
  
  # Check that new elements are added
  expect_true(!is.null(clust_colored$colors))
  expect_true(!is.null(clust_colored$clusters_colors))
  
  # Check that class is preserved
  expect_true(inherits(clust_colored, "bioregion.clusters"))
})

test_that("color assignment is deterministic", {
  # Run twice and check results are identical
  clust_colored1 <- bioregion_colors(clust_fish)
  clust_colored2 <- bioregion_colors(clust_fish)
  
  expect_equal(clust_colored1$colors, clust_colored2$colors)
  expect_equal(clust_colored1$clusters_colors, clust_colored2$clusters_colors)
})

test_that("handles clusters with NA values", {
  # Create a cluster object with some NAs
  clust_na <- clust_fish
  clust_na$clusters[[2]][1:2] <- NA
  
  expect_error({
    clust_colored <- bioregion_colors(clust_na)
    expect_true(!is.null(clust_colored$colors))
  }, NA)
})
