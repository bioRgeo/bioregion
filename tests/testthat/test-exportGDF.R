# Inputs -----------------------------------------------------------------------

# Unipartite network without weights
net_uni <- data.frame(
  Node1 = c("A", "A", "B", "C", "D"),
  Node2 = c("B", "C", "C", "D", "E")
)

# Unipartite network with weights
net_uni_w <- data.frame(
  Node1 = c("A", "A", "B", "C", "D"),
  Node2 = c("B", "C", "C", "D", "E"),
  Weight = c(1.5, 2.0, 1.0, 3.5, 0.5)
)

# Bipartite network without weights
net_bip <- data.frame(
  Site = c(rep("Site1", 2), rep("Site2", 3), rep("Site3", 2)),
  Species = c("Sp_a", "Sp_b", "Sp_a", "Sp_c", "Sp_d", "Sp_b", "Sp_d")
)

# Bipartite network with weights
net_bip_w <- data.frame(
  Site = c(rep("Site1", 2), rep("Site2", 3), rep("Site3", 2)),
  Species = c("Sp_a", "Sp_b", "Sp_a", "Sp_c", "Sp_d", "Sp_b", "Sp_d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

# Node characteristics without colors
node_chars_no_color <- data.frame(
  node_id = c("A", "B", "C", "D", "E"),
  category = c("Type1", "Type1", "Type2", "Type2", "Type1"),
  value = c(10, 20, 15, 25, 30)
)

# Node characteristics with colors
node_chars_color <- data.frame(
  node_id = c("A", "B", "C", "D", "E"),
  category = c("Type1", "Type1", "Type2", "Type2", "Type1"),
  node_color = c("#FF5733", "#33FF57", "#3357FF", "#FF33F5", "#F5FF33")
)

# Invalid inputs for testing
net_invalid1 <- data.frame(
  Node1 = c("A", "A", "B"),
  Node2 = c("B", "C", "C"),
  Weight = c("a", "b", "c")  # Non-numeric weight
)

net_invalid2 <- data.frame(
  Node1 = c("A", "A", "B"),
  Node2 = c("B", "C", "C"),
  Weight = c(1.5, NA, 1.0)  # NA in weight
)

net_invalid3 <- data.frame(
  Node1 = c("A", "A", "A"),
  Node2 = c("B", "C", "B")  # Duplicated pairs
)

net_invalid4 <- data.frame(
  Node1 = c("A", NA, "B"),
  Node2 = c("B", "C", "C")  # NA in nodes
)

node_chars_invalid <- data.frame(
  node_id = c("A", NA, "C"),  # NA in data.frame
  category = c("Type1", "Type1", "Type2")
)

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs - file creation", {
  
  # Test unipartite network without weights
  temp_file1 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni, file = temp_file1)
  expect_true(file.exists(temp_file1))
  unlink(temp_file1)
  
  # Test unipartite network with weights
  temp_file2 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", file = temp_file2)
  expect_true(file.exists(temp_file2))
  unlink(temp_file2)
  
  # Test bipartite network without weights
  temp_file3 <- tempfile(fileext = ".gdf")
  exportGDF(net_bip, col1 = "Site", col2 = "Species", file = temp_file3)
  expect_true(file.exists(temp_file3))
  unlink(temp_file3)
  
  # Test bipartite network with weights
  temp_file4 <- tempfile(fileext = ".gdf")
  exportGDF(net_bip_w, col1 = "Site", col2 = "Species", 
            weight = "Weight", file = temp_file4)
  expect_true(file.exists(temp_file4))
  unlink(temp_file4)
  
  # Test with node characteristics (no colors)
  temp_file5 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", 
            bioregions = node_chars_no_color,
            bioregionalization = "node_id",
            file = temp_file5)
  expect_true(file.exists(temp_file5))
  unlink(temp_file5)
  
  # Test with node characteristics (with colors)
  temp_file6 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", 
            bioregions = node_chars_color,
            bioregionalization = "node_id",
            color_column = "node_color",
            file = temp_file6)
  expect_true(file.exists(temp_file6))
  unlink(temp_file6)
})

test_that("valid outputs - file content structure", {
  
  # Test basic structure
  temp_file <- tempfile(fileext = ".gdf")
  exportGDF(net_uni, file = temp_file)
  content <- readLines(temp_file)
  
  # Check for nodedef and edgedef headers
  expect_true(any(grepl("^nodedef>", content)))
  expect_true(any(grepl("^edgedef>", content)))
  
  # Check number of nodes (5 unique nodes: A, B, C, D, E)
  nodedef_idx <- which(grepl("^nodedef>", content))
  edgedef_idx <- which(grepl("^edgedef>", content))
  num_nodes <- edgedef_idx - nodedef_idx - 1
  expect_equal(num_nodes, 5)
  
  # Check number of edges (5 edges)
  num_edges <- length(content) - edgedef_idx
  expect_equal(num_edges, 5)
  
  unlink(temp_file)
})

test_that("valid outputs - weighted vs unweighted", {
  
  # Unweighted network
  temp_file1 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni, file = temp_file1)
  content1 <- readLines(temp_file1)
  edgedef_line1 <- content1[grepl("^edgedef>", content1)]
  expect_false(grepl("weight", edgedef_line1))
  unlink(temp_file1)
  
  # Weighted network
  temp_file2 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", file = temp_file2)
  content2 <- readLines(temp_file2)
  edgedef_line2 <- content2[grepl("^edgedef>", content2)]
  expect_true(grepl("weight", edgedef_line2))
  unlink(temp_file2)
})

test_that("valid outputs - node attributes", {
  
  # Without node characteristics
  temp_file1 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni, file = temp_file1)
  content1 <- readLines(temp_file1)
  nodedef_line1 <- content1[grepl("^nodedef>", content1)]
  # Should only have name and label
  expect_true(grepl("name VARCHAR,label VARCHAR", nodedef_line1))
  expect_false(grepl("category", nodedef_line1))
  unlink(temp_file1)
  
  # With node characteristics (no colors)
  temp_file2 <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", 
            bioregions = node_chars_no_color,
            bioregionalization = "node_id",
            file = temp_file2)
  content2 <- readLines(temp_file2)
  nodedef_line2 <- content2[grepl("^nodedef>", content2)]
  # Should have category and value attributes
  expect_true(grepl("category", nodedef_line2))
  expect_true(grepl("value", nodedef_line2))
  unlink(temp_file2)
})

test_that("valid outputs - color conversion", {
  
  # With colors
  temp_file <- tempfile(fileext = ".gdf")
  exportGDF(net_uni_w, weight = "Weight", 
            bioregions = node_chars_color,
            bioregionalization = "node_id",
            color_column = "node_color",
            file = temp_file)
  content <- readLines(temp_file)
  
  # Check for color in header
  nodedef_line <- content[grepl("^nodedef>", content)]
  expect_true(grepl("color", nodedef_line))
  
  # Check that hex colors are converted to RGB
  # #FF5733 should become 255,87,51
  expect_true(any(grepl("255,87,51", content)))
  
  unlink(temp_file)
})

test_that("valid outputs - bipartite networks", {
  
  # Bipartite network without weights
  temp_file1 <- tempfile(fileext = ".gdf")
  exportGDF(net_bip, col1 = "Site", col2 = "Species", file = temp_file1)
  content1 <- readLines(temp_file1)
  
  # Check for both site and species nodes (3 sites + 4 species = 7 nodes)
  nodedef_idx <- which(grepl("^nodedef>", content1))
  edgedef_idx <- which(grepl("^edgedef>", content1))
  num_nodes <- edgedef_idx - nodedef_idx - 1
  expect_equal(num_nodes, 7)
  
  # Check number of edges (7 edges)
  num_edges <- length(content1) - edgedef_idx
  expect_equal(num_edges, 7)
  
  unlink(temp_file1)
  
  # Bipartite network with weights
  temp_file2 <- tempfile(fileext = ".gdf")
  exportGDF(net_bip_w, col1 = "Site", col2 = "Species", 
            weight = "Weight", file = temp_file2)
  content2 <- readLines(temp_file2)
  edgedef_line2 <- content2[grepl("^edgedef>", content2)]
  expect_true(grepl("weight", edgedef_line2))
  
  unlink(temp_file2)
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs - df parameter", {
  
  expect_error(
    exportGDF("not a dataframe"),
    "must be a data.frame",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(data.frame(x = 1:3)),
    "must be a data.frame with at least two columns",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_invalid3),
    "contain duplicated pairs of nodes",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_invalid4),
    "NA(s) detected",
    fixed = TRUE
  )
})

test_that("invalid inputs - col1 and col2 parameters", {
  
  expect_error(
    exportGDF(net_uni, col1 = 123),
    "col1 must be a character",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, col2 = TRUE),
    "col2 must be a character",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, col1 = "NonExistent"),
    "col1 ('NonExistent') is not a column name in df",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, col2 = "NonExistent"),
    "col2 ('NonExistent') is not a column name in df",
    fixed = TRUE
  )
})

test_that("invalid inputs - weight parameter", {
  
  expect_error(
    exportGDF(net_uni_w, weight = 123),
    "weight must be a character",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni_w, weight = "NonExistent"),
    "weight ('NonExistent') is not a column name in df",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_invalid1, weight = "Weight"),
    "The weight column must be numeric",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_invalid2, weight = "Weight"),
    "NA(s) detected in the weight column",
    fixed = TRUE
  )
})

test_that("invalid inputs - bioregions parameter", {
  
  expect_error(
    exportGDF(net_uni, bioregions = "not a dataframe"),
    "must be a data.frame",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, bioregions = node_chars_invalid),
    "NA(s) detected in the data.frame",
    fixed = TRUE
  )
})

test_that("invalid inputs - bioregionalization parameter", {
  
  expect_error(
    exportGDF(net_uni, bioregions = node_chars_no_color, 
              bioregionalization = 123),
    "bioregionalization must be a character",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, bioregions = node_chars_no_color, 
              bioregionalization = "NonExistent"),
    "bioregionalization ('NonExistent') is not a column name in bioregions",
    fixed = TRUE
  )
})

test_that("invalid inputs - color_column parameter", {
  
  expect_error(
    exportGDF(net_uni, bioregions = node_chars_color,
              bioregionalization = "node_id",
              color_column = 123),
    "color_column must be a character",
    fixed = TRUE
  )
  
  expect_error(
    exportGDF(net_uni, bioregions = node_chars_no_color,
              bioregionalization = "node_id",
              color_column = "NonExistent"),
    "Color column 'NonExistent' not found in node characteristics",
    fixed = TRUE
  )
})

test_that("invalid inputs - file parameter", {
  
  expect_error(
    exportGDF(net_uni, file = 123),
    "file must be a character",
    fixed = TRUE
  )
})

# Tests for zero-weight filtering ----------------------------------------------
test_that("zero-weight edges are filtered with warning", {
  
  # Create a network with zero-weight edges
  net_with_zeros <- data.frame(
    Node1 = c("A", "A", "B", "C", "D"),
    Node2 = c("B", "C", "C", "D", "E"),
    Weight = c(1.5, 0, 1.0, 0, 0.5)  # 2 edges with weight = 0
  )
  
  temp_file <- tempfile(fileext = ".gdf")
  
  # Should produce a warning
  expect_warning(
    exportGDF(net_with_zeros, weight = "Weight", file = temp_file),
    "2 edge\\(s\\) with weight = 0 detected and removed"
  )
  
  # Check that file was created
  expect_true(file.exists(temp_file))
  
  # Read the file and check that only 3 edges remain (out of 5 original)
  content <- readLines(temp_file)
  edgedef_idx <- which(grepl("^edgedef>", content))
  num_edges <- length(content) - edgedef_idx
  expect_equal(num_edges, 3)
  
  # Check that zero-weight edges are not in the file
  expect_false(any(grepl(",0$", content)))
  
  unlink(temp_file)
})

test_that("no warning when no zero-weight edges", {
  
  temp_file <- tempfile(fileext = ".gdf")
  
  # Should NOT produce a warning
  expect_no_warning(
    exportGDF(net_uni_w, weight = "Weight", file = temp_file)
  )
  
  unlink(temp_file)
})

# Tests for bioregion.clusters integration -------------------------------------
test_that("bioregion.clusters integration - basic functionality", {
  
  # Create a simple similarity network and clustering
  comat_small <- matrix(sample(0:100, 40, replace = TRUE), 5, 8)
  rownames(comat_small) <- c("A", "B", "C", "D", "E")
  colnames(comat_small) <- paste0("Species", 1:8)
  
  net <- similarity(comat_small, metric = "Simpson")
  clust <- netclu_greedy(net)
  clust_colored <- bioregion_colors(clust)
  
  # Create network data frame
  net_df <- data.frame(
    Node1 = c("A", "A", "B", "C", "D"),
    Node2 = c("B", "C", "C", "D", "E"),
    Weight = c(1.5, 2.0, 1.0, 3.5, 0.5)
  )
  
  temp_file <- tempfile(fileext = ".gdf")
  
  # Export with bioregion.clusters object
  expect_message(
    exportGDF(net_df, weight = "Weight", 
              bioregions = clust_colored,
              file = temp_file),
    "Using partition"
  )
  
  expect_true(file.exists(temp_file))
  
  # Read content and verify structure
  content <- readLines(temp_file)
  
  # Check for cluster attribute
  nodedef_line <- content[grepl("^nodedef>", content)]
  expect_true(grepl("cluster", nodedef_line))
  
  # Check for color attribute
  expect_true(grepl("color", nodedef_line))
  
  # Verify that RGB colors are present (format: '255,87,51')
  node_lines <- content[2:(which(grepl("^edgedef>", content)) - 1)]
  expect_true(any(grepl("'[0-9]+,[0-9]+,[0-9]+'", node_lines)))
  
  unlink(temp_file)
})

test_that("bioregion.clusters integration - multiple partitions", {
  
  # Create clustering with multiple partitions
  dissim <- similarity_to_dissimilarity(
    similarity(fishmat, metric = "Simpson")
  )
  # Use n_clust values that can actually be achieved
  clust_hier <- hclu_hierarclust(dissim, n_clust = c(3, 4, 8))
  clust_colored <- bioregion_colors(clust_hier)
  
  # Create simple network from fishmat
  net_df <- mat_to_net(fishmat[1:10, 1:20], weight = TRUE)
  
  temp_file1 <- tempfile(fileext = ".gdf")
  temp_file2 <- tempfile(fileext = ".gdf")
  
  # Export with default partition (first one)
  expect_message(
    exportGDF(net_df, weight = "Weight",
              bioregions = clust_colored,
              file = temp_file1),
    "Using partition"
  )
  
  # Export with specific partition
  expect_no_message(
    exportGDF(net_df, weight = "Weight",
              bioregions = clust_colored,
              cluster_column = "K_8",
              file = temp_file2)
  )
  
  expect_true(file.exists(temp_file1))
  expect_true(file.exists(temp_file2))
  
  # Verify different partitions produce different cluster assignments
  content1 <- readLines(temp_file1)
  content2 <- readLines(temp_file2)
  
  # Files should be different (different cluster assignments)
  expect_false(identical(content1, content2))
  
  unlink(temp_file1)
  unlink(temp_file2)
})

test_that("bioregion.clusters integration - without colors", {
  
  # Create clustering without colors
  comat_small <- matrix(sample(0:100, 40, replace = TRUE), 5, 8)
  rownames(comat_small) <- c("A", "B", "C", "D", "E")
  colnames(comat_small) <- paste0("Species", 1:8)
  
  net <- similarity(comat_small, metric = "Simpson")
  clust <- netclu_greedy(net)
  # Don't add colors
  
  # Create network data frame
  net_df <- data.frame(
    Node1 = c("A", "A", "B", "C", "D"),
    Node2 = c("B", "C", "C", "D", "E"),
    Weight = c(1.5, 2.0, 1.0, 3.5, 0.5)
  )
  
  temp_file <- tempfile(fileext = ".gdf")
  
  # Export with bioregion.clusters object without colors
  expect_message(
    exportGDF(net_df, weight = "Weight",
              bioregions = clust,
              file = temp_file),
    "Using partition"
  )
  
  expect_true(file.exists(temp_file))
  
  # Read content and verify structure
  content <- readLines(temp_file)
  
  # Check for cluster attribute
  nodedef_line <- content[grepl("^nodedef>", content)]
  expect_true(grepl("cluster", nodedef_line))
  
  # Should NOT have color attribute
  expect_false(grepl("color", nodedef_line))
  
  unlink(temp_file)
})

test_that("bioregion.clusters integration - invalid cluster_column", {
  
  # Create clustering
  comat_small <- matrix(sample(0:100, 40, replace = TRUE), 5, 8)
  rownames(comat_small) <- c("A", "B", "C", "D", "E")
  colnames(comat_small) <- paste0("Species", 1:8)
  
  net <- similarity(comat_small, metric = "Simpson")
  clust <- netclu_greedy(net)
  clust_colored <- bioregion_colors(clust)
  
  # Create network data frame
  net_df <- data.frame(
    Node1 = c("A", "A", "B", "C", "D"),
    Node2 = c("B", "C", "C", "D", "E")
  )
  
  temp_file <- tempfile(fileext = ".gdf")
  
  # Try to use non-existent partition
  expect_error(
    exportGDF(net_df,
              bioregions = clust_colored,
              cluster_column = "NonExistent",
              file = temp_file),
    "cluster_column 'NonExistent' not found in available partitions"
  )
})

test_that("bioregion.clusters integration - color consistency", {
  
  # Create clustering with colors
  net <- similarity(fishmat, metric = "Simpson")
  clust <- netclu_greedy(net)
  clust_colored <- bioregion_colors(clust)
  
  # Get the first partition name
  partition_name <- names(clust_colored$colors)[1]
  
  # Create network from fishmat
  net_df <- mat_to_net(fishmat[1:10, 1:20], weight = TRUE)
  
  temp_file <- tempfile(fileext = ".gdf")
  
  suppressMessages(
    exportGDF(net_df, weight = "Weight",
              bioregions = clust_colored,
              file = temp_file)
  )
  
  content <- readLines(temp_file)
  
  # Extract node lines
  nodedef_idx <- which(grepl("^nodedef>", content))
  edgedef_idx <- which(grepl("^edgedef>", content))
  node_lines <- content[(nodedef_idx + 1):(edgedef_idx - 1)]
  
  # Each node line should have a color in RGB format
  # Count lines with color patterns
  lines_with_colors <- sum(grepl("'[0-9]+,[0-9]+,[0-9]+'", node_lines))
  
  # All nodes should have colors
  expect_equal(lines_with_colors, length(node_lines))
  
  unlink(temp_file)
})
