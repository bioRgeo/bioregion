IHCT <- function(dist_mat,
                 sites = rownames(dist_mat),
                 method = "average",
                 n_runs = 100,
                 n_splits = 2,
                 depth = 1, 
                 previous_height = Inf, 
                 max_remaining_size = length(sites),
                 monotonicity_direction = c("top-down", "bottom-up"),
                 verbose = TRUE) {
  
  # n_runs & n_splits
  if(n_splits > n_runs){
    n_splits <- n_runs
  }

  # Progress info because the algorithm is quite long to run
  if (verbose && interactive()) {
    if(max_remaining_size > 2) {
      cat(sprintf(
        "\rProcessing height: %.3f | Current branch size: %d | Max cluster size remaining: %d     ", 
        previous_height,
        length(sites),
        max(unlist(max_remaining_size))
      ))
    } else if(max_remaining_size == 2) {
      cat(sprintf(
        "\rProcessing height: %.3f | Current branch size: %d | All clusters successfully processed     ", 
        previous_height,
        length(sites)
      ))
    }
    
    utils::flush.console()  # Ensures the message is displayed in real-time
  }
  
  # Base case: If the subcluster is a single-site cluster, return it as a leaf node
  if (length(sites) == 1) {
    return(list(name = sites, height = 0, depth = depth, children = NULL))
  }
  
  # Step 1: subset dissimilarity matrix to current cluster only 
  dist_mat_d <- dist_mat[sites,sites]

  # Step 2: find the majority vote for a binary split among many trees
  if(nrow(dist_mat_d) > 2) { # Only if we have more than 2 sites
    subclusters <- coassign_binary_split(dist_mat_d, 
                                         method = method,
                                         n_runs = n_runs,
                                         n_splits = n_splits,
                                         binsplit = "tree")
  } else { # if 2 sites only we already have the subclusters
    subclusters <- list(rownames(dist_mat_d)[1],
                        rownames(dist_mat_d)[2])
  }
  
  # Ensure the smallest subcluster is the first (same rules as for hclust)
  if (length(subclusters[[1]]) > length(subclusters[[2]])) {
    subclusters <- list(subclusters[[2]], subclusters[[1]])
  }
  
  # Calculate height based on the specified method
  cluster1_sites <- subclusters[[1]]
  cluster2_sites <- subclusters[[2]]
  pairwise_distances <- dist_mat_d[cluster1_sites, cluster2_sites]
  
  # Calculate centroids for ward.D and ward.D2
  centroid1 <- colMeans(dist_mat_d[cluster1_sites, , drop = FALSE])
  centroid2 <- colMeans(dist_mat_d[cluster2_sites, , drop = FALSE])
  
  # Compute squared Euclidean distance between centroids
  centroid_distance <- sum((centroid1 - centroid2)^2)

  # Compute height based on the selected method
  calculated_height <- switch(
    method,
    "single" = min(pairwise_distances),
    "complete" = max(pairwise_distances),
    "average" = mean(pairwise_distances),
    "mcquitty" = mean(pairwise_distances),  # Note: This is acceptable for initial height
    "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * sqrt(centroid_distance),
    "ward.D2" = (length(cluster1_sites) * length(cluster2_sites)) /  
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "centroid" = centroid_distance, 
    "median" = 0.5 * centroid_distance,  
    stop("method argument is not valid")
  )
  
  # Step 3: height constraint (monotonic constraint of trees)
  if (monotonicity_direction == "top-down") {
    height <- min(calculated_height, previous_height)
  } else if (monotonicity_direction == "bottom-up") {
    # Height will be adjusted after processing subclusters
    height <- calculated_height
  }
  if (height < 0) {
    height <- 0
  } 
  
  # If we are at the iteration where we work on the max cluster size,
  # then we update max_remaining_size based on the largest subcluster size
  # (always in position 2)
  if (length(sites) == max_remaining_size) {
    max_remaining_size <- length(subclusters[[2]])
  }
  
  # Process each subcluster and collect their structures
  children <- list()
  for (i in 1:2) {
    subcluster <- subclusters[[i]]
    child_structure <- IHCT(dist_mat_d,
                            sites = subcluster,
                            method = method,
                            n_runs = n_runs,
                            n_splits = n_splits,
                            depth = depth + 1,
                            previous_height = height,
                            max_remaining_size = max_remaining_size,
                            monotonicity_direction = monotonicity_direction,
                            verbose = verbose)
    children[[i]] <- child_structure
  }
  
  # After processing subclusters, adjust height if monotonicity is bottom-up
  if (monotonicity_direction == "bottom-up") {
    max_child_height <- max(sapply(children, function(x) x$height))
    if (max_child_height > calculated_height) {
      height <- max_child_height
    } else {
      height <- calculated_height
    }
  }
  
  # Return the tree structure for the current split
  return(list(name = paste(sites, collapse = ","),
              height = height,
              depth = depth,
              children = children))
}

# Function to find the majority vote for a binary split among many trees
coassign_binary_split <- function(dist_mat,
                                  method = "average", 
                                  n_runs = 100,
                                  n_splits = 2,
                                  binsplit = "tree") {
  
  # Number of sites
  n <- dim(dist_mat)[1]
  
  # Step 1: randomize distance matrices and generate trees
  trees <- list()
  ccc <- NULL
  for (run in 1:n_runs) {
    trees[[run]] <- list()
    rand_dist_mat <- stats::as.dist(randomize_dist(dist_mat))
    trees[[run]] <- fastcluster::hclust(rand_dist_mat, 
                                        method = method)
    ccc <- c(ccc, suppressWarnings(tree_eval(trees[[run]], 
                                             dist_mat)$cophcor))
  }
  trees <- trees[order(ccc, decreasing = TRUE)[1:n_splits]]
  
  # Step 2: cut trees at 2 clusters for each run
  clusts <- lapply(trees, function(tree) {cut_tree(stats::as.hclust(tree), 
                                                   n_clust = 2, 
                                                   find_h = FALSE,
                                                   verbose = FALSE)})
  
  # Step 3: make a data.frame with all partitions based on sorted ID
  fixed_order <- sort(clusts[[1]]$ID) 
  clusts <- do.call(cbind,
                    lapply(clusts, function(df) {df[match(fixed_order, df$ID), 
                                                    2, 
                                                    drop = FALSE]}))
  
  # Step 4: compute pairwise dissimilarity based on memberships
  coassign <- matrix(0, n, n)
  for (r in 1:dim(clusts)[2]) {
    coassign <- coassign + (outer(clusts[, r], clusts[, r], "=="))
  }
  coassign <- (1 - coassign / n_runs)
  rownames(coassign) <- fixed_order
  colnames(coassign) <- fixed_order
  
  coassign <- stats::as.dist(coassign)
  
  # Step 5: Split into two clusters form the coassign matrix
  if (binsplit == "tree") {
    pw_tree <- fastcluster::hclust(coassign, method = method)
    groups <- stats::cutree(pw_tree, k = 2)
  } else if (binsplit == "pam") {
    groups <- cluster::pam(coassign, k = 2)$clustering
  }
  groups <- data.frame(Site = names(groups), cluster = groups)

  # Return the two groups of sites
  return(split(groups$Site, groups$cluster))
}


reconstruct_hclust_bis <- function(tree) {
  # Step 1: Assign negative IDs to leaf nodes based on alphabetical order
  get_leaves <- function(node) {
    if (is.null(node$children)) {
      return(node$name)
    } else {
      return(unlist(lapply(node$children, get_leaves)))
    }
  }
  
  leaf_labels <- sort(unique(get_leaves(tree)))
  n_leaves <- length(leaf_labels)
  
  # Map leaf labels to negative IDs based on alphabetical order
  leaf_ids <- stats::setNames(-seq_along(leaf_labels), leaf_labels)
  
  # Step 2: Traverse the tree and collect internal nodes
  internal_nodes <- list()
  node_counter <- 0
  
  collect_internal_nodes <- function(node) {
    if (is.null(node$children)) {
      # Leaf node
      node_id <- leaf_ids[node$name]
      return(node_id)
    } else {
      # Internal node
      left_id <- collect_internal_nodes(node$children[[1]])
      right_id <- collect_internal_nodes(node$children[[2]])
      
      node_counter <<- node_counter + 1
      temp_id <- paste0("n", node_counter)
      internal_nodes[[temp_id]] <<- list(
        left = left_id,
        right = right_id,
        height = node$height
      )
      return(temp_id)
    }
  }
  
  # Collect all internal nodes
  root_id <- collect_internal_nodes(tree)
  
  # Now, we have internal_nodes as a list with temporary IDs, and their left, right, height
  # We need to assign actual IDs to internal nodes, in order of increasing height
  
  # First, create a data frame of internal nodes
  internal_nodes_df <- data.frame(
    temp_id = names(internal_nodes),
    left = sapply(internal_nodes, function(x) x$left),
    right = sapply(internal_nodes, function(x) x$right),
    height = sapply(internal_nodes, function(x) x$height),
    stringsAsFactors = FALSE
  )
  
  # Sort internal nodes by height (increasing order)
  internal_nodes_df <- internal_nodes_df[order(internal_nodes_df$height, decreasing = FALSE), ]
  
  # Assign IDs to internal nodes starting from 1
  internal_node_ids <- seq_len(nrow(internal_nodes_df))
  names(internal_node_ids) <- internal_nodes_df$temp_id
  
  # Function to map IDs
  map_id <- function(x) {
    if (x %in% names(internal_node_ids)) {
      return(internal_node_ids[x])
    } else {
      return(as.integer(x))  # Should be negative IDs for leaves
    }
  }
  
  # Now, update left and right IDs in internal_nodes_df using the mapping
  internal_nodes_df$left <- sapply(internal_nodes_df$left, map_id)
  internal_nodes_df$right <- sapply(internal_nodes_df$right, map_id)
  
  # Build the merge matrix and height vector
  merge <- as.matrix(internal_nodes_df[, c("left", "right")])
  height <- internal_nodes_df$height
  
  # Compute order vector via depth-first traversal
  order <- integer(0)
  compute_order <- function(node) {
    if (is.null(node$children)) {
      leaf_index <- which(leaf_labels == node$name)
      order <<- c(order, leaf_index)
    } else {
      compute_order(node$children[[1]])
      compute_order(node$children[[2]])
    }
  }
  compute_order(tree)
  
  # Assemble the hclust object
  hclust_obj <- list(
    merge = merge,
    height = height,
    order = order,
    labels = leaf_labels,
    method = "Iterative Hierarchical Consensus Tree"
  )
  class(hclust_obj) <- "hclust"
  return(hclust_obj)
}
