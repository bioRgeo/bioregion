IHCT <- function(dist_mat,
                 sites = rownames(dist_mat),
                 method = "average",
                 n_runs = 100,
                 top_n_trees = 2,
                 depth = 1,
                 previous_height = Inf,
                 max_remaining_size = length(sites),
                 monotonicity_direction = c("top-down", "bottom-up"),
                 stable_shortcircuit = FALSE,
                 stability_check = "top_n_trees",
                 verbose = TRUE) {
  
  # n_runs & top_n_trees
  if(top_n_trees > n_runs){
    top_n_trees <- n_runs
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
  split_stable <- FALSE
  split_trees <- NULL
  if(nrow(dist_mat_d) > 2) { # Only if we have more than 2 sites
    subclusters <- coassign_binary_split(dist_mat_d,
                                         method = method,
                                         n_runs = n_runs,
                                         top_n_trees = top_n_trees,
                                         binsplit = "tree",
                                         stability_check = if (stable_shortcircuit) stability_check else "none")
    if (stable_shortcircuit) {
      split_stable <- isTRUE(attr(subclusters, "stable"))
      split_trees <- attr(subclusters, "trees")
    }
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

    # Check if we can short-circuit this subcluster
    use_shortcircuit <- FALSE
    if (stable_shortcircuit && split_stable &&
        length(subcluster) > 2 && !is.null(split_trees)) {
      if (subtrees_identical(split_trees, subcluster)) {
        use_shortcircuit <- TRUE
      }
    }

    if (use_shortcircuit) {
      # Extract the sub-tree from the best tree (first in sorted list,
      # as they are all identical anyway)
      dend <- stats::as.dendrogram(split_trees[[1]])
      branch <- find_branch_for_sites(dend, subcluster)
      child_structure <- dend_to_IHCT_tree(branch, depth + 1,
                                           dist_mat_d, method)
      child_structure <- apply_monotonicity(child_structure, height,
                                            monotonicity_direction)
    } else {
      child_structure <- IHCT(dist_mat_d,
                              sites = subcluster,
                              method = method,
                              n_runs = n_runs,
                              top_n_trees = top_n_trees,
                              depth = depth + 1,
                              previous_height = height,
                              max_remaining_size = max_remaining_size,
                              monotonicity_direction = monotonicity_direction,
                              stable_shortcircuit = stable_shortcircuit,
                              stability_check = stability_check,
                              verbose = verbose)
    }
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
                                  top_n_trees = 2,
                                  binsplit = "tree",
                                  stability_check = "none") {

  # Number of sites
  n <- dim(dist_mat)[1]

  # Step 1: randomize distance matrices and generate trees
  all_trees <- list()
  ccc <- NULL
  for (run in 1:n_runs) {
    all_trees[[run]] <- list()
    rand_dist_mat <- stats::as.dist(randomize_dist(dist_mat))
    all_trees[[run]] <- fastcluster::hclust(rand_dist_mat,
                                            method = method)
    ccc <- c(ccc, suppressWarnings(tree_eval(all_trees[[run]],
                                             dist_mat)$cophcor))
  }

  # Sort by cophenetic correlation and select best trees
  sort_order <- order(ccc, decreasing = TRUE)
  trees <- all_trees[sort_order[1:top_n_trees]]

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

  # Step 5: Split into two clusters from the coassign matrix
  if (binsplit == "tree") {
    pw_tree <- fastcluster::hclust(coassign, method = method)
    groups <- stats::cutree(pw_tree, k = 2)
  } else if (binsplit == "pam") {
    groups <- cluster::pam(coassign, k = 2)$clustering
  }
  groups <- data.frame(Site = names(groups), cluster = groups)

  result <- split(groups$Site, groups$cluster)

  # Step 6: Stability detection (when requested)
  # This is a first filter only // gate 1:
  # we first check if the splits appear stable (based on clusters)
  # if not stable, there is no reason to test tree topology
  # Basically, the question asked in this first gate is:
  # "do all trees agree on which sites go into which group?"
  if (stability_check != "none") {
    if (stability_check == "all") {
      # Check stability across ALL n_runs trees
      all_cuts <- lapply(all_trees, function(tree) {
        stats::cutree(tree, k = 2)
      })
      # Normalize labels: first site always gets label 1
      ref_site <- names(all_cuts[[1]])[1]
      all_assignments <- vapply(all_cuts, function(cut) {
        if (cut[ref_site] != 1L) 3L - cut else cut
      }, integer(n))
      ref_col <- all_assignments[, 1]
      split_stable <- all(vapply(2:ncol(all_assignments), function(j) {
        identical(all_assignments[, j], ref_col)
      }, logical(1)))
    } else {
      # stability_check == "top_n_trees": check only the best trees
      # cut_tree/knbclu returns character cluster labels (e.g. "1", "2"),
      # so we normalize by swapping labels when the first site disagrees,
      # using character replacement rather than integer arithmetic.
      cluster_labels <- sort(unique(clusts[, 1]))
      norm_clusts <- clusts
      ref_label <- norm_clusts[1, 1]
      for (col_idx in seq_len(ncol(norm_clusts))) {
        if (norm_clusts[1, col_idx] != ref_label) {
          col_vals <- norm_clusts[, col_idx]
          norm_clusts[, col_idx] <- ifelse(
            col_vals == cluster_labels[1],
            cluster_labels[2],
            cluster_labels[1])
        }
      }
      ref_col <- norm_clusts[, 1]
      split_stable <- if (ncol(norm_clusts) > 1) {
        all(vapply(2:ncol(norm_clusts), function(j) {
          identical(norm_clusts[, j], ref_col)
        }, logical(1)))
      } else {
        TRUE
      }
    }
    attr(result, "stable") <- split_stable
    attr(result, "trees") <- trees
  }

  return(result)
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


# Helper functions for stable sub-tree short-circuit optimization 

# We ask: 
# for a given subcluster, do all trees have the same branching structure?
# Process:
# 1. convert to dendrogram
# 2. extract branch from the dendrogram containing only our sites
# 3. sort tree children alphabetically so that equivalent trees will
# always provide similar results, even if order of branch is not correct



# Check whether all hclust trees have identical sub-tree topology
# for a given subcluster
subtrees_identical <- function(trees, subcluster_sites) {
  canonicals <- vapply(trees, function(hc) {
    # 1. 
    dend <- stats::as.dendrogram(hc)
    # 2.
    branch <- find_branch_for_sites(dend, subcluster_sites)
    # 3.
    canonical_subtree(branch)
  }, character(1))
  length(unique(canonicals)) == 1
}

# Given a dendrogram (from an hclust tree cut at k=2), find the branch
# whose leaves match the specified sites
find_branch_for_sites <- function(dend, sites) {
  b1_labels <- sort(get_dend_leaves(dend[[1]]))
  target <- sort(sites)
  if (identical(b1_labels, target)) return(dend[[1]])
  return(dend[[2]])
}

# Get all leaf labels from a dendrogram
get_dend_leaves <- function(dend) {
  if (is.leaf(dend)) return(attr(dend, "label"))
  unlist(lapply(seq_along(dend), function(i) get_dend_leaves(dend[[i]])))
}



# Returns a canonical string representation of a dendrogram topology.
# Children are sorted by their canonical string so that two topologically
# identical trees produce the same output regardless of internal child order.
canonical_subtree <- function(dend) {
  if (is.leaf(dend)) return(attr(dend, "label"))
  child_strs <- vapply(seq_along(dend), function(i) {
    canonical_subtree(dend[[i]])
  }, character(1))
  child_strs <- sort(child_strs)
  paste0("(", paste(child_strs, collapse = ","), ")")
}



# Count leaves in an IHCT tree node
count_IHCT_leaves <- function(node) {
  if (is.null(node$children)) return(1L)
  sum(vapply(node$children, count_IHCT_leaves, integer(1)))
}

# Compute IHCT-style height for a binary split
compute_IHCT_height <- function(cluster1_sites, cluster2_sites,
                                dist_mat, method) {
  pairwise_distances <- dist_mat[cluster1_sites, cluster2_sites]

  centroid1 <- colMeans(dist_mat[cluster1_sites, , drop = FALSE])
  centroid2 <- colMeans(dist_mat[cluster2_sites, , drop = FALSE])
  centroid_distance <- sum((centroid1 - centroid2)^2)

  switch(
    method,
    "single" = min(pairwise_distances),
    "complete" = max(pairwise_distances),
    "average" = mean(pairwise_distances),
    "mcquitty" = mean(pairwise_distances),
    "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) /
      (length(cluster1_sites) + length(cluster2_sites)) * sqrt(centroid_distance),
    "ward.D2" = (length(cluster1_sites) * length(cluster2_sites)) /
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "centroid" = centroid_distance,
    "median" = 0.5 * centroid_distance,
    stop("method argument is not valid")
  )
}

# Convert a dendrogram into the IHCT list format (name, height, depth, children).
# Recalculates heights using IHCT's formulas rather than hclust's stored heights.
dend_to_IHCT_tree <- function(dend, depth, dist_mat, method) {
  if (is.leaf(dend)) {
    return(list(name = attr(dend, "label"), height = 0,
                depth = depth, children = NULL))
  }

  child1 <- dend_to_IHCT_tree(dend[[1]], depth + 1, dist_mat, method)
  child2 <- dend_to_IHCT_tree(dend[[2]], depth + 1, dist_mat, method)

  # Ensure smallest subcluster first (IHCT convention)
  n1 <- count_IHCT_leaves(child1)
  n2 <- count_IHCT_leaves(child2)
  if (n1 > n2) {
    children <- list(child2, child1)
  } else {
    children <- list(child1, child2)
  }

  # Get leaf labels for each branch to compute height
  c1_leaves <- get_dend_leaves(dend[[1]])
  c2_leaves <- get_dend_leaves(dend[[2]])
  all_leaves <- c(c1_leaves, c2_leaves)

  sub_dist <- dist_mat[all_leaves, all_leaves]
  calculated_height <- compute_IHCT_height(c1_leaves, c2_leaves, sub_dist, method)
  if (calculated_height < 0) calculated_height <- 0

  return(list(
    name = paste(all_leaves, collapse = ","),
    height = calculated_height,
    depth = depth,
    children = children
  ))
}

# Apply monotonicity constraints to a converted sub-tree
apply_monotonicity <- function(node, previous_height, direction) {
  if (is.null(node$children)) return(node)

  if (direction == "top-down") {
    node$height <- min(node$height, previous_height)
    if (node$height < 0) node$height <- 0
    node$children <- lapply(node$children, function(child) {
      apply_monotonicity(child, node$height, direction)
    })
  } else if (direction == "bottom-up") {
    # First recurse into children
    node$children <- lapply(node$children, function(child) {
      apply_monotonicity(child, Inf, direction)
    })
    # Then adjust height upward if needed
    max_child_height <- max(vapply(node$children,
                                   function(x) x$height, numeric(1)))
    if (max_child_height > node$height) {
      node$height <- max_child_height
    }
  }

  return(node)
}
