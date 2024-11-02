iterative_consensus_tree <- function(net, 
                                     sites = unique(net[, 1]), 
                                     index = NULL,
                                     method = "average",
                                     depth = 1, 
                                     tree_structure = list(), 
                                     previous_height = Inf, 
                                     verbose = TRUE,
                                     n_runs = 100,
                                     max_remaining_size = length(sites)) {


  # Progress info because the algorithm is quite long to run
  if (verbose && interactive()) {
    cat(sprintf(
      "\rProcessing height: %.3f | Current branch size: %d | Max cluster size remaining: %d     ", 
      previous_height,
      length(sites),
      max(unlist(max_remaining_size))
    ))
    utils::flush.console()  # Ensures the message is displayed in real-time
  }

  
  # Step 1: species composition matrix for the current cluster
  current_net <- net[net[, 1] %in% sites, ]
  comat <- net_to_mat(current_net, weight = FALSE)
  
  # Step 2: dissimilarity matrix
  # TODO: implement the possibility of using other indices than the ones we have 
  # in the package already, e.g. for phylobeta diversity
  # TODO: implement correctly the formulas for abc and ABC
  dissim <- dissimilarity(comat, metric = index)
  dist_matrix <- stats::as.dist(
    net_to_mat(dissim, weight = TRUE, squared = TRUE, symmetrical = TRUE)
  )
  
  # Step 3: find the majority vote for a binary split among many trees
  subclusters <- stable_binary_split(dist_matrix, 
                                     method = method,
                                     n_runs = n_runs)
  
  # Ensure the smallest subcluster is the first (same rules as for hclust)
  if (length(subclusters[[1]]) > length(subclusters[[2]])) {
    subclusters <- list(subclusters[[2]], subclusters[[1]])
  }
  
  # Calculate height based on the specified method
  cluster1_sites <- subclusters[[1]]
  cluster2_sites <- subclusters[[2]]
  pairwise_distances <- as.matrix(dist_matrix)[cluster1_sites, cluster2_sites]
  
  # Calculate centroids for ward.D and ward.D2
  centroid1 <- colMeans(comat[cluster1_sites, , drop = FALSE])
  centroid2 <- colMeans(comat[cluster2_sites, , drop = FALSE])
  centroid_distance <- sum((centroid1 - centroid2)^2)
  
  # Compute height based on the selected method
  # TODO: check if ward calculation align with calculations in hclust &
  # fastclust
  calculated_height <- switch(
    method,
    "single" = min(pairwise_distances),
    "complete" = max(pairwise_distances),
    "average" = mean(pairwise_distances),
    "mcquitty" = mean(c(
      mean(as.matrix(dist_matrix)[cluster1_sites, ]),
      mean(as.matrix(dist_matrix)[cluster2_sites, ])
    )),
    "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "ward.D2" = 0.5 * (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "centroid" = centroid_distance,
    "median" = stats::median(pairwise_distances),
    stop("method argument is not valid")
  )
  
  # Step 4: height constraint (monotonic constraint of trees)
  height <- min(calculated_height, previous_height - 0.001)
  if (height < 0) height <- 0
  
  # update tree_structure for the current split
  tree_structure <- c(tree_structure, list(
    list(
      cluster = sites,
      subclusters = subclusters,
      height = height,
      depth = depth
    )
  ))
  
  # If we are the iteration where we work on the max cluster size,
  # then we update max_remaining_size based on the largest subcluster size
  # (always in position 2)
  if (length(sites) == max_remaining_size) {
    max_remaining_size <- length(subclusters[[2]])
  }
  
  # Next we process each subcluster
  for (subcluster in subclusters) {
    if (length(subcluster) == 1) {
      # If the subcluster is a single-site cluster, add it as a leaf node
      tree_structure <- c(tree_structure, list(list(leaf = subcluster)))
    } else {
      # Else if it is a normal cluster, recursively process it, passing the 
      # current height as the previous height
      tree_structure <- iterative_consensus_tree(
        net = net, 
        sites = subcluster, 
        index = index,
        method = method,  
        depth = depth + 1, 
        tree_structure = tree_structure, 
        previous_height = height, 
        verbose = verbose,
        max_remaining_size = max_remaining_size)
    }
  }




  # It works!! :)
  return(tree_structure)
}


# Function to find the majority vote for a binary split among many trees
stable_binary_split <- function(dist.obj,
                                method = "average", 
                                n_runs = 100) {
  
  # Step 1: randomize distance matrices and generate trees
  randomtrees <- list()
  for (run in 1:n_runs) {
    randomtrees[[run]] <- list()
    randomtrees[[run]]$dist.matrix <- .randomizeDistance(dist.obj)
    randomtrees[[run]]$hierartree <-  
      fastcluster::hclust(randomtrees[[run]]$dist.matrix, 
                          method = method)
  }
  
  # Step 2: cut trees at 2 clusters for each run
  trees <- lapply(randomtrees, function(trial) trial$hierartree)
  clusts <- lapply(trees, function(tree) {
    suppressMessages(cut_tree(
      stats::as.hclust(tree), n_clust = 2, find_h = FALSE))
  })
  
  # Step 3: ensure we have the same order of sites in each partition
  fixed_order <- sort(clusts[[1]]$ID) 
  df_list_ordered <- lapply(clusts, function(df) {
    df[match(fixed_order, df$ID), ]
  })
  
  # Step 4: make a data frame with all partitions to use the same functions as
  # in compare_partitions()
  partitions <- data.frame(ID = fixed_order)
  partitions <- cbind(partitions, lapply(df_list_ordered, `[[`, 2))
  colnames(partitions)[-1] <- paste0("run", seq_along(clusts))
  
  # Step 5: Calculate pairwise memberships for clustering consistency
  site_pw_memberships <- get_pairwise_membership(partitions[, -1])
  site_names <- do.call(rbind, strsplit(rownames(site_pw_memberships), "_"))
  
  # Step 6: Compute pairwise dissimilarity based on memberships
  # First as a network format as it is easier to handle
  pairwise_proportion <- data.frame(
    Site1 = fixed_order[as.numeric(site_names[, 1])],
    Site2 = fixed_order[as.numeric(site_names[, 2])],
    Weight = rowSums(!site_pw_memberships) / ncol(site_pw_memberships)
  ) # Note with use !site_pw_memberships, because we need dissimilarity here,
  # not similarity
  
  # Then we convert to dist object
  dist_pw_prop <- stats::as.dist(
    net_to_mat(pairwise_proportion, weight = TRUE, squared = TRUE, symmetrical = TRUE)
  )
  
  # Step 7: Create hierarchical tree from the pairwise dissimilarity matrix
  pw_tree <- fastcluster::hclust(dist_pw_prop, method = method)
  # Note we do not randomize here, but the tree topology is very conspicuous 
  # and so we probably do not need to randomize it
  groups <- stats::cutree(pw_tree, k = 2)
  groups <- data.frame(Site = names(groups), cluster = groups)
  
  return(split(groups$Site, groups$cluster))
}

# Function to make hclust object based on the output list we prepared above
# Function to make hclust object based on the output list we prepared above
reconstruct_hclust <- function(tree_structure) {
  
  # Separate internal nodes from leaf nodes
  internal_nodes <- Filter(function(x) !is.null(x$depth), tree_structure)
  
  # Initialize merge and height vectors + labels
  n_clusters <- length(unique(unlist(lapply(tree_structure, 
                                            function(x) if (!is.null(x$leaf)) x$leaf else
                                              x$cluster))))
  merge <- matrix(0, nrow = n_clusters - 1, ncol = 2)  # Correct matrix size
  height <- numeric(n_clusters - 1)
  labels <- unique(unlist(lapply(tree_structure, 
                                 function(x) if (!is.null(x$leaf)) x$leaf else
                                   x$cluster)))
  
  # Map leaf nodes to negative IDs likewise to hclust
  leaf_map <- stats::setNames(-seq_along(labels), labels)
  cluster_id <- 1  # Start IDs for clusters at 1
  
  # Map to track each cluster ID
  cluster_map <- leaf_map  # Start with the map containing leaves
  
  # Sort internal nodes by increasing height to ensure monotonic height order
  internal_nodes <- internal_nodes[order(sapply(internal_nodes, function(x) x$height))]
  
 # Process each internal node bottom up
  for (node in internal_nodes) {
    
    # Process left and right subclusters
    left <- node$subclusters[[1]]
    right <- node$subclusters[[2]]
    
    # Ensure IDs for left and right subclusters are available in cluster_map
    left_key <- paste(left, collapse = "-")
    right_key <- paste(right, collapse = "-")
    
    # Assign a new cluster ID if the subcluster key is missing
    if (!left_key %in% names(cluster_map)) {
      cluster_map[[left_key]] <- cluster_id
      cluster_id <- cluster_id + 1
    }
    if (!right_key %in% names(cluster_map)) {
      cluster_map[[right_key]] <- cluster_id
      cluster_id <- cluster_id + 1
    }
    
    # Get IDs for left and right subclusters
    left_id <- cluster_map[[left_key]]
    right_id <- cluster_map[[right_key]]
    
    # Assign a new cluster ID for the current node if it’s not already assigned
    current_key <- paste(node$cluster, collapse = "-")
    if (!current_key %in% names(cluster_map)) {
      cluster_map[[current_key]] <- cluster_id
      cluster_id <- cluster_id + 1
    }
    
    # Add to merge and height
    merge[cluster_map[[current_key]], ] <- c(left_id, right_id)
    height[cluster_map[[current_key]]] <- node$height
  }
  
  # Compute order iteratively based on left and right structure
  order <- c()
  stack <- list(tree_structure[[1]])  # Start with the root node
  
  while (length(stack) > 0) {
    node <- stack[[1]]
    stack <- stack[-1]  # Remove first node at every iteration
    
    if (!is.null(node$leaf)) {
      # If it's a leaf, add it to the order
      order <- c(order, match(node$leaf, labels))
    } else if (!is.null(node$cluster)) {
      # If it's a cluster, determine if left and right are leaves or clusters
      if (length(node$subclusters[[2]]) == 1) {
        # Right subcluster is a leaf
        right_id <- match(node$subclusters[[2]], labels)
        order <- c(order, right_id)
      } else {
        # Right subcluster is a cluster, add to stack
        right_index <- which(
          sapply(internal_nodes,
                 function(x) identical(x$cluster, 
                                       node$subclusters[[2]])))
        if (length(right_index) == 1) stack <-
            c(list(internal_nodes[[right_index]]), stack)
      }
      
      if (length(node$subclusters[[1]]) == 1) {
        # Left subcluster is a leaf
        left_id <- match(node$subclusters[[1]], labels)
        order <- c(order, left_id)
      } else {
        # Left subcluster is a cluster, add to stack
        left_index <- which(
          sapply(internal_nodes, function(x) identical(x$cluster, 
                                                       node$subclusters[[1]])))
        if (length(left_index) == 1) stack <- 
            c(list(internal_nodes[[left_index]]), stack)
      }
    }
  }
  
  # Construct our beautiful hclust object
  hc <- list(
    merge = merge,
    height = height,
    order = order,
    labels = labels,
    method = "Iterative Hierarchical Consensus Tree"
  )
  class(hc) <- "hclust"
  
  # tadam!
  return(hc)
}


